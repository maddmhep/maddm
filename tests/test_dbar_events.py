"""Dbar forecasts: energy/flux Jacobians, sensitivity normalization and run hooks."""
import ast
import io
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from indirect import dbar_interface as dbar
from indirect import pbar_interface
from plotting import plot_ams02_pbar_fluxes as pbar_plot
from test_pbar_integration import Harness, record, ROOT


from flux_test_data import requires_dbar


class ForecastTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = dbar.AntideuteronEvents()

    def evaluate(self, records=None, card=None, xsi=1., md=4, pbar=None):
        from indirect import cr_config
        try:
            cr_config.table_path('dbar', card or {})
        except ValueError as error:
            self.skipTest(str(error))
        rows = records if records is not None else [record([5, -5], 3.e-26)]
        return self.engine.evaluate(100., rows, 1000022, md, xsi, card or {}, pbar)

    def test_constant_sensitivity_flux_produces_two_events_in_each_window(self):
        for bounds in dbar.WINDOWS.values():
            low, high, sensitivity = bounds
            result = dbar.window_forecast(lambda t: sensitivity, low, high, sensitivity)
            self.assertAlmostEqual(result['expected_events'], 2., places=13)
            self.assertAlmostEqual(result['integrated_flux']/sensitivity, high-low, places=13)

    def test_total_energy_to_per_nucleon_jacobian_and_force_field(self):
        t = np.array([.18, .7, 2., 4.2])
        template = lambda k: 1./np.asarray(k)**2
        # With phi=0, dPhi/dT=2*dPhi/dK at K=2T, so 1/(2T^2).
        np.testing.assert_allclose(dbar.flux_per_nucleon(template, t, 3.e-26, 0.), 1./(2*t*t))
        phi = .7
        k, u = 2*t, 2*t+phi
        expected = 2./u**2*k*(k+2*1.875613)/(u*(u+2*1.875613))
        np.testing.assert_allclose(dbar.flux_per_nucleon(template, t, 3.e-26, phi), expected)
        analytic = dbar.window_forecast(lambda t: 1./(2*t*t), .18, .7, 8.e-7)
        expected_N = (1/.18-1/.7)/((.7-.18)*8.e-7)
        self.assertAlmostEqual(analytic['expected_events']/expected_N, 1., places=12)

    def test_native_table_and_independent_dense_window_integration(self):
        result = self.evaluate(card={'dbar_phi': .7})
        table = self.engine.table
        self.assertEqual((len(table.masses), table.flux.shape[2]), (64, 19))
        provenance = json.loads((dbar.DATA_DIR/'provenance.json').read_text())
        self.assertIn('total kinetic', provenance['energy'])
        template = table.template(100., {'bb': 1.})
        for name, (low, high, sensitivity) in dbar.WINDOWS.items():
            grid = np.linspace(low, high, 20001)
            flux = dbar.flux_per_nucleon(template, grid, 3.e-26, .7)
            trapezoid = getattr(np, 'trapezoid', None) or np.trapz
            independent = 2*trapezoid(flux, grid)/((high-low)*sensitivity)
            self.assertLess(abs(result['cases']['all_dm']['windows'][name]['expected_events']/independent-1), 1.e-7)

    def test_branching_rate_density_and_thermal_linearity(self):
        bb = self.evaluate()
        ww = self.evaluate([record([24, -24], 3.e-26)])
        mix = self.evaluate([record([5, -5], .9e-26), record([24, -24], 2.1e-26)])
        scale = self.evaluate([record([5, -5], 6.e-26)], {'rhoDM': .8}, xsi=.5)
        for name in dbar.WINDOWS:
            value = lambda r, case='all_dm': r['cases'][case]['windows'][name]['expected_events']
            self.assertAlmostEqual(value(mix)/(0.3*value(bb)+.7*value(ww)), 1., places=11)
            self.assertAlmostEqual(value(scale)/value(bb), 8., places=11)
            self.assertAlmostEqual(value(scale, 'thermal')/value(scale), .25, places=11)

    def test_dirac_source_factor_is_the_same_as_pbar(self):
        majorana = self.evaluate()
        rows = [record([5,-5], 3.e-26, initial=(-1000022,1000022))]
        dirac = self.evaluate(rows, md=8)
        alternate = self.evaluate(rows, md=8, card={'pbar_table_convention':'particle_antiparticle'})
        m, d, a = [r['cases']['all_dm']['AMS02_total_events'] for r in (majorana,dirac,alternate)]
        self.assertAlmostEqual(d/m, .5, places=12)
        self.assertAlmostEqual(a/m, 1., places=12)

    def test_modulation_uses_model_point_fits_and_explicit_fixed_fallback(self):
        pb = dict(status='ok', model_points={'all_dm':{'phi_GV':1.1},'thermal':{'phi_GV':.9}},
                  upper_limit={'phi_GV':0.}, best_fit={'phi_GV':.2})
        result = self.evaluate(xsi=.5, pbar=pb)
        self.assertEqual(result['cases']['all_dm']['phi_GV'], 1.1)
        self.assertEqual(result['cases']['thermal']['phi_GV'], .9)
        fixed = self.evaluate(card={'dbar_phi':.3}, pbar={'status':'tabulated_channels'})
        self.assertEqual(fixed['cases']['all_dm']['phi_GV'], .3)
        self.assertIn('fixed', fixed['solar_modulation_source'])

    def test_unsupported_rate_is_never_silently_dropped_and_zero_rate_is_valid(self):
        result = self.evaluate([record([5,-5], 3.e-26),record([9999,9999], 1.e-25)])
        self.assertEqual(result['status'], 'unsupported_channels')
        self.assertNotIn('cases', result)
        zero = self.evaluate([record([5,-5], 0.)])
        self.assertTrue(zero['zero_model_rate'])
        self.assertEqual(zero['cases']['all_dm']['AMS02_total_events'], 0.)
        for mass in (4., 100001.):
            with self.assertRaises(ValueError):
                self.engine.evaluate(mass,[record([5,-5],3.e-26)],1000022,4,1.,{})

    def test_saved_forecast_and_flux_units_and_stale_files(self):
        result = self.evaluate(xsi=.5)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            dbar.write_outputs(out,result)
            saved = json.loads((out/'AMS02_dbar_events.json').read_text())
            self.assertNotIn('_flux_columns', saved)
            self.assertIn('(GeV/n)^-1',saved['flux_units'])
            rows = np.loadtxt(out/'AMS02_dbar_events.dat',dtype=str)
            self.assertEqual(rows.shape,(6,7))
            values = dbar.observables(result)
            self.assertAlmostEqual(values['AMS02_dbar_events_total'], values['AMS02_dbar_events_low']+values['AMS02_dbar_events_high'])
            text = io.StringIO(); dbar.save_summary(text,result)
            self.assertIn('AMS02_dbar_events_total',text.getvalue())
            dbar.write_outputs(out,{'status':'disabled'})
            self.assertFalse((out/'AMS02_dbar_fluxes.dat').exists())
            self.assertFalse((out/'AMS02_dbar_events.dat').exists())


class HookTests(unittest.TestCase):
    @requires_dbar
    def test_dbar_reuses_pbar_rates_and_does_not_fit_a_sigmav_flux(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory)
            h._milky_way_rates = [record([5,-5],3.e-26)]
            h._pbar_result = dict(status='tabulated_channels')
            h._collect_pbar_rates = lambda *a: (_ for _ in ()).throw(AssertionError('duplicate rate run'))
            h.calculate_dbar_events(100.)
            self.assertEqual(h.last_results['AMS02_dbar_status'],'ok')
            self.assertEqual(h.last_results['AMS02_dbar_phi'],.7)
            self.assertTrue((Path(directory)/'output/run_01/AMS02_dbar_events.json').exists())

    @requires_dbar
    def test_dbar_works_with_pbar_limits_disabled(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory)
            h.maddm_card['pbar_limits'] = False
            h.calculate_pbar_limits(100.)
            h.calculate_dbar_events(100.)
            self.assertEqual(h.last_results['AMS02_pbar_status'],'disabled')
            self.assertEqual(h.last_results['AMS02_dbar_status'],'ok')
            self.assertEqual(h._dbar_result['rates'][0]['source'],'inclusive_Milky_Way')

    def test_supplied_plot_settings_and_pdf_save_hook(self):
        import matplotlib.figure
        original = matplotlib.figure.Figure.savefig
        observed = {}
        def capture(fig,*args,**options):
            top,bottom = fig.axes
            observed.update(x=top.get_xscale(), y=top.get_yscale(), residual=bottom.get_yscale(),
                            xmax=bottom.get_xlim()[1], dpi=options.get('dpi'))
            return original(fig,*args,**options)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            shutil.copy2(ROOT/'tests/fixtures/pbar/AMS02_pbar_fluxes.dat',out)
            with patch.object(matplotlib.figure.Figure,'savefig',capture):
                saved = pbar_interface.save_plot(out)
            self.assertTrue(Path(saved).read_bytes().startswith(b'%PDF'))
            self.assertEqual(observed,dict(x='log',y='log',residual='linear',xmax=1000.,dpi=600))
            self.assertIsNone(pbar_interface.save_plot(out,enabled=False))
            self.assertFalse(Path(saved).exists())


if __name__ == '__main__':
    unittest.main()
