"""GCE normalization/statistics and shared propagation choices."""
import io
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from indirect import cr_config
from indirect import dbar_interface
from indirect import gce_interface
from indirect import pbar_interface
from plotting import plot_ams02_pbar_fluxes as pbar_plot
from indirect.pbar.core import FluxTable
from test_pbar_integration import Harness,record,ROOT,ENV,RUN,SOURCE
import ast

node=next(n for n in RUN.body if isinstance(n,ast.FunctionDef) and n.name=='calculate_gce_fit')
exec(compile(ast.Module(body=[node],type_ignores=[]),str(SOURCE),'exec'),ENV)
Harness.calculate_gce_fit=ENV['calculate_gce_fit']


from flux_test_data import requires_both, requires_dbar


class GCETests(unittest.TestCase):
    def test_first_block_selection_units_and_upper_limit_sentinels(self):
        energy,flux,used=gce_interface.read_data()
        self.assertEqual((used.sum(),(~used).sum()),(11,15))
        self.assertEqual(energy.shape,(26,3))
        self.assertAlmostEqual(flux[0,0],3.860e-4*1e-3/.4775494185)
        self.assertAlmostEqual(energy[0,0],.594607)
        self.assertTrue(np.all(flux[~used,0]<0))
        self.assertTrue(np.all(flux[used,3]<flux[used,0]))

    def test_split_normal_fit_and_delta_chi2_one_errors(self):
        p=np.array([1.,2.,3.]);y=2*p
        fit=gce_interface.fit_amplitude(p,y,np.full(3,.1),np.full(3,.2))
        self.assertAlmostEqual(fit['amplitude'],2.,places=11)
        self.assertAlmostEqual(2-fit['amplitude_low'],.1/np.sqrt(14),places=11)
        self.assertAlmostEqual(fit['amplitude_high']-2,.2/np.sqrt(14),places=11)
        self.assertAlmostEqual(fit['chi2'],0.,places=18)
        self.assertEqual(fit['ndof'],2)

    def test_geometry_and_dirac_scaling_of_best_fit_and_errors(self):
        engine=gce_interface.GCEFit()
        majorana=engine.evaluate(50.,[record([5,-5],3.e-26)],1000022,4,{})
        dirac=engine.evaluate(50.,[record([5,-5],3.e-26,initial=[-1000022,1000022])],1000022,8,{})
        doubled=engine.evaluate(50.,[record([5,-5],3.e-26)],1000022,4,{'gce_jbar':2*majorana['geometry']['jbar']})
        self.assertEqual(majorana['ndof'],10)
        self.assertEqual(majorana['geometry']['denominator_pi'],8)
        self.assertEqual(dirac['geometry']['denominator_pi'],16)
        for key in ('sigmav_best_cm3s','error_minus_cm3s','error_plus_cm3s'):
            self.assertAlmostEqual(dirac[key]/majorana[key],2.,places=9)
            self.assertAlmostEqual(doubled[key]/majorana[key],.5,places=9)
        self.assertAlmostEqual(majorana['chi2'],dirac['chi2'],places=10)
        self.assertEqual(len(majorana['upper_limit_checks']),15)

    def test_mass_interpolation_and_kinematic_support(self):
        engine=gce_interface.GCEFit()
        mass=55.;index=engine.channels.index('bb')
        lo=np.flatnonzero(engine.masses==50.)[0];hi=np.flatnonzero(engine.masses==60.)[0]
        native=.5*(engine.yields[lo,:,index]+engine.yields[hi,:,index])
        f=engine.spectrum(mass,{'bb':1.})
        e=mass*engine.x
        np.testing.assert_allclose(f(e),native/(e*np.log(10.)),rtol=1e-12,atol=1e-100)
        np.testing.assert_array_equal(f([mass*1.01,mass*2]),[0.,0.])

    def test_actual_multibody_spectrum_supported_without_fake_channel_mapping(self):
        engine=gce_interface.GCEFit();rows=[record([11,-11,13,-13],3.e-26,decay=True)]
        missing=engine.evaluate(100.,rows,1000022,4,{})
        self.assertEqual(missing['status'],'unsupported_channels')
        valid=engine.evaluate(100.,rows,1000022,4,{},runtime_spectrum=lambda e:np.exp(-np.asarray(e)/5)/np.asarray(e))
        self.assertEqual(valid['status'],'ok')
        self.assertIn('runtime',valid['source'])

    def test_actual_hook_reuses_mw_rates_and_writes_complete_products(self):
        with tempfile.TemporaryDirectory() as directory:
            h=Harness(directory);h.mode={'indirect':'sigmav'}
            h._milky_way_rates=[record([5,-5],3.e-26)]
            h.maddm_card['gce_plot']=False
            h._collect_pbar_rates=lambda *a:(_ for _ in ()).throw(AssertionError('extra rate run'))
            h.calculate_gce_fit(50.)
            self.assertEqual(h.last_results['GCE_status'],'ok')
            out=Path(directory)/'output/run_01'
            self.assertEqual(np.loadtxt(out/'GCE_spectrum.dat').shape,(26,11))
            self.assertTrue((out/'GCE_fit.txt').exists())
            saved=json.loads((out/'GCE_fit.json').read_text())
            self.assertNotIn('_bins',saved)


class CRTests(unittest.TestCase):
    @requires_both
    def test_all_54_tables_and_metadata_are_resolved(self):
        for scheme,setup,halo,sp in itertools.product(['QUAINT','BIG','SLIM'],['MIN','MED','MAX'],['NFW','Einasto','Burkert'],['pbar','dbar']):
            card=dict(cr_propagation_scheme=scheme,cr_propagation_set=setup,cr_halo_profile=halo)
            table=FluxTable(cr_config.table_path(sp,card))
            self.assertEqual(table.propagation,scheme+'-'+setup)
            self.assertEqual(table.halo,halo)
            self.assertEqual(len(table.masses),64)

    def test_secondary_selection_shares_only_identical_big_and_slim_inputs(self):
        q=cr_config.secondary_path({})
        big=cr_config.secondary_path({'cr_propagation_scheme':'BIG'})
        self.assertNotEqual(q,big)
        slim=cr_config.secondary_path({'cr_propagation_scheme':'SLIM'})
        self.assertTrue(slim.is_file())
        self.assertEqual(slim,big)
        self.assertEqual(cr_config.secondary_path({'cr_propagation_scheme':'SLIM','pbar_secondary_file':str(q)}),q)

    def test_native_limits_only_match_the_original_configuration(self):
        self.assertTrue(cr_config.native_channel_limits({}))
        for card in [{'cr_halo_profile':'Einasto'},{'cr_propagation_set':'MIN'},
                     {'cr_propagation_scheme':'BIG'},{'pbar_phi_max':1.}]:
            self.assertFalse(cr_config.native_channel_limits(card))

    @requires_dbar
    def test_dbar_forecast_and_metadata_follow_changed_settings(self):
        engine=dbar_interface.AntideuteronEvents();rows=[record([5,-5],3.e-26)]
        default=engine.evaluate(100.,rows,1000022,4,1.,{})
        other=engine.evaluate(100.,rows,1000022,4,1.,dict(cr_halo_profile='Burkert',cr_propagation_scheme='SLIM',cr_propagation_set='MAX'))
        self.assertEqual(other['propagation'],'SLIM-MAX');self.assertEqual(other['halo'],'Burkert')
        self.assertIn('secondary_BIG_IS.txt',other['secondary_path'])
        self.assertNotAlmostEqual(other['cases']['all_dm']['AMS02_total_events']/default['cases']['all_dm']['AMS02_total_events'],1.,places=4)

    def test_upper_configuration_plot_and_separate_model_curve(self):
        fixture=ROOT/'tests/fixtures/pbar/AMS02_pbar_fluxes.dat'
        import shutil
        with tempfile.TemporaryDirectory() as directory:
            shutil.copy2(fixture,directory)
            with patch.object(pbar_plot,'make_plot') as plot:
                pbar_interface.save_plot(directory)
            args,kwargs=plot.call_args
            a=np.loadtxt(fixture)
            np.testing.assert_array_equal(args[0]['total'],a[:,14])
            np.testing.assert_array_equal(kwargs['model_dm'],a[:,4])
            self.assertEqual(kwargs['prediction'],'upper')


if __name__=='__main__':unittest.main()
