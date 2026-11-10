"""Physics checks for the shared halo and the likelihood-independent Earth fluxes."""
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad, trapezoid

from indirect import cr_config
from indirect import cr_earth
from indirect import dbar_interface
from indirect import halo_model
from indirect.pbar.core import FluxTable, Likelihood, SIGMAV_REF


from flux_test_data import requires_both, requires_dbar, requires_pbar


class HaloTests(unittest.TestCase):
    def test_pair_average_has_correct_p_wave_moment(self):
        for v0, escape in ((238., 544.), (175., 400.), (220., 125.)):
            with self.subTest(v0=v0, escape=escape):
                rule = np.array(halo_model.relative_speed_rule(v0, escape))
                u, w = rule.T
                z = escape/v0
                # Isotropy implies <|v1-v2|²> = 2 <v²> for independent draws.
                moment = 2*v0*v0*quad(lambda x: x**4*np.exp(-x*x), 0., z)[0]/quad(
                    lambda x: x*x*np.exp(-x*x), 0., z)[0]
                self.assertAlmostEqual(w.sum(), 1., places=15)
                self.assertTrue(np.all((u > 0.) & (u*halo_model.C_KMS < 2*escape)))
                self.assertAlmostEqual(np.dot(w, (u*halo_model.C_KMS)**2)/moment, 1., places=12)

    def test_inverse_speed_is_normalized_with_and_without_lab_boost(self):
        v0, escape = 238., 544.
        for boost in (0., 250., 600.):
            with self.subTest(boost=boost):
                eta = lambda v: halo_model.mean_inverse_speed(v, v0, escape, boost)
                self.assertAlmostEqual(quad(eta, 0., escape+boost, points=[abs(escape-boost)])[0], 1., places=11)
                self.assertEqual(eta(escape+boost+1.), 0.)

    def test_rapidd_table_uses_configured_motion_and_vendor_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            default = Path(tmp)/'vendor.dat'
            default.write_text('200 2\n')
            path = Path(tmp)/'SHM.dat'
            halo_model.write_recoil_halo(path, {}, default)
            self.assertEqual(path.read_text().splitlines()[0], '200 0')
            data = np.loadtxt(path, skiprows=1)
            self.assertEqual(data.shape, (200, 2))
            self.assertAlmostEqual(trapezoid(data[:, 1], data[:, 0]), 1., places=4)
            halo_model.write_recoil_halo(path, {'vEarth_mod': 0.}, default)
            self.assertFalse(np.array_equal(data, np.loadtxt(path, skiprows=1)))


class EarthFluxTests(unittest.TestCase):
    mixture = dict(status='ok', total_sigmav_cm3s=SIGMAV_REF, branching={'bb': 1.})
    fitted = dict(status='ok', model_points=dict(all_dm=dict(phi_GV=.5, norm=.83, slope=.07),
                                               thermal=dict(phi_GV=.9, norm=1.05, slope=-.07)))

    @requires_pbar
    def test_pbar_components_agree_with_likelihood_at_measured_points(self):
        engine = cr_earth.EarthFluxes()
        # reference checks at the table normalization density (0.4 GeV/cm^3)
        result = engine.evaluate('pbar', 100., self.mixture, 4, .5, {'rhoDM': .4}, self.fitted)
        grid, sec_is, dm_is, total_is, sec_toa, dm_toa, total_toa = result['_columns'][:, :7].T
        table = FluxTable(cr_config.table_path('pbar', {}))
        likelihood = Likelihood(secondary_path=cr_config.secondary_path({}))
        sec, dm = likelihood.components(table.template(100., rho_local_GeVcm3=.4), SIGMAV_REF, .83, .07, .5)
        # The Earth file is per energy; the likelihood uses per rigidity.
        for expected, actual in ((sec, sec_toa), (dm, dm_toa)):
            select = (likelihood.t > grid[0]) & (likelihood.t < min(grid[-1], 30.))
            calculated = np.interp(np.log(likelihood.t[select]), np.log(grid), actual)
            np.testing.assert_allclose(calculated*likelihood.jacobian[select], expected[select], rtol=.002)
        np.testing.assert_array_equal(total_is, sec_is+dm_is)
        np.testing.assert_array_equal(total_toa, sec_toa+dm_toa)

    @requires_dbar
    def test_dbar_per_nucleon_and_shared_density(self):
        engine = cr_earth.EarthFluxes()
        settings = dict(cr_propagation_scheme='SLIM', cr_propagation_set='MAX', dbar_phi=.7, rhoDM=.4)
        result = engine.evaluate('dbar', 100., self.mixture, 4, .5, settings)
        grid, dm_is, dm_toa = result['_columns'][:, [0, 2, 5]].T
        template = FluxTable(cr_config.table_path('dbar', settings)).template(100., rho_local_GeVcm3=.4)
        np.testing.assert_allclose(dm_is, dbar_interface.flux_per_nucleon(template, grid, SIGMAV_REF, 0.))
        np.testing.assert_allclose(dm_toa, dbar_interface.flux_per_nucleon(template, grid, SIGMAV_REF, .7))
        scaled = engine.evaluate('dbar', 100., self.mixture, 8, .5, dict(settings, rhoDM=.8))
        np.testing.assert_allclose(scaled['_columns'][:, 5], 2*dm_toa)
        np.testing.assert_allclose(result['_columns'][:, 11], .25*dm_toa)
        self.assertIn('secondary_BIG_IS.txt', result['secondary_table'])

    @requires_both
    def test_missing_positron_tables_do_not_enter_active_outputs(self):
        engine = cr_earth.EarthFluxes()
        with patch.object(cr_config, 'table_path', wraps=cr_config.table_path) as paths:
            results = {s: engine.evaluate(s, 100., self.mixture, 4, 1., {}) for s in cr_earth.ACTIVE}
        self.assertEqual([call.args[0] for call in paths.call_args_list], ['pbar', 'dbar'])
        with tempfile.TemporaryDirectory() as tmp:
            summary = cr_earth.write_outputs(tmp, results)
            self.assertEqual(summary['positron']['status'], 'pending_tables')
            self.assertFalse((Path(tmp)/'CosmiXsPPPC_positron_flux_earth.dat').exists())
            for species in cr_earth.ACTIVE:
                self.assertTrue((Path(tmp)/summary[species]['file']).is_file())


if __name__ == '__main__':
    unittest.main()
