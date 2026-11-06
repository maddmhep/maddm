"""Numerical and physics validation. Run: python -m unittest discover -s tests -v."""
from dataclasses import replace
import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.integrate import quad
from scipy.optimize import differential_evolution

from pbar.core import (DMTemplate, FitConfig, FluxTable, Likelihood, SIGMAV_REF,
                        Spectrum, dtdr, force_field, kinetic_energy)
from pbar.core import MP


class PhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = FluxTable()
        cls.likelihood = Likelihood()

    def test_force_field_identity_and_charge(self):
        f = lambda t: np.asarray(t)**-2
        t = np.array([0.5, 5., 50.])
        np.testing.assert_allclose(force_field(f, t, 0.), f(t), rtol=1.e-14)
        phi = .7
        expected = (t+phi)**-2*t*(t+2*MP)/((t+phi)*(t+phi+2*MP))
        np.testing.assert_allclose(force_field(f, t, phi), expected, rtol=1.e-14)
        self.assertFalse(np.allclose(expected, force_field(f, t, phi/2), rtol=1.e-3))

    def test_rigidity_jacobian_conserves_flux(self):
        a, b = 1.08, 427.5
        f = lambda t: (1+t)**-2.4
        rigidity_integral = quad(lambda r: f(kinetic_energy(r))*dtdr(r), a, b, epsabs=1.e-12)[0]
        energy_integral = quad(f, kinetic_energy(a), kinetic_energy(b), epsabs=1.e-12)[0]
        self.assertAlmostEqual(rigidity_integral/energy_integral, 1., places=11)

    def test_interpolation_nodes_and_zero_handling(self):
        a = np.loadtxt(self.table.path)
        node = a[a[:, 0] == 100.]
        f = self.table.template(100., {'bb': 1.})
        np.testing.assert_allclose(f(10.**node[:, 1]), node[:, 12], rtol=5.e-13)
        zero = Spectrum([1., 2., 4., 8.], [0., 0., 1., 2.])
        self.assertEqual(float(zero(1.5)), 0.)
        self.assertGreater(float(zero(3.)), 0.)
        np.testing.assert_allclose(zero([1., 2., 4., 8.]), [0., 0., 1., 2.], atol=1.e-14)
        with self.assertRaises(ValueError):
            zero(.9)
        with self.assertRaises(ValueError):
            self.table.template(1.e6)

    def test_interpolation_derivative_and_interval_lower_bound(self):
        f = Spectrum([1., 2., 4., 8., 16.], [0., 0., 1., .3, .1])
        x = np.array([1.2, 2.3, 3.3, 5.1, 10.7])
        step = 1.e-5
        numerical = (f(x+step)-f(x-step))/(2*step)
        np.testing.assert_allclose(f.derivative(x), numerical, rtol=2.e-6, atol=1.e-10)
        for low, high in [(1.1, 3.9), (3., 10.), (4.5, 15.)]:
            envelope = f.minimum(low, high)
            self.assertLessEqual(float(envelope), float(np.min(f(np.geomspace(low, high, 1000))))+1.e-13)

    def test_linear_mass_and_branching_and_density(self):
        low = self.table.template(100., {'bb': 1.})
        high = self.table.template(110., {'bb': 1.})
        mid = self.table.template(105., {'bb': 1.})
        e = self.table.energy
        np.testing.assert_allclose(mid(e), .5*(low(e)+high(e)), rtol=5.e-13)
        mixed = self.table.template(100., {'bb': .7, 'WW': .3})
        ww = self.table.template(100., {'WW': 1.})
        x = np.geomspace(.2, 500., 37)
        np.testing.assert_allclose(mixed(x), .7*low(x)+.3*ww(x), rtol=1.e-14)
        denser = self.table.template(100., {'bb': 1.}, rho_local_GeVcm3=.8, source_factor=.5)
        np.testing.assert_allclose(denser(x), 2*low(x), rtol=1.e-14)
        with self.assertRaises(ValueError):
            self.table.template(100., {'bb': .5})

    def test_custom_template_normalization(self):
        x = np.geomspace(.1, 1.e4, 100)
        a = DMTemplate.from_propagated_flux(x, 1/x, sigmav_ref_cm3s=2*SIGMAV_REF)
        b = Spectrum(x, .5/x)
        np.testing.assert_allclose(a(x), b(x), rtol=1.e-14)

    def test_no_spurious_limit_from_empty_channel(self):
        empty = self.table.template(10., {'nu_e_nu_e': 1.})
        self.assertTrue(np.all(empty([1., 10., 100.]) == 0.))
        result = self.likelihood.upper_limit(empty)
        self.assertEqual(result['status'], 'no_signal_in_table')
        self.assertIsNone(result['sigmav95_cm3s'])

    def test_profile_against_independent_full_optimizer(self):
        f = self.table.template(100., {'bb': 1.})
        for rate in [0., 3.e-26, 1.e-25]:
            ours = self.likelihood.profile(f, rate)
            def full(x):
                sec, dm = self.likelihood.components(f, rate, *x)
                return np.sum(((self.likelihood.observed-sec-dm)/self.likelihood.errors)**2)
            c = self.likelihood.config
            independent = differential_evolution(full, [c.norm_bounds, c.slope_bounds, c.phi_bounds_GV],
                                seed=12, tol=1.e-10, atol=1.e-10, polish=True)
            self.assertTrue(independent.success)
            self.assertLess(abs(ours.chi2-independent.fun), 2.e-5)

    def test_secondary_tilt_before_modulation(self):
        l = self.likelihood
        phi, slope, norm = .6, .1, 1.1
        sec, dm = l.components(None, 0., norm, slope, phi)
        deformed = lambda t: norm*l.secondary(t)*(t/10.)**slope
        expected = force_field(deformed, l.t, phi)*l.jacobian
        np.testing.assert_allclose(sec, expected, rtol=1.e-14)
        self.assertTrue(np.all(dm == 0))

    def test_solar_modulation_recomputed_for_dm(self):
        f = self.table.template(100., {'bb': 1.})
        l = self.likelihood
        for phi in [.2, 1.2]:
            _, dm = l.components(f, SIGMAV_REF, 1., 0., phi)
            expected = force_field(f, l.t, phi)*l.jacobian
            np.testing.assert_allclose(dm, expected, rtol=1.e-14)

    def test_asimov_nuisance_recovery(self):
        l = self.likelihood
        sec, _ = l.components(None, 0., 1.07, -.034, .73)
        synthetic = Likelihood(observed_flux=sec)
        best = synthetic.profile(None, 0.)
        self.assertLess(best.chi2, 1.e-9)
        np.testing.assert_allclose([best.norm, best.slope, best.phi_GV], [1.07, -.034, .73], atol=2.e-5)


class LimitTests(unittest.TestCase):
    def test_asimov_upper_against_known_gaussian_answer(self):
        # Fix nuisance freedom to a negligible interval. An Asimov Gaussian
        # upper limit then has a closed-form answer sqrt(q)/||signal/error||.
        c = FitConfig(norm_bounds=(1.-1.e-9, 1.+1.e-9), slope_bounds=(-1.e-10, 1.e-10),
                      phi_bounds_GV=(.7-1.e-9, .7+1.e-9))
        l0 = Likelihood(c)
        template = FluxTable().template(100., {'bb': 1.})
        secondary, signal = l0.components(template, SIGMAV_REF, 1., 0., .7)
        l = Likelihood(c, observed_flux=secondary)
        expected = SIGMAV_REF*np.sqrt(c.delta_chi2)/np.linalg.norm(signal/l.errors)
        r = l.upper_limit(template)
        self.assertLess(abs(r['sigmav95_cm3s']/expected-1), 3.e-5)
        self.assertLess(r['best_fit']['sigmav_cm3s'], 1.e-30)
        self.assertAlmostEqual(r['delta_chi2_at_limit'], c.delta_chi2, places=4)
        self.assertTrue(all(.8 <= p['norm'] <= 1.2 for p in r['profile']))


if __name__ == '__main__':
    unittest.main()
