"""Shared Galactic distributions: recoil integrals, pair moments and limits."""
from pathlib import Path
import ctypes
import json
import math
import os
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad, trapezoid

import halo_model as halo


class VelocityProfilesTests(unittest.TestCase):
    def test_default_preserves_native_shm_and_numeric_agrees(self):
        card = {'vmp': 202., 'vescape': 480.}
        default = halo.velocity_profile(card)
        numeric = halo.velocity_profile(dict(card, halo_velocity_profile='shm_numeric'))
        self.assertEqual(default.name, 'SHM')
        self.assertEqual(halo.galactic_rule(default), halo.relative_speed_rule(202., 480.))
        np.testing.assert_allclose(halo.recoil_halo_table(card, 200, default),
                                   halo.recoil_halo_table(card, 200, numeric), rtol=1e-8, atol=1e-13)
        self.assertAlmostEqual(halo.rms_speed(halo.galactic_rule(default)),
                               halo.rms_speed(halo.galactic_rule(numeric)), places=12)

    def test_alternative_isotropic_profiles_match_radial_moments(self):
        v0, cut = 185., 470.
        laws = [({'halo_velocity_profile': 'SHM_beta', 'halo_beta': .8},
                 lambda r: math.exp(-(r/v0)**2)-.8*math.exp(-(cut/v0)**2)),
                ({'halo_velocity_profile': 'Lisanti', 'halo_lisanti_k': 1.5},
                 lambda r: (math.exp(-(r/v0)**2/1.5)-math.exp(-(cut/v0)**2/1.5))**1.5)]
        for settings, f in laws:
            with self.subTest(profile=settings['halo_velocity_profile']):
                profile = halo.velocity_profile(dict(settings, vmp=v0, vescape=cut))
                q, w = np.array(halo.galactic_rule(profile)).T
                u = q*halo.C_KMS
                r2, r4 = self.moments(f, cut)
                self.assertAlmostEqual(w.sum(), 1., places=14)
                self.assertTrue(np.all((u > 0) & (u < 2*cut)))
                self.assertAlmostEqual(np.dot(w, u**2)/(2*r2), 1., places=8)
                self.assertAlmostEqual(np.dot(w, u**4)/(2*r4+10*r2*r2/3), 1., places=8)
                # At zero boost eta is a direct integral of the one-particle PDF.
                norm = quad(lambda r: r*r*f(r), 0, cut)[0]
                for threshold in (0., 220., 450., 469.9, 490.):
                    want = quad(lambda r: r*f(r), threshold, cut)[0]/norm if threshold < cut else 0.
                    self.assertAlmostEqual(float(profile.primary.eta(threshold, 0.)), want, delta=2.e-11)

    def test_lisanti_k_one_is_beta_one(self):
        beta = halo.velocity_profile({'halo_velocity_profile': 'SHM_beta'})
        lisanti = halo.velocity_profile({'halo_velocity_profile': 'Lisanti', 'halo_lisanti_k': 1.})
        speeds = np.linspace(0., 810., 65)
        np.testing.assert_allclose(beta.primary.eta(speeds, 250.), lisanti.primary.eta(speeds, 250.),
                                   atol=1.e-13, rtol=1.e-10)
        self.assertAlmostEqual(halo.rms_speed(halo.galactic_rule(beta)),
                               halo.rms_speed(halo.galactic_rule(lisanti)), places=12)

    @staticmethod
    def moments(f, cut):
        norm = quad(lambda r: r**2*f(r), 0., cut)[0]
        return [quad(lambda r: r**k*f(r), 0., cut)[0]/norm for k in (4, 6)]

    def test_lmc_pair_weights_and_cross_term(self):
        w, bulk = .18, 640.
        card = dict(halo_velocity_profile='SHM_wLMC', halo_lmc_fraction=w, halo_lmc_vb=bulk,
                    vmp=210., vescape=490., halo_lmc_sigb=110., halo_lmc_vcut=155.)
        q, weight = np.array(halo.galactic_rule(halo.velocity_profile(card))).T
        u = q*halo.C_KMS
        am = self.moments(lambda r: math.exp(-(r/210.)**2), 490.)
        bm = self.moments(lambda r: math.exp(-(r/110.)**2), 155.)
        def pair(a, b, boost):
            s2 = a[0]+b[0]
            return np.array([s2+boost**2, a[1]+b[1]+10/3*(a[0]*b[0]+boost**2*s2)+boost**4])
        expected = (1-w)**2*pair(am, am, 0.)+2*w*(1-w)*pair(am, bm, bulk)+w*w*pair(bm, bm, 0.)
        np.testing.assert_allclose(np.array([u**2, u**4])@weight, expected, rtol=2.e-9)
        self.assertGreater(max(u), 2*490.)
        self.assertLess(max(u), bulk+490.+155.)
        self.assertEqual(halo.galactic_rule(halo.velocity_profile(dict(card, halo_lmc_fraction=0.))),
                         halo.relative_speed_rule(210., 490.))
        self.assertEqual(halo.galactic_rule(halo.velocity_profile(dict(card, halo_lmc_fraction=1.))),
                         halo.relative_speed_rule(110., 155.))

    def test_lab_direction_changes_recoil_but_not_annihilation(self):
        card = {'halo_velocity_profile':'SHM_wLMC', 'halo_lmc_fraction':.2}
        p = halo.velocity_profile(card)
        vb = halo.lmc_bulk_velocity(p, card)
        ref = np.array([11.1+29.2, 238.+12.2-.1, 7.3+5.9])
        self.assertAlmostEqual(np.dot(vb,ref)/(np.linalg.norm(vb)*np.linalg.norm(ref)), p.cosb)
        other = dict(card, halo_lmc_azimuth_deg=90.)
        self.assertEqual(halo.galactic_rule(p), halo.galactic_rule(halo.velocity_profile(other)))
        a, b = halo.recoil_halo_table(card, 200), halo.recoil_halo_table(other, 200)
        self.assertFalse(np.array_equal(a, b))
        for x, eta in (a.T, b.T):
            self.assertAlmostEqual(trapezoid(eta, x), 1., delta=2.e-5)
            self.assertTrue(np.all(np.diff(eta) <= 1.e-13))

    def test_reject_invalid_profiles_and_fixed_speed_combination(self):
        for opts in ({'halo_velocity_profile':'unknown'}, {'halo_velocity_profile':'SHM_beta','halo_beta':1.2},
                     {'halo_velocity_profile':'Lisanti','halo_lisanti_k':.05},
                     {'halo_velocity_profile':'SHM_wLMC','halo_lmc_fraction':1.1},
                     {'halo_velocity_profile':'SHM_wLMC','halo_lmc_cosb':1.5},
                     {'halo_velocity_profile':'Lisanti','mw_velocity_model':'fixed'}):
            with self.subTest(card=opts), self.assertRaises(ValueError):
                halo.velocity_profile(opts)

    def test_saved_tables_and_metadata_use_selected_profile(self):
        card = {'halo_velocity_profile': 'Lisanti', 'halo_lisanti_k': 1.7, 'rholocal': .47}
        profile = halo.velocity_profile(card)
        rule = halo.galactic_rule(profile)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            default = root/'vendor.dat'
            default.write_text('200 2\n')
            halo.write_recoil_halo(root/'DD_halo.dat', card, default, profile)
            halo.write_velocity_diagnostics(root, card, profile, rule, 'DD_halo.dat')
            data = json.loads((root/halo.VELOCITY_OUTPUT).read_text())
            self.assertEqual((data['halo_velocity_profile'], data['halo_lisanti_k']), ('Lisanti', 1.7))
            self.assertEqual(data['rapidd_recoil_halo_file'], 'DD_halo.dat')
            np.testing.assert_allclose(np.loadtxt(root/data['mw_quadrature_file']), np.array(rule)*[halo.C_KMS, 1.])
            np.testing.assert_allclose(np.loadtxt(root/'DD_halo.dat', skiprows=1), halo.recoil_halo_table(card, 200))

    @unittest.skipUnless(os.environ.get('MADDM_RAPIDD_TEST_LIB'), 'requires the optional RAPIDD shared library')
    def test_fixed_lab_eta_matches_rapidd_c(self):
        # Standalone comparison is at one lab speed, so annual averaging and
        # the azimuth convention do not enter the comparison.
        lib = ctypes.CDLL(os.environ['MADDM_RAPIDD_TEST_LIB'])
        f = lib.halo_f
        f.argtypes = [ctypes.c_char_p]+[ctypes.c_double]*6+[ctypes.c_int]+[ctypes.c_double]*5
        f.restype = ctypes.c_double
        for profile in ('SHM', 'SHM_numeric', 'Lisanti', 'SHM_wLMC'):
            for vmin in (80., 350., 680., 850.):
                with self.subTest(profile=profile, threshold=vmin):
                    expected = f(profile.encode(), vmin, 544., 238., 0., 250., 1.5, 0,
                                 .006, 570., -.71, 100., 200.)
                    if profile == 'Lisanti':
                        p = halo.velocity_profile({'halo_velocity_profile': profile})
                        actual = p.primary.eta(vmin, 250.)
                    else:
                        actual = halo.mean_inverse_speed(vmin, 238., 544., 250.)
                    if profile == 'SHM_wLMC':
                        drift = math.sqrt(570.**2+250.**2-2*570.*250.*(-.71))
                        actual = .994*actual+.006*halo.mean_inverse_speed(vmin, 100., 200., drift)
                    self.assertAlmostEqual(float(actual), expected, delta=max(1.e-10, 2.e-5*expected))


if __name__ == '__main__':
    unittest.main()
