"""Shared Galactic distributions: recoil integrals, pair moments and limits."""
from pathlib import Path
import ctypes
import json
import math
import os
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad, trapezoid

from indirect import halo_model as halo


class VelocityProfilesTests(unittest.TestCase):
    def test_default_preserves_native_shm_and_numeric_agrees(self):
        card = {'vmp': 202., 'vescape': 480.}
        default = halo.velocity_profile(card)
        numeric = halo.velocity_profile(dict(card, halo_velocity_profile='shm_numeric'))
        self.assertEqual(default.name, 'SHM')
        self.assertEqual(halo.galactic_rule(default), halo.relative_speed_rule(202., 480.))
        speeds = np.linspace(0., 760., 39)
        np.testing.assert_allclose(halo.mean_inverse_speed(speeds, 202., 480., 250.),
                                   numeric.primary.eta(speeds, 250.), rtol=1e-8, atol=1e-13)
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

    def test_lmc_direction_does_not_change_annihilation(self):
        card = {'halo_velocity_profile': 'SHM_wLMC', 'halo_lmc_fraction': .2}
        rule = halo.galactic_rule(halo.velocity_profile(card))
        self.assertEqual(rule, halo.galactic_rule(halo.velocity_profile(dict(card, halo_lmc_cosb=.3))))

    def test_reject_invalid_profiles(self):
        for opts in ({'halo_velocity_profile':'unknown'}, {'halo_velocity_profile':'SHM_beta','halo_beta':1.2},
                     {'halo_velocity_profile':'Lisanti','halo_lisanti_k':.05},
                     {'halo_velocity_profile':'SHM_wLMC','halo_lmc_fraction':1.1},
                     {'halo_velocity_profile':'SHM_wLMC','halo_lmc_cosb':1.5}):
            with self.subTest(card=opts), self.assertRaises(ValueError):
                halo.velocity_profile(opts)
        # with fixed-speed annihilation the profile is still read (it enters the recoil rates)
        fixed = halo.velocity_profile({'halo_velocity_profile': 'Lisanti', 'mw_velocity_model': 'fixed'})
        self.assertEqual(fixed.name, 'Lisanti')

    def test_saved_tables_and_metadata_use_selected_profile(self):
        card = {'halo_velocity_profile': 'Lisanti', 'halo_lisanti_k': 1.7, 'rhoDM': .47}
        profile = halo.velocity_profile(card)
        rule = halo.galactic_rule(profile)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            halo.write_velocity_diagnostics(root, card, profile, rule)
            data = json.loads((root/halo.VELOCITY_OUTPUT).read_text())
            self.assertEqual((data['halo_velocity_profile'], data['halo_lisanti_k']), ('Lisanti', 1.7))
            self.assertEqual(data['rhoDM_gev_cm3'], .47)
            np.testing.assert_allclose(np.loadtxt(root/data['mw_quadrature_file']), np.array(rule)*[halo.C_KMS, 1.])

    @unittest.skipUnless(os.environ.get('MADDM_RAPIDD_TEST_LIB'), 'requires the optional RAPIDD shared library')
    def test_fixed_lab_eta_matches_rapidd_c(self):
        # Standalone comparison at one lab speed.
        lib = ctypes.CDLL(os.environ['MADDM_RAPIDD_TEST_LIB'])
        f = lib.halo_f
        f.argtypes = [ctypes.c_char_p]+[ctypes.c_double]*6+[ctypes.c_int]+[ctypes.c_double]*5
        f.restype = ctypes.c_double
        for profile in ('SHM', 'SHM_numeric', 'SHM_beta', 'Lisanti', 'SHM_wLMC'):
            for vmin in (80., 350., 680., 850.):
                with self.subTest(profile=profile, threshold=vmin):
                    expected = f(profile.encode(), vmin, 544., 238., 0., 250., 1.5, 0,
                                 .006, 570., -.71, 100., 200.)
                    if profile in ('Lisanti', 'SHM_beta'):
                        # SHM_beta with the beta = 0 of the C call: the hard-cutoff SHM
                        p = halo.velocity_profile({'halo_velocity_profile': profile, 'halo_beta': 0.})
                        actual = p.primary.eta(vmin, 250.)
                    else:
                        actual = halo.mean_inverse_speed(vmin, 238., 544., 250.)
                    if profile == 'SHM_wLMC':
                        drift = math.sqrt(570.**2+250.**2-2*570.*250.*(-.71))
                        actual = .994*actual+.006*halo.mean_inverse_speed(vmin, 100., 200., drift)
                    self.assertAlmostEqual(float(actual), expected, delta=max(1.e-10, 2.e-5*expected))


RAPIDD_DIR = Path(__file__).resolve().parents[1]/'vendor'/'RAPIDD_for_DM'
RAPIDD_LIB = RAPIDD_DIR/'lib'/'build'/('libRAPIDD.dylib' if sys.platform == 'darwin' else 'libRAPIDD.so')

# card entries of the MadDM defaults plus a non-default value of every profile parameter
SHARED_CARD = {'vmp': 225., 'vescape': 520., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3,
               'vEarth_mod': 29.8, 'halo_beta': .6, 'halo_lisanti_k': 2.2, 'halo_lmc_fraction': .15,
               'halo_lmc_vb': 540., 'halo_lmc_sigb': 90., 'halo_lmc_vcut': 180., 'halo_lmc_cosb': -.5}


@unittest.skipUnless(RAPIDD_LIB.is_file(), 'requires the RAPIDD shared library')
class SharedDirectIndirectHaloTests(unittest.TestCase):
    """Direct detection (RAPIDD tables) and MW annihilation read the same card selection."""

    @classmethod
    def setUpClass(cls):
        if str(RAPIDD_DIR) not in sys.path:
            sys.path.insert(0, str(RAPIDD_DIR))
        import rapidd.core
        cls.core = rapidd.core

    def test_rapidd_reads_the_profile_and_parameters_of_annihilation(self):
        for name in halo.PROFILE_NAMES:
            with self.subTest(profile=name):
                card = dict(SHARED_CARD, halo_velocity_profile=name.lower())
                mw = halo.velocity_profile(card)
                dd = self.core.halo_kwargs_from_card(card)
                self.assertEqual(dd['profile'], mw.name)
                self.assertEqual((dd['v0'], dd['vesc']), (mw.primary.scale, mw.primary.cutoff))
                if name == 'SHM_beta':
                    self.assertEqual(dd['beta'], mw.primary.parameter)
                elif name == 'Lisanti':
                    self.assertEqual(dd['k_lisanti'], mw.primary.parameter)
                elif name == 'SHM_wLMC':
                    self.assertEqual((dd['w'], dd['vb'], dd['sigb'], dd['vcut'], dd['cosb']),
                                     (mw.fraction, mw.bulk_speed, mw.stream.scale, mw.stream.cutoff, mw.cosb))

    def test_rapidd_table_is_the_lab_eta_of_the_annihilation_law(self):
        """The RAPIDD annual-average table (eta_0 and eta_2) of every profile is the lab-frame boost of
        the Galactic law averaged by the annihilation rates."""
        ve = self.core.lab_frame_speed_annual_avg(SHARED_CARD['vmp'], (11.1, 12.2, 7.3))
        with tempfile.TemporaryDirectory() as tmp:
            for name in halo.PROFILE_NAMES:
                with self.subTest(profile=name):
                    card = dict(SHARED_CARD, halo_velocity_profile=name)
                    path = self.core.halo_path_from_card(card, tmp)
                    self.assertEqual(Path(path).name, name+'.dat')
                    table = np.loadtxt(path, skiprows=1)
                    v, eta0, eta2 = table[:, 0], table[:, 1], table[:, 3]
                    p = halo.velocity_profile(card)
                    laws = [(1-p.fraction, p.primary, ve)]
                    if p.stream is not None:
                        drift = math.sqrt(p.bulk_speed**2+ve**2-2*p.bulk_speed*ve*p.cosb)
                        laws.append((p.fraction, p.stream, drift))
                    want0 = sum(f*law.eta(v, b) for f, law, b in laws)
                    np.testing.assert_allclose(eta0, want0, rtol=2e-5, atol=1e-5*want0.max())
                    # eta_2 = int_{v>vmin} v f d^3v = vmin^2 eta_0(vmin) + 2 int_{vmin} x eta_0(x) dx
                    for vmin, value in zip(v[::20], eta2[::20]):
                        want2 = sum(f*(vmin**2*float(law.eta(vmin, b))
                                       + 2*quad(lambda x: x*law.eta(x, b), vmin, law.cutoff+b, limit=200)[0])
                                    for f, law, b in laws)
                        self.assertAlmostEqual(value, want2, delta=2e-4*eta2[0])

    def test_equivalent_reference_uses_the_standard_shm(self):
        with tempfile.TemporaryDirectory() as tmp:
            card = dict(SHARED_CARD, halo_velocity_profile='SHM_beta', halo_beta=1.)
            self.assertNotEqual(self.core.halo_path_from_card(card, tmp), self.core.halo_path)
            default = dict(SHARED_CARD, vmp=238., vescape=544.)
            self.assertEqual(self.core.halo_path_from_card(default, tmp), self.core.halo_path)


if __name__ == '__main__':
    unittest.main()
