import math
import unittest

from direct_recoil import native_nucleon_cross_sections


class NativeNucleonTests(unittest.TestCase):
    def test_analytic_vector_and_axial_projections(self):
        card = dict(SPd=0., SNd=0., VPd=1., VNd=2.,
                    AVPd=-.42, AVNd=.84, SigPd=0., SigNd=0.)
        xsec = native_nucleon_cross_sections(card, {'alpha_d': (0., 2.e-7, 3.e-7, 0.)}, 10.)
        for target, mn, vector, axial in [('p', .93827208943, 1., -.42),
                                         ('n', .93956542194, 2., .84)]:
            common = 3.893793721e-28 * (10.*mn/(10.+mn))**2 / math.pi
            self.assertAlmostEqual(xsec['sigmaN_SI_'+target] / (common*(2.e-7*vector)**2), 1.)
            self.assertAlmostEqual(xsec['sigmaN_SD_'+target] / (3.*common*(3.e-7*axial)**2), 1.)

    def test_missing_or_nonfinite_projections_are_not_a_zero_prediction(self):
        xsec = native_nucleon_cross_sections({}, {'alpha_d': [-1.]*4,
                                                  'alpha_u': [math.nan, 0., 0., 0.]}, 10.)
        self.assertEqual(list(xsec.values()), [-1.]*4)
