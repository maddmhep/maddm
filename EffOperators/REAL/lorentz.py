# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 14:25:58


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot


FFS1 = Lorentz(name = 'FFS1',
               spins = [ 2, 2, 1 ],
               structure = 'Identity(2,1)')

FFV1 = Lorentz(name = 'FFV1',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma5(-1,1)*Gamma(3,2,-1)')

FFV2 = Lorentz(name = 'FFV2',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFSS1 = Lorentz(name = 'FFSS1',
                spins = [ 2, 2, 1, 1 ],
                structure = 'Identity(2,1)')

FFVV1 = Lorentz(name = 'FFVV1',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Identity(2,1)*Metric(3,4)')

FFVV2 = Lorentz(name = 'FFVV2',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Epsilon(3,4,-1,-2)*P(-1,3)*Gamma(-2,2,-3)*ProjM(-3,1) - Epsilon(3,4,-1,-2)*P(-1,4)*Gamma(-2,2,-3)*ProjM(-3,1) - Epsilon(3,4,-1,-2)*P(-1,3)*Gamma(-2,2,-3)*ProjP(-3,1) + Epsilon(3,4,-1,-2)*P(-1,4)*Gamma(-2,2,-3)*ProjP(-3,1)')

