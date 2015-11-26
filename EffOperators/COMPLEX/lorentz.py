# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 15:06:46


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot


FFSS1 = Lorentz(name = 'FFSS1',
                spins = [ 2, 2, 1, 1 ],
                structure = 'Identity(2,1)')

FFFF1 = Lorentz(name = 'FFFF1',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Identity(2,1)*Identity(4,3)')

FFFF2 = Lorentz(name = 'FFFF2',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-3,3)*ProjM(-2,1)')

FFFF3 = Lorentz(name = 'FFFF3',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-3)*Gamma(-1,4,-2)*ProjM(-2,3)*ProjP(-3,1)')

FFFF4 = Lorentz(name = 'FFFF4',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-3)*Gamma(-1,4,-2)*ProjP(-3,1)*ProjP(-2,3)')

FFFF5 = Lorentz(name = 'FFFF5',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-2,1)*ProjP(-3,3)')

FFFF6 = Lorentz(name = 'FFFF6',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjP(-3,3)*ProjP(-2,1)')

FFVV1 = Lorentz(name = 'FFVV1',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Identity(2,1)*Metric(3,4)')

FFVV2 = Lorentz(name = 'FFVV2',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Epsilon(3,4,-1,-2)*P(-1,3)*Gamma(-2,2,-3)*ProjM(-3,1) - Epsilon(3,4,-1,-2)*P(-1,4)*Gamma(-2,2,-3)*ProjM(-3,1) - Epsilon(3,4,-1,-2)*P(-1,3)*Gamma(-2,2,-3)*ProjP(-3,1) + Epsilon(3,4,-1,-2)*P(-1,4)*Gamma(-2,2,-3)*ProjP(-3,1)')

