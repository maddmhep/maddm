# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 13.1.0 for Microsoft Windows (64-bit) (August 14, 2022)
# Date: Fri 21 Oct 2022 11:22:11


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot
try:
   import form_factors as ForFac 
except ImportError:
   pass


FFFF1 = Lorentz(name = 'FFFF1',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,1)*Gamma(-1,4,3)')

FFFF2 = Lorentz(name = 'FFFF2',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma5(-3,3)*Gamma5(-2,1)*Gamma(-1,2,-2)*Gamma(-1,4,-3)')

FFFF3 = Lorentz(name = 'FFFF3',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-2,-3,1)*Gamma(-2,4,-4)*Gamma(-1,-4,3)*Gamma(-1,2,-3)')

FFFF4 = Lorentz(name = 'FFFF4',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-2,-3,3)*Gamma(-2,2,-4)*Gamma(-1,-4,1)*Gamma(-1,4,-3)')

FFFF5 = Lorentz(name = 'FFFF5',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-2,-4,3)*Gamma(-2,-3,1)*Gamma(-1,2,-3)*Gamma(-1,4,-4)')

FFFF6 = Lorentz(name = 'FFFF6',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Identity(2,1)*Identity(4,3)')

