# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 13.1.0 for Microsoft Windows (64-bit) (August 14, 2022)
# Date: Fri 4 Nov 2022 18:20:51


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.e__plus__, P.e__minus__, P.Xd__tilde__, P.Xd ],
             color = [ '1' ],
             lorentz = [ L.FFFF1, L.FFFF2 ],
             couplings = {(0,0):C.GC_1,(0,1):C.GC_2})

