# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 13.1.0 for Microsoft Windows (64-bit) (August 14, 2022)
# Date: Fri 21 Oct 2022 11:22:11


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.e__plus__, P.e__minus__, P.Xd__tilde__, P.Xd ],
             color = [ '1' ],
             lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF4, L.FFFF5, L.FFFF6 ],
             couplings = {(0,1):C.GC_1,(0,0):C.GC_5,(0,4):C.GC_4,(0,3):C.GC_3,(0,2):C.GC_3,(0,5):C.GC_2})

