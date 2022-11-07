# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 15:06:46


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L



V_1 = Vertex(name = 'V_1',
             particles = [ P.e__plus__, P.e__minus__, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ '1' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})


V_2 = Vertex(name = 'V_2',
              particles = [ P.fdm__tilde__, P.P__tilde__fdm, P.e__plus__, P.e__minus__],
              color = [ '1' ],
              lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF5, L.FFFF6 ],
              couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,3):C.GC_1,(0,4):C.GC_2})


V_3 = Vertex(name = 'V_3',
              particles = [ P.e__plus__, P.e__minus__, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ '1' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})


