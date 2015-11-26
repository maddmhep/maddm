# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 15:06:46


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.d__tilde__, P.d, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})

V_2 = Vertex(name = 'V_2',
             particles = [ P.s__tilde__, P.s, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})

V_3 = Vertex(name = 'V_3',
             particles = [ P.b__tilde__, P.b, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})

V_4 = Vertex(name = 'V_4',
             particles = [ P.d__tilde__, P.d, P.fdm__tilde__, P.P__tilde__fdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF4, L.FFFF5 ],
             couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,4):C.GC_1,(0,3):C.GC_2})

V_5 = Vertex(name = 'V_5',
             particles = [ P.s__tilde__, P.s, P.fdm__tilde__, P.P__tilde__fdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF5, L.FFFF6 ],
             couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,3):C.GC_1,(0,4):C.GC_2})

V_6 = Vertex(name = 'V_6',
             particles = [ P.b__tilde__, P.b, P.fdm__tilde__, P.P__tilde__fdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF4, L.FFFF5 ],
             couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,4):C.GC_1,(0,3):C.GC_2})

V_7 = Vertex(name = 'V_7',
             particles = [ P.u__tilde__, P.u, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})

V_8 = Vertex(name = 'V_8',
             particles = [ P.c__tilde__, P.c, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})

V_9 = Vertex(name = 'V_9',
             particles = [ P.t__tilde__, P.t, P.sdm__tilde__, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_4})

V_10 = Vertex(name = 'V_10',
              particles = [ P.fdm__tilde__, P.P__tilde__fdm, P.u__tilde__, P.u ],
              color = [ 'Identity(3,4)' ],
              lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF5, L.FFFF6 ],
              couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,3):C.GC_1,(0,4):C.GC_2})

V_11 = Vertex(name = 'V_11',
              particles = [ P.fdm__tilde__, P.P__tilde__fdm, P.c__tilde__, P.c ],
              color = [ 'Identity(3,4)' ],
              lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF4, L.FFFF5 ],
              couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,4):C.GC_1,(0,3):C.GC_2})

V_12 = Vertex(name = 'V_12',
              particles = [ P.fdm__tilde__, P.P__tilde__fdm, P.t__tilde__, P.t ],
              color = [ 'Identity(3,4)' ],
              lorentz = [ L.FFFF1, L.FFFF2, L.FFFF3, L.FFFF5, L.FFFF6 ],
              couplings = {(0,0):C.GC_3,(0,1):C.GC_2,(0,2):C.GC_1,(0,3):C.GC_1,(0,4):C.GC_2})

V_13 = Vertex(name = 'V_13',
              particles = [ P.d__tilde__, P.d, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})

V_14 = Vertex(name = 'V_14',
              particles = [ P.s__tilde__, P.s, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})

V_15 = Vertex(name = 'V_15',
              particles = [ P.b__tilde__, P.b, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})

V_16 = Vertex(name = 'V_16',
              particles = [ P.u__tilde__, P.u, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})

V_17 = Vertex(name = 'V_17',
              particles = [ P.c__tilde__, P.c, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})

V_18 = Vertex(name = 'V_18',
              particles = [ P.t__tilde__, P.t, P.vdm__tilde__, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_6,(0,1):C.GC_5})

