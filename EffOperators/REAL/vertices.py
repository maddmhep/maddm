# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 14:25:58


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.d__tilde__, P.d, P.Fsdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_2 = Vertex(name = 'V_2',
             particles = [ P.s__tilde__, P.s, P.Fsdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_3 = Vertex(name = 'V_3',
             particles = [ P.b__tilde__, P.b, P.Fsdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_4 = Vertex(name = 'V_4',
             particles = [ P.d__tilde__, P.d, P.P__tilde__sdm, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_5})

V_5 = Vertex(name = 'V_5',
             particles = [ P.s__tilde__, P.s, P.P__tilde__sdm, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_5})

V_6 = Vertex(name = 'V_6',
             particles = [ P.b__tilde__, P.b, P.P__tilde__sdm, P.P__tilde__sdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFSS1 ],
             couplings = {(0,0):C.GC_5})

V_7 = Vertex(name = 'V_7',
             particles = [ P.P__tilde__fdm, P.P__tilde__fdm, P.Fsdm ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_4})

V_8 = Vertex(name = 'V_8',
             particles = [ P.u__tilde__, P.u, P.Fsdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_9 = Vertex(name = 'V_9',
             particles = [ P.c__tilde__, P.c, P.Fsdm ],
             color = [ 'Identity(1,2)' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_10 = Vertex(name = 'V_10',
              particles = [ P.t__tilde__, P.t, P.Fsdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1 ],
              couplings = {(0,0):C.GC_3})

V_11 = Vertex(name = 'V_11',
              particles = [ P.u__tilde__, P.u, P.P__tilde__sdm, P.P__tilde__sdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFSS1 ],
              couplings = {(0,0):C.GC_5})

V_12 = Vertex(name = 'V_12',
              particles = [ P.c__tilde__, P.c, P.P__tilde__sdm, P.P__tilde__sdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFSS1 ],
              couplings = {(0,0):C.GC_5})

V_13 = Vertex(name = 'V_13',
              particles = [ P.t__tilde__, P.t, P.P__tilde__sdm, P.P__tilde__sdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFSS1 ],
              couplings = {(0,0):C.GC_5})

V_14 = Vertex(name = 'V_14',
              particles = [ P.d__tilde__, P.d, P.P__tilde__vdm, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_7,(0,1):C.GC_6})

V_15 = Vertex(name = 'V_15',
              particles = [ P.s__tilde__, P.s, P.P__tilde__vdm, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_7,(0,1):C.GC_6})

V_16 = Vertex(name = 'V_16',
              particles = [ P.b__tilde__, P.b, P.P__tilde__vdm, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_7,(0,1):C.GC_6})

V_17 = Vertex(name = 'V_17',
              particles = [ P.u__tilde__, P.u, P.P__tilde__vdm, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_7,(0,1):C.GC_6})

V_18 = Vertex(name = 'V_18',
              particles = [ P.c__tilde__, P.c, P.P__tilde__vdm, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_7,(0,1):C.GC_6})

V_19 = Vertex(name = 'V_19',
              particles = [ P.t__tilde__, P.t, P.P__tilde__vdm, P.P__tilde__vdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFVV1, L.FFVV2 ],
              couplings = {(0,0):C.GC_7,(0,1):C.GC_6})

V_20 = Vertex(name = 'V_20',
              particles = [ P.d__tilde__, P.d, P.Fvdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV2 ],
              couplings = {(0,0):C.GC_2})

V_21 = Vertex(name = 'V_21',
              particles = [ P.s__tilde__, P.s, P.Fvdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV2 ],
              couplings = {(0,0):C.GC_2})

V_22 = Vertex(name = 'V_22',
              particles = [ P.b__tilde__, P.b, P.Fvdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV2 ],
              couplings = {(0,0):C.GC_2})

V_23 = Vertex(name = 'V_23',
              particles = [ P.P__tilde__fdm, P.P__tilde__fdm, P.Fvdm ],
              color = [ '1' ],
              lorentz = [ L.FFV1 ],
              couplings = {(0,0):C.GC_1})

V_24 = Vertex(name = 'V_24',
              particles = [ P.u__tilde__, P.u, P.Fvdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV2 ],
              couplings = {(0,0):C.GC_2})

V_25 = Vertex(name = 'V_25',
              particles = [ P.c__tilde__, P.c, P.Fvdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV2 ],
              couplings = {(0,0):C.GC_2})

V_26 = Vertex(name = 'V_26',
              particles = [ P.t__tilde__, P.t, P.Fvdm ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV2 ],
              couplings = {(0,0):C.GC_2})

