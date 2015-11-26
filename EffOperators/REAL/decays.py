# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 14:25:58


from object_library import all_decays, Decay
import particles as P


Decay_Fsdm = Decay(name = 'Decay_Fsdm',
                   particle = P.Fsdm,
                   partial_widths = {(P.d,P.d__tilde__):'((-24*fSIe**2*MD**2 + 6*fSIe**2*MFsdm**2)*cmath.sqrt(-4*MD**2*MFsdm**2 + MFsdm**4))/(16.*cmath.pi*abs(MFsdm)**3)',
                                     (P.s,P.s__tilde__):'((6*fSIe**2*MFsdm**2 - 24*fSIe**2*MS**2)*cmath.sqrt(MFsdm**4 - 4*MFsdm**2*MS**2))/(16.*cmath.pi*abs(MFsdm)**3)',
                                     (P.b,P.b__tilde__):'((-24*fSIe**2*MB**2 + 6*fSIe**2*MFsdm**2)*cmath.sqrt(-4*MB**2*MFsdm**2 + MFsdm**4))/(16.*cmath.pi*abs(MFsdm)**3)',
                                     (P.P__tilde__fdm,P.P__tilde__fdm):'((-32*fdmm**2*fSIe**2 + 8*fSIe**2*MFsdm**2)*cmath.sqrt(-4*fdmm**2*MFsdm**2 + MFsdm**4))/(32.*cmath.pi*abs(MFsdm)**3)',
                                     (P.u,P.u__tilde__):'((6*fSIe**2*MFsdm**2 - 24*fSIe**2*MU**2)*cmath.sqrt(MFsdm**4 - 4*MFsdm**2*MU**2))/(16.*cmath.pi*abs(MFsdm)**3)',
                                     (P.c,P.c__tilde__):'((-24*fSIe**2*MC**2 + 6*fSIe**2*MFsdm**2)*cmath.sqrt(-4*MC**2*MFsdm**2 + MFsdm**4))/(16.*cmath.pi*abs(MFsdm)**3)',
                                     (P.t,P.t__tilde__):'((6*fSIe**2*MFsdm**2 - 24*fSIe**2*MT**2)*cmath.sqrt(MFsdm**4 - 4*MFsdm**2*MT**2))/(16.*cmath.pi*abs(MFsdm)**3)'})

Decay_Fvdm = Decay(name = 'Decay_Fvdm',
                   particle = P.Fvdm,
                   partial_widths = {(P.d,P.d__tilde__):'((-24*fSDe**2*MD**2 + 24*fSDe**2*MFvdm**2)*cmath.sqrt(-4*MD**2*MFvdm**2 + MFvdm**4))/(48.*cmath.pi*abs(MFvdm)**3)',
                                     (P.s,P.s__tilde__):'((24*fSDe**2*MFvdm**2 - 24*fSDe**2*MS**2)*cmath.sqrt(MFvdm**4 - 4*MFvdm**2*MS**2))/(48.*cmath.pi*abs(MFvdm)**3)',
                                     (P.b,P.b__tilde__):'((-24*fSDe**2*MB**2 + 24*fSDe**2*MFvdm**2)*cmath.sqrt(-4*MB**2*MFvdm**2 + MFvdm**4))/(48.*cmath.pi*abs(MFvdm)**3)',
                                     (P.P__tilde__fdm,P.P__tilde__fdm):'((-64*fdmm**2*fSDe**2 + 16*fSDe**2*MFvdm**2)*cmath.sqrt(-4*fdmm**2*MFvdm**2 + MFvdm**4))/(96.*cmath.pi*abs(MFvdm)**3)',
                                     (P.u,P.u__tilde__):'((24*fSDe**2*MFvdm**2 - 24*fSDe**2*MU**2)*cmath.sqrt(MFvdm**4 - 4*MFvdm**2*MU**2))/(48.*cmath.pi*abs(MFvdm)**3)',
                                     (P.c,P.c__tilde__):'((-24*fSDe**2*MC**2 + 24*fSDe**2*MFvdm**2)*cmath.sqrt(-4*MC**2*MFvdm**2 + MFvdm**4))/(48.*cmath.pi*abs(MFvdm)**3)',
                                     (P.t,P.t__tilde__):'((24*fSDe**2*MFvdm**2 - 24*fSDe**2*MT**2)*cmath.sqrt(MFvdm**4 - 4*MFvdm**2*MT**2))/(48.*cmath.pi*abs(MFvdm)**3)'})

