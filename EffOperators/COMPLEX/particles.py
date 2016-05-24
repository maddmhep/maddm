# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 15:06:46


from __future__ import division
from object_library import all_particles, Particle
import parameters as Param

import propagators as Prop

u = Particle(pdg_code = 2,
             name = 'u',
             antiname = 'u~',
             spin = 2,
             color = 3,
             mass = Param.MU,
             width = Param.ZERO,
             texname = 'u',
             antitexname = 'u~',
             charge = 2/3,
             GhostNumber = 0,
             LeptonNumber = 0,
             Y = 0)

u__tilde__ = u.anti()

c = Particle(pdg_code = 4,
             name = 'c',
             antiname = 'c~',
             spin = 2,
             color = 3,
             mass = Param.MC,
             width = Param.ZERO,
             texname = 'c',
             antitexname = 'c~',
             charge = 2/3,
             GhostNumber = 0,
             LeptonNumber = 0,
             Y = 0)

c__tilde__ = c.anti()

t = Particle(pdg_code = 6,
             name = 't',
             antiname = 't~',
             spin = 2,
             color = 3,
             mass = Param.MT,
             width = Param.WT,
             texname = 't',
             antitexname = 't~',
             charge = 2/3,
             GhostNumber = 0,
             LeptonNumber = 0,
             Y = 0)

t__tilde__ = t.anti()

d = Particle(pdg_code = 1,
             name = 'd',
             antiname = 'd~',
             spin = 2,
             color = 3,
             mass = Param.MD,
             width = Param.ZERO,
             texname = 'd',
             antitexname = 'd~',
             charge = -1/3,
             GhostNumber = 0,
             LeptonNumber = 0,
             Y = 0)

d__tilde__ = d.anti()

s = Particle(pdg_code = 3,
             name = 's',
             antiname = 's~',
             spin = 2,
             color = 3,
             mass = Param.MS,
             width = Param.ZERO,
             texname = 's',
             antitexname = 's~',
             charge = -1/3,
             GhostNumber = 0,
             LeptonNumber = 0,
             Y = 0)

s__tilde__ = s.anti()

b = Particle(pdg_code = 5,
             name = 'b',
             antiname = 'b~',
             spin = 2,
             color = 3,
             mass = Param.MB,
             width = Param.ZERO,
             texname = 'b',
             antitexname = 'b~',
             charge = -1/3,
             GhostNumber = 0,
             LeptonNumber = 0,
             Y = 0)

b__tilde__ = b.anti()

P__tilde__sdm = Particle(pdg_code = 999000006,
                         name = '~sdm',
                         antiname = 'sdm~',
                         spin = 1,
                         color = 1,
                         mass = Param.sdmm,
                         width = Param.Wsdm,
                         texname = '~sdm',
                         antitexname = 'sdm~',
                         charge = 0,
                         GhostNumber = 0,
                         LeptonNumber = 0,
                         Y = 0)

sdm__tilde__ = P__tilde__sdm.anti()

P__tilde__fdm = Particle(pdg_code = 999000007,
                         name = '~fdm',
                         antiname = 'fdm~',
                         spin = 2,
                         color = 1,
                         mass = Param.fdmm,
                         width = Param.Wfdm,
                         texname = '~fdm',
                         antitexname = 'fdm~',
                         charge = 0,
                         GhostNumber = 0,
                         LeptonNumber = 0,
                         Y = 0)

fdm__tilde__ = P__tilde__fdm.anti()

P__tilde__vdm = Particle(pdg_code = 999000008,
                         name = '~vdm',
                         antiname = 'vdm~',
                         spin = 3,
                         color = 1,
                         mass = Param.vdmm,
                         width = Param.Wvdm,
                         texname = '~vdm',
                         antitexname = 'vdm~',
                         charge = 0,
                         GhostNumber = 0,
                         LeptonNumber = 0,
                         Y = 0)

vdm__tilde__ = P__tilde__vdm.anti()

