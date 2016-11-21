# This file was automatically created by FeynRules 2.0.20
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 20 Jan 2015 14:25:58



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

sSIe = Parameter(name = 'sSIe',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\text{sSIe}',
                 lhablock = 'EFFBlock',
                 lhacode = [ 1 ])

fSIe = Parameter(name = 'fSIe',
                 nature = 'external',
                 type = 'real',
                 value = 10000000000,
                 texname = '\\text{fSIe}',
                 lhablock = 'EFFBlock',
                 lhacode = [ 2 ])

vSIe = Parameter(name = 'vSIe',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\text{vSIe}',
                 lhablock = 'EFFBlock',
                 lhacode = [ 3 ])

fSDe = Parameter(name = 'fSDe',
                 nature = 'external',
                 type = 'real',
                 value = 10000000000,
                 texname = '\\text{fSDe}',
                 lhablock = 'EFFBlock',
                 lhacode = [ 4 ])

vSDe = Parameter(name = 'vSDe',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\text{vSDe}',
                 lhablock = 'EFFBlock',
                 lhacode = [ 5 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.00255,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.27,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.00504,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.101,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

sdmm = Parameter(name = 'sdmm',
                 nature = 'external',
                 type = 'real',
                 value = 70,
                 texname = '\\text{sdmm}',
                 lhablock = 'MASS',
                 lhacode = [ 999000006 ])

MFsdm = Parameter(name = 'MFsdm',
                  nature = 'external',
                  type = 'real',
                  value = 10000000000,
                  texname = '\\text{MFsdm}',
                  lhablock = 'MASS',
                  lhacode = [ 999000007 ])

fdmm = Parameter(name = 'fdmm',
                 nature = 'external',
                 type = 'real',
                 value = 70,
                 texname = '\\text{fdmm}',
                 lhablock = 'MASS',
                 lhacode = [ 999000008 ])

vdmm = Parameter(name = 'vdmm',
                 nature = 'external',
                 type = 'real',
                 value = 70,
                 texname = '\\text{vdmm}',
                 lhablock = 'MASS',
                 lhacode = [ 999000009 ])

MFvdm = Parameter(name = 'MFvdm',
                  nature = 'external',
                  type = 'real',
                  value = 10000000000,
                  texname = '\\text{MFvdm}',
                  lhablock = 'MASS',
                  lhacode = [ 999000010 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

Wsdm = Parameter(name = 'Wsdm',
                 nature = 'external',
                 type = 'real',
                 value = 0.,
                 texname = '\\text{Wsdm}',
                 lhablock = 'DECAY',
                 lhacode = [ 999000006 ])

WFsdm = Parameter(name = 'WFsdm',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = '\\text{WFsdm}',
                  lhablock = 'DECAY',
                  lhacode = [ 999000007 ])

Wfdm = Parameter(name = 'Wfdm',
                 nature = 'external',
                 type = 'real',
                 value = 0.,
                 texname = '\\text{Wfdm}',
                 lhablock = 'DECAY',
                 lhacode = [ 999000008 ])

Wvdm = Parameter(name = 'Wvdm',
                 nature = 'external',
                 type = 'real',
                 value = 0.,
                 texname = '\\text{Wvdm}',
                 lhablock = 'DECAY',
                 lhacode = [ 999000009 ])

WFvdm = Parameter(name = 'WFvdm',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = '\\text{WFvdm}',
                  lhablock = 'DECAY',
                  lhacode = [ 999000010 ])

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')
