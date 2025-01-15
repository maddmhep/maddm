ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_114 = -((MDL_GPU33*MDL_YT)/MDL_SQRT__2)
      GC_115 = (MDL_COMPLEXI*MDL_GSU33*MDL_YT)/MDL_SQRT__2
      GC_13 = -MDL_GPXD
      GC_14 = MDL_COMPLEXI*MDL_GSXD
      GC_1__1 = -(MDL_COMPLEXI*MDL_FSDE)
      GC_2__1 = MDL_COMPLEXI*MDL_FSDE
      GC_3__1 = MDL_COMPLEXI*MDL_FSIE
      GC_7__1 = MDL_COMPLEXI*MDL_GDMA
      GC_8__1 = MDL_COMPLEXI*MDL_GDMS
      END
