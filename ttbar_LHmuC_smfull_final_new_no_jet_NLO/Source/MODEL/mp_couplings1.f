ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MP_COUP1( )

      IMPLICIT NONE

      INCLUDE 'model_functions.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      MP__GC_1 = -(MP__MDL_EE*MP__MDL_COMPLEXI)/3.000000E+00_16
      MP__GC_2 = (2.000000E+00_16*MP__MDL_EE*MP__MDL_COMPLEXI)
     $ /3.000000E+00_16
      MP__GC_3 = -(MP__MDL_EE*MP__MDL_COMPLEXI)
      MP__GC_21 = -(MP__MDL_CW*MP__MDL_EE*MP__MDL_COMPLEXI)/(2.000000E
     $ +00_16*MP__MDL_SW)
      MP__GC_22 = (MP__MDL_CW*MP__MDL_EE*MP__MDL_COMPLEXI)/(2.000000E
     $ +00_16*MP__MDL_SW)
      MP__GC_23 = -(MP__MDL_EE*MP__MDL_COMPLEXI*MP__MDL_SW)/(6.000000E
     $ +00_16*MP__MDL_CW)
      MP__GC_24 = (MP__MDL_EE*MP__MDL_COMPLEXI*MP__MDL_SW)/(2.000000E
     $ +00_16*MP__MDL_CW)
      MP__GC_33 = -((MP__MDL_COMPLEXI*MP__MDL_YB)/MP__MDL_SQRT__2)
      MP__GC_34 = -((MP__MDL_COMPLEXI*MP__MDL_YC)/MP__MDL_SQRT__2)
      MP__GC_36 = -((MP__MDL_COMPLEXI*MP__MDL_YM)/MP__MDL_SQRT__2)
      MP__GC_37 = -((MP__MDL_COMPLEXI*MP__MDL_YT)/MP__MDL_SQRT__2)
      END
