ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1( )

      IMPLICIT NONE

      INCLUDE 'model_functions.inc'
      INCLUDE '../vector.inc'


      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_40 = (MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_106 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM3X1)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_107 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM3X2)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_108 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM3X3)/(MDL_SW
     $ *MDL_SQRT__2)
      END
