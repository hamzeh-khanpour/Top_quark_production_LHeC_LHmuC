ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP3( )

      IMPLICIT NONE

      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_5 = MDL_COMPLEXI*G
      UVWFCT_B_0_1EPS = COND(DCMPLX(MDL_MB),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))
     $ *3.000000D+00*MDL_CF))
      UVWFCT_T_0_1EPS = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))
     $ *3.000000D+00*MDL_CF))
      R2_DXTW = ((MDL_CKM31*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2))
     $ *MDL_R2MIXEDFACTOR_FIN_
      R2_SXTW = ((MDL_CKM32*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2))
     $ *MDL_R2MIXEDFACTOR_FIN_
      R2_BXTW = ((MDL_CKM33*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2))
     $ *MDL_R2MIXEDFACTOR_FIN_
      UVWFCT_B_0 = COND(DCMPLX(MDL_MB),DCMPLX(0.000000D+00),DCMPLX(
     $ -((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))*MDL_CF
     $ *(4.000000D+00-3.000000D+00*REGLOG(DCMPLX(MDL_MB__EXP__2
     $ /MDL_MU_R__EXP__2)))))
      UVWFCT_T_0 = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00),DCMPLX(
     $ -((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))*MDL_CF
     $ *(4.000000D+00-3.000000D+00*REGLOG(DCMPLX(MDL_MT__EXP__2
     $ /MDL_MU_R__EXP__2)))))
      END
