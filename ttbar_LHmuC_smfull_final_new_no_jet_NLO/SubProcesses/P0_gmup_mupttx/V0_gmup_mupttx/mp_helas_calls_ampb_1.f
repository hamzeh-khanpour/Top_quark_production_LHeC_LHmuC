      SUBROUTINE MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=32)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=6)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=49, NLOOPGROUPS=15, NCTAMPS=115)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=164)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=17,NLOOPWAVEFUNCS=95)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=2, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     

      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_OXXXXX(P(0,2),MDL_MM,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_IXXXXX(P(0,3),MDL_MM,NHEL(3),-1*IC(3),W(1,3))
      CALL MP_OXXXXX(P(0,4),MDL_MT,NHEL(4),+1*IC(4),W(1,4))
      CALL MP_IXXXXX(P(0,5),MDL_MT,NHEL(5),-1*IC(5),W(1,5))
      CALL MP_FFV1_1(W(1,4),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,6))
      CALL MP_FFV1P0_3(W(1,3),W(1,2),GC_3,ZERO,ZERO,W(1,7))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_FFV1_0(W(1,5),W(1,6),W(1,7),GC_2,AMP(1))
      CALL MP_FFS1_3(W(1,3),W(1,2),GC_36,MDL_MH,MDL_WH,W(1,8))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_FFS1_0(W(1,5),W(1,6),W(1,8),GC_37,AMP(2))
      CALL MP_FFV2_4_3(W(1,3),W(1,2),GC_21,GC_24,MDL_MZ,MDL_WZ,W(1,9))
C     Amplitude(s) for born diagram with ID 3
      CALL MP_FFV2_5_0(W(1,5),W(1,6),W(1,9),GC_22,GC_23,AMP(3))
      CALL MP_FFV1_2(W(1,5),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,10))
C     Amplitude(s) for born diagram with ID 4
      CALL MP_FFV1_0(W(1,10),W(1,4),W(1,7),GC_2,AMP(4))
C     Amplitude(s) for born diagram with ID 5
      CALL MP_FFS1_0(W(1,10),W(1,4),W(1,8),GC_37,AMP(5))
C     Amplitude(s) for born diagram with ID 6
      CALL MP_FFV2_5_0(W(1,10),W(1,4),W(1,9),GC_22,GC_23,AMP(6))
      CALL MP_FFV1_2(W(1,5),W(1,7),GC_2,MDL_MT,MDL_WT,W(1,11))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,11),W(1,6),R2_QQQ,R2_QQT,AMPL(1,1))
      CALL MP_R2_QQ_2_0(W(1,11),W(1,6),UV_TMASS_1EPS,AMPL(2,2))
      CALL MP_R2_QQ_2_0(W(1,11),W(1,6),UV_TMASS,AMPL(1,3))
      CALL MP_FFS1_2(W(1,5),W(1,8),GC_37,MDL_MT,MDL_WT,W(1,12))
C     Counter-term amplitude(s) for loop diagram number 8
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,12),W(1,6),R2_QQQ,R2_QQT,AMPL(1,4))
      CALL MP_R2_QQ_2_0(W(1,12),W(1,6),UV_TMASS_1EPS,AMPL(2,5))
      CALL MP_R2_QQ_2_0(W(1,12),W(1,6),UV_TMASS,AMPL(1,6))
      CALL MP_FFV2_5_2(W(1,5),W(1,9),GC_22,GC_23,MDL_MT,MDL_WT,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,13),W(1,6),R2_QQQ,R2_QQT,AMPL(1,7))
      CALL MP_R2_QQ_2_0(W(1,13),W(1,6),UV_TMASS_1EPS,AMPL(2,8))
      CALL MP_R2_QQ_2_0(W(1,13),W(1,6),UV_TMASS,AMPL(1,9))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL MP_FFV1_0(W(1,5),W(1,6),W(1,7),R2_UUA,AMPL(1,10))
C     Counter-term amplitude(s) for loop diagram number 11
      CALL MP_FFS1_0(W(1,5),W(1,6),W(1,8),R2_TTH,AMPL(1,11))
      CALL MP_FFS1_0(W(1,5),W(1,6),W(1,8),UV_HTT_1EPS,AMPL(2,12))
      CALL MP_FFS1_0(W(1,5),W(1,6),W(1,8),UV_HTT,AMPL(1,13))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL MP_FFV2_5_0(W(1,5),W(1,6),W(1,9),R2_UUZ_V2,R2_UUZ_V5,AMPL(1
     $ ,14))
      CALL MP_FFV1_1(W(1,4),W(1,7),GC_2,MDL_MT,MDL_WT,W(1,14))
C     Counter-term amplitude(s) for loop diagram number 13
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,14),R2_QQQ,R2_QQT,AMPL(1
     $ ,15))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,14),UV_TMASS_1EPS,AMPL(2,16))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,14),UV_TMASS,AMPL(1,17))
      CALL MP_FFS1_1(W(1,4),W(1,8),GC_37,MDL_MT,MDL_WT,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,15),R2_QQQ,R2_QQT,AMPL(1
     $ ,18))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,15),UV_TMASS_1EPS,AMPL(2,19))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,15),UV_TMASS,AMPL(1,20))
      CALL MP_FFV2_5_1(W(1,4),W(1,9),GC_22,GC_23,MDL_MT,MDL_WT,W(1,16))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,16),R2_QQQ,R2_QQT,AMPL(1
     $ ,21))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,16),UV_TMASS_1EPS,AMPL(2,22))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,16),UV_TMASS,AMPL(1,23))
C     Counter-term amplitude(s) for loop diagram number 16
      CALL MP_FFV1_0(W(1,10),W(1,4),W(1,7),R2_UUA,AMPL(1,24))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL MP_FFS1_0(W(1,10),W(1,4),W(1,8),R2_TTH,AMPL(1,25))
      CALL MP_FFS1_0(W(1,10),W(1,4),W(1,8),UV_HTT_1EPS,AMPL(2,26))
      CALL MP_FFS1_0(W(1,10),W(1,4),W(1,8),UV_HTT,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL MP_FFV2_5_0(W(1,10),W(1,4),W(1,9),R2_UUZ_V2,R2_UUZ_V5
     $ ,AMPL(1,28))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),R2_GQQ,AMPL(1,29))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQB_1EPS,AMPL(2,30))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQC_1EPS,AMPL(2,31))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,32))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,33))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,34))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,35))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQT_1EPS,AMPL(2,36))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQB,AMPL(1,37))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQC,AMPL(1,38))
      CALL MP_FFV1_0(W(1,11),W(1,4),W(1,1),UV_GQQT,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 20
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),R2_GQQ,AMPL(1,40))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQB_1EPS,AMPL(2,41))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQC_1EPS,AMPL(2,42))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQG_1EPS,AMPL(2,43))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQQ_1EPS,AMPL(2,44))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQQ_1EPS,AMPL(2,45))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQQ_1EPS,AMPL(2,46))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQT_1EPS,AMPL(2,47))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQB,AMPL(1,48))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQC,AMPL(1,49))
      CALL MP_FFV1_0(W(1,5),W(1,14),W(1,1),UV_GQQT,AMPL(1,50))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),R2_GQQ,AMPL(1,51))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQB_1EPS,AMPL(2,52))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQC_1EPS,AMPL(2,53))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,54))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,55))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,56))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,57))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQT_1EPS,AMPL(2,58))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQB,AMPL(1,59))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQC,AMPL(1,60))
      CALL MP_FFV1_0(W(1,12),W(1,4),W(1,1),UV_GQQT,AMPL(1,61))
C     Counter-term amplitude(s) for loop diagram number 22
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),R2_GQQ,AMPL(1,62))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQB_1EPS,AMPL(2,63))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQC_1EPS,AMPL(2,64))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQG_1EPS,AMPL(2,65))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQQ_1EPS,AMPL(2,66))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQQ_1EPS,AMPL(2,67))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQQ_1EPS,AMPL(2,68))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQT_1EPS,AMPL(2,69))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQB,AMPL(1,70))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQC,AMPL(1,71))
      CALL MP_FFV1_0(W(1,5),W(1,15),W(1,1),UV_GQQT,AMPL(1,72))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),R2_GQQ,AMPL(1,73))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQB_1EPS,AMPL(2,74))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQC_1EPS,AMPL(2,75))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,76))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,77))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,78))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,79))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQT_1EPS,AMPL(2,80))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQB,AMPL(1,81))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQC,AMPL(1,82))
      CALL MP_FFV1_0(W(1,13),W(1,4),W(1,1),UV_GQQT,AMPL(1,83))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),R2_GQQ,AMPL(1,84))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQB_1EPS,AMPL(2,85))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQC_1EPS,AMPL(2,86))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQG_1EPS,AMPL(2,87))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQQ_1EPS,AMPL(2,88))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQQ_1EPS,AMPL(2,89))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQQ_1EPS,AMPL(2,90))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQT_1EPS,AMPL(2,91))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQB,AMPL(1,92))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQC,AMPL(1,93))
      CALL MP_FFV1_0(W(1,5),W(1,16),W(1,1),UV_GQQT,AMPL(1,94))
      CALL MP_FFV1P0_3(W(1,5),W(1,4),GC_5,ZERO,ZERO,W(1,17))
C     Counter-term amplitude(s) for loop diagram number 40
      CALL MP_R2_GGZ_0(W(1,1),W(1,17),W(1,9),R2_GGZDOWN,AMPL(1,95))
      CALL MP_R2_GGZ_0(W(1,1),W(1,17),W(1,9),R2_GGZDOWN,AMPL(1,96))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL MP_R2_GGZ_0(W(1,1),W(1,17),W(1,9),R2_GGZUP,AMPL(1,97))
C     Counter-term amplitude(s) for loop diagram number 44
      CALL MP_VVS1_0(W(1,1),W(1,17),W(1,8),R2_GGHC,AMPL(1,98))
C     Counter-term amplitude(s) for loop diagram number 46
      CALL MP_R2_GGZ_0(W(1,1),W(1,17),W(1,9),R2_GGZUP,AMPL(1,99))
C     Counter-term amplitude(s) for loop diagram number 48
      CALL MP_VVS1_0(W(1,1),W(1,17),W(1,8),R2_GGHB,AMPL(1,100))
C     Counter-term amplitude(s) for loop diagram number 50
      CALL MP_R2_GGZ_0(W(1,1),W(1,17),W(1,9),R2_GGZDOWN,AMPL(1,101))
C     Counter-term amplitude(s) for loop diagram number 52
      CALL MP_VVS1_0(W(1,1),W(1,17),W(1,8),R2_GGHT,AMPL(1,102))
C     Counter-term amplitude(s) for loop diagram number 54
      CALL MP_R2_GGZ_0(W(1,1),W(1,17),W(1,9),R2_GGZUP,AMPL(1,103))
C     At this point, all CT amps needed for (QCD=4 QED=4), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

