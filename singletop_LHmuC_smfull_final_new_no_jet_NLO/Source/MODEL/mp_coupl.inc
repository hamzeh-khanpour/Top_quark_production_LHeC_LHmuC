ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*16 MP__G
      COMMON/MP_STRONG/ MP__G

      COMPLEX*32 MP__GAL(2)
      COMMON/MP_WEAK/ MP__GAL

      COMPLEX*32 MP__MU_R
      COMMON/MP_RSCALE/ MP__MU_R


      REAL*16 MP__MDL_MB,MP__MDL_MC,MP__MDL_MH,MP__MDL_MM,MP__MDL_MT
     $ ,MP__MDL_MTA,MP__MDL_MW,MP__MDL_MZ,MP__MDL_ME

      COMMON/MP_MASSES/ MP__MDL_MB,MP__MDL_MC,MP__MDL_MH,MP__MDL_MM
     $ ,MP__MDL_MT,MP__MDL_MTA,MP__MDL_MW,MP__MDL_MZ,MP__MDL_ME


      REAL*16 MP__MDL_WH,MP__MDL_WT,MP__MDL_WTAU,MP__MDL_WW,MP__MDL_WZ

      COMMON/MP_WIDTHS/ MP__MDL_WH,MP__MDL_WT,MP__MDL_WTAU,MP__MDL_WW
     $ ,MP__MDL_WZ


      COMPLEX*32 MP__GC_11,MP__GC_18,MP__GC_19,MP__GC_20

      COMPLEX*32 MP__GC_5,MP__UVWFCT_B_0_1EPS,MP__UVWFCT_T_0_1EPS
     $ ,MP__R2_DXTW,MP__R2_SXTW,MP__R2_BXTW,MP__UVWFCT_B_0
     $ ,MP__UVWFCT_T_0

      COMMON/MP_COUPLINGS/ MP__GC_5,MP__UVWFCT_B_0_1EPS
     $ ,MP__UVWFCT_T_0_1EPS,MP__R2_DXTW,MP__R2_SXTW,MP__R2_BXTW
     $ ,MP__UVWFCT_B_0,MP__UVWFCT_T_0,MP__GC_11,MP__GC_18,MP__GC_19
     $ ,MP__GC_20


