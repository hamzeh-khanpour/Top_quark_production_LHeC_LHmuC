ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*16 G
      COMMON/MP_STRONG/ G

      COMPLEX*32 GAL(2)
      COMMON/MP_WEAK/ GAL

      COMPLEX*32 MU_R
      COMMON/MP_RSCALE/ MU_R


      REAL*16 MDL_MB,MDL_MC,MDL_MH,MDL_MM,MDL_MT,MDL_MTA,MDL_MW,MDL_MZ
     $ ,MDL_ME

      COMMON/MP_MASSES/ MDL_MB,MDL_MC,MDL_MH,MDL_MM,MDL_MT,MDL_MTA
     $ ,MDL_MW,MDL_MZ,MDL_ME


      REAL*16 MDL_WH,MDL_WT,MDL_WTAU,MDL_WW,MDL_WZ

      COMMON/MP_WIDTHS/ MDL_WH,MDL_WT,MDL_WTAU,MDL_WW,MDL_WZ


      COMPLEX*32 GC_11,GC_18,GC_19,GC_20

      COMPLEX*32 GC_5,UVWFCT_B_0_1EPS,UVWFCT_T_0_1EPS,R2_DXTW,R2_SXTW
     $ ,R2_BXTW,UVWFCT_B_0,UVWFCT_T_0

      COMMON/MP_COUPLINGS/ GC_5,UVWFCT_B_0_1EPS,UVWFCT_T_0_1EPS
     $ ,R2_DXTW,R2_SXTW,R2_BXTW,UVWFCT_B_0,UVWFCT_T_0,GC_11,GC_18
     $ ,GC_19,GC_20


