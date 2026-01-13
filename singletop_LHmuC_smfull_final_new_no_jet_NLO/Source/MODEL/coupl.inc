ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     
C     NB: VECSIZE_MEMMAX is defined in vector.inc
C     NB: vector.inc must be included before coupl.inc
C     

      DOUBLE PRECISION G, ALL_G
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R
      COMMON/RSCALE/ MU_R, ALL_MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=3D0)
      DOUBLE PRECISION NL
      PARAMETER(NL=0D0)

      DOUBLE PRECISION MDL_MB,MDL_MC,MDL_MH,MDL_MM,MDL_MT,MDL_MTA
     $ ,MDL_MW,MDL_MZ,MDL_ME

      COMMON/MASSES/ MDL_MB,MDL_MC,MDL_MH,MDL_MM,MDL_MT,MDL_MTA,MDL_MW
     $ ,MDL_MZ,MDL_ME


      DOUBLE PRECISION MDL_WH,MDL_WT,MDL_WTAU,MDL_WW,MDL_WZ

      COMMON/WIDTHS/ MDL_WH,MDL_WT,MDL_WTAU,MDL_WW,MDL_WZ


      DOUBLE COMPLEX GC_11, GC_18, GC_19, GC_20

      DOUBLE COMPLEX GC_5, UVWFCT_B_0_1EPS, UVWFCT_T_0_1EPS, R2_DXTW,
     $  R2_SXTW, R2_BXTW, UVWFCT_B_0, UVWFCT_T_0

      COMMON/COUPLINGS/ GC_5, UVWFCT_B_0_1EPS, UVWFCT_T_0_1EPS,
     $  R2_DXTW, R2_SXTW, R2_BXTW, UVWFCT_B_0, UVWFCT_T_0, GC_11,
     $  GC_18, GC_19, GC_20

