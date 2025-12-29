ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     
C     NB: VECSIZE_MEMMAX is defined in vector.inc
C     NB: vector.inc must be included before coupl.inc
C     

      DOUBLE PRECISION G, ALL_G(1)
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R(1)
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


      DOUBLE COMPLEX GC_2, GC_3, GC_50, GC_51, GC_58, GC_59, GC_93,
     $  GC_94

      DOUBLE COMPLEX GC_11(VECSIZE_MEMMAX)

      COMMON/COUPLINGS/ GC_11, GC_2, GC_3, GC_50, GC_51, GC_58, GC_59,
     $  GC_93, GC_94

