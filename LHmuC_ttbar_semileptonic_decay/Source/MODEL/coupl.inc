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


      DOUBLE COMPLEX GC_2, GC_3, GC_40, GC_41, GC_42, GC_43, GC_44,
     $  GC_45, GC_46, GC_49, GC_50, GC_51, GC_58, GC_59, GC_93, GC_94,
     $  GC_100, GC_101, GC_102, GC_103, GC_104, GC_105, GC_108

      DOUBLE COMPLEX GC_11(VECSIZE_MEMMAX)

      COMMON/COUPLINGS/ GC_11, GC_2, GC_3, GC_40, GC_41, GC_42, GC_43,
     $  GC_44, GC_45, GC_46, GC_49, GC_50, GC_51, GC_58, GC_59, GC_93,
     $  GC_94, GC_100, GC_101, GC_102, GC_103, GC_104, GC_105, GC_108

