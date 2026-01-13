      INTEGER IFR
      INTEGER LMAXCONFIGSB_USED
      PARAMETER (LMAXCONFIGSB_USED=6)
      INTEGER MAX_BRANCHB_USED
      PARAMETER (MAX_BRANCHB_USED=3)
      INTEGER MAPCONFIG(0 : LMAXCONFIGSB_USED)
      INTEGER IFOREST(2, -MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER SPROP(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER TPRID(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      LOGICAL GFORCEBW(-MAX_BRANCHB_USED : -1, LMAXCONFIGSB_USED)

C     Here are the congifurations
C     Process: g mu+ > mu+ t t~ [ all = QCD QED ] QCD^2=4 QED^2=4
C     Diagram 1, Amplitude 1
      DATA MAPCONFIG(   1)/   1/
      DATA (IFOREST(IFR, -1,   1),IFR=1,2)/  1,  4/
      DATA TPRID(  -1,   1)/       6/
      DATA (IFOREST(IFR, -2,   1),IFR=1,2)/ -1,  5/
      DATA TPRID(  -2,   1)/      22/
      DATA (IFOREST(IFR, -3,   1),IFR=1,2)/ -2,  3/
C     Diagram 2, Amplitude 2
      DATA MAPCONFIG(   2)/   2/
      DATA (IFOREST(IFR, -1,   2),IFR=1,2)/  1,  4/
      DATA TPRID(  -1,   2)/       6/
      DATA (IFOREST(IFR, -2,   2),IFR=1,2)/ -1,  5/
      DATA TPRID(  -2,   2)/      25/
      DATA (IFOREST(IFR, -3,   2),IFR=1,2)/ -2,  3/
C     Diagram 3, Amplitude 3
      DATA MAPCONFIG(   3)/   3/
      DATA (IFOREST(IFR, -1,   3),IFR=1,2)/  1,  4/
      DATA TPRID(  -1,   3)/       6/
      DATA (IFOREST(IFR, -2,   3),IFR=1,2)/ -1,  5/
      DATA TPRID(  -2,   3)/      23/
      DATA (IFOREST(IFR, -3,   3),IFR=1,2)/ -2,  3/
C     Diagram 4, Amplitude 4
      DATA MAPCONFIG(   4)/   4/
      DATA (IFOREST(IFR, -1,   4),IFR=1,2)/  1,  5/
      DATA TPRID(  -1,   4)/       6/
      DATA (IFOREST(IFR, -2,   4),IFR=1,2)/ -1,  4/
      DATA TPRID(  -2,   4)/      22/
      DATA (IFOREST(IFR, -3,   4),IFR=1,2)/ -2,  3/
C     Diagram 5, Amplitude 5
      DATA MAPCONFIG(   5)/   5/
      DATA (IFOREST(IFR, -1,   5),IFR=1,2)/  1,  5/
      DATA TPRID(  -1,   5)/       6/
      DATA (IFOREST(IFR, -2,   5),IFR=1,2)/ -1,  4/
      DATA TPRID(  -2,   5)/      25/
      DATA (IFOREST(IFR, -3,   5),IFR=1,2)/ -2,  3/
C     Diagram 6, Amplitude 6
      DATA MAPCONFIG(   6)/   6/
      DATA (IFOREST(IFR, -1,   6),IFR=1,2)/  1,  5/
      DATA TPRID(  -1,   6)/       6/
      DATA (IFOREST(IFR, -2,   6),IFR=1,2)/ -1,  4/
      DATA TPRID(  -2,   6)/      23/
      DATA (IFOREST(IFR, -3,   6),IFR=1,2)/ -2,  3/
C     Number of configs
      DATA MAPCONFIG(0)/   6/

C     Here are the BWs
