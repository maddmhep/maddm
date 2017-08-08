      SUBROUTINE CACAI(Z, FNU, KODE, MR, N, Y, NZ, RL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CACAI
C***REFER TO  CAIRY
C
C     CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO CACON CAN RESULT IF CACON
C     IS CALLED FROM CAIRY.
C
C***ROUTINES CALLED  CASYI,CBKNU,CMLRI,CSERI,CS1S2,R1MACH
C***END PROLOGUE  CACAI
      COMPLEX CSGN, CSPN, C1, C2, Y, Z, ZN, CY
      REAL ALIM, ARG, ASCLE, AZ, CPN, DFNU, ELIM, FMR, FNU, PI, RL,
     * SGN, SPN, TOL, YY, R1MACH
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION Y(N), CY(2)
      DATA PI / 3.14159265358979324E0 /
      NZ = 0
      ZN = -Z
      AZ = CABS(Z)
      NN = N
      DFNU = FNU + FLOAT(N-1)
      IF (AZ.LE.2.0E0) GO TO 10
      IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL CSERI(ZN, FNU, KODE, NN, Y, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL CASYI(ZN, FNU, KODE, NN, Y, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 70
      GO TO 40
   30 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL CMLRI(ZN, FNU, KODE, NN, Y, NW, TOL)
      IF(NW.LT.0) GO TO 70
   40 CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      CALL CBKNU(ZN, FNU, KODE, 1, CY, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 70
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
      CSGN = CMPLX(0.0E0,SGN)
      IF (KODE.EQ.1) GO TO 50
      YY = -AIMAG(ZN)
      CPN = COS(YY)
      SPN = SIN(YY)
      CSGN = CSGN*CMPLX(CPN,SPN)
   50 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(INU,2).EQ.1) CSPN = -CSPN
      C1 = CY(1)
      C2 = Y(1)
      IF (KODE.EQ.1) GO TO 60
      IUF = 0
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   60 CONTINUE
      Y(1) = CSPN*C1 + CSGN*C2
      RETURN
   70 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE CACON(Z, FNU, KODE, MR, N, Y, NZ, RL, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CACON
C***REFER TO  CBESK,CBESH
C
C     CACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***ROUTINES CALLED  CBINU,CBKNU,CS1S2,R1MACH
C***END PROLOGUE  CACON
      COMPLEX CK, CONE, CS, CSCL, CSCR, CSGN, CSPN, CSS, CSR, C1, C2,
     * RZ, SC1, SC2, ST, S1, S2, Y, Z, ZN, CY
      REAL ALIM, ARG, ASCLE, AS2, BSCLE, BRY, CPN, C1I, C1M, C1R, ELIM,
     * FMR, FNU, FNUL, PI, RL, SGN, SPN, TOL, YY, R1MACH
      INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
      DIMENSION Y(N), CY(2), CSS(3), CSR(3), BRY(3)
      DATA PI / 3.14159265358979324E0 /
      DATA CONE / (1.0E0,0.0E0) /
      NZ = 0
      ZN = -Z
      NN = N
      CALL CBINU(ZN, FNU, KODE, NN, Y, NW, RL, FNUL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 80
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      NN = MIN0(2,N)
      CALL CBKNU(ZN, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      S1 = CY(1)
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
      CSGN = CMPLX(0.0E0,SGN)
      IF (KODE.EQ.1) GO TO 10
      YY = -AIMAG(ZN)
      CPN = COS(YY)
      SPN = SIN(YY)
      CSGN = CSGN*CMPLX(CPN,SPN)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(INU,2).EQ.1) CSPN = -CSPN
      IUF = 0
      C1 = S1
      C2 = Y(1)
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      IF (KODE.EQ.1) GO TO 20
      CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC1 = C1
   20 CONTINUE
      Y(1) = CSPN*C1 + CSGN*C2
      IF (N.EQ.1) RETURN
      CSPN = -CSPN
      S2 = CY(2)
      C1 = S2
      C2 = Y(2)
      IF (KODE.EQ.1) GO TO 30
      CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC2 = C1
   30 CONTINUE
      Y(2) = CSPN*C1 + CSGN*C2
      IF (N.EQ.2) RETURN
      CSPN = -CSPN
      RZ = CMPLX(2.0E0,0.0E0)/ZN
      CK = CMPLX(FNU+1.0E0,0.0E0)*RZ
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CSCR = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CSCR
      CSR(1) = CSCR
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = ASCLE
      BRY(2) = 1.0E0/ASCLE
      BRY(3) = R1MACH(2)
      AS2 = CABS(S2)
      KFLAG = 2
      IF (AS2.GT.BRY(1)) GO TO 40
      KFLAG = 1
      GO TO 50
   40 CONTINUE
      IF (AS2.LT.BRY(2)) GO TO 50
      KFLAG = 3
   50 CONTINUE
      BSCLE = BRY(KFLAG)
      S1 = S1*CSS(KFLAG)
      S2 = S2*CSS(KFLAG)
      CS = CSR(KFLAG)
      DO 70 I=3,N
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        C1 = S2*CS
        ST = C1
        C2 = Y(I)
        IF (KODE.EQ.1) GO TO 60
        IF (IUF.LT.0) GO TO 60
        CALL CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
        NZ = NZ + NW
        SC1 = SC2
        SC2 = C1
        IF (IUF.NE.3) GO TO 60
        IUF = -4
        S1 = SC1*CSS(KFLAG)
        S2 = SC2*CSS(KFLAG)
        ST = SC2
   60   CONTINUE
        Y(I) = CSPN*C1 + CSGN*C2
        CK = CK + RZ
        CSPN = -CSPN
        IF (KFLAG.GE.3) GO TO 70
        C1R = REAL(C1)
        C1I = AIMAG(C1)
        C1R = ABS(C1R)
        C1I = ABS(C1I)
        C1M = AMAX1(C1R,C1I)
        IF (C1M.LE.BSCLE) GO TO 70
        KFLAG = KFLAG + 1
        BSCLE = BRY(KFLAG)
        S1 = S1*CS
        S2 = ST
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        CS = CSR(KFLAG)
   70 CONTINUE
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE CAIRY(Z, ID, KODE, AI, NZ, IERR)
C***BEGIN PROLOGUE  CAIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C         ON KODE=1, CAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
C         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
C         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
C         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
C         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z)
C
C         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
C         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
C         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
C         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             AI=AI(Z)                ON ID=0 OR
C                             AI=DAI(Z)/DZ            ON ID=1
C                        = 2  RETURNS
C                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
C                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)
C
C         OUTPUT
C           AI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           NZ     - UNDERFLOW INDICATOR
C                    NZ= 0   , NORMAL RETURN
C                    NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
C                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
C                            TOO LARGE WITH KODE=1.
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C
C***LONG DESCRIPTION
C
C         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
C         FUNCTIONS BY
C
C            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
C                           C=1.0/(PI*SQRT(3.0))
C                           ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 1986
C
C***ROUTINES CALLED  CACAI,CBKNU,I1MACH,R1MACH
C***END PROLOGUE  CAIRY
      COMPLEX AI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
      REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CK, COEF, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, RL, R1M5, SFAC, TOL, TTH, ZI, ZR,
     * Z3I, Z3R, R1MACH, BB, ALAZ
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CY(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667E-01,
     * 3.55028053887817240E-01,2.58819403792806799E-01,
     * 1.83776298473930683E-01/
      DATA  CONE / (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = CABS(Z)
      TOL = AMAX1(R1MACH(4),1.0E-18)
      FID = FLOAT(ID)
      IF (AZ.GT.1.0E0) GO TO 60
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1 = CONE
      S2 = CONE
      IF (AZ.LT.TOL) GO TO 160
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1 = CONE
      TRM2 = CONE
      ATRM = 1.0E0
      Z3 = Z*Z*Z
      AZ3 = AZ*AA
      AK = 2.0E0 + FID
      BK = 3.0E0 - FID - FID
      CK = 4.0E0 - FID
      DK = 3.0E0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = AMIN1(D1,D2)
      AK = 24.0E0 + 9.0E0*FID
      BK = 30.0E0 - 9.0E0*FID
      Z3R = REAL(Z3)
      Z3I = AIMAG(Z3)
      DO 30 K=1,25
        TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
        S1 = S1 + TRM1
        TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
        S2 = S2 + TRM2
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = AMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0E0
        BK = BK + 18.0E0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AI = AI*CEXP(ZTA)
      RETURN
   50 CONTINUE
      AI = -S2*CMPLX(C2,0.0E0)
      IF (AZ.GT.TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AI = AI*CEXP(ZTA)
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   60 CONTINUE
      FNU = (1.0E0+FID)/3.0E0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      ALAZ=ALOG(AZ)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA=0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CSQ=CSQRT(Z)
      ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0E0
      ZI = AIMAG(Z)
      ZR = REAL(Z)
      AK = AIMAG(ZTA)
      IF (ZR.GE.0.0E0) GO TO 70
      BK = REAL(ZTA)
      CK = -ABS(BK)
      ZTA = CMPLX(CK,AK)
   70 CONTINUE
      IF (ZI.NE.0.0E0) GO TO 80
      IF (ZR.GT.0.0E0) GO TO 80
      ZTA = CMPLX(0.0E0,AK)
   80 CONTINUE
      AA = REAL(ZTA)
      IF (AA.GE.0.0E0 .AND. ZR.GT.0.0E0) GO TO 100
      IF (KODE.EQ.2) GO TO 90
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 90
      AA = -AA + 0.25E0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 240
   90 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0E0) MR = -1
      CALL CACAI(ZTA, FNU, KODE, MR, 1, CY, NN, RL, TOL, ELIM, ALIM)
      IF (NN.LT.0) GO TO 250
      NZ = NZ + NN
      GO TO 120
  100 CONTINUE
      IF (KODE.EQ.2) GO TO 110
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 110
      AA = -AA - 0.25E0*ALAZ
      IFLAG = 2
      SFAC = 1.0E0/TOL
      IF (AA.LT.(-ELIM)) GO TO 180
  110 CONTINUE
      CALL CBKNU(ZTA, FNU, KODE, 1, CY, NZ, TOL, ELIM, ALIM)
  120 CONTINUE
      S1 = CY(1)*CMPLX(COEF,0.0E0)
      IF (IFLAG.NE.0) GO TO 140
      IF (ID.EQ.1) GO TO 130
      AI = CSQ*S1
      RETURN
  130 AI = -Z*S1
      RETURN
  140 CONTINUE
      S1 = S1*CMPLX(SFAC,0.0E0)
      IF (ID.EQ.1) GO TO 150
      S1 = S1*CSQ
      AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  150 CONTINUE
      S1 = -S1*Z
      AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  160 CONTINUE
      AA = 1.0E+3*R1MACH(1)
      S1 = CMPLX(0.0E0,0.0E0)
      IF (ID.EQ.1) GO TO 170
      IF (AZ.GT.AA) S1 = CMPLX(C2,0.0E0)*Z
      AI = CMPLX(C1,0.0E0) - S1
      RETURN
  170 CONTINUE
      AI = -CMPLX(C2,0.0E0)
      AA = SQRT(AA)
      IF (AZ.GT.AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
      AI = AI + S1*CMPLX(C1,0.0E0)
      RETURN
  180 CONTINUE
      NZ = 1
      AI = CMPLX(0.0E0,0.0E0)
      RETURN
  240 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  250 CONTINUE
      IF(NN.EQ.(-1)) GO TO 240
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
      SUBROUTINE CASYI(Z, FNU, KODE, N, Y, NZ, RL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CASYI
C***REFER TO  CBESI,CBESK
C
C     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***ROUTINES CALLED  R1MACH
C***END PROLOGUE  CASYI
      COMPLEX AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1, RZ, S2,
     * Y, Z
      REAL AA, ACZ, AEZ, AK, ALIM, ARG, ARM, ATOL, AZ, BB, BK, DFNU,
     * DNU2, ELIM, FDN, FNU, PI, RL, RTPI, RTR1, S, SGN, SQK, TOL, X,
     * YY, R1MACH
      INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
      DIMENSION Y(N)
      DATA PI, RTPI  /3.14159265358979324E0 , 0.159154943091895336E0 /
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      NZ = 0
      AZ = CABS(Z)
      X = REAL(Z)
      ARM = 1.0E+3*R1MACH(1)
      RTR1 = SQRT(ARM)
      IL = MIN0(2,N)
      DFNU = FNU + FLOAT(N-IL)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      AK1 = CMPLX(RTPI,0.0E0)/Z
      AK1 = CSQRT(AK1)
      CZ = Z
      IF (KODE.EQ.2) CZ = Z - CMPLX(X,0.0E0)
      ACZ = REAL(CZ)
      IF (ABS(ACZ).GT.ELIM) GO TO 80
      DNU2 = DFNU + DFNU
      KODED = 1
      IF ((ABS(ACZ).GT.ALIM) .AND. (N.GT.2)) GO TO 10
      KODED = 0
      AK1 = AK1*CEXP(CZ)
   10 CONTINUE
      FDN = 0.0E0
      IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
      EZ = Z*CMPLX(8.0E0,0.0E0)
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
      AEZ = 8.0E0*AZ
      S = TOL/AEZ
      JL = INT(RL+RL) + 2
      YY = AIMAG(Z)
      P1 = CZERO
      IF (YY.EQ.0.0E0) GO TO 20
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-FLOAT(INU))*PI
      INU = INU + N - IL
      AK = -SIN(ARG)
      BK = COS(ARG)
      IF (YY.LT.0.0E0) BK = -BK
      P1 = CMPLX(AK,BK)
      IF (MOD(INU,2).EQ.1) P1 = -P1
   20 CONTINUE
      DO 50 K=1,IL
        SQK = FDN - 1.0E0
        ATOL = S*ABS(SQK)
        SGN = 1.0E0
        CS1 = CONE
        CS2 = CONE
        CK = CONE
        AK = 0.0E0
        AA = 1.0E0
        BB = AEZ
        DK = EZ
        DO 30 J=1,JL
          CK = CK*CMPLX(SQK,0.0E0)/DK
          CS2 = CS2 + CK
          SGN = -SGN
          CS1 = CS1 + CK*CMPLX(SGN,0.0E0)
          DK = DK + EZ
          AA = AA*ABS(SQK)/BB
          BB = BB + AEZ
          AK = AK + 8.0E0
          SQK = SQK - AK
          IF (AA.LE.ATOL) GO TO 40
   30   CONTINUE
        GO TO 90
   40   CONTINUE
        S2 = CS1
        IF (X+X.LT.ELIM) S2 = S2 + P1*CS2*CEXP(-Z-Z)
        FDN = FDN + 8.0E0*DFNU + 4.0E0
        P1 = -P1
        M = N - IL + K
        Y(M) = S2*AK1
   50 CONTINUE
      IF (N.LE.2) RETURN
      NN = N
      K = NN - 2
      AK = FLOAT(K)
      RZ = (CONE+CONE)/Z
      IB = 3
      DO 60 I=IB,NN
        Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
        AK = AK - 1.0E0
        K = K - 1
   60 CONTINUE
      IF (KODED.EQ.0) RETURN
      CK = CEXP(CZ)
      DO 70 I=1,NN
        Y(I) = Y(I)*CK
   70 CONTINUE
      RETURN
   80 CONTINUE
      NZ = -1
      RETURN
   90 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE CBESK(Z, FNU, KODE, N, CY, NZ, IERR)
C***BEGIN PROLOGUE  CBESK
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
C             BESSEL FUNCTION OF THE THIRD KIND
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C***DESCRIPTION
C
C         ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
C         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESK
C         RETURNS THE SCALED K FUNCTIONS,
C
C         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
C
C         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0E0
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(I)=K(FNU+I-1,Z), I=1,...,N
C                        = 2  RETURNS
C                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
C                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C                    DEPENDING ON KODE
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
C                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
C                              IN THE SEQUENCE.
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
C         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
C         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
C         HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
C         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
C
C         FOR NEGATIVE ORDERS, THE FORMULA
C
C                       K(-FNU,Z) = K(FNU,Z)
C
C         CAN BE USED.
C
C         CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
C         AVAILABLE.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 1986
C
C***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
C***END PROLOGUE  CBESK
C
      COMPLEX CY, Z
      REAL AA, ALIM, ALN, ARG, AZ, DIG, ELIM, FN, FNU, FNUL, RL, R1M5,
     * TOL, UFL, XX, YY, R1MACH, BB
      INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CY(N)
C***FIRST EXECUTABLE STATEMENT  CBESK
      IERR = 0
      NZ=0
      XX = REAL(Z)
      YY = AIMAG(Z)
      IF (YY.EQ.0.0E0 .AND. XX.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0E0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = AMAX1(R1MACH(4),1.0E-18)
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303E0*(FLOAT(K)*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*FLOAT(K1)
      DIG = AMIN1(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + AMAX1(-AA,-41.45E0)
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      AZ = CABS(Z)
      FN = FNU + FLOAT(NN-1)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA = 0.5E0/TOL
      BB=FLOAT(I1MACH(9))*0.5E0
      AA=AMIN1(AA,BB)
      IF(AZ.GT.AA) GO TO 210
      IF(FN.GT.AA) GO TO 210
      AA=SQRT(AA)
      IF(AZ.GT.AA) IERR=3
      IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
C     UFL = EXP(-ELIM)
      UFL = R1MACH(1)*1.0E+3
      IF (AZ.LT.UFL) GO TO 180
      IF (FNU.GT.FNUL) GO TO 80
      IF (FN.LE.1.0E0) GO TO 60
      IF (FN.GT.2.0E0) GO TO 50
      IF (AZ.GT.TOL) GO TO 60
      ARG = 0.5E0*AZ
      ALN = -FN*ALOG(ARG)
      IF (ALN.GT.ELIM) GO TO 180
      GO TO 60
   50 CONTINUE
      CALL CUOIK(Z, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 180
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 100
   60 CONTINUE
      IF (XX.LT.0.0E0) GO TO 70
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
C-----------------------------------------------------------------------
      CALL CBKNU(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 200
      NZ=NW
      RETURN
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
C-----------------------------------------------------------------------
   70 CONTINUE
      IF (NZ.NE.0) GO TO 180
      MR = 1
      IF (YY.LT.0.0E0) MR = -1
      CALL CACON(Z, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 200
      NZ=NW
      RETURN
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = 0
      IF (XX.GE.0.0E0) GO TO 90
      MR = 1
      IF (YY.LT.0.0E0) MR = -1
   90 CONTINUE
      CALL CBUNK(Z, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 200
      NZ = NZ + NW
      RETURN
  100 CONTINUE
      IF (XX.LT.0.0E0) GO TO 180
      RETURN
  180 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  200 CONTINUE
      IF(NW.EQ.(-1)) GO TO 180
      NZ=0
      IERR=5
      RETURN
  210 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
      SUBROUTINE CBINU(Z, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CBINU
C***REFER TO  CBESH,CBESI,CBESJ,CBESK,CAIRY,CBIRY
C
C     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  CASYI,CBUNI,CMLRI,CSERI,CUOIK,CWRSK
C***END PROLOGUE  CBINU
      COMPLEX CW, CY, CZERO, Z
      REAL ALIM, AZ, DFNU, ELIM, FNU, FNUL, RL, TOL
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CY(N), CW(2)
      DATA CZERO / (0.0E0,0.0E0) /
C
      NZ = 0
      AZ = CABS(Z)
      NN = N
      DFNU = FNU + FLOAT(N-1)
      IF (AZ.LE.2.0E0) GO TO 10
      IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
      CALL CSERI(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      INW = IABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + FLOAT(NN-1)
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0E0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30 CONTINUE
      CALL CASYI(Z, FNU, KODE, NN, CY, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0E0) GO TO 70
   50 CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      CALL CUOIK(Z, FNU, KODE, 1, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+FLOAT(NN-1)
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
      CALL CMLRI(Z, FNU, KODE, NN, CY, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
      CALL CUOIK(Z, FNU, KODE, 2, 2, CW, NW, TOL, ELIM, ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CY(I) = CZERO
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL CWRSK(Z, FNU, KODE, NN, CY, NW, CW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
      NUI = INT(FNUL-DFNU) + 1
      NUI = MAX0(NUI,0)
      CALL CBUNI(Z, FNU, KODE, NN, CY, NW, NUI, NLAST, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      SUBROUTINE CBKNU(Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CBKNU
C***REFER TO  CBESI,CBESK,CAIRY,CBESH
C
C     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
C
C***ROUTINES CALLED  CKSCL,CSHCH,GAMLN,I1MACH,R1MACH,CUCHK
C***END PROLOGUE  CBKNU
C
      COMPLEX CCH, CK, COEF, CONE, CRSC, CS, CSCL, CSH, CSR, CSS, CTWO,
     * CZ, CZERO, F, FMU, P, PT, P1, P2, Q, RZ, SMU, ST, S1, S2, Y, Z,
     * ZD, CELM, CY
      REAL AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, CC, DNU,
     * DNU2, ELIM, ETEST, FC, FHS, FK, FKS, FNU, FPI, G1, G2, HPI, PI,
     * P2I, P2M, P2R, RK, RTHPI, R1, S, SPI, TM, TOL, TTH, T1, T2, XX,
     * YY, GAMLN, R1MACH, HELIM, ELM, XD, YD, ALAS, AS
      INTEGER I, IDUM, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N,
     * NZ, I1MACH, NW, J, IC, INUB
      DIMENSION BRY(3), CC(8), CSS(3), CSR(3), Y(N), CY(2)
C
      DATA KMAX / 30 /
      DATA R1 / 2.0E0 /
      DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
C
      DATA PI, RTHPI, SPI ,HPI, FPI, TTH /
     1     3.14159265358979324E0,       1.25331413731550025E0,
     2     1.90985931710274403E0,       1.57079632679489662E0,
     3     1.89769999331517738E0,       6.66666666666666666E-01/
C
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1     5.77215664901532861E-01,    -4.20026350340952355E-02,
     2    -4.21977345555443367E-02,     7.21894324666309954E-03,
     3    -2.15241674114950973E-04,    -2.01348547807882387E-05,
     4     1.13302723198169588E-06,     6.11609510448141582E-09/
C
      XX = REAL(Z)
      YY = AIMAG(Z)
      CAZ = CABS(Z)
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RZ = CTWO/Z
      INU = INT(FNU+0.5E0)
      DNU = FNU - FLOAT(INU)
      IF (ABS(DNU).EQ.0.5E0) GO TO 110
      DNU2 = 0.0E0
      IF (ABS(DNU).GT.TOL) DNU2 = DNU*DNU
      IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR CABS(Z).LE.R1
C-----------------------------------------------------------------------
      FC = 1.0E0
      SMU = CLOG(RZ)
      FMU = SMU*CMPLX(DNU,0.0E0)
      CALL CSHCH(FMU, CSH, CCH)
      IF (DNU.EQ.0.0E0) GO TO 10
      FC = DNU*PI
      FC = FC/SIN(FC)
      SMU = CSH*CMPLX(1.0E0/DNU,0.0E0)
   10 CONTINUE
      A2 = 1.0E0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
      T2 = EXP(-GAMLN(A2,IDUM))
      T1 = 1.0E0/(T2*FC)
      IF (ABS(DNU).GT.0.1E0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
      AK = 1.0E0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (ABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = 0.5E0*(T1+T2)*FC
      G1 = G1*FC
      F = CMPLX(G1,0.0E0)*CCH + SMU*CMPLX(G2,0.0E0)
      PT = CEXP(FMU)
      P = CMPLX(0.5E0/T2,0.0E0)*PT
      Q = CMPLX(0.5E0/T1,0.0E0)/PT
      S1 = F
      S2 = P
      AK = 1.0E0
      A1 = 1.0E0
      CK = CONE
      BK = 1.0E0 - DNU2
      IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
      IF (CAZ.LT.TOL) GO TO 70
      CZ = Z*Z*CMPLX(0.25E0,0.0E0)
      T1 = 0.25E0*CAZ*CAZ
   60 CONTINUE
      F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
      P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
      Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
      RK = 1.0E0/AK
      CK = CK*CZ*CMPLX(RK,0.0)
      S1 = S1 + CK*F
      A1 = A1*T1*RK
      BK = BK + AK + AK + 1.0E0
      AK = AK + 1.0E0
      IF (A1.GT.TOL) GO TO 60
   70 CONTINUE
      Y(1) = S1
      IF (KODED.EQ.1) RETURN
      Y(1) = S1*CEXP(Z)
      RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80 CONTINUE
      IF (CAZ.LT.TOL) GO TO 100
      CZ = Z*Z*CMPLX(0.25E0,0.0E0)
      T1 = 0.25E0*CAZ*CAZ
   90 CONTINUE
      F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
      P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
      Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
      RK = 1.0E0/AK
      CK = CK*CZ*CMPLX(RK,0.0E0)
      S1 = S1 + CK*F
      S2 = S2 + CK*(P-F*CMPLX(AK,0.0E0))
      A1 = A1*T1*RK
      BK = BK + AK + AK + 1.0E0
      AK = AK + 1.0E0
      IF (A1.GT.TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      BK = REAL(SMU)
      A1 = FNU + 1.0E0
      AK = A1*ABS(BK)
      IF (AK.GT.ALIM) KFLAG = 3
      P2 = S2*CSS(KFLAG)
      S2 = P2*RZ
      S1 = S1*CSS(KFLAG)
      IF (KODED.EQ.1) GO TO 210
      F = CEXP(Z)
      S1 = S1*F
      S2 = S2*F
      GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110 CONTINUE
      COEF = CMPLX(RTHPI,0.0E0)/CSQRT(Z)
      KFLAG = 2
      IF (KODED.EQ.2) GO TO 120
      IF (XX.GT.ALIM) GO TO 290
C     BLANK LINE
      A1 = EXP(-XX)*REAL(CSS(KFLAG))
      PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
      COEF = COEF*PT
  120 CONTINUE
      IF (ABS(DNU).EQ.0.5E0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR CABS(Z).GT.R1
C-----------------------------------------------------------------------
      AK = COS(PI*DNU)
      AK = ABS(AK)
      IF (AK.EQ.0.0E0) GO TO 300
      FHS = ABS(0.25E0-DNU2)
      IF (FHS.EQ.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(11))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
      T1 = FLOAT(I1MACH(11)-1)*R1MACH(5)*3.321928094E0
      T1 = AMAX1(T1,12.0E0)
      T1 = AMIN1(T1,60.0E0)
      T2 = TTH*T1 - 6.0E0
      IF (XX.NE.0.0E0) GO TO 130
      T1 = HPI
      GO TO 140
  130 CONTINUE
      T1 = ATAN(YY/XX)
      T1 = ABS(T1)
  140 CONTINUE
      IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
C-----------------------------------------------------------------------
      ETEST = AK/(PI*CAZ*TOL)
      FK = 1.0E0
      IF (ETEST.LT.1.0E0) GO TO 180
      FKS = 2.0E0
      RK = CAZ + CAZ + 2.0E0
      A1 = 0.0E0
      A2 = 1.0E0
      DO 150 I=1,KMAX
        AK = FHS/FKS
        BK = RK/(FK+1.0E0)
        TM = A2
        A2 = BK*A2 - AK*A1
        A1 = TM
        RK = RK + 2.0E0
        FKS = FKS + FK + FK + 2.0E0
        FHS = FHS + FK + FK
        FK = FK + 1.0E0
        TM = ABS(A2)*FK
        IF (ETEST.LT.TM) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*SQRT(T2/CAZ)
      FHS = ABS(0.25E0-DNU2)
      GO TO 180
  170 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
C-----------------------------------------------------------------------
      A2 = SQRT(CAZ)
      AK = FPI*AK/(TOL*SQRT(A2))
      AA = 3.0E0*T1/(1.0E0+CAZ)
      BB = 14.7E0*T1/(28.0E0+CAZ)
      AK = (ALOG(AK)+CAZ*COS(AA)/(1.0E0+0.008E0*CAZ))/COS(BB)
      FK = 0.12125E0*AK*AK/CAZ + 1.5E0
  180 CONTINUE
      K = INT(FK)
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      FK = FLOAT(K)
      FKS = FK*FK
      P1 = CZERO
      P2 = CMPLX(TOL,0.0E0)
      CS = P2
      DO 190 I=1,K
        A1 = FKS - FK
        A2 = (FKS+FK)/(A1+FHS)
        RK = 2.0E0/(FK+1.0E0)
        T1 = (FK+XX)*RK
        T2 = YY*RK
        PT = P2
        P2 = (P2*CMPLX(T1,T2)-P1)*CMPLX(A2,0.0E0)
        P1 = PT
        CS = CS + P2
        FKS = A1 - FK + 1.0E0
        FK = FK - 1.0E0
  190 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
      TM = CABS(CS)
      PT = CMPLX(1.0E0/TM,0.0E0)
      S1 = PT*P2
      CS = CONJG(CS)*PT
      S1 = COEF*S1*CS
      IF (INU.GT.0 .OR. N.GT.1) GO TO 200
      ZD = Z
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  200 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
C-----------------------------------------------------------------------
      TM = CABS(P2)
      PT = CMPLX(1.0E0/TM,0.0E0)
      P1 = PT*P1
      P2 = CONJG(P2)*PT
      PT = P1*P2
      S2 = S1*(CONE+(CMPLX(DNU+0.5E0,0.0E0)-PT)/Z)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210 CONTINUE
      CK = CMPLX(DNU+1.0E0,0.0E0)*RZ
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 220
      IF (N.EQ.1) S1=S2
      ZD = Z
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      IF (IFLAG.EQ.1) GO TO 261
  225 CONTINUE
      P1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        CK = CK + RZ
        IF (KFLAG.GE.3) GO TO 230
        P2 = S2*P1
        P2R = REAL(P2)
        P2I = AIMAG(P2)
        P2R = ABS(P2R)
        P2I = ABS(P2I)
        P2M = AMAX1(P2R,P2I)
        IF (P2M.LE.ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*P1
        S2 = P2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        P1 = CSR(KFLAG)
  230 CONTINUE
      IF (N.EQ.1) S1 = S2
  240 CONTINUE
      Y(1) = S1*CSR(KFLAG)
      IF (N.EQ.1) RETURN
      Y(2) = S2*CSR(KFLAG)
      IF (N.EQ.2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      IF (KK.GT.N) RETURN
      P1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2 = S2
        S2 = CK*S2 + S1
        S1 = P2
        CK = CK + RZ
        P2 = S2*P1
        Y(I) = P2
        IF (KFLAG.GE.3) GO TO 260
        P2R = REAL(P2)
        P2I = AIMAG(P2)
        P2R = ABS(P2R)
        P2I = ABS(P2I)
        P2M = AMAX1(P2R,P2I)
        IF (P2M.LE.ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*P1
        S2 = P2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        P1 = CSR(KFLAG)
  260 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5E0*ELIM
      ELM = EXP(-ELIM)
      CELM = CMPLX(ELM,0.0)
      ASCLE = BRY(1)
      ZD = Z
      XD = XX
      YD = YY
      IC = -1
      J = 2
      DO 262 I=1,INU
        ST = S2
        S2 = CK*S2+S1
        S1 = ST
        CK = CK+RZ
        AS = CABS(S2)
        ALAS = ALOG(AS)
        P2R = -XD+ALAS
        IF(P2R.LT.(-ELIM)) GO TO 263
        P2 = -ZD+CLOG(S2)
        P2R = REAL(P2)
        P2I = AIMAG(P2)
        P2M = EXP(P2R)/TOL
        P1 = CMPLX(P2M,0.0E0)*CMPLX(COS(P2I),SIN(P2I))
        CALL CUCHK(P1,NW,ASCLE,TOL)
        IF(NW.NE.0) GO TO 263
        J=3-J
        CY(J) = P1
        IF(IC.EQ.(I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 262
        XD = XD-ELIM
        S1 = S1*CELM
        S2 = S2*CELM
        ZD = CMPLX(XD,YD)
  262 CONTINUE
      IF(N.EQ.1) S1 = S2
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2 = CY(J)
      J = 3 - J
      S1 = CY(J)
      IF(INUB.LE.INU) GO TO 225
      IF(N.EQ.1) S1 = S2
      GO TO 240
  270 CONTINUE
      Y(1) = S1
      IF (N.EQ.1) GO TO 280
      Y(2) = S2
  280 CONTINUE
      ASCLE = BRY(1)
      CALL CKSCL(ZD, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
      INU = N - NZ
      IF (INU.LE.0) RETURN
      KK = NZ + 1
      S1 = Y(KK)
      Y(KK) = S1*CSR(1)
      IF (INU.EQ.1) RETURN
      KK = NZ + 2
      S2 = Y(KK)
      Y(KK) = S2*CSR(1)
      IF (INU.EQ.2) RETURN
      T2 = FNU + FLOAT(KK-1)
      CK = CMPLX(T2,0.0E0)*RZ
      KFLAG = 1
      GO TO 250
  290 CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY EXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300 CONTINUE
      S1 = COEF
      S2 = COEF
      GO TO 210
  310 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE CBUNI(Z, FNU, KODE, N, Y, NZ, NUI, NLAST, FNUL, TOL,
     * ELIM, ALIM)
C***BEGIN PROLOGUE  CBUNI
C***REFER TO  CBESI,CBESK
C
C     CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***ROUTINES CALLED  CUNI1,CUNI2,R1MACH
C***END PROLOGUE  CBUNI
      COMPLEX CSCL, CSCR, CY, RZ, ST, S1, S2, Y, Z
      REAL ALIM, AX, AY, DFNU, ELIM, FNU, FNUI, FNUL, GNU, TOL, XX, YY,
     * ASCLE, BRY, STR, STI, STM, R1MACH
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION Y(N), CY(2), BRY(3)
      NZ = 0
      XX = REAL(Z)
      YY = AIMAG(Z)
      AX = ABS(XX)*1.7321E0
      AY = ABS(YY)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) GO TO 60
      FNUI = FLOAT(NUI)
      DFNU = FNU + FLOAT(N-1)
      GNU = DFNU + FNUI
      IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNI1(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNI2(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   20 CONTINUE
      IF (NW.LT.0) GO TO 50
      IF (NW.NE.0) GO TO 90
      AY = CABS(CY(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      AX = 1.0E0
      CSCL = CMPLX(AX,0.0E0)
      IF (AY.GT.BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      AX = 1.0E0/TOL
      CSCL = CMPLX(AX,0.0E0)
      GO TO 25
   21 CONTINUE
      IF (AY.LT.BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE = BRY(3)
      AX = TOL
      CSCL = CMPLX(AX,0.0E0)
   25 CONTINUE
      AY = 1.0E0/AX
      CSCR = CMPLX(AY,0.0E0)
      S1 = CY(2)*CSCL
      S2 = CY(1)*CSCL
      RZ = CMPLX(2.0E0,0.0E0)/Z
      DO 30 I=1,NUI
        ST = S2
        S2 = CMPLX(DFNU+FNUI,0.0E0)*RZ*S2 + S1
        S1 = ST
        FNUI = FNUI - 1.0E0
        IF (IFLAG.GE.3) GO TO 30
        ST = S2*CSCR
        STR = REAL(ST)
        STI = AIMAG(ST)
        STR = ABS(STR)
        STI = ABS(STI)
        STM = AMAX1(STR,STI)
        IF (STM.LE.ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1 = S1*CSCR
        S2 = ST
        AX = AX*TOL
        AY = 1.0E0/AX
        CSCL = CMPLX(AX,0.0E0)
        CSCR = CMPLX(AY,0.0E0)
        S1 = S1*CSCL
        S2 = S2*CSCL
   30 CONTINUE
      Y(N) = S2*CSCR
      IF (N.EQ.1) RETURN
      NL = N - 1
      FNUI = FLOAT(NL)
      K = NL
      DO 40 I=1,NL
        ST = S2
        S2 = CMPLX(FNU+FNUI,0.0E0)*RZ*S2 + S1
        S1 = ST
        ST = S2*CSCR
        Y(K) = ST
        FNUI = FNUI - 1.0E0
        K = K - 1
        IF (IFLAG.GE.3) GO TO 40
        STR = REAL(ST)
        STI = AIMAG(ST)
        STR = ABS(STR)
        STI = ABS(STI)
        STM = AMAX1(STR,STI)
        IF (STM.LE.ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1 = S1*CSCR
        S2 = ST
        AX = AX*TOL
        AY = 1.0E0/AX
        CSCL = CMPLX(AX,0.0E0)
        CSCR = CMPLX(AY,0.0E0)
        S1 = S1*CSCL
        S2 = S2*CSCL
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
   60 CONTINUE
      IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNI1(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
      GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNI2(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   80 CONTINUE
      IF (NW.LT.0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END
      SUBROUTINE CBUNK(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CBUNK
C***REFER TO  CBESK,CBESH
C
C     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
C
C***ROUTINES CALLED  CUNK1,CUNK2
C***END PROLOGUE  CBUNK
      COMPLEX Y, Z
      REAL ALIM, AX, AY, ELIM, FNU, TOL, XX, YY
      INTEGER KODE, MR, N, NZ
      DIMENSION Y(N)
      NZ = 0
      XX = REAL(Z)
      YY = AIMAG(Z)
      AX = ABS(XX)*1.7321E0
      AY = ABS(YY)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CKSCL(ZR, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
C***BEGIN PROLOGUE  CKSCL
C***REFER TO  CBKNU,CUNK1,CUNK2
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***ROUTINES CALLED  CUCHK
C***END PROLOGUE  CKSCL
      COMPLEX CK, CS, CY, CZERO, RZ, S1, S2, Y, ZR, ZD, CELM
      REAL AA, ASCLE, ACS, AS, CSI, CSR, ELIM, FN, FNU, TOL, XX, ZRI,
     * ELM, ALAS, HELIM
      INTEGER I, IC, K, KK, N, NN, NW, NZ
      DIMENSION Y(N), CY(2)
      DATA CZERO / (0.0E0,0.0E0) /
C
      NZ = 0
      IC = 0
      XX = REAL(ZR)
      NN = MIN0(2,N)
      DO 10 I=1,NN
        S1 = Y(I)
        CY(I) = S1
        AS = CABS(S1)
        ACS = -XX + ALOG(AS)
        NZ = NZ + 1
        Y(I) = CZERO
        IF (ACS.LT.(-ELIM)) GO TO 10
        CS = -ZR + CLOG(S1)
        CSR = REAL(CS)
        CSI = AIMAG(CS)
        AA = EXP(CSR)/TOL
        CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
        CALL CUCHK(CS, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        Y(I) = CS
        NZ = NZ - 1
        IC = I
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      Y(1) = CZERO
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0E0
      CK = CMPLX(FN,0.0E0)*RZ
      S1 = CY(1)
      S2 = CY(2)
      HELIM = 0.5E0*ELIM
      ELM = EXP(-ELIM)
      CELM = CMPLX(ELM,0.0E0)
      ZRI =AIMAG(ZR)
      ZD = ZR
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
      DO 30 I=3,N
        KK = I
        CS = S2
        S2 = CK*S2 + S1
        S1 = CS
        CK = CK + RZ
        AS = CABS(S2)
        ALAS = ALOG(AS)
        ACS = -XX + ALAS
        NZ = NZ + 1
        Y(I) = CZERO
        IF (ACS.LT.(-ELIM)) GO TO 25
        CS = -ZD + CLOG(S2)
        CSR = REAL(CS)
        CSI = AIMAG(CS)
        AA = EXP(CSR)/TOL
        CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
        CALL CUCHK(CS, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        Y(I) = CS
        NZ = NZ - 1
        IF (IC.EQ.(KK-1)) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        XX = XX-ELIM
        S1 = S1*CELM
        S2 = S2*CELM
        ZD = CMPLX(XX,ZRI)
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 K=1,NZ
        Y(K) = CZERO
   50 CONTINUE
      RETURN
      END
      SUBROUTINE CMLRI(Z, FNU, KODE, N, Y, NZ, TOL)
C***BEGIN PROLOGUE  CMLRI
C***REFER TO  CBESI,CBESK
C
C     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***ROUTINES CALLED  GAMLN,R1MACH
C***END PROLOGUE  CMLRI
      COMPLEX CK, CNORM, CONE, CTWO, CZERO, PT, P1, P2, RZ, SUM, Y, Z
      REAL ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF, FNU, RHO,
     * RHO2, SCLE, TFNF, TOL, TST, X, GAMLN, R1MACH
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N
      DIMENSION Y(N)
      DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
      SCLE = 1.0E+3*R1MACH(1)/TOL
      NZ=0
      AZ = CABS(Z)
      X = REAL(Z)
      IAZ = INT(AZ)
      IFNU = INT(FNU)
      INU = IFNU + N - 1
      AT = FLOAT(IAZ) + 1.0E0
      CK = CMPLX(AT,0.0E0)/Z
      RZ = CTWO/Z
      P1 = CZERO
      P2 = CONE
      ACK = (AT+1.0E0)/AZ
      RHO = ACK + SQRT(ACK*ACK-1.0E0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0E0)*(RHO-1.0E0))
      TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PT = P2
        P2 = P1 - CK*P2
        P1 = PT
        CK = CK + RZ
        AP = CABS(P2)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0E0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
      P1 = CZERO
      P2 = CONE
      AT = FLOAT(INU) + 1.0E0
      CK = CMPLX(AT,0.0E0)/Z
      ACK = AT/AZ
      TST = SQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PT = P2
        P2 = P1 - CK*P2
        P1 = PT
        CK = CK + RZ
        AP = CABS(P2)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = CABS(CK)
        FLAM = ACK + SQRT(ACK*ACK-1.0E0)
        FKAP = AP/CABS(P1)
        RHO = AMIN1(FLAM,FKAP)
        TST = TST*SQRT(RHO/(RHO*RHO-1.0E0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
      K = K + 1
      KK = MAX0(I+IAZ,K+INU)
      FKK = FLOAT(KK)
      P1 = CZERO
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
      P2 = CMPLX(SCLE,0.0E0)
      FNF = FNU - FLOAT(IFNU)
      TFNF = FNF + FNF
      BK = GAMLN(FKK+TFNF+1.0E0,IDUM) - GAMLN(FKK+1.0E0,IDUM)
     *     -GAMLN(TFNF+1.0E0,IDUM)
      BK = EXP(BK)
      SUM = CZERO
      KM = KK - INU
      DO 50 I=1,KM
        PT = P2
        P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
        P1 = PT
        AK = 1.0E0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
        BK = ACK
        FKK = FKK - 1.0E0
   50 CONTINUE
      Y(N) = P2
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PT = P2
        P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
        P1 = PT
        AK = 1.0E0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
        BK = ACK
        FKK = FKK - 1.0E0
        M = N - I + 1
        Y(M) = P2
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PT = P2
        P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
        P1 = PT
        AK = 1.0E0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
        BK = ACK
        FKK = FKK - 1.0E0
   80 CONTINUE
   90 CONTINUE
      PT = Z
      IF (KODE.EQ.2) PT = PT - CMPLX(X,0.0E0)
      P1 = -CMPLX(FNF,0.0E0)*CLOG(RZ) + PT
      AP = GAMLN(1.0E0+FNF,IDUM)
      PT = P1 - CMPLX(AP,0.0E0)
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
      P2 = P2 + SUM
      AP = CABS(P2)
      P1 = CMPLX(1.0E0/AP,0.0E0)
      CK = CEXP(PT)*P1
      PT = CONJG(P2)*P1
      CNORM = CK*PT
      DO 100 I=1,N
        Y(I) = Y(I)*CNORM
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE CRATI(Z, FNU, N, CY, TOL)
C***BEGIN PROLOGUE  CRATI
C***REFER TO  CBESI,CBESK,CBESH
C
C     CRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CRATI
      COMPLEX CDFNU, CONE, CY, CZERO, PT, P1, P2, RZ, T1, Z
      REAL AK, AMAGZ, AP1, AP2, ARG, AZ, DFNU, FDNU, FLAM, FNU, FNUP,
     * RAP1, RHO, TEST, TEST1, TOL
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CY(N)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
      AZ = CABS(Z)
      INU = INT(FNU)
      IDNU = INU + N - 1
      FDNU = FLOAT(IDNU)
      MAGZ = INT(AZ)
      AMAGZ = FLOAT(MAGZ+1)
      FNUP = AMAX1(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      RZ = (CONE+CONE)/Z
      T1 = CMPLX(FNUP,0.0E0)*RZ
      P2 = -T1
      P1 = CONE
      T1 = T1 + RZ
      IF (ID.GT.0) ID = 0
      AP2 = CABS(P2)
      AP1 = CABS(P1)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = SQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0E0/AP1
      P1 = P1*CMPLX(RAP1,0.0E0)
      P2 = P2*CMPLX(RAP1,0.0E0)
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PT = P2
      P2 = P1 - T1*P2
      P1 = PT
      T1 = T1 + RZ
      AP2 = CABS(P2)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = CABS(T1)*0.5E0
      FLAM = AK + SQRT(AK*AK-1.0E0)
      RHO = AMIN1(AP2/AP1,FLAM)
      TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0E0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = FLOAT(KK)
      DFNU = FNU + FLOAT(N-1)
      CDFNU = CMPLX(DFNU,0.0E0)
      T1 = CMPLX(AK,0.0E0)
      P1 = CMPLX(1.0E0/AP2,0.0E0)
      P2 = CZERO
      DO 30 I=1,KK
        PT = P1
        P1 = RZ*(CDFNU+T1)*P1 + P2
        P2 = PT
        T1 = T1 - CONE
   30 CONTINUE
      IF (REAL(P1).NE.0.0E0 .OR. AIMAG(P1).NE.0.0E0) GO TO 40
      P1 = CMPLX(TOL,TOL)
   40 CONTINUE
      CY(N) = P2/P1
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = FLOAT(K)
      T1 = CMPLX(AK,0.0E0)
      CDFNU = CMPLX(FNU,0.0E0)*RZ
      DO 60 I=2,N
        PT = CDFNU + T1*RZ + CY(K+1)
        IF (REAL(PT).NE.0.0E0 .OR. AIMAG(PT).NE.0.0E0) GO TO 50
        PT = CMPLX(TOL,TOL)
   50   CONTINUE
        CY(K) = CONE/PT
        T1 = T1 - CONE
        K = K - 1
   60 CONTINUE
      RETURN
      END
      SUBROUTINE CS1S2(ZR, S1, S2, NZ, ASCLE, ALIM, IUF)
C***BEGIN PROLOGUE  CS1S2
C***REFER TO  CBESK,CAIRY
C
C     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CS1S2
      COMPLEX CZERO, C1, S1, S1D, S2, ZR
      REAL AA, ALIM, ALN, ASCLE, AS1, AS2, XX
      INTEGER IUF, NZ
      DATA CZERO / (0.0E0,0.0E0) /
      NZ = 0
      AS1 = CABS(S1)
      AS2 = CABS(S2)
      AA = REAL(S1)
      ALN = AIMAG(S1)
      IF (AA.EQ.0.0E0 .AND. ALN.EQ.0.0E0) GO TO 10
      IF (AS1.EQ.0.0E0) GO TO 10
      XX = REAL(ZR)
      ALN = -XX - XX + ALOG(AS1)
      S1D = S1
      S1 = CZERO
      AS1 = 0.0E0
      IF (ALN.LT.(-ALIM)) GO TO 10
      C1 = CLOG(S1D) - ZR - ZR
      S1 = CEXP(C1)
      AS1 = CABS(S1)
      IUF = IUF + 1
   10 CONTINUE
      AA = AMAX1(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1 = CZERO
      S2 = CZERO
      NZ = 1
      IUF = 0
      RETURN
      END
      SUBROUTINE CSERI(Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CSERI
C***REFER TO  CBESI,CBESK
C
C     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***ROUTINES CALLED  CUCHK,GAMLN,R1MACH
C***END PROLOGUE  CSERI
      COMPLEX AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ, S1, S2, W,
     * Y, Z
      REAL AA, ACZ, AK, ALIM, ARM, ASCLE, ATOL, AZ, DFNU, ELIM, FNU,
     * FNUP, RAK1, RS, RTR1, S, SS, TOL, X, GAMLN, R1MACH
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NW, NZ
      DIMENSION Y(N), W(2)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      NZ = 0
      AZ = CABS(Z)
      IF (AZ.EQ.0.0E0) GO TO 150
      X = REAL(Z)
      ARM = 1.0E+3*R1MACH(1)
      RTR1 = SQRT(ARM)
      CRSC = CMPLX(1.0E0,0.0E0)
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 140
      HZ = Z*CMPLX(0.5E0,0.0E0)
      CZ = CZERO
      IF (AZ.GT.RTR1) CZ = HZ*HZ
      ACZ = CABS(CZ)
      NN = N
      CK = CLOG(HZ)
   10 CONTINUE
      DFNU = FNU + FLOAT(NN-1)
      FNUP = DFNU + 1.0E0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      AK1 = CK*CMPLX(DFNU,0.0E0)
      AK = GAMLN(FNUP,IDUM)
      AK1 = AK1 - CMPLX(AK,0.0E0)
      IF (KODE.EQ.2) AK1 = AK1 - CMPLX(X,0.0E0)
      RAK1 = REAL(AK1)
      IF (RAK1.GT.(-ELIM)) GO TO 30
   20 CONTINUE
      NZ = NZ + 1
      Y(NN) = CZERO
      IF (ACZ.GT.DFNU) GO TO 170
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 10
   30 CONTINUE
      IF (RAK1.GT.(-ALIM)) GO TO 40
      IFLAG = 1
      SS = 1.0E0/TOL
      CRSC = CMPLX(TOL,0.0E0)
      ASCLE = ARM*SS
   40 CONTINUE
      AK = AIMAG(AK1)
      AA = EXP(RAK1)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEF = CMPLX(AA,0.0E0)*CMPLX(COS(AK),SIN(AK))
      ATOL = TOL*ACZ/FNUP
      IL = MIN0(2,NN)
      DO 80 I=1,IL
        DFNU = FNU + FLOAT(NN-I)
        FNUP = DFNU + 1.0E0
        S1 = CONE
        IF (ACZ.LT.TOL*FNUP) GO TO 60
        AK1 = CONE
        AK = FNUP + 2.0E0
        S = FNUP
        AA = 2.0E0
   50   CONTINUE
        RS = 1.0E0/S
        AK1 = AK1*CZ*CMPLX(RS,0.0E0)
        S1 = S1 + AK1
        S = S + AK
        AK = AK + 2.0E0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 50
   60   CONTINUE
        M = NN - I + 1
        S2 = S1*COEF
        W(I) = S2
        IF (IFLAG.EQ.0) GO TO 70
        CALL CUCHK(S2, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 20
   70   CONTINUE
        Y(M) = S2*CRSC
        IF (I.NE.IL) COEF = COEF*CMPLX(DFNU,0.0E0)/HZ
   80 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = FLOAT(K)
      RZ = (CONE+CONE)/Z
      IF (IFLAG.EQ.1) GO TO 110
      IB = 3
   90 CONTINUE
      DO 100 I=IB,NN
        Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
        AK = AK - 1.0E0
        K = K - 1
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  110 CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3
C-----------------------------------------------------------------------
      S1 = W(1)
      S2 = W(2)
      DO 120 L=3,NN
        CK = S2
        S2 = S1 + CMPLX(AK+FNU,0.0E0)*RZ*S2
        S1 = CK
        CK = S2*CRSC
        Y(K) = CK
        AK = AK - 1.0E0
        K = K - 1
        IF (CABS(CK).GT.ASCLE) GO TO 130
  120 CONTINUE
      RETURN
  130 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 90
  140 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0E0) NZ = NZ - 1
  150 CONTINUE
      Y(1) = CZERO
      IF (FNU.EQ.0.0E0) Y(1) = CONE
      IF (N.EQ.1) RETURN
      DO 160 I=2,N
        Y(I) = CZERO
  160 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
C-----------------------------------------------------------------------
  170 CONTINUE
      NZ = -NZ
      RETURN
      END
      SUBROUTINE CSHCH(Z, CSH, CCH)
C***BEGIN PROLOGUE  CSHCH
C***REFER TO  CBESK,CBESH
C
C     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CSHCH
      COMPLEX CCH, CSH, Z
      REAL CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y, COSH, SINH
      X = REAL(Z)
      Y = AIMAG(Z)
      SH = SINH(X)
      CH = COSH(X)
      SN = SIN(Y)
      CN = COS(Y)
      CSHR = SH*CN
      CSHI = CH*SN
      CSH = CMPLX(CSHR,CSHI)
      CCHR = CH*CN
      CCHI = SH*SN
      CCH = CMPLX(CCHR,CCHI)
      RETURN
      END
      SUBROUTINE CUCHK(Y, NZ, ASCLE, TOL)
C***BEGIN PROLOGUE  CUCHK
C***REFER TO CSERI,CUOIK,CUNK1,CUNK2,CUNI1,CUNI2,CKSCL
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*R1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDER FLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUCHK
C
      COMPLEX Y
      REAL ASCLE, SS, ST, TOL, YR, YI
      INTEGER NZ
      NZ = 0
      YR = REAL(Y)
      YI = AIMAG(Y)
      YR = ABS(YR)
      YI = ABS(YI)
      ST = AMIN1(YR,YI)
      IF (ST.GT.ASCLE) RETURN
      SS = AMAX1(YR,YI)
      ST=ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
      END
      SUBROUTINE CUNHJ(Z, FNU, IPMTR, TOL, PHI, ARG, ZETA1, ZETA2,
     * ASUM, BSUM)
C***BEGIN PROLOGUE  CUNHJ
C***REFER TO  CBESI,CBESK
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         CUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUNHJ
      COMPLEX ARG, ASUM, BSUM, CFNU, CONE, CR, CZERO, DR, P, PHI,
     * PRZTH, PTFN, RFN13, RTZTA, RZTH, SUMA, SUMB, TFN, T2, UP, W, W2,
     * Z, ZA, ZB, ZC, ZETA, ZETA1, ZETA2, ZTH
      REAL ALFA, ANG, AP, AR, ATOL, AW2, AZTH, BETA, BR, BTOL, C, EX1,
     * EX2, FNU, FN13, FN23, GAMA, HPI, PI, PP, RFNU, RFNU2, THPI, TOL,
     * WI, WR, ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR, ASUMR, ASUMI, BSUMR,
     * BSUMI, TEST, TSTR, TSTI, AC
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), P(30), UP(14), CR(14), DR(14)
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000E+00,     1.04166666666666667E-01,
     3     8.35503472222222222E-02,     1.28226574556327160E-01,
     4     2.91849026464140464E-01,     8.81627267443757652E-01,
     5     3.32140828186276754E+00,     1.49957629868625547E+01,
     6     7.89230130115865181E+01,     4.74451538868264323E+02,
     7     3.20749009089066193E+03,     2.40865496408740049E+04,
     8     1.98923119169509794E+05,     1.79190200777534383E+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000E+00,    -1.45833333333333333E-01,
     3    -9.87413194444444444E-02,    -1.43312053915895062E-01,
     4    -3.17227202678413548E-01,    -9.42429147957120249E-01,
     5    -3.51120304082635426E+00,    -1.57272636203680451E+01,
     6    -8.22814390971859444E+01,    -4.92355370523670524E+02,
     7    -3.31621856854797251E+03,    -2.48276742452085896E+04,
     8    -2.04526587315129788E+05,    -1.83844491706820990E+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000E+00,    -2.08333333333333333E-01,
     4     1.25000000000000000E-01,     3.34201388888888889E-01,
     5    -4.01041666666666667E-01,     7.03125000000000000E-02,
     6    -1.02581259645061728E+00,     1.84646267361111111E+00,
     7    -8.91210937500000000E-01,     7.32421875000000000E-02,
     8     4.66958442342624743E+00,    -1.12070026162229938E+01,
     9     8.78912353515625000E+00,    -2.36408691406250000E+00,
     A     1.12152099609375000E-01,    -2.82120725582002449E+01,
     B     8.46362176746007346E+01,    -9.18182415432400174E+01,
     C     4.25349987453884549E+01,    -7.36879435947963170E+00,
     D     2.27108001708984375E-01,     2.12570130039217123E+02,
     E    -7.65252468141181642E+02,     1.05999045252799988E+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541E+02,     2.18190511744211590E+02,
     4    -2.64914304869515555E+01,     5.72501420974731445E-01,
     5    -1.91945766231840700E+03,     8.06172218173730938E+03,
     6    -1.35865500064341374E+04,     1.16553933368645332E+04,
     7    -5.30564697861340311E+03,     1.20090291321635246E+03,
     8    -1.08090919788394656E+02,     1.72772750258445740E+00,
     9     2.02042913309661486E+04,    -9.69805983886375135E+04,
     A     1.92547001232531532E+05,    -2.03400177280415534E+05,
     B     1.22200464983017460E+05,    -4.11926549688975513E+04,
     C     7.10951430248936372E+03,    -4.93915304773088012E+02,
     D     6.07404200127348304E+00,    -2.42919187900551333E+05,
     E     1.31176361466297720E+06,    -2.99801591853810675E+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400E+06,    -2.81356322658653411E+06,
     4     1.26836527332162478E+06,    -3.31645172484563578E+05,
     5     4.52187689813627263E+04,    -2.49983048181120962E+03,
     6     2.43805296995560639E+01,     3.28446985307203782E+06,
     7    -1.97068191184322269E+07,     5.09526024926646422E+07,
     8    -7.41051482115326577E+07,     6.63445122747290267E+07,
     9    -3.75671766607633513E+07,     1.32887671664218183E+07,
     A    -2.78561812808645469E+06,     3.08186404612662398E+05,
     B    -1.38860897537170405E+04,     1.10017140269246738E+02,
     C    -4.93292536645099620E+07,     3.25573074185765749E+08,
     D    -9.39462359681578403E+08,     1.55359689957058006E+09,
     E    -1.62108055210833708E+09,     1.10684281682301447E+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309E+08,     1.42062907797533095E+08,
     4    -2.44740627257387285E+07,     2.24376817792244943E+06,
     5    -8.40054336030240853E+04,     5.51335896122020586E+02,
     6     8.14789096118312115E+08,    -5.86648149205184723E+09,
     7     1.86882075092958249E+10,    -3.46320433881587779E+10,
     8     4.12801855797539740E+10,    -3.30265997498007231E+10,
     9     1.79542137311556001E+10,    -6.56329379261928433E+09,
     A     1.55927986487925751E+09,    -2.25105661889415278E+08,
     B     1.73951075539781645E+07,    -5.49842327572288687E+05,
     C     3.03809051092238427E+03,    -1.46792612476956167E+10,
     D     1.14498237732025810E+11,    -3.99096175224466498E+11,
     E     8.19218669548577329E+11,    -1.09837515608122331E+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209E+12,    -6.45364869245376503E+11,
     3     2.87900649906150589E+11,    -8.78670721780232657E+10,
     4     1.76347306068349694E+10,    -2.16716498322379509E+09,
     5     1.43157876718888981E+08,    -3.87183344257261262E+06,
     6     1.82577554742931747E+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444E-03,    -9.22077922077922078E-04,
     5    -8.84892884892884893E-05,     1.65927687832449737E-04,
     6     2.46691372741792910E-04,     2.65995589346254780E-04,
     7     2.61824297061500945E-04,     2.48730437344655609E-04,
     8     2.32721040083232098E-04,     2.16362485712365082E-04,
     9     2.00738858762752355E-04,     1.86267636637545172E-04,
     A     1.73060775917876493E-04,     1.61091705929015752E-04,
     B     1.50274774160908134E-04,     1.40503497391269794E-04,
     C     1.31668816545922806E-04,     1.23667445598253261E-04,
     D     1.16405271474737902E-04,     1.09798298372713369E-04,
     E     1.03772410422992823E-04,     9.82626078369363448E-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256E-05,     8.85710852478711718E-05,
     5     8.42963105715700223E-05,     8.03497548407791151E-05,
     6     7.66981345359207388E-05,     7.33122157481777809E-05,
     7     7.01662625163141333E-05,     6.72375633790160292E-05,
     8     6.93735541354588974E-04,     2.32241745182921654E-04,
     9    -1.41986273556691197E-05,    -1.16444931672048640E-04,
     A    -1.50803558053048762E-04,    -1.55121924918096223E-04,
     B    -1.46809756646465549E-04,    -1.33815503867491367E-04,
     C    -1.19744975684254051E-04,    -1.06184319207974020E-04,
     D    -9.37699549891194492E-05,    -8.26923045588193274E-05,
     E    -7.29374348155221211E-05,    -6.44042357721016283E-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048E-05,    -5.04731044303561628E-05,
     5    -4.48134868008882786E-05,    -3.98688727717598864E-05,
     6    -3.55400532972042498E-05,    -3.17414256609022480E-05,
     7    -2.83996793904174811E-05,    -2.54522720634870566E-05,
     8    -2.28459297164724555E-05,    -2.05352753106480604E-05,
     9    -1.84816217627666085E-05,    -1.66519330021393806E-05,
     A    -1.50179412980119482E-05,    -1.35554031379040526E-05,
     B    -1.22434746473858131E-05,    -1.10641884811308169E-05,
     C    -3.54211971457743841E-04,    -1.56161263945159416E-04,
     D     3.04465503594936410E-05,     1.30198655773242693E-04,
     E     1.67471106699712269E-04,     1.70222587683592569E-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704E-04,     1.36339170977445120E-04,
     5     1.14886692029825128E-04,     9.45869093034688111E-05,
     6     7.64498419250898258E-05,     6.07570334965197354E-05,
     7     4.74394299290508799E-05,     3.62757512005344297E-05,
     8     2.69939714979224901E-05,     1.93210938247939253E-05,
     9     1.30056674793963203E-05,     7.82620866744496661E-06,
     A     3.59257485819351583E-06,     1.44040049814251817E-07,
     B    -2.65396769697939116E-06,    -4.91346867098485910E-06,
     C    -6.72739296091248287E-06,    -8.17269379678657923E-06,
     D    -9.31304715093561232E-06,    -1.02011418798016441E-05,
     E    -1.08805962510592880E-05,    -1.13875481509603555E-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414E-05,    -1.19987364870944141E-05,
     5     3.78194199201772914E-04,     2.02471952761816167E-04,
     6    -6.37938506318862408E-05,    -2.38598230603005903E-04,
     7    -3.10916256027361568E-04,    -3.13680115247576316E-04,
     8    -2.78950273791323387E-04,    -2.28564082619141374E-04,
     9    -1.75245280340846749E-04,    -1.25544063060690348E-04,
     A    -8.22982872820208365E-05,    -4.62860730588116458E-05,
     B    -1.72334302366962267E-05,     5.60690482304602267E-06,
     C     2.31395443148286800E-05,     3.62642745856793957E-05,
     D     4.58006124490188752E-05,     5.24595294959114050E-05,
     E     5.68396208545815266E-05,     5.94349820393104052E-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742E-05,     6.08023907788436497E-05,
     5     6.01577894539460388E-05,     5.89199657344698500E-05,
     6     5.72515823777593053E-05,     5.52804375585852577E-05,
     7     5.31063773802880170E-05,     5.08069302012325706E-05,
     8     4.84418647620094842E-05,     4.60568581607475370E-05,
     9    -6.91141397288294174E-04,    -4.29976633058871912E-04,
     A     1.83067735980039018E-04,     6.60088147542014144E-04,
     B     8.75964969951185931E-04,     8.77335235958235514E-04,
     C     7.49369585378990637E-04,     5.63832329756980918E-04,
     D     3.68059319971443156E-04,     1.88464535514455599E-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149E-05,    -8.28520220232137023E-05,
     5    -1.72751952869172998E-04,    -2.36314873605872983E-04,
     6    -2.77966150694906658E-04,    -3.02079514155456919E-04,
     7    -3.12594712643820127E-04,    -3.12872558758067163E-04,
     8    -3.05678038466324377E-04,    -2.93226470614557331E-04,
     9    -2.77255655582934777E-04,    -2.59103928467031709E-04,
     A    -2.39784014396480342E-04,    -2.20048260045422848E-04,
     B    -2.00443911094971498E-04,    -1.81358692210970687E-04,
     C    -1.63057674478657464E-04,    -1.45712672175205844E-04,
     D    -1.29425421983924587E-04,    -1.14245691942445952E-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885E-03,     1.35592576302022234E-03,
     5    -7.17858090421302995E-04,    -2.58084802575270346E-03,
     6    -3.49271130826168475E-03,    -3.46986299340960628E-03,
     7    -2.82285233351310182E-03,    -1.88103076404891354E-03,
     8    -8.89531718383947600E-04,     3.87912102631035228E-06,
     9     7.28688540119691412E-04,     1.26566373053457758E-03,
     A     1.62518158372674427E-03,     1.83203153216373172E-03,
     B     1.91588388990527909E-03,     1.90588846755546138E-03,
     C     1.82798982421825727E-03,     1.70389506421121530E-03,
     D     1.55097127171097686E-03,     1.38261421852276159E-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774E-03,     1.03676532638344962E-03,
     3     8.71437918068619115E-04,     7.16080155297701002E-04,
     4     5.72637002558129372E-04,     4.42089819465802277E-04,
     5     3.24724948503090564E-04,     2.20342042730246599E-04,
     6     1.28412898401353882E-04,     4.82005924552095464E-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309E-02,     5.59964911064388073E-03,
     5     2.88501402231132779E-03,     1.80096606761053941E-03,
     6     1.24753110589199202E-03,     9.22878876572938311E-04,
     7     7.14430421727287357E-04,     5.71787281789704872E-04,
     8     4.69431007606481533E-04,     3.93232835462916638E-04,
     9     3.34818889318297664E-04,     2.88952148495751517E-04,
     A     2.52211615549573284E-04,     2.22280580798883327E-04,
     B     1.97541838033062524E-04,     1.76836855019718004E-04,
     C     1.59316899661821081E-04,     1.44347930197333986E-04,
     D     1.31448068119965379E-04,     1.20245444949302884E-04,
     E     1.10449144504599392E-04,     1.01828770740567258E-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509E-05,     8.74130545753834437E-05,
     5     8.13466262162801467E-05,     7.59002269646219339E-05,
     6     7.09906300634153481E-05,     6.65482874842468183E-05,
     7     6.25146958969275078E-05,     5.88403394426251749E-05,
     8    -1.49282953213429172E-03,    -8.78204709546389328E-04,
     9    -5.02916549572034614E-04,    -2.94822138512746025E-04,
     A    -1.75463996970782828E-04,    -1.04008550460816434E-04,
     B    -5.96141953046457895E-05,    -3.12038929076098340E-05,
     C    -1.26089735980230047E-05,    -2.42892608575730389E-07,
     D     8.05996165414273571E-06,     1.36507009262147391E-05,
     E     1.73964125472926261E-05,     1.98672978842133780E-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639E-05,     2.23954659232456514E-05,
     5     2.28967783814712629E-05,     2.30785389811177817E-05,
     6     2.30321976080909144E-05,     2.28236073720348722E-05,
     7     2.25005881105292418E-05,     2.20981015361991429E-05,
     8     2.16418427448103905E-05,     2.11507649256220843E-05,
     9     2.06388749782170737E-05,     2.01165241997081666E-05,
     A     1.95913450141179244E-05,     1.90689367910436740E-05,
     B     1.85533719641636667E-05,     1.80475722259674218E-05,
     C     5.52213076721292790E-04,     4.47932581552384646E-04,
     D     2.79520653992020589E-04,     1.52468156198446602E-04,
     E     6.93271105657043598E-05,     1.76258683069991397E-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136E-05,    -3.17972413350427135E-05,
     5    -4.18861861696693365E-05,    -4.69004889379141029E-05,
     6    -4.87665447413787352E-05,    -4.87010031186735069E-05,
     7    -4.74755620890086638E-05,    -4.55813058138628452E-05,
     8    -4.33309644511266036E-05,    -4.09230193157750364E-05,
     9    -3.84822638603221274E-05,    -3.60857167535410501E-05,
     A    -3.37793306123367417E-05,    -3.15888560772109621E-05,
     B    -2.95269561750807315E-05,    -2.75978914828335759E-05,
     C    -2.58006174666883713E-05,    -2.41308356761280200E-05,
     D    -2.25823509518346033E-05,    -2.11479656768912971E-05,
     E    -1.98200638885294927E-05,    -1.85909870801065077E-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224E-05,    -1.63997823854497997E-05,
     5    -4.74617796559959808E-04,    -4.77864567147321487E-04,
     6    -3.20390228067037603E-04,    -1.61105016119962282E-04,
     7    -4.25778101285435204E-05,     3.44571294294967503E-05,
     8     7.97092684075674924E-05,     1.03138236708272200E-04,
     9     1.12466775262204158E-04,     1.13103642108481389E-04,
     A     1.08651634848774268E-04,     1.01437951597661973E-04,
     B     9.29298396593363896E-05,     8.40293133016089978E-05,
     C     7.52727991349134062E-05,     6.69632521975730872E-05,
     D     5.92564547323194704E-05,     5.22169308826975567E-05,
     E     4.58539485165360646E-05,     4.01445513891486808E-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081E-05,     3.05157995034346659E-05,
     5     2.64956119950516039E-05,     2.29363633690998152E-05,
     6     1.97893056664021636E-05,     1.70091984636412623E-05,
     7     1.45547428261524004E-05,     1.23886640995878413E-05,
     8     1.04775876076583236E-05,     8.79179954978479373E-06,
     9     7.36465810572578444E-04,     8.72790805146193976E-04,
     A     6.22614862573135066E-04,     2.85998154194304147E-04,
     B     3.84737672879366102E-06,    -1.87906003636971558E-04,
     C    -2.97603646594554535E-04,    -3.45998126832656348E-04,
     D    -3.53382470916037712E-04,    -3.35715635775048757E-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809E-04,    -2.66722723047612821E-04,
     5    -2.27654214122819527E-04,    -1.89922611854562356E-04,
     6    -1.55058918599093870E-04,    -1.23778240761873630E-04,
     7    -9.62926147717644187E-05,    -7.25178327714425337E-05,
     8    -5.22070028895633801E-05,    -3.50347750511900522E-05,
     9    -2.06489761035551757E-05,    -8.70106096849767054E-06,
     A     1.13698686675100290E-06,     9.16426474122778849E-06,
     B     1.56477785428872620E-05,     2.08223629482466847E-05,
     C     2.48923381004595156E-05,     2.80340509574146325E-05,
     D     3.03987774629861915E-05,     3.21156731406700616E-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708E-03,    -2.43402962938042533E-03,
     5    -1.83422663549856802E-03,    -7.62204596354009765E-04,
     6     2.39079475256927218E-04,     9.49266117176881141E-04,
     7     1.34467449701540359E-03,     1.48457495259449178E-03,
     8     1.44732339830617591E-03,     1.30268261285657186E-03,
     9     1.10351597375642682E-03,     8.86047440419791759E-04,
     A     6.73073208165665473E-04,     4.77603872856582378E-04,
     B     3.05991926358789362E-04,     1.60315694594721630E-04,
     C     4.00749555270613286E-05,    -5.66607461635251611E-05,
     D    -1.32506186772982638E-04,    -1.90296187989614057E-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408E-04,    -2.62628811464668841E-04,
     5    -2.82050469867598672E-04,    -2.93081563192861167E-04,
     6    -2.97435962176316616E-04,    -2.96557334239348078E-04,
     7    -2.91647363312090861E-04,    -2.83696203837734166E-04,
     8    -2.73512317095673346E-04,    -2.61750155806768580E-04,
     9     6.38585891212050914E-03,     9.62374215806377941E-03,
     A     7.61878061207001043E-03,     2.83219055545628054E-03,
     B    -2.09841352012720090E-03,    -5.73826764216626498E-03,
     C    -7.70804244495414620E-03,    -8.21011692264844401E-03,
     D    -7.65824520346905413E-03,    -6.47209729391045177E-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473E-03,    -3.45612289713133280E-03,
     5    -2.01785580014170775E-03,    -7.59430686781961401E-04,
     6     2.84173631523859138E-04,     1.10891667586337403E-03,
     7     1.72901493872728771E-03,     2.16812590802684701E-03,
     8     2.45357710494539735E-03,     2.61281821058334862E-03,
     9     2.67141039656276912E-03,     2.65203073395980430E-03,
     A     2.57411652877287315E-03,     2.45389126236094427E-03,
     B     2.30460058071795494E-03,     2.13684837686712662E-03,
     C     1.95896528478870911E-03,     1.77737008679454412E-03,
     D     1.59690280765839059E-03,     1.42111975664438546E-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582E-01,     2.51984209978974633E-01,
     5     1.54790300415655846E-01,     1.10713062416159013E-01,
     6     8.57309395527394825E-02,     6.97161316958684292E-02,
     7     5.86085671893713576E-02,     5.04698873536310685E-02,
     8     4.42600580689154809E-02,     3.93720661543509966E-02,
     9     3.54283195924455368E-02,     3.21818857502098231E-02,
     A     2.94646240791157679E-02,     2.71581677112934479E-02,
     B     2.51768272973861779E-02,     2.34570755306078891E-02,
     C     2.19508390134907203E-02,     2.06210828235646240E-02,
     D     1.94388240897880846E-02,     1.83810633800683158E-02,
     E     1.74293213231963172E-02,     1.65685837786612353E-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445E-02,     1.50729501494095594E-02,
     3     1.44193250839954639E-02,     1.38184805735341786E-02,
     4     1.32643378994276568E-02,     1.27517121970498651E-02,
     5     1.22761545318762767E-02,     1.18338262398482403E-02/
      DATA EX1, EX2, HPI, PI, THPI /
     1     3.33333333333333333E-01,     6.66666666666666667E-01,
     2     1.57079632679489662E+00,     3.14159265358979324E+00,
     3     4.71238898038468986E+00/
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      RFNU = 1.0E0/FNU
C     ZB = Z*CMPLX(RFNU,0.0E0)
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TSTR = REAL(Z)
      TSTI = AIMAG(Z)
      TEST = R1MACH(1)*1.0E+3
      AC = FNU*TEST
      IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
      AC = 2.0E0*ABS(ALOG(TEST))+FNU
      ZETA1 = CMPLX(AC,0.0E0)
      ZETA2 = CMPLX(FNU,0.0E0)
      PHI=CONE
      ARG=CONE
      RETURN
   15 CONTINUE
      ZB = Z*CMPLX(RFNU,0.0E0)
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = CMPLX(1.0E0/FN13,0.0E0)
      W2 = CONE - ZB*ZB
      AW2 = CABS(W2)
      IF (AW2.GT.0.25E0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(W2).LE.0.25E0
C-----------------------------------------------------------------------
      K = 1
      P(1) = CONE
      SUMA = CMPLX(GAMA(1),0.0E0)
      AP(1) = 1.0E0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        P(K) = P(K-1)*W2
        SUMA = SUMA + P(K)*CMPLX(GAMA(K),0.0E0)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETA = W2*SUMA
      ARG = ZETA*CMPLX(FN23,0.0E0)
      ZA = CSQRT(SUMA)
      ZETA2 = CSQRT(W2)*CMPLX(FNU,0.0E0)
      ZETA1 = ZETA2*(CONE+ZETA*ZA*CMPLX(EX2,0.0E0))
      ZA = ZA + ZA
      PHI = CSQRT(ZA)*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMB = CZERO
      DO 30 K=1,KMAX
        SUMB = SUMB + P(K)*CMPLX(BETA(K),0.0E0)
   30 CONTINUE
      ASUM = CZERO
      BSUM = SUMB
      L1 = 0
      L2 = 30
      BTOL = TOL*CABS(BSUM)
      ATOL = TOL
      PP = 1.0E0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMA = CZERO
        DO 40 K=1,KMAX
          M = L1 + K
          SUMA = SUMA + P(K)*CMPLX(ALFA(M),0.0E0)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUM = ASUM + SUMA*CMPLX(PP,0.0E0)
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMB = CZERO
        DO 70 K=1,KMAX
          M = L2 + K
          SUMB = SUMB + P(K)*CMPLX(BETA(M),0.0E0)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUM = BSUM + SUMB*CMPLX(PP,0.0E0)
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUM = ASUM + CONE
      PP = RFNU*REAL(RFN13)
      BSUM = BSUM*CMPLX(PP,0.0E0)
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     CABS(W2).GT.0.25E0
C-----------------------------------------------------------------------
  130 CONTINUE
      W = CSQRT(W2)
      WR = REAL(W)
      WI = AIMAG(W)
      IF (WR.LT.0.0E0) WR = 0.0E0
      IF (WI.LT.0.0E0) WI = 0.0E0
      W = CMPLX(WR,WI)
      ZA = (CONE+W)/ZB
      ZC = CLOG(ZA)
      ZCR = REAL(ZC)
      ZCI = AIMAG(ZC)
      IF (ZCI.LT.0.0E0) ZCI = 0.0E0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0E0) ZCR = 0.0E0
      ZC = CMPLX(ZCR,ZCI)
      ZTH = (ZC-W)*CMPLX(1.5E0,0.0E0)
      CFNU = CMPLX(FNU,0.0E0)
      ZETA1 = ZC*CFNU
      ZETA2 = W*CFNU
      AZTH = CABS(ZTH)
      ZTHR = REAL(ZTH)
      ZTHI = AIMAG(ZTH)
      ANG = THPI
      IF (ZTHR.GE.0.0E0 .AND. ZTHI.LT.0.0E0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0E0) GO TO 140
      ANG = ATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0E0) ANG = ANG + PI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*COS(ANG)
      ZETAI = PP*SIN(ANG)
      IF (ZETAI.LT.0.0E0) ZETAI = 0.0E0
      ZETA = CMPLX(ZETAR,ZETAI)
      ARG = ZETA*CMPLX(FN23,0.0E0)
      RTZTA = ZTH/ZETA
      ZA = RTZTA/W
      PHI = CSQRT(ZA+ZA)*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      TFN = CMPLX(RFNU,0.0E0)/W
      RZTH = CMPLX(RFNU,0.0E0)/ZTH
      ZC = RZTH*CMPLX(AR(2),0.0E0)
      T2 = CONE/W2
      UP(2) = (T2*CMPLX(C(2),0.0E0)+CMPLX(C(3),0.0E0))*TFN
      BSUM = UP(2) + ZC
      ASUM = CZERO
      IF (RFNU.LT.TOL) GO TO 220
      PRZTH = RZTH
      PTFN = TFN
      UP(1) = CONE
      PP = 1.0E0
      BSUMR = REAL(BSUM)
      BSUMI = AIMAG(BSUM)
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZA = CMPLX(C(L),0.0E0)
          DO 150 J=2,KP1
            L = L + 1
            ZA = ZA*T2 + CMPLX(C(L),0.0E0)
  150     CONTINUE
          PTFN = PTFN*TFN
          UP(KP1) = PTFN*ZA
          CR(KS) = PRZTH*CMPLX(BR(KS+1),0.0E0)
          PRZTH = PRZTH*RZTH
          DR(KS) = PRZTH*CMPLX(AR(KS+2),0.0E0)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMA = UP(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMA = SUMA + CR(JR)*UP(JU)
  170   CONTINUE
        ASUM = ASUM + SUMA
        ASUMR = REAL(ASUM)
        ASUMI = AIMAG(ASUM)
        TEST = ABS(ASUMR) + ABS(ASUMI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMB = UP(LR+2) + UP(LRP1)*ZC
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMB = SUMB + DR(JR)*UP(JU)
  190   CONTINUE
        BSUM = BSUM + SUMB
        BSUMR = REAL(BSUM)
        BSUMI = AIMAG(BSUM)
        TEST = ABS(BSUMR) + ABS(BSUMI)
        IF (PP.LT.BTOL .AND. TEST.LT.TOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUM = ASUM + CONE
      BSUM = -BSUM*RFN13/RTZTA
      GO TO 120
      END
      SUBROUTINE CUNI1(Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CUNI1
C***REFER TO  CBESI,CBESK
C
C     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  CUCHK,CUNIK,CUOIK,R1MACH
C***END PROLOGUE  CUNI1
      COMPLEX CFN, CONE, CRSC, CSCL, CSR, CSS, CWRK, CZERO, C1, C2,
     * PHI, RZ, SUM, S1, S2, Y, Z, ZETA1, ZETA2, CY
      REAL ALIM, APHI, ASCLE, BRY, C2I, C2M, C2R, ELIM, FN, FNU, FNUL,
     * RS1, TOL, YY, R1MACH
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), Y(N), CWRK(16), CSS(3), CSR(3), CY(2)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = AMAX1(FNU,1.0E0)
      INIT = 0
      CALL CUNIK(Z, FN, 1, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
      IF (KODE.EQ.1) GO TO 10
      CFN = CMPLX(FN,0.0E0)
      S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
      GO TO 20
   10 CONTINUE
      S1 = -ZETA1 + ZETA2
   20 CONTINUE
      RS1 = REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 130
   30 CONTINUE
      NN = MIN0(2,ND)
      DO 80 I=1,NN
        FN = FNU + FLOAT(ND-I)
        INIT = 0
        CALL CUNIK(Z, FN, 1, 0, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
        IF (KODE.EQ.1) GO TO 40
        CFN = CMPLX(FN,0.0E0)
        YY = AIMAG(Z)
        S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + CMPLX(0.0E0,YY)
        GO TO 50
   40   CONTINUE
        S1 = -ZETA1 + ZETA2
   50   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHI)
        RS1 = RS1 + ALOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 60
        IF (I.EQ.1) IFLAG = 3
   60   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF CABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
        S2 = PHI*SUM
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 70
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 110
   70   CONTINUE
        M = ND - I + 1
        CY(I) = S2
        Y(M) = S2*CSR(IFLAG)
   80 CONTINUE
      IF (ND.LE.2) GO TO 100
      RZ = CMPLX(2.0E0,0.0E0)/Z
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = FLOAT(K)
      DO 90 I=3,ND
        C2 = S2
        S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
        S1 = C2
        C2 = S2*C1
        Y(K) = C2
        K = K - 1
        FN = FN - 1.0E0
        IF (IFLAG.GE.3) GO TO 90
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        C1 = CSR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 120
      Y(ND) = CZERO
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 100
      CALL CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 100
      FN = FNU + FLOAT(ND-1)
      IF (FN.GE.FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 120
      NZ = N
      DO 140 I=1,N
        Y(I) = CZERO
  140 CONTINUE
      RETURN
      END
      SUBROUTINE CUNI2(Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM,
     * ALIM)
C***BEGIN PROLOGUE  CUNI2
C***REFER TO  CBESI,CBESK
C
C     CUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***ROUTINES CALLED  CAIRY,CUCHK,CUNHJ,CUOIK,R1MACH
C***END PROLOGUE  CUNI2
      COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CID, CIP, CONE, CRSC, CSCL,
     * CSR, CSS, CY, CZERO, C1, C2, DAI, PHI, RZ, S1, S2, Y, Z, ZB,
     * ZETA1, ZETA2, ZN, ZAR
      REAL AARG, AIC, ALIM, ANG, APHI, ASCLE, AY, BRY, CAR, C2I, C2M,
     * C2R, ELIM, FN, FNU, FNUL, HPI, RS1, SAR, TOL, YY, R1MACH
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     * NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), Y(N), CIP(4), CSS(3), CSR(3), CY(2)
      DATA CZERO,CONE,CI/(0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0)/
      DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     1 (1.0E0,0.0E0), (0.0E0,1.0E0), (-1.0E0,0.0E0), (0.0E0,-1.0E0)/
      DATA HPI, AIC  /
     1      1.57079632679489662E+00,     1.265512123484645396E+00/
C
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      YY = AIMAG(Z)
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
      ZN = -Z*CI
      ZB = Z
      CID = -CI
      INU = INT(FNU)
      ANG = HPI*(FNU-FLOAT(INU))
      CAR = COS(ANG)
      SAR = SIN(ANG)
      C2 = CMPLX(CAR,SAR)
      ZAR = C2
      IN = INU + N - 1
      IN = MOD(IN,4)
      C2 = C2*CIP(IN+1)
      IF (YY.GT.0.0E0) GO TO 10
      ZN = CONJG(-ZN)
      ZB = CONJG(ZB)
      CID = -CID
      C2 = CONJG(C2)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = AMAX1(FNU,1.0E0)
      CALL CUNHJ(ZN, FN, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
      IF (KODE.EQ.1) GO TO 20
      CFN = CMPLX(FNU,0.0E0)
      S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2))
      GO TO 30
   20 CONTINUE
      S1 = -ZETA1 + ZETA2
   30 CONTINUE
      RS1 = REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 150
   40 CONTINUE
      NN = MIN0(2,ND)
      DO 90 I=1,NN
        FN = FNU + FLOAT(ND-I)
        CALL CUNHJ(ZN, FN, 0, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
        IF (KODE.EQ.1) GO TO 50
        CFN = CMPLX(FN,0.0E0)
        AY = ABS(YY)
        S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2)) + CMPLX(0.0E0,AY)
        GO TO 60
   50   CONTINUE
        S1 = -ZETA1 + ZETA2
   60   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        APHI = CABS(PHI)
        AARG = CABS(ARG)
        RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 70
        IF (I.EQ.1) IFLAG = 3
   70   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        CALL CAIRY(ARG, 0, 2, AI, NAI, IDUM)
        CALL CAIRY(ARG, 1, 2, DAI, NDAI, IDUM)
        S2 = PHI*(AI*ASUM+DAI*BSUM)
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 80
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 120
   80   CONTINUE
        IF (YY.LE.0.0E0) S2 = CONJG(S2)
        J = ND - I + 1
        S2 = S2*C2
        CY(I) = S2
        Y(J) = S2*CSR(IFLAG)
        C2 = C2*CID
   90 CONTINUE
      IF (ND.LE.2) GO TO 110
      RZ = CMPLX(2.0E0,0.0E0)/Z
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = FLOAT(K)
      DO 100 I=3,ND
        C2 = S2
        S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
        S1 = C2
        C2 = S2*C1
        Y(K) = C2
        K = K - 1
        FN = FN - 1.0E0
        IF (IFLAG.GE.3) GO TO 100
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        C1 = CSR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
      Y(ND) = CZERO
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 110
      CALL CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 110
      FN = FNU + FLOAT(ND-1)
      IF (FN.LT.FNUL) GO TO 130
C      FN = AIMAG(CID)
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1 = CIP(K)
C      IF (FN.LT.0.0E0) S1 = CONJG(S1)
C      C2 = C2*S1
      IN = INU + ND - 1
      IN = MOD(IN,4) + 1
      C2 = ZAR*CIP(IN)
      IF (YY.LE.0.0E0)C2=CONJG(C2)
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 140
      NZ = N
      DO 160 I=1,N
        Y(I) = CZERO
  160 CONTINUE
      RETURN
      END
      SUBROUTINE CUNIK(ZR, FNU, IKFLG, IPMTR, TOL, INIT, PHI, ZETA1,
     * ZETA2, SUM, CWRK)
C***BEGIN PROLOGUE  CUNIK
C***REFER TO  CBESI,CBESK
C
C        CUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CUNIK
      COMPLEX CFN, CON, CONE, CRFN, CWRK, CZERO, PHI, S, SR, SUM, T,
     * T2, ZETA1, ZETA2, ZN, ZR
      REAL AC, C, FNU, RFN, TEST, TOL, TSTR, TSTI
      INTEGER I, IKFLG, INIT, IPMTR, J, K, L
      DIMENSION C(120), CWRK(16), CON(2)
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
      DATA CON(1), CON(2)  /
     1(3.98942280401432678E-01,0.0E0),(1.25331413731550025E+00,0.0E0)/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000E+00,    -2.08333333333333333E-01,
     4     1.25000000000000000E-01,     3.34201388888888889E-01,
     5    -4.01041666666666667E-01,     7.03125000000000000E-02,
     6    -1.02581259645061728E+00,     1.84646267361111111E+00,
     7    -8.91210937500000000E-01,     7.32421875000000000E-02,
     8     4.66958442342624743E+00,    -1.12070026162229938E+01,
     9     8.78912353515625000E+00,    -2.36408691406250000E+00,
     A     1.12152099609375000E-01,    -2.82120725582002449E+01,
     B     8.46362176746007346E+01,    -9.18182415432400174E+01,
     C     4.25349987453884549E+01,    -7.36879435947963170E+00,
     D     2.27108001708984375E-01,     2.12570130039217123E+02,
     E    -7.65252468141181642E+02,     1.05999045252799988E+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541E+02,     2.18190511744211590E+02,
     4    -2.64914304869515555E+01,     5.72501420974731445E-01,
     5    -1.91945766231840700E+03,     8.06172218173730938E+03,
     6    -1.35865500064341374E+04,     1.16553933368645332E+04,
     7    -5.30564697861340311E+03,     1.20090291321635246E+03,
     8    -1.08090919788394656E+02,     1.72772750258445740E+00,
     9     2.02042913309661486E+04,    -9.69805983886375135E+04,
     A     1.92547001232531532E+05,    -2.03400177280415534E+05,
     B     1.22200464983017460E+05,    -4.11926549688975513E+04,
     C     7.10951430248936372E+03,    -4.93915304773088012E+02,
     D     6.07404200127348304E+00,    -2.42919187900551333E+05,
     E     1.31176361466297720E+06,    -2.99801591853810675E+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400E+06,    -2.81356322658653411E+06,
     4     1.26836527332162478E+06,    -3.31645172484563578E+05,
     5     4.52187689813627263E+04,    -2.49983048181120962E+03,
     6     2.43805296995560639E+01,     3.28446985307203782E+06,
     7    -1.97068191184322269E+07,     5.09526024926646422E+07,
     8    -7.41051482115326577E+07,     6.63445122747290267E+07,
     9    -3.75671766607633513E+07,     1.32887671664218183E+07,
     A    -2.78561812808645469E+06,     3.08186404612662398E+05,
     B    -1.38860897537170405E+04,     1.10017140269246738E+02,
     C    -4.93292536645099620E+07,     3.25573074185765749E+08,
     D    -9.39462359681578403E+08,     1.55359689957058006E+09,
     E    -1.62108055210833708E+09,     1.10684281682301447E+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309E+08,     1.42062907797533095E+08,
     4    -2.44740627257387285E+07,     2.24376817792244943E+06,
     5    -8.40054336030240853E+04,     5.51335896122020586E+02,
     6     8.14789096118312115E+08,    -5.86648149205184723E+09,
     7     1.86882075092958249E+10,    -3.46320433881587779E+10,
     8     4.12801855797539740E+10,    -3.30265997498007231E+10,
     9     1.79542137311556001E+10,    -6.56329379261928433E+09,
     A     1.55927986487925751E+09,    -2.25105661889415278E+08,
     B     1.73951075539781645E+07,    -5.49842327572288687E+05,
     C     3.03809051092238427E+03,    -1.46792612476956167E+10,
     D     1.14498237732025810E+11,    -3.99096175224466498E+11,
     E     8.19218669548577329E+11,    -1.09837515608122331E+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2     C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3     1.00815810686538209E+12,    -6.45364869245376503E+11,
     4     2.87900649906150589E+11,    -8.78670721780232657E+10,
     5     1.76347306068349694E+10,    -2.16716498322379509E+09,
     6     1.43157876718888981E+08,    -3.87183344257261262E+06,
     7     1.82577554742931747E+04,     2.86464035717679043E+11,
     8    -2.40629790002850396E+12,     9.10934118523989896E+12,
     9    -2.05168994109344374E+13,     3.05651255199353206E+13,
     A    -3.16670885847851584E+13,     2.33483640445818409E+13,
     B    -1.23204913055982872E+13,     4.61272578084913197E+12,
     C    -1.19655288019618160E+12,     2.05914503232410016E+11,
     D    -2.18229277575292237E+10,     1.24700929351271032E+09/
      DATA C(119), C(120)/
     1    -2.91883881222208134E+07,     1.18838426256783253E+05/
C
      IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
      RFN = 1.0E0/FNU
      CRFN = CMPLX(RFN,0.0E0)
C     T = ZR*CRFN
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TSTR = REAL(ZR)
      TSTI = AIMAG(ZR)
      TEST = R1MACH(1)*1.0E+3
      AC = FNU*TEST
      IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
      AC = 2.0E0*ABS(ALOG(TEST))+FNU
      ZETA1 = CMPLX(AC,0.0E0)
      ZETA2 = CMPLX(FNU,0.0E0)
      PHI=CONE
      RETURN
   15 CONTINUE
      T=ZR*CRFN
      S = CONE + T*T
      SR = CSQRT(S)
      CFN = CMPLX(FNU,0.0E0)
      ZN = (CONE+SR)/T
      ZETA1 = CFN*CLOG(ZN)
      ZETA2 = CFN*SR
      T = CONE/SR
      SR = T*CRFN
      CWRK(16) = CSQRT(SR)
      PHI = CWRK(16)*CON(IKFLG)
      IF (IPMTR.NE.0) RETURN
      T2 = CONE/S
      CWRK(1) = CONE
      CRFN = CONE
      AC = 1.0E0
      L = 1
      DO 20 K=2,15
        S = CZERO
        DO 10 J=1,K
          L = L + 1
          S = S*T2 + CMPLX(C(L),0.0E0)
   10   CONTINUE
        CRFN = CRFN*SR
        CWRK(K) = CRFN*S
        AC = AC*RFN
        TSTR = REAL(CWRK(K))
        TSTI = AIMAG(CWRK(K))
        TEST = ABS(TSTR) + ABS(TSTI)
        IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20 CONTINUE
      K = 15
   30 CONTINUE
      INIT = K
   40 CONTINUE
      IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
      S = CZERO
      DO 50 I=1,INIT
        S = S + CWRK(I)
   50 CONTINUE
      SUM = S
      PHI = CWRK(16)*CON(1)
      RETURN
   60 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
      S = CZERO
      T = CONE
      DO 70 I=1,INIT
        S = S + T*CWRK(I)
        T = -T
   70 CONTINUE
      SUM = S
      PHI = CWRK(16)*CON(2)
      RETURN
      END
      SUBROUTINE CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CUNK1
C***REFER TO  CBESK
C
C     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  CS1S2,CUCHK,CUNIK,R1MACH
C***END PROLOGUE  CUNK1
      COMPLEX CFN, CK, CONE, CRSC, CS, CSCL, CSGN, CSPN, CSR, CSS,
     * CWRK, CY, CZERO, C1, C2, PHI,  RZ, SUM,  S1, S2, Y, Z,
     * ZETA1,  ZETA2,  ZR, PHID, ZETA1D, ZETA2D, SUMD
      REAL ALIM, ANG, APHI, ASC, ASCLE, BRY, CPN, C2I, C2M, C2R, ELIM,
     * FMR, FN, FNF, FNU, PI, RS1, SGN, SPN, TOL, X, R1MACH
      INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     * KK, KODE, MR, N, NW, NZ, J, IPARD, INITD, IC
      DIMENSION BRY(3), INIT(2), Y(N), SUM(2), PHI(2), ZETA1(2),
     * ZETA2(2), CY(2), CWRK(16,3), CSS(3), CSR(3)
      DATA CZERO, CONE / (0.0E0,0.0E0) , (1.0E0,0.0E0) /
      DATA PI / 3.14159265358979324E0 /
C
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      X = REAL(Z)
      ZR = Z
      IF (X.LT.0.0E0) ZR = -Z
      J=2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + FLOAT(I-1)
        INIT(J) = 0
        CALL CUNIK(ZR, FN, 2, 0, TOL, INIT(J), PHI(J), ZETA1(J),
     *   ZETA2(J), SUM(J), CWRK(1,J))
        IF (KODE.EQ.1) GO TO 20
        CFN = CMPLX(FN,0.0E0)
        S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
        GO TO 30
   20   CONTINUE
        S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHI(J))
        RS1 = RS1 + ALOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        S2 = PHI(J)*SUM(J)
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(KFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (KFLAG.NE.1) GO TO 50
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        CY(KDFLG) = S2
        Y(I) = S2*CSR(KFLAG)
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 290
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (X.LT.0.0E0) GO TO 290
        KDFLG = 1
        Y(I) = CZERO
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF (Y(I-1).EQ.CZERO) GO TO 70
        Y(I-1) = CZERO
        NZ=NZ+1
   70 CONTINUE
      I=N
   75 CONTINUE
      RZ = CMPLX(2.0E0,0.0E0)/ZR
      CK = CMPLX(FN,0.0E0)*RZ
      IB = I+1
      IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
C     ON UNDERFLOW
C-----------------------------------------------------------------------
      FN = FNU+FLOAT(N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      INITD = 0
      CALL CUNIK(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD,
     *CWRK(1,3))
      IF (KODE.EQ.1) GO TO 80
      CFN=CMPLX(FN,0.0E0)
      S1=ZETA1D-CFN*(CFN/(ZR+ZETA2D))
      GO TO 90
   80 CONTINUE
      S1=ZETA1D-ZETA2D
   90 CONTINUE
      RS1=REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 95
      IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI=CABS(PHID)
      RS1=RS1+ALOG(APHI)
      IF (ABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 290
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (X.LT.0.0E0) GO TO 290
      NZ=N
      DO 96 I=1,N
        Y(I) = CZERO
   96 CONTINUE
      RETURN
  100 CONTINUE
C-----------------------------------------------------------------------
C     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2 = S2
        S2 = CK*S2 + S1
        S1 = C2
        CK = CK + RZ
        C2 = S2*C1
        Y(I) = C2
        IF (KFLAG.GE.3) GO TO 120
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        C1 = CSR(KFLAG)
  120 CONTINUE
  160 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
C-----------------------------------------------------------------------
      CSGN = CMPLX(0.0E0,SGN)
      INU = INT(FNU)
      FNF = FNU - FLOAT(INU)
      IFN = INU + N - 1
      ANG = FNF*SGN
      CPN = COS(ANG)
      SPN = SIN(ANG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
      ASC = BRY(1)
      KK = N
      IUF = 0
      KDFLG = 1
      IB = IB-1
      IC = IB-1
      DO 260 K=1,N
        FN = FNU + FLOAT(KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        M=3
        IF (N.GT.2) GO TO 175
  170   CONTINUE
        INITD = INIT(J)
        PHID = PHI(J)
        ZETA1D = ZETA1(J)
        ZETA2D = ZETA2(J)
        SUMD = SUM(J)
        M = J
        J = 3 - J
        GO TO 180
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 170
        INITD = 0
  180   CONTINUE
        CALL CUNIK(ZR, FN, 1, 0, TOL, INITD, PHID, ZETA1D,
     *   ZETA2D, SUMD, CWRK(1,M))
        IF (KODE.EQ.1) GO TO 190
        CFN = CMPLX(FN,0.0E0)
        S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
        GO TO 200
  190   CONTINUE
        S1 = -ZETA1D + ZETA2D
  200   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 250
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 210
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHID)
        RS1 = RS1 + ALOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 250
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 210
        IF (KDFLG.EQ.1) IFLAG = 3
  210   CONTINUE
        S2 = CSGN*PHID*SUMD
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 220
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) S2 = CMPLX(0.0E0,0.0E0)
  220   CONTINUE
        CY(KDFLG) = S2
        C2 = S2
        S2 = S2*CSR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1 = Y(KK)
        IF (KODE.EQ.1) GO TO 240
        CALL CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  240   CONTINUE
        Y(KK) = S1*CSPN + S2
        KK = KK - 1
        CSPN = -CSPN
        IF (C2.NE.CZERO) GO TO 245
        KDFLG = 1
        GO TO 260
  245   CONTINUE
        IF (KDFLG.EQ.2) GO TO 265
        KDFLG = 2
        GO TO 260
  250   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 290
        S2 = CZERO
        GO TO 220
  260 CONTINUE
      K = N
  265 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      CS = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = FLOAT(INU+IL)
      DO 280 I=1,IL
        C2 = S2
        S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
        S1 = C2
        FN = FN - 1.0E0
        C2 = S2*CS
        CK = C2
        C1 = Y(KK)
        IF (KODE.EQ.1) GO TO 270
        CALL CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  270   CONTINUE
        Y(KK) = C1*CSPN + C2
        KK = KK - 1
        CSPN = -CSPN
        IF (IFLAG.GE.3) GO TO 280
        C2R = REAL(CK)
        C2I = AIMAG(CK)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 280
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*CS
        S2 = CK
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        CS = CSR(IFLAG)
  280 CONTINUE
      RETURN
  290 CONTINUE
      NZ = -1
      RETURN
      END
      SUBROUTINE CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CUNK2
C***REFER TO  CBESK
C
C     CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***ROUTINES CALLED  CAIRY,CS1S2,CUCHK,CUNHJ,R1MACH
C***END PROLOGUE  CUNK2
      COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CIP,
     * CK, CONE, CRSC, CR1, CR2, CS, CSCL, CSGN, CSPN, CSR, CSS, CY,
     * CZERO, C1, C2, DAI, PHI,  RZ, S1, S2, Y, Z, ZB, ZETA1,
     * ZETA2, ZN, ZR, PHID, ARGD, ZETA1D, ZETA2D, ASUMD, BSUMD
      REAL AARG, AIC, ALIM, ANG, APHI, ASC, ASCLE, BRY, CAR, CPN, C2I,
     * C2M, C2R, ELIM, FMR, FN, FNF, FNU, HPI, PI, RS1, SAR, SGN, SPN,
     * TOL, X, YY, R1MACH
      INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
      DIMENSION BRY(3), Y(N), ASUM(2), BSUM(2), PHI(2), ARG(2),
     * ZETA1(2), ZETA2(2), CY(2), CIP(4), CSS(3), CSR(3)
      DATA CZERO, CONE, CI, CR1, CR2 /
     1         (0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0),
     1(1.0E0,1.73205080756887729E0),(-0.5E0,-8.66025403784438647E-01)/
      DATA HPI, PI, AIC /
     1     1.57079632679489662E+00,     3.14159265358979324E+00,
     1     1.26551212348464539E+00/
      DATA CIP(1),CIP(2),CIP(3),CIP(4)/
     1 (1.0E0,0.0E0), (0.0E0,-1.0E0), (-1.0E0,0.0E0), (0.0E0,1.0E0)/
C
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      CRSC = CMPLX(TOL,0.0E0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = 1.0E+3*R1MACH(1)/TOL
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = R1MACH(2)
      X = REAL(Z)
      ZR = Z
      IF (X.LT.0.0E0) ZR = -Z
      YY = AIMAG(ZR)
      ZN = -ZR*CI
      ZB = ZR
      INU = INT(FNU)
      FNF = FNU - FLOAT(INU)
      ANG = -HPI*FNF
      CAR = COS(ANG)
      SAR = SIN(ANG)
      CPN = -HPI*CAR
      SPN = -HPI*SAR
      C2 = CMPLX(-SPN,CPN)
      KK = MOD(INU,4) + 1
      CS = CR1*C2*CIP(KK)
      IF (YY.GT.0.0E0) GO TO 10
      ZN = CONJG(-ZN)
      ZB = CONJG(ZB)
   10 CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      J = 2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + FLOAT(I-1)
        CALL CUNHJ(ZN, FN, 0, TOL, PHI(J), ARG(J), ZETA1(J), ZETA2(J),
     *   ASUM(J), BSUM(J))
        IF (KODE.EQ.1) GO TO 20
        CFN = CMPLX(FN,0.0E0)
        S1 = ZETA1(J) - CFN*(CFN/(ZB+ZETA2(J)))
        GO TO 30
   20   CONTINUE
        S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHI(J))
        AARG = CABS(ARG(J))
        RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        C2 = ARG(J)*CR2
        CALL CAIRY(C2, 0, 2, AI, NAI, IDUM)
        CALL CAIRY(C2, 1, 2, DAI, NDAI, IDUM)
        S2 = CS*PHI(J)*(AI*ASUM(J)+CR2*DAI*BSUM(J))
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(KFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (KFLAG.NE.1) GO TO 50
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        IF (YY.LE.0.0E0) S2 = CONJG(S2)
        CY(KDFLG) = S2
        Y(I) = S2*CSR(KFLAG)
        CS = -CI*CS
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (X.LT.0.0E0) GO TO 300
        KDFLG = 1
        Y(I) = CZERO
        CS = -CI*CS
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF (Y(I-1).EQ.CZERO) GO TO 70
        Y(I-1) = CZERO
        NZ=NZ+1
   70 CONTINUE
      I=N
   75 CONTINUE
      RZ = CMPLX(2.0E0,0.0E0)/ZR
      CK = CMPLX(FN,0.0E0)*RZ
      IB = I + 1
      IF (N.LT.IB) GO TO 170
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
C     ON UNDERFLOW
C-----------------------------------------------------------------------
      FN = FNU+FLOAT(N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      CALL CUNHJ(ZN,FN,IPARD,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD,BSUMD)
      IF (KODE.EQ.1) GO TO 80
      CFN=CMPLX(FN,0.0E0)
      S1=ZETA1D-CFN*(CFN/(ZB+ZETA2D))
      GO TO 90
   80 CONTINUE
      S1=ZETA1D-ZETA2D
   90 CONTINUE
      RS1=REAL(S1)
      IF (ABS(RS1).GT.ELIM) GO TO 95
      IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI=CABS(PHID)
      AARG = CABS(ARGD)
      RS1=RS1+ALOG(APHI)-0.25E0*ALOG(AARG)-AIC
      IF (ABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (RS1.GT.0.0E0) GO TO 300
C-----------------------------------------------------------------------
C     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (X.LT.0.0E0) GO TO 300
      NZ=N
      DO 96 I=1,N
        Y(I) = CZERO
   96 CONTINUE
      RETURN
  100 CONTINUE
C-----------------------------------------------------------------------
C     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2 = S2
        S2 = CK*S2 + S1
        S1 = C2
        CK = CK + RZ
        C2 = S2*C1
        Y(I) = C2
        IF (KFLAG.GE.3) GO TO 120
        C2R = REAL(C2)
        C2I = AIMAG(C2)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1 = S1*C1
        S2 = C2
        S1 = S1*CSS(KFLAG)
        S2 = S2*CSS(KFLAG)
        C1 = CSR(KFLAG)
  120 CONTINUE
  170 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = FLOAT(MR)
      SGN = -SIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGN = CMPLX(0.0E0,SGN)
      IF (YY.LE.0.0E0) CSGN = CONJG(CSGN)
      IFN = INU + N - 1
      ANG = FNF*SGN
      CPN = COS(ANG)
      SPN = SIN(ANG)
      CSPN = CMPLX(CPN,SPN)
      IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      CS = CMPLX(CAR,-SAR)*CSGN
      IN = MOD(IFN,4) + 1
      C2 = CIP(IN)
      CS = CS*CONJG(C2)
      ASC = BRY(1)
      KK = N
      KDFLG = 1
      IB = IB-1
      IC = IB-1
      IUF = 0
      DO 270 K=1,N
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        FN = FNU+FLOAT(KK-1)
        IF (N.GT.2) GO TO 180
  175   CONTINUE
        PHID = PHI(J)
        ARGD = ARG(J)
        ZETA1D = ZETA1(J)
        ZETA2D = ZETA2(J)
        ASUMD = ASUM(J)
        BSUMD = BSUM(J)
        J = 3 - J
        GO TO 190
  180   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 190
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 175
        CALL CUNHJ(ZN, FN, 0, TOL, PHID, ARGD, ZETA1D, ZETA2D,
     *   ASUMD, BSUMD)
  190   CONTINUE
        IF (KODE.EQ.1) GO TO 200
        CFN = CMPLX(FN,0.0E0)
        S1 = -ZETA1D + CFN*(CFN/(ZB+ZETA2D))
        GO TO 210
  200   CONTINUE
        S1 = -ZETA1D + ZETA2D
  210   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = REAL(S1)
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = CABS(PHID)
        AARG = CABS(ARGD)
        RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0E0) GO TO 220
        IF (KDFLG.EQ.1) IFLAG = 3
  220   CONTINUE
        CALL CAIRY(ARGD, 0, 2, AI, NAI, IDUM)
        CALL CAIRY(ARGD, 1, 2, DAI, NDAI, IDUM)
        S2 = CS*PHID*(AI*ASUMD+DAI*BSUMD)
        C2R = REAL(S1)
        C2I = AIMAG(S1)
        C2M = EXP(C2R)*REAL(CSS(IFLAG))
        S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
        S2 = S2*S1
        IF (IFLAG.NE.1) GO TO 230
        CALL CUCHK(S2, NW, BRY(1), TOL)
        IF (NW.NE.0) S2 = CMPLX(0.0E0,0.0E0)
  230   CONTINUE
        IF (YY.LE.0.0E0) S2 = CONJG(S2)
        CY(KDFLG) = S2
        C2 = S2
        S2 = S2*CSR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1 = Y(KK)
        IF (KODE.EQ.1) GO TO 250
        CALL CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  250   CONTINUE
        Y(KK) = S1*CSPN + S2
        KK = KK - 1
        CSPN = -CSPN
        CS = -CS*CI
        IF (C2.NE.CZERO) GO TO 255
        KDFLG = 1
        GO TO 270
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 275
        KDFLG = 2
        GO TO 270
  260   CONTINUE
        IF (RS1.GT.0.0E0) GO TO 300
        S2 = CZERO
        GO TO 230
  270 CONTINUE
      K = N
  275 CONTINUE
      IL = N-K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      CS = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = FLOAT(INU+IL)
      DO 290 I=1,IL
        C2 = S2
        S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
        S1 = C2
        FN = FN - 1.0E0
        C2 = S2*CS
        CK = C2
        C1 = Y(KK)
        IF (KODE.EQ.1) GO TO 280
        CALL CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  280   CONTINUE
        Y(KK) = C1*CSPN + C2
        KK = KK - 1
        CSPN = -CSPN
        IF (IFLAG.GE.3) GO TO 290
        C2R = REAL(CK)
        C2I = AIMAG(CK)
        C2R = ABS(C2R)
        C2I = ABS(C2I)
        C2M = AMAX1(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 290
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1 = S1*CS
        S2 = CK
        S1 = S1*CSS(IFLAG)
        S2 = S2*CSS(IFLAG)
        CS = CSR(IFLAG)
  290 CONTINUE
      RETURN
  300 CONTINUE
      NZ = -1
      RETURN
      END
      SUBROUTINE CUOIK(Z, FNU, KODE, IKFLG, N, Y, NUF, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CUOIK
C***REFER TO  CBESI,CBESK,CBESH
C
C     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***ROUTINES CALLED  CUCHK,CUNHJ,CUNIK,R1MACH
C***END PROLOGUE  CUOIK
      COMPLEX ARG, ASUM, BSUM, CWRK, CZ, CZERO, PHI, SUM, Y, Z, ZB,
     * ZETA1, ZETA2, ZN, ZR
      REAL AARG, AIC, ALIM, APHI, ASCLE, AX, AY, ELIM, FNN, FNU, GNN,
     * GNU, RCZ, TOL, X, YY
      INTEGER I, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
      DIMENSION Y(N), CWRK(16)
      DATA CZERO / (0.0E0,0.0E0) /
      DATA AIC / 1.265512123484645396E+00 /
      NUF = 0
      NN = N
      X = REAL(Z)
      ZR = Z
      IF (X.LT.0.0E0) ZR = -Z
      ZB = ZR
      YY = AIMAG(ZR)
      AX = ABS(X)*1.7321E0
      AY = ABS(YY)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = AMAX1(FNU,1.0E0)
      IF (IKFLG.EQ.1) GO TO 10
      FNN = FLOAT(NN)
      GNN = FNU + FNN - 1.0E0
      GNU = AMAX1(GNN,FNN)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
      IF (IFORM.EQ.2) GO TO 20
      INIT = 0
      CALL CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM,
     * CWRK)
      CZ = -ZETA1 + ZETA2
      GO TO 40
   20 CONTINUE
      ZN = -ZR*CMPLX(0.0E0,1.0E0)
      IF (YY.GT.0.0E0) GO TO 30
      ZN = CONJG(-ZN)
   30 CONTINUE
      CALL CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
      CZ = -ZETA1 + ZETA2
      AARG = CABS(ARG)
   40 CONTINUE
      IF (KODE.EQ.2) CZ = CZ - ZB
      IF (IKFLG.EQ.2) CZ = -CZ
      APHI = CABS(PHI)
      RCZ = REAL(CZ)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.GT.ELIM) GO TO 170
      IF (RCZ.LT.ALIM) GO TO 50
      RCZ = RCZ + ALOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
      IF (RCZ.GT.ELIM) GO TO 170
      GO TO 100
   50 CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.LT.(-ELIM)) GO TO 60
      IF (RCZ.GT.(-ALIM)) GO TO 100
      RCZ = RCZ + ALOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 80
   60 CONTINUE
      DO 70 I=1,NN
        Y(I) = CZERO
   70 CONTINUE
      NUF = NN
      RETURN
   80 CONTINUE
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CZ = CZ + CLOG(PHI)
      IF (IFORM.EQ.1) GO TO 90
      CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
   90 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = AIMAG(CZ)
      CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
      CALL CUCHK(CZ, NW, ASCLE, TOL)
      IF (NW.EQ.1) GO TO 60
  100 CONTINUE
      IF (IKFLG.EQ.2) RETURN
      IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  110 CONTINUE
      GNU = FNU + FLOAT(NN-1)
      IF (IFORM.EQ.2) GO TO 120
      INIT = 0
      CALL CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM,
     * CWRK)
      CZ = -ZETA1 + ZETA2
      GO TO 130
  120 CONTINUE
      CALL CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
      CZ = -ZETA1 + ZETA2
      AARG = CABS(ARG)
  130 CONTINUE
      IF (KODE.EQ.2) CZ = CZ - ZB
      APHI = CABS(PHI)
      RCZ = REAL(CZ)
      IF (RCZ.LT.(-ELIM)) GO TO 140
      IF (RCZ.GT.(-ALIM)) RETURN
      RCZ = RCZ + ALOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 150
  140 CONTINUE
      Y(NN) = CZERO
      NN = NN - 1
      NUF = NUF + 1
      IF (NN.EQ.0) RETURN
      GO TO 110
  150 CONTINUE
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CZ = CZ + CLOG(PHI)
      IF (IFORM.EQ.1) GO TO 160
      CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
  160 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = AIMAG(CZ)
      CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
      CALL CUCHK(CZ, NW, ASCLE, TOL)
      IF (NW.EQ.1) GO TO 140
      RETURN
  170 CONTINUE
      NUF = -1
      RETURN
      END
      SUBROUTINE CWRSK(ZR, FNU, KODE, N, Y, NZ, CW, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CWRSK
C***REFER TO  CBESI,CBESK
C
C     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN
C
C***ROUTINES CALLED  CBKNU,CRATI,R1MACH
C***END PROLOGUE  CWRSK
      COMPLEX CINU, CSCL, CT, CW, C1, C2, RCT, ST, Y, ZR
      REAL ACT, ACW, ALIM, ASCLE, ELIM, FNU, S1, S2, TOL, YY
      INTEGER I, KODE, N, NW, NZ
      DIMENSION Y(N), CW(2)
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
      NZ = 0
      CALL CBKNU(ZR, FNU, KODE, 2, CW, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 50
      CALL CRATI(ZR, FNU, N, Y, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
      CINU = CMPLX(1.0E0,0.0E0)
      IF (KODE.EQ.1) GO TO 10
      YY = AIMAG(ZR)
      S1 = COS(YY)
      S2 = SIN(YY)
      CINU = CMPLX(S1,S2)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
      ACW = CABS(CW(2))
      ASCLE = 1.0E+3*R1MACH(1)/TOL
      CSCL = CMPLX(1.0E0,0.0E0)
      IF (ACW.GT.ASCLE) GO TO 20
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0E0/ASCLE
      IF (ACW.LT.ASCLE) GO TO 30
      CSCL = CMPLX(TOL,0.0E0)
   30 CONTINUE
      C1 = CW(1)*CSCL
      C2 = CW(2)*CSCL
      ST = Y(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0E0/CABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
C-----------------------------------------------------------------------
      CT = ZR*(C2+ST*C1)
      ACT = CABS(CT)
      RCT = CMPLX(1.0E0/ACT,0.0E0)
      CT = CONJG(CT)*RCT
      CINU = CINU*RCT*CT
      Y(1) = CINU*CSCL
      IF (N.EQ.1) RETURN
      DO 40 I=2,N
        CINU = ST*CINU
        ST = Y(I)
        Y(I) = CINU*CSCL
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      FUNCTION GAMLN(Z,IERR)
C***BEGIN PROLOGUE  GAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT
C           Z      - REAL ARGUMENT, Z.GT.0.0E0
C
C         OUTPUT
C           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                    IERR=1, Z.LE.0.0E0,    NO COMPUTATION
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,R1MACH
C***END PROLOGUE  GAMLN
C
      INTEGER I, I1M, K, MZ, NZ, IERR, I1MACH
      REAL CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST, T1, WDTOL, Z,
     * ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ
      REAL R1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4     0.00000000000000000E+00,     0.00000000000000000E+00,
     5     6.93147180559945309E-01,     1.79175946922805500E+00,
     6     3.17805383034794562E+00,     4.78749174278204599E+00,
     7     6.57925121201010100E+00,     8.52516136106541430E+00,
     8     1.06046029027452502E+01,     1.28018274800814696E+01,
     9     1.51044125730755153E+01,     1.75023078458738858E+01,
     A     1.99872144956618861E+01,     2.25521638531234229E+01,
     B     2.51912211827386815E+01,     2.78992713838408916E+01,
     C     3.06718601060806728E+01,     3.35050734501368889E+01,
     D     3.63954452080330536E+01,     3.93398841871994940E+01,
     E     4.23356164607534850E+01,     4.53801388984769080E+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4     4.84711813518352239E+01,     5.16066755677643736E+01,
     5     5.47847293981123192E+01,     5.80036052229805199E+01,
     6     6.12617017610020020E+01,     6.45575386270063311E+01,
     7     6.78897431371815350E+01,     7.12570389671680090E+01,
     8     7.46582363488301644E+01,     7.80922235533153106E+01,
     9     8.15579594561150372E+01,     8.50544670175815174E+01,
     A     8.85808275421976788E+01,     9.21361756036870925E+01,
     B     9.57196945421432025E+01,     9.93306124547874269E+01,
     C     1.02968198614513813E+02,     1.06631760260643459E+02,
     D     1.10320639714757395E+02,     1.14034211781461703E+02,
     E     1.17771881399745072E+02,     1.21533081515438634E+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4     1.25317271149356895E+02,     1.29123933639127215E+02,
     5     1.32952575035616310E+02,     1.36802722637326368E+02,
     6     1.40673923648234259E+02,     1.44565743946344886E+02,
     7     1.48477766951773032E+02,     1.52409592584497358E+02,
     8     1.56360836303078785E+02,     1.60331128216630907E+02,
     9     1.64320112263195181E+02,     1.68327445448427652E+02,
     A     1.72352797139162802E+02,     1.76395848406997352E+02,
     B     1.80456291417543771E+02,     1.84533828861449491E+02,
     C     1.88628173423671591E+02,     1.92739047287844902E+02,
     D     1.96866181672889994E+02,     2.01009316399281527E+02,
     E     2.05168199482641199E+02,     2.09342586752536836E+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4     2.13532241494563261E+02,     2.17736934113954227E+02,
     5     2.21956441819130334E+02,     2.26190548323727593E+02,
     6     2.30439043565776952E+02,     2.34701723442818268E+02,
     7     2.38978389561834323E+02,     2.43268849002982714E+02,
     8     2.47572914096186884E+02,     2.51890402209723194E+02,
     9     2.56221135550009525E+02,     2.60564940971863209E+02,
     A     2.64921649798552801E+02,     2.69291097651019823E+02,
     B     2.73673124285693704E+02,     2.78067573440366143E+02,
     C     2.82474292687630396E+02,     2.86893133295426994E+02,
     D     2.91323950094270308E+02,     2.95766601350760624E+02,
     E     3.00220948647014132E+02,     3.04686856765668715E+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2     3.09164193580146922E+02,     3.13652829949879062E+02,
     3     3.18152639620209327E+02,     3.22663499126726177E+02,
     4     3.27185287703775217E+02,     3.31717887196928473E+02,
     5     3.36261181979198477E+02,     3.40815058870799018E+02,
     6     3.45379407062266854E+02,     3.49954118040770237E+02,
     7     3.54539085519440809E+02,     3.59134205369575399E+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3     8.33333333333333333E-02,    -2.77777777777777778E-03,
     4     7.93650793650793651E-04,    -5.95238095238095238E-04,
     5     8.41750841750841751E-04,    -1.91752691752691753E-03,
     6     6.41025641025641026E-03,    -2.95506535947712418E-02,
     7     1.79644372368830573E-01,    -1.39243221690590112E+00,
     8     1.34028640441683920E+01,    -1.56848284626002017E+02,
     9     2.19310333333333333E+03,    -3.61087712537249894E+04,
     A     6.91472268851313067E+05,    -1.52382215394074162E+07,
     B     3.82900751391414141E+08,    -1.08822660357843911E+10,
     C     3.47320283765002252E+11,    -1.23696021422692745E+13,
     D     4.88788064793079335E+14,    -2.13203339609193739E+16/
C
C             LN(2*PI)
      DATA CON                    /     1.83787706640934548E+00/
C
C***FIRST EXECUTABLE STATEMENT  GAMLN
      IERR=0
      IF (Z.LE.0.0E0) GO TO 70
      IF (Z.GT.101.0E0) GO TO 10
      NZ = INT(Z)
      FZ = Z - FLOAT(NZ)
      IF (FZ.GT.0.0E0) GO TO 10
      IF (NZ.GT.100) GO TO 10
      GAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = R1MACH(4)
      WDTOL = AMAX1(WDTOL,0.5E-18)
      I1M = I1MACH(11)
      RLN = R1MACH(5)*FLOAT(I1M)
      FLN = AMIN1(RLN,20.0E0)
      FLN = AMAX1(FLN,3.0E0)
      FLN = FLN - 3.0E0
      ZM = 1.8000E0 + 0.3875E0*FLN
      MZ = INT(ZM) + 1
      ZMIN = FLOAT(MZ)
      ZDMY = Z
      ZINC = 0.0E0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - FLOAT(NZ)
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0E0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (ABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0E0) GO TO 50
      TLG = ALOG(Z)
      GAMLN = Z*(TLG-1.0E0) + 0.5E0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0E0
      NZ = INT(ZINC)
      DO 60 I=1,NZ
        ZP = ZP*(Z+FLOAT(I-1))
   60 CONTINUE
      TLG = ALOG(ZDMY)
      GAMLN = ZDMY*(TLG-1.0E0) - ALOG(ZP) + 0.5E0*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      IERR=1
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / 2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -126 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / 2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  129 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /    7 /
C     DATA IMACH( 2) /    2 /
C     DATA IMACH( 3) /    2 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -256 /
C     DATA IMACH(13) /  255 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) / -256 /
C     DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -50 /
C     DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -32754 /
C     DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     7 /
C     DATA IMACH( 4) /     6 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -4095 /
C     DATA IMACH(13) /  4094 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -4095 /
C     DATA IMACH(16) /  4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   47 /
C     DATA IMACH(12) / -929 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   94 /
C     DATA IMACH(15) / -929 /
C     DATA IMACH(16) / 1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    0 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) / Z'7FFFFFFF' /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -126 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     7/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    31/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -128/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1024/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /   11 /
C     DATA IMACH( 2) /   12 /
C     DATA IMACH( 3) /    8 /
C     DATA IMACH( 4) /   10 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) /32767 /
C     DATA IMACH(10) /   16 /
C     DATA IMACH(11) /    6 /
C     DATA IMACH(12) /  -64 /
C     DATA IMACH(13) /   63 /
C     DATA IMACH(14) /   14 /
C     DATA IMACH(15) /  -64 /
C     DATA IMACH(16) /   63 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     6/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    32/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -126/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1022/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /       5 /
C     DATA IMACH( 2) /       6 /
C     DATA IMACH( 3) /       0 /
C     DATA IMACH( 4) /       6 /
C     DATA IMACH( 5) /      24 /
C     DATA IMACH( 6) /       3 /
C     DATA IMACH( 7) /       2 /
C     DATA IMACH( 8) /      23 /
C     DATA IMACH( 9) / 8388607 /
C     DATA IMACH(10) /       2 /
C     DATA IMACH(11) /      23 /
C     DATA IMACH(12) /    -127 /
C     DATA IMACH(13) /     127 /
C     DATA IMACH(14) /      38 /
C     DATA IMACH(15) /    -127 /
C     DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /   43 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   63 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5/
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     39 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5 /
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     55 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH(1)  /    5 /
C     DATA IMACH(2)  /    6 /
C     DATA IMACH(3)  /    6 /
C     DATA IMACH(3)  /    7 /
C     DATA IMACH(5)  /   32 /
C     DATA IMACH(6)  /    4 /
C     DATA IMACH(7)  /    2 /
C     DATA IMACH(8)  /   32 /
C     DATA IMACH(9)  /2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -126 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   53 /
C     DATA IMACH(15) /-1015 /
C     DATA IMACH(16) / 1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     0 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -125 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1021 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   32 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     6 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    23 /
C     DATA IMACH(12) /  -126 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    52 /
C     DATA IMACH(15) / -1022 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  128 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    1 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) /-1024 /
C     DATA IMACH(16) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780, G-FLOAT OPTION
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1022 /
C     DATA IMACH(16)/  1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /     1/
C     DATA IMACH( 2) /     1/
C     DATA IMACH( 3) /     0/
C     DATA IMACH( 4) /     1/
C     DATA IMACH( 5) /    16/
C     DATA IMACH( 6) /     2/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    15/
C     DATA IMACH( 9) / 32767/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -127/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    56/
C     DATA IMACH(15) /  -127/
C     DATA IMACH(16) /   127/
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
C
      STOP
      END
*DECK R1MACH
      REAL FUNCTION R1MACH(I)
C***BEGIN PROLOGUE  R1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  890213   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=SINGLE PRECISION(R1MACH-S D1MACH-D),
C             MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns single precision machine dependent constants
C***DESCRIPTION
C
C   R1MACH can be used to obtain machine-dependent parameters
C   for the local machine environment.  It is a function
C   subroutine with one (input) argument, and can be called
C   as follows, for example
C
C        A = R1MACH(I)
C
C   where I=1,...,5.  The (output) value of A above is
C   determined by the (input) value of I.  The results for
C   various values of I are discussed below.
C
C   Single-Precision Machine Constants
C   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   R1MACH(3) = B**(-T), the smallest relative spacing.
C   R1MACH(4) = B**(1-T), the largest relative spacing.
C   R1MACH(5) = LOG10(B)
C
C   Assume single precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(11) = T, the number of base-B digits.
C   I1MACH(12) = EMIN, the smallest exponent E.
C   I1MACH(13) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment,
C   the desired set of DATA statements should be activated by
C   removing the C from column 1.  Also, the values of
C   R1MACH(1) - R1MACH(4) should be checked for consistency
C   with the local operating system.
C
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
      SAVE RMACH
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7EFFFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1) / 16#00800000 /
C     DATA LARGE(1) / 16#7FFFFFFF /
C     DATA RIGHT(1) / 16#33800000 /
C     DATA DIVER(1) / 16#34000000 /
C     DATA LOG10(1) / 16#3E9A209B /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA RMACH(1) / Z"3001800000000000" /
C     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
C     DATA RMACH(3) / Z"3FD2800000000000" /
C     DATA RMACH(4) / Z"3FD3800000000000" /
C     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX C-1
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7FFFFFFF'X /
C     DATA RIGHT(1) / '34800000'X /
C     DATA DIVER(1) / '35000000'X /
C     DATA LOG10(1) / '3F9A209B'X /
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL /    20K,       0 /
C     DATA LARGE / 77777K, 177777K /
C     DATA RIGHT / 35420K,       0 /
C     DATA DIVER / 36020K,       0 /
C     DATA LOG10 / 40423K,  42023K /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE91), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1) / 00004000000B /
C     DATA LARGE(1) / 17677777777B /
C     DATA RIGHT(1) / 06340000000B /
C     DATA DIVER(1) / 06400000000B /
C     DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C       ASSUMING REAL*4 IS THE DEFAULT REAL
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7F7FFFFF'X /
C     DATA RIGHT(1) / '33800000'X /
C     DATA DIVER(1) / '34000000'X /
C     DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA SMALL(1) /     8420761 /
C     DATA LARGE(1) /  2139081118 /
C     DATA RIGHT(1) /   863997169 /
C     DATA DIVER(1) /   872385777 /
C     DATA LOG10(1) /  1050288283 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /   128,     0 /
C     DATA LARGE(1), LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
C     DATA DIVER(1), DIVER(2) / 13568,     0 /
C     DATA LOG10(1), LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1), DIVER(2) / O032400, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS
C
c     data rmach(1) / 1.17549 424 e-38 /
c     data rmach(2) / 3.40282 356 e+38 /
c     data rmach(3) / 1.19209 290 e-07 /
c     data rmach(4) / 2.38418 579 e-07 /
c     data rmach(5) / 0.30103 001 /
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'34000000' /
C     DATA DIVER(1) / Z'34800000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1), SMALL(2) /     0,    256/
C     DATA LARGE(1), LARGE(2) /    -1,   -129/
C     DATA RIGHT(1), RIGHT(2) /     0,  26880/
C     DATA DIVER(1), DIVER(2) /     0,  27136/
C     DATA LOG10(1), LOG10(2) /  8347,  32538/
C
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ('R1MACH -- I OUT OF BOUNDS', 25, 1, 2)
C
      R1MACH = RMACH(I)
      RETURN
C
      END
      SUBROUTINE XERROR(MESS,NMESS,L1,L2)
C
C     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
C     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
C     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
C     ROUTINE.
C
      CHARACTER*(*) MESS
      NN=NMESS/70
      NR=NMESS-70*NN
      IF(NR.NE.0) NN=NN+1
      K=1
      PRINT 900
  900 FORMAT(/)
      DO 10 I=1,NN
        KMIN=MIN0(K+69,NMESS)
        PRINT *, MESS(K:KMIN)
        K=K+70
   10 CONTINUE
      PRINT 900
      RETURN
      END


      double precision function  besselk(order, x, J)
      double precision x
      integer order
      integer j
c     if j =1 returns the besselk
c     if j=2 returns the besselk * exp(x)
      COMPLEX CY, Z
      DIMENSION CY(1)
      INTEGER IERR, NZ
      Z = CMPLX(x, 0d0)
      call cbesk(Z, DBLE(order), J, 1, CY, NZ, IERR)
      if (IERR.eq.0)then
         besselk = DBLE(CY(1))
      else
         besselk = 0d0
      endif
      return
      end
      
