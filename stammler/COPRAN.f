CCOPRAN        COPRAN/COLPROB, BARILOCHE 1979.
      SUBROUTINE COPRAN(MMAX, IG, R, SIG, P, GAM)
C
C COPRAN = COLLISION PROBABILITIES IN ANNULAR GEOMETRY.
C          BY INGVAR CARLVIK, 21 NOVEMBER 1968.
C
C MMAX   = NUMBER OF ANNULAR REGIONS.             -- MMAX.LE.20 --
C IG     = NUMBER OF POINTS IN INTEGRATION         -- IG.LE.5 --
C R()    = RADII IN INCREASING ORDER.
C P(I,J) = VOL(I)*SIG(I)*(COLPROB FROM REGION I TO J).
C GAM()  = PARTIAL BLACKNESS.
C TAU()  = OPTICAL DISTANCES.
C SIGV() = SIGMA*VOLUME.
      REAL KI3
      DIMENSION R(1), SIG(1), P(MMAX, 1), GAM(1),
     1          GJC(36), GJP(6,6), GJW(6,6),
     2          R2(20), SIGV(20), TAU(20), TX(21)
      EQUIVALENCE (GJC(1), GJP(1,1), GJW(1,1)), (TX(2), TAU(1))
C
C CONSTANTS FOR GAUSS-JACOBI INTEGRATION.
C ABRAMOWITZ AND STEGUN, PAGE 921, FOR K=1.
C X(I)=1.0-X(I)**2,  W(I)=4.0*W(I),  GJP(I,N).
C W(I) ARE EXTRA DOUBLED FOR INTEGRATION IN TWO DIRECTIONS.
      DATA GJC
     1 /.00000000,.55995556,.87393877,.95491150,.98046718,.99029084,
     2 2.00000000,.00000000,.28606124,.65127016,.82660307,.90725799,
     3  .72783448,1.27216552,.0000000,.16932809,.47704397,.68412769, 
     4  .27930792,.91696442,.80372766,.00000000,.11094751,.35681753,
     5  .12472388,.51939018,.81385828,.54202764, .00000000,.07803490,
     6  .06299166,.29563548,.58554794,.66869856,.38712636,.00000000/
C
      DATA PI /3.1415927/
C
C PREPARE.
      N = MMAX
      Y = 0.0
      DO 10 I = 1, N
      R2(I) = R(I) * R(I)
      SIGV(I) = PI * SIG(I) * (R2(I) - Y)
      Y = R2(I)
      DO 10 J = I, N
   10 P(I,J) = 0.0
C
C EVALUATE FIRST P(I,J) AS INTEGRAL (KI3(TAU(IJ+)) - KI3(TAU(IJ-))).
      DO 20 I = 1, N
      RSTART = 0.0
      IF(I.GT.1) RSTART = R(I-1)
      DR = R(I) - RSTART
      DO 20 K = 1, IG
      Y = RSTART + DR*GJP(IG+1, K)
      FAC = DR*GJW(K, IG+1)
      TAU(I-1) = 0.0
      TPLUS = 0.0
      Y2 = Y*Y
      DO 20 J = I, N
      TMINUS = TPLUS
      TPLUS = SQRT(R2(J) - Y2)
      TAU(J) = TAU(J-1) + SIG(J)*(TPLUS - TMINUS)
      DO 20 L = I, J
   20 P(L, J) = P(L, J) + FAC*(KI3(TAU(J)+TAU(L)) - KI3(TAU(J)-TAU(L)))
C
C COMPOSE P(I, J) FOR BLACK BOUNDARIES.
      DO 40 I = 1, N
      J = N - I + 1
      DO 30 K = J, N
      L = J + N - K
      IF(L.EQ.1) GO TO 30
      IF(L.EQ.J) P(J, L-1) = P(L-1, J)
      P(J, L) = P(J, L) - P(J, L-1)
      IF(J.EQ.1) GO TO 30
      P(J, L) = P(J, L) - P(J-1, L) + P(J-1, L-1)
   30 P(L, J) = P(J, L)
   40 P(J, J) = P(J, J) + SIGV(J)
C
C FIRST-FLIGHT BLACKNESSES.
      S = 2.0 / (PI*R(N))
      DO 60 I = 1, N
      SUM = SIGV(I)
      DO 50 J = 1, N
   50 SUM = SUM - P(I, J)
   60 GAM(I) = S*SUM
C
      RETURN
      END