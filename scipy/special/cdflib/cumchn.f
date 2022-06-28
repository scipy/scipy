      SUBROUTINE cumchn(x,df,pnonc,cum,ccum)
C***********************************************************************
C
C     SUBROUTINE CUMCHN(X,DF,PNONC,CUM,CCUM)
C             CUMulative of the Non-central CHi-square distribution
C
C                               Function
C
C     Calculates     the       cumulative      non-central    chi-square
C     distribution, i.e.,  the probability   that  a   random   variable
C     which    follows  the  non-central chi-square  distribution,  with
C     non-centrality  parameter    PNONC  and   continuous  degrees   of
C     freedom DF, is less than or equal to X.
C
C                              Arguments
C
C     X       --> Upper limit of integration of the non-central
C                 chi-square distribution.
C                                                 X is DOUBLE PRECISION
C
C     DF      --> Degrees of freedom of the non-central
C                 chi-square distribution.
C                                                 DF is DOUBLE PRECISION
C
C     PNONC   --> Non-centrality parameter of the non-central
C                 chi-square distribution.
C                                                 PNONC is DOUBLE PRECIS
C
C     CUM <-- Cumulative non-central chi-square distribution.
C                                                 CUM is DOUBLE PRECISIO
C
C     CCUM <-- Compliment of Cumulative non-central chi-square distribut
C                                                 CCUM is DOUBLE PRECISI
C
C
C                                Method
C
C     Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of
C     Mathematical    Functions,  US   NBS   (1966)    to calculate  the
C     non-central chi-square.
C
C                                Variables
C
C     EPS     --- Convergence criterion.  The sum stops when a
C                 term is less than EPS*SUM.
C                                                 EPS is DOUBLE PRECISIO
C
C     CCUM <-- Compliment of Cumulative non-central
C              chi-square distribution.
C                                                 CCUM is DOUBLE PRECISI
C
C***********************************************************************
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION ccum,cum,df,pnonc,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION adj,centaj,centwt,chid2,dfd2,eps,lcntaj,lcntwt,
     +                 lfact,pcent,pterm,sum,sumadj,term,wt,xnonc,xx,
     +                 abstol
      INTEGER i,icent
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alngam
      EXTERNAL alngam
C     ..
C     .. External Subroutines ..
      EXTERNAL cumchi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,exp,int,log
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION dg
      LOGICAL qsmall
C     ..
C     .. Data statements ..
      DATA eps/1.0D-15/
      DATA abstol/1.0D-300/
C     ..
C     .. Statement Function definitions ..
      qsmall(xx) = .NOT. (sum .GE. abstol .AND. xx .GE. eps*sum)
      dg(i) = df + 2.0D0*dble(i)
C     ..
C
      IF (.NOT. (x.LE.0.0D0)) GO TO 10
      cum = 0.0D0
      ccum = 1.0D0
      RETURN

   10 IF (.NOT. (pnonc.LE.1.0D-10)) GO TO 20
C
C
C     When non-centrality parameter is (essentially) zero,
C     use cumulative chi-square distribution
C
C
      CALL cumchi(x,df,cum,ccum)
      RETURN

   20 xnonc = pnonc/2.0D0
C***********************************************************************
C
C     The following code calculates the weight, chi-square, and
C     adjustment term for the central term in the infinite series.
C     The central term is the one in which the poisson weight is
C     greatest.  The adjustment term is the amount that must
C     be subtracted from the chi-square to move up two degrees
C     of freedom.
C
C***********************************************************************
      icent = int(xnonc)
      IF (icent.EQ.0) icent = 1
      chid2 = x/2.0D0
C
C
C     Calculate central weight term
C
C
      lfact = alngam(dble(icent+1))
      lcntwt = -xnonc + icent*log(xnonc) - lfact
      centwt = exp(lcntwt)
C
C
C     Calculate central chi-square
C
C
      CALL cumchi(x,dg(icent),pcent,ccum)
C
C
C     Calculate central adjustment term
C
C
      dfd2 = dg(icent)/2.0D0
      lfact = alngam(1.0D0+dfd2)
      lcntaj = dfd2*log(chid2) - chid2 - lfact
      centaj = exp(lcntaj)
      sum = centwt*pcent
C***********************************************************************
C
C     Sum backwards from the central term towards zero.
C     Quit whenever either
C     (1) the zero term is reached, or
C     (2) the term gets small relative to the sum, or
C
C***********************************************************************
      sumadj = 0.0D0
      adj = centaj
      wt = centwt
      i = icent
C
      GO TO 40

   30 IF (qsmall(term) .OR. i.EQ.0) GO TO 50
   40 dfd2 = dg(i)/2.0D0
C
C
C     Adjust chi-square for two fewer degrees of freedom.
C     The adjusted value ends up in PTERM.
C
C
      adj = adj*dfd2/chid2
      sumadj = sumadj + adj
      pterm = pcent + sumadj
C
C
C     Adjust poisson weight for J decreased by one
C
C
      wt = wt* (i/xnonc)
      term = wt*pterm
      sum = sum + term
      i = i - 1
      GO TO 30

   50 sumadj = centaj
C***********************************************************************
C
C     Now sum forward from the central term towards infinity.
C     Quit when either
C     (1) the term gets small relative to the sum, or
C
C***********************************************************************
      adj = centaj
      wt = centwt
      i = icent
C
      GO TO 70

   60 IF (qsmall(term)) GO TO 80
C
C
C     Update weights for next higher J
C
C
   70 wt = wt* (xnonc/ (i+1))
C
C
C     Calculate PTERM and add term to sum
C
C
      pterm = pcent - sumadj
      term = wt*pterm
      sum = sum + term
C
C
C     Update adjustment term for DF for next iteration
C
C
      i = i + 1
      dfd2 = dg(i)/2.0D0
      adj = adj*chid2/dfd2
      sumadj = sumadj + adj
      GO TO 60

   80 cum = sum
      ccum = 0.5D0 + (0.5D0-cum)
C
      RETURN

      END
