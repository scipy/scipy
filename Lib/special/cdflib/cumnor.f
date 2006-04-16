      SUBROUTINE cumnor(arg,result,ccum)
C**********************************************************************
C
C     SUBROUINE CUMNOR(X,RESULT,CCUM)
C
C
C                              Function
C
C
C     Computes the cumulative  of    the  normal   distribution,   i.e.,
C     the integral from -infinity to x of
C          (1/sqrt(2*pi)) exp(-u*u/2) du
C
C     X --> Upper limit of integration.
C                                        X is DOUBLE PRECISION
C
C     RESULT <-- Cumulative normal distribution.
C                                        RESULT is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative normal distribution.
C                                        CCUM is DOUBLE PRECISION
C
C
C     Renaming of function ANORM from:
C
C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
C     Package of Special Function Routines and Test Drivers"
C     acm Transactions on Mathematical Software. 19, 22-32.
C
C     with slight modifications to return ccum and to deal with
C     machine constants.
C
C**********************************************************************
C
C
C Original Comments:
C------------------------------------------------------------------
C
C This function evaluates the normal distribution function:
C
C                              / x
C                     1       |       -t*t/2
C          P(x) = ----------- |      e       dt
C                 sqrt(2 pi)  |
C                             /-oo
C
C   The main computation evaluates near-minimax approximations
C   derived from those in "Rational Chebyshev approximations for
C   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
C   This transportable program uses rational functions that
C   theoretically approximate the normal distribution function to
C   at least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.
C
C   MIN   = smallest machine representable number.
C
C   EPS   = argument below which anorm(x) may be represented by
C           0.5  and above which  x*x  will not underflow.
C           A conservative value is the largest machine number X
C           such that   1.0 + X = 1.0   to machine precision.
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ANORM = 0     for  ARG .LE. XLOW.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 15, 1992
C
C------------------------------------------------------------------
      INTEGER i
      DOUBLE PRECISION a,arg,b,c,d,del,eps,half,p,one,q,result,sixten,
     +                 temp,sqrpi,thrsh,root32,x,xden,xnum,y,xsq,zero,
     +                 min,ccum
      DIMENSION a(5),b(4),c(9),d(8),p(6),q(5)
C------------------------------------------------------------------
C  External Function
C------------------------------------------------------------------
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C------------------------------------------------------------------
C  Mathematical constants
C
C  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
C  THRSH is the argument for which anorm = 0.75.
C------------------------------------------------------------------
      DATA one,half,zero,sixten/1.0D0,0.5D0,0.0D0,1.60D0/,
     +     sqrpi/3.9894228040143267794D-1/,thrsh/0.66291D0/,
     +     root32/5.656854248D0/
C------------------------------------------------------------------
C  Coefficients for approximation in first interval
C------------------------------------------------------------------
      DATA a/2.2352520354606839287D00,1.6102823106855587881D02,
     +     1.0676894854603709582D03,1.8154981253343561249D04,
     +     6.5682337918207449113D-2/
      DATA b/4.7202581904688241870D01,9.7609855173777669322D02,
     +     1.0260932208618978205D04,4.5507789335026729956D04/
C------------------------------------------------------------------
C  Coefficients for approximation in second interval
C------------------------------------------------------------------
      DATA c/3.9894151208813466764D-1,8.8831497943883759412D00,
     +     9.3506656132177855979D01,5.9727027639480026226D02,
     +     2.4945375852903726711D03,6.8481904505362823326D03,
     +     1.1602651437647350124D04,9.8427148383839780218D03,
     +     1.0765576773720192317D-8/
      DATA d/2.2266688044328115691D01,2.3538790178262499861D02,
     +     1.5193775994075548050D03,6.4855582982667607550D03,
     +     1.8615571640885098091D04,3.4900952721145977266D04,
     +     3.8912003286093271411D04,1.9685429676859990727D04/
C------------------------------------------------------------------
C  Coefficients for approximation in third interval
C------------------------------------------------------------------
      DATA p/2.1589853405795699D-1,1.274011611602473639D-1,
     +     2.2235277870649807D-2,1.421619193227893466D-3,
     +     2.9112874951168792D-5,2.307344176494017303D-2/
      DATA q/1.28426009614491121D00,4.68238212480865118D-1,
     +     6.59881378689285515D-2,3.78239633202758244D-3,
     +     7.29751555083966205D-5/
C------------------------------------------------------------------
C  Machine dependent constants
C------------------------------------------------------------------
      eps = spmpar(1)*0.5D0
      min = spmpar(2)
C------------------------------------------------------------------
      x = arg
      y = abs(x)
      IF (y.LE.thrsh) THEN
C------------------------------------------------------------------
C  Evaluate  anorm  for  |X| <= 0.66291
C------------------------------------------------------------------
          xsq = zero
          IF (y.GT.eps) xsq = x*x
          xnum = a(5)*xsq
          xden = xsq
          DO 10 i = 1,3
              xnum = (xnum+a(i))*xsq
              xden = (xden+b(i))*xsq
   10     CONTINUE
          result = x* (xnum+a(4))/ (xden+b(4))
          temp = result
          result = half + temp
          ccum = half - temp
C------------------------------------------------------------------
C  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
C------------------------------------------------------------------
      ELSE IF (y.LE.root32) THEN
          xnum = c(9)*y
          xden = y
          DO 20 i = 1,7
              xnum = (xnum+c(i))*y
              xden = (xden+d(i))*y
   20     CONTINUE
          result = (xnum+c(8))/ (xden+d(8))
          xsq = aint(y*sixten)/sixten
          del = (y-xsq)* (y+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF
C------------------------------------------------------------------
C  Evaluate  anorm  for |X| > sqrt(32)
C------------------------------------------------------------------
      ELSE
          result = zero
          xsq = one/ (x*x)
          xnum = p(6)*xsq
          xden = xsq
          DO 30 i = 1,4
              xnum = (xnum+p(i))*xsq
              xden = (xden+q(i))*xsq
   30     CONTINUE
          result = xsq* (xnum+p(5))/ (xden+q(5))
          result = (sqrpi-result)/y
          xsq = aint(x*sixten)/sixten
          del = (x-xsq)* (x+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF

      END IF

      IF (result.LT.min) result = 0.0D0
      IF (ccum.LT.min) ccum = 0.0D0
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
C----------Last card of ANORM ----------
      END
