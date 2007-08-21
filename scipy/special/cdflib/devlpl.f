      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
C              Double precision EVALuate a PoLynomial at X
C
C
C                              Function
C
C
C     returns
C          A(1) + A(2)*X + ... + A(N)*X**(N-1)
C
C
C                              Arguments
C
C
C     A --> Array of coefficients of the polynomial.
C                                        A is DOUBLE PRECISION(N)
C
C     N --> Length of A, also degree of polynomial - 1.
C                                        N is INTEGER
C
C     X --> Point at which the polynomial is to be evaluated.
C                                        X is DOUBLE PRECISION
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION a(n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION term
      INTEGER i
C     ..
C     .. Executable Statements ..
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
   10 CONTINUE
      devlpl = term
      RETURN

      END
