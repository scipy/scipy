      SUBROUTINE cumf(f,dfn,dfd,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMF(F,DFN,DFD,CUM,CCUM)
C                    CUMulative F distribution
C
C
C                              Function
C
C
C     Computes  the  integral from  0  to  F of  the f-density  with DFN
C     and DFD degrees of freedom.
C
C
C                              Arguments
C
C
C     F --> Upper limit of integration of the f-density.
C                                                  F is DOUBLE PRECISION
C
C     DFN --> Degrees of freedom of the numerator sum of squares.
C                                                  DFN is DOUBLE PRECISI
C
C     DFD --> Degrees of freedom of the denominator sum of squares.
C                                                  DFD is DOUBLE PRECISI
C
C     CUM <-- Cumulative f distribution.
C                                                  CUM is DOUBLE PRECISI
C
C     CCUM <-- Compliment of Cumulative f distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C
C                              Method
C
C
C     Formula  26.5.28 of  Abramowitz and   Stegun   is  used to  reduce
C     the cumulative F to a cumulative beta distribution.
C
C
C                              Note
C
C
C     If F is less than or equal to 0, 0 is returned.
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION dfd,dfn,f,cum,ccum
C     ..
C     .. Local Scalars ..

      DOUBLE PRECISION dsum,prod,xx,yy
      INTEGER ierr
C     ..
C     .. Parameters ..
      DOUBLE PRECISION half
      PARAMETER (half=0.5D0)
      DOUBLE PRECISION done
      PARAMETER (done=1.0D0)
C     ..
C     .. External Subroutines ..
      EXTERNAL bratio
C     ..
C     .. Executable Statements ..

      IF (.NOT. (f.LE.0.0D0)) GO TO 10
      cum = 0.0D0
      ccum = 1.0D0
      RETURN

   10 prod = dfn*f
C
C     XX is such that the incomplete beta with parameters
C     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
C
C     YY is 1 - XX
C
C     Calculate the smaller of XX and YY accurately
C
      dsum = dfd + prod
      xx = dfd/dsum
      IF (xx.GT.half) THEN
          yy = prod/dsum
          xx = done - yy

      ELSE
          yy = done - xx
      END IF

      CALL bratio(dfd*half,dfn*half,xx,yy,ccum,cum,ierr)
      RETURN

      END
