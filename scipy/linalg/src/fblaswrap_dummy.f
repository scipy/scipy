      subroutine wcdotu (r, n, cx, incx, cy, incy)
      external cdotu
      complex cdotu, r
      integer n
      complex cx (*)
      integer incx
      complex cy (*)
      integer incy
      r = cdotu (n, cx, incx, cy, incy)
      end

      subroutine wzdotu (r, n, zx, incx, zy, incy)
      external zdotu
      double complex zdotu, r
      integer n
      double complex zx (*)
      integer incx
      double complex zy (*)
      integer incy
      r = zdotu (n, zx, incx, zy, incy)
      end

      subroutine wcdotc (r, n, cx, incx, cy, incy)
      external cdotc
      complex cdotc, r
      integer n
      complex cx (*)
      integer incx
      complex cy (*)
      integer incy
      r = cdotc (n, cx, incx, cy, incy)
      end

      subroutine wzdotc (r, n, zx, incx, zy, incy)
      external zdotc
      double complex zdotc, r
      integer n
      double complex zx (*)
      integer incx
      double complex zy (*)
      integer incy
      r = zdotc (n, zx, incx, zy, incy)
      end

      REAL FUNCTION WSDOT( N, SX, INCX, SY, INCY )
      INTEGER INCX, INCY, N
      REAL SX(*), SY(*)
      EXTERNAL SDOT
      REAL SDOT
      WSDOT = SDOT( N, SX, INCX, SY, INCY )
      END FUNCTION
      
      REAL FUNCTION WSASUM( N, SX, INCX )
      INTEGER INCX, N
      REAL SX(*)
      EXTERNAL SASUM
      REAL SASUM
      WSASUM = SASUM( N, SX, INCX )
      END FUNCTION

      REAL FUNCTION WSNRM2( N, SX, INCX )
      INTEGER INCX, N
      REAL SX(*)
      EXTERNAL SNRM2
      REAL SNRM2
      WSNRM2 = SNRM2( N, SX, INCX )
      END FUNCTION

      REAL FUNCTION WSCASUM( N, CX, INCX )
      INTEGER INCX, N
      COMPLEX CX(*)
      EXTERNAL SCASUM
      REAL SCASUM
      WSCASUM = SCASUM( N, CX, INCX )
      END FUNCTION

      REAL FUNCTION WSCNRM2( N, CX, INCX )
      INTEGER INCX, N
      COMPLEX CX(*)
      EXTERNAL SCNRM2
      REAL SCNRM2
      WSCNRM2 = SCNRM2( N, CX, INCX )
      END FUNCTION
