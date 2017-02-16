c Wrappers allowing to link MacOSX's Accelerate framework to
c gfortran compiled code

c Accelerate BLAS is cblas (http://www.netlib.org/blas/blast-forum/cblas.tgz);
c these wrappers call the cblas functions via the C-functions defined
c in veclib_cabi.c

      REAL FUNCTION WSDOT( N, SX, INCX, SY, INCY )
      INTEGER INCX, INCY, N
      REAL SX(*), SY(*)
      REAL RESULT
      EXTERNAL ACC_SDOT
      CALL ACC_SDOT( N, SX, INCX, SY, INCY, RESULT )
      WSDOT = RESULT
      END FUNCTION

      REAL FUNCTION WSDSDOT( N, SB, SX, INCX, SY, INCY )
      REAL SB
      INTEGER INCX, INCY, N
      REAL SX(*), SY(*)
      REAL RESULT
      EXTERNAL ACC_SDSDOT
      CALL ACC_SDSDOT( N, SB, SX, INCX, SY, INCY, RESULT )
      WSDSDOT = RESULT
      END FUNCTION

      REAL FUNCTION WSASUM( N, SX, INCX )
      INTEGER INCX, N
      REAL SX(*)
      REAL RESULT
      EXTERNAL ACC_SASUM
      CALL ACC_SASUM( N, SX, INCX, RESULT )
      WSASUM = RESULT
      END FUNCTION

      REAL FUNCTION WSNRM2( N, SX, INCX )
      INTEGER INCX, N
      REAL SX(*)
      REAL RESULT
      EXTERNAL ACC_SNRM2
      CALL ACC_SNRM2( N, SX, INCX, RESULT )
      WSNRM2 = RESULT
      END FUNCTION

      REAL FUNCTION WSCASUM( N, CX, INCX )
      INTEGER INCX, N
      COMPLEX CX(*)
      REAL RESULT
      EXTERNAL ACC_SCASUM
      CALL ACC_SCASUM( N, CX, INCX, RESULT )
      WSCASUM = RESULT
      END FUNCTION

      REAL FUNCTION WSCNRM2( N, CX, INCX )
      INTEGER INCX, N
      COMPLEX CX(*)
      REAL RESULT
      EXTERNAL ACC_SCNRM2
      CALL ACC_SCNRM2( N, CX, INCX, RESULT )
      WSCNRM2 = RESULT
      END FUNCTION

c The LAPACK in the Accelerate framework is a CLAPACK
c (www.netlib.org/clapack) and has hence a different interface than the
c modern Fortran LAPACK libraries. These wrappers here help to link
c Fortran code to Accelerate.
c This wrapper files covers all Lapack functions that are in all versions
c before Lapack 3.2 (Lapack 3.2 adds CLANHF and SLANSF that would be
c problematic, but those do not exist in OSX <= 10.6, and are actually not
c used in scipy)

      REAL FUNCTION WCLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      CHARACTER          NORM
      INTEGER            KL, KU, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANGB
      DOUBLE PRECISION   CLANGB
      WCLANGB = REAL(CLANGB( NORM, N, KL, KU, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANGE( NORM, M, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, M, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANGE
      DOUBLE PRECISION   CLANGE
      WCLANGE = REAL(CLANGE( NORM, M, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANGT( NORM, N, DL, D, DU )
      CHARACTER          NORM
      INTEGER            N
      COMPLEX            D( * ), DL( * ), DU( * )
      EXTERNAL           CLANGT
      DOUBLE PRECISION   CLANGT
      WCLANGT = REAL(CLANGT( NORM, N, DL, D, DU ))
      END FUNCTION

      REAL FUNCTION WCLANHB( NORM, UPLO, N, K, AB, LDAB, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANHB
      DOUBLE PRECISION   CLANHB
      WCLANHB = REAL(CLANHB( NORM, UPLO, N, K, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANHE( NORM, UPLO, N, A, LDA, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANHE
      DOUBLE PRECISION   CLANHE
      WCLANHE = REAL(CLANHE( NORM, UPLO, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANHP( NORM, UPLO, N, AP, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            N
      REAL               WORK( * )
      COMPLEX            AP( * )
      EXTERNAL           CLANHP
      DOUBLE PRECISION   CLANHP
      WCLANHP = REAL(CLANHP( NORM, UPLO, N, AP, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANHS( NORM, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANHS
      DOUBLE PRECISION   CLANHS
      WCLANHS = REAL(CLANHS( NORM, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANHT( NORM, N, D, E )
      CHARACTER          NORM
      INTEGER            N
      REAL               D( * )
      COMPLEX            E( * )
      EXTERNAL           CLANHT
      DOUBLE PRECISION   CLANHT
      WCLANHT = REAL(CLANHT( NORM, N, D, E ))
      END FUNCTION

      REAL FUNCTION WCLANSB( NORM, UPLO, N, K, AB, LDAB, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANSB
      DOUBLE PRECISION   CLANSB
      WCLANSB = REAL(CLANSB( NORM, UPLO, N, K, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANSP( NORM, UPLO, N, AP, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            N
      REAL               WORK( * )
      COMPLEX            AP( * )
      EXTERNAL           CLANSP
      DOUBLE PRECISION   CLANSP
      WCLANSP = REAL(CLANSP( NORM, UPLO, N, AP, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANSY( NORM, UPLO, N, A, LDA, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANSY
      DOUBLE PRECISION   CLANSY
      WCLANSY = REAL(CLANSY( NORM, UPLO, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANTB
      DOUBLE PRECISION   CLANTB
      WCLANTB = REAL(CLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANTP( NORM, UPLO, DIAG, N, AP, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            N
      REAL               WORK( * )
      COMPLEX            AP( * )
      EXTERNAL           CLANTP
      DOUBLE PRECISION   CLANTP
      WCLANTP = REAL(CLANTP( NORM, UPLO, DIAG, N, AP, WORK ))
      END FUNCTION

      REAL FUNCTION WCLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            LDA, M, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANTR
      DOUBLE PRECISION   CLANTR
      WCLANTR = REAL(CLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WSCSUM1( N, CX, INCX )
      INTEGER            INCX, N
      COMPLEX            CX( * )
      EXTERNAL           SCSUM1
      DOUBLE PRECISION   SCSUM1
      WSCSUM1 = REAL(SCSUM1( N, CX, INCX ))
      END FUNCTION

      REAL FUNCTION WSLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      CHARACTER          NORM
      INTEGER            KL, KU, LDAB, N
      REAL               AB( LDAB, * ), WORK( * )
      EXTERNAL           SLANGB
      DOUBLE PRECISION   SLANGB
      WSLANGB = REAL(SLANGB( NORM, N, KL, KU, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANGE( NORM, M, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, M, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANGE
      DOUBLE PRECISION   SLANGE
      WSLANGE = REAL(SLANGE( NORM, M, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANGT( NORM, N, DL, D, DU )
      CHARACTER          NORM
      INTEGER            N
      REAL               D( * ), DL( * ), DU( * )
      EXTERNAL           SLANGT
      DOUBLE PRECISION   SLANGT
      WSLANGT = REAL(SLANGT( NORM, N, DL, D, DU ))
      END FUNCTION

      REAL FUNCTION WSLANHS( NORM, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANHS
      DOUBLE PRECISION   SLANHS
      WSLANHS = REAL(SLANHS( NORM, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANSB( NORM, UPLO, N, K, AB, LDAB, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               AB( LDAB, * ), WORK( * )
      EXTERNAL           SLANSB
      DOUBLE PRECISION   SLANSB
      WSLANSB = REAL(SLANSB( NORM, UPLO, N, K, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANSP( NORM, UPLO, N, AP, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            N
      REAL               AP( * ), WORK( * )
      EXTERNAL           SLANSP
      DOUBLE PRECISION   SLANSP
      WSLANSP = REAL(SLANSP( NORM, UPLO, N, AP, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANST( NORM, N, D, E )
      CHARACTER          NORM
      INTEGER            N
      REAL               D( * ), E( * )
      EXTERNAL           SLANST
      DOUBLE PRECISION   SLANST
      WSLANST = REAL(SLANST( NORM, N, D, E ))
      END FUNCTION

      REAL FUNCTION WSLANSY( NORM, UPLO, N, A, LDA, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANSY
      DOUBLE PRECISION   SLANSY
      WSLANSY = REAL(SLANSY( NORM, UPLO, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               AB( LDAB, * ), WORK( * )
      EXTERNAL           SLANTB
      DOUBLE PRECISION   SLANTB
      WSLANTB = REAL(SLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANTP( NORM, UPLO, DIAG, N, AP, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            N
      REAL               AP( * ), WORK( * )
      EXTERNAL           SLANTP
      DOUBLE PRECISION   SLANTP
      WSLANTP = REAL(SLANTP( NORM, UPLO, DIAG, N, AP, WORK ))
      END FUNCTION

      REAL FUNCTION WSLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            LDA, M, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANTR
      DOUBLE PRECISION   SLANTR
      WSLANTR = REAL(SLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK ))
      END FUNCTION

      REAL FUNCTION WSLAPY2( X, Y )
      REAL               X, Y
      EXTERNAL           SLAPY2
      DOUBLE PRECISION   SLAPY2
      WSLAPY2 = REAL(SLAPY2( X, Y ))
      END FUNCTION

      REAL FUNCTION WSLAPY3( X, Y, Z )
      REAL               X, Y, Z
      EXTERNAL           SLAPY3
      DOUBLE PRECISION   SLAPY3
      WSLAPY3 = REAL(SLAPY3( X, Y, Z ))
      END FUNCTION

      REAL FUNCTION WSLAMCH( CMACH )
      CHARACTER          CMACH
      EXTERNAL           SLAMCH
      DOUBLE PRECISION   SLAMCH
      WSLAMCH = REAL(SLAMCH( CMACH ))
      END FUNCTION

      REAL FUNCTION WSLAMC3( A, B )
      REAL               A, B
      EXTERNAL           SLAMC3
      DOUBLE PRECISION   SLAMC3
      WSLAMC3 = REAL(SLAMC3( A, B ))
      END FUNCTION
