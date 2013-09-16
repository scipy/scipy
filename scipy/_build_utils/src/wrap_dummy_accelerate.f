      REAL FUNCTION WSDOT( N, SX, INCX, SY, INCY )
      INTEGER INCX, INCY, N
      REAL SX(*), SY(*)
      EXTERNAL SDOT
      REAL SDOT
      WSDOT = SDOT( N, SX, INCX, SY, INCY )
      END FUNCTION

      REAL FUNCTION WSDSDOT( N, SB, SX, INCX, SY, INCY )
      REAL SB
      INTEGER INCX, INCY, N
      REAL SX(*), SY(*)
      EXTERNAL SDSDOT
      REAL SDSDOT
      WSDSDOT = SDSDOT( N, SB, SX, INCX, SY, INCY )
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

c The LAPACK in the Accelerate framework is a CLAPACK
c (www.netlib.org/clapack) and has hence a different interface than the
c modern Fortran LAPACK libraries. These wrappers here help to link
c Fortran code to Accelerate.

      REAL FUNCTION WCLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      CHARACTER          NORM
      INTEGER            KL, KU, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANGB
      REAL               CLANGB
      WCLANGB = CLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WCLANGE( NORM, M, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, M, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANGE
      REAL               CLANGE
      WCLANGE = CLANGE( NORM, M, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WCLANGT( NORM, N, DL, D, DU )
      CHARACTER          NORM
      INTEGER            N
      COMPLEX            D( * ), DL( * ), DU( * )
      EXTERNAL           CLANGT
      REAL               CLANGT
      WCLANGT = CLANGT( NORM, N, DL, D, DU )
      END FUNCTION

      REAL FUNCTION WCLANHB( NORM, UPLO, N, K, AB, LDAB, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANHB
      REAL               CLANHB
      WCLANHB = CLANHB( NORM, UPLO, N, K, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WCLANHE( NORM, UPLO, N, A, LDA, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANHE
      REAL               CLANHE
      WCLANHE = CLANHE( NORM, UPLO, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WCLANHP( NORM, UPLO, N, AP, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            N
      REAL               WORK( * )
      COMPLEX            AP( * )
      EXTERNAL           CLANHP
      REAL               CLANHP
      WCLANHP = CLANHP( NORM, UPLO, N, AP, WORK )
      END FUNCTION

      REAL FUNCTION WCLANHS( NORM, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANHS
      REAL               CLANHS
      WCLANHS = CLANHS( NORM, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WCLANHT( NORM, N, D, E )
      CHARACTER          NORM
      INTEGER            N
      REAL               D( * )
      COMPLEX            E( * )
      EXTERNAL           CLANHT
      REAL               CLANHT
      WCLANHT = CLANHT( NORM, N, D, E )
      END FUNCTION

      REAL FUNCTION WCLANSB( NORM, UPLO, N, K, AB, LDAB, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANSB
      REAL               CLANSB
      WCLANSB = CLANSB( NORM, UPLO, N, K, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WCLANSP( NORM, UPLO, N, AP, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            N
      REAL               WORK( * )
      COMPLEX            AP( * )
      EXTERNAL           CLANSP
      REAL               CLANSP
      WCLANSP = CLANSP( NORM, UPLO, N, AP, WORK )
      END FUNCTION

      REAL FUNCTION WCLANSY( NORM, UPLO, N, A, LDA, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANSY
      REAL               CLANSY
      WCLANSY = CLANSY( NORM, UPLO, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WCLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      EXTERNAL           CLANTB
      REAL               CLANTB
      WCLANTB = CLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WCLANTP( NORM, UPLO, DIAG, N, AP, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            N
      REAL               WORK( * )
      COMPLEX            AP( * )
      EXTERNAL           CLANTP
      REAL               CLANTP
      WCLANTP = CLANTP( NORM, UPLO, DIAG, N, AP, WORK )
      END FUNCTION

      REAL FUNCTION WCLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            LDA, M, N
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      EXTERNAL           CLANTR
      REAL               CLANTR
      WCLANTR = CLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WSCSUM1( N, CX, INCX )
      INTEGER            INCX, N
      COMPLEX            CX( * )
      EXTERNAL           SCSUM1
      REAL               SCSUM1
      WSCSUM1 = SCSUM1( N, CX, INCX )
      END FUNCTION

      REAL FUNCTION WSLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      CHARACTER          NORM
      INTEGER            KL, KU, LDAB, N
      REAL               AB( LDAB, * ), WORK( * )
      EXTERNAL           SLANGB
      REAL               SLANGB
      WSLANGB = SLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WSLANGE( NORM, M, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, M, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANGE
      REAL               SLANGE
      WSLANGE = SLANGE( NORM, M, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WSLANGT( NORM, N, DL, D, DU )
      CHARACTER          NORM
      INTEGER            N
      REAL               D( * ), DL( * ), DU( * )
      EXTERNAL           SLANGT
      REAL               SLANGT
      WSLANGT = SLANGT( NORM, N, DL, D, DU )
      END FUNCTION

      REAL FUNCTION WSLANHS( NORM, N, A, LDA, WORK )
      CHARACTER          NORM
      INTEGER            LDA, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANHS
      REAL               SLANHS
      WSLANHS = SLANHS( NORM, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WSLANSB( NORM, UPLO, N, K, AB, LDAB, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               AB( LDAB, * ), WORK( * )
      EXTERNAL           SLANSB
      REAL               SLANSB
      WSLANSB = SLANSB( NORM, UPLO, N, K, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WSLANSP( NORM, UPLO, N, AP, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            N
      REAL               AP( * ), WORK( * )
      EXTERNAL           SLANSP
      REAL               SLANSP
      WSLANSP = SLANSP( NORM, UPLO, N, AP, WORK )
      END FUNCTION

      REAL FUNCTION WSLANST( NORM, N, D, E )
      CHARACTER          NORM
      INTEGER            N
      REAL               D( * ), E( * )
      EXTERNAL           SLANST
      REAL               SLANST
      WSLANST = SLANST( NORM, N, D, E )
      END FUNCTION

      REAL FUNCTION WSLANSY( NORM, UPLO, N, A, LDA, WORK )
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANSY
      REAL               SLANSY
      WSLANSY = SLANSY( NORM, UPLO, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WSLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            K, LDAB, N
      REAL               AB( LDAB, * ), WORK( * )
      EXTERNAL           SLANTB
      REAL               SLANTB
      WSLANTB = SLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )
      END FUNCTION

      REAL FUNCTION WSLANTP( NORM, UPLO, DIAG, N, AP, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            N
      REAL               AP( * ), WORK( * )
      EXTERNAL           SLANTP
      REAL               SLANTP
      WSLANTP = SLANTP( NORM, UPLO, DIAG, N, AP, WORK )
      END FUNCTION

      REAL FUNCTION WSLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            LDA, M, N
      REAL               A( LDA, * ), WORK( * )
      EXTERNAL           SLANTR
      REAL               SLANTR
      WSLANTR = SLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )
      END FUNCTION

      REAL FUNCTION WSLAPY2( X, Y )
      REAL               X, Y
      EXTERNAL           SLAPY2
      REAL               SLAPY2
      WSLAPY2 = SLAPY2( X, Y )
      END FUNCTION

      REAL FUNCTION WSLAPY3( X, Y, Z )
      REAL               X, Y, Z
      EXTERNAL           SLAPY3
      REAL               SLAPY3
      WSLAPY3 = SLAPY3( X, Y, Z )
      END FUNCTION

      REAL FUNCTION WSLAMCH( CMACH )
      CHARACTER          CMACH
      EXTERNAL           SLAMCH
      REAL               SLAMCH
      WSLAMCH = SLAMCH( CMACH )
      END FUNCTION

      REAL FUNCTION WSLAMC3( A, B )
      REAL               A, B
      EXTERNAL           SLAMC3
      REAL               SLAMC3
      WSLAMC3 = SLAMC3( A, B )
      END FUNCTION
