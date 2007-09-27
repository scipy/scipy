C*******************************************************************
C RETURN THE HISTOGRAM OF ARRAY X, THAT IS, THE NUMBER OF ELEMENTS
C IN X FALLING INTO EACH BIN.
C THE BIN ARRAY CONSISTS IN N BINS STARTING AT BIN0 WITH WIDTH DELTA.
C HISTO H : | LOWER OUTLIERS | 1 | 2 | 3 | ... |  N  | UPPER OUTLIERS |
C INDEX i : |        1       | 2 | 3 | 4 | ... | N+1 |      N+2       |

      SUBROUTINE FIXED_BINSIZE(X, BIN0, DELTA, N, NX, H)

C PARAMETERS
C ----------
C X : ARRAY 
C BIN0 : LEFT BIN EDGE
C DELTA : BIN WIDTH
C N : NUMBER OF BINS
C H : HISTOGRAM

      IMPLICIT NONE
      INTEGER :: N, NX, i, K
      DOUBLE PRECISION ::  X(NX), BIN0, DELTA
      INTEGER :: H(N+2), UP, LOW

CF2PY INTEGER INTENT(IN) :: N
CF2PY INTEGER INTENT(HIDE) :: NX = LEN(X)
CF2PY DOUBLE PRECISION DIMENSION(NX), INTENT(IN) :: X
CF2PY DOUBLE PRECISION INTENT(IN) :: BIN0, DELTA
CF2PY INTEGER DIMENSION(N+2), INTENT(OUT) :: H


      DO i=1,N+2
        H(i) = 0
      ENDDO
      
C     OUTLIERS INDICES
      UP = N+2
      LOW = 1

      DO i=1,NX
        IF (X(i) >= BIN0) THEN
          K = INT((X(i)-BIN0)/DELTA)+1
          IF (K <= N) THEN
            H(K+1) = H(K+1) + 1
          ELSE 
            H(UP) = H(UP) + 1
          ENDIF
        ELSE 
          H(LOW) = H(LOW) + 1
        ENDIF
      ENDDO

      END SUBROUTINE



C*******************************************************************
C RETURN THE WEIGHTED HISTOGRAM OF ARRAY X, THAT IS, THE SUM OF THE 
C WEIGHTS OF THE ELEMENTS OF X FALLING INTO EACH BIN.
C THE BIN ARRAY CONSISTS IN N BINS STARTING AT BIN0 WITH WIDTH DELTA.
C HISTO H : | LOWER OUTLIERS | 1 | 2 | 3 | ... |  N  | UPPER OUTLIERS |
C INDEX i : |        1       | 2 | 3 | 4 | ... | N+1 |      N+2       |

      SUBROUTINE WEIGHTED_FIXED_BINSIZE(X, W, BIN0, DELTA, N, NX, H)

C PARAMETERS
C ----------
C X : ARRAY 
C W : WEIGHTS
C BIN0 : LEFT BIN EDGE
C DELTA : BIN WIDTH
C N : NUMBER OF BINS
C H : HISTOGRAM

      IMPLICIT NONE
      INTEGER :: N, NX, i, K
      DOUBLE PRECISION ::  X(NX), W(NX), BIN0, DELTA, H(N+2)
      INTEGER :: UP, LOW

CF2PY INTEGER INTENT(IN) :: N
CF2PY INTEGER INTENT(HIDE) :: NX = LEN(X)
CF2PY DOUBLE PRECISION DIMENSION(NX), INTENT(IN) :: X, W
CF2PY DOUBLE PRECISION INTENT(IN) :: BIN0, DELTA
CF2PY DOUBLE PRECISION DIMENSION(N+2), INTENT(OUT) :: H


      DO i=1,N+2
        H(i) = 0.D0
      ENDDO
      
C     OUTLIERS INDICES
      UP = N+2
      LOW = 1

      DO i=1,NX
        IF (X(i) >= BIN0) THEN
          K = INT((X(i)-BIN0)/DELTA)+1
          IF (K <= N) THEN
            H(K+1) = H(K+1) + W(i)
          ELSE 
            H(UP) = H(UP) + W(i)
          ENDIF
        ELSE 
          H(LOW) = H(LOW) + W(i)
        ENDIF
      ENDDO

      END SUBROUTINE


C*****************************************************************************
C COMPUTE N DIMENSIONAL FLATTENED HISTOGRAM

      SUBROUTINE FIXED_BINSIZE_ND(X, BIN0, DELTA, N, COUNT, NX,D,NC)

C PARAMETERS
C ----------
C X : ARRAY (NXD)
C BIN0 : LEFT BIN EDGES (D)      
C DELTA : BIN WIDTH (D)
C N : NUMBER OF BINS (D)
C COUNT : FLATTENED HISTOGRAM (NC)
C NC : PROD(N(:)+2)

      IMPLICIT NONE
      INTEGER :: NX, D, NC,N(D), i, j, k, T
      DOUBLE PRECISION :: X(NX,D), BIN0(D), DELTA(D)
      INTEGER :: INDEX(NX), ORDER(D), MULT, COUNT(NC)


CF2PY DOUBLE PRECISION DIMENSION(NX,D), INTENT(IN) :: X
CF2PY DOUBLE PRECISION DIMENSION(D) :: BIN0, DELTA
CF2PY INTEGER INTENT(IN) :: N
CF2PY INTEGER DIMENSION(NC), INTENT(OUT) :: COUNT
CF2PY INTEGER INTENT(HIDE) :: NX=SHAPE(X,1)
CF2PY INTEGER INTENT(HIDE) :: D=SHAPE(X,2)


C     INITIALIZE INDEX
      DO i=1, NX
        INDEX(i) = 0
      ENDDO

C     INITIALIZE COUNT
      DO i=1,NC
        COUNT(i) = 0
      ENDDO

C     ORDER THE BIN SIZE ARRAY N(D)
      CALL QSORTI(ORDER, D, N)

C     INITIALIZE THE DIMENSIONAL MULTIPLIER
      MULT=1

C     FIND THE FLATTENED INDEX OF EACH SAMPLE
      DO j=1, D
        k = ORDER(j)
        MULT=MULT*N(k)

        DO i=1, NX 
          IF (X(i,k) >= BIN0(k)) THEN
            T = INT((X(i, k)-BIN0(k))/DELTA(k))+1
            IF (T <= N(k)) THEN
              T = T+1
            ELSE
              T = N(k)+2
            ENDIF
          ELSE
            T = 1
          ENDIF

          INDEX(i) = INDEX(I) + T*MULT
        ENDDO
      ENDDO

C     COUNT THE NUMBER OF SAMPLES FALLING INTO EACH BIN
      DO i=1,NX
        COUNT(INDEX(i)) =  COUNT(INDEX(i)) + 1
      ENDDO

      END SUBROUTINE 


C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C 
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.


      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      INTEGER X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      INTEGER A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) RETURN
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END
