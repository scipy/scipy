        SUBROUTINE UMD2P1 (WHO, WHERE,
     $          N, NE, JOB, TRANS, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, LX, W, LW)
        INTEGER WHO, WHERE, N, NE, JOB, LVALUE, LINDEX, INDEX (LINDEX),
     $          KEEP (20), ICNTL (20), INFO (40), LX, LW
        DOUBLE PRECISION
     $          VALUE (LVALUE), B (LX), X (LX), W (LW)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        LOGICAL TRANS
        
C=== UMD2P1 ============================================================
C
C  Unsymmetric-pattern MultiFrontal Package (UMFPACK). Version 2.2d
C  Copyright (C) 1997, Timothy A. Davis, University of Florida, USA.
C  ALL RIGHTS RESERVED.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  July 7, 1997. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C***********************************************************************
C* NOTICE:  "The UMFPACK Package may be used SOLELY for educational,   *
C* research, and benchmarking purposes by non-profit organizations and *
C* the U.S. government.  Commercial and other organizations may make   *
C* use of UMFPACK SOLELY for benchmarking purposes only.  UMFPACK may  *
C* be modified by or on behalf of the User for such use but at no time *
C* shall UMFPACK or any such modified version of UMFPACK become the    *
C* property of the User.  UMFPACK is provided without warranty of any  *
C* kind, either expressed or implied.  Neither the Authors nor their   *
C* employers shall be liable for any direct or consequential loss or   *
C* damage whatsoever arising out of the use or misuse of UMFPACK by    *
C* the User.  UMFPACK must not be sold.  You may make copies of        *
C* UMFPACK, but this NOTICE and the Copyright notice must appear in    *
C* all copies.  Any other use of UMFPACK requires written permission.  *
C* Your use of UMFPACK is an implicit agreement to these conditions."  *
C*                                                                     *
C* The MA38 Package in Release 12 of the Harwell Subroutine Library    *
C* (HSL) has equivalent functionality (and identical calling interface)*
C* as UMFPACK (the HSL has single and double precision versions only,  *
C* however).  It is available for commercial use.   Technical reports, *
C* information on HSL, and matrices are available via the World Wide   *
C* Web at http://www.cis.rl.ac.uk/struct/ARCD/NUM.html, or by          *
C* anonymous ftp at seamus.cc.rl.ac.uk/pub.  Also contact Dr. Scott    *
C* Roberts, Harwell Subroutine Library, B 552, AEA Technology,         *
C* Harwell, Didcot, Oxon OX11 0RA, England.                            *
C* telephone (44) 1235 434988, fax (44) 1235 434136                    *
C* email Scott.Roberts@aeat.co.uk, who will provide details of price   *
C* and conditions of use.                                              *
C***********************************************************************

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  print input/output arguments for UMD2FA, UMD2RF, and UMD2SO

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be deleted on installation (replaced with a dummy
C  routine that just returns without printing) in order to completely
C  disable the printing of all input/output parameters.  To completely
C  disable all I/O, you can also replace the UMD2P2 routine with a
C  dummy subroutine.  If you make this modification, please do
C  not delete any original code - just comment it out instead.  Add a
C  comment and date to your modifications.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who:            what routine called UMD2P1:
C                       1: UMD2FA, 2: UMD2RF, 3: UMD2SO
C       where:          called from where:
C                       1: entry of routine, else exit of routine
C       Icntl (3):      if < 3 then print nothing, if 3 then print
C                       terse output, if 4 print info and matrix
C                       values.  If 5, print everything.
C       Icntl (2):      I/O unit on which to print.  No printing
C                       occurs if < 0.
C
C       Parameters to print, see UMD2FA, UMD2RF, or UMD2SO for
C       descriptions:
C
C           n, ne, job, trans, lvalue, lindex, Value, Index, Keep,
C           Icntl, Info, Rinfo, B, X, lx, W, lw

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       on Icntl (2) I/O unit only

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  UMD2FA, UMD2RF, UMD2SO
C       functions called:       MIN
        INTRINSIC MIN

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        LOGICAL TRANSA, TRANSC, PRLU, BADLU, SGLTON, PRESRV, SYMBOL
        INTEGER IO, PRL, PRN, K, LUI1, LUI2, LUX1, LUX2, ROW, COL,
     $          FACNE, FACN, NZ, FACJOB, NBLKS, NZOFF, FACTRA, CPERMP,
     $          RPERMP, APP, AXP, AIP, OFFIP, OFFXP, LUBLPP, OFFPP,
     $          BLKPP, P1, P2, P, BLK, K1, K2, KN, LUIIP, LUXXP, NPIV,
     $          NLU, E, LUK, LUPP, LUIP, LUXP, LUDEGR, LUDEGC, LUNSON,
     $          LUSONP, LUCP, LURP, I, J, NZCOL, NZROW, UXP, SON,
     $          PRMAX, LUDIMR, LUDIMC, MAXDR, MAXDC, LUIR1, IP1, IP2,
     $          XP1
        DOUBLE PRECISION
     $          ONE
        PARAMETER (PRMAX = 10)

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C  prn:     number of entries printed so far
C  prmax:   maximum number of entries to print if prl = 3
C  prlu:    true if printing LU factors
C
C  Location and status of LU factors:
C  ----------------------------------
C  transc:  TRANSC argument in UMD2SO
C  transa:  TRANSA argument in UMD2FA or UMD2RF when matrix factorized
C  badlu:   true if LU factors uncomputed or corrupted
C  presrv:  true if original matrix was preserved when factorized
C  symbol:  true if only symbolic part of LU factors needed on input
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for a call to UMD2RF
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors start in Value (lux1...)
C  lux2:    real part of LU factors end in Value (...lux1)
C  ip1:     pointer into leading part of LU factors in Index
C  ip2:     pointer into trailing part of LU factors in Index
C  xp1:     pointer into leading part of LU factors in Value
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  app:     Ap (1..n+1) array located in Index (app...app+n)
C  axp:     Ax (1..nz) array located in Value (axp...axp+nz-1)
C  aip:     Ai (1..nz) array located in Index (aip...aip+nz-1)
C  offip:   Offi (1..nzoff) array loc. in Index (offip...offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array loc. in Value (offxp...offxp+nzoff-1)
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  ...      seven scalars in Index (lui2-6...lui2):
C  factra:  0/1 if TRANSA argument was false/true in UMD2FA or UMD2RF
C  nzoff:   number of entries in off-diagonal part
C  nblks:   number of diagonal blocks
C  facjob:  JOB argument in UMD2FA or UMD2RF when matrix factorized 
C  nz:      entries in A
C  facn:    N argument in UMD2FA or UMD2RF when matrix factorized 
C  facne:   NE argument in UMD2FA or UMD2RF when matrix factorized 
C
C  A single diagonal block and its LU factors:
C  -------------------------------------------
C  blk:     current diagonal block
C  k1,k2:   current diagonal is A (k1..k2, k1..k2)
C  kn:      order of current diagonal block (= k2-k1+1)
C  sglton:  true if current diagonal block is 1-by-1 (a singleton)
C  luiip:   LU factors of a diagonal block start in Index (luiip...)
C  luxxp:   LU factors of a diagonal block start in Value (luxxp...)
C  npiv:    number of pivots in a diagonal block (0 <= npiv <= kn)
C  nlu:     number of elements in a diagonal block
C  lupp:    LUp (1..nlu) array located in Index (lupp...lupp+nlu-1)
C
C  An element in the LU factors of a single diagonal block:
C  --------------------------------------------------------
C  e:       element
C  luk:     number of pivots in element e
C  luip:    integer part of element is in Index (luip...)
C  luxp:    real part of element e is in Value (luxp...)
C  ludegr:  row degree (number of columns) of U2 block in element e
C  ludegc:  column degree (number of rows) of L2 block in element e
C  lunson:  number of sons of element e in the assembly DAG
C  lusonp:  list of sons of element e in Index(lusonp...lusonp+lunson-1)
C  lucp:    column pattern (row indices) of L2 block in Index (lucp..)
C  lurp:    row pattern (column indices) of U2 block in Index (lurp..)
C  nzcol:   entries in a column of L, including unit diagonal
C  nzrow:   entries in a row of U, including non-unit diagonal
C  uxp:     a row of the U2 block located in Value (uxp...)
C  son:     a son of the element e
C  ludimr:  row dimension (number of columns) in frontal matrix
C  ludimc:  column dimension (number of rows) in frontal matrix
C  maxdr:   largest ludimr for this block
C  maxdc:   largest ludimc for this block
C
C  Other:
C  ------
C  row:     row index
C  col:     column index
C  k:       kth pivot, and general loop index
C  i, j:    loop indices
C  p:       pointer
C  p1:      column of A starts Ai/Ax (p1...), or row Offi/x (p1...)
C  p2:      column of A ends in Ai/Ax (...p2), or row Offi/x (...p2)

C=======================================================================
C  EXECUTABLE STATEMENTS:
C       if (printing disabled on installation) return
C=======================================================================

C-----------------------------------------------------------------------
C  get printing control parameters
C-----------------------------------------------------------------------

        ONE = 1
        IO = ICNTL (2)
        PRL = ICNTL (3)
        IF (PRL .LT. 3 .OR. IO .LT. 0) THEN 
C          printing has not been requested
           RETURN
        ENDIF 

C-----------------------------------------------------------------------
C  who is this, and where.  Determine if LU factors are to be printed
C-----------------------------------------------------------------------

        IF (WHO .EQ. 1) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 6) 'UMD2FA input:'
              PRLU = .FALSE.
           ELSE 
              WRITE (IO, 6) 'UMD2FA output:'
              PRLU = .TRUE.
           ENDIF 
        ELSE IF (WHO .EQ. 2) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 6) 'UMD2RF input:'
              PRLU = .TRUE.
           ELSE 
              WRITE (IO, 6) 'UMD2RF output:'
              PRLU = .TRUE.
           ENDIF 
        ELSE IF (WHO .EQ. 3) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 6) 'UMD2SO input:'
              PRLU = .TRUE.
           ELSE 
              WRITE (IO, 6) 'UMD2SO output:'
              PRLU = .FALSE.
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print scalar input arguments: n, ne, job, trans, lvalue, lindex
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1) THEN 
           WRITE (IO, 1)  'Scalar arguments:'
           WRITE (IO, 1)  '   N:         ', N, ' : order of matrix A'
           IF (WHO .EQ. 3) THEN 
C             UMD2SO:
              LUI2 = KEEP (5)

C             was A or A^T factorized?
              TRANSA = .FALSE.
              IF (LUI2-6 .GE. 1 .AND. LUI2-6 .LE. LINDEX) THEN 
                 TRANSA = INDEX (LUI2-6) .NE. 0
              ENDIF 
              TRANSC = TRANS
              IF (.NOT. TRANSC) THEN 
                 IF (JOB .EQ. 1) THEN 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve P''Lx=b'
                 ELSE IF (JOB .EQ. 2) THEN 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve UQ''x=b'
                 ELSE IF (.NOT. TRANSA) THEN 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve Ax=b (PAQ=LU was factorized)'
                 ELSE 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve A''x=b (PA''Q=LU was factorized)'
                 ENDIF 
              ELSE 
                 IF (JOB .EQ. 1) THEN 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve L''Px=b'
                 ELSE IF (JOB .EQ. 2) THEN 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve QU''x=b'
                 ELSE IF (.NOT. TRANSA) THEN 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve A''x=b (PAQ=LU was factorized)'
                 ELSE 
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve Ax=b (PA''Q=LU was factorized)'
                 ENDIF 
              ENDIF 
              IF (TRANSC) THEN 
                 WRITE (IO, 1)
     $           '   TRANSC:          .true. : see JOB above '
              ELSE 
                 WRITE (IO, 1)
     $           '   TRANSC:         .false. : see JOB above '
              ENDIF 

           ELSE 
C             UMD2FA or UMD2RF:
              WRITE (IO, 1) '   NE:        ', NE,
     $        ' : entries in matrix A'
              IF (JOB .EQ. 1) THEN 
                 WRITE (IO, 1) '   JOB:       ', JOB,
     $           ' : matrix A preserved'
              ELSE 
                 WRITE (IO, 1) '   JOB:       ', JOB,
     $           ' : matrix A not preserved'
              ENDIF 

              TRANSA = TRANS
              IF (TRANSA) THEN 
                 WRITE (IO, 1)
     $           '   TRANSA:          .true. : factorize A transpose'
              ELSE 
                 WRITE (IO, 1)
     $           '   TRANSA:         .false. : factorize A'
              ENDIF 

           ENDIF 
           WRITE (IO, 1) '   LVALUE:    ',LVALUE,
     $     ' : size of VALUE array'
           WRITE (IO, 1) '   LINDEX:    ',LINDEX,
     $     ' : size of INDEX array'
        ENDIF 

C-----------------------------------------------------------------------
C  print control parameters:  Icntl, Cntl, and Keep (6..8)
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1) THEN 
           WRITE (IO, 1)
     $     'Control parameters, normally initialized by UMD21I:'
           WRITE (IO, 1) '   ICNTL (1): ', ICNTL (1),
     $     ' : I/O unit for error and warning messages'
           WRITE (IO, 1) '   ICNTL (2): ', IO,
     $     ' : I/O unit for diagnostics'
           WRITE (IO, 1) '   ICNTL (3): ', PRL,
     $     ' : printing control'
           IF (WHO .EQ. 1) THEN 
              IF (ICNTL (4) .EQ. 1) THEN 
                 WRITE (IO, 1) '   ICNTL (4): ', ICNTL (4),
     $           ' : use block triangular form (BTF)'
              ELSE 
                 WRITE (IO, 1) '   ICNTL (4): ', ICNTL (4),
     $           ' : do not permute to block triangular form (BTF)'
              ENDIF 
              WRITE (IO, 1) '   ICNTL (5): ', ICNTL (5),
     $        ' : columns examined during pivot search'
              IF (ICNTL (6) .NE. 0) THEN 
                 WRITE (IO, 1) '   ICNTL (6): ', ICNTL (6),
     $           ' : preserve symmetry'
              ELSE 
                 WRITE (IO, 1) '   ICNTL (6): ', ICNTL (6),
     $           ' : do not preserve symmetry'
              ENDIF 
           ENDIF 
           IF (WHO .NE. 3) THEN 
              WRITE (IO, 1) '   ICNTL (7): ', ICNTL (7),
     $        ' : block size for dense matrix multiply'
           ELSE 
              WRITE (IO, 1) '   ICNTL (8): ', ICNTL (8),
     $        ' : maximum number of iterative refinement steps'
           ENDIF 
           IF (WHO .EQ. 1) THEN 
              WRITE (IO, 3) '   CNTL (1):   ',CNTL (1),
     $        ' : relative pivot tolerance'
              WRITE (IO, 3) '   CNTL (2):   ',CNTL (2),
     $        ' : frontal matrix growth factor'
              WRITE (IO, 1) '   KEEP (6):  ',KEEP(6),
     $        ' : largest positive integer'
              WRITE (IO, 1) '   KEEP (7):  ',KEEP(7),
     $        ' : dense row/col control, d1'
              WRITE (IO, 1) '   KEEP (8):  ',KEEP(8),
     $        ' : dense row/col control, d2'
           ELSE IF (WHO .EQ. 3) THEN 
              WRITE (IO, 3) '   CNTL (3):   ',CNTL(3),
     $        ' : machine epsilon'
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print the informational output
C-----------------------------------------------------------------------

        IF (WHERE .NE. 1) THEN 
           WRITE (IO, 1) 'Output information:'
           IF (INFO (1) .LT. 0) THEN 
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : error occurred!'
           ELSE IF (INFO (1) .GT. 0) THEN 
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : warning occurred'
           ELSE 
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : no error or warning occurred'
           ENDIF 
           IF (WHO .NE. 3) THEN 
              WRITE (IO, 1) '   INFO (2):  ', INFO (2),
     $        ' : duplicate entries in A'
              WRITE (IO, 1) '   INFO (3):  ', INFO (3),
     $        ' : invalid entries in A (indices not in 1..N)'
              WRITE (IO, 1) '   INFO (4):  ', INFO (4),
     $        ' : invalid entries in A (not in prior pattern)'
              WRITE (IO, 1) '   INFO (5):  ', INFO (5),
     $        ' : entries in A after summing duplicates'
              WRITE (IO, 1)
     $  '                             and removing invalid entries'
              WRITE (IO, 1) '   INFO (6):  ', INFO (6),
     $        ' : entries in diagonal blocks of A'
              WRITE (IO, 1) '   INFO (7):  ', INFO (7),
     $        ' : entries in off-diagonal blocks of A'
              WRITE (IO, 1) '   INFO (8):  ', INFO (8),
     $        ' : 1-by-1 diagonal blocks in A'
              WRITE (IO, 1) '   INFO (9):  ', INFO (9),
     $        ' : diagonal blocks in A (>1 only if BTF used)'
              WRITE (IO, 1) '   INFO (10): ', INFO (10),
     $        ' : entries below diagonal in L'
              WRITE (IO, 1) '   INFO (11): ', INFO (11),
     $        ' : entries above diagonal in U'
              WRITE (IO, 1) '   INFO (12): ', INFO (12),
     $        ' : entries in L + U + offdiagonal blocks of A'
              WRITE (IO, 1) '   INFO (13): ', INFO (13),
     $        ' : frontal matrices'
              WRITE (IO, 1) '   INFO (14): ', INFO (14),
     $        ' : integer garbage collections'
              WRITE (IO, 1) '   INFO (15): ', INFO (15),
     $        ' : real garbage collections'
              WRITE (IO, 1) '   INFO (16): ', INFO (16),
     $        ' : diagonal pivots chosen'
              WRITE (IO, 1) '   INFO (17): ', INFO (17),
     $        ' : numerically valid pivots found in A'
              WRITE (IO, 1) '   INFO (18): ', INFO (18),
     $        ' : memory used in INDEX'
              WRITE (IO, 1) '   INFO (19): ', INFO (19),
     $        ' : minimum memory needed in INDEX'
              WRITE (IO, 1) '   INFO (20): ', INFO (20),
     $        ' : memory used in VALUE'
              WRITE (IO, 1) '   INFO (21): ', INFO (21),
     $        ' : minimum memory needed in VALUE'
              WRITE (IO, 1) '   INFO (22): ', INFO (22),
     $        ' : memory needed in INDEX for next call to UMD2RF'
              WRITE (IO, 1) '   INFO (23): ', INFO (23),
     $        ' : memory needed in VALUE for next call to UMD2RF'
           ELSE 
              WRITE (IO, 1) '   INFO (24): ', INFO (24),
     $        ' : steps of iterative refinement taken'
           ENDIF 
           IF (WHO .NE. 3) THEN 
              WRITE (IO, 3) '   RINFO (1):  ', RINFO (1),
     $        ' : total BLAS flop count'
              WRITE (IO, 3) '   RINFO (2):  ', RINFO (2),
     $        ' : assembly flop count'
              WRITE (IO, 3) '   RINFO (3):  ', RINFO (3),
     $        ' : pivot search flop count'
              WRITE (IO, 3) '   RINFO (4):  ', RINFO (4),
     $        ' : Level-1 BLAS flop count'
              WRITE (IO, 3) '   RINFO (5):  ', RINFO (5),
     $        ' : Level-2 BLAS flop count'
              WRITE (IO, 3) '   RINFO (6):  ', RINFO (6),
     $        ' : Level-3 BLAS flop count'
           ELSE IF (LW .EQ. 4*N) THEN 
              WRITE (IO, 3) '   RINFO (7):  ', RINFO (7),
     $        ' : sparse error estimate omega1'
              WRITE (IO, 3) '   RINFO (8):  ', RINFO (8),
     $        ' : sparse error estimate omega2'
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print input matrix A, in triplet form, for UMD2FA and UMD2RF
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1 .AND. WHO .NE. 3) THEN 

           IF (TRANSA) THEN 
              IF (PRL .GE. 5) THEN 
                 WRITE (IO, 1) 'The input matrix A transpose:'
                 WRITE (IO, 1) '   VALUE (1 ... ',NE,
     $           ' ): numerical values'
                 WRITE (IO, 1) '   INDEX (1 ... ',NE,
     $           ' ): column indices'
                 WRITE (IO, 1) '   INDEX (',NE+1,' ... ',2*NE,
     $           ' ): row indices'
              ENDIF 
              WRITE (IO, 1)
     $        'Input matrix A transpose (entry: row, column, value):'
           ELSE 

              IF (PRL .GE. 5) THEN 
                 WRITE (IO, 1) 'The input matrix A:'
                 WRITE (IO, 1) '   VALUE (1 ... ',NE,
     $           ' ): numerical values'
                 WRITE (IO, 1) '   INDEX (1 ... ',NE,
     $           ' ): row indices'
                 WRITE (IO, 1) '   INDEX (',NE+1,' ... ',2*NE,
     $           ' ): column indices'
              ENDIF 
              WRITE (IO, 1)
     $        'Input matrix A (entry: row, column, value):'

           ENDIF 

           PRN = MIN (PRMAX, NE)
           IF (PRL .GE. 4) THEN 
              PRN = NE
           ENDIF 
           DO 20 K = 1, PRN 

              IF (TRANSA) THEN 
                 ROW = INDEX (K+NE)
                 COL = INDEX (K)
              ELSE 
                 ROW = INDEX (K)
                 COL = INDEX (K+NE)
              ENDIF 

              WRITE (IO, 2) K, ROW, COL, VALUE (K)
20         CONTINUE 
           IF (PRN .LT. NE) THEN 
              WRITE (IO, 7)
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print the LU factors:  UMD2FA output, UMD2RF input/output,
C                         and UMD2SO input
C-----------------------------------------------------------------------

        IF (PRLU .AND. INFO (1) .LT. 0) THEN 
           WRITE (IO, 1)
     $     'LU factors not printed because of error flag, INFO (1) ='
     $     , INFO (1)
           PRLU = .FALSE.
        ENDIF 

        IF (PRLU) THEN 

C          -------------------------------------------------------------
C          description of what must be preserved between calls
C          -------------------------------------------------------------

           LUX1 = KEEP (1)
           LUX2 = KEEP (2)
           LUI1 = KEEP (3)
           LUIR1 = KEEP (4)
           LUI2 = KEEP (5)

           XP1 = LUX1
           IP1 = LUI1
           IP2 = LUI2

C          -------------------------------------------------------------
C          on input to UMD2RF, only the symbol information is used
C          -------------------------------------------------------------

           SYMBOL = WHO .EQ. 2 .AND. WHERE .EQ. 1

           IF (PRL .GE. 5) THEN 
              IF (SYMBOL) THEN 
                 WRITE (IO, 1)
     $           'KEEP (4...5) gives the location of LU factors'
                 WRITE (IO, 1)
     $           '   which must be preserved for calls to UMD2RF: '
              ELSE 
                 WRITE (IO, 1)
     $           'KEEP (1...5) gives the location of LU factors'
                 WRITE (IO, 1)
     $           '   which must be preserved for calls to UMD2SO: '
                 WRITE (IO, 1) '      VALUE ( KEEP (1): ', LUX1,
     $           ' ... KEEP (2): ', LUX2,' )'
                 WRITE (IO, 1) '      INDEX ( KEEP (3): ', LUI1,
     $           ' ... KEEP (5): ', LUI2,' )'
                 WRITE (IO, 1) '   and for calls to UMD2RF: '
              ENDIF 
              WRITE (IO, 1) '      INDEX ( KEEP (4): ',LUIR1,
     $        ' ... KEEP (5): ', LUI2,' )'
           ENDIF 

           BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR.
     $        LUI2 .GT. LINDEX
           IF (.NOT. SYMBOL) THEN 
              BADLU = BADLU .OR. LUX1 .LE. 0 .OR.
     $        LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE .OR. LUI1 .LE. 0 .OR.
     $        LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
           ENDIF 

C          -------------------------------------------------------------
C          get the 7 scalars, and location of permutation vectors
C          -------------------------------------------------------------

           IF (BADLU) THEN 
C             pointers are bad, so these values cannot be obtained
              FACNE  = 0
              FACN   = 0
              NZ     = 0
              FACJOB = 0
              NBLKS  = 0
              NZOFF  = 0
              FACTRA = 0
           ELSE 
              FACNE  = INDEX (LUI2)
              FACN   = INDEX (LUI2-1)
              NZ     = INDEX (LUI2-2)
              FACJOB = INDEX (LUI2-3)
              NBLKS  = INDEX (LUI2-4)
              NZOFF  = INDEX (LUI2-5)
              FACTRA = INDEX (LUI2-6)
           ENDIF 

           PRESRV = FACJOB .NE. 0
           TRANSA = FACTRA .NE. 0
           RPERMP = (LUI2-6) - (FACN)
           CPERMP = RPERMP - (FACN)
           IP2 = CPERMP - 1

           IF (PRL .GE. 5) THEN 
              IF (SYMBOL) THEN 
                 WRITE (IO, 1) 'Layout of LU factors in INDEX:'
              ELSE 
                 WRITE (IO, 1)
     $           'Layout of LU factors in VALUE and INDEX:'
              ENDIF 
           ENDIF 

C          -------------------------------------------------------------
C          get location of preserved input matrix
C          -------------------------------------------------------------

           IF (PRESRV) THEN 
C             preserved column-form of original matrix
              APP = IP1
              AIP = APP + (FACN+1)
              IP1 = AIP + (NZ)
              AXP = XP1
              XP1 = XP1 + (NZ)
              IF (PRL .GE. 5 .AND. .NOT. SYMBOL) THEN 
                 WRITE (IO, 1)'   preserved copy of original matrix:'
                 WRITE (IO, 1)'      INDEX ( ',APP,' ... ', AIP-1,
     $           ' ): column pointers'
                 WRITE (IO, 1)'      INDEX ( ',AIP,' ... ', IP1-1,
     $           ' ): row indices'
                 WRITE (IO, 1)'      VALUE ( ',AXP,' ... ', XP1-1,
     $           ' ): numerical values'
              ENDIF 
           ELSE 
              IF (PRL .GE. 5 .AND. .NOT. SYMBOL) THEN 
                 WRITE (IO, 1) '   original matrix not preserved.'
              ENDIF 
           ENDIF 

           BADLU = BADLU .OR.
     $          N .NE. FACN .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
           IF (.NOT. SYMBOL) THEN 
              BADLU = BADLU .OR. XP1 .GT. LUX2 .OR. NZOFF .LT. 0
           ENDIF 
           IF (BADLU) THEN 
              NBLKS = 0
           ENDIF 

           IF (NBLKS .LE. 1) THEN 

C             ----------------------------------------------------------
C             single block (or block triangular form not used),
C             or LU factors are corrupted
C             ----------------------------------------------------------

              IF (PRL .GE. 5) THEN 
                 WRITE (IO, 1)
     $           '   collection of elements in LU factors:'
                 WRITE (IO, 1) '      INDEX ( ',LUIR1,' ... ', IP2,
     $           ' ): integer data'
                 IF (.NOT. SYMBOL) THEN 
                    WRITE (IO, 1) '      VALUE ( ',XP1,' ... ', LUX2,
     $              ' ): numerical values'
                 ENDIF 
              ENDIF 

           ELSE 

C             ----------------------------------------------------------
C             block triangular form with more than one block
C             ----------------------------------------------------------

              OFFIP = IP1
              IP1 = IP1 + (NZOFF)
              OFFXP = XP1
              XP1 = XP1 + (NZOFF)
              OFFPP = CPERMP - (N+1)
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              IP2 = LUBLPP - 1
              BADLU = BADLU .OR. LUIR1 .GT. IP2
              IF (.NOT. SYMBOL) THEN 
                 BADLU = BADLU .OR. IP1 .GT. IP2 .OR.
     $           XP1 .GT. LUX2 .OR. LUIR1 .NE. IP1
              ENDIF 

              IF (PRL .GE. 5) THEN 
                 WRITE (IO, 1)
     $           '   matrix permuted to upper block triangular form.'
                 IF (NZOFF .NE. 0 .AND. .NOT. SYMBOL) THEN 
                    WRITE (IO, 1)'   entries not in diagonal blocks:'
                    WRITE (IO, 1)'      INDEX ( ',OFFIP,' ... ',
     $              LUIR1-1, ' ): row indices'
                    WRITE (IO, 1)'      VALUE ( ',OFFXP,' ... ',
     $              XP1-1, ' ): numerical values'
                 ENDIF 
                 WRITE (IO, 1)
     $  '   collection of elements in LU factors of diagonal blocks:'
                 IF (LUIR1 .LE. LUBLPP-1) THEN 
                    WRITE (IO, 1) '      INDEX ( ',LUIR1,' ... ',
     $              IP2, ' ): integer data'
                 ENDIF 
                 IF (XP1 .LE. LUX2 .AND. .NOT. SYMBOL) THEN 
                    WRITE (IO, 1) '      VALUE ( ',XP1,' ... ', LUX2,
     $              ' ): numerical values'
                 ENDIF 
                 WRITE (IO, 1) '   other block triangular data:'
                 WRITE (IO, 1) '      INDEX ( ',LUBLPP,' ... ',
     $           BLKPP-1, ' ): pointers to block factors' 
                 WRITE (IO, 1) '      INDEX ( ', BLKPP,' ... ',
     $           OFFPP-1, ' ): index range of blocks'
                 IF (.NOT. SYMBOL) THEN 
                    WRITE (IO, 1) '      INDEX ( ', OFFPP,' ... ',
     $              LUI2-7,' ): off-diagonal row pointers'
                 ENDIF 
              ENDIF 

           ENDIF 

C          -------------------------------------------------------------
C          print location of permutation vectors and 7 scalars at tail
C          -------------------------------------------------------------

           IF (PRL .GE. 5) THEN 
              WRITE (IO, 1)
     $        '   permutation vectors (start at KEEP(4)-2*N-6):'
              WRITE (IO, 1) '      INDEX ( ',CPERMP,' ... ',RPERMP-1,
     $        ' ): column permutations'
              WRITE (IO, 1) '      INDEX ( ',RPERMP,' ... ',LUI2-7,
     $        ' ): row permutations'
              WRITE (IO, 1) '   other data in INDEX: '
              WRITE (IO, 1) '      INDEX ( ',LUI2-6,' ): ', FACTRA,
     $        ' : TRANSA UMD2FA/UMD2RF argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2-5,' ): ', NZOFF,
     $        ' : entries in off-diagonal part'
              WRITE (IO, 1) '      INDEX ( ',LUI2-4,' ): ', NBLKS,
     $        ' : number of diagonal blocks'
              WRITE (IO, 1) '      INDEX ( ',LUI2-3,' ): ', FACJOB,
     $        ' : JOB UMD2FA/UMD2RF argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2-2,' ): ', NZ,
     $        ' : entries in original matrix'
              WRITE (IO, 1) '      INDEX ( ',LUI2-1,' ): ', FACN,
     $        ' : N UMD2FA/UMD2RF argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2  ,' ): ', FACNE,
     $        ' : NE UMD2FA/UMD2RF argument'
           ENDIF 

           IF (.NOT. SYMBOL) THEN 
              BADLU = BADLU .OR. IP1 .NE. LUIR1
           ENDIF 
           IP1 = LUIR1
           IF (BADLU) THEN 
              WRITE (IO, 1) 'LU factors uncomputed or corrupted!'
              PRESRV = .FALSE.
              NBLKS = 0
           ENDIF 

C          -------------------------------------------------------------
C          copy of original matrix in column-oriented form
C          -------------------------------------------------------------

           IF (PRESRV .AND. .NOT. SYMBOL) THEN 
              WRITE (IO, 8)
              WRITE (IO, 1)
     $        'Preserved copy of original matrix (stored by column),'
              WRITE (IO, 1) 'one entry per line (row index, value):'
              DO 40 COL = 1, N 
                 P1 = INDEX (APP-1 + COL)
                 P2 = INDEX (APP-1 + COL+1) - 1
                 WRITE (IO, 1) '   col: ', COL
                 IF (PRL .EQ. 3) THEN 
                    P2 = MIN (PRMAX, P2)
                 ENDIF 
                 DO 30 P = P1, P2 
                    WRITE (IO, 5) INDEX (AIP-1 + P), VALUE (AXP-1 + P)
30               CONTINUE 
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN 
C                   exit out of loop if done printing:
                    WRITE (IO, 7)
                    GO TO 50
                 ENDIF 
40            CONTINUE 
C             loop exit label:
50            CONTINUE
           ENDIF 

C          -------------------------------------------------------------
C          entries in off-diagonal blocks, in row-oriented form
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1 .AND. .NOT. SYMBOL) THEN 
              WRITE (IO, 8)
              WRITE (IO, 1)
     $        'Entries not in diagonal blocks (stored by row):'
              WRITE (IO, 1) 'one entry per line (column index, value):'
              IF (NZOFF .EQ. 0) THEN 
                 WRITE (IO, 1) '   (none)'
              ENDIF 
              DO 70 ROW = 1, N 
                 P1 = INDEX (OFFPP-1 + ROW)
                 P2 = INDEX (OFFPP-1 + ROW+1) - 1
                 IF (P2 .GE. P1) THEN 
                    WRITE (IO, 1) '   row: ', ROW
                    IF (PRL .EQ. 3) THEN 
                       P2 = MIN (PRMAX, P2)
                    ENDIF 
                    DO 60 P = P1, P2 
                       WRITE (IO, 5)
     $                 INDEX (OFFIP-1 + P), VALUE (OFFXP-1 + P)
60                  CONTINUE 
                 ENDIF 
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN 
C                   exit out of loop if done printing:
                    WRITE (IO, 7)
                    GO TO 80
                 ENDIF 
70            CONTINUE 
C             loop exit label:
80            CONTINUE
           ENDIF 

C          -------------------------------------------------------------
C          LU factors of each diagonal block
C          -------------------------------------------------------------

           WRITE (IO, 8)
           IF (NBLKS .GT. 0) THEN 
              IF (SYMBOL) THEN 
                 WRITE (IO, 1) 'Nonzero pattern of prior LU factors:'
              ELSE 
                 WRITE (IO, 1) 'LU factors:'
              ENDIF 
           ENDIF 
           PRN = 0
           DO 200 BLK = 1, NBLKS 

C             ----------------------------------------------------------
C             print the factors of a single diagonal block
C             ----------------------------------------------------------

              IF (NBLKS .GT. 1) THEN 
                 K1 = INDEX (BLKPP-1 + BLK)
                 K2 = INDEX (BLKPP-1 + BLK+1) - 1
                 KN = K2-K1+1
                 SGLTON = KN .EQ. 1
                 IF (SGLTON) THEN 
C                   this is a singleton
                    LUXXP = XP1-1 + INDEX (LUBLPP-1 + BLK)
                 ELSE 
                    LUIIP = IP1-1 + INDEX (LUBLPP-1 + BLK)
                 ENDIF 
                 IF (BLK .GT. 1) THEN 
                    WRITE (IO, 9)
                 ENDIF 
              ELSE 
                 SGLTON = .FALSE.
                 K1 = 1
                 K2 = N
                 KN = N
                 LUIIP = IP1
              ENDIF 

              IF (SGLTON) THEN 

C                -------------------------------------------------------
C                this is a singleton
C                -------------------------------------------------------

                 IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN 
C                   exit out of loop if done printing:
                    WRITE (IO, 7)
                    GO TO 210
                 ENDIF 
                 PRN = PRN + 1
                 IF (SYMBOL) THEN 
                    WRITE (IO, 1) 'Block: ', BLK,
     $              ' (singleton) at index : ', K1
                 ELSE 
                    WRITE (IO, 4) 'Block: ', BLK,
     $              ' (singleton) at index : ', K1,'       value: ',
     $              VALUE (LUXXP)
                 ENDIF 
                 IF (PRL .GE. 5) THEN 
                    WRITE (IO, 1) 'located in VALUE ( ', LUXXP,' )'
                 ENDIF 

              ELSE 

C                -------------------------------------------------------
C                this block is larger than 1-by-1
C                -------------------------------------------------------

                 LUXXP = XP1-1 + INDEX (LUIIP)
                 NLU = INDEX (LUIIP+1)
                 NPIV = INDEX (LUIIP+2)
                 MAXDC = INDEX (LUIIP+3)
                 MAXDR = INDEX (LUIIP+4)
                 LUPP = LUIIP+5
                 IF (NBLKS .GT. 1) THEN 
                    WRITE (IO, 1) 'Block: ',BLK,' first index: ',K1,
     $              ' last index: ',K2
                 ENDIF 
                 IF (PRL .GE. 5) THEN 
                    WRITE (IO, 1) 'elements: ', NLU, ' pivots: ', NPIV
                    WRITE (IO, 1) 'largest contribution block: ',
     $                         MAXDC, ' by ', MAXDR
                    WRITE (IO, 1)'located in INDEX ( ',LUIIP,' ... )'
                    IF (.NOT. SYMBOL) THEN 
                       WRITE (IO, 1) 'and in VALUE ( ',LUXXP,' ... )'
                    ENDIF 
                 ENDIF 
                 LUIIP = LUPP + NLU

C                Note: the indices of the LU factors of the block range
C                from 1 to kn, even though the kn-by-kn block resides in
C                A (k1 ... k2, k1 ... k2).
                 K = 0

                 DO 190 E = 1, NLU 

C                   ----------------------------------------------------
C                   print a single element
C                   ----------------------------------------------------

                    LUIP = LUIIP-1 + INDEX (LUPP-1 + E)
                    LUXP = LUXXP-1 + INDEX (LUIP)
                    LUK  = INDEX (LUIP+1)
                    LUDEGR = INDEX (LUIP+2)
                    LUDEGC = INDEX (LUIP+3)
                    LUNSON = INDEX (LUIP+4)
                    LUDIMR = INDEX (LUIP+5)
                    LUDIMC = INDEX (LUIP+6)
                    LUCP = LUIP + 7
                    LURP = LUCP + LUDEGC
                    LUSONP = LURP + LUDEGR
                    IF (PRL .GE. 5) THEN 
                       WRITE (IO, 1) '   e: ', E, ' pivots: ', LUK
                       WRITE (IO, 1) '   children in dag: ', LUNSON,
     $                 ' frontal matrix: ', LUDIMR, ' by ', LUDIMC
                       ENDIF 

C                   ----------------------------------------------------
C                   print the columns of L
C                   ----------------------------------------------------

                    P = LUXP
                    DO 140 J = 1, LUK 
                       COL = K+J
                       NZCOL = LUK-J+1+LUDEGC
                       WRITE (IO, 1) '      L, col: ', COL
                       PRN = PRN + 1
                       ROW = COL
                       IF (SYMBOL) THEN 
                          WRITE (IO, 5) ROW
                       ELSE 
C                         L is unit diagonal:
                          WRITE (IO, 5) ROW, ONE
                       ENDIF 
                       P = P + 1
C                      pivot block
                       DO 120 I = J+1, LUK 
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN 
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF 
                          PRN = PRN + 1
                          ROW = K+I
                          IF (SYMBOL) THEN 
                             WRITE (IO, 5) ROW
                          ELSE 
                             WRITE (IO, 5) ROW, VALUE (P)
                          ENDIF 
                          P = P + 1
120                    CONTINUE 
C                      L block
                       DO 130 I = 1, LUDEGC 
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN  
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF 
                          PRN = PRN + 1
                          ROW = INDEX (LUCP-1+I)
                          IF (SYMBOL) THEN 
                             WRITE (IO, 5) ROW
                          ELSE 
                             WRITE (IO, 5) ROW, VALUE (P)
                          ENDIF 
                          P = P + 1
130                    CONTINUE 
                       P = P + J
140                 CONTINUE 

C                   ----------------------------------------------------
C                   print the rows of U
C                   ----------------------------------------------------

                    UXP = LUXP + LUK*(LUDEGC+LUK)
                    DO 170 I = 1, LUK 
                       ROW = K+I
                       NZROW = LUK-I+1+LUDEGR
                       WRITE (IO, 1) '      U, row: ', ROW
                       P = LUXP + (I-1) + (I-1) * (LUDEGC+LUK)
C                      pivot block
                       DO 150 J = I, LUK 
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN 
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF 
                          PRN = PRN + 1
                          COL = K+J
                          IF (SYMBOL) THEN 
                             WRITE (IO, 5) COL
                          ELSE 
                             WRITE (IO, 5) COL, VALUE (P)
                          ENDIF 
                          P = P + (LUDEGC+LUK)
150                    CONTINUE 
                       P = UXP
C                      U block
                       DO 160 J = 1, LUDEGR 
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN 
C                            exit out of loop if done printing:
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF 
                          PRN = PRN + 1
                          COL = INDEX (LURP-1+J)
                          IF (SYMBOL) THEN 
                             WRITE (IO, 5) COL
                          ELSE 
                             WRITE (IO, 5) COL, VALUE (P)
                          ENDIF 
                          P = P + LUK
160                    CONTINUE 
                       UXP = UXP + 1
170                 CONTINUE 

C                   ----------------------------------------------------
C                   print the sons of the element in the assembly DAG
C                   ----------------------------------------------------

                    IF (PRL .GE. 5) THEN 
                       DO 180 I = 1, LUNSON 
                          PRN = PRN + 1
                          SON = INDEX (LUSONP-1+I)
                          IF (SON .LE. KN) THEN 
C                            an LUson
                             WRITE (IO, 1) '      LUson: ', SON
                          ELSE IF (SON .LE. 2*KN) THEN 
C                            a Uson
                             WRITE (IO, 1) '      Uson:  ', SON-KN
                          ELSE 
C                            an Lson
                             WRITE (IO, 1) '      Lson:  ', SON-2*KN
                          ENDIF 
180                    CONTINUE 
                    ENDIF 

C                   ----------------------------------------------------
C                   increment count of pivots within this block
C                   ----------------------------------------------------

                    K = K + LUK
190              CONTINUE 
              ENDIF 
200        CONTINUE 
C          loop exit label:
210        CONTINUE

C          -------------------------------------------------------------
C          row and column permutations
C          -------------------------------------------------------------

           IF (.NOT. BADLU) THEN 
              PRN = MIN (PRMAX, N)
              IF (PRL .GE. 4) THEN 
C                print all of Cperm and Rperm
                 PRN = N
              ENDIF 
              WRITE (IO, 8)
              WRITE (IO, 1) 'Column permutations'
              DO 220 I = 1, PRN 
                 WRITE (IO, 5) INDEX (CPERMP+I-1)
220           CONTINUE 
              IF (PRN .LT. N) THEN 
                 WRITE (IO, 7)
              ENDIF 
              WRITE (IO, 8)
              WRITE (IO, 1) 'Row permutations'
              DO 230 I = 1, PRN 
                 WRITE (IO, 5) INDEX (RPERMP+I-1)
230           CONTINUE 
              IF (PRN .LT. N) THEN 
                 WRITE (IO, 7)
              ENDIF 
           ENDIF 

        ENDIF 

C-----------------------------------------------------------------------
C  print B (on input) or W and X (on output) for UMD2SO
C-----------------------------------------------------------------------

        IF (WHO .EQ. 3) THEN 
           WRITE (IO, 8)
           PRN = MIN (PRMAX, N)
           IF (PRL .GE. 4) THEN 
C             print all of B, or W and X
              PRN = N
           ENDIF 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 1) 'W (1 ... ',LW,
     $        ' ), work vector: not printed'
              WRITE (IO, 1) 'B (1 ... ',N,' ), right-hand side: '
              DO 240 I = 1, PRN 
                 WRITE (IO, 5) I, B (I)
240           CONTINUE 
              IF (PRN .LT. N) THEN 
                 WRITE (IO, 7)
              ENDIF 
           ELSE 
              IF (INFO (1) .LT. 0) THEN 
                 WRITE (IO, 1) 'W (1 ... ',LW,' ), work vector, and'
                 WRITE (IO, 1) 'X (1 ... ',N, ' ), solution,'
                 WRITE (IO, 1)
     $           '   not printed because of error flag, INFO (1) = ',
     $           INFO (1)
              ELSE 
                 IF (LW .EQ. 4*N) THEN 
C                   UMD2SO did iterative refinement
                    WRITE (IO, 1) 'W (1 ... ',N,' ), residual: '
                    DO 250 I = 1, PRN 
                       WRITE (IO, 5) I, W (I)
250                 CONTINUE 
                    IF (PRN .LT. N) THEN 
                       WRITE (IO, 7)
                    ENDIF 
                    WRITE (IO, 1) 'W (',N+1,' ... ',LW,
     $              ' ), work vector: not printed'
                 ELSE 
C                   no iterative refinement
                    WRITE (IO, 1) 'W (1 ... ',LW,
     $              ' ), work vector: not printed'
                 ENDIF 
                 WRITE (IO, 1) 'X (1 ... ',N,' ), solution: '
                 DO 260 I = 1, PRN 
                    WRITE (IO, 5) I, X (I)
260              CONTINUE 
                 IF (PRN .LT. N) THEN 
                    WRITE (IO, 7)
                 ENDIF 
              ENDIF 
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  who is this, and where:
C-----------------------------------------------------------------------

        IF (WHO .EQ. 1) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 6) 'end of UMD2FA input '
           ELSE 
              WRITE (IO, 6) 'end of UMD2FA output'
           ENDIF 
        ELSE IF (WHO .EQ. 2) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 6) 'end of UMD2RF input '
           ELSE 
              WRITE (IO, 6) 'end of UMD2RF output'
           ENDIF 
        ELSE IF (WHO .EQ. 3) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 6) 'end of UMD2SO input '
           ELSE 
              WRITE (IO, 6) 'end of UMD2SO output'
           ENDIF 
        ENDIF 

        RETURN

C=======================================================================
C  FORMAT STATMENTS
C=======================================================================

1       FORMAT (' ', A, :, I12, :, A, :, I12, :,
     $               A, :, I12, :, A, :, I12, :, A, :, I12)
2       FORMAT (' ', I12, ': ', I12, ' ', I12, ' ',
     $          D11.4)
3       FORMAT (' ', A, D11.4, A)
4       FORMAT (' ', A, I12, A, I12, /, A,
     $          D11.4)
5       FORMAT (' ', I12, :, ': ',
     $          D11.4)
6       FORMAT (' ', 59('='), A)
7       FORMAT ('    ...')
8       FORMAT (' ', 79 ('-'))
9       FORMAT (' ', 79 ('.'))
        END 
