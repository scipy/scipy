        SUBROUTINE UMD2SO (N, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)
        INTEGER N, JOB, LVALUE, LINDEX, INDEX (LINDEX), KEEP (20),
     $          ICNTL (20), INFO (40)
        DOUBLE PRECISION
     $          VALUE (LVALUE), B (N), X (N), W (*)
        DOUBLE PRECISION
     $          CNTL (10), RINFO (20)
        LOGICAL TRANSC
        
C=== UMD2SO ============================================================
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
C  HSL Compatibility:  this routine has the same arguments as MA38C/CD. 

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Given LU factors computed by UMD2FA or UMD2RF, and the
C  right-hand-side, B, solve a linear system for the solution X.
C
C  This routine handles all permutations, so that B and X are in terms
C  of the original order of the matrix, A, and not in terms of the
C  permuted matrix.
C
C  If iterative refinement is done, then the residual is returned in W,
C  and the sparse backward error estimates are returned in
C  Rinfo (7) and Rinfo (8).  The computed solution X is the
C  exact solution of the equation (A + dA)x = (b + db), where
C    dA (i,j)  <= max (Rinfo (7), Rinfo (8)) * abs (A(i,j))
C  and
C    db (i) <= max (Rinfo (7) * abs (b (i)),
C                   Rinfo (8) * maxnorm (A) * maxnorm (x computed))
C  Note that dA has the same sparsity pattern as A.
C  The method used to compute the sparse backward error estimate is
C  described in M. Arioli, J. W. Demmel, and I. S. Duff, "Solving
C  sparse linear systems with sparse backward error," SIAM J. Matrix
C  Analysis and Applications, vol 10, 1989, pp. 165-190.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C           ------------------------------------------------------------
C  n:       An integer variable.
C           Must be set by caller on input (not modified).
C           Must be the same as passed to UMD2FA or UMD2RF.

C           ------------------------------------------------------------
C  job:     An integer variable.
C           Must be set by caller on input (not modified).
C           What system to solve (see the transc argument below).
C           Iterative refinement is only performed if job = 0,
C           Icntl (8) > 0, and only if the original matrix was
C           preserved (job = 1 in UMD2FA or UMD2RF).

C           ------------------------------------------------------------
C  transc:  A logical variable.
C           Must be set by caller on input (not modified).
C           solve with L and U factors or with L' and U', where
C           transa was passed to UMD2FA or UMD2RF.
C
C           If transa = false, then PAQ = LU was performed,
C           and the following systems are solved:
C
C                               transc = false          transc = true
C                               ----------------        ----------------
C                  job = 0      solve Ax = b            solve A'x = b
C                  job = 1      solve P'Lx = b          solve L'Px = b
C                  job = 2      solve UQ'x = b          solve QU'x = b
C
C           If transa = true, then A was transformed prior to LU
C           factorization, and P(A')Q = LU
C
C                               transc = false          transc = true
C                               ----------------        ----------------
C                  job = 0      solve A'x = b           solve Ax = b
C                  job = 1      solve P'Lx = b          solve L'Px = b
C                  job = 2      solve UQ'x = b          solve QU'x = b
C
C           Other values of job are treated as zero.  Iterative
C           refinement can be done only when solving Ax=b or A'x=b.
C
C           The comments below use Matlab notation, where
C           x = L \ b means x = (L^(-1)) * b, premultiplication by
C           the inverse of L.

C           ------------------------------------------------------------
C  lvalue:  An integer variable.
C           Must be set by caller on input (not modified).
C           The size of Value.

C           ------------------------------------------------------------
C  lindex:  An integer variable.
C           Must be set by caller on input (not modified).
C           The size of Index.

C           ------------------------------------------------------------
C  Value:   A double precision array of size lvalue.
C           Must be set by caller on input (normally from last call to
C           UMD2FA or UMD2RF) (not modified).
C           The LU factors, in Value (Keep (1) ... Keep (2)).
C           The entries in Value (1 ... Keep (1) - 1) and in
C           Value (Keep (2) + 1 ... lvalue) are not accessed.

C           ------------------------------------------------------------
C  Index:   An integer array of size lindex.
C           Must be set by caller on input (normally from last call to
C           UMD2FA or UMD2RF) (not modified).
C           The LU factors, in Index (Keep (3) ... Keep (5)).
C           The entries in Index (1 ... Keep (3) - 1) and in
C           Index (Keep (5) + 1 ... lindex) are not accessed.

C           ------------------------------------------------------------
C  Keep:    An integer array of size 20.
C
C           Keep (1..5): Must be set by caller on input (normally from
C               last call to UMD2FA or UMD2RF) (not modified).
C               Layout of the LU factors in Value and Index

C           ------------------------------------------------------------
C  B:       A double precision array of size n.
C           Must be set by caller on input (not modified).
C           The right hand side, b, of the system to solve.

C           ------------------------------------------------------------
C  W:       A double precision array of size 2*n or 4*n.
C           Need not be set by caller on input.  Modified on output.
C           Workspace of size W (1..2*n) if Icntl (8) = 0, which
C           is the default value.  If iterative refinement is
C           performed, and W must be of size W (1..4*n) and the
C           residual b-Ax (or b-A'x) is returned in W (1..n).

C           ------------------------------------------------------------
C  X:       A double precision array of size n.
C           Need not be set by caller on input.  Modified on output.
C           The solution, x, of the system that was solved.  Valid only
C           if Info (1) is greater than or equal to 0.

C           ------------------------------------------------------------
C  Cntl:    A double precision array of size 10.
C           Must be set by caller on input (not modified).
C           real control parameters, see UMD21I for a description,
C           which sets the defaults.  

C           ------------------------------------------------------------
C  Icntl:   An integer array of size 20.
C           Must be set by caller on input (not modified).
C           Integer control parameters, see UMD21I for a description,
C           which sets the defaults.  In particular, Icntl (8) is
C           the maximum number of steps of iterative refinement to be
C           performed.

C           ------------------------------------------------------------
C  Info:    An integer array of size 40.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of UMD2SO.
C
C           Info (1) is the error flag.  If Info (1) is -7, then
C           the LU factors are uncomputed, or have been corrupted since
C           the last call to UMD2FA or UMD2RF.  No system is solved,
C           and X (1..n) is not valid on output.  If Info (1) is 8,
C           then iterative refinement was requested but cannot be done.
C           To perform iterative refinement, the original matrix must be
C           preserved (job = 1 in UMD2FA or UMD2RF) and Ax=b or A'x=b
C           must be solved (job = 0 in UMD2SO).  Info (24) is the
C           steps of iterative refinement actually taken.

C           ------------------------------------------------------------
C  Rinfo:   A double precision array of size 20.
C           Need not be set by caller on input.  Modified on output.
C           It contains information about the execution of UMD2SO.
C
C           If iterative refinement was performed then
C           Rinfo (7) is the sparse error estimate, omega1, and
C           Rinfo (8) is the sparse error estimate, omega2.

C=======================================================================
C  TO BE PRESERVED BETWEEN CALLS TO UMD2FA, UMD2RF, UMD2SO:
C=======================================================================
C
C  The following must be unchanged since the call to UMD2FA or UMD2RF
C  that computed the LU factors:
C
C       n
C       Value (Keep (1) ... Keep (2))
C       Index (Keep (3) ... Keep (5))
C       Keep (1 ... 20)

C## End of user documentation for UMD2SO ###############################

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine
C       subroutines called:     UMD2ER, UMD2P1, UMD2S2
C       functions called:       MAX
        INTRINSIC MAX

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER NBLKS, OFFIP, OFFXP, N1, NZ, NE, OFFPP, BLKPP, LUBLPP,
     $          APP, AN, ANZ, ON, LUI1, LUI2, LUX1, LUX2, AIP, AXP,
     $          CPERMP, RPERMP, NZOFF, IRSTEP, YP, LY, LW, SP, IP1, IP2,
     $          XP1, LUIR1, IO, PRL
        LOGICAL PRESRV, BADLU

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C
C  Location and status of LU factors:
C  ----------------------------------
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for a call to UMD2RF
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors start in Value (lux1...)
C  lux2:    real part of LU factors end in Value (...lux1)
C  ip1:     pointer into leading part of LU factors in Index
C  ip2:     pointer into trailing part of LU factors in Index
C  xp1:     pointer into leading part of LU factors in Value
C  badlu:   if true, then LU factors are corrupted or not computed
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  app:     Ap (1..n+1) array located in Index (app...app+n)
C  axp:     Ax (1..nz) array located in Value (axp...axp+nz-1)
C  aip:     Ai (1..nz) array located in Index (aip...aip+nz-1)
C  an:      n if A is preserved, 1 otherwise
C  anz:     nz if A is preserved, 1 otherwise
C  offip:   Offi (1..nzoff) array loc. in Index (offip...offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array loc. in Value (offxp...offxp+nzoff-1)
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  on:      size of Offp (1..n+1):  n if nblks > 1, 1 otherwise
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  ...      seven scalars in Index (lui2-6...lui2):
C  nzoff:   number of entries in off-diagonal part
C  nblks:   number of diagonal blocks
C  presrv:  true if original matrix was preserved when factorized
C  nz:      entries in A
C  n1:      N argument in UMD2FA or UMD2RF when matrix factorized
C  ne:      NE argument in UMD2FA or UMD2RF when matrix factorized
C
C  Arrays allocated from W work array:
C  -----------------------------------
C  lw:      size of W
C  yp:      Y (1..n) located in W (yp...yp+n-1) 
C  sp:      S (1..n) located in W (sp...sp+n-1) 
C  ly:      size of Y and S
C
C  Other:
C  ------
C  irstep:  maximum number of iterative refinement steps to take

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IO = ICNTL (2)
        PRL = ICNTL (3)

C-----------------------------------------------------------------------
C  clear informational output
C-----------------------------------------------------------------------

        INFO (1) = 0
        INFO (24) = 0
        RINFO (7) = 0
        RINFO (8) = 0

C-----------------------------------------------------------------------
C  print input arguments if requested
C-----------------------------------------------------------------------

        IRSTEP = MAX (0, ICNTL (8))
        IF (IRSTEP .EQ. 0) THEN 
           LW = 2*N
        ELSE 
           LW = 4*N
        ENDIF 
        CALL UMD2P1 (3, 1,
     $          N, NE, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, N, W, LW)

C-----------------------------------------------------------------------
C  get pointers to LU factors
C-----------------------------------------------------------------------

        LUX1 = KEEP (1)
        LUX2 = KEEP (2)
        LUI1 = KEEP (3)
        LUIR1 = KEEP (4)
        LUI2 = KEEP (5)
        BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1
     $     .OR. LUI2 .GT. LINDEX
     $     .OR. LUX1 .LE. 0 .OR. LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE
     $     .OR. LUI1 .LE. 0 .OR. LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
        IF (BADLU) THEN 
           CALL UMD2ER (3, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  get seven scalars (transa, nzoff, nblks, presrv, nz, n, ne) from LU
C-----------------------------------------------------------------------

        NE = INDEX (LUI2)
        N1 = INDEX (LUI2-1)
        NZ = INDEX (LUI2-2)
        PRESRV = INDEX (LUI2-3) .NE. 0
        NBLKS = INDEX (LUI2-4)
        NZOFF = INDEX (LUI2-5)
C       transa = Index (lui2-6) .ne. 0, we don't actually need this here

C-----------------------------------------------------------------------
C  get pointers to permutation vectors
C-----------------------------------------------------------------------

        RPERMP = (LUI2-6) - N
        CPERMP = RPERMP - N
        IP2 = CPERMP - 1
        XP1 = LUX1
        IP1 = LUI1

C-----------------------------------------------------------------------
C  get pointers to preserved column-oriented copy of input matrix
C-----------------------------------------------------------------------

        IF (PRESRV) THEN 

C          -------------------------------------------------------------
C          original matrix preserved in Index (lui1..lui1+nz+n) and
C          Value (lux1..lux1+nz-1)
C          -------------------------------------------------------------

           APP = IP1
           AIP = APP + N+1
           IP1 = AIP + NZ
           AXP = XP1
           XP1 = AXP + NZ
           AN = N
           ANZ = NZ

        ELSE 

C          -------------------------------------------------------------
C          original matrix not preserved, pass dummy argument to UMD2S2
C          -------------------------------------------------------------

           APP = 1
           AIP = 1
           AXP = 1
           AN = 1
           ANZ = 1

        ENDIF 

C-----------------------------------------------------------------------
C  get pointers to block-triangular information, if BTF was used
C-----------------------------------------------------------------------

        IF (NBLKS .GT. 1) THEN 

C          -------------------------------------------------------------
C          get pointers to off-diagonal nonzeros, and BTF arrays
C          -------------------------------------------------------------

           OFFIP = IP1
           IP1 = IP1 + NZOFF
           OFFXP = XP1
           XP1 = XP1 + NZOFF
           OFFPP = CPERMP - (N+1)
           BLKPP = OFFPP - (NBLKS+1)
           LUBLPP = BLKPP - (NBLKS)
           IP2 = LUBLPP - 1
           ON = N

        ELSE 

C          -------------------------------------------------------------
C          matrix was factorized as a single block, pass dummy arg.
C          -------------------------------------------------------------

           OFFIP = 1
           OFFXP = 1
           OFFPP = 1
           BLKPP = 1
           LUBLPP = 1
           ON = 1

        ENDIF 

        BADLU = N .NE. N1 .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $     NBLKS .LE. 0 .OR. NBLKS .GT. N .OR.
     $     XP1 .GT. LUX2 .OR. NZOFF .LT. 0 .OR. IP1 .NE. LUIR1
        IF (BADLU) THEN 
           CALL UMD2ER (3, ICNTL, INFO, -7, 0)
C          error return, LU factors are corrupted:
           GO TO 9000
        ENDIF 

C-----------------------------------------------------------------------
C  get the number of steps of iterative refinement
C-----------------------------------------------------------------------

        IF (IRSTEP .GT. 0 .AND. .NOT. PRESRV) THEN 
C          original matrix not preserved (UMD2FA/UMD2RF job .ne. 1)
           CALL UMD2ER (3, ICNTL, INFO, 8, 0)
           IRSTEP = 0
        ENDIF 
        IF (IRSTEP .GT. 0 .AND. (JOB .EQ. 1 .OR. JOB .EQ. 2)) THEN 
C          iterative refinement for Ax=b and A'x=b only (job = 0)
           CALL UMD2ER (3, ICNTL, INFO, 8, 1)
           IRSTEP = 0
        ENDIF 
        IF (IRSTEP .EQ. 0) THEN 
C          pass a dummy argument as Y, which is not accessed in UMD2S2
           YP = 1
           LY = 1
           SP = 1
           LW = 2*N
        ELSE 
C          pass W (yp ... yp+n-1) as Y (1..n) to UMD2S2
           YP = 2*N+1
           LY = N
           SP = 3*N+1
           LW = 4*N
        ENDIF 

C-----------------------------------------------------------------------
C  solve; optional iterative refinement and sparse backward error
C-----------------------------------------------------------------------

        CALL UMD2S2 (N, JOB, TRANSC, LUX2-XP1+1, VALUE (XP1),
     $     IP2-LUIR1+1, INDEX (LUIR1), B, X,
     $     W, W (N+1), LY, W (YP), W (SP),
     $     CNTL, ICNTL, INFO, RINFO, INDEX (CPERMP), INDEX (RPERMP),
     $     PRESRV, AN, ANZ, INDEX (APP), INDEX (AIP), VALUE (AXP),
     $     ON, MAX (1, NZOFF), INDEX (OFFPP), INDEX (OFFIP),
     $     VALUE (OFFXP), NBLKS, INDEX (LUBLPP), INDEX (BLKPP), IRSTEP)

C-----------------------------------------------------------------------
C  print output arguments if requested
C-----------------------------------------------------------------------

C       error return label:
9000    CONTINUE

        CALL UMD2P1 (3, 2,
     $          N, NE, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, N, W, LW)
        RETURN
        END 
