        SUBROUTINE UMD21I (KEEP, CNTL, ICNTL)
        INTEGER ICNTL (20), KEEP (20)
        DOUBLE PRECISION
     $          CNTL (10)
        
C=== UMD21I ============================================================
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
C  HSL Compatibility:  this routine has the same arguments as MA38I/ID. 

C=======================================================================
C  USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Initialize user-controllable parameters to default values, and
C  non-user controllable parameters.  This routine is normally
C  called once prior to any call to UMD2FA.
C
C  This routine sets the default control parameters.  We recommend
C  changing these defaults under certain circumstances:
C
C  (1) If you know that your matrix has nearly symmetric nonzero
C       pattern, then we recommend setting Icntl (6) to 1 so that
C       diagonal pivoting is preferred.  This can have a significant
C       impact on the performance for matrices that are essentially
C       symmetric in pattern.
C
C   (2) If you know that your matrix is not reducible to block
C       triangular form, then we recommend setting Icntl (4) to 0
C       so that UMFPACK does not try to permute the matrix to block
C       triangular form (it will not do any useful work and will
C       leave the matrix in its irreducible form).  The work saved
C       is typically small, however.
C
C   The other control parameters typically have less effect on overall
C   performance.

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be modified on installation to reflect the computer
C  or environment in which this package is installed (printing control,
C  maximum integer, block size, and machine epsilon in particular).  If
C  you, the installer, modify this routine, please comment out the
C  original code, and add comments (with date) to describe the
C  installation.  Do not delete any original code.

C=======================================================================
C  ARGUMENTS:
C=======================================================================

C               --------------------------------------------------------
C  Icntl:       An integer array of size 20.  Need not be set by
C               caller on input.  On output, it contains default
C               integer control parameters.
C
C  Icntl (1):   Fortran output unit for error messages.
C               Default: 6
C
C  Icntl (2):   Fortran output unit for diagnostic messages.
C               Default: 6
C
C  Icntl (3):   printing-level.
C               0 or less: no output
C               1: error messages only
C               2: error messages and terse diagnostics
C               3: as 2, and print first few entries of all input and
C                       output arguments.  Invalid and duplicate entries
C                       are printed.
C               4: as 2, and print all entries of all input and
C                       output arguments.  Invalid and duplicate entries
C                       are printed.  The entire input matrix and its
C                       factors are printed.
C               5: as 4, and print out information on the data
C                       structures used to represent the LU factors,
C                       the assembly DAG, etc.
C               Default: 2
C
C  Icntl (4):   whether or not to attempt a permutation to block
C               triangular form.  If equal to one, then attempt the
C               permutation.  If you know the matrix is not reducible
C               to block triangular form, then setting Icntl (4) to
C               zero can save a small amount of computing time.
C               Default: 1 (attempt the permutation)
C
C  Icntl (5):   the number of columns to examine during the global
C               pivot search.  A value less than one is treated as one.
C               Default: 4
C
C  Icntl (6):   if not equal to zero, then pivots from the diagonal
C               of A (or the diagonal of the block-triangular form) are
C               preferred.  If the nonzero pattern of the matrix is
C               basically symmetric, we recommend that you change this
C               default value to 1 so that pivots on the diagonal
C               are preferred.
C               Default: 0 (do not prefer the diagonal)
C
C  Icntl (7):   block size for the BLAS, controlling the tradeoff
C               between the Level-2 and Level-3 BLAS.  Values less than
C               one are treated as one.
C               Default: 16, which is suitable for the CRAY YMP.
C
C  Icntl (8):   number of steps of iterative refinement to perform.
C               Values less than zero are treated as zero.  The matrix
C               must be preserved for iterative refinement to be done
C               (job=1 in UMD2FA or UMD2RF).
C               Default: 0  (no iterative refinement)
C
C  Icntl (9 ... 20):  set to zero.  Reserved for future releases.

C               --------------------------------------------------------
C  Cntl:        A double precision array of size 10.
C               Need not be set by caller on input.  On output, contains
C               default control parameters.
C
C  Cntl (1):    pivoting tradeoff between sparsity-preservation
C               and numerical stability.  An entry A(k,k) is numerically
C               acceptable if:
C                  abs (A(k,k)) >= Cntl (1) * max (abs (A(*,k)))
C               Values less than zero are treated as zero (no numerical
C               constraints).  Values greater than one are treated as
C               one (partial pivoting with row interchanges).
C               Default: 0.1
C
C  Cntl (2):    amalgamation parameter.  If the first pivot in a
C               frontal matrix has row degree r and column degree c,
C               then a working array of size
C                  (Cntl (2) * c) - by - (Cntl (2) * r)
C               is allocated for the frontal matrix.  Subsequent pivots
C               within the same frontal matrix must fit within this
C               working array, or they are not selected for this frontal
C               matrix.  Values less than one are treated as one (no
C               fill-in due to amalgamation).  Some fill-in due to
C               amalgamation is necessary for efficient use of the BLAS
C               and to reduce the assembly operations required.
C               Default: 2.0
C
C  Cntl (3):    Normally not modified by the user.
C               Defines the smallest positive number,
C               epsilon = Cntl (3), such that fl (1.0 + epsilon)
C               is greater than 1.0 (fl (x) is the floating-point
C               representation of x).  If the floating-point mantissa
C               is binary, then Cntl (3) is 2 ** (-b+1), where b
C               is the number of bits in the mantissa (including the
C               implied bit, if applicable).
C
C               Typical defaults:
C               For IEEE double precision, Cntl (3) = 2 ** (-53+1)
C               For IEEE single precision, Cntl (3) = 2 ** (-24+1)
C               For CRAY double precision, Cntl (3) = 2 ** (-96+1)
C               For CRAY single precision, Cntl (3) = 2 ** (-48+1)
C
C               A value of Cntl (3) less than or equal to zero
C               or greater than 2 ** (-15) is treated as 2 ** (-15),
C               which assumes that any floating point representation
C               has at least a 16-bit mantissa.  Cntl (3) is only
C               used in UMD2S2 to compute the sparse backward error
C               estimates, Rinfo (7) and Rinfo (8), when
C               Icntl (8) > 0 (the default is Icntl (8) = 0,
C               so by default, Cntl (3) is not used).
C
C  Cntl (4 ... 10):  set to zero.  Reserved for future releases.

C               --------------------------------------------------------
C  Keep:        An integer array of size 20.
C               Need not be set by the caller.  On output, contains
C               integer control parameters that are (normally) non-user
C               controllable (but can of course be modified by the
C               "expert" user or library installer).
C
C  Keep (1 ... 5):  unmodified (see UMD2FA or UMD2RF for a description).
C
C  Keep (6):    Largest representable positive integer.  Set to
C               2^31 - 1 = 2147483647 for 32-bit machines with 2's
C               complement arithmetic (the usual case).
C               Default: 2147483647
C
C  Keep (7) and Keep (8): A column is treated as "dense" if
C               it has more than
C               max (0, Keep(7), Keep(8)*int(sqrt(float(n))))
C               original entries.  "Dense" columns are treated
C               differently that "sparse" rows and columns.  Dense
C               columns are transformed into a priori contribution
C               blocks of dimension cdeg-by-1, where cdeg is the number
C               of original entries in the column.  Modifying these two
C               parameters can change the pivot order.
C               Default:  Keep (7) = 64
C               Default:  Keep (8) = 1
C
C  Keep (9 ... 20):  set to zero.  Reserved for future releases.

C## End of user documentation for UMD21I ###############################

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   user routine

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER I

C  i:       loop index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

C       ----------------------------------------------------------------
C       integer control parameters:
C       ----------------------------------------------------------------

        ICNTL (1) = 6
        ICNTL (2) = 6
        ICNTL (3) = 2
        ICNTL (4) = 1
        ICNTL (5) = 4
        ICNTL (6) = 0
        ICNTL (7) = 16
        ICNTL (8) = 0

C       Icntl (9 ... 20) is reserved for future releases:
        DO 10 I = 9, 20 
           ICNTL (I) = 0
10      CONTINUE 

C       ----------------------------------------------------------------
C       control parameters:
C       ----------------------------------------------------------------

        CNTL (1) = 0.1
        CNTL (2) = 2

C       IEEE double precision:  epsilon = 2 ** (-53)
        CNTL (3) = 2.0 ** (-52)

C       Cntl (4 ... 10) is reserved for future releases:
        DO 30 I = 4, 10 
           CNTL (I) = 0
30      CONTINUE 

C       ----------------------------------------------------------------
C       integer control parameters in Keep:
C       ----------------------------------------------------------------

        KEEP (6) = 2147483647
        KEEP (7) = 64
        KEEP (8) = 1

C       Keep (9 ... 20) is reserved for future releases:
        DO 20 I = 9, 20 
           KEEP (I) = 0
20      CONTINUE 

        RETURN
        END 
