      SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      INTEGER            LVL, MSUB, N, ND
*     ..
*     .. Array Arguments ..
      INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
*     ..
*
*  Purpose
*  =======
*
*  SLASDT creates a tree of subproblems for bidiagonal divide and
*  conquer.
*
*  Arguments
*  =========
*
*   N      (input) INTEGER
*          On entry, the number of diagonal elements of the
*          bidiagonal matrix.
*
*   LVL    (output) INTEGER
*          On exit, the number of levels on the computation tree.
*
*   ND     (output) INTEGER
*          On exit, the number of nodes on the tree.
*
*   INODE  (output) INTEGER array, dimension ( N )
*          On exit, centers of subproblems.
*
*   NDIML  (output) INTEGER array, dimension ( N )
*          On exit, row dimensions of left children.
*
*   NDIMR  (output) INTEGER array, dimension ( N )
*          On exit, row dimensions of right children.
*
*   MSUB   (input) INTEGER.
*          On entry, the maximum row dimension each subproblem at the
*          bottom of the tree can be of.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ming Gu and Huan Ren, Computer Science Division, University of
*     California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IR, LLST, MAXN, NCRNT, NLVL
      REAL               TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, LOG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Find the number of levels on the tree.
*
      MAXN = MAX( 1, N )
      TEMP = LOG( REAL( MAXN ) / REAL( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP ) + 1
*
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - 1
      IL = 0
      IR = 1
      LLST = 1
      DO 20 NLVL = 1, LVL - 1
*
*        Constructing the tree at (NLVL+1)-st level. The number of
*        nodes created on this level is LLST * 2.
*
         DO 10 I = 0, LLST - 1
            IL = IL + 2
            IR = IR + 2
            NCRNT = LLST + I
            NDIML( IL ) = NDIML( NCRNT ) / 2
            NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - 1
            INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - 1
            NDIML( IR ) = NDIMR( NCRNT ) / 2
            NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - 1
            INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + 1
   10    CONTINUE
         LLST = LLST*2
   20 CONTINUE
      ND = LLST*2 - 1
*
      RETURN
*
*     End of SLASDT
*
      END
