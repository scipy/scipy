/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* umf_dump.h: debugging definitions. */

#ifndef NDEBUG

GLOBAL void UMF_dump_dense
(
    Entry *C,
    Int dim,
    Int m,
    Int n
) ;

GLOBAL void UMF_dump_element
(
    NumericType *Numeric,
    WorkType *Work,
    Int e,
    Int clean
) ;

GLOBAL void UMF_dump_rowcol
(
    Int dump_which,
    NumericType *Numeric,
    WorkType *Work,
    Int dump_index,
    Int check_degree
) ;

GLOBAL void UMF_dump_matrix
(
    NumericType *Numeric,
    WorkType *Work,
    Int check_degree
) ;

GLOBAL void UMF_dump_current_front
(
    NumericType *Numeric,
    WorkType *Work,
    Int check
) ;

GLOBAL void UMF_dump_lu
(
    NumericType *Numeric
) ;

GLOBAL void UMF_dump_memory
(
    NumericType *Numeric
) ;

GLOBAL void UMF_dump_packed_memory
(
    NumericType *Numeric,
    WorkType *Work
) ;

GLOBAL void UMF_dump_col_matrix
(
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    const Int Ai [ ],
    const Int Ap [ ],
    Int n_row,
    Int n_col,
    Int nz
) ;

GLOBAL void UMF_dump_chain
(
    Int frontid,
    Int Front_parent [ ],
    Int Front_npivcol [ ],
    Int Front_nrows [ ],
    Int Front_ncols [ ],
    Int nfr
) ;

GLOBAL void UMF_dump_rowmerge
(
    NumericType *Numeric,
    SymbolicType *Symbolic,
    WorkType *Work
) ;

GLOBAL void UMF_dump_start
(
    void
) ;


GLOBAL void UMF_dump_diagonal_map
(
    Int Diagonal_map [ ],
    Int Diagonal_imap [ ],
    Int n1,
    Int nn,
    Int nempty
) ;

#define UMF_DBMAX 50000
GLOBAL extern Int UMF_debug ;
GLOBAL extern Int UMF_allocfail ;
GLOBAL extern double UMF_gprob ;

#define DEBUGk(k,params) { if (UMF_debug >= (k)) { PRINTF (params) ; } }

#define DEBUGm4(params) DEBUGk (-4, params)
#define DEBUGm3(params) DEBUGk (-3, params)
#define DEBUGm2(params) DEBUGk (-2, params)
#define DEBUGm1(params) DEBUGk (-1, params)
#define DEBUG0(params) DEBUGk (0, params)
#define DEBUG1(params) DEBUGk (1, params)
#define DEBUG2(params) DEBUGk (2, params)
#define DEBUG3(params) DEBUGk (3, params)
#define DEBUG4(params) DEBUGk (4, params)
#define DEBUG5(params) DEBUGk (5, params)
#define DEBUG6(params) DEBUGk (6, params)
#define DEBUG7(params) DEBUGk (7, params)
#define DEBUG8(params) DEBUGk (8, params)
#define DEBUG9(params) DEBUGk (9, params)

#define EDEBUGk(k,a) { if (UMF_debug >= (k)) { PRINT_ENTRY (a) ; } }

#define EDEBUG0(a) EDEBUGk (0, a)
#define EDEBUG1(a) EDEBUGk (1, a)
#define EDEBUG2(a) EDEBUGk (2, a)
#define EDEBUG3(a) EDEBUGk (3, a)
#define EDEBUG4(a) EDEBUGk (4, a)
#define EDEBUG5(a) EDEBUGk (5, a)
#define EDEBUG6(a) EDEBUGk (6, a)
#define EDEBUG7(a) EDEBUGk (7, a)
#define EDEBUG8(a) EDEBUGk (8, a)
#define EDEBUG9(a) EDEBUGk (9, a)

/* ASSERT defined in amd_dump.h */

#else

/* ========================================================================== */
/* === No debugging ========================================================= */
/* ========================================================================== */

/* turn off all debugging macros */

#define DEBUGk(k,params)

#define DEBUGm4(params)
#define DEBUGm3(params)
#define DEBUGm2(params)
#define DEBUGm1(params)
#define DEBUG0(params)
#define DEBUG1(params)
#define DEBUG2(params)
#define DEBUG3(params)
#define DEBUG4(params)
#define DEBUG5(params)
#define DEBUG6(params)
#define DEBUG7(params)
#define DEBUG8(params)
#define DEBUG9(params)

#define EDEBUGk(k,a)

#define EDEBUG0(a)
#define EDEBUG1(a)
#define EDEBUG2(a)
#define EDEBUG3(a)
#define EDEBUG4(a)
#define EDEBUG5(a)
#define EDEBUG6(a)
#define EDEBUG7(a)
#define EDEBUG8(a)
#define EDEBUG9(a)

#endif /* NDEBUG */
