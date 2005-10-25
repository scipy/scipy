/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_row_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic,
    Int cdeg0,
    Int cdeg1,
    const Int Pattern [ ],
    const Int Pos [ ],
    Int pivrow [2],
    Int rdeg [2],
    Int W_i [ ],
    Int W_o [ ],
    Int prior_pivrow [2],
    const Entry Wxy [ ],
    Int pivcol,
    Int freebie [2]
) ;

#define IN 0
#define OUT 1

#define IN_IN 0
#define IN_OUT 1
#define OUT_IN 2
#define OUT_OUT 3
