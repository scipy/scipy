/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_transpose
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    const Int P [ ],
    const Int Q [ ],
    Int nq,
    Int Rp [ ],
    Int Ri [ ],
    double Rx [ ],
    Int W [ ],
    Int check
#ifdef COMPLEX
    , const double Az [ ]
    , double Rz [ ]
    , Int do_conjugate
#endif
) ;
