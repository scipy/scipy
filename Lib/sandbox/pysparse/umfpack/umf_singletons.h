/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_singletons
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const Int Quser [ ],
    Int Cdeg [ ],
    Int Cperm [ ],
    Int Rdeg [ ],
    Int Rperm [ ],
    Int InvRperm [ ],
    Int *n1,
    Int *n1c,
    Int *n1r,
    Int *nempty_col,
    Int *nempty_row,
    Int *is_sym,
    Int *max_rdeg,
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Next [ ]
) ;
