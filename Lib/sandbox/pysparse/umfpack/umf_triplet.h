/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_triplet_map_x
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],
    const Int Tj [ ],
    Int Ap [ ],
    Int Ai [ ],
    Int Rp [ ],
    Int Rj [ ],
    Int W [ ],
    Int RowCount [ ]
    , const double Tx [ ]
    , double Ax [ ]
    , double Rx [ ]
#ifdef COMPLEX
    , const double Tz [ ]
    , double Az [ ]
    , double Rz [ ]
#endif
    , Int Map [ ]
    , Int Map2 [ ]
) ;

GLOBAL Int UMF_triplet_map_nox
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],
    const Int Tj [ ],
    Int Ap [ ],
    Int Ai [ ],
    Int Rp [ ],
    Int Rj [ ],
    Int W [ ],
    Int RowCount [ ]
    , Int Map [ ]
    , Int Map2 [ ]
) ;

GLOBAL Int UMF_triplet_nomap_x
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],
    const Int Tj [ ],
    Int Ap [ ],
    Int Ai [ ],
    Int Rp [ ],
    Int Rj [ ],
    Int W [ ],
    Int RowCount [ ]
    , const double Tx [ ]
    , double Ax [ ]
    , double Rx [ ]
#ifdef COMPLEX
    , const double Tz [ ]
    , double Az [ ]
    , double Rz [ ]
#endif
) ;

GLOBAL Int UMF_triplet_nomap_nox
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],
    const Int Tj [ ],
    Int Ap [ ],
    Int Ai [ ],
    Int Rp [ ],
    Int Rj [ ],
    Int W [ ],
    Int RowCount [ ]
) ;
