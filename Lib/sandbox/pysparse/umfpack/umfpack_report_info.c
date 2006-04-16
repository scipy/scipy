/* ========================================================================== */
/* === UMFPACK_report_info ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints the Info array.  See umfpack_report_info.h for
    details.
*/

#include "umf_internal.h"

#define PRINT_INFO(format,x) \
{ \
    if (SCALAR_IS_NAN (x) || (!SCALAR_IS_LTZERO (x))) \
    { \
	PRINTF ((format, x)) ; \
    } \
}

/* RATIO macro uses a double relop, but ignore NaN case: */
#define RATIO(a,b,c) (((b) == 0) ? (c) : (((double) a)/((double) b)))

/* ========================================================================== */
/* === print_ratio ========================================================== */
/* ========================================================================== */

PRIVATE void print_ratio
(
    char *what,
    char *format,
    double estimate,
    double actual
)
{
    if (estimate < 0 && actual < 0)	/* double relop, but ignore Nan case */
    {
	return ;
    }
    PRINTF (("    %-27s", what)) ;
    if (estimate >= 0)			/* double relop, but ignore Nan case */
    {
	PRINTF ((format, estimate)) ;
    }
    else
    {
	PRINTF (("                    -")) ;
    }
    if (actual >= 0)			/* double relop, but ignore Nan case */
    {
	PRINTF ((format, actual)) ;
    }
    else
    {
	PRINTF (("                    -")) ;
    }
    if (estimate >= 0 && actual >= 0)	/* double relop, but ignore Nan case */
    {
	PRINTF ((" %5.0f%%\n", 100 * RATIO (actual, estimate, 1))) ;
    }
    else
    {
	PRINTF (("      -\n")) ;
    }
}

/* ========================================================================== */
/* === UMFPACK_report_info ================================================== */
/* ========================================================================== */

GLOBAL void UMFPACK_report_info
(
    const double Control [UMFPACK_CONTROL],
    const double Info [UMFPACK_INFO]
)
{

    double lnz_est, unz_est, lunz_est, lnz, unz, lunz, tsym, tnum, fnum, tsolve,
	fsolve, ttot, ftot, twsym, twnum, twsolve, twtot, n2 ;
    Int n_row, n_col, n_inner, prl, is_sym ;

    /* ---------------------------------------------------------------------- */
    /* get control settings and status to determine what to print */
    /* ---------------------------------------------------------------------- */

    prl = GET_CONTROL (UMFPACK_PRL, UMFPACK_DEFAULT_PRL) ;

    if (!Info || prl < 2)
    {
	/* no output generated if Info is (double *) NULL */
	/* or if prl is less than 2 */
	return ;
    }

    /* ---------------------------------------------------------------------- */
    /* print umfpack version */
    /* ---------------------------------------------------------------------- */

    PRINTF (("\n%s, Info:\n", UMFPACK_VERSION)) ;

#ifndef NDEBUG
    PRINTF ((
"**** Debugging enabled (UMFPACK will be exceedingly slow!) *****************\n"
    )) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* print run-time options */
    /* ---------------------------------------------------------------------- */

#ifdef DINT
    PRINTF (("    matrix entry defined as:          double\n")) ;
    PRINTF (("    Int (generic integer) defined as: int\n")) ;
#endif
#ifdef DLONG
    PRINTF (("    matrix entry defined as:          double\n")) ;
    PRINTF (("    Int (generic integer) defined as: long\n")) ;
#endif
#ifdef ZINT
    PRINTF (("    matrix entry defined as:          double complex\n")) ;
    PRINTF (("    Int (generic integer) defined as: int\n")) ;
#endif
#ifdef ZLONG
    PRINTF (("    matrix entry defined as:          double complex\n")) ;
    PRINTF (("    Int (generic integer) defined as: long\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* print compile-time options */
    /* ---------------------------------------------------------------------- */

    PRINTF (("    BLAS library used:                ")) ;

#if defined (USE_NO_BLAS)
    PRINTF (("none.  UMFPACK will be slow.\n")) ;
#elif defined (USE_C_BLAS)
    PRINTF (("C-BLAS.\n")) ;
#elif defined (USE_MATLAB_BLAS)
    PRINTF (("built-in MATLAB BLAS.\n")) ;
#elif defined (USE_SUNPERF_BLAS)
    PRINTF (("Sun Performance Library BLAS.\n")) ;
#elif defined (USE_SCSL_BLAS)
    PRINTF (("SGI SCSL BLAS.\n")) ;
#elif defined (USE_FORTRAN_BLAS)
    PRINTF (("Fortran BLAS.\n")) ;
#endif

    PRINTF (("    MATLAB:                           ")) ;
#ifdef MATLAB_MEX_FILE
    PRINTF (("yes.\n")) ;
#else
#ifdef MATHWORKS
    PRINTF (("yes (using internal ut* routines).\n")) ;
#else
    PRINTF (("no.\n")) ;
#endif
#endif

    PRINTF (("    CPU timer:                        ")) ;
#ifdef GETRUSAGE
    PRINTF (("getrusage.\n")) ;
#else
    PRINTF (("ANSI C clock (may wrap-around).\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* print n and nz */
    /* ---------------------------------------------------------------------- */

    n_row = (Int) Info [UMFPACK_NROW] ;
    n_col = (Int) Info [UMFPACK_NCOL] ;
    n_inner = MIN (n_row, n_col) ;

    PRINT_INFO ("    number of rows in matrix A:       "ID"\n", n_row) ;
    PRINT_INFO ("    number of columns in matrix A:    "ID"\n", n_col) ;
    PRINT_INFO ("    entries in matrix A:              "ID"\n",
	(Int) Info [UMFPACK_NZ]) ;
    PRINT_INFO ("    memory usage reported in:         "ID"-byte Units\n",
	(Int) Info [UMFPACK_SIZE_OF_UNIT]) ;

    PRINT_INFO ("    size of int:                      "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_INT]) ;
    PRINT_INFO ("    size of long:                     "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_LONG]) ;
    PRINT_INFO ("    size of pointer:                  "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_POINTER]) ;
    PRINT_INFO ("    size of numerical entry:          "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_ENTRY]) ;

    /* ---------------------------------------------------------------------- */
    /* symbolic parameters */
    /* ---------------------------------------------------------------------- */

    if (Info [UMFPACK_STRATEGY_USED] == UMFPACK_STRATEGY_SYMMETRIC)
    {
	PRINTF (("\n    strategy used:                    symmetric\n")) ;
    }
    else if (Info [UMFPACK_STRATEGY_USED] == UMFPACK_STRATEGY_UNSYMMETRIC)
    {
	PRINTF (("\n    strategy used:                    unsymmetric\n")) ;
    }
    else if (Info [UMFPACK_STRATEGY_USED] == UMFPACK_STRATEGY_2BY2)
    {
	PRINTF (("\n    strategy used:                    symmetric 2-by-2\n"));
    }

    if (Info [UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_AMD)
    {
	PRINTF (("    ordering used:                    amd on A+A'\n")) ;
    }
    else if (Info [UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_COLAMD)
    {
	PRINTF (("    ordering used:                    colamd on A\n")) ;
    }
    else if (Info [UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_GIVEN)
    {
	PRINTF (("    ordering used:                    provided by user\n")) ;
    }

    if (Info [UMFPACK_QFIXED] == 1)
    {
	PRINTF (("    modify Q during factorization:    no\n")) ;
    }
    else if (Info [UMFPACK_QFIXED] == 0)
    {
	PRINTF (("    modify Q during factorization:    yes\n")) ;
    }

    if (Info [UMFPACK_DIAG_PREFERRED] == 0)
    {
	PRINTF (("    prefer diagonal pivoting:         no\n")) ;
    }
    else if (Info [UMFPACK_DIAG_PREFERRED] == 1)
    {
	PRINTF (("    prefer diagonal pivoting:         yes\n")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* singleton statistics */
    /* ---------------------------------------------------------------------- */

    PRINT_INFO ("    pivots with zero Markowitz cost:               %0.f\n",
	Info [UMFPACK_COL_SINGLETONS] + Info [UMFPACK_ROW_SINGLETONS]) ;
    PRINT_INFO ("    submatrix S after removing zero-cost pivots:\n"
		"        number of \"dense\" rows:                    %.0f\n",
	Info [UMFPACK_NDENSE_ROW]) ;
    PRINT_INFO ("        number of \"dense\" columns:                 %.0f\n",
	Info [UMFPACK_NDENSE_COL]) ;
    PRINT_INFO ("        number of empty rows:                      %.0f\n",
	Info [UMFPACK_NEMPTY_ROW]) ;
    PRINT_INFO ("        number of empty columns                    %.0f\n",
	Info [UMFPACK_NEMPTY_COL]) ;
    is_sym = Info [UMFPACK_S_SYMMETRIC] ;
    if (is_sym > 0)
    {
	PRINTF (("        submatrix S square and diagonal preserved\n")) ;
    }
    else if (is_sym == 0)
    {
	PRINTF (("        submatrix S not square or diagonal not preserved\n"));
    }

    /* ---------------------------------------------------------------------- */
    /* statistics from amd_aat */
    /* ---------------------------------------------------------------------- */

    n2 = Info [UMFPACK_N2] ;
    if (n2 >= 0)
    {
	PRINTF (("    pattern of square submatrix S:\n")) ;
    }
    PRINT_INFO ("        number rows and columns                    %.0f\n",
	n2) ;
    PRINT_INFO ("        symmetry of nonzero pattern:               %.6f\n",
	Info [UMFPACK_PATTERN_SYMMETRY]) ;
    PRINT_INFO ("        nz in S+S' (excl. diagonal):               %.0f\n",
	Info [UMFPACK_NZ_A_PLUS_AT]) ;
    PRINT_INFO ("        nz on diagonal of matrix S:                %.0f\n",
	Info [UMFPACK_NZDIAG]) ;
    if (Info [UMFPACK_NZDIAG] >= 0 && n2 > 0)
    {
	PRINTF (("        fraction of nz on diagonal:                %.6f\n",
	Info [UMFPACK_NZDIAG] / n2)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* statistics from 2-by-2 permutation */
    /* ---------------------------------------------------------------------- */

    PRINT_INFO ("    2-by-2 pivoting to place large entries on diagonal:\n"
		"        # of small diagonal entries of S:          %.0f\n",
	Info [UMFPACK_2BY2_NWEAK]) ;
    PRINT_INFO ("        # unmatched:                               %.0f\n",
	Info [UMFPACK_2BY2_UNMATCHED]) ;
    PRINT_INFO ("        symmetry of P2*S:                          %.6f\n",
	Info [UMFPACK_2BY2_PATTERN_SYMMETRY]) ;
    PRINT_INFO ("        nz in P2*S+(P2*S)' (excl. diag.):          %.0f\n",
	Info [UMFPACK_2BY2_NZ_PA_PLUS_PAT]) ;
    PRINT_INFO ("        nz on diagonal of P2*S:                    %.0f\n",
	Info [UMFPACK_2BY2_NZDIAG]) ;
    if (Info [UMFPACK_2BY2_NZDIAG] >= 0 && n2 > 0)
    {
	PRINTF (("        fraction of nz on diag of P2*S:            %.6f\n",
	Info [UMFPACK_2BY2_NZDIAG] / n2)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* statistics from AMD */
    /* ---------------------------------------------------------------------- */

    if (Info [UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_AMD)
    {
	double dmax = Info [UMFPACK_SYMMETRIC_DMAX] ;
	PRINTF (("    AMD statistics, for strict diagonal pivoting:\n")) ;
	PRINT_INFO ("        est. flops for LU factorization:           %.5e\n",
	    Info [UMFPACK_SYMMETRIC_FLOPS]) ;
	PRINT_INFO ("        est. nz in L+U (incl. diagonal):           %.0f\n",
	    Info [UMFPACK_SYMMETRIC_LUNZ]) ;
	PRINT_INFO ("        est. largest front (# entries):            %.0f\n",
	    dmax*dmax) ;
	PRINT_INFO ("        est. max nz in any column of L:            %.0f\n",
	    dmax) ;
	PRINT_INFO (
	    "        number of \"dense\" rows/columns in S+S':    %.0f\n",
	    Info [UMFPACK_SYMMETRIC_NDENSE]) ;
    }

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization */
    /* ---------------------------------------------------------------------- */

    tsym = Info [UMFPACK_SYMBOLIC_TIME] ;
    twsym = Info [UMFPACK_SYMBOLIC_WALLTIME] ;

    PRINT_INFO ("    symbolic factorization defragmentations:       %.0f\n",
	Info [UMFPACK_SYMBOLIC_DEFRAG]) ;
    PRINT_INFO ("    symbolic memory usage (Units):                 %.0f\n",
	Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]) ;
    PRINT_INFO ("    symbolic memory usage (MBytes):                %.1f\n",
	MBYTES (Info [UMFPACK_SYMBOLIC_PEAK_MEMORY])) ;
    PRINT_INFO ("    Symbolic size (Units):                         %.0f\n",
	Info [UMFPACK_SYMBOLIC_SIZE]) ;
    PRINT_INFO ("    Symbolic size (MBytes):                        %.0f\n",
	MBYTES (Info [UMFPACK_SYMBOLIC_SIZE])) ;
    PRINT_INFO ("    symbolic factorization CPU time (sec):         %.2f\n",
	tsym) ;
    PRINT_INFO ("    symbolic factorization wallclock time(sec):    %.2f\n",
	twsym) ;

    /* ---------------------------------------------------------------------- */
    /* scaling, from numerical factorization */
    /* ---------------------------------------------------------------------- */

    if (Info [UMFPACK_WAS_SCALED] == UMFPACK_SCALE_NONE)
    {
	PRINTF (("\n    matrix scaled: no\n")) ;
    }
    else if (Info [UMFPACK_WAS_SCALED] == UMFPACK_SCALE_SUM)
    {
	PRINTF (("\n    matrix scaled: yes ")) ;
	PRINTF (("(divided each row by sum of abs values in each row)\n")) ;
	PRINTF (("    minimum sum (abs (rows of A)):              %.5e\n",
	    Info [UMFPACK_RSMIN])) ;
	PRINTF (("    maximum sum (abs (rows of A)):              %.5e\n",
	    Info [UMFPACK_RSMAX])) ;
    }
    else if (Info [UMFPACK_WAS_SCALED] == UMFPACK_SCALE_MAX)
    {
	PRINTF (("\n    matrix scaled: yes ")) ;
	PRINTF (("(divided each row by max abs value in each row)\n")) ;
	PRINTF (("    minimum max (abs (rows of A)):              %.5e\n",
	    Info [UMFPACK_RSMIN])) ;
	PRINTF (("    maximum max (abs (rows of A)):              %.5e\n",
	    Info [UMFPACK_RSMAX])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* estimate/actual in symbolic/numeric factorization */
    /* ---------------------------------------------------------------------- */

    /* double relop, but ignore NaN case: */
    if (Info [UMFPACK_SYMBOLIC_DEFRAG] >= 0	/* UMFPACK_*symbolic called */
    ||  Info [UMFPACK_NUMERIC_DEFRAG] >= 0)	/* UMFPACK_numeric called */
    {
	PRINTF (("\n    symbolic/numeric factorization:      upper bound")) ;
	PRINTF (("               actual      %%\n")) ;
	PRINTF (("    variable-sized part of Numeric object:\n")) ;
    }
    print_ratio ("    initial size (Units)", " %20.0f",
	Info [UMFPACK_VARIABLE_INIT_ESTIMATE], Info [UMFPACK_VARIABLE_INIT]) ;
    print_ratio ("    peak size (Units)", " %20.0f",
	Info [UMFPACK_VARIABLE_PEAK_ESTIMATE], Info [UMFPACK_VARIABLE_PEAK]) ;
    print_ratio ("    final size (Units)", " %20.0f",
	Info [UMFPACK_VARIABLE_FINAL_ESTIMATE], Info [UMFPACK_VARIABLE_FINAL]) ;
    print_ratio ("Numeric final size (Units)", " %20.0f",
	Info [UMFPACK_NUMERIC_SIZE_ESTIMATE], Info [UMFPACK_NUMERIC_SIZE]) ;
    print_ratio ("Numeric final size (MBytes)", " %20.1f",
	MBYTES (Info [UMFPACK_NUMERIC_SIZE_ESTIMATE]),
	MBYTES (Info [UMFPACK_NUMERIC_SIZE])) ;
    print_ratio ("peak memory usage (Units)", " %20.0f",
	Info [UMFPACK_PEAK_MEMORY_ESTIMATE], Info [UMFPACK_PEAK_MEMORY]) ;
    print_ratio ("peak memory usage (MBytes)", " %20.1f",
	MBYTES (Info [UMFPACK_PEAK_MEMORY_ESTIMATE]),
	MBYTES (Info [UMFPACK_PEAK_MEMORY])) ;
    print_ratio ("numeric factorization flops", " %20.5e",
	Info [UMFPACK_FLOPS_ESTIMATE], Info [UMFPACK_FLOPS]) ;

    lnz_est = Info [UMFPACK_LNZ_ESTIMATE] ;
    unz_est = Info [UMFPACK_UNZ_ESTIMATE] ;
    if (lnz_est >= 0 && unz_est >= 0)	/* double relop, but ignore NaN case */
    {
	lunz_est = lnz_est + unz_est - n_inner ;
    }
    else
    {
	lunz_est = EMPTY ;
    }
    lnz = Info [UMFPACK_LNZ] ;
    unz = Info [UMFPACK_UNZ] ;
    if (lnz >= 0 && unz >= 0)		/* double relop, but ignore NaN case */
    {
	lunz = lnz + unz - n_inner ;
    }
    else
    {
	lunz = EMPTY ;
    }
    print_ratio ("nz in L (incl diagonal)", " %20.0f", lnz_est, lnz) ;
    print_ratio ("nz in U (incl diagonal)", " %20.0f", unz_est, unz) ;
    print_ratio ("nz in L+U (incl diagonal)", " %20.0f", lunz_est, lunz) ;

    print_ratio ("largest front (# entries)", " %20.0f",
	Info [UMFPACK_MAX_FRONT_SIZE_ESTIMATE], Info [UMFPACK_MAX_FRONT_SIZE]) ;
    print_ratio ("largest # rows in front", " %20.0f",
	Info [UMFPACK_MAX_FRONT_NROWS_ESTIMATE],
	Info [UMFPACK_MAX_FRONT_NROWS]) ;
    print_ratio ("largest # columns in front", " %20.0f",
	Info [UMFPACK_MAX_FRONT_NCOLS_ESTIMATE],
	Info [UMFPACK_MAX_FRONT_NCOLS]) ;

    /* ---------------------------------------------------------------------- */
    /* numeric factorization */
    /* ---------------------------------------------------------------------- */

    tnum = Info [UMFPACK_NUMERIC_TIME] ;
    twnum = Info [UMFPACK_NUMERIC_WALLTIME] ;
    fnum = Info [UMFPACK_FLOPS] ;

    PRINT_INFO ("\n    initial allocation ratio used:                 %0.3g\n",
	Info [UMFPACK_ALLOC_INIT_USED]) ;
    PRINT_INFO ("    # of forced updates due to frontal growth:     %.0f\n",
	Info [UMFPACK_FORCED_UPDATES]) ;
    PRINT_INFO ("    number of off-diagonal pivots:                 %.0f\n",
	Info [UMFPACK_NOFF_DIAG]) ;
    PRINT_INFO ("    nonzeros on diagonal of U:                     %.0f\n",
	Info [UMFPACK_UDIAG_NZ]) ;
    PRINT_INFO ("    min abs. value on diagonal of U:               %.2e\n",
	Info [UMFPACK_UMIN]) ;
    PRINT_INFO ("    max abs. value on diagonal of U:               %.2e\n",
	Info [UMFPACK_UMAX]) ;
    PRINT_INFO ("    estimate of reciprocal of condition number:    %.2e\n",
	Info [UMFPACK_RCOND]) ;
    PRINT_INFO ("    indices in compressed pattern:                 %.0f\n",
	Info [UMFPACK_COMPRESSED_PATTERN]) ;
    PRINT_INFO ("    numerical values stored in Numeric object:     %.0f\n",
	Info [UMFPACK_LU_ENTRIES]) ;
    PRINT_INFO ("    numeric factorization defragmentations:        %.0f\n",
	Info [UMFPACK_NUMERIC_DEFRAG]) ;
    PRINT_INFO ("    numeric factorization reallocations:           %.0f\n",
	Info [UMFPACK_NUMERIC_REALLOC]) ;
    PRINT_INFO ("    costly numeric factorization reallocations:    %.0f\n",
	Info [UMFPACK_NUMERIC_COSTLY_REALLOC]) ;
    PRINT_INFO ("    numeric factorization CPU time (sec):          %.2f\n",
	tnum) ;
    PRINT_INFO ("    numeric factorization wallclock time (sec):    %.2f\n",
	twnum) ;

    if (tnum > 0 && fnum > 0)
    {
	PRINT_INFO (
	   "    numeric factorization mflops (CPU time):       %.2f\n",
	   1e-6 * fnum / tnum) ;
    }
    if (twnum > 0 && fnum > 0)
    {
	PRINT_INFO (
	   "    numeric factorization mflops (wallclock):      %.2f\n",
	   1e-6 * fnum / twnum) ;
    }

    ttot = EMPTY ;
    ftot = fnum ;
    if (tsym >= 0 && tnum >= 0)
    {
	ttot = tsym + tnum ;
	PRINT_INFO ("    symbolic + numeric CPU time (sec):             %.2f\n",
	    ttot) ;
	if (ftot > 0 && ttot > 0)
	{
	    PRINT_INFO (
		"    symbolic + numeric mflops (CPU time):          %.2f\n",
		1e-6 * ftot / ttot) ;
	}
    }

    twtot = EMPTY ;
    if (twsym >= 0 && twnum >= 0)
    {
	twtot = twsym + twnum ;
	PRINT_INFO ("    symbolic + numeric wall clock time (sec):      %.2f\n",
	    twtot) ;
	if (ftot > 0 && twtot > 0)
	{
	    PRINT_INFO (
		"    symbolic + numeric mflops (wall clock):        %.2f\n",
		1e-6 * ftot / twtot) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* solve */
    /* ---------------------------------------------------------------------- */

    tsolve = Info [UMFPACK_SOLVE_TIME] ;
    twsolve = Info [UMFPACK_SOLVE_WALLTIME] ;
    fsolve = Info [UMFPACK_SOLVE_FLOPS] ;

    PRINT_INFO ("\n    solve flops:                                   %.5e\n",
	fsolve) ;
    PRINT_INFO ("    iterative refinement steps taken:              %.0f\n",
	Info [UMFPACK_IR_TAKEN]) ;
    PRINT_INFO ("    iterative refinement steps attempted:          %.0f\n",
	Info [UMFPACK_IR_ATTEMPTED]) ;
    PRINT_INFO ("    sparse backward error omega1:                  %.2e\n",
	Info [UMFPACK_OMEGA1]) ;
    PRINT_INFO ("    sparse backward error omega2:                  %.2e\n",
	Info [UMFPACK_OMEGA2]) ;
    PRINT_INFO ("    solve CPU time (sec):                          %.2f\n",
	tsolve) ;
    PRINT_INFO ("    solve wall clock time (sec):                   %.2f\n",
	twsolve) ;
    if (fsolve > 0 && tsolve > 0)
    {
	PRINT_INFO (
	    "    solve mflops (CPU time):                       %.2f\n",
	    1e-6 * fsolve / tsolve) ;
    }
    if (fsolve > 0 && twsolve > 0)
    {
	PRINT_INFO (
	    "    solve mflops (wall clock time):                %.2f\n",
	    1e-6 * fsolve / twsolve) ;
    }

    if (ftot >= 0 && fsolve >= 0)
    {
	ftot += fsolve ;
	PRINT_INFO (
	"\n    total symbolic + numeric + solve flops:        %.5e\n", ftot) ;
    }

    if (tsolve >= 0)
    {
	if (ttot >= 0 && ftot >= 0)
	{
	    ttot += tsolve ;
	    PRINT_INFO (
		"    total symbolic + numeric + solve CPU time:     %.2f\n",
		ttot) ;
	    if (ftot > 0 && ttot > 0)
	    {
		PRINT_INFO (
		"    total symbolic + numeric + solve mflops (CPU): %.2f\n",
		1e-6 * ftot / ttot) ;
	    }
	}
    }

    if (twsolve >= 0)
    {
	if (twtot >= 0 && ftot >= 0)
	{
	    twtot += tsolve ;
	    PRINT_INFO (
		"    total symbolic+numeric+solve wall clock time:  %.2f\n",
		twtot) ;
	    if (ftot > 0 && twtot > 0)
	    {
		PRINT_INFO (
		"    total symbolic+numeric+solve mflops(wallclock) %.2f\n",
		1e-6 * ftot / twtot) ;
	    }
	}
    }
    PRINTF (("\n")) ;
}
