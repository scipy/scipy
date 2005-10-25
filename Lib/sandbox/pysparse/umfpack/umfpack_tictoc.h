/* ========================================================================== */
/* === umfpack_tictoc ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

void umfpack_tic (double stats [2]) ;

void umfpack_toc (double stats [2]) ;


/*
Syntax (for all versions: di, dl, zi, and zl):

    #include "umfpack.h"
    double stats [2] ;
    umfpack_tic (stats) ;
    ...
    umfpack_toc (stats) ;

Purpose:

    umfpack_tic returns the CPU time and wall clock time used by the process.
    The CPU time includes both "user" and "system" time (the latter is time
    spent by the system on behalf of the process, and is thus charged to the
    process).  umfpack_toc returns the CPU time and wall clock time since the
    last call to umfpack_tic with the same stats array.

    Typical usage:

	umfpack_tic (stats) ;
	... do some work ...
	umfpack_toc (stats) ;

    then stats [1] contains the time in seconds used by the code between
    umfpack_tic and umfpack_toc, and stats [0] contains the wall clock time
    elapsed between the umfpack_tic and umfpack_toc.  These two routines act
    just like tic and toc in MATLAB, except that the both process time and
    wall clock time are returned.

    This routine normally uses the sysconf and times routines in the POSIX
    standard.  If -DNPOSIX is defined at compile time, then the ANSI C clock
    routine is used instead, and only the CPU time is returned (stats [0]
    is set to zero).

    umfpack_tic and umfpack_toc are the routines used internally in UMFPACK
    to time the symbolic analysis, numerical factorization, and the forward/
    backward solve.

Arguments:

    double stats [2]:

	stats [0]:  wall clock time, in seconds
	stats [1]:  CPU time, in seconds
*/
