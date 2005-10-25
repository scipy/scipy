/* ========================================================================== */
/* === umfpack_tictoc ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Returns the time in seconds used by the process, and
    the current wall clock time.  BE CAREFUL:  if you compare the run time of
    UMFPACK with other sparse matrix packages, be sure to use the same timer.
    See umfpack_tictoc.h for details.

    These routines conform to the POSIX standard.  See umf_config.h for
    more details.
*/

#include "umf_internal.h"

#define TINY_TIME 1e-4

#ifndef NPOSIX

#include <unistd.h>
#include <sys/times.h>

void umfpack_tic (double stats [2])
{
    /* Return the current time */
    /* stats [0]: current wallclock time, in seconds */
    /* stats [1]: user + system time for the process, in seconds */

    struct tms t ;
    double ticks ;

    ticks = (double) sysconf (_SC_CLK_TCK) ;
    stats [0] = (double) times (&t) / ticks ;
    stats [1] = (double) (t.tms_utime + t.tms_stime) / ticks ;

    /* if time is tiny, just return zero */
    if (stats [0] < TINY_TIME) stats [0] = 0 ;
    if (stats [1] < TINY_TIME) stats [1] = 0 ;
}

#else

/* Generic ANSI C: use the ANSI clock function.  No wallclock time. */

#include <time.h>

void umfpack_tic (double stats [2])
{
    stats [0] = 0 ;
    stats [1] = ((double) (clock ( ))) / ((double) (CLOCKS_PER_SEC)) ;
    if (stats [1] < TINY_TIME) stats [1] = 0 ;
}

#endif

/* -------------------------------------------------------------------------- */

void umfpack_toc (double stats [2])
{
    /* Return the current time since the last call to umfpack_tic. */
    /* On input, stats holds the values returned by umfpack_tic. */
    /* On ouput, stats holds the time since the last umfpack_tic. */

    double done [2] ;
    umfpack_tic (done) ;

    stats [0] = done [0] - stats [0] ;
    stats [1] = done [1] - stats [1] ;

    if (stats [0] < 0) stats [0] = 0 ;
    if (stats [1] < 0) stats [1] = 0 ;

}
