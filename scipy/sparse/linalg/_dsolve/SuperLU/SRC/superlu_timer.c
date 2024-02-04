/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file superlu_timer.c
 * \brief Returns the time used
 *
 * <pre>
 * Purpose
 * ======= 
 * 
 * Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 * </pre>
 */

#undef USE_TIMES

#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>

double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#elif _WIN32

#include <time.h>

double SuperLU_timer_()
{
    clock_t t;
    t=clock();

    return ((double)t)/CLOCKS_PER_SEC;
}

#elif defined( USE_TIMES )

#ifndef NO_TIMER
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>
#endif

/*! \brief Timer function
 */ 
double SuperLU_timer_()
{
#ifdef NO_TIMER
    /* no sys/times.h on WIN32 */
    double tmp;
    tmp = 0.0;
#else
    struct tms use;
    double tmp;
    int clocks_per_sec = sysconf(_SC_CLK_TCK);

    times ( &use );
    tmp = use.tms_utime;
    tmp += use.tms_stime;
#endif
    return (double)(tmp) / clocks_per_sec;
}

#else

#include <sys/time.h>
#include <stdlib.h>

/*! \brief Timer function
 */ 
double SuperLU_timer_(void)
{
    struct timeval tp;
    double tmp;
    gettimeofday(&tp, NULL);
    tmp = tp.tv_sec + tp.tv_usec/1000000.0;
    return (tmp);
}

#endif

