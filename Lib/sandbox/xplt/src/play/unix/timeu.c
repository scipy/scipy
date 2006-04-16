/*
 * timeu.c -- $Id$
 * p_cpu_secs for UNIX machines (ANSI C version not UNIX specific)
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "play.h"

/* if getrusage exists, it is likely best choice (eg- Digital UNIX) */
#ifdef USE_GETRUSAGE

/* this is BSD way to get user and system time */
#undef _POSIX_SOURCE
#include <sys/time.h>
#include <sys/resource.h>
double
p_cpu_secs(double *sys)
{
  struct rusage cpu;
  getrusage(RUSAGE_SELF, &cpu);
  if (sys) *sys = cpu.ru_stime.tv_sec + 1.0e-6*cpu.ru_stime.tv_usec;
  return cpu.ru_utime.tv_sec + 1.0e-6*cpu.ru_utime.tv_usec;
}

#else
# ifdef USE_TIMES

/* this is POSIX 1003.1-1990 standard timing interface */
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#include <time.h>
#include <sys/times.h>
/* try to handle modest deviations from POSIX standard (e.g.- Sun) */
#  ifndef CLK_TCK
#   include <unistd.h>
#   ifndef CLK_TCK
#    define CLK_TCK sysconf(_SC_CLK_TCK)
#   endif
#  endif
static double secs_per_tick = 0.0;
double
p_cpu_secs(double *sys)
{
  struct tms cpu;
  times(&cpu);
  if (secs_per_tick==0.0) secs_per_tick = 1./((double)CLK_TCK);
  if (sys) *sys = cpu.tms_stime*secs_per_tick;
  return cpu.tms_utime*secs_per_tick;
}

# else

/* ANSI C standard should at least compile anywhere */
#include "time.h"
static double secs_per_tick = 0.0;
double
p_cpu_secs(double *sys)
{
  if (secs_per_tick==0.0) secs_per_tick = 1./((double)CLOCKS_PER_SEC);
  if (sys) *sys = 0.0;
  return clock()*secs_per_tick;
}

# endif
#endif
