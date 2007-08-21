/*
 * timeu.c -- $Id$
 * p_cpu_secs for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "play.h"

/* functions in kernel32.lib */
#include <windows.h>

/* note: Windows NT provides a GetProcessTimes function that
 * returns system time as well as user time (also GetCurrentProcess) */

static double p_cpu_unit = 0.;
double
p_cpu_secs(double *sys)
{
  LARGE_INTEGER t;
  if (sys) *sys = 0.;
  if (p_cpu_unit == 0.) {
    if (QueryPerformanceFrequency(&t))
      p_cpu_unit = 1./(t.LowPart + 4294967296.*t.HighPart);
    else
      p_cpu_unit = -1.;
  }
  if (p_cpu_unit == -1. ||
      !QueryPerformanceCounter(&t))
    return p_wall_secs();
  return p_cpu_unit*(t.LowPart + 4294967296.*t.HighPart);
}
