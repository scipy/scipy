/*
 * timew.c -- $Id$
 * p_wall_secs for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "play.h"

/* GetTickCount in kernel32.lib */
#include <windows.h>

static int p_wall_init = 0;

static double p_wall_0 = 0.0;
static DWORD p_wall_1 = 0;
double
p_wall_secs(void)
{
  DWORD wall = GetTickCount();
  if (!p_wall_init) {
    p_wall_1 = wall;
    p_wall_init = 1;
  } else if (wall < p_wall_1) {
    p_wall_0 += 1.e-3*((~(DWORD)0) - p_wall_1);
    p_wall_0 += 1.e-3;
    p_wall_1 = 0;
  }
  p_wall_0 += 1.e-3*(wall-p_wall_1);
  p_wall_1 = wall;
  return p_wall_0;
}
