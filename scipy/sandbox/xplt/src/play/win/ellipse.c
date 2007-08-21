/*
 * ellipse.c -- $Id$
 * p_ellipse for MS Windows
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playw.h"

void
p_ellipse(p_win *w, int x0, int y0, int x1, int y1, int border)
{
  HDC dc = w_getdc(w, border? 18 : 12);
  if (dc)
    Ellipse(dc, x0, y0, x1+1+(!border), y1+1+(!border));
}
