/*
 * pen.m
 * Choosing line width and style for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

static float dashed[] = { 5, 5 };
static float dotted[] = { 1, 3 };
static float dashdot[] = { 5, 2, 1, 2 };
static float dashdotdot[] = { 5, 2, 1, 2, 1, 2 };
static float *x_dash[] = { 0, dashed, dotted, dashdot, dashdotdot };
static int x_ndash[] = { 0, 2, 2, 4, 6 };

void
p_pen(p_win *w, int width, int type)
{ CGContextRef cr = w->cr;
  int disjoint = (type & P_SQUARE);
  type ^= disjoint;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (width<1) width = 1;
  else if (width>100) width = 100;
  CGContextSetLineWidth(cr, (float)width);
  if (type != P_SOLID)
  { int n = x_ndash[type];
    if (width<2)
      CGContextSetLineDash(cr, 0.0, x_dash[type], n);
    else
    { int i;
      float dash[6];
      for (i=0; i<n; i++)
      { float x = x_dash[type][i];
        dash[i] = x > 1.0 ? width * x : 1;
      }
      CGContextSetLineDash(cr, 0.0, dash, n);
    }
  }
  if (disjoint)
  { CGContextSetLineCap(cr, kCGLineCapSquare);
    CGContextSetLineJoin(cr, kCGLineJoinMiter);
  }
  else
  { CGContextSetLineCap(cr, kCGLineCapRound);
    CGContextSetLineJoin(cr, kCGLineJoinRound);
  }
}
