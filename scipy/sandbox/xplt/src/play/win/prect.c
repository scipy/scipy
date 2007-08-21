/*
 * prect.c -- $Id$
 * p_rect for MS Windows
 *
 * Copyright (c) 2001.  See accompanying LEGAL file for details.
 */

#include "playw.h"

void
p_rect(p_win *w, int x0, int y0, int x1, int y1, int border)
{
  if (border) {
    HDC dc = w_getdc(w, 18);
    if (dc)
      Rectangle(dc, x0, y0, x1+1, y1+1);
  } else {
    HDC dc = w_getdc(w, 4);
    if (dc && w->brush) {
      RECT r;
      r.left = x0;
      r.top = y0;
      r.right = x1;
      r.bottom = y1;
      FillRect(dc, &r, w->brush);
    }
  }
}
