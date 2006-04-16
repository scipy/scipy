/*
 * clips.c -- $Id$
 * p_clip for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"

void
p_clip(p_win *w, int x0, int y0, int x1, int y1)
{
  HDC dc = w_getdc(w, 0);
  if (dc) {
    HRGN rgn = (x1<=x0 || y1<=y0)? 0 : CreateRectRgn(x0, y0, x1, y1);
    SelectClipRgn(dc, rgn);
    if (rgn) DeleteObject(rgn);
  }
}
