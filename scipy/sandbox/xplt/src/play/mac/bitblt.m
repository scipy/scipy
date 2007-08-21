/*
 * bitblt.m
 * p_bitblt for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

void
p_bitblt(p_win *w, int x, int y, p_win *offscreen,
         int x0, int y0, int x1, int y1)
{
  View* view = w->view;
  View* source = offscreen->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view)
  { [view lockFocus];
    if (source)
    { int gs = [source gState];
      if (gs) NSCopyBits(gs, NSMakeRect(x0,y0,x1-x0,y1-y0), NSMakePoint(x,y+y1-y0));
    }
    [view unlockFocus];
    [w->w flushWindow];
  }
}
