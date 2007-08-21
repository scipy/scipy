/*
 * rect.c -- $Id$
 * p_rect for X11
 *
 * Copyright (c) 2001.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

void
p_rect(p_win *w, int x0, int y0, int x1, int y1, int border)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  GC gc = x_getgc(s, w, FillSolid);
  int tmp;
  if (x1 > x0) x1 -= x0;
  else tmp = x0-x1, x0 = x1, x1 = tmp;
  if (y1 > y0) y1 -= y0;
  else tmp = y0-y1, y0 = y1, y1 = tmp;
  if (border)
    XDrawRectangle(dpy, w->d, gc, x0, y0, x1, y1);
  else
    XFillRectangle(dpy, w->d, gc, x0, y0, x1, y1);
  if (p_signalling) p_abort();
}
