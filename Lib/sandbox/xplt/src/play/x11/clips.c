/*
 * clips.c -- $Id$
 * p_clip for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

void
p_clip(p_win *w, int x0, int y0, int x1, int y1)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  GC gc = x_getgc(s, (p_win *)0, FillSolid);
  w->xyclip[0] = x0;
  w->xyclip[1] = y0;
  w->xyclip[2] = x1;
  w->xyclip[3] = y1;
  x_clip(dpy, gc, x0, y0, x1, y1);
  s->gc_w_clip = w;
}

void
x_clip(Display *dpy, GC gc, int x0, int y0, int x1, int y1)
{
  XRectangle xr;
  if (x1>x0 && y1>y0) {
    xr.width = x1 - x0;
    xr.height = y1 - y0;
    xr.x = x0;
    xr.y = y0;
    XSetClipRectangles(dpy, gc, 0,0, &xr, 1, YXBanded);
  } else {
    XSetClipMask(dpy, gc, None);
  }
  if (p_signalling) p_abort();
}
