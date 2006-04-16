/*
 * lines.c -- $Id$
 * p_lines and p_disjoint for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

void
p_dots(p_win *w)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  GC gc = x_getgc(s, w, FillSolid);
  int nmx = XMaxRequestSize(dpy)-3;
  int n = x_pt_count;
  x_pt_count = 0;
  while (n>0) {
    if (n<nmx) nmx = n;
    XDrawPoints(dpy, w->d, gc, x_pt_list, nmx, CoordModeOrigin);
    n -= nmx;
  }
  if (p_signalling) p_abort();
}

void
p_lines(p_win *w)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  GC gc = x_getgc(s, w, FillSolid);
  int nmx = XMaxRequestSize(dpy)-3;
  int n = x_pt_count;
  x_pt_count = 0;
  while (n>1) {
    if (n<nmx) nmx = n;
    XDrawLines(dpy, w->d, gc, x_pt_list, nmx, CoordModeOrigin);
    n -= nmx;
  }
  if (p_signalling) p_abort();
}

void
p_segments(p_win *w)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  GC gc = x_getgc(s, w, FillSolid);
  int nmx = (XMaxRequestSize(dpy)-3)/2;
  int n = x_pt_count / 2;
  x_pt_count = 0;
  while (n>0) {
    if (n<nmx) nmx = n;
    /* note: assume here that XPoint[2] identical to XSegment */
    XDrawSegments(dpy, w->d, gc, (XSegment *)x_pt_list, nmx);
    n -= nmx;
  }
  if (p_signalling) p_abort();
}

static char dashed[] = { 5, 5 };
static char dotted[] = { 1, 3 };
static char dashdot[] = { 5, 2, 1, 2 };
static char dashdotdot[] = { 5, 2, 1, 2, 1, 2 };
static char *x_dash[] = { 0, dashed, dotted, dashdot, dashdotdot };
static int x_ndash[] = { 0, 2, 2, 4, 6 };

void
p_pen(p_win *w, int width, int type)
{
  p_scr *s = w->s;
  GC gc = s->gc;
  int disjoint = (type & P_SQUARE);
  int same_type = (s->gc_type == type);

  if (width<2) width = 0;
  else if (width>100) width = 100;

  if (s->gc_width==width && same_type) return;

  type ^= disjoint;
  if (type>4 || type<0) type = 0;
  XSetLineAttributes(s->xdpy->dpy, gc, width,
                     type? LineOnOffDash : LineSolid,
                     disjoint? CapProjecting : CapRound,
                     disjoint? JoinMiter : JoinRound);
  if (!same_type) s->gc_type = (type | disjoint);
  s->gc_width = width;

  if (type) {
    /* dash pattern depends on linestyle */
    int n = x_ndash[type];
    if (width<2) {
      XSetDashes(s->xdpy->dpy, gc, 0, x_dash[type], n);
    } else {
      /* dash pattern must scale with line thickness */
      int i;
      char dash[6];
      for (i=0 ; i<n ; i++)
        dash[i] = x_dash[type][i]>1? width*x_dash[type][i] : 1;
      XSetDashes(s->xdpy->dpy, gc, 0, dash, n);
    }
  }
}
