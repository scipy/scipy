/*
 * feep.c -- $Id$
 * p_feep for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

void
p_feep(p_win *w)
{
  if (w->s->xdpy->dpy) XBell(w->s->xdpy->dpy, 100);
}
