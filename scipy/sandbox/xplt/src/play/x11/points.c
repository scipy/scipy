/*
 * points.c -- $Id$
 * p_i_pnts, p_d_pnts, p_d_map for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

XPoint x_pt_list[2050];
int x_pt_count = 0;
static double x_pt_xa=1., x_pt_xb=0., x_pt_ya=1., x_pt_yb=0.;

/* ARGSUSED */
void
p_d_map(p_win *w, double xt[], double yt[], int set)
{
  if (set) {
    x_pt_xa = xt[0];
    x_pt_xb = xt[1];
    x_pt_ya = yt[0];
    x_pt_yb = yt[1];
  } else {
    xt[0] = x_pt_xa;
    xt[1] = x_pt_xb;
    yt[0] = x_pt_ya;
    yt[1] = x_pt_yb;
  }
}

void
p_i_pnts(p_win *w, const int *x, const int *y, int n)
{
  if (n == -1) {
    if (x_pt_count < 2048) {
      n = x_pt_count++;
      x_pt_list[n].x = x[0];
      x_pt_list[n].y = y[0];
    } else {
      x_pt_count = 0;
    }
  } else {
    XPoint *wrk = x_pt_list;
    if (n >= 0) {
      x_pt_count = n;
    } else {
      wrk += x_pt_count;
      n = -n;
      x_pt_count += n;
    }
    if (x_pt_count <= 2048) {
      while (n--) {
        wrk[0].x = *x++;
        wrk[0].y = *y++;
        wrk++;
      }
    } else {
      x_pt_count = 0;
    }
  }
}

/* ARGSUSED */
void
p_d_pnts(p_win *w, const double *x, const double *y, int n)
{
  if (n == -1) {
    if (x_pt_count < 2048) {
      n = x_pt_count++;
      x_pt_list[n].x = (short)(x_pt_xa*x[0] + x_pt_xb);
      x_pt_list[n].y = (short)(x_pt_ya*y[0] + x_pt_yb);
    } else {
      x_pt_count = 0;
    }
  } else {
    XPoint *wrk = x_pt_list;
    if (n >= 0) {
      x_pt_count = n;
    } else {
      wrk += x_pt_count;
      n = -n;
      x_pt_count += n;
    }
    if (x_pt_count <= 2048) {
      while (n--) {
        wrk[0].x = (short)(x_pt_xa*(*x++) + x_pt_xb);
        wrk[0].y = (short)(x_pt_ya*(*y++) + x_pt_yb);
        wrk++;
      }
    } else {
      x_pt_count = 0;
    }
  }
}
