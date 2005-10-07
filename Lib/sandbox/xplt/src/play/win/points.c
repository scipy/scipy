/*
 * points.c -- $Id$
 * p_i_pnts, p_d_pnts, p_d_map for MS Windows
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "playw.h"

POINT w_pt_list[2050];
int w_pt_count = 0;
static double w_pt_xa=1., w_pt_xb=0., w_pt_ya=1., w_pt_yb=0.;

/* ARGSUSED */
void
p_d_map(p_win *w, double xt[], double yt[], int set)
{
  if (set) {
    w_pt_xa = xt[0];
    w_pt_xb = xt[1];
    w_pt_ya = yt[0];
    w_pt_yb = yt[1];
  } else {
    xt[0] = w_pt_xa;
    xt[1] = w_pt_xb;
    yt[0] = w_pt_ya;
    yt[1] = w_pt_yb;
  }
}

void
p_i_pnts(p_win *w, const int *x, const int *y, int n)
{
  if (n == -1) {
    if (w_pt_count < 2048) {
      n = w_pt_count++;
      w_pt_list[n].x = x[0];
      w_pt_list[n].y = y[0];
    } else {
      w_pt_count = 0;
    }
  } else {
    POINT *wrk = w_pt_list;
    if (n >= 0) {
      w_pt_count = n;
    } else {
      wrk += w_pt_count;
      n = -n;
      w_pt_count += n;
    }
    if (w_pt_count <= 2048) {
      while (n--) {
        wrk[0].x = *x++;
        wrk[0].y = *y++;
        wrk++;
      }
    } else {
      w_pt_count = 0;
    }
  }
}

/* ARGSUSED */
void
p_d_pnts(p_win *w, const double *x, const double *y, int n)
{
  if (n == -1) {
    if (w_pt_count < 2048) {
      n = w_pt_count++;
      w_pt_list[n].x = (long)(w_pt_xa*x[0] + w_pt_xb);
      w_pt_list[n].y = (long)(w_pt_ya*y[0] + w_pt_yb);
    } else {
      w_pt_count = 0;
    }
  } else {
    POINT *wrk = w_pt_list;
    if (n >= 0) {
      w_pt_count = n;
    } else {
      wrk += w_pt_count;
      n = -n;
      w_pt_count += n;
    }
    if (w_pt_count <= 2048) {
      while (n--) {
        wrk[0].x = (long)(w_pt_xa*(*x++) + w_pt_xb);
        wrk[0].y = (long)(w_pt_ya*(*y++) + w_pt_yb);
        wrk++;
      }
    } else {
      w_pt_count = 0;
    }
  }
}
