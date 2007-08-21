/*
 * points.m
 * p_i_pnts, p_d_pnts, p_d_map, p_lines, p_segments, p_dots, p_rect, and p_fill
 * for Mac OS X.
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "playm.h"

static CGPoint m_pt_list[2050];
static int m_pt_count = 0;
static double m_pt_xa=1., m_pt_xb=0., m_pt_ya=1., m_pt_yb=0.;

void
p_d_map(p_win *w, double xt[], double yt[], int set)
{
  if (set) {
    m_pt_xa = xt[0];
    m_pt_xb = xt[1];
    m_pt_ya = yt[0];
    m_pt_yb = yt[1];
  } else {
    xt[0] = m_pt_xa;
    xt[1] = m_pt_xb;
    yt[0] = m_pt_ya;
    yt[1] = m_pt_yb;
  }
}

void
p_i_pnts(p_win *w, const int *x, const int *y, int n)
{
  if (n == -1) {
    if (m_pt_count < 2048) {
      n = m_pt_count++;
      m_pt_list[n].x = x[0];
      m_pt_list[n].y = y[0];
    } else {
      m_pt_count = 0;
    }
  } else {
    CGPoint *wrk = m_pt_list;
    if (n >= 0) {
      m_pt_count = n;
    } else {
      wrk += m_pt_count;
      n = -n;
      m_pt_count += n;
    }
    if (m_pt_count <= 2048) {
      while (n--) {
        wrk[0].x = *x++;
        wrk[0].y = *y++;
        wrk++;
      }
    } else {
      m_pt_count = 0;
    }
  }
}

void
p_d_pnts(p_win *w, const double *x, const double *y, int n)
{
  if (n == -1) {
    if (m_pt_count < 2048) {
      n = m_pt_count++;
      m_pt_list[n].x = (long)(m_pt_xa*x[0] + m_pt_xb);
      m_pt_list[n].y = (long)(m_pt_ya*y[0] + m_pt_yb);
    } else {
      m_pt_count = 0;
    }
  } else {
    CGPoint *wrk = m_pt_list;
    if (n >= 0) {
      m_pt_count = n;
    } else {
      wrk += m_pt_count;
      n = -n;
      m_pt_count += n;
    }
    if (m_pt_count <= 2048) {
      while (n--) {
        wrk[0].x = (long)(m_pt_xa*(*x++) + m_pt_xb);
        wrk[0].y = (long)(m_pt_ya*(*y++) + m_pt_yb);
        wrk++;
      }
    } else {
      m_pt_count = 0;
    }
  }
}

void
p_lines(p_win *w)
{
  const int n = m_pt_count;
  View* view = w->view;
  CGContextRef cr = w->cr;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view && n>0) {
    CGPoint* pxy = m_pt_list;
    CGContextBeginPath(cr);
    if (n>1) {
      const int closed = (pxy[0].x==pxy[n-1].x && pxy[0].y==pxy[n-1].y);
      if (closed)
      { CGContextAddLines(cr, pxy, n-1);
        CGContextClosePath(cr);
      }
      else CGContextAddLines(cr, pxy, n);
    } else {
      CGContextMoveToPoint(cr, pxy[0].x, pxy[0].y);
      CGContextAddLineToPoint(cr, pxy[0].x+1,pxy[0].y);
    }
    CGContextStrokePath(cr);
  }
}

void
p_segments(p_win *w)
{
  int n = m_pt_count / 2;
  View* view = w->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view && n>0) {
    CGContextRef cr = w->cr;
    CGPoint* pxy = m_pt_list;
    CGContextBeginPath(cr);
    while (n-- > 0) {
      CGContextMoveToPoint(cr, pxy[0].x, pxy[0].y);
      if (pxy[1].x!=pxy[0].x || pxy[1].y!=pxy[0].y) {
        CGContextAddLineToPoint(cr, pxy[1].x,pxy[1].y);
      } else {
        CGContextAddLineToPoint(cr, pxy[1].x+1,pxy[1].y);
      }
      pxy += 2;
    }
    CGContextStrokePath(cr);
  }
}

void
p_dots(p_win *w)
{
  int n = m_pt_count;
  View* view = w->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view) {
    CGContextRef cr = w->cr;
    CGPoint *pxy = m_pt_list;
    CGRect rect;
    rect.size.width = 1.0;
    rect.size.height = 1.0;
    while (n-- > 0) {
      rect.origin.x = pxy->x - 0.5;
      rect.origin.y = pxy->y - 0.5;
      CGContextFillRect(cr, rect);
      pxy++;
      /* Due to anti-aliasing, this effectively fills four pixels
       * in a two-by-two square. */
    }
  }
}

void
p_rect(p_win *w, int x0, int y0, int x1, int y1, int border)
{
  View* view = w->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view) {
    CGContextRef cr = w->cr;
    CGRect r = CGRectMake(x0,y0,x1-x0,y1-y0);
    if (border) CGContextStrokeRect(cr, r);
    else CGContextFillRect(cr, r);
  }
}

void
p_fill(p_win *w, int convexity)
{
  View* view = w->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view) {
    CGContextRef cr = w->cr;
    CGContextAddLines(cr, m_pt_list, m_pt_count);
    CGContextClosePath(cr);
    CGContextFillPath(cr);
  }
}
