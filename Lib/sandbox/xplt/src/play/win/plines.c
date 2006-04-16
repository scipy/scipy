/*
 * plines.c -- $Id$
 * p_lines for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"

/* Following advantages of paths would be nice, BUT paths cost
 *   a factor of two or so in performance, so forget it:
 * Reasons to use paths for all lines and fills:
 * (1) Win9x join and cap styles only work for paths
 * (2) Other fill primitives use pen, requiring it to be set NULL
 *     affects prect.c, ellipse.c
 */

/* MS Windows unclear on what happens when begin and end points
 * of a line are identical, making the inconsistent claim that
 * the first point on a line is always included and the last point
 * on a line is always excluded
 * -- experimentally, the last point can be included by adding
 *    a single segment one pixel long to the end (last segment
 *    cannot be zero pixels long) -- Munro's Win95 box
 * -- but on Langer's WinNT box, last point *is* included,
 * the two boxes also differ in does_linetypes, so the hack here
 * uses that to distinguish -- I dont know what else to do... */

static int w_polyline_to(HDC dc, int type, int width, POINT *pxy, int n,
                         int closed);
static void w_final_pt(HDC dc, POINT *pxy);
static float w_hypot(float dx, float dy);

void
p_lines(p_win *w)
{
  int n = w_pt_count;
  HDC dc = w_getdc(w, 2);
  if (dc && n>0) {
    POINT *pxy = w_pt_list;
    if (n>1) {
      int closed = (pxy[0].x==pxy[n-1].x && pxy[0].y==pxy[n-1].y);
      if ((w->pen_type&7)==P_SOLID || w->s->does_linetypes) {
        /* NOTE: PolylineTo is MUCH slower than Polyline */
        Polyline(dc, pxy, n);
        if (!closed && !w->s->does_linetypes)
           MoveToEx(dc, pxy[n-1].x, pxy[n-1].y, 0),
             w_final_pt(dc, &pxy[n-2]);
      } else {
        MoveToEx(dc, pxy[0].x, pxy[0].y, 0);
        w_polyline_to(dc, w->pen_type&7, w->pen_width, &pxy[0], n-1, closed);
      }
    } else {
      MoveToEx(dc, pxy[0].x, pxy[0].y, 0);
      LineTo(dc, pxy[0].x+1, pxy[0].y);
    }
  }
}

void
p_segments(p_win *w)
{
  int n = w_pt_count / 2;
  HDC dc = w_getdc(w, 2);
  if (dc && n>0) {
    POINT *pxy = w_pt_list;
    while (n-- > 0) {
      MoveToEx(dc, pxy[0].x, pxy[0].y, 0);
      if (pxy[1].x!=pxy[0].x || pxy[1].y!=pxy[0].y) {
        if ((w->pen_type&7)==P_SOLID || w->s->does_linetypes) {
          LineTo(dc, pxy[1].x, pxy[1].y);
          if (!w->s->does_linetypes) w_final_pt(dc, pxy);
        } else {
          w_polyline_to(dc, w->pen_type&7, w->pen_width, &pxy[0], 1, 0);
        }
      } else {
        LineTo(dc, pxy[1].x+1, pxy[1].y);
      }
      pxy += 2;
    }
  }
}

void
p_dots(p_win *w)
{
  int n = w_pt_count;
  HDC dc = w_getdc(w, 0);
  if (dc) {
    POINT *pxy = w_pt_list;
    COLORREF color = w_color(w, w->color);
    while (n-- > 0) {
      SetPixelV(dc, pxy[0].x, pxy[0].y, color);
      pxy++;
    }
  }
}

static void
w_final_pt(HDC dc, POINT *pxy)
{
  /* Windows maddeningly wont draw final pixel on a polyline
   * -- this hack works pretty well onscreen, worse in metafiles */
  long dx = pxy[1].x - pxy[0].x;
  long dy = pxy[1].y - pxy[0].y;
  int ix = 1, iy = 0;
  if (dx < 0) {
    if (dy < 0) {
      ix = -((dy > dx) || (dx < 2*dy));
      iy = -((dx > dy) || (dy < 2*dx));
    } else {
      ix = -((dy < -dx) || (-dx > 2*dy));
      iy = ((-dx < dy) || (dy > -2*dx));
    }
  } else {
    if (dy < 0) {
      ix = ((-dy < dx) || (dx > -2*dy));
      iy = -((dx < -dy) || (-dy > 2*dx));
    } else if (dy > 0) {
      ix = ((dy < dx) || (dx > 2*dy));
      iy = ((dx < dy) || (dy > 2*dx));
    }
  }
  LineTo(dc, pxy[1].x+ix, pxy[1].y+iy);
}

static float dashed[] = { 5, 5 };
static float dotted[] = { 1, 3 };
static float dashdot[] = { 5, 2, 1, 2 };
static float dashdotdot[] = { 5, 2, 1, 2, 1, 2 };
static float *x_dash[] = { 0, dashed, dotted, dashdot, dashdotdot };
static int x_ndash[] = { 0, 2, 2, 4, 6 };

static int
w_polyline_to(HDC dc, int type, int width, POINT *pxy, int n, int closed)
{
  int ndash = x_ndash[type];
  int i, x, y, xmv=0, ymv=0;
  int moved = 0;
  float dash[6], len, dx, dy, d=0.0;

  for (i=0 ; i<ndash ; i++) {
    dash[i] = x_dash[type][i];
    if (width>1 && dash[i]>1) dash[i] *= width;
  }

  i = 0;
  while (n-- > 0) {
    x = pxy[0].x;
    y = pxy[0].y;
    pxy++;
    dx = (float)(pxy[0].x - x);
    dy = (float)(pxy[0].y - y);
    len = w_hypot(dx, dy);
    if (d+len > dash[i]) {
      float s = -d;
      float m = 1.0f/len;
      float f;
      int xx, yy;
      do {
        s += dash[i];
        f = (s+0.5f)*m;
        xx = x + (int)(f*dx);
        yy = y + (int)(f*dy);
        if (i&1) {
          xmv = xx, ymv = yy, moved = 1;
        } else {
          if (moved) MoveToEx(dc, xmv, ymv, 0), moved = 0;
          LineTo(dc, xx, yy);
        }
        if ((++i) == ndash) i = 0;
      } while (s+dash[i] < len);
      d = len - s;  /* 0 < d <= dash[i] */
    } else {
      d += len;
    }
    if (i&1) {
      xmv = pxy[0].x, ymv = pxy[0].y, moved = 1;
    } else {
      if (moved) MoveToEx(dc, xmv, ymv, 0), moved = 0;
      LineTo(dc, pxy[0].x, pxy[0].y);
    }
  }
  if (!moved && !closed) LineTo(dc, pxy[0].x+1, pxy[0].y);
  return 1;
}

static float
w_hypot(float dx, float dy)
{
  if (dx<0.0) dx = -dx;
  if (dy<0.0) dy = -dy;
  if (dy > dx) {
    float tmp = dy;
    dy = dx;
    dx = tmp;
  }
  /* biggest error is when dx==dy -- off by 6 percent */
  if (dx) dx += 0.5f*dy*dy/dx;
  return dx;
}
