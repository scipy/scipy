/*
 * colors.c -- $Id$
 * color handling for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

GC
x_getgc(p_scr *s, p_win *w, int fillstyle)
{
  GC gc = s->gc;
  if (w && w != s->gc_w_clip) {
    x_clip(s->xdpy->dpy, gc,
           w->xyclip[0], w->xyclip[1], w->xyclip[2], w->xyclip[3]);
    s->gc_w_clip = w;
  }
  if (fillstyle != s->gc_fillstyle) {
    /* note: if gc_color == P_GRAYB--P_GRAYC, then gc_fillstyle
     *       will not reflect the actual fillstyle in gc
     * the fillstyle would be unnecessary but for rotated text,
     * which requires stippled fills, and therefore may not work
     * properly with P_GRAYB--P_GRAYC colors -- see textout.c */
    XSetFillStyle(s->xdpy->dpy, gc, fillstyle);
    s->gc_fillstyle = fillstyle;
  }
  return gc;
}

void
p_color(p_win *w, p_col_t color)
{
  p_scr *s = w->s;
  GC gc = s->gc;
  p_col_t old_color = s->gc_color;

  if (color != old_color) {
    p_col_t pixel;

    s->gc_color = -1; /* invalidate; also when w->pixels changes */
    pixel = x_getpixel(w, color);

    if (color==P_XOR) XSetFunction(s->xdpy->dpy, gc, GXxor);
    else if (old_color==P_XOR ||
             old_color==-1) XSetFunction(s->xdpy->dpy, gc, GXcopy);

    if ((color==P_GRAYC || color==P_GRAYB) && s->gui_flags) {
      XSetFillStyle(s->xdpy->dpy, gc, FillOpaqueStippled);
      XSetStipple(s->xdpy->dpy, gc, s->gray);
      XSetBackground(s->xdpy->dpy, gc, s->colors[3].pixel);
      /* note that s->gc_fillstyle cannot be set here */
    } else if ((old_color==P_GRAYC || old_color==P_GRAYB) && s->gui_flags) {
      XSetFillStyle(s->xdpy->dpy, gc, FillSolid);
      XSetBackground(s->xdpy->dpy, gc, s->colors[0].pixel);
    }

    XSetForeground(s->xdpy->dpy, gc, pixel);
    s->gc_color = color;
  }
}

p_col_t
x_getpixel(p_win *w, p_col_t color)
{
  p_scr *s = w->s;
  p_col_t pixel;
  if (w->parent) w = w->parent;  /* offscreens use parent pixels */
  if (!P_IS_RGB(color)) {       /* standard and indexed color models */
    pixel = w->pixels[color];

  } else {                      /* rgb color model */
    unsigned int r = P_R(color);
    unsigned int g = P_G(color);
    unsigned int b = P_B(color);
    if (s->vclass==TrueColor || s->vclass==DirectColor) {
      pixel = (s->pixels[r]&s->rmask) |
        (s->pixels[g]&s->gmask) | (s->pixels[b]&s->bmask);
    } else if (s->vclass!=PseudoColor) {
      pixel = s->pixels[(r+g+b)/3];
    } else if (w->rgb_pixels || x_rgb_palette(w)) {
      /* use precomputed 5-9-5 color brick */
      r = (r+32)>>6;
      g = (g+16)>>5;
      b = (b+32)>>6;
      g += b+(b<<3);
      pixel = w->rgb_pixels[r+g+(g<<2)];  /* r + 5*g * 45*b */
    } else {
      pixel = s->colors[1].pixel;
    }
  }
  return pixel;
}
