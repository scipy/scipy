/*
 * rgbread.c -- $Id$
 * p_rgbread for X11, reads window contents into pixmap
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

#include <X11/Xutil.h>

static int rgb_find_shift(unsigned long mask);

void
p_rgb_read(p_win *w, unsigned char *rgbs, int x0, int y0, int x1, int y1)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  XImage *image;
  unsigned long pixel;

  if (s->image) x_imzap(s);
  x1 -= x0;
  y1 -= y0;
  s->own_image_data = 0;
  image = s->image = XGetImage(dpy, w->d, x0, y0, x1, y1, AllPlanes, ZPixmap);

  if (s->vclass==TrueColor || s->vclass==DirectColor) {
    unsigned long rmask = s->rmask;
    unsigned long gmask = s->gmask;
    unsigned long bmask = s->bmask;
    int rshr = rgb_find_shift(rmask);
    int gshr = rgb_find_shift(gmask);
    int bshr = rgb_find_shift(bmask);
    for (y0=0 ; y0<y1 ; y0++)
      for (x0=0 ; x0<x1 ; x0++) {
        pixel = XGetPixel(image, x0, y0);
        rgbs[0] = rshr>=0? ((pixel&rmask) >> rshr) : ((pixel&rmask) << -rshr);
        rgbs[1] = gshr>=0? ((pixel&gmask) >> gshr) : ((pixel&gmask) << -gshr);
        rgbs[2] = bshr>=0? ((pixel&bmask) >> bshr) : ((pixel&bmask) << -bshr);
        rgbs += 3;
      }

  } else {
    Colormap cmap = (w->cmap==None? DefaultColormap(dpy,s->scr_num) : w->cmap);
    int map_size = DefaultVisual(dpy,s->scr_num)->map_entries;
    XColor map[256];
    if (map_size>256) map_size = 256;
    for (y0=0 ; y0<map_size ; y0++) map[y0].pixel = y0;
    XQueryColors(dpy, cmap, map, map_size);
    for (y0=0 ; y0<y1 ; y0++)
      for (x0=0 ; x0<x1 ; x0++) {
        pixel = XGetPixel(image, x0, y0);
        if (pixel<256) {
          rgbs[0] = map[pixel].red>>8;
          rgbs[1] = map[pixel].green>>8;
          rgbs[2] = map[pixel].blue>>8;
        } else {
          rgbs[0] = rgbs[1] = rgbs[2] = 0;
        }
        rgbs += 3;
      }
  }

  x_imzap(s);
  if (p_signalling) p_abort();
}

static int
rgb_find_shift(unsigned long mask)
{
  int s;
  for (s=0 ; s<8*sizeof(mask)-1 ; s++) if ((mask>>s)&1) break;
  for (s++ ; s<8*sizeof(mask) ; s++) if (!((mask>>s)&1)) break;
  return s-8;
}
