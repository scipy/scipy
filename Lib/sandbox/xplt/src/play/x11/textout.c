/*
 * textout.c -- $Id$
 * p_text for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

#include "pstdlib.h"

#include <string.h>

int
p_txheight(p_scr *s, int font, int pixsize, int *baseline)
{
  XFontStruct *f = x_font(s->xdpy, font, pixsize);
  if (baseline) *baseline = f->ascent;
  return f->ascent + f->descent;  /* hopefully close to pixsize */
}

int
p_txwidth(p_scr *s, const char *text, int n, int font, int pixsize)
{
  XFontStruct *f = x_font(s->xdpy, font, pixsize);
  int len = strlen(text);
  if (n<=0 || n>len) n = len;
  return XTextWidth(f, (char *)text, n);
}

void
p_font(p_win *w, int font, int pixsize, int orient)
{
  p_scr *s = w->s;

  if (s->rotgc || s->tmp || s->image || s->pixmap!=None) x_rotzap(s);

  if (!orient) {
    s->rotgc_orient = 0;
    if (font!=s->gc_font || pixsize!=s->gc_pixsize) {
      XFontStruct *f = x_font(s->xdpy, font, pixsize);
      XSetFont(s->xdpy->dpy, s->gc, f->fid);
      s->gc_font = font;
      s->gc_pixsize = pixsize;
    }

  } else {
    /* defer rotated font selection to p_text call */
    s->rotgc_font = font;
    s->rotgc_pixsize = pixsize;
    s->rotgc_orient = orient;
  }
}

void
p_text(p_win *w, int x0, int y0, const char *text, int n)
{
  p_scr *s = w->s;
  x_display *xdpy = s->xdpy;
  Display *dpy = xdpy->dpy;
  int orient = s->rotgc_orient;
  GC gc = x_getgc(s, w, orient? FillStippled : FillSolid);
  Drawable d = w->d;
  int i;

  if (s->rotgc || s->tmp || s->image || s->pixmap!=None) x_rotzap(s);

  if (n<=0) n = 16350;
  for (i=0 ; i<n ; i++) if (!text[i]) break;
  n = i;

  if (!orient) {
    XDrawString(dpy, d, gc, x0, y0, (char *)text, n);

  } else {
    XFontStruct *f = x_font(xdpy, s->rotgc_font, s->rotgc_pixsize);
    int width = XTextWidth(f, (char *)text, n);
    int height = f->ascent + f->descent;
    XGCValues values;

    char *data = p_malloc(((width-1)/8+1)*height);
    if (!data) { x_rotzap(s);  return; }

    /* create a bitmap and a special rotgc,
     * draw text to bitmap and destroy the rotgc */
    s->pixmap = XCreatePixmap(dpy, s->root, width, height, 1);
    values.foreground = 1;
    values.background = 0;
    values.font = f->fid;
    s->rotgc = XCreateGC(dpy, s->pixmap,
                         GCForeground|GCBackground|GCFont, &values);
    XDrawImageString(dpy, s->pixmap, s->rotgc, 0, f->ascent, (char *)text, n);

    /* create an image in the format used by XCreateBitmapFromData */
    s->own_image_data = 1;
    s->image = XCreateImage(dpy, (Visual *)0, 1, XYBitmap, 0,
                            data, width, height, 8, 0);
    s->image->byte_order = s->image->bitmap_bit_order = LSBFirst;

    /* read the bitmap containing the text, then destroy it
     * - in current XFree86 source, this is spectacularly inefficient
     *   using XGetPixel and XPutPixel
     * - XPutImage is currently the only efficient image routine, but
     *   the labor of dealing with arbitrary bitmaps is more fairly
     *   theirs than mine; hopefully by selecting the image format
     *   used by their own bitmap data files, I get optimized first */
    XGetSubImage(dpy, s->pixmap, 0,0, width, height, 1, XYPixmap,
                 s->image, 0,0);
    x_pxzap(dpy, &s->pixmap);
    s->image->data = 0;
    s->tmp = data;

    if (orient==2) {
      s->image->data = p_malloc(((width-1)/8+1)*height);
      if (!s->image->data) { x_rotzap(s);  return; }
      s->pixmap = XCreatePixmap(dpy, s->root, width, height, 1);
      p_lrot180((void *)data, (void *)s->image->data, width, height);
      XPutImage(dpy, s->pixmap, s->rotgc, s->image, 0,0, 0,0, width,height);
      x0 -= width-1;
      y0 -= f->descent-1;
      /* do switcharoo for this case, since upside down is unusual */
      i=width; width=height; height=i;

    } else {
      /* stupid to destroy image, then recreate a new one, but the
       * interface is botched to make it tedious to recompute all
       * the image members that change with the shape transpose */
      x_imzap(s);
      s->image = XCreateImage(dpy, (Visual *)0, 1, XYBitmap, 0,
                              p_malloc(((height-1)/8+1)*width),
                              height, width, 8, 0);
      if (!s->image->data) { x_rotzap(s);  return; }
      s->image->byte_order = s->image->bitmap_bit_order = LSBFirst;
      s->pixmap = XCreatePixmap(dpy, s->root, height, width, 1);
      if (orient==1) {
        p_lrot090((void *)data, (void *)s->image->data, width, height);
        x0 -= f->ascent;
        y0 -= width-1;
      } else {
        p_lrot270((void *)data, (void *)s->image->data, width, height);
        x0 -= f->descent-1;
      }
      XPutImage(dpy, s->pixmap, s->rotgc, s->image, 0,0, 0,0, height,width);
    }

    /* stipple fillstyle already set by x_getgc,
     * draw the text as a stippled rectangle */
    XSetStipple(dpy, gc, s->pixmap);
    XSetTSOrigin(dpy, gc, x0, y0);
    XFillRectangle(dpy, d, gc, x0, y0, height, width);

    x_rotzap(s);
  }
  if (p_signalling) p_abort();
}
