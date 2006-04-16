/*
 * images.c -- $Id$
 * p_pixmap and p_rgbmap for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

#include "pstdlib.h"

static void x_image(p_win *w, unsigned char *bytes, int rgb, int ncols,
                    int nrows, int x0, int y0, int x1, int y1);

static int x_put1(XImage *im, int i, int j, p_col_t pixel);
static int x_put2b(XImage *im, int i, int j, p_col_t pixel);
static int x_put2l(XImage *im, int i, int j, p_col_t pixel);
static int x_put3b(XImage *im, int i, int j, p_col_t pixel);
static int x_put3l(XImage *im, int i, int j, p_col_t pixel);
static int x_put4b(XImage *im, int i, int j, p_col_t pixel);
static int x_put4l(XImage *im, int i, int j, p_col_t pixel);

void
p_ndx_cell(p_win *w, unsigned char *ndxs, int ncols, int nrows,
           int x0, int y0, int x1, int y1)
{
  x_image(w, ndxs, 0, ncols, nrows, x0, y0, x1, y1);
}

void
p_rgb_cell(p_win *w, unsigned char *rgbs, int ncols, int nrows,
           int x0, int y0, int x1, int y1)
{
  x_image(w, rgbs, 1, ncols, nrows, x0, y0, x1, y1);
}

static void
x_image(p_win *w, unsigned char *bytes, int rgb, int ncols,
        int nrows, int x0, int y0, int x1, int y1)
{
  unsigned char *cell;
  int bpp, xa, ya, xb, yb, r, dr, dy, c, dc, dx;
  int width = x1-x0;
  int height = y1-y0;
  int yexpand = (height>=nrows);
  int xexpand = (width>=ncols);
  int drow = rgb? 3 : 1;
  int dcol = rgb? 3*ncols : ncols;
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  GC gc = x_getgc(s, w, FillSolid);
  XImage *im;

  if (s->image) x_imzap(s);
  s->own_image_data = 1;
  im = s->image = XCreateImage(dpy, DefaultVisual(dpy, s->scr_num),
                               s->depth, ZPixmap, 0, (char *)0,
                               width, height, BitmapPad(dpy), 0);
  bpp = im->bits_per_pixel;

  if (yexpand) {
    r = ((unsigned int)nrows)>>1;
    dr = height%nrows - nrows;
    dy = height/nrows;
  } else {
    r = ((unsigned int)height)>>1;
    dr = nrows%height - height;
    dy = (nrows/height)*dcol;
  }
  if (xexpand) {
    c = ((unsigned int)ncols)>>1;
    dc = width%ncols - ncols;
    dx = width/ncols;
  } else {
    c = ((unsigned int)width)>>1;
    dc = ncols%width - width;
    dx = (ncols/width)*drow;
  }

  if ((long)bpp*(long)width*(long)height > 288L*(long)ncols*(long)nrows) {
    /* less X protocol traffic in a series of filled rectangles */
    int c0 = ((unsigned int)c)>>1;
    x_imzap(s);  /* wasteful to have created image,
                  * but easiest way to get bpp... */
    for (ya=yb=y0 ; ya<y1 ; ya=yb) {
      r += dr;
      if (yexpand) {
        yb += dy;
        if (r>=0) yb++;
        else r += nrows;
      } else {
        yb++;
      }
      cell = bytes;
      c = c0;
      for (xa=xb=x0 ; xa<x1 ; xa=xb) {
        c += dc;
        if (xexpand) {
          xb += dx;
          if (c>=0) xb++;
          else c += ncols;
        } else {
          xb++;
        }
        p_color(w, rgb? P_RGB(cell[0], cell[1], cell[2]) : cell[0]);
        XFillRectangle(dpy, w->d, gc, xa, ya, xb-xa, yb-ya);
        if (xexpand) {
          cell += drow;
        } else {
          cell += dx;
          if (c>=0) cell += drow;
          else c += width;
        }
      }
      if (yexpand) {
        bytes += dcol;
      } else {
        bytes += dy;
        if (r>=0) bytes += dcol;
        else r += height;
      }
    }

  } else {
    /* less X protocol traffic in one image */
    int (*put_pixel)(XImage *, int, int, p_col_t)= 0;
    int bpl = im->bytes_per_line;
    p_col_t pixel;

    im->data = p_malloc(bpl*height);
    if (!im->data) { x_imzap(s); return; }

    if (bpp==8) {
      put_pixel = x_put1, bpp = 1;
    } else if (im->byte_order==MSBFirst) {
      if      (bpp==32) put_pixel = x_put4b, bpp = 4;
      else if (bpp==16) put_pixel = x_put2b, bpp = 2;
      else if (bpp==24) put_pixel = x_put3b, bpp = 3;
    } else {
      if      (bpp==32) put_pixel = x_put4l, bpp = 4;
      else if (bpp==16) put_pixel = x_put2l, bpp = 2;
      else if (bpp==24) put_pixel = x_put3l, bpp = 3;
    }

    if (put_pixel) {
      y1 = height*bpl;
      if (yexpand) dy *= bpl;
      x1 = width*bpp;
      if (xexpand) dx *= bpp;
    } else {
      /* generic case is slow, but still more efficient than setting
       * the image parameters and letting XPutImage do the work --
       * XPutImage uses XGetPixel as well as XPutPixel for this case
       * (this is same as XPutPixel macro in Xutil.h) */
      put_pixel = (int (*)(XImage *, int, int, p_col_t))im->f.put_pixel;
      if (bpp==1) p_color(w, P_FG);
      y1 = height;
      bpl = 1;
      x1 = width;
      bpp = 1;
    }

    if (ncols==width && nrows==height) {
      for (r=0 ; r<y1 ; r+=bpl)
        for (c=0 ; c<x1 ; c+=bpp) {
          pixel = rgb? P_RGB(bytes[0], bytes[1], bytes[2]) : bytes[0];
          bytes += drow;
          put_pixel(im, c, r, x_getpixel(w, pixel));
        }

    } else {
      int i, j;
      int c0 = ((unsigned int)c)>>1;
      for (ya=yb=0 ; yb<y1 ; ya=yb) {
        r += dr;
        if (yexpand) {
          yb += dy;
          if (r>=0) yb += bpl;
          else r += nrows;
        } else {
          yb += bpl;
        }
        cell = bytes;
        c = c0;
        for (xa=xb=0 ; xb<x1 ; xa=xb) {
          c += dc;
          if (xexpand) {
            xb += dx;
            if (c>=0) xb += bpp;
            else c += ncols;
          } else {
            xb += bpp;
          }
          pixel = rgb? P_RGB(cell[0], cell[1], cell[2]) : cell[0];
          pixel = x_getpixel(w, pixel);
          for (j=ya ; j<yb ; j+=bpl)
            for (i=xa ; i<xb ; i+=bpp)
              put_pixel(im, i, j, pixel);
          if (xexpand) {
            cell += drow;
          } else {
            cell += dx;
            if (c>=0) cell += drow;
            else c += width;
          }
        }
        if (yexpand) {
          bytes += dcol;
        } else {
          bytes += dy;
          if (r>=0) bytes += dcol;
          else r += height;
        }
      }
    }

    XPutImage(dpy, w->d, gc, im, 0,0, x0, y0, width, height);
    x_imzap(s);
  }
  if (p_signalling) p_abort();
}

static int
x_put1(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i] = pixel;
  return 0;
}

static int
x_put2b(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i+0] = pixel>>8;
  data[j+i+1] = pixel;
  return 0;
}

static int
x_put2l(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i+0] = pixel;
  data[j+i+1] = pixel>>8;
  return 0;
}

static int
x_put3b(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i+0] = pixel>>16;
  data[j+i+1] = pixel>>8;
  data[j+i+2] = pixel;
  return 0;
}

static int
x_put3l(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i+0] = pixel;
  data[j+i+1] = pixel>>8;
  data[j+i+2] = pixel>>16;
  return 0;
}

static int
x_put4b(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i+0] = pixel>>24;
  data[j+i+1] = pixel>>16;
  data[j+i+2] = pixel>>8;
  data[j+i+3] = pixel;
  return 0;
}

static int
x_put4l(XImage *im, int i, int j, p_col_t pixel)
{
  unsigned char *data = (unsigned char *)im->data;
  data[j+i+0] = pixel;
  data[j+i+1] = pixel>>8;
  data[j+i+2] = pixel>>16;
  data[j+i+3] = pixel>>24;
  return 0;
}
