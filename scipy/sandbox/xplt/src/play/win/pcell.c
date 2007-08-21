/*
 * pcell.c -- $Id$
 * p_ndx_cell, p_rgb_cell for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"

static void w_cell(p_win *w, unsigned char *ndxs, unsigned char *rgbs,
                   int ncols, int nrows, int x0, int y0, int x1, int y1);

static void *wbm_alloc(HDC dc, BITMAPINFO *bmi, HDC *pbmdc, HBITMAP *pbm);
static void wbm_free(HDC bmdc, HBITMAP bm);
static BITMAPINFO *wbm_head(BITMAPINFO *bmi, int ncols, int nrows, int nbits);

void
p_ndx_cell(p_win *w, unsigned char *ndxs, int ncols, int nrows,
           int x0, int y0, int x1, int y1)
{
  w_cell(w, ndxs, 0, ncols, nrows, x0, y0, x1, y1);
}

void
p_rgb_cell(p_win *w, unsigned char *rgbs, int ncols, int nrows,
           int x0, int y0, int x1, int y1)
{
  if (w->s->sys_pal && !w->rgb_mode) {
    p_palette(w, p_595, 225);
    w->rgb_mode = 1;
  }
  w_cell(w, 0, rgbs, ncols, nrows, x0, y0, x1, y1);
}

void
p_bitblt(p_win *w, int x, int y, p_win *offscreen,
         int x0, int y0, int x1, int y1)
{
  HDC dc = w_getdc(w, 0);
  HDC dc2 = w_getdc(offscreen, 0);
  if (dc && dc2)
    BitBlt(dc, x, y, x1-x0, y1-y0, dc2, x0, y0, SRCCOPY);
}

void
p_rgb_read(p_win *w, unsigned char *rgbs,
           int x0, int y0, int x1, int y1)
{
  HDC dc = w_getdc(w, 0);
  x1 -= x0;
  y1 -= y0;
  if (dc) {
    BITMAPINFO bmi;
    HDC bmdc;
    HBITMAP bm;
    DWORD *bmrgb = wbm_alloc(dc, wbm_head(&bmi, x1, y1, 32), &bmdc, &bm);
    if (bmrgb) {
      int i;
      long j, npix = (long)x1 * (long)y1;
      BitBlt(bmdc, 0, 0, x1, y1, dc, x0, y0, SRCCOPY);
      for (j=0 ; j<npix ; j+=x1) {
        for (i=0 ; i<x1 ; i++,rgbs+=3) {
          rgbs[0] = GetBValue(bmrgb[i+j]);                /* yes, backwards */
          rgbs[1] = GetGValue(bmrgb[i+j]);
          rgbs[2] = GetRValue(bmrgb[i+j]);
        }
      }
      wbm_free(bmdc, bm);
    }
  }
}

static void
w_cell(p_win *w, unsigned char *ndxs, unsigned char *rgbs,
       int ncols, int nrows, int x0, int y0, int x1, int y1)
{
  /* huge effort here failed to make 8-bit DIBs work when the
   * system is set to 256-color mode (w->s->sys_pal!=0)
   * -- expect this to perform poorly in that mode */
  HDC dc = w_getdc(w, 0);
  if (dc) {
    struct {
      BITMAPINFOHEADER bmih;
      WORD ndx[256];
    } h;
    BITMAPINFO *bmi = (BITMAPINFO *)&h;
    int nbits = (w->s->sys_pal && w->palette)? 8 : 32;
    HDC bmdc;
    HBITMAP bm;
    DWORD *bmrgb = wbm_alloc(dc, wbm_head(bmi, ncols, nrows, nbits),
                            &bmdc, &bm);
    if (bmrgb) {
      long i, j;
      if (nbits == 8) {
        unsigned int offset = w->s->sys_offset;
        COLORREF cr, *sys_index = w->s->sys_index;
        int n = w->rgb_mode? 225 :
          (w->parent? w->parent->n_pixels : w->n_pixels);
        for (i=0 ; i<n ; i++) h.ndx[i] = (WORD)(offset + i);
        for (; i<=P_EXTRA ; i++) h.ndx[i] = 0;  /* black */
        for (; i<256 ; i++) {
          cr = w->s->sys_index[255-i] & 0xff;
          if (cr < offset) h.ndx[i] = (WORD)cr;
          else  h.ndx[i] = (WORD)(n+(cr & 0xff));
        }
      }
      if (ndxs) {
        if (nbits != 8) {
          long npix = (long)ncols * (long)nrows;
          COLORREF pxl, *pixels = w->parent? w->parent->pixels : w->pixels;
          for (j=0 ; j<npix ; j+=ncols) {
            for (i=0 ; i<ncols ; i++,ndxs++) {
              pxl = pixels[ndxs[0]];
              bmrgb[i+j] = W_SWAPRGB(pxl);                /* yes, backwards */
            }
          }
        } else {
          int nc = (ncols+3) & ~3;
          long npix = (long)nc * (long)nrows;
          unsigned char *bmbyt = (unsigned char *)bmrgb;
          for (j=0 ; j<npix ; j+=nc) {
            for (i=0 ; i<ncols ; i++,ndxs++)
              bmbyt[i+j] = ndxs[0];
          }
        }
      } else {
        if (nbits != 8) {
          long npix = (long)ncols * (long)nrows;
          for (j=0 ; j<npix ; j+=ncols) {
            for (i=0 ; i<ncols ; i++,rgbs+=3)
              bmrgb[i+j] = RGB(rgbs[2],rgbs[1],rgbs[0]);  /* yes, backwards */
          }
        } else {
          int nc = (ncols+3) & ~3;
          long npix = (long)nc * (long)nrows;
          unsigned char *bmbyt = (unsigned char *)bmrgb;
          unsigned int r, g, b;
          for (j=0 ; j<npix ; j+=nc) {
            for (i=0 ; i<ncols ; i++,rgbs+=3) {
              r = rgbs[0], g = rgbs[1], b = rgbs[2];
              r = (r+32)>>6;
              g = (g+16)>>5;
              b = (b+32)>>6;
              g += b+(b<<3);
              bmbyt[i+j] = (unsigned char)(r+g+(g<<2));
            }
          }
        }
      }
      StretchDIBits(dc, x0, y0, x1-x0, y1-y0, 0, 0, ncols, nrows,
                    bmrgb, bmi, nbits!=8? DIB_RGB_COLORS : DIB_PAL_COLORS,
                    SRCCOPY);
      wbm_free(bmdc, bm);
    }
  }
}

static void *
wbm_alloc(HDC dc, BITMAPINFO *bmi, HDC *pbmdc, HBITMAP *pbm)
{
  int nbits = bmi->bmiHeader.biBitCount;
  void *bmrgb = 0;
  HBITMAP bm = CreateDIBSection(dc, bmi,
                               nbits!=8? DIB_RGB_COLORS : DIB_RGB_COLORS,
                               &bmrgb, 0, 0);
  if (bm) {
    HDC bmdc = CreateCompatibleDC(dc);
    if (bmdc) {
      SelectObject(bmdc, bm);
      *pbmdc = bmdc;
      *pbm = bm;
    } else {
      DeleteObject(bm);
    }
  }
  return bmrgb;
}

static void
wbm_free(HDC bmdc, HBITMAP bm)
{
  if (bmdc) DeleteDC(bmdc);
  if (bm) DeleteObject(bm);
}

static BITMAPINFO *
wbm_head(BITMAPINFO *bmi, int ncols, int nrows, int nbits)
{
  bmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
  bmi->bmiHeader.biWidth = ncols;
  bmi->bmiHeader.biHeight = -nrows;  /* top-down order */
  bmi->bmiHeader.biPlanes = 1;
  bmi->bmiHeader.biBitCount = nbits;
  bmi->bmiHeader.biCompression = BI_RGB;
  bmi->bmiHeader.biSizeImage = 0;
  bmi->bmiHeader.biXPelsPerMeter = 0;
  bmi->bmiHeader.biYPelsPerMeter = 0;
  bmi->bmiHeader.biClrUsed = 0;
  bmi->bmiHeader.biClrImportant = 0;
  return bmi;
}
