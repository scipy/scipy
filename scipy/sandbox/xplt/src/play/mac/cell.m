/*
 * cell.m
 * p_ndx_cell, p_rgb_cell for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

static void m_cell(p_win *w, unsigned char *ndxs, unsigned char *rgbs,
                   int ncols, int nrows, int x0, int y0, int x1, int y1);
static void m_release_data(void *info, const void *data, size_t size);


void
p_ndx_cell(p_win *w, unsigned char *ndxs, int ncols, int nrows,
           int x0, int y0, int x1, int y1)
{
  m_cell(w, ndxs, 0, ncols, nrows, x0, y0, x1, y1);
}

void
p_rgb_cell(p_win *w, unsigned char *rgbs, int ncols, int nrows,
           int x0, int y0, int x1, int y1)
{
  m_cell(w, 0, rgbs, ncols, nrows, x0, y0, x1, y1);
}

void
p_rgb_read(p_win *w, unsigned char *rgbs,
           int x0, int y0, int x1, int y1)
{
  printf ("p_rgb_read\n");
  /* This function is needed for g_rgb_read, which is never called. */
}

static void
m_cell(p_win *w, unsigned char *ndxs, unsigned char *rgbs,
       int ncols, int nrows, int x0, int y0, int x1, int y1)
{
  View* view = w->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view) {
    CGContextRef cr = w->cr;
    const size_t nComponents = 3; /* red, green, blue */
    const size_t bytesPerComponent = 1;
    const size_t bitsPerComponent = 8 * bytesPerComponent;
    const size_t bitsPerPixel = bitsPerComponent * nComponents;
    const size_t bytesPerRow = nComponents * bytesPerComponent * ncols;
    const size_t size = bytesPerRow * nrows;
    unsigned char* data = (unsigned char*)malloc(size*sizeof(unsigned char));
    if (data)
    { int ii, jj;
      CGColorSpaceRef colorspace = CGColorSpaceCreateDeviceRGB();
      if (ndxs) {
        int kk = nrows*ncols;
        unsigned char* index = ndxs;
        for (ii = 0; ii < nrows; ii++) {
          kk -= ncols;
          for (jj = 0; jj < ncols; jj++, index++) {
            float* color = w->pixels[*index];
            data[3*(kk+jj)]   = 256 * color[0];
            data[3*(kk+jj)+1] = 256 * color[1];
            data[3*(kk+jj)+2] = 256 * color[2];
          }
        }
      } else {
        int kk = size;
        for (ii = 0; ii < size; ii+=bytesPerRow) {
          kk -= bytesPerRow;
          for (jj = 0; jj < bytesPerRow; jj++) {
            data[kk+jj] = rgbs[ii+jj]; 
          }
        }
      }
      CGDataProviderRef provider = CGDataProviderCreateWithData (NULL,
                                                                 data,
                                                                 size, 
                                                                 m_release_data); 
      CGImageRef bitmap = CGImageCreate (ncols,
                                         nrows,
                                         bitsPerComponent,
                                         bitsPerPixel,
                                         bytesPerRow,
                                         colorspace,
                                         kCGImageAlphaNone,
                                         provider,
                                         NULL,
                                         false,
                                         kCGRenderingIntentDefault);
      CGColorSpaceRelease(colorspace);

      if(bitmap)
      { CGContextDrawImage(cr, CGRectMake(x0,y0,x1-x0,y1-y0), bitmap);
        CGImageRelease(bitmap);
      }
      CGDataProviderRelease(provider);
    }
  }
}

static void
m_release_data(void *info, const void *data, size_t size)
{ unsigned char* p = (unsigned char*) data;
  free(p);
}
