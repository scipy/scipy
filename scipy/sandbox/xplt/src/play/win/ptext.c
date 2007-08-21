/*
 * ptext.c -- $Id$
 * p_text for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

void
p_text(p_win *w, int x0, int y0, const char *text, int n)
{
  HDC dc = w_getdc(w, 1);
  if (dc) {
    int i;
    if (n <= 0) n = 16350;
    for (i=0 ; i<n ; i++) if (!text[i]) break;
    n = i;
    if (!w->orient) {
      TextOut(dc, x0, y0, text, n);
    } else {
      /* only do this with opaque background,
       * since only 8x8 patterned brushes supported under Win95
       * (could read screen, rotate, TextOut, rotate back) */
      TEXTMETRIC tm;
      SIZE sz;
      if (n>0 && GetTextMetrics(dc, &tm) &&
          GetTextExtentPoint32(dc, text, n, &sz) &&
          sz.cx>0) {
        int width = sz.cx;
        int height = tm.tmHeight + 1;
        int base = tm.tmAscent;
        int scandx = ((width-1)/(8*sizeof(DWORD)) + 1)*(8*sizeof(DWORD));
        int scandy = ((height-1)/(8*sizeof(DWORD)) + 1)*(8*sizeof(DWORD));
        HDC offdc = CreateCompatibleDC(dc);
        HBITMAP offbm = CreateCompatibleBitmap(dc, scandx, scandy);
        HFONT font = GetCurrentObject(dc, OBJ_FONT);
        void *monodat = p_malloc(scandx*scandy/4);
        void *monorot = (char *)monodat + (scandx*scandy/8);
        struct {
          BITMAPINFOHEADER bmiHeader;
          DWORD bmiColors[2];
        } mono;
        mono.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        mono.bmiHeader.biWidth = scandx;
        mono.bmiHeader.biHeight = -scandy;
        mono.bmiHeader.biPlanes = 1;
        mono.bmiHeader.biBitCount = 1;
        mono.bmiHeader.biCompression = BI_RGB;
        mono.bmiHeader.biSizeImage = 0;
        mono.bmiHeader.biXPelsPerMeter = 0;
        mono.bmiHeader.biYPelsPerMeter = 0;
        mono.bmiHeader.biClrUsed = 0;
        mono.bmiHeader.biClrImportant = 0;
        /* apparently, bmiColors ignored on read
         * -white bits are set, black bits reset no matter what */
        mono.bmiColors[0] = RGB(0,0,0);
        mono.bmiColors[1] = RGB(255,255,255);
        if (offdc && offbm && font && monodat) {
          RECT r;
          SelectObject(offdc, offbm);
          SelectObject(offdc, font);
          SetBkMode(offdc, TRANSPARENT);
          SetTextColor(offdc, RGB(255,255,255));
          SetTextAlign(offdc, TA_LEFT | TA_BASELINE | TA_NOUPDATECP);
          r.left = r.top = 0;
          r.right = scandx;
          r.bottom = scandy;
          FillRect(offdc, &r, GetStockObject(BLACK_BRUSH));
          if (w->orient == 2) {
            TextOut(offdc, scandx-width, base, text, n);
          } else if (w->orient == 1) {
            TextOut(offdc, scandx-width, base, text, n);
          } else {
            TextOut(offdc, 0, scandy-height+base, text, n);
          }
          GetDIBits(offdc, offbm, 0, scandy, monodat,
                    (BITMAPINFO *)&mono, DIB_RGB_COLORS);
        }
        if (offdc) DeleteDC(offdc);
        if (offbm) DeleteObject(offbm);

        if (monodat) {
          COLORREF bg = w_color(w,w->bg);
          COLORREF fg = w_color(w,w->color);
          if (w->orient == 2) {
            p_mrot180(monodat, monorot, scandx, height);
            x0 -= width-1;
            y0 -= (height-base)-1;
            scandx = width, width = height, height = scandx;
          } else if (w->orient == 1) {
            p_mrot090(monodat, monorot, scandx, scandy);
            x0 -= base;
            y0 -= width-1;
          } else {
            p_mrot270(monodat, monorot, scandx, scandy);
            x0 -= (height-base)-1;
          }
          mono.bmiColors[0] = W_SWAPRGB(bg);
          mono.bmiColors[1] = W_SWAPRGB(fg);
          mono.bmiHeader.biWidth = height;
          mono.bmiHeader.biHeight = -width;
          /* prefer this to SetDIBitsToDevice? */
          StretchDIBits(dc, x0, y0,
                        height, width,
                        0, 0, height, width, monorot,
                        (BITMAPINFO *)&mono, DIB_RGB_COLORS, SRCCOPY);

          p_free(monodat);
        }
      }
    }
  }
}
