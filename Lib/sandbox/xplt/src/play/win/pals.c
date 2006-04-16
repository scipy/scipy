/*
 * pals.c -- $Id$
 * p_palette for MS Windows
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

void
p_palette(p_win *w, unsigned long *colors, int n)
{
  p_scr *s = w->s;
  int i;
  unsigned int r, g, b;
  PALETTEENTRY space[257];
  LOGPALETTE *pal = (LOGPALETTE *)space;
  PALETTEENTRY *entry = pal->palPalEntry;
  PALETTEENTRY *sys_pal = s->sys_pal;
  int offset = s->sys_offset;
  int change_pal = (sys_pal != 0);
  int shallow = change_pal;

  if (n+2*offset > 256) n = 256-2*offset;
  if (n > 240) n = 240;    /* system reserves 20 colors -- may only get 236 */

  /* but look out for metafile... it should always copy parent
   * window palette at creation time, never call p_palette */
  if (w->parent) w = w->parent;
  if (w->menu) shallow = change_pal = 0;
  else if (w->rgb_mode && change_pal) shallow = 2, change_pal = 0;

  if (!p_signalling) {
    for (i=0 ; i<n ; i++) {
      r = P_R(colors[i]);
      g = P_G(colors[i]);
      b = P_B(colors[i]);
      if (!shallow) {
        w->pixels[i] = RGB(r, g, b);
      } else if (shallow==1) {
        w->pixels[i] = PALETTEINDEX(offset+i);
      } else {  /* must be consistent with w_color in getdc.c */
        r = (r+32)>>6;
        g = (g+16)>>5;
        b = (b+32)>>6;
        g += b+(b<<3);
        w->pixels[i] = PALETTEINDEX(offset+r+g+(g<<2));  /* r+5*g+45*b */
      }
      if (!change_pal) continue;
      entry[offset+i].peRed = (unsigned char)r;
      entry[offset+i].peGreen = (unsigned char)g;
      entry[offset+i].peBlue = (unsigned char)b;
      entry[offset+i].peFlags = 0;
    }
    for (; i<w->n_pixels ; i++)
      w->pixels[i] = s->sys_colors[255-P_FG];
    w->n_pixels = n;

    if (change_pal) {
      pal->palVersion = 0x300;     /* apparently this means Windows 3.0 */
      pal->palNumEntries = 2*offset+n;
      if (w->palette) {
        UnrealizeObject(w->palette);
        DeleteObject(w->palette);
        w->palette = 0;
      }
      if (n) {
        for (i=0 ; i<offset ; i++) entry[i] = s->sys_pal[i];
        for (; i<2*offset ; i++) entry[n+i] = s->sys_pal[i];
        w->palette = CreatePalette(pal);
        if (w->palette) {
          if (s->active==w || !s->active) {
            if (!s->active) s->active = w;
            SelectPalette(w->dc, w->palette, 0);
            /* no InvalidateRect here -- arrange to only call
             * p_palette before a redraw display list sequence */
          }
          SelectPalette(w->dc, w->palette, 1);
          RealizePalette(w->dc);
        }
      }
      if (!w->palette && s->active==w) s->active = 0;
    }
  }
}
