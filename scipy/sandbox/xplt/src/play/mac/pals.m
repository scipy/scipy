/*
 * pals.m
 * p_palette for Mac OS X.
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playm.h"
#include "pstdlib.h"

void
p_palette(p_win *w, unsigned long *colors, int n)
{
  p_scr *s = w->s;
  int i;

  if (n > 256) n = 256;
  if (n > 240) n = 240;    /* system reserves 20 colors -- may only get 236 */

  if (w->parent) w = w->parent;

  if (!p_signalling) {
    for (i=0 ; i<n ; i++) {
      w->pixels[i][0] = P_R(colors[i]) / 255.0;
      w->pixels[i][1] = P_G(colors[i]) / 255.0;
      w->pixels[i][2] = P_B(colors[i]) / 255.0;
      w->pixels[i][3] = 1.0;
    }
    for (; i<w->n_pixels ; i++)
    { w->pixels[i][0] = s->sys_colors[255-P_FG][0];
      w->pixels[i][1] = s->sys_colors[255-P_FG][1];
      w->pixels[i][2] = s->sys_colors[255-P_FG][2];
      w->pixels[i][3] = s->sys_colors[255-P_FG][3];
    }
    w->n_pixels = n;
  }
}
