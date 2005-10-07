/*
 * color.m
 * Color handling for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

void
p_color(p_win *w, unsigned long color)
{ int change = 0;
  p_scr* s = w->s;
  View* view = w->view;
  CGContextRef cr = w->cr;
  if (w->parent) w = w->parent;
  if (p_signalling)
  { p_abort();
    return;
  }
  if (w->color != color) change = 1;
  /* p_color may be called before p_rect, without p_clip getting
   * called first. In that case, the focus is not yet locked on the
   * view. */
  if (s->lockedView != view)
  { CGAffineTransform tf = CGContextGetTextMatrix(cr);
    tf.a = +1.0;
    tf.d = -1.0; /* -1.0 so that text appears upright in a flipped view */
    CGContextSetTextMatrix(cr, tf);
    if (s->lockedView) [s->lockedView unlockFocus];
    [view lockFocus];
    s->lockedView = view;
    change = 1;
  }
  if (change) {
/* It seems that Cocoa doesn't support XOR drawing at this point.
    const int de_xor = (w->color == P_XOR);
    if (de_xor || color==P_XOR) {
      if (w->view)
        SetROP2(w->dc, de_xor? R2_COPYPEN : R2_NOT);
    }
*/
    w->color = color;
    if (!P_IS_RGB(color)) {
      w->components[0] = w->pixels[color][0];
      w->components[1] = w->pixels[color][1];
      w->components[2] = w->pixels[color][2];
      w->components[3] = w->pixels[color][3];
    } else {
      w->components[0] = P_R(color) / 255.0;
      w->components[1] = P_G(color) / 255.0;
      w->components[2] = P_B(color) / 255.0;
      w->components[3] = 1.0;
    }
    CGContextSetStrokeColorSpace(cr, s->colorspace);
    CGContextSetStrokeColor(cr, w->components);
    CGContextSetFillColorSpace(cr, s->colorspace);
    CGContextSetFillColor(cr, w->components);
  }
}
