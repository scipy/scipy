/*
 * clips.m
 * p_clip for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

void
p_clip(p_win *w, int x0, int y0, int x1, int y1)
{
  p_scr* s = w->s;
  View* view = w->view;
  CGContextRef cr = w->cr;
  static CGRect rect;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (s->lockedView) [s->lockedView unlockFocus];
  [view lockFocus];
  s->lockedView = view;
  if (x1 > x0 || y1 > y0)
  { rect.origin.x = x0;
    rect.origin.y = y0;
    rect.size.width = x1-x0;
    rect.size.height = y1-y0;
    CGContextClipToRect(cr, rect);
  }
  /* Note that lockFocus resets the colors in the graphics state */
  CGContextSetStrokeColorSpace(cr, s->colorspace);
  CGContextSetStrokeColor(cr, w->components);
  CGContextSetFillColorSpace(cr, s->colorspace);
  CGContextSetFillColor(cr, w->components);
  /* tf.a and tf.d are reset to 13.0 after the window moves to the front
   * (no clue why). Reset them here. */
  CGAffineTransform tf = CGContextGetTextMatrix(cr);
  tf.a = +1.0;
  tf.d = -1.0; /* -1.0 so that text appears upright in a flipped view */
  CGContextSetTextMatrix(cr, tf);
}
