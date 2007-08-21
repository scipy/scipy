/*
 * text.m
 * p_text for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

void
p_text(p_win *w, int x0, int y0, const char *text, int n)
{
  View* view = w->view;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (view) {
    int i;
    CGContextRef cr = w->cr;
    if (n <= 0) n = 16350;
    for (i=0 ; i<n ; i++) if (!text[i]) break;
    n = i;
    if (!w->orient) {
      CGContextSetTextPosition (cr, (float)x0, (float)y0);
      for (i = 0; i < n; i++)
      { CGGlyph c = (unsigned short)(text[i]) - 29;
        CGContextShowGlyphs(cr, &c, 1);
      }
      // The following would have been easier, but unfortunately it only works
      // if the font is selected using CGContextSelectFont. That, however, would
      // mean that we need to create the font each time we want to draw some
      // text. Hence, use CGContextSetFont with CGContextShowGlyphs instead.
      // CGContextShowTextAtPoint(cr, (float)x0, (float)y0, text, n);
    } else {
      const float angle = -M_PI/2;
      CGContextRotateCTM(cr, angle);
      CGContextSetTextPosition (cr, -(float)y0, (float)x0);
      for (i = 0; i < n; i++)
      { CGGlyph c = (unsigned short)(text[i]) - 29;
        CGContextShowGlyphs(cr, &c, 1);
      }
      // Same comment as above.
      // CGContextShowTextAtPoint(cr, -(float)y0, (float)x0, text, n);
      CGContextRotateCTM(cr, -angle);
    }
  }
}
