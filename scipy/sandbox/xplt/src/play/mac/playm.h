/*
 * playm.h
 * Mac OS X-private portability layer declarations
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "play.h"
#include <Cocoa/Cocoa.h>

/* ------------------------------------------------------------------------ */

@interface View: NSView
{ int mousebutton; // Mouse button state
  p_win* pw;
}
- (void)expose:(id)dummy;
@end

/* ------------------------------------------------------------------------ */


struct p_scr {
  int width, height, depth;
  int x0, y0;                  /* usable area may not start at (0,0) */

  float sys_colors[15][4];
  CGColorSpaceRef colorspace;

  ATSFontMetrics metric; /* Standard gui font */
  CGFontRef font;        /* Standard gui font */
  ATSFontMetrics metrics[5][4]; /* Font cache */
  CGFontRef fonts[5][4];        /* Font cache */
  CGContextRef scratch;
  /* scratch space for drawing strings to determine their width */

  View* lockedView; /* Currently locked view */

  NSCursor* cursors[P_NONE]; /* Will be loaded as needed */
};

struct p_win {
  void *ctx;
  p_scr *s;

  NSWindow* w;
  View* view;
  CGContextRef cr;

  p_win *parent;
  /* DC keeps three separate colors, only set when actually used */
  float pixels[256][4];
  int n_pixels;
  float components[4];
  /* (r,g,b,alpha) components of the color to be used */
  unsigned long color;
  /* identifier of the color currently loaded in components */

  int width; /* line width */
  int dash;  /* dashed line style */
  int squarelinestyle; /* Square/round caps and miter/round joining */

  CGFontRef fontref;      /* Currently loaded font */
  int fontid;
  int orient;      /* font orientation */

  unsigned long bg;          /* for p_clear */
};
