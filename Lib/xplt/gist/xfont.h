/*
 * XFONT.H
 *
 * $Id$
 *
 * Declare X windows font utilities for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef XFONT_H
#define XFONT_H

#include "gist.h"
#include "xicky.h"

/* Get the closest thing to (font,height) known to the given display
   at the given resolution.  May be wildly off if fonts don't exist,
   but guaranteed to return SOMETHING.  Up to five fonts are maintained
   without having to query the server.  Return loadedID.  */
extern int GxIDFont(GxDisplay *xdpy, int dpi, int font, GpReal height);
/* Build the name of the font having the given loadedID-- should be
   returned by GxIDFont to be sure it actually exists.  */
extern char *GxNameFont(int id);

/* ifont==font>>2 or -1, slant, bold==0, 1, or -1, isize=0..5 or -1,
   where -1 means not specified.  The nsize mask bits are set if
   the given ifont, slant, and bold are available in that size.
   Returned char** is static in all cases.  Don't modify or free it.  */
extern char **GxFontFaces(GxFontProps *fontProps, int isize,
			  int slant, int bold, int *nfonts, int *mask);
extern char **GxFontSizes(GxFontProps *fontProps, int ifont,
			  int slant, int bold, int *nsizes, int *mask);
extern char **GxFontSlants(GxFontProps *fontProps, int ifont,
			   int isize, int bold, int *nslants, int *mask);
extern char **GxFontWeights(GxFontProps *fontProps, int ifont,
			    int isize, int slant, int *nweights, int *mask);
extern int gxFontSize[FONT_SIZES];  /* 8, 10, 12, 14, 18, 24 */

extern void GxGrabFonts(GxDisplay *xdpy, char *permFont);

#endif
