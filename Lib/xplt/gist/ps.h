/*
 * PS.H
 *
 * $Id$
 *
 * Declare the PostScript engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef PS_H
#define PS_H

#include "gist.h"
#include "engine.h"

#include <stdio.h>

typedef struct GpsBBox GpsBBox;
struct GpsBBox {
  int xll, yll, xur, yur;
};

typedef struct PSEngine PSEngine;
struct PSEngine {
  Engine e;

  /* --------------- Specific to PSEngine ------------------- */

  char *filename;
  FILE *file;     /* 0 until file is actually written into */
  int closed;     /* if file==0 and closed!=0, there was a write error */

  /* Page orientation and color table can only be changed at the beginning
     of each page, so the entries in the Engine base class are repeated
     here as the "currently in effect" values.  When a new page begins,
     these values are brought into agreement with those in the base
     class.  The ChangePalette virtual function temporarily resets colorMode
     to 0, to avoid any references to the new palette until the
     next page begins.  */
  int landscape;
  int colorMode;
  int nColors;

  GpsBBox pageBB;   /* bounding box for current page */
  GpsBBox docBB;    /* bounding box for entire document */
  int currentPage;  /* current page number, incremented by EndPage */
  long fonts;       /* bits correspond to Gist fonts (set when used) */

  /* The ps.ps GistPrimitives assume knowledge of the several
     graphical state parameters, mostly from gistA.  These are
     reset at the beginning of each page.  */
  GpBox clipBox;
  int curClip;
  int curColor;
  int curType;
  GpReal curWidth;
  int curFont;
  GpReal curHeight;
  int curAlignH, curAlignV;
  int curOpaque;

  /* When clipping is turned off, the ps.ps GistPrimitives state
     partially reverts to its condition when clipping was turned on.  */
  int clipColor;
  int clipType;
  GpReal clipWidth;
  int clipFont;
  GpReal clipHeight;

  char line[80];   /* buffer in which to build current output line */
  int nchars;      /* current number of characters in line */
};

extern PSEngine *GisPSEngine(Engine *engine);

#endif
