/*
 * CGM.H
 *
 * $Id$
 *
 * Declare the CGM binary metafile engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef CGM_H
#define CGM_H

#include "gist.h"
#include "engine.h"

#include <stdio.h>
#ifndef SEEK_CUR
/* Sun strikes again...  This may break other non-standard machines */
#define SEEK_CUR 1
#endif

/* This Engine is based on the ANSI CGM Standard, ANSI X3.122 - 1986
   Parts 1 and 3.  It is fully conforming with the following caveats:
   1. The font names in the FONT LIST metafile descriptor element are
      standard PostScript font names, rather than ISO standard names
      (which I couldn't find).  The font numbers match the numbering
      scheme used by the GPLOT program.
   2. A Gist smooth polyline is output using the ordinary POLYLINE
      primitive, preceded by an APPLICATION DATA comment.  This gives
      reasonable output (identical to that produced by the Gist X
      engine, in fact), but not as nice as the Gist PostScript engine.
 */

typedef unsigned char Octet;  /* as defined in the CGM standard */
#define MAX_PARTITION 0x7ffc  /* maximum length of a CELL ARRAY partition */

typedef struct CGMEngine CGMEngine;
struct CGMEngine {
  Engine e;

  /* --------------- Specific to CGMEngine ------------------- */

  char *filename;
  GpReal scale;   /* VDC units per Gist NDC unit, 25545.2 by default */
  long fileSize;  /* approximate maximum in bytes, default is 1 Meg */
  void (*IncrementName)(char *filename);
                  /* function to increment filename IN PLACE
                     (a reasonable default is supplied) */
  p_file *file;   /* 0 until file is actually written into */
  int state;      /* CGM state as described in fig. 12 of ANSI X3.122 */

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

  int currentPage;  /* current page number, incremented by EndPage */

  /* The CGM engine must keep track of several graphical state
     parameters, mostly from gistA.  These are reset at the
     beginning of each page.  */
  GpBox clipBox;
  int curClip;
  unsigned long curColor[5];
  int curType;
  GpReal curWidth;
  int curMark;
  GpReal curSize;
  int curFont;
  GpReal curHeight;
  int curAlignH, curAlignV;
  int curPath;
  int curOpaque;  /* unlike Gist, this applies to more than text... */
  int curEtype;
  GpReal curEwidth;
};

extern CGMEngine *GisCGMEngine(Engine *engine);

/* Default CGM scale factor and maximum file size */
extern GpReal gCGMScale;
extern long gCGMFileSize;

/* To change away from the default fileSize, just change the
   CGMEngine fileSize member.  To change the CGM scale factor
   (the initial default value is 25545.2, which corresponds to
   2400 dpi resolution), use GcgmSetScale after creating the
   engine but BEFORE DRAWING ANYTHING WITH IT: */
extern void GcgmSetScale(CGMEngine *cgmEngine, GpReal scale);

#endif
