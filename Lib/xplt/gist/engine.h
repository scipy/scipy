/*
 * ENGINE.H
 *
 * $Id$
 *
 * Declare common properties of all GIST engines
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef ENGINE_H
#define ENGINE_H

#include "gist.h"

/* ------------------------------------------------------------------------ */

struct Engine {
  Engine *next;
  Engine *nextActive;
  char *name;
  char *type;    /* same for every instance of a given Engine type */

  int active;
  int marked;    /* set if any marks have been made on current page */

  int landscape;          /* non-0 if page is wider than tall */
  GpTransform transform;    /* viewport in NDC, window in VDC */
  GpXYMap devMap; /* actual coefficients for NDC->VDC mapping */
  GpXYMap map;     /* actual coefficients for WC->VDC mapping */

  /* Current color palette (see gist.h) */
  int colorChange;  /* set this to alert Engine that palette has changed */
  int colorMode;    /* for hardcopy devices, set to dump palette to file */
  int nColors;
  GpColorCell *palette;  /* pointer is NOT owned by this engine */

  /* Handling incremental changes and damage to a drawing.  */
  Drawing *drawing;
  int lastDrawn;
  long systemsSeen[2];
  int inhibit;     /* set by GdAlert if next primitive is to be ignored */
  int damaged;
  GpBox damage;    /* Never used if no ClearArea function provided */

  /* --------------- Virtual function table ------------------- */

  /* Close device and free all related memory (GpKillEngine)  */
  void (*Kill)(Engine *engine);

  /* New page (frame advance) */
  int (*Clear)(Engine *engine, int always);

  /* Flush output buffers */
  int (*Flush)(Engine *engine);

  /* Change coordinate transformation (GpSetTrans)
     If engine damaged, must set clipping to damage.  */
  void (*ChangeMap)(Engine *engine);

  /* Change color palette (GpSetPalette) */
  int (*ChangePalette)(Engine *engine); /* returns maximum palette size */

  /* Polyline primitive.  Segments have already been clipped.  */
  int (*DrawLines)(Engine *engine, long n, const GpReal *px,
		   const GpReal *py, int closed, int smooth);
  /* Polymarker primitive.  Points have already been clipped.
     Can use GpPseudoMark if no intrinsic polymarker primitive.  */
  int (*DrawMarkers)(Engine *engine, long n,
		     const GpReal *px, const GpReal *py);
  /* Text primitive.  No clipping has yet been done.  */
  int (*DrawText)(Engine *engine, GpReal x0, GpReal y0, const char *text);
  /* Filled area primitive.  Polygon has already been clipped.  */
  int (*DrawFill)(Engine *engine, long n, const GpReal *px, const GpReal *py);
  /* Cell array primitive.  No clipping has yet been done.  */
  int (*DrawCells)(Engine *engine, GpReal px, GpReal py, GpReal qx, GpReal qy,
		   long width, long height, long nColumns,
		   const GpColor *colors);
  /* Disjoint line primitive.  Segments have already been clipped.  */
  int (*DrawDisjoint)(Engine *engine, long n, const GpReal *px,
		      const GpReal *py, const GpReal *qx, const GpReal *qy);

  /* Clear a damaged area (box in NDC units) - set to 0 by GpNewEngine */
  void (*ClearArea)(Engine *engine, GpBox *box);
};

/* Linked lists of all engines and of all active engines */
extern Engine *gistEngines;
extern Engine *gistActive;

/* Generic Engine constructor and destructor */
extern Engine *GpNewEngine(long size, char *name, char *type,
			   GpTransform *transform, int landscape,
  void (*Kill)(Engine*), int (*Clear)(Engine*,int), int (*Flush)(Engine*),
  void (*ChangeMap)(Engine*), int (*ChangePalette)(Engine*),
  int (*DrawLines)(Engine*,long,const GpReal*,const GpReal*,int,int),
  int (*DrawMarkers)(Engine*,long,const GpReal*,const GpReal *),
  int (*DrawText)(Engine*e,GpReal,GpReal,const char*),
  int (*DrawFill)(Engine*,long,const GpReal*,const GpReal*),
  int (*DrawCells)(Engine*,GpReal,GpReal,GpReal,GpReal,
		   long,long,long,const GpColor*),
  int (*DrawDisjoint)(Engine*,long,const GpReal*,const GpReal*,
		      const GpReal*,const GpReal*));
extern void GpDelEngine(Engine *engine);

/* ------------------------------------------------------------------------ */
/* Coordinate mapping */

/* Set engine->devMap from engine->transform */
extern void GpDeviceMap(Engine *engine);

/* Compose WC->VDC engine->map coefficients given gistT (WC->NDC) and
   engine->devMap (NDC->VDC).  */
extern void GpComposeMap(Engine *engine);

/* The X window, CGM, and PostScript devices can all be based on
   integer coordinate systems using points and segments which are
   short integers.  Here are some common routines for forming and
   manipulating these.  */

typedef struct GpPoint GpPoint;
struct GpPoint {  /* same as X windows XPoint */
  short x, y;
};

typedef struct GpSegment GpSegment;
struct GpSegment {  /* same as X windows XSegment */
  short x1, y1, x2, y2;
};

/* Returns number of points processed this pass (<=maxPoints) */
extern long GpIntPoints(const GpXYMap *map, long maxPoints, long n,
			const GpReal *x, const GpReal *y, GpPoint **result);

/* Returns number of segments processed this pass (<=maxSegs) */
extern long GpIntSegs(const GpXYMap *map, long maxSegs, long n,
		      const GpReal *x1, const GpReal *y1,
		      const GpReal *x2, const GpReal *y2, GpSegment **result);

/* ------------------------------------------------------------------------ */

/* DrawMarkers based on DrawText.  Uses Helvetica font, for use when:
   (1) A device has no polymarker primitive, or
   (2) You want the marker to be a character (used by GaLines)   */
/* Note: GpPseudoMark is defined in gist.c to share static functions.  */
extern int GpPseudoMark(Engine *engine, long n,
			const GpReal *px, const GpReal *py);

/* Routine to clip cell array (one axis at a time) */
extern long GpClipCells(GpMap *map, GpReal *px, GpReal *qx,
			GpReal xmin, GpReal xmax, long ncells, long *off);

/* Raw routine to inflict damage */
extern void GpDamage(Engine *eng, Drawing *drawing, GpBox *box);

/* ------------------------------------------------------------------------ */

#endif
