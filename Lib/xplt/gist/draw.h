/*
 * DRAW.H
 *
 * $Id$
 *
 * Declare display list structures of GIST C interface
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef DRAW_H
#define DRAW_H

#include "gist.h"

/* Virtual function table for display list objects */
typedef struct GdOpTable GdOpTable;
struct GdOpTable {
  int type;    /* element type-- E_LINES, ... E_SYSTEM */

  /* Kill is virtual destructor (constructors are GdLines, etc.)
     If all is non-0, then every element in this ring is killed */
  void (*Kill)(void *el);

  /* GetProps sets gistD and gistA from element, returns element type */
  int (*GetProps)(void *el);

  /* SetProps sets element properties from gistA, gistD
     If xyzChanged has CHANGE_XY or CHANGE_Z or both set, then
     log x and y arrays are freed or contours are recomputed or both.
     Returns 0 on success or 1 if memory error recomputing contours.  */
  int (*SetProps)(void *el, int xyzChanged);

  /* Draw generates appropriate drawing primitives for this element.
     Returns 0 if successful.  */
  int (*Draw)(void *el, int xIsLog, int yIsLog);

  /* Scan rescans the (x,y) extent of this element according to the
     flags (extreme values, log scaling) and limits, setting extreme
     limits as appropriate for this curve.  If the limits are restricted
     and if no points will be plotted, limits.xmin=limits.xmax or
     limits.ymin=limits.ymax on return.  The el->box is also set
     appropriately.  */
  int (*Scan)(void *el, int flags, GpBox *limits);

  /* Margin returns a box representing the amount of NDC space which
     should be left around the element's box to allow for finite line
     thickness, projecting text, and the like. */
  void (*Margin)(void *el, GpBox *margin);
};

/* Generic display list element */
typedef struct GdElement GdElement;
struct GdElement {
  GdOpTable *ops;  /* virtual function table */
  GdElement *next, *prev;  /* elements form doubly linked rings */
  GpBox box;       /* extreme values of coordinates for this object */
  int hidden;      /* hidden flag */
  char *legend;    /* descriptive text */
  int number;      /* drawing->nElements when this element added */
};

/* GdLines element */
typedef struct GeLines GeLines;
struct GeLines {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  int n;                        /* number of points */
  GpReal *x, *y, *xlog, *ylog;  /* (x,y) points */
  GpLineAttribs l;              /* line attributes */
  GaLineAttribs dl;             /* decorated line attributes */
  GpMarkerAttribs m;            /* marker attributes */
};

/* GdDisjoint element */
typedef struct GeDisjoint GeDisjoint;
struct GeDisjoint {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  int n;                        /* number of segments */
  GpReal *x, *y, *xlog, *ylog;  /* (px,py) points */
  GpReal *xq, *yq, *xqlog, *yqlog;  /* (qx,qy) points */
  GpLineAttribs l;              /* line attributes */
};

/* GdText element */
typedef struct GeText GeText;
struct GeText {
  GdElement el;
  GpReal x0, y0;                /* text position */
  char *text;                   /* text string */
  GpTextAttribs t;              /* text attributes */
};

/* GdCells element */
typedef struct GeCells GeCells;
struct GeCells {
  GdElement el;
  GpReal px, py, qx, qy;        /* colors[0][width-1] at (qx,py) */
  long width, height;           /* width and height of image */
  GpColor *colors;              /* image array */
};

/* GdFill element */
typedef struct GePolys GePolys;
struct GePolys {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  GpReal *x, *y, *xlog, *ylog;  /* (px,py) points */
  long n, *pn;                  /* pn[n] list of polygon lengths */
  GpColor *colors;              /* n fill colors */
  GpLineAttribs e;              /* edge attributes */
};

/* GdMesh element */
typedef struct GeMesh GeMesh;
struct GeMesh {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  int noCopy;                   /* also has bit for mesh.reg */
  GaQuadMesh mesh;              /* mesh coordinates */
  GpReal *xlog, *ylog;          /* 0 unless log axes */
  int region, boundary;         /* region number, boundary flag */
  GpLineAttribs l;              /* line attributes */
  int inhibit;                  /* +1 to inhibit j=const lines,
				   +2 to inhibit i=const lines */
};

/* GdFillMesh element */
typedef struct GeFill GeFill;
struct GeFill {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  int noCopy;                   /* also has bit for mesh.reg */
  GaQuadMesh mesh;              /* mesh coordinates */
  GpReal *xlog, *ylog;          /* 0 unless log axes */
  int region;                   /* this and above must match GeMesh */

  GpColor *colors;              /* color array */
  long nColumns;                /* row dimension of colors */
  GpLineAttribs e;              /* edge attributes */
};

/* GdVectors element */
typedef struct GeVectors GeVectors;
struct GeVectors {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  int noCopy;                   /* also has bit for mesh.reg */
  GaQuadMesh mesh;              /* mesh coordinates */
  GpReal *xlog, *ylog;          /* 0 unless log axes */
  int region;                   /* this and above must match GeMesh */

  GpReal *u, *v;                /* iMax-by-jMax vector components (u,v) */
  GpReal scale;                 /* vectors plot as (dx,dy)=scale*(u,v) */
  GpLineAttribs l;              /* line attributes */
  GpFillAttribs f;              /* filled area attributes */
  GaVectAttribs vect;           /* vector attributes */
};

/* GdContours element */
typedef struct GeContours GeContours;
struct GeContours {
  GdElement el;
  GpBox linBox, logBox;         /* extreme values */
  int noCopy;                   /* also has bit for mesh.reg */
  GaQuadMesh mesh;              /* mesh coordinates */
  GpReal *xlog, *ylog;          /* 0 unless log axes */
  int region;                   /* this and above must match GeMesh */

  GpReal *z;                    /* iMax-by-jMax contour array */
  int nLevels;                  /* number of contour levels */
  GpReal *levels;               /* list of contour levels */
  GeLines **groups;             /* groups of contour curves at each level */
  GpLineAttribs l;              /* line attributes */
  GaLineAttribs dl;             /* decorated line attributes */
  GpMarkerAttribs m;            /* marker attributes */
};

/* Coordinate system-- a GpTransform with a tick style */
typedef struct GeSystem GeSystem;
struct GeSystem {
  GdElement el;
  GaTickStyle ticks;      /* tick style */
  GpTransform trans;      /* viewport and window */
  int flags;              /* for computing the limits (see GdProperties) */
  int rescan;             /* set if must be rescanned before redraw */
  int unscanned;          /* element number of 1st unscanned element */
  GdElement *elements;    /* elements in this coordinate system */
  GpBox savedWindow;      /* saved window limits for GdSaveLimits */
  int savedFlags;         /* saved flags for GdSaveLimits */
  GaAltTicks *xtick, *ytick;    /* alternative tick and label generators */
  GaAltLabel *xlabel, *ylabel;
};

typedef struct GeLegendBox GeLegendBox;
struct GeLegendBox {
  GpReal x, y;              /* NDC location of this legend box */
  GpReal dx, dy;            /* if non-zero, offset to 2nd column */
  GpTextAttribs textStyle;  /* font, size, etc. of these legends */
  int nchars, nlines;       /* max number of characters per line, lines */
  int nwrap;                /* max number of lines to wrap long legends */
};

/* Drawing-- the complete display list */
struct Drawing {
  Drawing *next;          /* drawings kept in a list */
  int cleared;            /* next element added will delete everything */
  int nSystems;           /* total number of systems (now) in drawing */
  int nElements;          /* total number of elements (ever) in drawing */
  GeSystem *systems;      /* coordinate systems */
  GdElement *elements;    /* elements not belonging to any system */
  int damaged;            /* set if damage box meaningful */
  GpBox damage;           /* region damaged by GdEdit, etc. */
  int landscape;          /* non-0 for landscape, 0 for portrait */
  GeLegendBox legends[2]; /* 0 is normal legend box, 1 is contours */
};

/* The list of GIST drawings */
extern Drawing *gistDrawList;

/* The following functions are intended to assist in writing the
   constructors for new types of Drawing Elements */

/* generic function for adding elements to current system--
   note that you must reset el->ops by hand, since only the predefined
   types are treated properly */
extern void GeAddElement(int type, GdElement *element);

/* generic function to mark element as unscanned if in a coordinate
   system, else set its box */
extern void GeMarkForScan(GdElement *el, GpBox *linBox);

/* generic function to get the current mesh, returning a GeMeshXY*,
   and optionally not copying the gistA.mesh arrays.  The return
   value is iMax*jMax (0 on failure) */
extern long GeGetMesh(int noCopy, GaQuadMesh *meshin, int region,
		      void *vMeshEl);

/* Updates the system limits if necessary -- do not use lightly.  */
extern int GdScan(GeSystem *system);

/* Defined in engine.c.  These routines communicate state information
   from the drawing to the engine.  They are not intended for external
   use.  */
extern int GdBeginDr(Drawing *drawing, GpBox *damage, int landscape);
   /* GdBeginSy returns 1 bit to draw elements, 2 bit to draw ticks */
extern int GdBeginSy(GpBox *tickOut, GpBox *tickIn,
		     GpBox *viewport, int number, int sysIndex);
   /* GdBeginEl returns 1 if elements should be drawn, else 0 */
extern int GdBeginEl(GpBox *box, int number);
extern void GdEndDr(void);

#endif
