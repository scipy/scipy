/*
 * GIST.H
 *
 * $Id$
 *
 * Declare GIST interface for C programs
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef GIST_H
#define GIST_H

/* current code breaks if GpReal is float instead of double */
typedef double GpReal;
typedef unsigned char GpColor;

/* Normalized device coordinates (NDC) have a definite meaning to GIST--
   PostScript and printing devices like units of 1/72.27 inch,
   X11 release 4 likes units of 1/75 inch or 1/100 inch (for its fonts).
   Next, 1.0 NDC should correspond roughly to 10.5 inches (the wide
   dimension of an ordinary sheet of paper.  This will happen if
   one point is 0.0013 NDC units (1.0 NDC unit is about 10.64 inch).
   Then NDC origin is approximately in the lower left hand corner of
   the page, although this can be shifted slightly to get nice margins.
 */
#define ONE_POINT 0.0013000
#define INCHES_PER_POINT (1./72.27)
#define ONE_INCH (72.27*ONE_POINT)

extern char gistError[128];  /* most recent error message */

extern char *gistPathDefault;  /* set in Makefile or gread.c, can be
                                  overridden by resetting directly */

/* ------------------------------------------------------------------------ */
/* Initialization Functions */

/*
  An engine in GIST is a set of device drivers.  A particular instance
  of an engine will also include a specific I/O connection for these
  drivers.  Instead of opening a GKS workstation, in GIST, you create
  a specific instance of a particular type of engine -- a PostScript
  engine, a CGM engine, or an X window engine.  Since there is a
  separate function for creating each different type of engine, GIST
  has the important advantage over GKS that only those device drivers
  actually used by your code are loaded.

   Four sets of device drivers are supplied with GIST:
     PostScript - hopefully can be fed to your laser printer
     CGM - binary format metafile, much more compact than PostScript,
           especially if there are many pages of graphics
     BX - basic play window
     FX - play window with odometer bar at top, mouse zooming, etc
 */

typedef struct Engine Engine;
extern Engine *GpPSEngine(char *name, int landscape, int mode, char *file);
extern Engine *GpCGMEngine(char *name, int landscape, int mode, char *file);
extern Engine *GpBXEngine(char *name, int landscape, int dpi, char *display);
extern Engine *GpFXEngine(char *name, int landscape, int dpi, char *display);

extern void GpKillEngine(Engine *engine);

extern int gist_input_hint, gist_private_map, gist_rgb_hint;
extern void g_initializer(int *pargc, char *argv[]);
extern void (*g_on_keyline)(char *msg);
extern void (*g_stdout)(char *output_line);

/* ------------------------------------------------------------------------ */
/* Control Functions */

/* Engines are kept in two lists - a list of all engines and
   a ring of active engines.  To ignore the active list and send
   only to a specific engine, call GpPreempt.  GpPreempt(0) turns
   off preemptive mode.  The preempting engine need not be active;
   it will still get all output until preempt mode is turned off, at
   which time it returns to its original state.  */
extern int GpActivate(Engine *engine);
extern int GpDeactivate(Engine *engine);
extern int GpPreempt(Engine *engine);
extern int GpActive(Engine *engine);  /* 1 if active or preempting, else 0 */

/* Pass engine==0 to GpClear or GpFlush to affect all active engines.  */
extern int GpClear(Engine *engine, int flag);
#define CONDITIONALLY 0
#define ALWAYS 1
extern int GpFlush(Engine *engine);

/* The Next routines provide a way to iterate throught the engine lists.
   Engines will be returned in order they are created or activated.
   Use engine==0 to get the first engine in the list; will return 0
   when there are no more engines in the list.  If preemptive mode has
   been set with GpPreempt, GpNextActive returns the preempting engine,
   whether or not it is in the active list.  */
extern Engine *GpNextEngine(Engine *engine);
extern Engine *GpNextActive(Engine *engine);

/* ------------------------------------------------------------------------ */
/* Transformations and Clipping */

typedef struct GpBox GpBox;
struct GpBox {
  GpReal xmin, xmax, ymin, ymax;
};

typedef struct GpTransform GpTransform;
struct GpTransform {
  GpBox viewport, window;
};

/* GpSetTrans sets the current transformation.  This involves notifying
   each active engine of the change, in turn, so that an appropriate
   transformation to device coordinates can be computed.  (GpActivate
   also sets the device coordinate transformation in the engine.)
   This transformation is also loaded into the gistA attribute list.  */
extern int GpSetTrans(const GpTransform *trans);

/* Although NDC space has a particular meaning, an 8.5x11 sheet of
   paper may have x along the 8.5 inch edge (portrait mode), or x along
   the 11 inch edge (landscape mode).  In either case, the origin is the
   lower left corner of the page.  This property can be set on a per
   Engine basis using GpLandscape, or for all active engines with
   GpLandscape(0).  If you use the D level routines, you shouldn't
   need this.  */
extern int GpLandscape(Engine *engine, int landscape);

/* current transformation */
extern GpTransform gistT;  /* use GpSetTrans to set this properly */

/* Turn clipping on or off by setting gistClip.  */
extern int gistClip;         /* 1 to clip to map.viewport, 0 to not clip */

/* Often, the linear transformation represented by a GpTransform is
   required in the simpler form dst= scale*src+offset.  The following
   convenience routines used internally by GIST are provided:  */
typedef struct GpMap GpMap;
struct GpMap {
  GpReal scale, offset;
};

typedef struct GpXYMap GpXYMap;
struct GpXYMap {
  GpMap x, y;
};

extern void GpSetMap(const GpBox *src, const GpBox *dst, GpXYMap *map);

/* gPortrait and gLandscape contain, respectively, boxes representing
   an 8.5-by-11 inch and an 11-by8.5 inch page in NDC */
extern GpBox gPortrait;
extern GpBox gLandscape;

/* Utilities for GpBox, assuming min<=max for both boxes */
extern int GpIntersect(const GpBox *box1, const GpBox *box2);
extern int GpContains(const GpBox *box1, const GpBox *box2); /* box1>=box2 */
extern void GpSwallow(GpBox *preditor, const GpBox *prey);

/* ------------------------------------------------------------------------ */
/* Output Primitives */

extern int GpLines(long n, const GpReal *px, const GpReal *py);
extern int GpMarkers(long n, const GpReal *px, const GpReal *py);
extern int GpText(GpReal x0, GpReal y0, const char *text);
  /* WARNING- for GpText, (x0,y0) are in WC, but the text size and
              orientation are specified in NDC (unlike GKS).  */
extern int GpFill(long n, const GpReal *px, const GpReal *py);
extern int GpCells(GpReal px, GpReal py, GpReal qx, GpReal qy,
                   long width, long height, long nColumns,
                   const GpColor *colors);

/* GKS seems to be missing a disjoint line primitive... */
extern int GpDisjoint(long n, const GpReal *px, const GpReal *py,
                      const GpReal *qx, const GpReal *qy);

/* ------------------------------------------------------------------------ */
/* Output attributes */

/* Instead of a set and get function pair for each attribute, the GIST
   C interface simply provides global access to its state table.  */
typedef struct GaAttributes GaAttributes;

typedef struct GpLineAttribs GpLineAttribs;
struct GpLineAttribs {
  unsigned long color;
  int type;       /* line types given by L_SOLID, etc.  */
  GpReal width;   /* default 1.0 is normal width of a line */

#define L_NONE 0
#define L_SOLID 1
#define L_DASH 2
#define L_DOT 3
#define L_DASHDOT 4
#define L_DASHDOTDOT 5

#define DEFAULT_LINE_WIDTH (0.5*ONE_POINT)
#define DEFAULT_LINE_INCHES (0.5*INCHES_PER_POINT)
};

typedef struct GpMarkerAttribs GpMarkerAttribs;
struct GpMarkerAttribs {
  unsigned long color;
  int type;       /* marker types given by M_ASTERISK, etc., or by ASCII  */
  GpReal size;    /* default 1.0 is normal size of a marker */

#define M_POINT 1
#define M_PLUS 2
#define M_ASTERISK 3
#define M_CIRCLE 4
#define M_CROSS 5

#define DEFAULT_MARKER_SIZE (10.0*ONE_POINT)
};

typedef struct GpFillAttribs GpFillAttribs;
struct GpFillAttribs {
  unsigned long color;
  int style;      /* fill area styles given by F_SOLID, etc.  */

#define F_HOLLOW 0
#define F_SOLID 1
#define F_PATTERN 2
#define F_HATCH 3
#define F_EMPTY 4
};

typedef struct GpTextAttribs GpTextAttribs;
struct GpTextAttribs {
  unsigned long color;
  int font;         /* text font (T_HELVETICA, etc.) */
  GpReal height;    /* character height in NDC, default 0.0156 (12pt)
                       UNLIKE GKS, GIST font sizes are always specified
                       in NDC.  This drastically simplifies coding for
                       devices like X windows, which must load a font
                       at each size.  It also conforms better with
                       a Mac-like user interface in which font size
                       in points is selected by the user.  */
  int orient;          /* text paths given by TX_RIGHT, etc.  */
  int alignH, alignV;  /* text alignments given by TH_NORMAL, etc.  */

  /* GKS is missing a text opacity flag.  */
  int opaque;

/* A font is a type face optionally ORed with T_BOLD and/or T_ITALIC. */
/* Available point sizes (for X) are 8, 10, 12, 14, 18, and 24 */
#define T_BOLD 1
#define T_ITALIC 2
#define T_COURIER 0
#define T_TIMES 4
#define T_HELVETICA 8
#define T_SYMBOL 12
#define T_NEWCENTURY 16

#define TX_RIGHT 0
#define TX_UP 1
#define TX_LEFT 2
#define TX_DOWN 3

#define TH_NORMAL 0
#define TH_LEFT 1
#define TH_CENTER 2
#define TH_RIGHT 3

#define TV_NORMAL 0
#define TV_TOP 1
#define TV_CAP 2
#define TV_HALF 3
#define TV_BASE 4
#define TV_BOTTOM 5
};

/* GaLines output function supports polylines with occasional marker
   characters and/or ray arrows.  */
typedef struct GaLineAttribs GaLineAttribs;
struct GaLineAttribs {
  int closed;   /* 0 for open curve, 1 for closed curve */
  int smooth;   /* 0 for no smoothing, 1, 2, or 3 for progresively
                   smoother splines (not implemented for all Engines)
                   Note: If marks or rays, smooth is ignored.  */

  /* Note: GaLineAttribs and GaMarkerAttribs determine the style, size,
     and color of the line and occasional markers.  If l.style is L_NONE,
     polymarkers will be plotted at every point rather than occasional
     markers.  */
  int marks;    /* 0 if no occasional markers, 1 if occasional markers */
  GpReal mSpace, mPhase;    /* occasional marker spacing and phase in NDC */

  int rays;     /* 0 if no ray arrows, 1 to get ray arrows */
  GpReal rSpace, rPhase;    /* ray spacing and phase in NDC */
  GpReal arrowL, arrowW;    /* ray arrowhead size, 1.0 is normal arrow */

#define DEFAULT_ARROW_LENGTH (10.0*ONE_POINT)
#define DEFAULT_ARROW_WIDTH (4.0*ONE_POINT)
};

typedef struct GaVectAttribs GaVectAttribs;
struct GaVectAttribs {
  int hollow;      /* 1 to outline darts, 0 to fill
                      (uses current line or filled area attibutes) */
  GpReal aspect;   /*  half-width/length of dart arrows */
};

struct GaAttributes {

  /* line attributes (GpLines, GpDisjoint) */
  GpLineAttribs l;

  /* marker attributes (GpMarkers) */
  GpMarkerAttribs m;

  /* filled area attributes (GpFill) */
  GpFillAttribs f;

  /* text attributes (GpText) */
  GpTextAttribs t;

  /* decorated line attributes (GaLines) */
  GaLineAttribs dl;

  /* vector attributes (GpVectors) */
  GaVectAttribs vect;

  /* edge attributes -- intended for 3D extensions (GpFill) */
  GpLineAttribs e;

  int rgb;  /* for GpCells */
};

extern GaAttributes gistA;

typedef struct GaAxisStyle GaAxisStyle;
struct GaAxisStyle {
#define TICK_LEVELS 5
  GpReal nMajor, nMinor, logAdjMajor, logAdjMinor;
  int nDigits, gridLevel;
  int flags;   /* TICK_L, ... LABEL_L, ... GRID_F below */

  GpReal tickOff, labelOff;  /* offsets in NDC from the edge of the
                                viewport to the ticks or labels */
  GpReal tickLen[TICK_LEVELS];  /* tick lengths in NDC */

  GpLineAttribs tickStyle, gridStyle;
  GpTextAttribs textStyle;   /* alignment ignored, set correctly */
  GpReal xOver, yOver;       /* position for overflow label */

/* Flags determine whether there are ticks at the lower or upper
   (left or bottom is lower) edges of the viewport, whether the ticks
   go inward or outward, whether the lower or upper edge ticks have
   labels, and whether there is a full grid, or a single grid line
   at the origin.
   Also whether an alternative tick and/or label generator is used.  */
#define TICK_L 0x001
#define TICK_U 0x002
#define TICK_C 0x004
#define TICK_IN 0x008
#define TICK_OUT 0x010
#define LABEL_L 0x020
#define LABEL_U 0x040
#define GRID_F 0x080
#define GRID_O 0x100
#define ALT_TICK 0x200
#define ALT_LABEL 0x400
};

typedef struct GaTickStyle GaTickStyle;
struct GaTickStyle {
  GaAxisStyle horiz, vert;
  int frame;
  GpLineAttribs frameStyle;
};

/* ------------------------------------------------------------------------ */

/*
  GIST includes a few output routines at a higher level than GKS.

  These cover the traditional types of plots available as output from
  the LASNEX laser fusion code, a large physics simulation code, plus
  a filled mesh plot.
 */

extern int GaLines(long n, const GpReal *px, const GpReal *py);
       /* Like GpLines, but includes GaAttributes fancy line attributes
          as well as the line attributes from GpAttributes.  */

/* You must load the components of a mesh into a GaMeshXY data structure
   before using the A level routines that require a mesh.  Only the
   contour routines use the triangle array.  If you don't supply a
   reg array (mesh->reg==0), GaMesh, GaFillMesh GaVectors, or
   GaContourInit will supply a default region number array for you.
   YOU ARE RESPONSIBLE FOR RELEASING THE ASSOCIATED STORAGE (using free).
   Failure to supply a triangle array may result in crossing contours
   in GaContours, but GaContourInit will not allocate one.  */
typedef struct GaQuadMesh GaQuadMesh;
struct GaQuadMesh {
  long iMax, jMax; /* mesh logical dimensions */
  GpReal *x, *y;   /* iMax-by-jMax mesh coordinate arrays */
  int *reg;        /* iMax-by-(jMax+1)+1 mesh region number array:
                      reg[j][i] non-zero means that the zone bounded by
                        [j-1][i-1], [j-1][i], [j][i], and [j][i-1]
                        exists,
                      reg[0][i]= reg[j][0]= reg[jMax][i]= reg[j][iMax]= 0,
                        so that the (iMax-1)-by-(jMax-1) possible zones
                        are surrounded by non-existent zones.  */
  short *triangle; /* iMax-by-jMax triangulation marker array
                      triangle[j][i]= 1, -1, or 0 as the zone bounded by
                        [j-1][i-1], [j-1][i], [j][i], and [j][i-1]
                      has been triangulated from (-1,-1) to (0,0),
                      from (-1,0) to (0,-1), or has not yet been
                      triangulated.  */
};

extern int GaMesh(GaQuadMesh *mesh, int region, int boundary, int inhibit);
       /* Plots the quadrilateral mesh.  If boundary==0,
          a mesh line is plotted whenever either zone it borders
          belongs to the given region.  An exception is made if
          region==0; in this case the entire mesh is plotted.
          If boundary!=0, a mesh line is plotted whenever one but
          not both of its bordering zones belong to the given region.
          Again, if region==0, the boundary for the whole mesh is
          plotted.
          If inhibit&1, then lines at constant j are not drawn;
          if inhibit&2, then lines at constant i are not drawn.  */

extern int GaFillMesh(GaQuadMesh *mesh, int region, const GpColor *colors,
                      long nColumns);
       /* Fills each zone of the quadrilateral mesh according to the
          (iMax-1)-by-(jMax-1) array of colors.  The colors array may
          be a subarray of a larger rectangular array; therefore, you
          must specify the actual length of a row of the colors array
          using the nColumns parameter (nColumns>=iMax-1 for sensible
          results).  Only the region specified is plotted, unless
          region==0, in which case all non-zero regions are filled.
          (As a special case, if colors==0, every zone is filled with
          gistA.f.color.)  If gistA.e.type!=L_NONE, an edge is drawn
          around each zone as it is drawn.  */

extern int GaFillMarker(long n, const GpReal *px, const GpReal *py,
                        GpReal x0, GpReal y0);
       /* Fills (px,py)[n] in NDC with origin at (x0,y0) in world.  */

extern int GaVectors(GaQuadMesh *mesh, int region,
                     const GpReal *u, const GpReal *v, GpReal scale);
       /* Plot a vector scale*(u,v) at each point of the current
          mesh (gistA.mesh).  If region==0, the entire mesh is
          drawn, otherwise only the specified region.  */

extern int GaContourInit(GaQuadMesh *mesh, int region,
                         const GpReal *z, GpReal level);
       /* Find the edges cut by the current contour, remembering z, mesh
          region, and level for the GaContour routine, which actually
          walks the contour.  The z array represents function values at
          the mesh points.  If region==0, contours are drawn on the
          entire mesh, otherwise only on the specified region.
          The mesh->triangle array is updated by GaContour as follows:
          If a contour passes through an untriangulated saddle zone,
          GaContour chooses a triangulation and marks the triangle
          array approriately, so that subsequent calls to GaContour
          will never produce intersecting contour curves.
          If mesh->triangle==0 on input to GaContourInit, the
          curves returned by GaContour may intersect the curves
          returned after a subsequent call to GaContourInit with a
          new contour level.  */

extern int GaContour(long *cn, GpReal **cx, GpReal **cy, int *closed);
       /* After a call to GaContourInit, GaContour must be called
          repeatedly to generate the sequence of curves obtained by
          walking the edges cut by the contour level plane.  GaContour
          returns 1 until there are no more contours to be plotted,
          when it returns 0.  GaContour signals an error by returning
          0, but setting *cn!=0.  The curve coordinates (*cx,*cy)
          use internal scratch space and the associated storage must
          not be freed.  (*cx, *cy) are valid only until the next
          GaContour or GaMesh call.  */

extern int GaTicks(GaTickStyle *ticks, int xIsLog, int yIsLog);
       /* Draws a system of tick marks and labels for the current
          transformation (gistT), according to the specified
          tick style.  */

extern int GaFreeScratch(void);
       /* Frees scratch space required by GaMesh, GaContour,
          and GaContourInit, and GaTicks.  Ordinarily, there is no reason
          to call this, as the amount of scratch space required is modest
          in comparison to the size of the mesh, and it is reused
          in subsequent calls.  */

typedef int GaAltTicks(GpReal lo, GpReal hi, GpReal nMajor, GpReal nMinor,
                       GpReal *ticks, int nlevel[TICK_LEVELS]);
      /* An alternative tick generator function accepts the lo and hi
         axis limits (lo<hi) and the maximum total number of ticks
         allowed for all levels of the hierarchy.  It is given space
         ticks[ceil(nMinor)] to return the tick values, which are given
         in order within each level -- that is, all level 0 ticks first,
         followed by all level 1 ticks, then level 2, and so on.  The
         elements of nlevel array should be set to the number of ticks
         ticks returned in ticks up to that level.  On entry, nlevel[i]
         is set to zero for each i.  On exit, max(nlevel)<=maxnum.
         The function should return 0 if successful.  If it returns
         non-zero, the default tick generation scheme will be used.  */
typedef int GaAltLabel(char *label, GpReal value);
      /* An alternative label generator stores label (<=31 characters)
         in label, which is to represent the given value.  The function
         should return 0 on success, non-zero to drop back to the default
         scheme.  If label==0 on input, the function should return the
         correct success value without actually computing the label;
         it will be called once in this mode for every value, and if it
         fails for any value, the default scheme will be used instead.  */

extern int Base60Ticks(GpReal lo, GpReal hi, GpReal nMajor, GpReal nMinor,
                       GpReal *ticks, int nlevel[TICK_LEVELS]);
     /* Puts ticks at multiples of 30, failing if lo<=-3600 or hi>=+3600.
        For ticks at multiples of 10 or less, the subdivisions are
        identical to the default decimal tick scheme.  */
extern int DegreeLabels(char *label, GpReal value);
     /* Prints (value+180)%360-180 instead of just value.  */
extern int HourLabels(char *label, GpReal value);
     /* Prints hh:mm (or mm:ss) for value 60*hh+mm (or 60*mm+ss).  */

extern int GaAltTick(GaTickStyle *ticks, int xIsLog, int yIsLog,
                      GaAltTicks *xtick, GaAltLabel *xlabel,
                      GaAltTicks *ytick, GaAltLabel *ylabel);
       /* Like GaTicks, but uses the specified ticker and labeler for
          non-log axes.  Log axes always use the defaults.  Either
          ticker or labeler or both may be zero to use the default.  */

/* ------------------------------------------------------------------------ */

/* These general purpose contouring routines actually have nothing to
 * do with graphics; I provide them as a convenience for computing
 * inputs to the GaLines and GpFill routines.
 * GcInit1 takes the same inputs as GaContourInit.  It returns nparts,
 * the number of disjoint contour curves, and its return value is the
 * total number of points on all these parts.  (On closed curves, the
 * start point will be repeated, and so has been counted twice.)
 * GcInit2 accepts two level values instead of one; it returns the
 * boundary of the region between these two contour levels, with
 * slits cut as necessary to render each disjoint part simply connected,
 * traced counterclockwise (in logical space) about the region between
 * the levels.  If you specify a non-zero value for nchunk, then the
 * mesh will be chunked into pieces of no more than nchunk^2 zones, which
 * limits the number of sides of any returned polygon to about 4*nchunk^2.
 * After you call GcInit1 or GcInit2, you must allocate px[ntotal],
 * py[ntotal], and n[nparts] arrays, then call GcTrace to fill them.
 */
extern long GcInit1(GaQuadMesh *mesh, int region, const GpReal *zz,
                    GpReal lev, long *nparts);
extern long GcInit2(GaQuadMesh *mesh, int region, const GpReal *zz,
                    GpReal levs[2], long nchunk, long *nparts);
extern long GcTrace(long *n, GpReal *px, GpReal *py);

/* ------------------------------------------------------------------------ */
/* Display list routines */

/* GIST maintains a list of drawings much like its list of engines.
   Unlike engines, only one drawing is active at a time.  The
   contents of the drawing can be rendered on all active engines by
   calling GdDraw.  This clears each engine before it begins to draw,
   so that each call to GdDraw results in one page of graphics output
   on all active engines.  X window engines may remember the most recent
   drawing rendered on them; they will then respond to the GdUpdate
   command by repairing any damage or making any additions to the
   window containing the drawing.  On metafile engines, GdUpdate is
   a noop.
   A drawing consists of a list of coordinate systems (a GKS normalization
   transformation plus GIST ticks and labels and a list of elements to be
   displayed), plus a list of graphical elements outside any coordinate
   systems (such as plot titles or legends).  The general layout of the
   coordinate systems (their viewport and tick and label style) can be
   specified using GdNewSystem, but the preferred technique is to put
   this data into a "style sheet", which is an ASCII file.  */

typedef struct Drauing Drauing;

extern Drauing *GdNewDrawing(char *gsFile);
extern void GdKillDrawing(Drauing *drawing);

/* After GdNewDrawing, the new drawing becomes the current drawing.
   GdSetDrawing can be used to change the current drawing.  This
   action saves the state of the previous drawing (the current system,
   element, etc.), which can be recalled using GdSetDrawing(0).
   Otherwise, GdSetDrawing scans the drawing to find the latest
   system, element, and mesh, and these become current.  */
extern int GdSetDrawing(Drauing *drawing);

/* GdClear marks the drawing as cleared, but no action is taken until
   something new is drawn, at which point the cleared display list
   elements are discarded.  GpClear is not sent until the next
   GdDraw.  If drawing is 0, the current drawing is assumed.  */
extern int GdClear(Drauing *drawing);

/* If changesOnly is non-zero, only new elements or damaged areas
   are redrawn.  The amount drawn may be recorded by each engine, and
   may vary from engine to engine.
   The Drauing* is also recorded on each engine.  */
extern int GdDraw(int changesOnly);

/* Graphical elements are added to the current drawing by setting
   attributes in gistA, then calling GdLines, GdDisjoint, GdCells, GdMesh,
   GdFillMesh, GdVectors, or GdContours.  GdText puts the text outside
   all coordinate systems unless toSys is non-0.  Each of these routines
   returns a unique identification number (equal to the total number
   of objects defined in the current drawing since the last GdClear).
   They return -1 if an error occurs.  */

extern int GdLines(long n, const GpReal *px, const GpReal *py);
extern int GdDisjoint(long n, const GpReal *px, const GpReal *py,
                      const GpReal *qx, const GpReal *qy);
extern int GdText(GpReal x0, GpReal y0, const char *text, int toSys);
extern int GdCells(GpReal px, GpReal py, GpReal qx, GpReal qy,
                   long width, long height, long nColumns,
                   const GpColor *colors);
extern int GdFill(long n, const GpColor *colors, const GpReal *px,
                  const GpReal *py, const long *pn);

/* The other D level primitives involve mesh arrays.  It is often
   convenient NOT to copy mesh data, so that several objects may
   point to a single mesh array (strictly speaking, it is the mesh->x,
   mesh->y, mesh->reg, and mesh->triangle arrays which will be copied
   or not copied).  Therefore, a noCopy flag is provided, which can be
   constructed by ORing together NOCOPY_MESH, etc.  If you set any of
   these bits, you are promising that you will not free the corresponding
   object for the lifetime of this display list object, conversely,
   when Gist discards the object, it will not free the pointers (the
   memory management responsibility is yours).  However, if you
   supply a GdFree routine, Gist will call that for uncopied objects.  */
extern int GdMesh(int noCopy, GaQuadMesh *mesh, int region, int boundary,
                  int inhibit);
extern int GdFillMesh(int noCopy, GaQuadMesh *mesh, int region,
                      GpColor *colors, long nColumns);
extern int GdVectors(int noCopy, GaQuadMesh *mesh, int region,
                     GpReal *u, GpReal *v, GpReal scale);
extern int GdContours(int noCopy, GaQuadMesh *mesh, int region,
                      GpReal *z, const GpReal *levels, int nLevels);
#define NOCOPY_MESH 1
#define NOCOPY_COLORS 2
#define NOCOPY_UV 4
#define NOCOPY_Z 8

/* graphical element types */
#define E_NONE 0
#define E_LINES 1
#define E_DISJOINT 2
#define E_TEXT 3
#define E_MESH 4
#define E_FILLED 5
#define E_VECTORS 6
#define E_CONTOURS 7
#define E_CELLS 8
#define E_POLYS 9
#define E_SYSTEM 10

/* Properties of drawing elements are similar to attributes of
   graphical primitives.  The properties include everything passed
   in the parameter lists to the primitives.  The global property
   list gistD is used to inquire about or change these values in
   the drawing elements-- i.e.- to edit the drawing.  */

typedef struct GdProperties GdProperties;
struct GdProperties {
  /* Properties of all elements */
  int hidden;           /* 1 if element should not be plotted */
  char *legend;         /* description of this element */

  /* Properties of coordinate systems */
  GaTickStyle ticks;    /* for GaTicks to produce ticks and labels */
  GpTransform trans;    /* transform for GpSetTrans for current system */
  int flags;            /* flags for computing the limits */
#define D_XMIN 0x001
#define D_XMAX 0x002
#define D_YMIN 0x004
#define D_YMAX 0x008
#define D_RESTRICT 0x010
#define D_NICE 0x020
#define D_SQUARE 0x040
#define D_LOGX 0x080
#define D_LOGY 0x100
#define D_ZOOMED 0x200
  GpBox limits;         /* current trans->window or exp10(trans->window) */

  /* Properties of elements */
  /* --- GdLines- also uses gistA.l, gistA.dl, gistA.m, GdFill */
  int n;
  GpReal *x, *y;
  /* --- GdDisjoint */
  GpReal *xq, *yq;
  /* --- GdText- also uses gistA.t */
  GpReal x0, y0;
  char *text;
  /* --- GdCells */
  GpReal px, py, qx, qy;
  long width, height;
  /* --- GdCells, GdFillMesh */
  long nColumns;         /* if colors copied, this is width or iMax-1 */
  GpColor *colors;       /* also GdFill */
  /* --- GdMesh, GdFillMesh, GdVectors, GdContours */
  int noCopy;  /* if bit is set, Gist will not free corresponding array
                  -- however, if non-0, GdFree will be called */
/* if NOCOPY_MESH was set, but reg or triangle was 0, then gist
   owns the reg or triangle array, and the following bits are NOT set: */
#define NOCOPY_REG 16
#define NOCOPY_TRI 32
  GaQuadMesh mesh;
  int region;
  /* --- GdMesh- also uses gistA.l */
  int boundary, inhibit;
  /* --- GdVectors- also uses gistA.l, gistA.f, gistA.vect */
  GpReal *u, *v, scale;
  /* --- GdContours- individual level curves are like GdLines
          contour itself uses gistA.l, gistA.dl, gistA.m as defaults */
  GpReal *z, *levels;
  int nLevels;
  /* --- GdFill */
  long *pn;
};

extern GdProperties gistD;

/* GdSetSystem, GdSetElement, and GdSetContour cause the contents
   of one drawing element to be copied into gistD and gistA for
   examination or editing (with GdEdit or GdRemove).
   GdNewSystem does an automatic GdSetSystem, and GdLines,
   GdText, GdCells, GdMesh, GdFillMesh, GdVectors, and GdContours
   do an automatic GdSetElement to the newly created element.
   The return values are E_LINES, ... E_SYSTEM; E_NONE on failure.  */
extern int GdSetSystem(int sysIndex);
extern int GdSetElement(int elIndex);
extern int GdSetContour(int levIndex);

/* return current sysIndex, or <0 if unavailable */
extern int GdGetSystem(void);

/* The elIndex required by GdSetElement is NOT the id number returned
   by GdLines, GdMesh, etc.  Instead, elIndex begins at 0 and goes
   up to 1 less than the number of elements in the current system.
   To get the elIndex corresponding to a given id, use GdFindIndex(id).
   The sysIndex varies from 0 (outside of all systems), up to the
   total number of coordinate systems defined in the current Drauing.
   To find the sysIndex for a given element id, use GdFindSystem(id).  */
extern int GdFindIndex(int id);  /* returns -1 if no such element
                                    if current system */
extern int GdFindSystem(int id); /* returns -1 if no such element
                                    anywhere in current Drauing */

/* After you call one of the 3 GdSet... routines, you may
   alter the contents of gistA and/or gistD, then call
   GdEdit to write the changes back into the element.
   If you change gistD.x, gistD.y, gistD.mesh->x, or gistD.mesh->y,
   use GdEdit(CHANGE_XY) (this recomputes log coordinates if necessary),
   otherwise, use GdEdit(0).  If you want to change any of the
   pointers, you should first free the existing pointer (with free),
   and be aware that GdKillDrawing, GdClear, or GdRemove will
   free the pointer you are providing.  If the element is a set of
   contours, and you change gistD.z or gistD.levels, gistD.nLevels, or
   gistD.mesh->triangle, use GdEdit(CHANGE_Z).
   Use GdEdit(CHANGE_XY | CHANGE_Z) if both apply.
   GdRemove completely removes the element from the drawing.  */
extern int GdEdit(int xyzChanged);
#define CHANGE_XY 1
#define CHANGE_Z 2
extern int GdRemove(void);

/* GdGetLimits copies the current coordinate systems limits into
   gistD.limits and gistD.flags.  GdSetLimits sets the limits in
   the current coordinate system to gistD.limits and gistD.flags
   (if the various adjustment flags are set, the window limits
   may be adjusted, but this will not happen until the drawing
   is rendered on at least one engine).  GdSetPort sets the tick
   style and viewport for the current coordinate system.  */
extern int GdGetLimits(void);
extern int GdSetLimits(void);
extern int GdSetPort(void);

/* Each coordinate system keeps a "saved" set of limits, which
   can be "reverted" to return the limits to the state at the time
   of the last save operation.  (gistD plays no role in this operation.)
   This is used by the fancy X Engine to save the limits prior to
   a series of mouse-driven zoom or pan operations. 
   GdRevertLimits(1) reverts to the save limits only if the current
   limits are marked D_ZOOMED, while GdRevertLimits(0) reverts
   unconditionally.  GdSaveLimits(1) assures that D_ZOOMED is not
   set in the saved flags.  */
extern int GdSaveLimits(int resetZoomed);
extern int GdRevertLimits(int ifZoomed);

/* You can register an alternative tick and/or label generator for
   each axis of the current corodinate system.  They will not be used
   unless the ALT_TICK or ALT_LABEL limits flag is set.  Passing 0
   to GdAltTick leaves the corresponding function unchanged -- there
   is no way to unregister it (just turn off the limits flag).  */
extern int GdAltTick(GaAltTicks *xtick, GaAltLabel *xlabel,
                     GaAltTicks *ytick, GaAltLabel *ylabel);

/* Unlike the other `D' level drawing routines, GdDrawLegends acts
   immediately.  If engine==0, all active engines are used.
   The intent is that legends not be rendered in an X window or
   other interactive device, but this optional routine can render
   them on a hardcopy device.  */
extern int GdDrawLegends(Engine *engine);

/* ------------------------------------------------------------------------ */
/* The following Gd routines are not intended for ordinary use:
   GdReadStyle, GdNewSystem, GdLandscape, and GdLegendBox are
   unnecessary if you use drawing style files (see gread.c).
   GdDetach is provided here for completeness; several of the
   other Gd routines use it (e.g.- GdClear).
   GdClearSystem is required for the specialized animation mode
   described in hlevel.c.  */

/* GdReadStyle clears the drawing and removes all its current systems
   as a side effect.  It is the worker for GdNewDrawing.  */
extern int GdReadStyle(Drauing *drawing, const char *gsFile);

/* Coordinate systems are usually defined in the styleFile, but
   a new coordinate system may be added with GdNewSystem, which
   returns the index of the newly created system.  */
extern int GdNewSystem(GpBox *viewport, GaTickStyle *ticks);

/* Set current drawing to landscape or portrait mode */
extern int GdLandscape(int landscape);

/* Set location and properties of legends for output.   */
extern int GdLegendBox(int which, GpReal x, GpReal y, GpReal dx, GpReal dy,
                       const GpTextAttribs *t, int nchars, int nlines,
                       int nwrap);

/* GdDetach detaches one or all (if drawing==0) drawings from one or all
   (if engine==0) engines.  Whether an engine is active or not is
   unimportant.  A drawing is attached to an engine by GdDraw.  */
extern void GdDetach(Drauing *drawing, Engine *engine);

/* Clear the current coordinate system-- returns 0 if there is no
   current system.  Damages viewport only if all limits fixed,
   otherwise damages both ticks and viewport.  Returns damaged box.  */
extern GpBox *GdClearSystem(void);

/* ------------------------------------------------------------------------ */
/* Color */

/* GIST cell arrays and filled meshes are aimed at displaying pseudocolor
   images, not true color images.  (Another primitive might be added to
   handle true color images.  This is probably most useful for shading
   projections of 3D objects, which requires additional primitives as
   anyway.)  The colors arrays in the G*Cells and G*FillMesh primitives
   are arrays of indices into a pseudocolor map.  The colors for the
   various attribute lists may come from this pseudocolor map, or
   they may be special values representing the foreground and background
   colors, or a primary color.  */

/* duplicate P_BG, P_FG, etc from play.h */
#define BG_COLOR (255UL)
#define FG_COLOR (254UL)
#define BLACK_COLOR (253UL)
#define WHITE_COLOR (252UL)
#define RED_COLOR (251UL)
#define GREEN_COLOR (250UL)
#define BLUE_COLOR (249UL)
#define CYAN_COLOR (248UL)
#define MAGENTA_COLOR (247UL)
#define YELLOW_COLOR (246UL)

/* Each Engine maintains a list of colors for the pseudocolor
   map.  These can be set or retrieved directly from the nColors and
   palette Engine members.  The pseudocolor map is in addition to
   the 10 standard colors defined in gist.h -- note that the meaning of
   a colorcell need not be the same from one engine to another, although
   the 4th component (gray) has been added to make it easier to control
   both gray and color devices using the same GpColorCell array.  Note
   that no Engine owns a GpColorCell array; allocating and freeing this
   storage is the responsibility of the application program.

   A change in the colormap for an Engine will generally have no effect
   until the next page of graphics output.  However, individual types
   of engines may implement a function to try to make immediate
   changes to the pseudocolor map.  */

/* duplicate p_col_t from play.h */
typedef unsigned long GpColorCell;

/* Two routines are provided to set the gray component of the cells in
   a colormap--  GpPutGray for gray=(red+green+blue)/3 and GpPutNTSC
   for gray=(30*red+59*green+11*blue).  The latter is the NTSC standard
   for converting color TV signals to monochrome.  GpPutRGB copies the
   gray entry into the red, green, and blue entries.  */
/* these all no-ops now */
extern void GpPutGray(int nColors, GpColorCell *palette);
extern void GpPutNTSC(int nColors, GpColorCell *palette);
extern void GpPutRGB(int nColors, GpColorCell *palette);

/* Set the palette; the change takes place as soon as possible for
   this engine and the effect on previously drawn objects is
   engine dependent.  The maximum usable palette size is returned.  */
extern int GpSetPalette(Engine *engine, GpColorCell *palette, int nColors);

/* GpGetPalette returns the number of colors in the palette-
   this is the nColors passed to GpSetPalette.  */
extern int GpGetPalette(Engine *engine, GpColorCell **palette);

/* GpReadPalette returns the number of colors found in the palette
   file (see gread.c for format), as well as the palette itself.
   If the number of colors in the palette exceeds maxColors, attempts
   to scale to maxColors.  */
extern int GpReadPalette(Engine *engine, const char *gpFile,
                         GpColorCell **palette, int maxColors);

/* The PostScript and CGM formats, unlike Xlib, guarantee that
   the color table from one page does not automatically carry over to
   the next.  Additionally, color images cannot be rendered (and may be
   inconvenient if they can be rendered) on most hardcopy devices.
   Consequently, by default, Gist produces a simple gray scale for
   PostScript or CGM output.  However, you can optionally force the
   current palette to be dumped on each output page by setting the
   color mode flag.  If any marks have been made on the current page,
   the color table cannot be dumped until the next page.
   For an X Engine, the colorMode is initially 0, which means to
   use the best available shared colors.  Setting colorMode to 1
   results in exact private colors, which may require a private
   colormap.  */
extern int GpDumpColors(Engine *engine, int colorMode);

/* ------------------------------------------------------------------------ */
/* Error handlers */

/* The Xlib functions invoked by X engines will call exit unless you
   set an alternative error recovery routine.  */
extern int GpSetXHandler(void (*ErrHandler)(char *errMsg));

/* ------------------------------------------------------------------------ */
/* Memory management */

/* GdFree, if non-0, will be called to free objects
   marked with the noCopy flag.  */
extern void (*GdFree)(void *);

#endif
