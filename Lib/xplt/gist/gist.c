/*
 * GIST.C
 *
 * $Id$
 *
 * Implement non-device specific portion of GIST C interface
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "gist.h"
#include "engine.h"
#include "clip.h"

/* Math functions (for determining lengths) */
extern double fabs(double);
extern double sqrt(double);
/* half 1.4142 13562 37309 50488, precise value is not critical here */
#define SQRT_HALF 0.7071067811865

/* Memory manager functions (for scratch space) */
extern void *malloc(long);  /* really size_t which is unsigned */
extern void free(void *);
extern char *strcpy(char *, const char *);

void *(*GmMalloc)(long)= &malloc;
void (*GmFree)(void *)= &free;
void (*GdFree)(void *)= 0;

/* Font for character occasional markers on curves */
#define MARKER_FONT T_HELVETICA

/* ------------------------------------------------------------------------ */

GaAttributes gistA= {
  { FG_COLOR, L_SOLID, 1.0 },                     /* line attributes */
  { FG_COLOR, M_ASTERISK, 1.0 },                /* marker attributes */
  { FG_COLOR, F_SOLID, 0, 0.01, 0.01, 0.0, 0.0 }, /* fill attributes */
  { FG_COLOR, 0, 0.0156,
      TX_RIGHT, TH_NORMAL, TV_NORMAL, 0 },        /* text attributes */
  { 0, 0, 0, 0.16, 0.14, 0,
    0.13, 0.11375, 1.0, 1.0 },          /* decorated line attributes */
  { 0, 0.125 },                                 /* vector attributes */
  { FG_COLOR, L_NONE, 1.0 },                      /* edge attributes */
  };

GpTransform gistT= {   /* default mapping */
  {0.0, 1.0, 0.0, 1.0},
  {0.0, 1.0, 0.0, 1.0}
  };

int gistClip= 0;         /* default clip flag */

GpBox gPortrait=  { 0.0, 0.798584, 0.0, 1.033461 };  /* 8.5 x 11 inches */
GpBox gLandscape= { 0.0, 1.033461, 0.0, 0.798584 };  /* 11 x 8.5 inches */

char gistError[128]= ""; /* most recent error message */

static void InitializeClip(void);
static void MMError(void);
static int GetScratch(long n);
static char SwapTextMark(void);
static void SwapMarkText(void);
static void SwapNormMap(GpReal *scalx, GpReal *offx,
			GpReal *scaly, GpReal *offy);
static void SwapMapNorm(void);
static void ClipArrow(GpReal x[3], GpReal y[3]);
static GpReal TrueNorm(GpReal dx, GpReal dy);
static GpReal OctagNorm(GpReal dx, GpReal dy);
static void DecorateLines(long n, const GpReal *px, const GpReal *py,
			  int closed);
static int MeshRowF(long iMax, long ijMax, int *ireg, int region,
		    long *ii, long *k);
static int MeshColF(long iMax, long ijMax, int *ireg, int region,
		    GpReal *x, GpReal *y, long *jj, long *kk);
static int MeshRowR(long iMax, long ijMax, int *ireg, int region,
		    long *ii, long *k);
static int MeshColR(long iMax, long ijMax, int *ireg, int region,
		    GpReal *x, GpReal *y, long *jj, long *kk);
static int MeshRowB(long iMax, long ijMax, int *ireg, int region,
		    long *ii, long *k);
static int MeshColB(long iMax, long ijMax, int *ireg, int region,
		    GpReal *x, GpReal *y, long *jj, long *kk);
static int *NewReg(long iMax, long ijMax);
static void FreeTmpReg(void);
static int DoSaddle(long zone, long step, long ij, long inc);
static void DoSmoothing(long *n, const GpReal **px, const GpReal **py,
			int closed, int smooth, GpReal scalx, GpReal offx,
			GpReal scaly, GpReal offy);
static int SmoothLines(long n, const GpReal *px, const GpReal *py,
		       int closed, int smooth, int clip);

/* ------------------------------------------------------------------------ */

/* GaLines communicates privately with GpLines via these flags.  */
static int gpCloseNext= 0;   /* (used by GaVectors also) */
static int gpSmoothNext= 0;
static int gpClipDone= 0;
static int gpClipInit= 0;    /* (used by GpFill also) */

static void InitializeClip(void)
{
  int already= gpClipInit;
  gpClipInit= 0;
  if (!already && gistClip) {
    ClipSetup(gistT.window.xmin, gistT.window.xmax,
	      gistT.window.ymin, gistT.window.ymax);
  }
}

int GpLines(long n, const GpReal *px, const GpReal *py)
{
  int value= 0;
  Engine *engine;
  int closed= gpCloseNext;
  int smooth= gpSmoothNext;
  int clip= gistClip && !gpClipDone;
  gpCloseNext= gpSmoothNext=  gpClipDone= 0;

  if (smooth) return SmoothLines(n, px, py, closed, smooth, clip);

  if (clip) InitializeClip();
  else gpClipInit= 0;

  if (!clip || ClipBegin(px, py, n, closed)) {
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      if (!engine->inhibit)
	value|= engine->DrawLines(engine, n, px, py, closed, smooth);
  } else while ((n=ClipMore())) {
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      if (!engine->inhibit)
	value|= engine->DrawLines(engine, n, xClip, yClip, 0, smooth);
  }

  return value;
}

int GpMarkers(long n, const GpReal *px, const GpReal *py)
{
  int value= 0;
  Engine *engine;

  if (gistClip) {
    InitializeClip();
    n= ClipPoints(px, py, n);
    px= xClip;
    py= yClip;
  }
  gpClipInit= 0;

  for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine)) {
    if (!engine->inhibit) {
      if (gistA.m.type<=32) value|= engine->DrawMarkers(engine, n, px, py);
      else value|= GpPseudoMark(engine, n, px, py);
    }
  }

  return value;
}

int GpText(GpReal x0, GpReal y0, const char *text)
{
  int value= 0;
  Engine *engine;

  for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
    if (!engine->inhibit)
      value|= engine->DrawText(engine, x0, y0, text);

  return value;
}

int GpFill(long n, const GpReal *px, const GpReal *py)
{
  int value= 0;
  Engine *engine;

  if (gistClip) {
    InitializeClip();
    n= ClipFilled(px, py, n);
    px= xClip;
    py= yClip;
  }
  gpClipInit= 0;
  if (n<2) return 0;

  for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
    if (!engine->inhibit)
      value|= engine->DrawFill(engine, n, px, py);

  return value;
}

int GpCells(GpReal px, GpReal py, GpReal qx, GpReal qy,
	    long width, long height, long nColumns,
	    const GpColor *colors)
{
  int value= 0;
  Engine *engine;

  for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
    if (!engine->inhibit)
      value|= engine->DrawCells(engine, px, py, qx, qy,
				width, height, nColumns, colors);

  return value;
}

int GpDisjoint(long n, const GpReal *px, const GpReal *py,
	       const GpReal *qx, const GpReal *qy)
{
  int value= 0;
  Engine *engine;

  if (gistClip) {
    InitializeClip();
    n= ClipDisjoint(px, py, qx, qy, n);
    px= xClip;
    py= yClip;
    qx= xClip1;
    qy= yClip1;
  }
  gpClipInit= 0;

  for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
    if (!engine->inhibit)
      value|= engine->DrawDisjoint(engine, n, px, py, qx, qy);

  return value;
}

/* ------------------------------------------------------------------------ */

/* Scratch space requirements--

     GaMesh - needs jMax GpReal points to gather columns
     GaContourInit - needs #edges GpReal points for contours
                     needs ijMax short array for edge markers
		     may need ijMax short array for triangulation
     GaTicks - needs a bunch
     GpLines - needs 3*n real points for smoothing (can't overlap
               with contour scratch space)
 */

extern int GaGetScratchP(long n);
extern int GaGetScratchS(long n);
extern GpReal *gaxScratch, *gayScratch;
extern short *gasScratch;
static long nScratchP= 0, nScratchS= 0;
GpReal *gaxScratch, *gayScratch;
short *gasScratch;

static GpReal *xScratch, *yScratch;
static long nScratch= 0;

static void MMError(void)
{
  strcpy(gistError, "memory manager failed in gist.c function");
}

int GaGetScratchP(long n)
{
  if (n<=nScratchP) return 0;
  if (nScratchP>0) { GmFree(gaxScratch);  GmFree(gayScratch); }
  gaxScratch= (GpReal *)GmMalloc(sizeof(GpReal)*n);
  gayScratch= (GpReal *)GmMalloc(sizeof(GpReal)*n);
  if (!gaxScratch || !gayScratch) {
    if (gaxScratch) GmFree(gaxScratch);
    if (gayScratch) GmFree(gayScratch);
    nScratchP= 0;
    MMError();
    return 1;
  }
  nScratchP= n;
  return 0;
}

static int GetScratch(long n)
{
  if (n<=nScratch) return 0;
  if (nScratch>0) { GmFree(xScratch);  GmFree(yScratch); }
  xScratch= (GpReal *)GmMalloc(sizeof(GpReal)*n);
  yScratch= (GpReal *)GmMalloc(sizeof(GpReal)*n);
  if (!xScratch || !yScratch) {
    if (xScratch) GmFree(xScratch);
    if (yScratch) GmFree(yScratch);
    nScratch= 0;
    MMError();
    return 1;
  }
  nScratch= n;
  return 0;
}

int GaGetScratchS(long n)
{
  if (n<=nScratchS) return 0;
  if (nScratchS>0) GmFree(gasScratch);
  gasScratch= (short *)GmMalloc(sizeof(short)*n);
  if (!gasScratch) {
    nScratchS= 0;
    MMError();
    return 1;
  }
  nScratchS= n;
  return 0;
}

int GaFreeScratch(void)
{
  if (nScratchP>0) { GmFree(gaxScratch);  GmFree(gayScratch); }
  if (nScratchS>0) GmFree(gasScratch);
  if (nScratch>0) { GmFree(xScratch);   GmFree(yScratch); }
  nScratchP= nScratchS= nScratch= 0;
  return 0;
}

/* ------------------------------------------------------------------------ */

static GpTextAttribs textSave;

static char SwapTextMark(void)
{
  /* Remember current text attributes */
  textSave= gistA.t;

  /* Swap in text attributes appropriate for markers */
  gistA.t.color= gistA.m.color;
  gistA.t.font= MARKER_FONT;
  gistA.t.height= gistA.m.size * DEFAULT_MARKER_SIZE;
  gistA.t.orient= TX_RIGHT;
  gistA.t.alignH= TH_CENTER;
  if (gistA.m.type != M_POINT) gistA.t.alignV= TV_HALF;
  else gistA.t.alignV= TV_BASE;
  gistA.t.opaque= 0;     /* later curves can cross opaque markers anyway */

  if (gistA.m.type>M_CROSS || gistA.m.type==0) return (char)gistA.m.type;
  else if (gistA.m.type==M_POINT) return '.';
  else if (gistA.m.type==M_PLUS) return '+';
  else if (gistA.m.type==M_ASTERISK) return '*';
  else if (gistA.m.type==M_CIRCLE) return 'O';
  else return 'X';
}

static void SwapMarkText(void)
{
  gistA.t= textSave;
}

static GpBox windowSave;

static void SwapNormMap(GpReal *scalx, GpReal *offx,
			GpReal *scaly, GpReal *offy)
{
  windowSave= gistT.window;

  *scalx= (gistT.viewport.xmax-gistT.viewport.xmin)/
    (windowSave.xmax-windowSave.xmin);
  *offx= gistT.viewport.xmin - windowSave.xmin*(*scalx);
  *scaly= (gistT.viewport.ymax-gistT.viewport.ymin)/
    (windowSave.ymax-windowSave.ymin);
  *offy= gistT.viewport.ymin - windowSave.ymin*(*scaly);

  gistT.window= gistT.viewport;
  GpSetTrans(&gistT);          /* GpSetTrans checks for this case... */
}

static void SwapMapNorm(void)
{
  gistT.window= windowSave;
  GpSetTrans(&gistT);
}

static void ClipArrow(GpReal x[3], GpReal y[3])
{
  GpReal xmin= gistT.window.xmin, xmax= gistT.window.xmax;
  GpReal ymin= gistT.window.ymin, ymax= gistT.window.ymax;

  /* NB- (x[1],y[1]) is guaranteed to be inside the window */

  if (x[0]<xmin) {
    y[0]+= (xmin-x[0])*(y[1]-y[0])/(x[1]-x[0]);
    x[0]= xmin;
  } else if (x[0]>xmax) {
    y[0]+= (xmax-x[0])*(y[1]-y[0])/(x[1]-x[0]);
    x[0]= xmax;
  }
  if (y[0]<ymin) {
    x[0]+= (ymin-y[0])*(x[1]-x[0])/(y[1]-y[0]);
    y[0]= ymin;
  } else if (y[0]>ymax) {
    x[0]+= (ymax-y[0])*(x[1]-x[0])/(y[1]-y[0]);
    y[0]= ymax;
  }

  if (x[2]<xmin) {
    y[2]+= (xmin-x[2])*(y[1]-y[2])/(x[1]-x[2]);
    x[2]= xmin;
  } else if (x[2]>xmax) {
    y[2]+= (xmax-x[2])*(y[1]-y[2])/(x[1]-x[2]);
    x[2]= xmax;
  }
  if (y[2]<ymin) {
    x[2]+= (ymin-y[2])*(x[1]-x[2])/(y[1]-y[2]);
    y[2]= ymin;
  } else if (y[2]>ymax) {
    x[2]+= (ymax-y[2])*(x[1]-x[2])/(y[1]-y[2]);
    y[2]= ymax;
  }
}

static GpReal TrueNorm(GpReal dx, GpReal dy)
{
  double x= fabs((double)dx);
  double y= fabs((double)dy);
  if (x>y) {
    y= y/x;
    return (GpReal)(x*sqrt(1.0+y*y));
  } else if (y) {
    x= x/y;
    return (GpReal)(y*sqrt(1.0+x*x));
  } else {
    return (GpReal)0.0;
  }
}

static GpReal OctagNorm(GpReal dx, GpReal dy)
{
  double x= fabs((double)dx);
  double y= fabs((double)dy);
  double z= (x+y)*SQRT_HALF;
  if (x>y) {
    return x>z? x : z;
  } else {
    return y>z? y : z;
  }
}

static void DecorateLines(long n, const GpReal *px, const GpReal *py,
			  int closed)
{
  GpReal dx, dy, x1, y1, x0, y0, x00= *px, y00= *py;
  GpReal len, trueLen, xL, yL, xW, yW, x[3], y[3];
  GpReal markSpacing, raySpacing, markPhase, rayPhase;
  GpReal arrowL, arrowW, frac;
  char markText[2];
  int marks= gistA.dl.marks;
  int rays= gistA.dl.rays;
  int type;

  /* Save transform map, set up for transform here */
  GpReal scalx, offx, scaly, offy;
  SwapNormMap(&scalx, &offx, &scaly, &offy);

  markSpacing= gistA.dl.mSpace;
  markPhase= gistA.dl.mPhase;
  if (marks) {
    markText[0]= SwapTextMark();
    markText[1]= '\0';
  }
  raySpacing= gistA.dl.rSpace;
  rayPhase= gistA.dl.rPhase;
  arrowL= gistA.dl.arrowL * DEFAULT_ARROW_LENGTH;
  arrowW= gistA.dl.arrowW * DEFAULT_ARROW_WIDTH;
  type= gistA.l.type;
  if (rays) {
    gistA.l.type= L_SOLID; /* dashed ray arrows can be confusing */
  }

  x0= *px++;
  y0= *py++;
  while (--n>0 || closed) {
    if (n<=0) {
      closed= 0;
      x1= x00;
      y1= y00;
    } else {
      x1= *px++;
      y1= *py++;
    }
    dx= scalx*(x1-x0);
    dy= scaly*(y1-y0);
    len= OctagNorm(dx, dy);

    if (marks) {
      markPhase+= len;
      while (markPhase>=markSpacing) {
	markPhase-= markSpacing;
	frac= markPhase/len;   /* 0<=frac<1, pt=pt0*frac+pt1*(1-frac) */
	/* Since the point is guaranteed to be inside current window,
	   no floating point clip is required here.  Rays are harder.  */
	GpText(x1*scalx+offx-dx*frac, y1*scaly+offy-dy*frac, markText);
      }
    }

    if (rays) {
      rayPhase+= len;
      if (rayPhase>=raySpacing) {
	trueLen= TrueNorm(dx, dy);
	xL= dx/trueLen;
	yW= -arrowW*xL;
	xL*= -arrowL;
	yL= dy/trueLen;
	xW= arrowW*yL;
	yL*= -arrowL;
	while (rayPhase>=raySpacing) {
	  rayPhase-= raySpacing;
	  frac= rayPhase/len;   /* 0<=frac<1, pt=pt0*frac+pt1*(1-frac) */
	  x[1]= x1*scalx+offx-dx*frac;
	  y[1]= y1*scaly+offy-dy*frac;
	  x[0]= x[1]+xL+xW;
	  y[0]= y[1]+yL+yW;
	  x[2]= x[1]+xL-xW;
	  y[2]= y[1]+yL-yW;
	  /* Guaranteed that point of arrow (x[1],y[1]) is inside window,
	     but points may hang outside.  Clipping is tricky because
	     this routine is called from inside a clipping loop.  */
	  ClipArrow(x,y);
	  gpClipDone= 1;
	  GpLines(3, x, y);
	}
      }
    }
    x0= x1;
    y0= y1;
  }

  SwapMapNorm();
  if (marks) SwapMarkText();
  if (rays) gistA.l.type= type;
}

/* GpPseudoMark is potentially useful as the DrawMarkers method
   for an Engine which has no polymarker primitive (e.g.- an X window).
   It is therefore declared as extern in engine.h, although it is
   defined here so it can share the Swap... functions.  */
int GpPseudoMark(Engine *engine, long n, const GpReal *px, const GpReal *py)
{
  int value= 0;
  char text[2];

  /* Set up marker for GpText, swap in appropriate text attributes */
  text[0]= SwapTextMark();
  text[1]= '\0';

  while (--n>=0) value|= engine->DrawText(engine, *px++, *py++, text);
  engine->marked= 1;

  /* Restore text attributes */
  SwapMarkText();
  return value;
}

int GaLines(long n, const GpReal *px, const GpReal *py)
     /* Like GpLines, but includes GaAttributes fancy line attributes
	as well as the line attributes from GpAttributes.  */
{
  int value= 0;

  /* Redirect to polymarker if no line type */
  if (gistA.l.type==L_NONE) {
    /* if (!gistA.dl.marks) return 0;    makes no sense */
    return GpMarkers(n, px, py);
  }

  /* Use the (potentially faster) GpLines routine to draw
     an undecorated polyline.  */
  if (!gistA.dl.marks && !gistA.dl.rays) {
    gpCloseNext= gistA.dl.closed;
    gpSmoothNext= gistA.dl.smooth;
    return GpLines(n, px, py);
  }

  if (gistClip) InitializeClip();
  gpClipInit= 0;

  /* Note that a decorated line cannot be smooth... */
  if (!gistClip || ClipBegin(px, py, n, gistA.dl.closed)) {
    gpCloseNext= gistA.dl.closed;
    gpClipDone= 1;
    value= GpLines(n, px, py);
    DecorateLines(n, px, py, gistA.dl.closed);
  } else while ((n= ClipMore())) {
    gpClipDone= 1;
    value|= GpLines(n, xClip, yClip);
    DecorateLines(n, xClip, yClip, 0);
  }

  return value;
}

/* ------------------------------------------------------------------------ */

/* ARGSUSED */
static int MeshRowF(long iMax, long ijMax, int *ireg, int region,
		    long *ii, long *k)
{
  long i= *ii;
  while (++i<ijMax)   /* scan till edge exists */
    if (ireg[i] || ireg[i+iMax]) break;
  if (i>=ijMax) return 1;
  *k= i-1;
  while (++i<ijMax)   /* scan till edge does not exist */
    if (!ireg[i] && !ireg[i+iMax]) break;
  *ii= i;
  return 0;
}

/* ARGSUSED */
static int MeshColF(long iMax, long ijMax, int *ireg, int region,
		    GpReal *x, GpReal *y, long *jj, long *kk)
{
  long k, j= *jj;
  while ((j+=iMax)<ijMax)   /* scan till edge exists */
    if (ireg[j] || ireg[j+1]) break;
  if (j>=ijMax) return 1;
  gaxScratch[0]= x[j-iMax];   /* gather 1st segment into scratch */
  gayScratch[0]= y[j-iMax];
  gaxScratch[1]= x[j];
  gayScratch[1]= y[j];
  k= 2;
  while ((j+=iMax)<ijMax) { /* scan till edge does not exist */
    if (!ireg[j] && !ireg[j+1]) break;
    gaxScratch[k]= x[j];      /* gather next segment into scratch */
    gayScratch[k]= y[j];
    k++;
  }
  *jj= j;
  *kk= k;
  return 0;
}

static int MeshRowR(long iMax, long ijMax, int *ireg, int region,
		    long *ii, long *k)
{
  long i= *ii;
  while (++i<ijMax)   /* scan till edge exists */
    if (ireg[i]==region || ireg[i+iMax]==region) break;
  if (i>=ijMax) return 1;
  *k= i-1;
  while (++i<ijMax)   /* scan till edge does not exist */
    if (ireg[i]!=region && ireg[i+iMax]!=region ) break;
  *ii= i;
  return 0;
}

static int MeshColR(long iMax, long ijMax, int *ireg, int region,
		    GpReal *x, GpReal *y, long *jj, long *kk)
{
  long k, j= *jj;
  while ((j+=iMax)<ijMax)   /* scan till edge exists */
    if (ireg[j]==region || ireg[j+1]==region) break;
  if (j>=ijMax) return 1;
  gaxScratch[0]= x[j-iMax];   /* gather 1st segment into scratch */
  gayScratch[0]= y[j-iMax];
  gaxScratch[1]= x[j];
  gayScratch[1]= y[j];
  k= 2;
  while ((j+=iMax)<ijMax) { /* scan till edge does not exist */
    if (ireg[j]!=region && ireg[j+1]!=region) break;
    gaxScratch[k]= x[j];      /* gather next segment into scratch */
    gayScratch[k]= y[j];
    k++;
  }
  *jj= j;
  *kk= k;
  return 0;
}

static int MeshRowB(long iMax, long ijMax, int *ireg, int region,
		    long *ii, long *k)
{
  long i= *ii;
  while (++i<ijMax)   /* scan till edge exists */
    if ((ireg[i]==region) ^ (ireg[i+iMax]==region)) break;
  if (i>=ijMax) return 1;
  *k= i-1;
  while (++i<ijMax)   /* scan till edge does not exist */
    if ((ireg[i]!=region) ^ (ireg[i+iMax]==region)) break;
  *ii= i;
  return 0;
}

static int MeshColB(long iMax, long ijMax, int *ireg, int region,
		    GpReal *x, GpReal *y, long *jj, long *kk)
{
  long k, j= *jj;
  while ((j+=iMax)<ijMax)   /* scan till edge exists */
    if ((ireg[j]==region) ^ (ireg[j+1]==region)) break;
  if (j>=ijMax) return 1;
  gaxScratch[0]= x[j-iMax];   /* gather 1st segment into scratch */
  gayScratch[0]= y[j-iMax];
  gaxScratch[1]= x[j];
  gayScratch[1]= y[j];
  k= 2;
  while ((j+=iMax)<ijMax) { /* scan till edge does not exist */
    if ((ireg[j]!=region) ^ (ireg[j+1]==region)) break;
    gaxScratch[k]= x[j];      /* gather next segment into scratch */
    gayScratch[k]= y[j];
    k++;
  }
  *jj= j;
  *kk= k;
  return 0;
}

static int *tmpReg= 0;

static void FreeTmpReg(void)
{
  int *reg= tmpReg;
  tmpReg= 0;
  GmFree(reg);
}

static int *NewReg(long iMax, long ijMax)
{
  if (tmpReg) FreeTmpReg();
  tmpReg= (int *)GmMalloc(sizeof(int)*(ijMax+iMax+1));
  if (tmpReg) {
    long i, j=0;
    for (i=0 ; i<ijMax+iMax+1 ; i++) {
      if (i<1 || i>=iMax || j<1) tmpReg[i]= 0;
      else tmpReg[i]= 1;
      j++;
      if (j==iMax) j= 0;
    }
  } else {
    MMError();
  }
  return tmpReg;
}

int GaMesh(GaQuadMesh *mesh, int region, int boundary, int inhibit)
{
  int value= 0;
  long iMax= mesh->iMax, jMax= mesh->jMax;
  long ijMax= iMax*jMax;
  GpReal *x= mesh->x, *y= mesh->y;
  int *ireg= mesh->reg;
  long i, j, k;
  int (*MeshRow)(long, long, int*, int, long*, long*);
  int (*MeshCol)(long, long, int*, int, GpReal*, GpReal*, long*, long*);

  /* Load appropriate edge existence scanners */
  if (!boundary) {
    if (region==0) {
      /* draw entire mesh */
      MeshRow= &MeshRowF;
      MeshCol= &MeshColF;
    } else {
      /* draw single region */
      MeshRow= &MeshRowR;
      MeshCol= &MeshColR;
    }
  } else {
    /* draw region boundary */
    MeshRow= &MeshRowB;
    MeshCol= &MeshColB;
  }

  /* Be sure there is enough scratch space to gather a column */
  if (!(inhibit&2) && GaGetScratchP(jMax)) return 1;

  /* Create default region array if none supplied */
  if (!ireg) {
    ireg= NewReg(iMax, ijMax);
    if (!ireg) return 1;
    mesh->reg= ireg;
  }

  /* Draw rows */
  if (!(inhibit&1)) {
    for (i=0 ; i<ijMax ; /* i incremented in MeshRow */) {
      if (MeshRow(iMax, ijMax, ireg, region, &i, &j)) break;
      value|= GpLines(i-j, x+j, y+j);
    }
  }

  /* Draw columns */
  if (!(inhibit&2)) {
    for (i=0 ; i<iMax ; i++) {
      j= i;
      for (;;) {
	if (MeshCol(iMax, ijMax, ireg, region, x, y, &j, &k)) break;
	value|= GpLines(k, gaxScratch, gayScratch);
	if (j>=ijMax) break;
      }
    }
  }

  if (tmpReg) FreeTmpReg();
  return value;
}

/* ------------------------------------------------------------------------ */

#define EXISTS(ij) (region? ireg[ij]==region : ireg[ij])

int GaFillMesh(GaQuadMesh *mesh, int region, const GpColor *colors,
	       long nColumns)
{
  int value= 0;
  long iMax= mesh->iMax;
  long ijMax= iMax*mesh->jMax;
  GpReal *x= mesh->x, *y= mesh->y;
  int *ireg= mesh->reg;
  GpReal qx[4], qy[4];
  long ij, row, col;

  /* Create default region array if none supplied */
  if (!ireg) {
    ireg= NewReg(iMax, ijMax);
    if (!ireg) return 1;
    mesh->reg= ireg;
  }

  InitializeClip();

  /* The only filled area attribute set is the color.  */
  row= col= 0;
  for (ij=iMax+1 ; ij<ijMax ; ij++) {
    if (EXISTS(ij)) {
      qx[0]= x[ij-iMax-1];  qy[0]= y[ij-iMax-1];
      qx[1]= x[ij-iMax  ];  qy[1]= y[ij-iMax  ];
      qx[2]= x[ij       ];  qy[2]= y[ij       ];
      qx[3]= x[ij     -1];  qy[3]= y[ij     -1];
      if (colors) gistA.f.color= colors[row+col];
      else gistA.f.color= BG_COLOR;
      gpClipInit= 1;
      value|= GpFill(4, qx, qy);
    }
    col++;
    if (col==iMax) {
      col= 0;
      row+= nColumns;
    }
  }

  if (tmpReg) FreeTmpReg();
  return value;
}

int GaFillMarker(long n, const GpReal *px, const GpReal *py,
		 GpReal x0, GpReal y0)
{
  int value= 0;
  Engine *engine;
  long i;

  /* Save transform map, set up for transform here */
  GpReal scalx, offx, scaly, offy;
  SwapNormMap(&scalx, &offx, &scaly, &offy);
  x0= x0*scalx + offx;
  y0= y0*scaly + offy;

  /* get scratch space, copy points to scratch, adding specified offsets */
  GaGetScratchP(n);
  for (i=0 ; i<n ; i++) {
    gaxScratch[i]= px[i] + x0;
    gayScratch[i]= py[i] + y0;
  }
  px= gaxScratch;
  py= gayScratch;

  if (gistClip) {
    GpReal xmin= gistT.viewport.xmin;
    GpReal xmax= gistT.viewport.xmax;
    GpReal ymin= gistT.viewport.ymin;
    GpReal ymax= gistT.viewport.ymax;
    if (xmin > xmax) { GpReal tmp= xmin; xmin= xmax; xmax= tmp; }
    if (ymin > ymax) { GpReal tmp= ymin; ymin= ymax; ymax= tmp; }
    ClipSetup(xmin, xmax, ymin, ymax);
    n= ClipFilled(px, py, n);
    px= xClip;
    py= yClip;
  }

  if (n>=2) {
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      if (!engine->inhibit)
	value|= engine->DrawFill(engine, n, px, py);
  }

  SwapMapNorm();
  return value;
}

/* ------------------------------------------------------------------------ */

int GaVectors(GaQuadMesh *mesh, int region,
	      const GpReal *u, const GpReal *v, GpReal scale)
{
  int value= 0;
  long iMax= mesh->iMax;
  long ijMax= iMax*mesh->jMax;
  GpReal *x= mesh->x, *y= mesh->y;
  int *ireg= mesh->reg;
  int hollow= gistA.vect.hollow;
  GpReal aspect= gistA.vect.aspect;
  int etype= gistA.e.type;
  GpReal xc, yc, dx, dy, vx[3], vy[3];
  long ij;
  GpReal scalx, offx, scaly, offy, dxscale, dyscale;

  /* Create default region array if none supplied */
  if (!ireg) {
    ireg= NewReg(iMax, ijMax);
    if (!ireg) return 1;
    mesh->reg= ireg;
  }

  /* Save transform map, set up for transform here */
  SwapNormMap(&scalx, &offx, &scaly, &offy);

  dxscale= 0.3333333333*scale*scalx;
  dyscale= 0.3333333333*scale*scaly;
  aspect*= 3.0;

  if (!hollow) gistA.e.type= L_NONE;

  InitializeClip();

  for (ij=0 ; ij<ijMax ; ij++) {
    if (EXISTS(ij) || EXISTS(ij+1) || EXISTS(ij+1+iMax) || EXISTS(ij+iMax)) {
      xc= scalx*x[ij]+offx;
      yc= scaly*y[ij]+offy;
      dx= dxscale*u[ij];
      dy= dyscale*v[ij];

      /* Dart has centroid at (xc,yc), length 3*(dx,dy), given aspect */
      vx[1]= xc + 2.0*dx;
      vy[1]= yc + 2.0*dy;
      vx[0]= xc - dx + aspect*dy;
      vx[2]= xc - dx - aspect*dy;
      vy[0]= yc - dy - aspect*dx;
      vy[2]= yc - dy + aspect*dx;

      if (hollow) {
	gpCloseNext= gpClipInit= 1;
	value|= GpLines(3, vx, vy);
      } else {
	gpClipInit= 1;
	value|= GpFill(3, vx, vy);
      }
    }
  }

  if (!hollow) gistA.e.type= etype;
  if (tmpReg) FreeTmpReg();
  SwapMapNorm();
  return value;
}

/* ------------------------------------------------------------------------ */

/* To use the contouring routines:

      mesh= ...;          **( iMax, jMax, x, y, reg, triangle )**
      for (i=0 ; i<nlevels ; i++) {        **( loop on levels )**
        gistA.dl= linestyle[i];
        if (GaContourInit(mesh, 0, z, level[i]))
          while (GaContour(&n, &px, &py, &gistA.dl.closed))
	    GaLines(n, px, py);    **( loop on disjoint parts )**
      }
 */

/* GaContour state information */
static GaQuadMesh *mesh;
static int region;
static const GpReal *z;
static GpReal level;
static long ni, nj, nib, njb, ibegin, jbegin;
static short *iedges, *jedges;
static int keepLeft;

int GaContourInit(GaQuadMesh *msh, int regn,
		  const GpReal *zz, GpReal lev)
       /* Find the edges cut by the current contour, remembering
	  z, triangle, and level for the GaContour routine, which
	  actually walks the contour.  The z array represents function
	  values on the mesh msh.  If triangle!=0,
	  it represents the initial value of the triangulation array,
	  which determines the behavior of GaContour in saddle zones.
	    triangle[j][i]= 1, -1, or 0 as the zone bounded by
	      [j-1][i-1], [j-1][i], [j][i], and [j][i-1]
	      has been triangulated from (-1,-1) to (0,0),
	      from (-1,0) to (0,-1), or has not yet been
	      triangulated.
	  The msh->triangle array is updated by GaContour.
	  If a contour passes through an untriangulated saddle zone,
	  GaContour chooses a triangulation and marks the triangle
	  array approriately, so that subsequent calls to GaContour
	  will never produce intersecting contour curves.  */
{
  long iMax= msh->iMax;
  long ijMax= iMax*msh->jMax;
  int *ireg= msh->reg;
  long ij;

  /* Create default region array if none supplied */
  if (!ireg) {
    ireg= NewReg(iMax, ijMax);
    if (!ireg) return 0;
    msh->reg= ireg;
  }

  /* Remember data for GaContour */
  mesh= msh;
  region= regn;
  z= zz;
  level= lev;

  /* Get scratch space to hold edges */
  if (GaGetScratchS(2*ijMax)) return 0;
  iedges= gasScratch;
  jedges= gasScratch+ijMax;

  /* Find all points above contour level */
  for (ij=0 ; ij<ijMax ; ij++) iedges[ij]= zz[ij]>lev;

  /* Find j=const edges which cut level plane */
  nj= njb= 0;
  for (ij=1 ; ij<ijMax ; ij++) {
    if ((iedges[ij]^iedges[ij-1]) && (EXISTS(ij)||EXISTS(ij+iMax))) {
      if ((ireg[ij]==region) ^ (ireg[ij+iMax]==region)) {
	jedges[ij]= 2;  /* contour enters mesh here */
	njb++;
      } else {
	jedges[ij]= 1;  /* interior edge */
	nj++;
      }
    } else {
      jedges[ij]= 0;
    }
  }
  jbegin= 1;

  /* Find i=const edges which cut level plane */
  ni= nib= 0;
  for (ij=ijMax-1 ; ij>=iMax ; ij--) {
    if ((iedges[ij]^iedges[ij-iMax]) && (EXISTS(ij)||EXISTS(ij+1))) {
      if ((ireg[ij]==region) ^ (ireg[ij+1]==region)) {
	iedges[ij]= 2;  /* contour enters mesh here */
	nib++;
      } else {
	iedges[ij]= 1;  /* interior edge */
	ni++;
      }
    } else {
      iedges[ij]= 0;
    }
  }
  ibegin= iMax;

  /* Set keepLeft to a known value (arbitrary but repeatable) */
  keepLeft= 0;

  /* Get scratch space for level curves */
  if (GaGetScratchP(ni+nib+nj+njb+1)) return 0;

  if (tmpReg) FreeTmpReg();
  return ni || nib || nj || njb;
}

static int DoSaddle(long zone, long step, long ij, long inc)
{
  /* The contour level plane intersects all four edges of this zone.
     DoSaddle examines the triangulation map to see if the decision
     has already been made.  If not, DoSaddle decides which way to
     turn on the basis of the keepLeft flag (set if the previous
     turn was to the right), the sets the triangulation map entry
     for this zone accordingly.  This ensures future triangulation
     decisions in this zone will be made consistently, so that
     contours at different levels never cross.  */
  short *triang= mesh->triangle;

  /* Return immediately if triangulation already decided.  */
  if (triang && triang[zone]) return triang[zone]>0? (step>0) : (step<0);

  if (inc==1) {  /* currently on j=const edge */
    if (triang) triang[zone]= keepLeft? -1 : 1;
    return keepLeft? (step<0) : (step>0);

  } else {       /* currently on i=const edge */
    if (triang) triang[zone]= keepLeft? 1 : -1;
    return keepLeft? (step>0) : (step<0);
  }
}

int GaContour(long *cn, GpReal **cx, GpReal **cy, int *closed)
       /* After a call to GaContourInit, GaContour must be called
	  repeatedly to generate the sequence of curves obtained by
	  walking the edges cut by the contour level plane.  GaContour
	  returns 1 until there are no more contours to be plotted,
	  when it returns 0.  GaContour signals an error by returning
	  0, but setting *cn!=0.  The curve coordinates (*cx,*cy)
	  use internal scratch space and the associated storage must
	  not be freed.  (*cx, *cy) are valid only until the next
	  GaContour or GaMesh call.  */
{
  long iMax= mesh->iMax;
  long ijMax= iMax*mesh->jMax;
  GpReal *x= mesh->x, *y= mesh->y;
  int *ireg= mesh->reg;
  long ij, zone, step, inc, n;
  GpReal frac;
  int isClosed;
  short mark;

  /* Find a starting point -- get current zone and edge */
  if (nib>0) {
    for (ij=ibegin ; ij<ijMax ; ij++) if (iedges[ij]==2) break;
    if (ij>=ijMax)  return 0; /* this is a bug */
    iedges[ij]= 0;
    zone= EXISTS(ij)? ij : ij+1;
    step= EXISTS(ij)? -1 : 1;
    if (--nib) ibegin= ij+1;
    else ibegin= iMax;
    inc= iMax;

  } else if (njb>0) {
    for (ij=jbegin ; ij<ijMax ; ij++) if (jedges[ij]==2) break;
    if (ij>=ijMax)  return 0; /* this is a bug */
    jedges[ij]= 0;
    zone= EXISTS(ij)? ij : ij+iMax;
    step= EXISTS(ij)? -iMax : iMax;
    if (--njb) jbegin= ij+1;
    else jbegin= 1;
    inc= 1;

  } else if (ni>0) {
    for (ij=ibegin ; ij<ijMax ; ij++) if (iedges[ij]) break;
    if (ij>=ijMax)  return 0; /* this is a bug */
    iedges[ij]= 0;
    zone= ij+1;    /* or ij, doesn't really matter... */
    step= 1;       /* ...this choice tends to go counterclockwise */
    ni--;
    ibegin= ij+1;
    inc= iMax;

  } else if (nj>0) {
    for (ij=jbegin ; ij<ijMax ; ij++) if (jedges[ij]) break;
    if (ij>=ijMax)  return 0; /* this is a bug */
    jedges[ij]= 0;
    zone= ij;      /* or ij+iMax, doesn't really matter... */
    step= -iMax;   /* ...this choice tends to go counterclockwise */
    nj--;
    jbegin= ij+1;
    inc= 1;

  } else {
    return 0;      /* no more edges, only correct way out */
  }

  /* Salt away first point */
  frac= (z[ij]-level)/(z[ij]-z[ij-inc]);
  gaxScratch[0]= frac*(x[ij-inc]-x[ij]) + x[ij];
  gayScratch[0]= frac*(y[ij-inc]-y[ij]) + y[ij];
  n= 1;

  /* Walk the contour */
  isClosed= 0;
  for(;;) {
    /* Find exit from current zone */
    if (inc==1) { /* this is a j=const edge */
      if (jedges[ij+step] && !iedges[zone]) {
	/* step, inc unchanged */
	ij+= step;
      } else if (iedges[zone] &&
		 (!iedges[zone-1] || DoSaddle(zone, step, ij, inc))) {
	ij= zone;
	if (step>0) keepLeft= 1;  /* just turned right */
	else keepLeft= 0;
	step= 1;
	inc= iMax;
      } else if (iedges[zone-1]) {
	ij= zone-1;
	if (step>0) keepLeft= 0;  /* just turned left */
	else keepLeft= 1;
	step= -1;
	inc= iMax;
      } else {
	isClosed= 1;  /* end of a closed contour */
	break;
      }
    } else {      /* this is a i=const edge */
      if (iedges[ij+step] && !jedges[zone]) {
	/* step, inc unchanged */
	ij+= step;
      } else if (jedges[zone] &&
		 (!jedges[zone-iMax] || DoSaddle(zone, step, ij, inc))) {
	ij= zone;
	if (step>0) keepLeft= 0;  /* just turned left */
	else keepLeft= 1;
	step= iMax;
	inc= 1;
      } else if (jedges[zone-iMax]) {
	ij= zone-iMax;
	if (step>0) keepLeft= 1;  /* just turned right */
	else keepLeft= 0;
	step= -iMax;
	inc= 1;
      } else {
	isClosed= 1;  /* end of a closed contour */
	break;
      }
    }

    /* Salt away current point */
    frac= (z[ij]-level)/(z[ij]-z[ij-inc]);
    gaxScratch[n]= frac*(x[ij-inc]-x[ij]) + x[ij];
    gayScratch[n]= frac*(y[ij-inc]-y[ij]) + y[ij];
    n++;

    /* Step into next zone */
    zone+= step;

    /* Erase edge marker for entry edge */
    if (inc==1) {
      mark= jedges[ij];
      jedges[ij]= 0;
      if (mark==2) { /* end of an open contour */
	njb--;
	if (!njb) jbegin= 1;
        break;
      } else nj--;
    } else {
      mark= iedges[ij];
      iedges[ij]= 0;
      if (mark==2) { /* end of an open contour */
	nib--;
	if (!nib) ibegin= iMax;
	break;
      }
      else ni--;
    }
  }

  *cn= n;
  *cx= gaxScratch;
  *cy= gayScratch;
  *closed= isClosed;

  /* Copy first point for closed curves (as a convenience) */
  if (isClosed) {
    gaxScratch[n]= gaxScratch[0];
    gayScratch[n]= gayScratch[0];
  }

  return 1;
}

/* ------------------------------------------------------------------------ */

static void DoSmoothing(long *n, const GpReal **px, const GpReal **py,
			int closed, int smooth, GpReal scalx, GpReal offx,
			GpReal scaly, GpReal offy)
     /* Each point is "split" into three colinear points:
	The central point of each triad is the original point, and
	the other two lie along the direction bisecting the original turn
	at the center point.  For open curves, the first and last points
	are simply duplicated.  Each group of 4 points is suitable as
	the initial, two control, and final points of a Bezier curve.  */
{
  long nn= *n;
  const GpReal *x= *px,  *y= *py;
  GpReal smoothness;

  /* temporaries required in loop */
  GpReal x0, y0, x1, y1, dx0, dy0, dx1, dy1, dsx, dsy, ds0, ds1;
  long i, j;

  if (GetScratch(3*nn+2)) {
    *n= 0;
    return;
  }

  /* support four graded amounts of smoothing--
     the triad of points is spread out more and more as smooth increases */
  if (smooth==1) smoothness= 0.5*0.25/3.0;
  else if (smooth==2) smoothness= 0.5*0.50/3.0;
  else if (smooth==3) smoothness= 0.5*0.75/3.0;
  else smoothness= 0.5*1.0/3.0;

  /* initialize loop on segments */
  x1= scalx*x[0]+offx;
  y1= scaly*y[0]+offy;
  if (closed) {
    /* Previous segment connects last point to first point */
    x0= scalx*x[nn-1]+offx;
    y0= scaly*y[nn-1]+offy;
    dx1= x1-x0;
    dy1= y1-y0;
    ds1= sqrt(dx1*dx1 + dy1*dy1);
    dx1= ds1!=0.0? dx1/ds1 : 0.0;
    dy1= ds1!=0.0? dy1/ds1 : 0.0;
  } else {
    /* Previous segment is zero */
    dx1= dy1= ds1= 0.0;
  }
  j= 1;

  /* do nn-1 segments- (dx1,dy1) is current segment, (dx0,dy0) previous */
  for (i=1 ; i<nn ; i++) {
    xScratch[j]= x0= x1;
    yScratch[j]= y0= y1;
    x1= scalx*x[i]+offx;
    y1= scaly*y[i]+offy;

    dx0= dx1;
    dx1= x1-x0;
    dy0= dy1;
    dy1= y1-y0;
    ds0= ds1;
    /* Note- clipped, normalized coordinates   */
    /* ==> there is no danger of overflow here */
    ds1= sqrt(dx1*dx1 + dy1*dy1);
    dx1= ds1!=0.0? dx1/ds1 : 0.0;
    dy1= ds1!=0.0? dy1/ds1 : 0.0;
    dsx= smoothness * (dx0 + dx1);
    dsy= smoothness * (dy0 + dy1);
    xScratch[j-1]= x0 - ds0*dsx;
    xScratch[j+1]= x0 + ds1*dsx;
    yScratch[j-1]= y0 - ds0*dsy;
    yScratch[j+1]= y0 + ds1*dsy;

    j+= 3;
  }
  /* now i= n and j= 3*n-2, xScratch[3*n-4] has been set */

  if (closed) {
    /* final segment connects last point to first */
    xScratch[j]= x0= x1;
    yScratch[j]= y0= y1;
    x1= scalx*x[0]+offx;
    y1= scaly*y[0]+offy;

    dx0= dx1;
    dx1= x1-x0;
    dy0= dy1;
    dy1= y1-y0;
    ds0= ds1;
    /* Note- clipped, normalized coordinates   */
    /* ==> there is no danger of overflow here */
    ds1= sqrt(dx1*dx1 + dy1*dy1);
    dx1= ds1!=0.0? dx1/ds1 : 0.0;
    dy1= ds1!=0.0? dy1/ds1 : 0.0;
    dsx= smoothness * (dx0 + dx1);
    dsy= smoothness * (dy0 + dy1);
    xScratch[j-1]= x0 - ds0*dsx;
    xScratch[j+1]= x0 + ds1*dsx;
    yScratch[j-1]= y0 - ds0*dsy;
    yScratch[j+1]= y0 + ds1*dsy;

    /* last control point was computed when i=1, and final knot point
       is first point */
    xScratch[j+2]= xScratch[0];
    yScratch[j+2]= yScratch[0];
    xScratch[j+3]= x1;  /* == xScratch[1] */
    yScratch[j+3]= y1;  /* == yScratch[1] */
    *n= j+3;  /* == 3*n+1 (counts first/last point twice) */

  } else {
    /* final control point and final knot point coincide */
    xScratch[j]= xScratch[j-1]= x1;
    yScratch[j]= yScratch[j-1]= y1;
    *n= j;  /* == 3*n-2 */
  }

  *px= xScratch+1;
  *py= yScratch+1;
}

static int SmoothLines(long n, const GpReal *px, const GpReal *py,
		       int closed, int smooth, int clip)
{
  int value= 0;
  Engine *engine;
  GpReal scalx, offx, scaly, offy;

  if (clip && !gpClipInit) InitializeClip();
  else gpClipInit= 0;

  /* Now that clip has been set up, can change coordiante transform--
     output from DoSmoothing is in normalized coordinates, while input
     is in world coordinates...  */
  SwapNormMap(&scalx, &offx, &scaly, &offy);

  if (!clip || ClipBegin(px, py, n, closed)) {
    DoSmoothing(&n, &px, &py, closed, smooth, scalx, offx, scaly, offy);
    /* Note: if closed, n= 3*n+1, last point same as first */
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      if (!engine->inhibit)
	value|= engine->DrawLines(engine, n, px, py, 0, smooth);
  } else while ((n=ClipMore())) {
    px= xClip;
    py= yClip;
    DoSmoothing(&n, &px, &py, 0, smooth, scalx, offx, scaly, offy);
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      if (!engine->inhibit)
	value|= engine->DrawLines(engine, n, px, py, 0, smooth);
  }

  SwapMapNorm();

  return value;
}

/* ------------------------------------------------------------------------ */
