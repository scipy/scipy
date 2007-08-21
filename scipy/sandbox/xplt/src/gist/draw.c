/*
 * DRAW.C
 *
 * $Id$
 *
 * Implement display list portion of GIST C interface
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "draw.h"
#include "gtext.h"
#include "pstdlib.h"

/* Generating default contour labels requires sprintf function */
#include <stdio.h>

#include <string.h>

extern double log10(double);
#define SAFELOG0 (-999.)
#define SAFELOG(x) ((x)>0? log10(x) : ((x)<0? log10(-(x)) : -999.))

/* In tick.c */
extern GpReal GpNiceUnit(GpReal finest, int *base, int *power);

extern double floor(double);
extern double ceil(double);
#ifndef NO_EXP10
  extern double exp10(double);
#else
# define exp10(x) pow(10.,x)
  extern double pow(double,double);
#endif
#define LOG2 0.301029996

/* ------------------------------------------------------------------------ */

static Drauing *currentDr= 0;  /* Drauing from GdNewDrawing or GdGetDrawing */
static GeSystem *currentSy;    /* System from GdNewSystem or GdSetSystem */
static GdElement *currentEl;   /* Element extracted with GdSetSystem, or
                                  GdGetElement */
static int currentCn;          /* level selected by GdGetContour */

/* Saved state information used by GdSetDrawing */
static Drauing *saveDr= 0;
static GeSystem *saveSy; 
static GdElement *saveEl;
static int saveCn;

GdProperties gistD= {
  0, 0,
  {
  {7.5, 50., 1.2, 1.2, 4, 1, TICK_L|TICK_U|TICK_OUT|LABEL_L,
   0.0, 14.0*ONE_POINT,
   {12.*ONE_POINT, 8.*ONE_POINT, 5.*ONE_POINT, 3.*ONE_POINT, 2.*ONE_POINT},
     {FG_COLOR, L_SOLID, 1.0},
     {FG_COLOR, L_DOT, 1.0},
     {FG_COLOR, T_HELVETICA, 14.*ONE_POINT,
        TX_RIGHT, TH_NORMAL, TV_NORMAL, 1},
     .425, .5-52.*ONE_POINT},
  {7.5, 50., 1.2, 1.2, 3, 1, TICK_L|TICK_U|TICK_OUT|LABEL_L,
   0.0, 14.0*ONE_POINT,
   {12.*ONE_POINT, 8.*ONE_POINT, 5.*ONE_POINT, 3.*ONE_POINT, 2.*ONE_POINT},
     {FG_COLOR, L_SOLID, 1.0},
     {FG_COLOR, L_DOT, 1.0},
     {FG_COLOR, T_HELVETICA, 14.*ONE_POINT,
        TX_RIGHT, TH_NORMAL, TV_NORMAL, 1},
     .25, .5-52.*ONE_POINT},
  0, {FG_COLOR, L_SOLID, 1.0}
  },
  {{ 0.0, 1.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0, 1.0 }},
  D_XMIN | D_XMAX | D_YMIN | D_YMAX,    /* flags */
  { 0.0, 1.0, 0.0, 1.0 },               /* limits */
  0, 0, 0, 0, 0,
  0.0, 0.0, 0,
  0.0, 0.0, 0.0, 0.0, 0, 0,             /* GdCells */
  0, 0,
  0, { 0, 0, 0, 0, 0, 0 }, 0,           /* noCopy, mesh, region */
  0, 0,
  0, 0, 0.0,
  0, 0, 0,  0 };

Drauing *gistDrawList= 0;

static void ClearDrawing(Drauing *drawing);
static void Damage(GeSystem *sys, GdElement *el);
static void SquareAdjust(GpReal *umin, GpReal *umax,
                         GpReal dv, int doMin, int doMax);
static void NiceAdjust(GpReal *umin, GpReal *umax, int isLog, int isMin);
static void EqAdjust(GpReal *umin, GpReal *umax);
static void EmptyAdjust(GpReal *umin, GpReal *umax, int doMin, int doMax);
static void EqualAdjust(GpReal *umin, GpReal *umax, int doMin, int doMax);
extern int Gd_DrawRing(void *elv, int xIsLog, int yIsLog,
                       GeSystem *sys, int t);
static void InitLegends(int contours, GeSystem *systems, GdElement *elements,
                        int size);
static void NextContours(void);
static int NextRing(void);
static int NextLegend(void);
static int BuildLegends(int more, int contours, GeSystem *systems,
                        GdElement *elements, GeLegendBox *lbox);
static int MemoryError(void);
static void *Copy1(const void *orig, long size);
static void *Copy2(void *x1, const void *orig1, const void *orig2, long size);
extern void Gd_ScanZ(long n, const GpReal *z, GpReal *zmin, GpReal *zmax);
static void ScanXY(long n, const GpReal *x, const GpReal *y, GpBox *extrema);

extern void Gd_NextMeshBlock(long *ii, long *jj, long len, long iMax,
                             int *reg, int region);
extern void Gd_MeshXYGet(void *vMeshEl);
static int AutoMarker(GaLineAttribs *dl, int number);
extern int Gd_MakeContours(GeContours *con);
static void GuessBox(GpBox *box, GpBox *viewport, GaTickStyle *ticks);
static GdElement *NextConCurve(GdElement *el);
static int GeFindIndex(int id, GeSystem *sys);

extern void Gd_KillRing(void *elv);
extern void Gd_KillMeshXY(void *vMeshEl);

static void (*DisjointKill)(void *el);
static void (*FilledKill)(void *el);
static void (*VectorsKill)(void *el);
static void (*ContoursKill)(void *el);
static void (*SystemKill)(void *el);
static int (*LinesGet)(void *el);
static int (*ContoursGet)(void *el);
extern void Gd_LinesSubSet(void *el);
static int (*SystemDraw)(void *el, int xIsLog, int yIsLog);

/* ------------------------------------------------------------------------ */
/* Set virtual function tables */

extern GdOpTable *GetDrawingOpTables(void);  /* in draw0.c */
static GdOpTable *opTables= 0;

/* ------------------------------------------------------------------------ */
/* Constructor and destructor for Drauing declared in gist.h */

Drauing *GdNewDrawing(char *gsFile)
{
  Drauing *drawing= p_malloc(sizeof(Drauing));
  if (!drawing) return 0;
  if (!opTables) {
    opTables= GetDrawingOpTables();
    DisjointKill= opTables[E_DISJOINT].Kill;
    FilledKill= opTables[E_FILLED].Kill;
    VectorsKill= opTables[E_VECTORS].Kill;
    ContoursKill= opTables[E_CONTOURS].Kill;
    SystemKill= opTables[E_SYSTEM].Kill;
    LinesGet= opTables[E_LINES].GetProps;
    ContoursGet= opTables[E_CONTOURS].GetProps;
    SystemDraw= opTables[E_SYSTEM].Draw;
  }

  drawing->next= gistDrawList;
  gistDrawList= drawing;
  drawing->cleared=
    drawing->nSystems= drawing->nElements= 0;
  drawing->systems= 0;
  drawing->elements= 0;
  drawing->damaged= 0;
  drawing->damage.xmin= drawing->damage.xmax=
    drawing->damage.ymin= drawing->damage.ymax= 0.0;
  drawing->landscape= 0;
  drawing->legends[0].nlines= drawing->legends[1].nlines= 0;

  GdSetDrawing(drawing);

  if (GdReadStyle(drawing, gsFile)) {
    GdSetDrawing(0);
    GdKillDrawing(drawing);
    return 0;
  }

  return drawing;
}

int GdLandscape(int landscape)
{
  if (!currentDr) return 1;
  if (landscape) landscape= 1;
  if (currentDr->landscape!=landscape) {
    currentDr->landscape= landscape;
    GdDetach(currentDr, 0);
  }
  return 0;
}

void GdKillDrawing(Drauing *drawing)
{
  if (!drawing) {
    drawing= currentDr;
    if (!drawing) return;
  }

  ClearDrawing(drawing);
  Gd_KillRing(drawing->systems);

  if (drawing==gistDrawList) gistDrawList= drawing->next;
  else {
    Drauing *draw= gistDrawList;
    while (draw->next!=drawing) draw= draw->next;
    draw->next= drawing->next;
  }

  if (drawing==currentDr) currentDr= 0;

  p_free(drawing);
}

extern void GdKillSystems(void);

void GdKillSystems(void)
{
  if (!currentDr) return;
  ClearDrawing(currentDr);
  Gd_KillRing(currentDr->systems);
  currentDr->systems= 0;
  currentDr->nSystems= 0;
}

int GdSetDrawing(Drauing *drawing)
{
  int nMax, sysIndex, i;
  GeSystem *sys;

  if (!drawing) {  /* swap current and saved state info */
    Drauing *tmpDr= currentDr;
    GeSystem *tmpSy= currentSy;
    GdElement *tmpEl= currentEl;
    int tmpCn= currentCn;
    currentDr= saveDr;   saveDr= tmpDr;
    currentSy= saveSy;   saveSy= tmpSy;
    currentEl= saveEl;   saveEl= tmpEl;
    currentCn= saveCn;   saveCn= tmpCn;
    return 0;
  }

  saveDr= currentDr;
  saveSy= currentSy;
  saveEl= currentEl;
  saveCn= currentCn;

  currentDr= drawing;

  /* Make a reasonable guess at current system and element */
  nMax= drawing->elements? drawing->elements->prev->number : -1;
  sysIndex= drawing->nSystems? 1 : 0;
  i= 0;
  if ((sys= drawing->systems)) do {
    i++;
    if (sys->el.number>nMax) { nMax= sys->el.number;  sysIndex= i; }
    sys= (GeSystem *)sys->el.next;
  } while (sys!=drawing->systems);
  GdSetSystem(sysIndex);
  if (nMax<0) {
    if (sysIndex<1) currentSy= 0;
    currentEl= 0;
  } else {
    GdElement *el= currentSy? currentSy->elements : drawing->elements;
    if (el) {
      currentEl= el->prev;
      currentEl->ops->GetProps(currentEl);
    } else {
      currentEl= 0;
    }
  }
  currentCn= -1;

  return 0;
}

int GdClear(Drauing *drawing)
{
  if (!drawing) drawing= currentDr;
  if (!drawing) return 1;
  drawing->cleared= 1;
  return 0;
}

static void ClearDrawing(Drauing *drawing)
{
  GeSystem *sys, *sys0= drawing->systems;
  int nSystems= 0;
  if ((sys= sys0)) do {
    Gd_KillRing(sys->elements);
    sys->elements= 0;
    sys->rescan= 0;
    sys->unscanned= -1;
    sys->el.number= -1;
    sys= (GeSystem *)sys->el.next;
    nSystems++;
  } while (sys!=sys0);
  Gd_KillRing(drawing->elements);
  drawing->elements= 0;
  drawing->nElements= 0;
  drawing->nSystems= nSystems;
  drawing->cleared= 2;

  if (drawing==currentDr) {
    currentSy= drawing->systems; /* as after GdSetDrawing */
    currentEl= 0;
    currentCn= -1;
  }

  /* Must detatch drawing from all engines, since even inactive ones
     need to know that the drawing has been erased.  */
  GdDetach(drawing, (Engine *)0);
}

/* ------------------------------------------------------------------------ */

/* The GIST display list is called a "drawing".  Any number of drawings
   may exist, but only one may be active at a time.

   The entire drawing is rendered on all active engines by calling
   GdDraw(0).  The drawing sequence is preceded by GpClear(0, CONDITIONALLY).

   GdDraw(-1) is like GdDraw(0) except the drawing is damaged to force
   all data to be rescanned as well.

   GdDraw(1) draws only the changes since the previous GdDraw.  Changes can
   be either destructive or non-destructive:

     If the engine->damage box is set, then some operation has
     damaged that part of the drawing since the prior GdDraw
     (such as changing a line style or plot limits).  The damaged
     section is cleared, then the display list is traversed, redrawing
     all elements which may intersect the damaged box, clipped to
     the damaged box.

     Then, any elements which were added since the prior GdDraw
     (and which caused no damage like a change in limits) are
     drawn as non-destructive updates.

 */

static void Damage(GeSystem *sys, GdElement *el)
{
  GpBox *box, adjustBox;
  if (!currentDr) return;
  if (!el) {
    if (!sys) return;
    /* If no element, damage the entire coordinate system, ticks and all */
    box= &sys->el.box;
  } else if (sys) {
    /* If element is in a coordinate system, damage the whole viewport */
    box= &sys->trans.viewport;
  } else {
    /* Elements not in a coordinate system already have NDC box--
       but may need to adjust it to allow for projecting junk.  */
    el->ops->Margin(el, &adjustBox);
    adjustBox.xmin+= el->box.xmin;
    adjustBox.xmax+= el->box.xmax;
    adjustBox.ymin+= el->box.ymin;
    adjustBox.ymax+= el->box.ymax;
    box= &adjustBox;
  }
  if (currentDr->damaged) {
    GpSwallow(&currentDr->damage, box);
  } else {
    currentDr->damage= *box;
    currentDr->damaged= 1;
  }
}

static void SquareAdjust(GpReal *umin, GpReal *umax,
                         GpReal dv, int doMin, int doMax)
{
  if (doMin) {
    if (doMax) umax[0]= 0.5*(umin[0]+umax[0]+dv);
    umin[0]= umax[0]-dv;
  } else if (doMax) {
    umax[0]= umin[0]+dv;
  }
}

static void NiceAdjust(GpReal *umin, GpReal *umax, int isLog, int isMin)
{
  GpReal un= *umin,   ux= *umax;
  GpReal unit,   du= ux-un;
  int base, power, reverted= 0;
  if (isLog) {
    if (du <= LOG2) {
      /* revert to linear scale */
      un= exp10(un);
      ux= exp10(ux);
      du= ux-un;
      reverted= 1;
    } else if (du>6.0 && isMin) {
      un= ux-6.0;
      du= 6.0;
    }
  }
  unit= GpNiceUnit(du/3.0, &base, &power);
  if (!isLog || reverted || unit>0.75) {
    un= unit*floor(un/unit);
    ux= unit*ceil(ux/unit);
    if (reverted) {
      un= log10(un);
      ux= log10(ux);
    }
  } else {
    /* subdecade log scale (sigh), use 2, 5, or 10 */
    GpReal dn= floor(un+0.0001),   dx= ceil(ux-0.0001);
    if (un<dn+(LOG2-0.0001)) un= dn;
    else if (un<dn+(0.9999-LOG2)) un= dn+LOG2;
    else un= dn+(1.0-LOG2);
    if (ux>dx-(LOG2+0.0001)) ux= dx;
    else if (ux>dx-(1.0001-LOG2)) ux= dx-LOG2;
    else ux= ux-(1.0-LOG2);
  }
  *umin= un;
  *umax= ux;
}

static void EqAdjust(GpReal *umin, GpReal *umax)
{
  GpReal nudge= *umin>0.0? 0.001*(*umin) : -0.001*(*umin);
  if (nudge==0.0) nudge= 1.0e-6;
  *umin-= nudge;
  *umax+= nudge;
}

static void EmptyAdjust(GpReal *umin, GpReal *umax, int doMin, int doMax)
{
  if (doMin) {
    if (doMax) { *umin= -1.0e-6; *umax= +1.0e-6; }
    else if (*umax>0.0) *umin= 0.999*(*umax);
    else if (*umax<0.0) *umin= 1.001*(*umax);
    else *umin= -1.0e-6;
  } else if (doMax) {
    if (*umin>0.0) *umax= 1.001*(*umin);
    else if (*umin<0.0) *umax= 0.999*(*umin);
    else *umax= +1.0e-6;
  } else if ((*umin)==(*umax)) {
    EqAdjust(umin, umax);
  }
}

static void EqualAdjust(GpReal *umin, GpReal *umax, int doMin, int doMax)
{
  if (doMin && doMax) EqAdjust(umin, umax);
  else EmptyAdjust(umin, umax, doMin, doMax);
}

int GdScan(GeSystem *sys)
{
  int flags= sys->flags;
  GpBox limits, tmp, *w= &sys->trans.window;
  GpReal xmin= w->xmin, xmax= w->xmax, ymin= w->ymin, ymax= w->ymax;
  int swapx, swapy;
  GdElement *el, *el0= sys->elements;
  int begin, damaged, first;

  /* Handle case of no curves (if, e.g., all elements removed) */
  if (!el0) {
    EmptyAdjust(&w->xmin, &w->xmax, flags&D_XMIN, flags&D_XMAX);
    EmptyAdjust(&w->ymin, &w->ymax, flags&D_YMIN, flags&D_YMAX);
    return 0;
  }

  /* Assure that limits are ordered even if window is not */
  swapx= (xmin > xmax) && !(flags&(D_XMIN|D_XMAX));
  swapy= (ymin > ymax) && !(flags&(D_YMIN|D_YMAX));
  if (!swapx) { limits.xmin= xmin;  limits.xmax= xmax; }
  else { limits.xmin= xmax;  limits.xmax= xmin; }
  if (!swapy) { limits.ymin= ymin;  limits.ymax= ymax; }
  else { limits.ymin= ymax;  limits.ymax= ymin; }
  tmp= limits;

  el= el0;
  begin= sys->rescan? -1 : sys->unscanned;

  /* Scan limits for each element */
  damaged= 0;
  first= 1;
  do {
    if (!el->hidden) {
      if (el->number>=begin) {
        /* Scan ensures log values present, sets box, scans xy values */
        if (el->ops->Scan(el, flags, &tmp)) return 1; /* mem failure */
        if (first) {
          /* first non-hidden element gives first cut at limits */
          limits= tmp;
          damaged= 1;
        } else {
          /* subsequent elements may cause limits to be adjusted */
          if (tmp.xmin<=tmp.xmax) {
            if (tmp.xmin<limits.xmin) limits.xmin= tmp.xmin;
            if (tmp.xmax>limits.xmax) limits.xmax= tmp.xmax;
          }
          if (tmp.ymin<=tmp.ymax) {
            if (tmp.ymin<limits.ymin) limits.ymin= tmp.ymin;
            if (tmp.ymax>limits.ymax) limits.ymax= tmp.ymax;
          }
        }
      }
      first= 0;
    }
    el= el->next;
  } while (el!=el0);

  /* #1- adjust if min==max */
  if (limits.xmin==limits.xmax)
    EqualAdjust(&limits.xmin, &limits.xmax, flags&D_XMIN, flags&D_XMAX);
  if (limits.ymin==limits.ymax)
    EqualAdjust(&limits.ymin, &limits.ymax, flags&D_XMIN, flags&D_XMAX);

  /* #2- adjust if log axis and minimum was SAFELOG(0) */
  if ((flags & D_LOGX) && (flags & D_XMIN) && limits.xmin==SAFELOG0
      && limits.xmax>SAFELOG0+10.0) limits.xmin= limits.xmax-10.0;
  if ((flags & D_LOGY) && (flags & D_YMIN) && limits.ymin==SAFELOG0
      && limits.ymax>SAFELOG0+10.0) limits.ymin= limits.ymax-10.0;

  /* #3- adjust if square limits specified and not semi-logarithmic */
  if ((flags & D_SQUARE) &&
      !(((flags&D_LOGX)!=0) ^ ((flags&D_LOGY)!=0))) {
    /* (Square axes don't make sense for semi-log scales) */
    GpReal dx= limits.xmax-limits.xmin;
    GpReal dy= limits.ymax-limits.ymin;
    GpReal dydx= (sys->trans.viewport.ymax-sys->trans.viewport.ymin)/
                 (sys->trans.viewport.xmax-sys->trans.viewport.xmin);
    /* Adjust y if either (1) dx>dy, or (2) x limits are both fixed
       (NB- SquareAdjust is a noop if both limits fixed) */
    if ((dx*dydx>dy && (flags&(D_YMIN|D_YMAX))) ||
        !(flags&(D_XMIN|D_XMAX)))
      SquareAdjust(&limits.ymin, &limits.ymax, dx*dydx,
                   flags&D_YMIN, flags&D_YMAX);
    else /* adjust x */
      SquareAdjust(&limits.xmin, &limits.xmax, dy/dydx,
                   flags&D_XMIN, flags&D_XMAX);
  }

  /* #4- adjust if nice limits specified */
  if (flags & D_NICE) {
    NiceAdjust(&limits.xmin, &limits.xmax, flags&D_LOGX, flags&D_XMIN);
    NiceAdjust(&limits.ymin, &limits.ymax, flags&D_LOGY, flags&D_YMIN);
  }

  if (swapx) {
    GpReal tmp= limits.xmin;  limits.xmin= limits.xmax;  limits.xmax= tmp;
  }
  if (swapy) {
    GpReal tmp= limits.ymin;  limits.ymin= limits.ymax;  limits.ymax= tmp;
  }
  if (damaged || limits.xmin!=xmin || limits.xmax!=xmax ||
      limits.ymin!=ymin || limits.ymax!=ymax)
    Damage(sys, (GdElement *)0);
  w->xmin= limits.xmin;
  w->xmax= limits.xmax;
  w->ymin= limits.ymin;
  w->ymax= limits.ymax;

  sys->rescan= 0;
  sys->unscanned= -1;

  return 0;
}

int Gd_DrawRing(void *elv, int xIsLog, int yIsLog, GeSystem *sys, int t)
{
  GdElement *el0, *el= elv;
  GpBox adjustBox, *box;
  int value= 0, drawIt= t;
  if ((el0= el)) {
    do {
      if (!t) {
        if (!sys) {
          el->ops->Margin(el, &adjustBox);
          adjustBox.xmin+= el->box.xmin;
          adjustBox.xmax+= el->box.xmax;
          adjustBox.ymin+= el->box.ymin;
          adjustBox.ymax+= el->box.ymax;
          box= &adjustBox;
        } else {
          box= &sys->trans.viewport;
        }
        drawIt= GdBeginEl(box, el->number);
      }
      if (drawIt) value|= el->ops->Draw(el, xIsLog, yIsLog);
      el= el->next;
    } while (el!=el0);
  }
  return value;
}

static GpTransform unitTrans= { {0., 2., 0., 2.}, {0., 2., 0., 2.} };

int GdDraw(int changesOnly)
{
  int value= 0;
  GpBox *damage;
  int systemCounter;
  int rescan= 0;

  if (!currentDr) return 1;

  if (changesOnly==-1) {
    rescan= 1;
    changesOnly= 0;
  }

  /* Take care of conditional clear */
  if (currentDr->cleared==1) {
    if (changesOnly) return 0;
    else ClearDrawing(currentDr);
  }
  if (!changesOnly || currentDr->cleared) {
    GpClear(0, CONDITIONALLY);
    currentDr->cleared= 0;
  }

  /* Check if any coordinate systems need to be rescanned */
  if (currentDr->systems) {
    int changed;
    GeSystem *sys, *sys0;
    sys= sys0= currentDr->systems;
    do {
      if (rescan) sys->rescan= 1;
      changed= (sys->rescan || sys->unscanned>=0);
      if (changed) changesOnly= 0;
      if (changed && GdScan(sys)) return 1;  /* memory manager failure */
      sys= (GeSystem *)sys->el.next;
    } while (sys!=sys0);
  }

  /* Give engines a chance to prepare for a drawing */
  if (currentDr->damaged) {
    damage= &currentDr->damage;
    currentDr->damaged= 0;
  } else {
    damage= 0;
  }
  /* GdBeginDr returns 1 if any active engine has been cleared or
     partially cleared.  */
  if (!GdBeginDr(currentDr, damage, currentDr->landscape) && changesOnly)
    return 0;

  /* Do coordinate systems */
  if (currentDr->systems) {
    GeSystem *sys, *sys0;
    sys= sys0= currentDr->systems;
    systemCounter= 0;
    do {
      value|= SystemDraw(sys, systemCounter, 0);
      systemCounter++;
      sys= (GeSystem *)sys->el.next;
    } while (sys!=sys0);
  }

  /* Do elements outside of coordinate systems */
  GpSetTrans(&unitTrans);
  gistClip= 0;
  value|= Gd_DrawRing(currentDr->elements, 0, 0, (GeSystem *)0, 0);

  /* Give engines a chance to clean up after a drawing */
  GdEndDr();

  return value;
}

/* ------------------------------------------------------------------------ */
/* Legend routines */

int GdLegendBox(int which, GpReal x, GpReal y, GpReal dx, GpReal dy,
                const GpTextAttribs *t, int nchars, int nlines, int nwrap)
{
  GeLegendBox *lbox;
  if (!currentDr || nchars<0) return 1;
  lbox= currentDr->legends;
  if (which) lbox++;
  lbox->x= x;     lbox->y= y;
  lbox->dx= dx;   lbox->dy= dy;
  lbox->textStyle= *t;
  lbox->nchars= nchars;
  lbox->nlines= nlines;
  lbox->nwrap= nwrap;
  return 0;
}

static char *legendText= 0;
static long lenLegends, maxLegends= 0;

static int nRemaining, curWrap;
static char *curLegend;
static int curMarker= 0;

static int doingContours, levelCurve, nLevels;
static GdElement *curElement, *cur0Element, *drElements, *curCon, *cur0Con;
static GeSystem *curSystem, *cur0System;
static GpReal *levelValue;
static GeLines **levelGroup;
static char levelLegend[32];

static void InitLegends(int contours, GeSystem *systems, GdElement *elements,
                        int size)
{
  doingContours= levelCurve= contours;
  if (doingContours) curCon= 0;
  curElement= 0;
  curSystem= cur0System= systems;
  drElements= elements;
  curLegend= 0;
  curMarker= 0;
  nRemaining= 0;

  if (size>maxLegends) {
    if (legendText) p_free(legendText);
    legendText= p_malloc((long)size);
  }
}

static void NextContours(void)
{
  if (!levelCurve) {
    /* Set up for the ring of level curves */
    GeContours *con= (GeContours *)curCon;
    nLevels= con->nLevels;
    levelValue= con->levels;
    levelGroup= con->groups;
    levelCurve= 1;
    if (levelGroup) {
      while (nLevels && !levelGroup[0]) {
        levelValue++;
        levelGroup++;
        nLevels--;
      }
    } else {
      nLevels= 0;
    }
    if (nLevels>0) curElement= (GdElement *)levelGroup[0];
    else curElement= 0;
    return;
  }

  levelCurve= 0;
  curElement= 0;
  if (curCon) {
    curCon= curCon->next;
    if (curCon==cur0Con) curCon= 0;
  }
  for (;;) {
    if (curCon) {
      do {
        if (curCon->ops->type==E_CONTOURS && !curCon->hidden) {
          /* Set up for contour element itself-- terminates immediately */
          curElement= curCon;
          cur0Element= curElement->next;
          return;
        }
        curCon= curCon->next;
      } while (curCon!=cur0Con);
    }

    if (curSystem) {
      curCon= cur0Con= curSystem->elements;
      curSystem= (GeSystem *)curSystem->el.next;
      if (curSystem==cur0System) curSystem= 0;
    } else if (drElements) {
      curCon= cur0Con= drElements;
      drElements= 0;
    } else {
      break;
    }
  }
}

static int NextRing(void)
{
  if (doingContours) {
    NextContours();
    if (!curElement) return 0;
  } else if (curSystem) {
    curElement= cur0Element= curSystem->elements;
    curSystem= (GeSystem *)curSystem->el.next;
    if (curSystem==cur0System) curSystem= 0;
  } else if (drElements) {
    curElement= cur0Element= drElements;
    drElements= 0;
  } else {
    return 0;
  }
  return 1;
}

static int specialMarks[5]= { '.', '+', '*', 'o', 'x' };

static int NextLegend(void)
{
  curLegend= 0;
  curMarker= 0;
  do {
    while (curElement) {
      if (!curElement->hidden) {
        int type= curElement->ops->type;
        if (curElement->legend) curLegend= curElement->legend;
        else if (levelCurve) {
          /* automatically generate level curve legend if not supplied */
          curLegend= levelLegend;
          sprintf(curLegend, "\001: %.4g", *levelValue);
        }
        if (curLegend) {
          nRemaining= strlen(curLegend);
          curWrap= 0;
          if ((type==E_LINES || type==E_CONTOURS) && curLegend[0]=='\001') {
            /* insert marker into E_LINES legend if so directed */
            curMarker= type==E_LINES? ((GeLines *)curElement)->m.type :
            ((GeContours *)curElement)->m.type;
            if (curMarker>=1 && curMarker<=5)
              curMarker= specialMarks[curMarker-1];
            else if (curMarker<' ' || curMarker>='\177')
              curMarker= ' ';
          }
        }
      }
      if (levelCurve) {
        do {
          levelValue++;
          levelGroup++;
          nLevels--;
        } while (nLevels && !levelGroup[0]);
        if (nLevels>0) curElement= (GdElement *)levelGroup[0];
        else curElement= 0;
      } else {
        curElement= curElement->next;
        if (curElement==cur0Element) curElement= 0;
      }
      if (curLegend) return 1;
    }
  } while (NextRing());
  return 0;
}

static int BuildLegends(int more, int contours, GeSystem *systems,
                        GdElement *elements, GeLegendBox *lbox)
{
  int firstLine= 1;
  int nlines= lbox->nlines;
  int nchars= lbox->nchars;
  int nwrap= lbox->nwrap;
  int nc;

  lenLegends= 0;
  if (!more) {
    if (nlines<=0 || nchars<=0) return 0;
    InitLegends(contours, systems, elements, (nchars+1)*nlines);
    if (!legendText) return 0;
  }

  for ( ; ; nlines--) {
    if (!curLegend && !NextLegend()) { more= 0;   break; }
    if (nlines<=0) { more= !more;   break; }
    if (firstLine) firstLine= 0;
    else legendText[lenLegends++]= '\n';
    nc= nRemaining>nchars? nchars : nRemaining;
    strncpy(legendText+lenLegends, curLegend, nc);
    if (curMarker) {
      legendText[lenLegends]= (char)curMarker;
      curMarker= 0;
    }
    lenLegends+= nc;
    nRemaining-= nc;
    if (nRemaining>0 && curWrap++<nwrap) curLegend+= nc;
    else { curLegend= 0; curMarker= 0; }
  }

  legendText[lenLegends]= '\0';
  return more;
}

int GdDrawLegends(Engine *engine)
{
  GpReal x, y;
  int type, more;
  GeLegendBox *lbox;
  if (!currentDr) return 1;

  if (engine) GpPreempt(engine);

  for (type=0 ; type<2 ; type++) {
    lbox= &currentDr->legends[type];
    x= lbox->x;
    y= lbox->y;
    gistA.t= lbox->textStyle;
    GpSetTrans(&unitTrans);
    gistClip= 0;
    if (lbox->nlines <= 0) continue;
    for (more=0 ; ; ) {
      more= BuildLegends(more, type, currentDr->systems, currentDr->elements,
                         lbox);
      if (!legendText) {
        /* memory error */
        if (engine) GpPreempt(0);
        return 1;
      }
      if (lenLegends>0) GpText(x, y, legendText);
      if (!more || (lbox->dx==0.0 && lbox->dy==0.0)) break;
      x+= lbox->dx;
      y+= lbox->dy;
    }
  }

  if (engine) GpPreempt(0);
  return 0;
}

/* ------------------------------------------------------------------------ */
/* Utility routines */

static int MemoryError(void)
{
  if (currentDr)
    strcpy(gistError, "memory manager failed in Gd function");
  else
    strcpy(gistError, "currentDr not set in Gd function");
  return -1;
}

static void *Copy1(const void *orig, long size)
{
  void *px;
  if (size<=0) return 0;
  px= p_malloc(size);
  if (!px) MemoryError();
  else if (orig) memcpy(px, orig, size);
  return px;
}

static void *Copy2(void *x1, const void *orig1, const void *orig2, long size)
{
  void *x2, **x1p= (void **)x1;
  *x1p= Copy1(orig1, size);
  if (!*x1p) return 0;
  x2= Copy1(orig2, size);
  if (!x2) { p_free(*x1p);  *x1p= 0; }
  return x2;
}

void Gd_ScanZ(long n, const GpReal *z, GpReal *zmin, GpReal *zmax)
{
  long i;
  GpReal zn, zx;
  zn= zx= z[0];
  for (i=1 ; i<n ; i++) {
    if (z[i]<zn) zn= z[i];
    else if (z[i]>zx) zx= z[i];
  }
  *zmin= zn;
  *zmax= zx;
}

static void ScanXY(long n, const GpReal *x, const GpReal *y, GpBox *extrema)
{
  Gd_ScanZ(n, x, &extrema->xmin, &extrema->xmax);
  Gd_ScanZ(n, y, &extrema->ymin, &extrema->ymax);
}

void GeAddElement(int type, GdElement *element)
{
  GdElement *old;
  Drauing *drawing= currentDr;
  GeSystem *sys;

  if (drawing->cleared==1) ClearDrawing(drawing);
  sys= currentSy;

  old= sys? sys->elements : drawing->elements;
  if (!old) {  /* this is first element */
    if (sys) sys->elements= element;
    else drawing->elements= element;
    element->prev= element->next= element;
  } else {     /* insert element at end of ring */
    element->prev= old->prev;
    element->next= old;
    old->prev= element->prev->next= element;
  }
  element->ops= opTables + type;
  element->hidden= gistD.hidden;
  if (gistD.legend) {
    element->legend= Copy1(gistD.legend, strlen(gistD.legend)+1);
    /* sigh. ignore memory error here */
  } else {
    element->legend= 0;
  }
  element->number= drawing->nElements++;
  /* System nust always have number of its largest element for
     GdBeginSy to work properly */
  if (sys) sys->el.number= element->number;
  else Damage((GeSystem *)0, element);
}

void Gd_NextMeshBlock(long *ii, long *jj, long len, long iMax,
                      int *reg, int region)
{   /* Find next contiguous run of mesh points in given region */
  long i= *ii;
  long j= *jj;
  if (region==0) {
    for (j=i ; j<len ; j++)
      if (reg[j] || reg[j+1] || reg[j+iMax] || reg[j+iMax+1]) break;
    i= j;
    for (j=i+1 ; j<len ; j++)
      if (!reg[j] && !reg[j+1] && !reg[j+iMax] && !reg[j+iMax+1]) break;
  } else {
    for (j=i ; j<len ; j++)
      if (reg[j]==region || reg[j+1]==region ||
          reg[j+iMax]==region || reg[j+iMax+1]==region) break;
    i= j;
    for (j=i+1 ; j<len ; j++)
      if (reg[j]!=region && reg[j+1]!=region &&
          reg[j+iMax]!=region && reg[j+iMax+1]!=region) break;
  }
  *ii= i;
  *jj= j;
}

long GeGetMesh(int noCopy, GaQuadMesh *meshin, int region, void *vMeshEl)
{
  GeMesh *meshEl= vMeshEl;
  GaQuadMesh *mesh= &meshEl->mesh;
  GpBox *linBox= &meshEl->linBox;
  long iMax, jMax, i, j, len;
  int *reg;

  if (currentDr->cleared==1) ClearDrawing(currentDr);

  /* retrieve mesh shape from meshin */
  mesh->iMax= iMax= meshin->iMax;
  mesh->jMax= jMax= meshin->jMax;
  len= iMax*jMax;

  mesh->reg= 0;       /* set up for possible error return */
  mesh->triangle= 0;  /* only needed by GdContours, so handled there */

  meshEl->xlog= meshEl->ylog= 0;
  meshEl->region= region;

  meshEl->noCopy= noCopy & (NOCOPY_MESH|NOCOPY_COLORS|NOCOPY_UV|NOCOPY_Z);
  if (noCopy&NOCOPY_MESH) {
    /* Just copy pointers to mesh arrays -- NOCOPY also means not
       to free the pointer later */
    mesh->x= meshin->x;
    mesh->y= meshin->y;
  } else {
    /* Copy the mesh arrays themselves */
    mesh->y= Copy2(&mesh->x, meshin->x, meshin->y, sizeof(GpReal)*len);
    if (!mesh->y) { p_free(vMeshEl);  return 0; }
  }

  if ((noCopy&NOCOPY_MESH) && meshin->reg) {
    mesh->reg= reg= meshin->reg;
    meshEl->noCopy|= NOCOPY_REG;
  } else {
    mesh->reg= reg= Copy1(meshin->reg, sizeof(int)*(len+iMax+1));
    if (!reg) { Gd_KillMeshXY(vMeshEl);  p_free(vMeshEl);  return 0; }
  }

  /* Be sure region array is legal */
  for (i=0 ; i<iMax ; i++) reg[i]= 0;
  if (!meshin->reg) for (i=iMax ; i<len ; i++) reg[i]= 1;
  for (i=len ; i<len+iMax+1 ; i++) reg[i]= 0;
  for (i=0 ; i<len ; i+=iMax) reg[i]= 0;

  /* Scan mesh for extreme values */
  if (meshin->reg) {
    GpBox box;
    int first= 1;

    /* Use ScanXY on the longest contiguous run(s) of points */
    for (i=0 ; i<len ; ) {
      Gd_NextMeshBlock(&i, &j, len, iMax, reg, region);
      if (i>=len) break;
      ScanXY(j-i, mesh->x+i, mesh->y+i, &box);
      if (first) { *linBox= box;  first= 0; }
      else GpSwallow(linBox, &box);
      i= j+1;
    }
    if (first)
      linBox->xmin= linBox->xmax= linBox->ymin= linBox->ymax= 0.0;

  } else {
    ScanXY(len, mesh->x, mesh->y, linBox);
  }
  if (!currentSy) meshEl->el.box= *linBox;  /* for GeAddElement */

  /* copy mesh properties to gistD */
  Gd_MeshXYGet(vMeshEl);

  return len;
}

void GeMarkForScan(GdElement *el, GpBox *linBox)
{
  if (currentSy) {
    if (currentSy->unscanned<0) currentSy->unscanned= el->number;
  } else {
    el->box= *linBox;
  }
}

static int AutoMarker(GaLineAttribs *dl, int number)
{
  int p;
  p= (number+3)%4;
  if (number>=26) number%= 26;
  dl->mPhase= 0.25*(0.5+p)*dl->mSpace;
  dl->rPhase= 0.25*(0.5+p)*dl->rSpace;
  return 'A'+number;
}

/* ------------------------------------------------------------------------ */
/* Constructors for drawing elements are public routines declared in gist.h */

int GdLines(long n, const GpReal *px, const GpReal *py)
{
  GeLines *el;
  if (n<=0) return -1;
  el= currentDr? p_malloc(sizeof(GeLines)) : 0;
  if (!el) return MemoryError();
  el->xlog= el->ylog= 0;

  /* make private copies of x and y arrays */
  el->y= Copy2(&el->x, px, py, sizeof(GpReal)*n);
  if (!el->y) { p_free(el);  return -1; }
  el->n= n;

  /* scan for min and max of x and y arrays */
  ScanXY(n, px, py, &el->linBox);
  if (!currentSy) el->el.box= el->linBox;  /* for GeAddElement */

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->l= gistA.l;
  el->dl= gistA.dl;
  el->m= gistA.m;

  /* set base class members */
  GeAddElement(E_LINES, &el->el);
  if (gistA.m.type==0) el->m.type= AutoMarker(&el->dl, el->el.number);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.n= n;
  gistD.x= el->x;
  gistD.y= el->y;

  return el->el.number;
}

int GdDisjoint(long n, const GpReal *px, const GpReal *py,
               const GpReal *qx, const GpReal *qy)
{
  GeDisjoint *el;
  GpBox box;
  if (n<=0) return -1;
  el= currentDr? p_malloc(sizeof(GeDisjoint)) : 0;
  if (!el) return MemoryError();
  el->el.next= el->el.prev= 0;
  el->xlog= el->ylog= el->xqlog= el->yqlog= 0;

  /* make private copies of x, y, xq, yq arrays */
  el->y= Copy2(&el->x, px, py, sizeof(GpReal)*n);
  if (!el->y) { p_free(el);  return -1; }
  el->yq= Copy2(&el->xq, qx, qy, sizeof(GpReal)*n);
  if (!el->yq) { DisjointKill(el);  return -1; }
  el->n= n;

  /* scan for min and max of x and y arrays */
  ScanXY(n, px, py, &box);
  ScanXY(n, qx, qy, &el->linBox);
  GpSwallow(&el->linBox, &box);
  if (!currentSy) el->el.box= el->linBox;  /* for GeAddElement */

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->l= gistA.l;

  /* set base class members */
  GeAddElement(E_DISJOINT, &el->el);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.n= n;
  gistD.x= el->x;
  gistD.y= el->y;
  gistD.xq= el->xq;
  gistD.yq= el->yq;

  return el->el.number;
}

int GdText(GpReal x0, GpReal y0, const char *text, int toSys)
{
  GeText *el= currentDr? p_malloc(sizeof(GeText)) : 0;
  GeSystem *sys= currentSy;
  if (!el) return MemoryError();

  /* make private copy of text string */
  el->text= Copy1(text, strlen(text)+1);
  if (!el->text) { p_free(el);  return -1; }
  el->x0= x0;
  el->y0= y0;

  /* Without some sort of common font metric, there is no way to
     know the box associated with text.  Even with such a metric,
     the box would have to change with the coordinate transform,
     unlike any other element.  For now, punt.  */
  el->el.box.xmin= el->el.box.xmax= x0;
  el->el.box.ymin= el->el.box.ymax= y0;

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->t= gistA.t;

  /* set base class members */
  if (currentDr->cleared==1) ClearDrawing(currentDr); /* else currentSy */
  if (!toSys) currentSy= 0;                      /* can be clobbered... */
  GeAddElement(E_TEXT, &el->el);
  if (currentSy && currentSy->unscanned<0)
    currentSy->unscanned= el->el.number;
  if (!toSys) currentSy= sys;

  /* copy properties to gistD */
  gistD.x0= el->x0;
  gistD.y0= el->y0;
  gistD.text= el->text;

  return el->el.number;
}

int GdCells(GpReal px, GpReal py, GpReal qx, GpReal qy,
            long width, long height, long nColumns, const GpColor *colors)
{
  GpReal x[2], y[2];
  GpBox linBox;
  long ncells= width*height;
  long len= sizeof(GpColor)*ncells;
  GeCells *el= currentDr? p_malloc(sizeof(GeCells)) : 0;
  if (!el) return MemoryError();

  /* make private copy of colors array */
  el->rgb = gistA.rgb;
  if (gistA.rgb) len *= 3;
  gistA.rgb = 0;
  el->colors= p_malloc(len);
  if (!el->colors) { p_free(el); return MemoryError(); }
  el->px= x[0]= px;
  el->py= y[0]= py;
  el->qx= x[1]= qx;
  el->qy= y[1]= qy;
  el->width= width;
  el->height= height;
  if (nColumns==width) {
    memcpy(el->colors, colors, len);
  } else {
    GpColor *newcols= el->colors;
    long i, rowSize= sizeof(GpColor)*width;
    for (i=0 ; i<height ; i++) {
      memcpy(newcols, colors, rowSize);
      newcols+= width;
      colors+= nColumns;
    }
  }
  ScanXY(2L, x, y, &linBox);
  if (!currentSy) el->el.box= linBox;  /* for GeAddElement */

  /* set base class members */
  GeAddElement(E_CELLS, &el->el);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &linBox);

  /* copy properties to gistD */
  gistD.px= el->px;
  gistD.py= el->py;
  gistD.qx= el->qx;
  gistD.qy= el->qy;
  gistD.width= el->width;
  gistD.height= el->height;
  gistD.colors= el->colors;

  return el->el.number;
}

int GdFill(long n, const GpColor *colors, const GpReal *px,
           const GpReal *py, const long *pn)
{
  GePolys *el;
  long i, ntot;
  if (n<=0) return -1;
  el= currentDr? p_malloc(sizeof(GePolys)) : 0;
  if (!el) return MemoryError();
  el->xlog= el->ylog= 0;

  /* make private copy of colors array */
  if (colors) {
    long ncol = (gistA.rgb? 3*n : n);
    el->rgb = gistA.rgb;
    el->colors= p_malloc(ncol);
    if (!el->colors) { p_free(el); return MemoryError(); }
    memcpy(el->colors, colors, ncol);
  } else {
    el->rgb = 0;
    el->colors= 0;
  }
  gistA.rgb = 0;

  /* make private copy of lengths array */
  el->pn= p_malloc(sizeof(long)*n);
  if (!el->pn) { p_free(el->colors); p_free(el); return MemoryError(); }
  for (ntot=i=0 ; i<n ; i++) {
    el->pn[i]= pn[i];
    ntot+= pn[i];
  }

  /* make private copies of x and y arrays */
  el->y= Copy2(&el->x, px, py, sizeof(GpReal)*ntot);
  if (!el->y) { p_free(el->pn); p_free(el->colors); p_free(el);  return -1; }
  el->n= n;

  /* scan for min and max of x and y arrays */
  if (n<2 || pn[1]>1) ScanXY(ntot, px, py, &el->linBox);
  else ScanXY(ntot-pn[0], px+pn[0], py+pn[0], &el->linBox);
  if (!currentSy) el->el.box= el->linBox;  /* for GeAddElement */

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->e= gistA.e;  /* for edges */

  /* set base class members */
  GeAddElement(E_POLYS, &el->el);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.n= n;
  gistD.x= el->x;
  gistD.y= el->y;
  gistD.pn= el->pn;
  gistD.colors= el->colors;

  return el->el.number;
}

int GdMesh(int noCopy, GaQuadMesh *mesh, int region, int boundary,
           int inhibit)
{
  long len;
  GeMesh *el= currentDr? p_malloc(sizeof(GeMesh)) : 0;
  if (!el) return MemoryError();
  el->el.next= el->el.prev= 0;

  /* get mesh */
  len= GeGetMesh(noCopy, mesh, region, el);
  if (!len) return -1;
  el->boundary= boundary;
  el->inhibit= inhibit;

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->l= gistA.l;

  /* set base class members */
  GeAddElement(E_MESH, &el->el);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.boundary= el->boundary;
  gistD.inhibit= el->inhibit;

  return el->el.number;
}

int GdFillMesh(int noCopy, GaQuadMesh *mesh, int region,
               GpColor *colors, long nColumns)
{
  long len;
  GeFill *el= currentDr? p_malloc(sizeof(GeFill)) : 0;
  if (!el) return MemoryError();
  el->el.next= el->el.prev= 0;

  /* get mesh */
  len= GeGetMesh(noCopy, mesh, region, el);
  if (!len) return -1;

  /* make private copy of colors array */
  el->rgb = gistA.rgb;
  if (noCopy&NOCOPY_COLORS || !colors) {
    el->colors= colors;
  } else {
    long iMax1= mesh->iMax-1;
    long len1= len - mesh->jMax - iMax1;
    int rgb = gistA.rgb;
    if (rgb) len1 *= 3;
    el->colors= Copy1(nColumns==iMax1?colors:0, sizeof(GpColor)*len1);
    if (!el->colors) { FilledKill(el);  return -1; }
    if (nColumns!=iMax1) {
      long i, j=0, k=0;
      for (i=0 ; i<len1 ; i++) {
        if (rgb) {
          el->colors[i++]= colors[3*(j+k)];
          el->colors[i++]= colors[3*(j+k)+1];
          el->colors[i]= colors[3*(j+k)+2];
        } else {
          el->colors[i]= colors[j+k];
        }
        j++;
        if (j==iMax1) { k+= nColumns; j= 0; }
      }
      nColumns= iMax1;
    }
  }
  gistA.rgb = 0;
  el->nColumns= nColumns;

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->e= gistA.e;  /* for edges */

  /* set base class members */
  GeAddElement(E_FILLED, &el->el);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.nColumns= nColumns;
  gistD.colors= el->colors;

  return el->el.number;
}

int GdVectors(int noCopy, GaQuadMesh *mesh, int region,
              GpReal *u, GpReal *v, GpReal scale)
{
  long len;
  GeVectors *el= currentDr? p_malloc(sizeof(GeVectors)) : 0;
  if (!el) return MemoryError();
  el->el.next= el->el.prev= 0;

  /* get mesh */
  len= GeGetMesh(noCopy, mesh, region, el);
  if (!len) return -1;

  /* make private copy of (u,v) arrays */
  if (noCopy&NOCOPY_UV) {
    el->u= u;
    el->v= v;
  } else {
    el->v= Copy2(&el->u, u, v, sizeof(GpReal)*len);
    if (!el->v) { VectorsKill(el);  return -1; }
  }
  el->scale= scale;

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->l= gistA.l;
  el->f= gistA.f;
  el->vect= gistA.vect;

  /* set base class members */
  GeAddElement(E_VECTORS, &el->el);

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.u= el->u;
  gistD.v= el->v;
  gistD.scale= el->scale;

  return el->el.number;
}

int Gd_MakeContours(GeContours *con)
{
  long n;
  GpReal dphase, *px, *py;
  int i;
  GeLines *group;
  GeLines *el, *prev;
  int marker;

  /* Generic properties copied to all contours (m.type reset below)  */
  gistA.l= con->l;
  gistA.dl= con->dl;
  gistA.m= con->m;
  marker= gistA.m.type>32? gistA.m.type : 'A';
  dphase= 0.25*con->dl.mSpace;

  for (i=0 ; i<con->nLevels ; i++) con->groups[i]= 0;

  for (i=0 ; i<con->nLevels ; i++) {
    gistA.m.type= marker++;
    if (marker=='Z'+1 || marker=='z'+1) marker= 'A';
    group= prev= 0;
    if (GaContourInit(&con->mesh, con->region, con->z, con->levels[i])) {
      while (GaContour(&n, &px, &py, &gistA.dl.closed)) {
        el= currentDr? p_malloc(sizeof(GeLines)) : 0;
        if (!el) return MemoryError();

        /* make private copies of x and y arrays */
        el->y= Copy2(&el->x, px, py, sizeof(GpReal)*n);
        if (!el->y) { p_free(el);  return -1; }
        el->n= n;
        el->xlog= el->ylog= 0;

        /* scan for min and max of x and y arrays */
        ScanXY(n, px, py, &el->linBox);
        if (!currentSy) el->el.box= el->linBox;  /* for GeAddElement */

        el->el.ops= opTables + E_LINES;
        el->el.hidden= 0;
        el->el.legend= 0;
        el->el.box= el->linBox;

        /* GeContour number always matches number of largest element
           for GdBeginSy to work properly */
        el->el.number= con->el.number= currentDr->nElements++;

        if (prev) {
          prev->el.next= group->el.prev= &el->el;
          el->el.prev= &prev->el;
          el->el.next= &group->el;
        } else {
          con->groups[i]= group= el;
          el->el.next= el->el.prev= &el->el;
        }
        prev= el;

        el->l= gistA.l;
        el->dl= gistA.dl;
        el->m= gistA.m;

        gistA.dl.mPhase+= dphase;
        if (gistA.dl.mPhase>gistA.dl.mSpace)
          gistA.dl.mPhase-= gistA.dl.mSpace;
      }
    }
  }

  return 0;
}

int GdContours(int noCopy, GaQuadMesh *mesh, int region,
               GpReal *z, const GpReal *levels, int nLevels)
{
  long len;
  GeContours *el= currentDr? p_malloc(sizeof(GeContours)) : 0;
  if (!el) return MemoryError();
  el->el.next= el->el.prev= 0;
  el->z= el->levels= 0;
  el->groups= 0;

  /* get mesh */
  len= GeGetMesh(noCopy, mesh, region, el);
  if (!len) return -1;

  /* make private copy of z and levels arrays */
  if (noCopy&NOCOPY_Z) {
    el->z= z;
  } else {
    el->z= Copy1(z, sizeof(GpReal)*len);
    if (!el->z) { ContoursKill(el);  return -1; }
  }

  /* create triangle array now if necessary */
  if (noCopy&NOCOPY_MESH && mesh->triangle) {
    el->mesh.triangle= mesh->triangle;
    el->noCopy|= NOCOPY_TRI;
    gistD.noCopy|= NOCOPY_TRI;
  } else {
    el->mesh.triangle= Copy1(mesh->triangle, sizeof(short)*len);
    if (!el->mesh.triangle) { ContoursKill(el);  return -1; }
  }

  /* copy relevant attributes from gistA
     -- This must be done BEFORE GeAddElement, since the damage
        calculation depends on the Margins! */
  el->l= gistA.l;
  el->dl= gistA.dl;
  el->m= gistA.m;

  /* set base class members */
  GeAddElement(E_CONTOURS, &el->el);
  if (gistA.m.type==0) el->m.type= AutoMarker(&el->dl, el->el.number);

  el->nLevels= nLevels;
  if (nLevels>0) {
    el->levels= Copy1(levels, sizeof(GpReal)*nLevels);
    if (!el->levels) { ContoursKill(el);  return -1; }
    el->groups= (GeLines **)p_malloc(sizeof(GeLines *)*nLevels);
    if (!el->groups || Gd_MakeContours(el)) { ContoursKill(el);  return -1; }
  } else {
    nLevels= 0;
    el->levels= 0;
    el->groups= 0;
  }

  /* current box not set, mark as unscanned if in system */
  GeMarkForScan(&el->el, &el->linBox);

  /* copy properties to gistD */
  gistD.z= el->z;
  gistD.nLevels= el->nLevels;
  gistD.levels= el->levels;

  return el->el.number;
}

/* ------------------------------------------------------------------------ */

static void GuessBox(GpBox *box, GpBox *viewport, GaTickStyle *ticks)
{
  GpReal dxmin= 0.0, xmin= viewport->xmin;
  GpReal dxmax= 0.0, xmax= viewport->xmax;
  GpReal dymin= 0.0, ymin= viewport->ymin;
  GpReal dymax= 0.0, ymax= viewport->ymax;
  int vf= ticks->vert.flags;
  int hf= ticks->horiz.flags;
  GpReal vlen= ((((vf&TICK_IN)&&(vf&TICK_OUT))||(vf&TICK_C))? 0.5 : 1.0) *
    ticks->vert.tickLen[0];
  GpReal hlen= ((((vf&TICK_IN)&&(vf&TICK_OUT))||(vf&TICK_C))? 0.5 : 1.0) *
    ticks->horiz.tickLen[0];
  GpReal cy= ticks->horiz.textStyle.height;
  GpReal cx= ticks->vert.textStyle.height*0.6;  /* guess at char width */
  GpReal hx= cy*0.6;                            /* guess at char width */
  GpBox overflow;

  /* Note-- extra 0.4 for nudged log decades (see DrawXLabels, tick.c) */
  cx*= ticks->vert.nDigits+2.4;         /* largest width of y label */
  hx*= 0.5*(ticks->horiz.nDigits+2.4);  /* maximum distance x label
                                           can project past xmin, xmax */

  if (((vf&TICK_L)&&(vf&TICK_OUT)) || (vf&TICK_C))
    dxmin= ticks->vert.tickOff+vlen;
  if (((vf&TICK_U)&&(vf&TICK_OUT)) || (vf&TICK_C))
    dxmax= ticks->vert.tickOff+vlen;
  if (((hf&TICK_L)&&(hf&TICK_OUT)) || (hf&TICK_C))
    dymin= ticks->horiz.tickOff+hlen;
  if (((hf&TICK_U)&&(hf&TICK_OUT)) || (hf&TICK_C))
    dymax= ticks->horiz.tickOff+hlen;

  if (vf & LABEL_L) xmin-= ticks->vert.labelOff+cx;
  else if ((hf&(LABEL_L|LABEL_U)) && hx>dxmin) xmin-= hx;
  else xmin-= dxmin;
  if (vf & LABEL_U) xmax+= ticks->vert.labelOff+cx;
  else if ((hf&(LABEL_L|LABEL_U)) && hx>dxmax) xmax+= hx;
  else xmax+= dxmax;

  if (hf & LABEL_L) ymin-= ticks->horiz.labelOff+2.0*cy;
  else if ((vf&(LABEL_L|LABEL_U)) && 0.5*cy>dymin) ymin-= 0.5*cy;
  else ymin-= dymin;
  if (hf & LABEL_U) xmax+= ticks->horiz.labelOff+2.0*cy;
  else if ((vf&(LABEL_L|LABEL_U)) && 0.5*cy>dymax) ymax+= 0.5*cy;
  else ymax+= dymax;

  if (vf & (TICK_L|TICK_U)) {
    xmin-= 0.5*ticks->vert.tickStyle.width*DEFAULT_LINE_WIDTH;
    xmax+= 0.5*ticks->vert.tickStyle.width*DEFAULT_LINE_WIDTH;
  }
  if (hf & (TICK_L|TICK_U)) {
    ymin-= 0.5*ticks->horiz.tickStyle.width*DEFAULT_LINE_WIDTH;
    ymax+= 0.5*ticks->horiz.tickStyle.width*DEFAULT_LINE_WIDTH;
  }

  box->xmin= xmin;
  box->xmax= xmax;
  box->ymin= ymin;
  box->ymax= ymax;

  /* Finally, swallow overflow boxes, assuming 22 characters max */
  overflow.xmin= ticks->horiz.xOver;
  overflow.ymin= ticks->horiz.yOver-ticks->horiz.textStyle.height*0.2;
  overflow.xmax= overflow.xmin+ticks->horiz.textStyle.height*(0.6*22.0);
  overflow.ymax= overflow.ymin+ticks->horiz.textStyle.height;
  GpSwallow(box, &overflow);
  overflow.xmin= ticks->vert.xOver;
  overflow.ymin= ticks->vert.yOver-ticks->vert.textStyle.height*0.2;
  overflow.xmax= overflow.xmin+ticks->vert.textStyle.height*(0.6*22.0);
  overflow.ymax= overflow.ymin+ticks->vert.textStyle.height;
  GpSwallow(box, &overflow);
}

int GdNewSystem(GpBox *viewport, GaTickStyle *ticks)
{
  GeSystem *sys;
  int sysIndex;

  if (!currentDr) return -1;

  /* Adding a new system clears the drawing */
  if (currentDr->cleared!=2) ClearDrawing(currentDr);
  sysIndex= currentDr->nSystems+1;

  sys= p_malloc(sizeof(GeSystem));
  if (!sys) return -1;
  sys->el.ops= opTables + E_SYSTEM;
  if (gistD.legend) {
    sys->el.legend= Copy1(gistD.legend, strlen(gistD.legend)+1);
    if (!sys->el.legend) { p_free(sys);  return -1; }
  } else sys->el.legend= 0;
  sys->el.hidden= gistD.hidden;

  if (sysIndex>1) {
    GdElement *prev= currentDr->systems->el.prev;
    prev->next= &sys->el;
    sys->el.prev= prev;
    sys->el.next= &currentDr->systems->el;
    currentDr->systems->el.prev= &sys->el;
  } else {
    sys->el.prev= sys->el.next= &sys->el;
    currentDr->systems= sys;
  }
  sys->el.number= -1;
  currentDr->nSystems++;
  sys->rescan= 0;
  sys->unscanned= -1;

  GuessBox(&sys->el.box, viewport, ticks);

  if (viewport->xmin<viewport->xmax) {
    sys->trans.viewport.xmin= viewport->xmin;
    sys->trans.viewport.xmax= viewport->xmax;
  } else {
    sys->trans.viewport.xmin= viewport->xmax;
    sys->trans.viewport.xmax= viewport->xmin;
  }
  if (viewport->ymin<viewport->ymax) {
    sys->trans.viewport.ymin= viewport->ymin;
    sys->trans.viewport.ymax= viewport->ymax;
  } else {
    sys->trans.viewport.ymin= viewport->ymax;
    sys->trans.viewport.ymax= viewport->ymin;
  }
  sys->trans.window.xmin= sys->trans.window.ymin= 0.0;
  sys->trans.window.xmax= sys->trans.window.ymax= 1.0;
  sys->ticks= *ticks;
  sys->flags= D_XMIN | D_XMAX | D_YMIN | D_YMAX;
  sys->elements= 0;
  sys->savedWindow.xmin= sys->savedWindow.ymin= 0.0;
  sys->savedWindow.xmax= sys->savedWindow.ymax= 1.0;
  sys->savedFlags= D_XMIN | D_XMAX | D_YMIN | D_YMAX;

  sys->xtick= sys->ytick= 0;
  sys->xlabel= sys->ylabel= 0;

  GdSetSystem(sysIndex);
  return sysIndex;
}

int GdGetSystem(void)
{
  GdElement *sys, *sys0;
  int sysIndex= 0;
  if (!currentDr) return -1;
  if (!currentDr->systems || !currentSy) return 0;

  /* ClearDrawing sets currentSy-- must not be pending */
  if (currentDr->cleared==1) ClearDrawing(currentDr);

  sys= sys0= (GdElement *)currentDr->systems;
  for (sysIndex=1 ; sys!=&currentSy->el ; sysIndex++) {
    sys= sys->next;
    if (sys==sys0) return -2;
  }

  return sysIndex;
}

int GdSetSystem(int sysIndex)
{
  GeSystem *sys;
  GdElement *sys0;
  if (!currentDr || !currentDr->systems) return E_NONE;

  /* ClearDrawing sets currentSy-- must not be pending */
  if (currentDr->cleared==1) ClearDrawing(currentDr);

  currentEl= 0;
  currentCn= -1;
  if (sysIndex<1) {  /* Set current system to none */
    currentSy= 0;
    gistD.trans.viewport.xmin= gistD.trans.viewport.xmax=
      gistD.trans.viewport.ymin= gistD.trans.viewport.ymax= 0.0;
    gistD.flags= 0;
    return E_NONE;
  }

  sys= currentDr->systems;
  sys0= &sys->el;
  while (--sysIndex && sys->el.next!=sys0)
    sys= (GeSystem *)sys->el.next;
  if (sysIndex>0) return E_NONE;

  currentSy= sys;
  gistD.hidden= sys->el.hidden;
  gistD.legend= sys->el.legend;
  gistD.ticks= sys->ticks;
  gistD.trans.viewport= sys->trans.viewport;
  if (GdGetLimits()) {
    SystemKill(sys);
    return E_NONE;
  }
  return E_SYSTEM;
}

int GdSetElement(int elIndex)
{
  GdElement *el, *el0;
  if (!currentDr) return E_NONE;

  el= currentSy? currentSy->elements : currentDr->elements;

  if (elIndex<0 || !el) {   /* set current element to none */
    currentEl= 0;
    currentCn= -1;
    return E_NONE;
  }

  el0= el;
  while (elIndex-- && el->next!=el0) el= el->next;
  if (elIndex>=0) return E_NONE;

  currentEl= el;
  currentCn= -1;
  return el->ops->GetProps(el);
}

static GdElement *NextConCurve(GdElement *el)
{
  GeContours *con= (GeContours *)currentEl;
  GdElement *el0= &con->groups[currentCn]->el;
  if (!el) el= el0;
  else if (el->next==el0) el= 0;
  else el= el->next;
  return el;
}

int GdSetContour(int levIndex)
{
  GeContours *con;
  GdElement *el;
  if (!currentDr || !currentEl || currentEl->ops->type!=E_CONTOURS)
    return E_NONE;
  con= (GeContours *)currentEl;
  if (levIndex>=con->nLevels) return E_NONE;
  if (levIndex<0) {
    currentCn= -1;
    return E_NONE;
  }
  currentCn= levIndex;
  el= NextConCurve((GdElement *)0);
  if (el) LinesGet(el);
  else ContoursGet(con);
  return E_LINES;
}

static int GeFindIndex(int id, GeSystem *sys)
{
  int elIndex;
  GdElement *el, *el0;
  if (!currentDr) return -1;

  el= sys? sys->elements : currentDr->elements;
  if (!el) return -1;

  el0= el;
  elIndex= 0;
  while (el->number != id) {
    if (el->next==el0) return -1;
    el= el->next;
    elIndex++;
  }

  return elIndex;
}

int GdFindIndex(int id)
{
  return GeFindIndex(id, currentSy);
}

int GdFindSystem(int id)
{
  int sysIndex;
  GeSystem *sys;
  GdElement *sys0;
  if (!currentDr) return -1;

  if (GeFindIndex(id, 0)>=0) return 0;
  sys= currentDr->systems;
  if (!sys) return -1;
  sys0= &sys->el;
  sysIndex= 1;
  while (GeFindIndex(id, sys)<0) {
    if (sys->el.next==sys0) return -1;
    sys= (GeSystem *)sys->el.next;
    sysIndex++;
  }
  return sysIndex;
}

int GdAltTick(GaAltTicks *xtick, GaAltLabel *xlabel,
              GaAltTicks *ytick, GaAltLabel *ylabel)
{
  if (!currentDr || !currentSy) return 1;
  if (xtick) currentSy->xtick= xtick;
  if (ytick) currentSy->ytick= ytick;
  if (xlabel) currentSy->xlabel= xlabel;
  if (ylabel) currentSy->ylabel= ylabel;
  return 0;
}

/* ------------------------------------------------------------------------ */

int GdEdit(int xyzChanged)
{
  GdElement *el= currentEl;
  if (!currentDr || !el) return 1;

  /* Changing linestyles or most other simple changes may incur damage
     in a way that is very difficult to anticipate, hence, must call
     Damage here.  On the other hand, the elements need to be rescanned
     only if CHANGE_XY or CHANGE_Z has been set, so only set rescan
     in this case.  If only linestyles and such have been changed, GdScan
     will not call Damage again, although if they have it will-- hence, if
     you are changing both coordinates and linestyles, there is no way to
     avoid "double damage".  */
  Damage(currentSy, el);
  if (currentSy && xyzChanged) currentSy->rescan= 1;

  if (currentCn>=0) {
    el= NextConCurve((GdElement *)0);
    if (el) {
      /* legend only changes on first piece of contour */
      el->legend= gistD.legend;
      Gd_LinesSubSet(el);
      /* other line properties propagate to all contour pieces
         -- but NEVER attempt to change the (x,y) values */
      while ((el= NextConCurve(el))) Gd_LinesSubSet(el);
    }
    return 0;
  }
  return el->ops->SetProps(el, xyzChanged);
}

int GdRemove(void)
{
  GdElement *el= currentEl;
  if (!currentDr || !el || currentCn>=0) return 1;

  /* Damage alert must take place here-- unfortunately, if this remove
     changes extreme values, a second call to Damage will be made in
     GdScan.  Hopefully, GdRemove is a rare enough operation that this
     inefficiency is negligible.  */
  Damage(currentSy, el);

  if (currentSy) {
    GdElement *prev= el->prev;
    if (el==prev) {
      currentSy->unscanned= -1;
      currentSy->rescan= 0;
      currentSy->el.number= -1;
    } else {
      if (el->number==currentSy->unscanned) {
        if (el->next != currentSy->elements)
          currentSy->unscanned= el->next->number;
        else currentSy->unscanned= -1;
      }
      if (el->number<currentSy->unscanned && !el->hidden)
        currentSy->rescan= 1;
      if (el->number==currentSy->el.number)
        currentSy->el.number= prev->number;
    }
  }

  if (currentSy && el==currentSy->elements) {
    if (el->next==el) currentSy->elements= 0;
    else currentSy->elements= el->next;
  } else if (el==currentDr->elements) {
    if (el->next==el) currentDr->elements= 0;
    else currentDr->elements= el->next;
  }

  el->ops->Kill(el);
  currentEl= 0;
  return 0;
}

int GdGetLimits(void)
{
  if (!currentDr || !currentSy) return 1;
  if ((currentSy->rescan || currentSy->unscanned>=0)
      && GdScan(currentSy)) return 1;  /* memory manager failure */
  gistD.trans.window= currentSy->trans.window;
  gistD.flags= currentSy->flags;

  if (gistD.flags & D_LOGX) {
    gistD.limits.xmin= exp10(gistD.trans.window.xmin);
    gistD.limits.xmax= exp10(gistD.trans.window.xmax);
  } else {
    gistD.limits.xmin= gistD.trans.window.xmin;
    gistD.limits.xmax= gistD.trans.window.xmax;
  }
  if (gistD.flags & D_LOGY) {
    gistD.limits.ymin= exp10(gistD.trans.window.ymin);
    gistD.limits.ymax= exp10(gistD.trans.window.ymax);
  } else {
    gistD.limits.ymin= gistD.trans.window.ymin;
    gistD.limits.ymax= gistD.trans.window.ymax;
  }
  return 0;
}

int GdSetLimits(void)
{
  int flags, rescan;
  if (!currentDr || !currentSy) return 1;

  if (gistD.flags & D_LOGX) {
    gistD.trans.window.xmin= SAFELOG(gistD.limits.xmin);
    gistD.trans.window.xmax= SAFELOG(gistD.limits.xmax);
  } else {
    gistD.trans.window.xmin= gistD.limits.xmin;
    gistD.trans.window.xmax= gistD.limits.xmax;
  }
  if (gistD.flags & D_LOGY) {
    gistD.trans.window.ymin= SAFELOG(gistD.limits.ymin);
    gistD.trans.window.ymax= SAFELOG(gistD.limits.ymax);
  } else {
    gistD.trans.window.ymin= gistD.limits.ymin;
    gistD.trans.window.ymax= gistD.limits.ymax;
  }

  flags= currentSy->flags;
  currentSy->flags= gistD.flags;

  /* Normally, setting the limits damages the entire system.
     However, would like to allow a special case for fixing limits to
     their existing extreme values.  */
  rescan= 1;
  if ( ! ( (flags^gistD.flags) & (~(D_XMIN|D_XMAX|D_YMIN|D_YMAX)) ) ) {
    if (((flags&D_XMIN)==(gistD.flags&D_XMIN) || (gistD.flags&D_XMIN)==0) &&
        ((flags&D_XMAX)==(gistD.flags&D_XMAX) || (gistD.flags&D_XMAX)==0) &&
        ((flags&D_YMIN)==(gistD.flags&D_YMIN) || (gistD.flags&D_YMIN)==0) &&
        ((flags&D_YMAX)==(gistD.flags&D_YMAX) || (gistD.flags&D_YMAX)==0)) {
      GpBox *w= &currentSy->trans.window;
      if (w->xmin==gistD.trans.window.xmin &&
          w->xmax==gistD.trans.window.xmax &&
          w->ymin==gistD.trans.window.ymin &&
          w->ymax==gistD.trans.window.ymax)
        rescan= 0;
    }
  }
  currentSy->trans.window= gistD.trans.window;
  currentSy->rescan|= rescan;

  /* damage alert takes place in GdScan just before rendering */
  return 0;
}

int GdSaveLimits(int resetZoomed)
{
  if (!currentDr || !currentSy) return 1;
  currentSy->savedWindow= currentSy->trans.window;
  currentSy->savedFlags= currentSy->flags;
  if (resetZoomed) currentSy->savedFlags&= ~D_ZOOMED;
  return 0;
}

int GdRevertLimits(int ifZoomed)
{
  if (!currentDr || !currentSy ||
      (ifZoomed && !(currentSy->flags&D_ZOOMED))) return 1;
  if (currentSy->savedFlags!=currentSy->flags ||
      currentSy->savedWindow.xmin!=currentSy->trans.window.xmin ||
      currentSy->savedWindow.xmax!=currentSy->trans.window.xmax ||
      currentSy->savedWindow.ymin!=currentSy->trans.window.ymin ||
      currentSy->savedWindow.ymax!=currentSy->trans.window.ymax) {
    currentSy->trans.window= currentSy->savedWindow;
    currentSy->flags= currentSy->savedFlags;
    currentSy->rescan= 1;
  }
  return 0;
}

int GdSetPort(void)
{
  GpBox *v, oldBox;
  if (!currentDr || !currentSy) return 1;

  currentSy->el.hidden= gistD.hidden;
  /*
  if (gistD.legend) {
    currentSy->el.legend= Copy1(gistD.legend, strlen(gistD.legend)+1);
  }
  */

  /* First, damage current coordinate system box.  */
  Damage(currentSy, (GdElement *)0);

  /* Save old box, set new ticks, viewport, and correponding box */
  oldBox= currentSy->el.box;
  v= &currentSy->trans.viewport;
  currentSy->ticks= gistD.ticks;
  *v= gistD.trans.viewport;
  GuessBox(&currentSy->el.box, &gistD.trans.viewport, &gistD.ticks);

  /* Since stacking order hasn't changed, new box must be damaged
     if it is not contained in the old box.  */
  v= &currentSy->el.box;
  if (v->xmin<oldBox.xmin || v->xmax>oldBox.xmax ||
      v->ymin<oldBox.ymin || v->ymax>oldBox.ymax)
    Damage(currentSy, (GdElement *)0);

  return 0;
}

GpBox *GdClearSystem(void)
{
  GpBox *dBox;
  int n, nel;
  GeSystem *sys, *sys0;
  GdElement *el, *el0;
  /* Intended for use with animation... */
  if (!currentDr || !currentSy) return 0;

  Gd_KillRing(currentSy->elements);
  currentSy->elements= 0;
  currentSy->el.number= currentSy->unscanned= -1;
  currentSy->rescan= 0;

  sys0= currentDr->systems;
  el0= currentDr->elements;
  nel= -1;
  if ((sys= sys0)) do {
    if (sys==currentSy) continue;
    n= currentSy->el.number;
    if (n>nel) nel= n;
    sys= (GeSystem *)sys->el.next;
  } while (sys!=sys0);
  if ((el= el0)) do {
    n= el->number;
    if (n>nel) nel= n;
    el= el->next;
  } while (el!=el0);
  currentDr->nElements= nel+1;

  if (currentSy->flags & (D_XMIN|D_XMAX|D_YMIN|D_YMAX)) {
    /* Some extreme value set, damage whole box */
    dBox= &currentSy->el.box;
    Damage(currentSy, (GdElement *)0);
  } else {
    /* All limits fixed, damage only viewport */
    dBox= &currentSy->trans.viewport;
    Damage(currentSy, &currentSy->el);
  }

  return dBox;
}

/* ------------------------------------------------------------------------ */
