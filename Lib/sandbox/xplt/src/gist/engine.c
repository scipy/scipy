/*
 * ENGINE.C
 *
 * $Id$
 *
 * Implement common properties of all GIST engines
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "gist.h"
#include "engine.h"
#include "draw.h"
#include "pstdlib.h"

Engine *gistEngines= 0;
Engine *gistActive= 0;
Engine *gistPreempt= 0;

#include <string.h>

static void DefaultClearArea(Engine *engine, GpBox *box);
static void MoreScratch(long np, long ns);

/* ------------------------------------------------------------------------ */

/* ARGSUSED */
static void DefaultClearArea(Engine *engine, GpBox *box)
{
  /* Default ClearArea triggers complete redraw */
  engine->Clear(engine, CONDITIONALLY);
  engine->lastDrawn= -1;
  engine->systemsSeen[0]= engine->systemsSeen[1]= 0;
  engine->damaged= engine->inhibit= 0;
}

Engine *GpNewEngine(long size, char *name, char *type,
                    GpTransform *transform, int landscape,
  void (*Kill)(Engine*), int (*Clear)(Engine*,int), int (*Flush)(Engine*),
  void (*ChangeMap)(Engine*), int (*ChangePalette)(Engine*),
  int (*DrawLines)(Engine*,long,const GpReal*,const GpReal*,int,int),
  int (*DrawMarkers)(Engine*,long,const GpReal*,const GpReal *),
  int (*DrwText)(Engine*e,GpReal,GpReal,const char*),
  int (*DrawFill)(Engine*,long,const GpReal*,const GpReal*),
  int (*DrawCells)(Engine*,GpReal,GpReal,GpReal,GpReal,
                   long,long,long,const GpColor*),
  int (*DrawDisjoint)(Engine*,long,const GpReal*,const GpReal*,
                      const GpReal*,const GpReal*))
{
  long lname= name? strlen(name) : 0;
  Engine *engine;
  /* For Electric Fence package and maybe others, it is nice to ensure
     that size of block allocated for Engine is a multiple of the size
     of the most restrictively aligned object which can be in any
     Engine; assume this is a double.  */
  lname= (lname/sizeof(double) + 1)*sizeof(double);  /* >= lname+1 */
  engine= (Engine *)p_malloc(size+lname);
  if (!engine) return 0;

  /* Fill in Engine properties, link into gistEngines list */
  engine->next= gistEngines;
  gistEngines= engine;
  engine->nextActive= 0;
  engine->name= (char *)engine + size;
  strcpy(name? engine->name : "", name);
  engine->type= type;
  engine->active= 0;
  engine->marked= 0;

  engine->transform= *transform;
  engine->landscape= landscape? 1 : 0;
  GpDeviceMap(engine);
  /* (a proper map will be installed when the engine is activated) */
  engine->map.x.scale= engine->map.y.scale= 1.0;
  engine->map.x.offset= engine->map.y.offset= 0.0;

  /* No pseudocolor map initially */
  engine->colorChange= 0;
  engine->colorMode= 0;
  engine->nColors= 0;
  engine->palette= 0;

  /* No associated drawing initially */
  engine->drawing= 0;
  engine->lastDrawn= -1;
  engine->systemsSeen[0]= engine->systemsSeen[1]= 0;
  engine->inhibit= 0;
  engine->damaged= 0;  /* causes Clear if no ClearArea virtual function */
  engine->damage.xmin= engine->damage.xmax=
    engine->damage.ymin= engine->damage.ymax= 0.0;

  /* Fill in virtual function table */
  engine->Kill= Kill;
  engine->Clear= Clear;
  engine->Flush= Flush;
  engine->ChangeMap= ChangeMap;
  engine->ChangePalette= ChangePalette;
  engine->DrawLines= DrawLines;
  engine->DrawMarkers= DrawMarkers;
  engine->DrwText= DrwText;
  engine->DrawFill= DrawFill;
  engine->DrawCells= DrawCells;
  engine->DrawDisjoint= DrawDisjoint;
  engine->ClearArea= &DefaultClearArea;   /* damage causes complete redraw */

  return engine;
}

void GpDelEngine(Engine *engine)
{
  Engine *eng= gistEngines;
  if (!engine) return;

  /* Unlink from gistEngines list */
  if (engine->active) GpDeactivate(engine);
  if (eng==engine) gistEngines= engine->next;
  else {
    /* Because of recursive deletes necessary to deal with X window
       deletions (see xbasic.c:ShutDown, hlevel.c:ShutDownDev), if
       the engine has already been removed from the list, it means that
       this routine is being called for the second time for this engine,
       and p_free must NOT be called.  Fix this someday.  */
    while (eng && eng->next!=engine) eng= eng->next;
    if (!eng) return;
    eng->next= engine->next;
  }

  p_free(engine);
}

/* ------------------------------------------------------------------------ */

void GpKillEngine(Engine *engine)
{
  if (engine) engine->Kill(engine);
}

int GpActivate(Engine *engine)
{
  if (!engine) return 1;
  if (!engine->active) {
    engine->active= 1;
    engine->nextActive= gistActive;
    gistActive= engine;
    engine->ChangeMap(engine);
  }
  return 0;
}

int GpDeactivate(Engine *engine)
{
  if (!engine) return 1;
  if (engine->active) {
    Engine *active= gistActive;
    engine->active= 0;
    if (active==engine) gistActive= engine->nextActive;
    else {
      while (active->nextActive!=engine) active= active->nextActive;
      active->nextActive= engine->nextActive;
    }
  }
  return 0;
}

int GpPreempt(Engine *engine)
{
  gistPreempt= engine;
  if (engine && !engine->active) engine->ChangeMap(engine);
  return 0;
}

int GpActive(Engine *engine)
{
  if (!engine) return 0;
  return engine==gistPreempt? 1 : engine->active;
}

int GpClear(Engine *engine, int flag)
{
  int value= 0;
  if (!engine) {
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine)) {
      engine->damaged= engine->inhibit= 0;
      engine->lastDrawn= -1;
      engine->systemsSeen[0]= engine->systemsSeen[1]= 0;
      value|= engine->Clear(engine, flag);
    }
  } else {
    engine->damaged= engine->inhibit= 0;
    engine->lastDrawn= -1;
    engine->systemsSeen[0]= engine->systemsSeen[1]= 0;
    value= engine->Clear(engine, flag);
  }
  return value;
}

int GpFlush(Engine *engine)
{
  if (!engine) {
    int value= 0;
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      value|= engine->Flush(engine);
    return value;
  }
  return engine->Flush(engine);
}

Engine *GpNextEngine(Engine *engine)
{
  return engine? engine->next : gistEngines;
}

Engine *GpNextActive(Engine *engine)
{
  if (gistPreempt) return engine? 0 : gistPreempt;
  else return engine? engine->nextActive : gistActive;
}

/* ------------------------------------------------------------------------ */

int GpSetTrans(const GpTransform *trans)
{
  Engine *engine;

  if (trans!=&gistT) gistT= *trans;

  for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
    engine->ChangeMap(engine);

  return 0;
}

int GpLandscape(Engine *engine, int landscape)
{
  if (!engine) {
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      engine->landscape= landscape;
  } else {
    engine->landscape= landscape;
  }
  return 0;
}

void GpSetMap(const GpBox *src, const GpBox *dst, GpXYMap *map)
{
  map->x.scale= (dst->xmax-dst->xmin)/(src->xmax-src->xmin);
  map->x.offset=  dst->xmin - map->x.scale*src->xmin;
  map->y.scale= (dst->ymax-dst->ymin)/(src->ymax-src->ymin);
  map->y.offset=  dst->ymin - map->y.scale*src->ymin;
}

void GpDeviceMap(Engine *engine)
{
  GpSetMap(&engine->transform.viewport, &engine->transform.window,
           &engine->devMap);
}

void GpComposeMap(Engine *engine)
{
  GpMap *devx= &engine->devMap.x;
  GpMap *devy= &engine->devMap.y;
  GpMap *mapx= &engine->map.x;
  GpMap *mapy= &engine->map.y;
  mapx->scale=
    devx->scale*(gistT.viewport.xmax-gistT.viewport.xmin)/
                (gistT.window.xmax-gistT.window.xmin);
  mapx->offset= devx->offset + devx->scale*gistT.viewport.xmin -
                mapx->scale*gistT.window.xmin;
  mapy->scale=
    devy->scale*(gistT.viewport.ymax-gistT.viewport.ymin)/
                (gistT.window.ymax-gistT.window.ymin);
  mapy->offset= devy->offset + devy->scale*gistT.viewport.ymin -
                mapy->scale*gistT.window.ymin;
}

/* ------------------------------------------------------------------------ */

/* Scratch space used by GpIntPoints and GpIntSegs */
static void *scratch= 0;
static long scratchPoints= 0, scratchSegs= 0;

static void MoreScratch(long np, long ns)
{
  if (scratch) p_free(scratch);
  if (np) {
    np+= 64;
    scratch= (void *)p_malloc(sizeof(GpPoint)*np);
    scratchPoints= np;
    scratchSegs= (sizeof(GpPoint)*np)/sizeof(GpSegment);
  } else {
    ns+= 32;
    scratch= (void *)p_malloc(sizeof(GpSegment)*ns);
    scratchSegs= ns;
    scratchPoints= (sizeof(GpSegment)*ns)/sizeof(GpPoint);
  }
}

long GpIntPoints(const GpXYMap *map, long maxPoints, long n,
                 const GpReal *x, const GpReal *y, GpPoint **result)
{
  GpReal scalx= map->x.scale, offx= map->x.offset;
  GpReal scaly= map->y.scale, offy= map->y.offset;
  long i, np= maxPoints<n? maxPoints : n;
  GpPoint *point;

  if (np+1>scratchPoints) MoreScratch(np+1, 0); /* allow for closure pt */
  *result= point= scratch;

  for (i=0 ; i<np ; i++) {
    point[i].x= (short)(scalx*x[i]+offx);
    point[i].y= (short)(scaly*y[i]+offy);
  }

  return np;
}

long GpIntSegs(const GpXYMap *map, long maxSegs, long n,
               const GpReal *x1, const GpReal *y1,
               const GpReal *x2, const GpReal *y2, GpSegment **result)
{
  GpReal scalx= map->x.scale, offx= map->x.offset;
  GpReal scaly= map->y.scale, offy= map->y.offset;
  long i, ns= maxSegs<n? maxSegs : n;
  GpSegment *seg;

  if (ns>scratchSegs) MoreScratch(0, ns);
  *result= seg= scratch;

  for (i=0 ; i<ns ; i++) {
    seg[i].x1= (short)(scalx*x1[i]+offx);
    seg[i].y1= (short)(scaly*y1[i]+offy);
    seg[i].x2= (short)(scalx*x2[i]+offx);
    seg[i].y2= (short)(scaly*y2[i]+offy);
  }

  return ns;
}

/* ------------------------------------------------------------------------ */

void GpPutGray(int nColors, GpColorCell *palette)
{
  /*
  while (nColors--) {
    palette->gray=
      ((int)palette->red+(int)palette->green+(int)palette->blue)/3;
    palette++;
  }
  */
}

void GpPutNTSC(int nColors, GpColorCell *palette)
{
  /*
  while (nColors--) {
    palette->gray=
      (30*(int)palette->red+59*(int)palette->green+11*(int)palette->blue)/100;
    palette++;
  }
  */
}

void GpPutRGB(int nColors, GpColorCell *palette)
{
  /*
  while (nColors--) {
    palette->red= palette->green= palette->blue= palette->gray;
    palette++;
  }
  */
}

int GpSetPalette(Engine *engine, GpColorCell *palette, int nColors)
{
  if (!engine) return 0;
  if (nColors<0) {
    palette= 0;
    nColors= 0;
  }
  engine->palette= palette;
  engine->nColors= nColors;
  engine->colorChange= 1;
  return engine->ChangePalette(engine);
}

int GpGetPalette(Engine *engine, GpColorCell **palette)
{
  *palette= engine? engine->palette : 0;
  return engine? engine->nColors : 0;
}

int GpDumpColors(Engine *engine, int colorMode)
{
  if (!engine) {
    for (engine=GpNextActive(0) ; engine ; engine=GpNextActive(engine))
      { engine->colorMode= colorMode;   engine->colorChange= 1; }
  } else {
    engine->colorMode= colorMode;   engine->colorChange= 1;
  }
  return 0;
}

/* ------------------------------------------------------------------------ */

long GpClipCells(GpMap *map, GpReal *px, GpReal *qx,
                 GpReal xmin, GpReal xmax, long ncells, long *off)
{
  long imin, imax;
  GpReal p, q, dx;
  GpReal scale= map->scale;
  GpReal offset= map->offset;

  xmin= xmin*scale+offset;
  xmax= xmax*scale+offset;
  if (xmin>xmax) {GpReal tmp=xmin; xmin=xmax; xmax=tmp;}
  p= (*px)*scale+offset;
  q= (*qx)*scale+offset;

  if (p<q && q>=xmin && p<=xmax) {
    dx= (q-p)/(GpReal)ncells;
    if (p<xmin) {
      imin= (long)((xmin-p)/dx);
      p+= dx*(GpReal)imin;
    } else {
      imin= 0;
    }
    if (q>xmax) {
      imax= (long)((q-xmax)/dx);
      q-= dx*(GpReal)imax;
      imax= ncells-imax;
    } else {
      imax= ncells;
    }
    if (imax-imin<2) {
      if (imax==imin) {
        if (p<xmin) p= xmin;
        if (q>xmax) q= xmax;
      } else {
        if (p<xmin && q>xmax) {
          if (q-xmax > xmin-p) { q-= xmin-p;  p= xmin; }
          else { p+= q-xmax;  q= xmax; }
        }
      }
    }
  } else if (p>q && p>=xmin && q<=xmax) {
    dx= (p-q)/(GpReal)ncells;
    if (q<xmin) {
      imax= (long)((xmin-q)/dx);
      q+= dx*(GpReal)imax;
      imax= ncells-imax;
    } else {
      imax= ncells;
    }
    if (p>xmax) {
      imin= (long)((p-xmax)/dx);
      p-= dx*(GpReal)imin;
    } else {
      imin= 0;
    }
    if (imax-imin<2) {
      if (imax==imin) {
        if (q<xmin) q= xmin;
        if (p>xmax) p= xmax;
      } else {
        if (q<xmin && p>xmax) {
          if (p-xmax > xmin-q) { p-= xmin-q;  q= xmin; }
          else { q+= p-xmax;  p= xmax; }
        }
      }
    }
  } else {
    imin= 0;
    imax= -1;
  }

  *px= p;
  *qx= q;
  *off= imin;

  return imax-imin;
}

/* ------------------------------------------------------------------------ */

int GpIntersect(const GpBox *box1, const GpBox *box2)
{
  /* Algorithm assumes min<max for x and y in both boxes */
  return box1->xmin<=box2->xmax && box1->xmax>=box2->xmin &&
         box1->ymin<=box2->ymax && box1->ymax>=box2->ymin;
}

int GpContains(const GpBox *box1, const GpBox *box2)
{
  /* Algorithm assumes min<max for x and y in both boxes */
  return box1->xmin<=box2->xmin && box1->xmax>=box2->xmax &&
         box1->ymin<=box2->ymin && box1->ymax>=box2->ymax;
}

void GpSwallow(GpBox *preditor, const GpBox *prey)
{
  /* Algorithm assumes min<max for x and y in both boxes */
  if (preditor->xmin>prey->xmin) preditor->xmin= prey->xmin;
  if (preditor->xmax<prey->xmax) preditor->xmax= prey->xmax;
  if (preditor->ymin>prey->ymin) preditor->ymin= prey->ymin;
  if (preditor->ymax<prey->ymax) preditor->ymax= prey->ymax;
}

/* ------------------------------------------------------------------------ */

/* These recondite routines are required to handle editing a drawing
   on one or more interactive engines.  The restriction to few
   routines builds in certain inefficiencies; if every drawing were always
   associated with one interactive engine some of the inefficiency could
   be reduced.  These are not intended for external use.  */

extern int gdNowRendering, gdMaxRendered;
int gdNowRendering= -1;
int gdMaxRendered= -1;

int GdBeginDr(Drauing *drawing, GpBox *damage, int landscape)
{
  int needToRedraw= 0;
  Engine *eng;

  if (damage) {
    /* If drawing has incurred damage, report damage to ALL engines
       interested in the drawing (not just active engines).  */
    for (eng=GpNextEngine(0) ; eng ; eng=GpNextEngine(eng))
      if (eng->drawing==drawing) GpDamage(eng, drawing, damage);
  }

  /* Loop on active engines to alert them that drawing is coming.  */
  for (eng=GpNextActive(0) ; eng ; eng=GpNextActive(eng)) {
    if (eng->drawing!=drawing) {
      /* This engine is not marked as interested in this drawing.
         Mark it, and reset damaged and lastDrawn flags so that no
         elements will be inhibited.  */
      eng->drawing= drawing;
      eng->lastDrawn= -1;
      eng->damaged= 0;
      if (landscape != eng->landscape) {
        eng->landscape= landscape;
        /* This change will be detected and acted upon by the first call
           to the ChangeMap method (GpSetTrans).  */
      }
      /* The semantics here are subtle --
         After a ClearDrawing, GdDetach zeroes eng->drawing in order to
         communicate that the drawing has been cleared.  Thus, the code
         gets here on a GdDraw after the drawing has been cleared, so
         the time has come to carry out the deferred clearing of this
         engine's plotting surface.  */
      GpClear(eng, CONDITIONALLY);
      needToRedraw= 1;

    } else if (eng->damaged) {
      /* This engine was interested in the drawing, which has been
         damaged.  Clear damaged area in preparation for repair work.
         (This is redundant if the damage was due to an X windows
          expose event, but the resulting inefficiency is very small.)  */
      eng->ClearArea(eng, &eng->damage);
      needToRedraw= 1;

    } else if (eng->lastDrawn<drawing->nElements-1) {
      needToRedraw= 1;
    }
  }

  gdNowRendering= gdMaxRendered= -1;
  return needToRedraw;
}

int GdBeginSy(GpBox *tickOut, GpBox *tickIn, GpBox *viewport,
              int number, int sysIndex)
{
  Engine *eng;
  int value= 0;
  long sysMask;

  /* Note that this is harmless if sysIndex>2*sizeof(long)--
     just slightly inefficient in that ticks and elements will ALWAYS
     be drawn...  This shouldn't be a practical problem.  */
  if (sysIndex>sizeof(long)) {
    sysMask= 1 << (sysIndex-sizeof(long));
    sysIndex= 1;
  } else {
    sysMask= 1 << sysIndex;
    sysIndex= 0;
  }

  /* Loop on active engines to determine whether any require ticks or
     elements to be drawn.  Set inhibit switches for ticks.  */
  for (eng=GpNextActive(0) ; eng ; eng=GpNextActive(eng)) {
    if ( ! (eng->systemsSeen[sysIndex] & sysMask) ) {
      /* this engine has never seen this system */
      value|= 3;
      eng->inhibit= 0;
      eng->systemsSeen[sysIndex]|= sysMask;

    } else if (eng->damaged && GpIntersect(tickOut, &eng->damage)) {
      /* engine damage touches this coordinate system--
         redraw ticks if region between tickIn and tickOut damaged,
         redraw elements if viewport damaged */
      if (!tickIn || !GpContains(tickIn, &eng->damage)) {
        value|= 2;
        eng->inhibit= 0;
      } else eng->inhibit= 1;
      if (number>eng->lastDrawn || GpIntersect(viewport, &eng->damage))
        value|= 1;

    } else {
      /* engine undamaged or damage doesn't touch this system--
         redraw elements if any new ones, don't redraw ticks */
      eng->inhibit= 1;
      if (number>eng->lastDrawn) value|= 1;
    }
  }

  return value;
}

int GdBeginEl(GpBox *box, int number)
{
  Engine *eng;
  int value= 0;

  /* Loop on active engines to determine whether any require this
     element to be drawn, and to set inhibit switches so that some
     may draw it and others not.  */
  for (eng=GpNextActive(0) ; eng ; eng=GpNextActive(eng)) {
    if (number>eng->lastDrawn) {
      /* this engine hasn't seen this element before */
      eng->inhibit= 0;
      value= 1;
      if (eng->damaged && gdMaxRendered<=eng->lastDrawn) {
        /* If this is the first new element, the damage flag
           must be reset, and ChangeMap must be called to set the
           clip rectangle back to its undamaged boundary.  */
        eng->damaged= 0;
        eng->ChangeMap(eng);
      }

    } else if (box && eng->damaged && GpIntersect(box, &eng->damage)) {
      /* engine damage touches this element */
      eng->inhibit= 0;
      value= 1;

    } else {
      /* this element has been seen before and hasn't been damaged */
      eng->inhibit= 1;
    }

    /* set number of element currently being drawn for GdEndDr */
    gdNowRendering= number;
    if (gdMaxRendered<gdNowRendering) gdMaxRendered= gdNowRendering;
  }

  return value;
}

void GdEndDr(void)
{
  Engine *eng;
  /* Done with this drawing- reset inhibit, damaged, and lastDrawn flags */
  for (eng=GpNextActive(0) ; eng ; eng=GpNextActive(eng)) {
    if (eng->lastDrawn<gdMaxRendered) eng->lastDrawn= gdMaxRendered;
    eng->inhibit= eng->damaged= 0;
  }
}

void GpDamage(Engine *eng, Drauing *drawing, GpBox *box)
{
  if (eng->drawing!=drawing || !eng->marked) return;
  if (eng->ClearArea==&DefaultClearArea) {
    /* This engine doesn't need to record the damage box */
    eng->damaged= 1;
  } else if (eng->damaged) {
    /* drawing is already damaged on this engine */
    if (eng->damage.xmin>box->xmin) eng->damage.xmin= box->xmin;
    if (eng->damage.xmax<box->xmax) eng->damage.xmax= box->xmax;
    if (eng->damage.ymin>box->ymin) eng->damage.ymin= box->ymin;
    if (eng->damage.ymax<box->ymax) eng->damage.ymax= box->ymax;
  } else {
    /* drawing is currently undamaged on this engine */
    eng->damaged= 1;
    eng->damage= *box;
  }
}

void GdDetach(Drauing *drawing, Engine *engine)
{
  Engine *eng;
  for (eng=GpNextEngine(0) ; eng ; eng=GpNextEngine(eng)) {
    if (!drawing || eng->drawing==drawing) {
      eng->drawing= 0;
      eng->inhibit= eng->damaged= 0;
      eng->lastDrawn= -1;
    }
  }
}

/* ------------------------------------------------------------------------ */
