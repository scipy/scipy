/*
   DRAW0.C
   Virtual functions for Drawing class.

   $Id$
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "draw.h"
#include "gtext.h"

extern char *strcpy(char *, const char *);

extern double log10(double);
#define SAFELOG0 (-999.)
#define SAFELOG(x) ((x)>0? log10(x) : ((x)<0? log10(-(x)) : -999.))

extern void Gd_KillRing(void *elv);

static void KillElement(GdElement *el);
static void LinesKill(void *el);
static void DisjointKill(void *el);
static void TextKill(void *el);
static void CellsKill(void *el);
static void PolysKill(void *el);
static void MeshKill(void *el);
static void FilledKill(void *el);
static void VectorsKill(void *el);
static void KillConGrps(GeLines **grps, int ngrps);
static void ContoursKill(void *el);
static void SystemKill(void *el);
extern void Gd_KillMeshXY(void *vMeshEl);

static int LinesGet(void *el);
static int DisjointGet(void *el);
static int TextGet(void *el);
static int CellsGet(void *el);
static int PolysGet(void *el);
static int MeshGet(void *el);
static int FilledGet(void *el);
static int VectorsGet(void *el);
static int ContoursGet(void *el);
static int SystemGet(void *el);
extern void Gd_MeshXYGet(void *vMeshEl);

static int LinesSet(void *el, int xyzChanged);
extern void Gd_LinesSubSet(void *el);
static int DisjointSet(void *el, int xyzChanged);
static int TextSet(void *el, int xyzChanged);
static int CellsSet(void *el, int xyzChanged);
static int PolysSet(void *el, int xyzChanged);
static void MeshXYSet(void *vMeshEl, int xyzChanged);
static int MeshSet(void *el, int xyzChanged);
static int FilledSet(void *el, int xyzChanged);
static int VectorsSet(void *el, int xyzChanged);
static int ContoursSet(void *el, int xyzChanged);
static int SystemSet(void *el, int xyzChanged);

static int LinesDraw(void *el, int xIsLog, int yIsLog);
static int DisjointDraw(void *el, int xIsLog, int yIsLog);
static int TextDraw(void *el, int xIsLog, int yIsLog);
static int CellsDraw(void *el, int xIsLog, int yIsLog);
static int PolysDraw(void *el, int xIsLog, int yIsLog);
static void MeshDrawSet(GaQuadMesh *mesh, void *vMeshEl,
			int xIsLog, int yIsLog);
static int MeshDraw(void *el, int xIsLog, int yIsLog);
static int FilledDraw(void *el, int xIsLog, int yIsLog);
static int VectorsDraw(void *el, int xIsLog, int yIsLog);
static int ContoursDraw(void *el, int xIsLog, int yIsLog);
static int SystemDraw(void *el, int xIsLog, int yIsLog);

static int LinesScan(void *el, int flags, GpBox *limits);
static int DisjointScan(void *el, int flags, GpBox *limits);
static int TextScan(void *el, int flags, GpBox *limits);
static int CellsScan(void *el, int flags, GpBox *limits);
static int PolysScan(void *el, int flags, GpBox *limits);
static int MeshXYScan(void *vMeshEl, int flags, GpBox *limits, GpBox *box);
static int MeshScan(void *el, int flags, GpBox *limits);
static int FilledScan(void *el, int flags, GpBox *limits);
static int VectorsScan(void *el, int flags, GpBox *limits);
static int ContoursScan(void *el, int flags, GpBox *limits);
static int SystemScan(void *el, int flags, GpBox *limits);

static void NoMargin(void *el, GpBox *margin);
static void LinesMargin(void *el, GpBox *margin);
static void DisjointMargin(void *el, GpBox *margin);
static void TextMargin(void *el, GpBox *margin);
static void MeshMargin(void *el, GpBox *margin);
static void VectorsMargin(void *el, GpBox *margin);
static void ContoursMargin(void *el, GpBox *margin);

extern void Gd_ScanZ(long n, const GpReal *z, GpReal *zmin, GpReal *zmax);
static int ScanMn(long n, GpReal *x, GpReal *y, GpReal ymin,
		  GpReal *xmin, GpReal *xmax);
static int ScanMx(long n, GpReal *x, GpReal *y, GpReal ymax,
		  GpReal *xmin, GpReal *xmax);
static int ScanMnMx(long n, GpReal *x, GpReal *y, GpReal ymin, GpReal ymax,
		    GpReal *xmin, GpReal *xmax);
static void ScanRXY(long n, GpReal *x, GpReal *y,
		    int flags, GpBox *limits, GpBox *box);
static int GetLogZ(long n, GpReal *z, GpReal **zlog,
		   GpReal *zmin, GpReal *zmax);
static int Get_LogZ(long n, long nndc, GpReal *z, GpReal **zlog,
		    GpReal *zmin, GpReal *zmax);

extern int Gd_MakeContours(GeContours *con);
extern int Gd_DrawRing(void *elv, int xIsLog, int yIsLog,
		       GeSystem *sys, int t);
extern void Gd_NextMeshBlock(long *ii, long *jj, long len, long iMax,
			     int *reg, int region);

/* ------------------------------------------------------------------------ */
/* Set virtual function tables */

extern GdOpTable *GetDrawingOpTables(void);
static GdOpTable opTables[E_SYSTEM+1]= {
  { E_NONE, 0, 0, 0, 0, 0, 0 },
  { E_LINES, &LinesKill, &LinesGet, &LinesSet,
      &LinesDraw, &LinesScan, &LinesMargin },
  { E_DISJOINT, &DisjointKill, &DisjointGet, &DisjointSet,
      &DisjointDraw, &DisjointScan, &DisjointMargin },
  { E_TEXT, &TextKill, &TextGet, &TextSet,
      &TextDraw, &TextScan, &TextMargin },
  { E_MESH, &MeshKill, &MeshGet, &MeshSet,
      &MeshDraw, &MeshScan, &MeshMargin },
  { E_FILLED, &FilledKill, &FilledGet, &FilledSet,
      &FilledDraw, &FilledScan, &NoMargin },
  { E_VECTORS, &VectorsKill, &VectorsGet, &VectorsSet,
      &VectorsDraw, &VectorsScan, &VectorsMargin },
  { E_CONTOURS, &ContoursKill, &ContoursGet, &ContoursSet,
      &ContoursDraw, &ContoursScan, &ContoursMargin },
  { E_CELLS, &CellsKill, &CellsGet, &CellsSet,
      &CellsDraw, &CellsScan, &NoMargin },
  { E_POLYS, &PolysKill, &PolysGet, &PolysSet,
      &PolysDraw, &PolysScan, &NoMargin },
  { E_SYSTEM, &SystemKill, &SystemGet, &SystemSet,
      &SystemDraw, &SystemScan, &NoMargin }
};

/* this is called at the first call to GdNewDrawing */
GdOpTable *GetDrawingOpTables(void)
{
  return opTables;
}

/* ------------------------------------------------------------------------ */
/* Destructors for drawing elements are private, accessed via the
   Kill virtual function */

void Gd_KillRing(void *elv)
{
  GdElement *el, *next= elv;
  while ((el= next)) {
    next= el->next;
    if (el == next) next= 0;
    el->ops->Kill(el);
  }
}

static void KillElement(GdElement *el)
{
  GdElement *next= el->next;
  if (el->legend) GmFree(el->legend);
  if (next) {
    if (next==el) next= 0;
    else { next->prev= el->prev;   el->prev->next= next; }
  }
  GmFree(el);
  return;
}

static void LinesKill(void *el)
{
  GeLines *lines= el;
  if (lines->x) GmFree(lines->x);
  if (lines->y) GmFree(lines->y);
  if (lines->xlog) GmFree(lines->xlog);
  if (lines->ylog) GmFree(lines->ylog);
  KillElement(el);
}

static void DisjointKill(void *el)
{
  GeDisjoint *lines= el;
  if (lines->x) GmFree(lines->x);
  if (lines->y) GmFree(lines->y);
  if (lines->xlog) GmFree(lines->xlog);
  if (lines->ylog) GmFree(lines->ylog);
  if (lines->xq) GmFree(lines->xq);
  if (lines->yq) GmFree(lines->yq);
  if (lines->xqlog) GmFree(lines->xqlog);
  if (lines->yqlog) GmFree(lines->yqlog);
  KillElement(el);
}

static void TextKill(void *el)
{
  GeText *text= el;
  if (text->text) GmFree(text->text);
  KillElement(el);
}

static void CellsKill(void *el)
{
  GeCells *cells= el;
  if (cells->colors) GmFree(cells->colors);
  KillElement(el);
}

static void PolysKill(void *el)
{
  GePolys *polys= el;
  if (polys->x) GmFree(polys->x);
  if (polys->y) GmFree(polys->y);
  if (polys->xlog) GmFree(polys->xlog);
  if (polys->ylog) GmFree(polys->ylog);
  if (polys->pn) GmFree(polys->pn);
  if (polys->colors) GmFree(polys->colors);
  KillElement(el);
}

static void MeshKill(void *el)
{
  Gd_KillMeshXY(el);
  KillElement(el);
}

static void FilledKill(void *el)
{
  GeFill *fill= el;
  Gd_KillMeshXY(el);
  if (fill->colors) {
    if (!(fill->noCopy&NOCOPY_COLORS)) GmFree(fill->colors);
    else if (GdFree) GdFree(fill->colors);
  }
  KillElement(el);
}

static void VectorsKill(void *el)
{
  GeVectors *vec= el;
  Gd_KillMeshXY(el);
  if (!(vec->noCopy&NOCOPY_UV)) {
    if (vec->u) GmFree(vec->u);
    if (vec->v) GmFree(vec->v);
  } else if (GdFree) {
    if (vec->u) GdFree(vec->u);
    if (vec->v) GdFree(vec->v);
  }
  KillElement(el);
}

static void KillConGrps(GeLines **grps, int ngrps)
{
  int i;
  for (i=0 ; i<ngrps ; i++) { Gd_KillRing(grps[i]);  grps[i]= 0; }
}

static void ContoursKill(void *el)
{
  GeContours *con= el;
  Gd_KillMeshXY(el);
  if (con->z) {
    if (!(con->noCopy&NOCOPY_Z)) GmFree(con->z);
    else if (GdFree) GdFree(con->z);
  }
  if (con->levels) GmFree(con->levels);
  if (con->groups) {
    KillConGrps(con->groups, con->nLevels);
    GmFree(con->groups);
  }
  KillElement(el);
}

static void SystemKill(void *el)
{
  GeSystem *sys= el;
  Gd_KillRing(sys->elements);
  KillElement(el);
}

void Gd_KillMeshXY(void *vMeshEl)
{
  GeMesh *meshEl= vMeshEl;
  GaQuadMesh *mesh= &meshEl->mesh;
  int noCopy= meshEl->noCopy;
  if (!(noCopy&NOCOPY_MESH)) {
    if (mesh->x) GmFree(mesh->x);
    if (mesh->y) GmFree(mesh->y);
  } else if (GdFree) {
    if (mesh->x) GdFree(mesh->x);
    if (mesh->y) GdFree(mesh->y);
  }
  if (mesh->reg) {
    if (!(noCopy&NOCOPY_REG)) GmFree(mesh->reg);
    else if (GdFree) GdFree(mesh->reg);
  }
  if (mesh->triangle) {
    if (!(noCopy&NOCOPY_TRI)) GmFree(mesh->triangle);
    else if (GdFree) GdFree(mesh->triangle);
  }
}

/* ------------------------------------------------------------------------ */
/* GetProps virtual function loads gistA, gistD from GdElement */

static int LinesGet(void *el)
{
  GeLines *e= el;
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.n= e->n;
  gistD.x= e->x;
  gistD.y= e->y;
  gistA.l= e->l;
  gistA.dl= e->dl;
  gistA.m= e->m;
  return E_LINES;
}

static int DisjointGet(void *el)
{
  GeDisjoint *e= el;
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.n= e->n;
  gistD.x= e->x;
  gistD.y= e->y;
  gistD.xq= e->xq;
  gistD.yq= e->yq;
  gistA.l= e->l;
  return E_DISJOINT;
}

static int TextGet(void *el)
{
  GeText *e= el;
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.x0= e->x0;
  gistD.y0= e->y0;
  gistD.text= e->text;
  gistA.t= e->t;
  return E_TEXT;
}

static int CellsGet(void *el)
{
  GeCells *e= el;
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.px= e->px;
  gistD.py= e->py;
  gistD.qx= e->qx;
  gistD.qy= e->qy;
  gistD.width= e->width;
  gistD.height= e->height;
  gistD.colors= e->colors;
  return E_CELLS;
}

static int PolysGet(void *el)
{
  GePolys *e= el;
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.n= e->n;
  gistD.x= e->x;
  gistD.y= e->y;
  gistD.pn= e->pn;
  gistD.colors= e->colors;
  return E_POLYS;
}

static int MeshGet(void *el)
{
  GeMesh *e= el;
  Gd_MeshXYGet(el);
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.boundary= e->boundary;
  gistD.inhibit= e->inhibit;
  gistA.l= e->l;
  return E_MESH;
}

static int FilledGet(void *el)
{
  GeFill *e= el;
  Gd_MeshXYGet(el);
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.nColumns= e->nColumns;
  gistD.colors= e->colors;
  gistA.e= e->e;
  return E_FILLED;
}

static int VectorsGet(void *el)
{
  GeVectors *e= el;
  Gd_MeshXYGet(el);
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.u= e->u;
  gistD.v= e->v;
  gistD.scale= e->scale;
  gistA.l= e->l;
  gistA.f= e->f;
  gistA.vect= e->vect;
  return E_VECTORS;
}

static int ContoursGet(void *el)
{
  GeContours *e= el;
  Gd_MeshXYGet(el);
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  gistD.z= e->z;
  gistD.nLevels= e->nLevels;
  gistD.levels= e->levels;
  gistA.l= e->l;
  gistA.dl= e->dl;
  gistA.m= e->m;
  return E_CONTOURS;
}

static int SystemGet(void *el)
{
  GeSystem *e= el;
  gistD.hidden= e->el.hidden;
  gistD.legend= e->el.legend;
  return E_SYSTEM;
}

void Gd_MeshXYGet(void *vMeshEl)
{
  GeMesh *meshEl= vMeshEl;
  GaQuadMesh *mesh= &meshEl->mesh;
  gistD.noCopy= meshEl->noCopy;
  gistD.mesh.iMax= mesh->iMax;
  gistD.mesh.jMax= mesh->jMax;
  gistD.mesh.x= mesh->x;
  gistD.mesh.y= mesh->y;
  gistD.mesh.reg= mesh->reg;
  gistD.mesh.triangle= mesh->triangle;
  gistD.region= meshEl->region;
}

/* ------------------------------------------------------------------------ */
/* SetProps virtual function loads GdElement from gistA, gistD */

static int LinesSet(void *el, int xyzChanged)
{
  GeLines *e= el;
  Gd_LinesSubSet(el);
  e->el.legend= gistD.legend;
  if (xyzChanged & CHANGE_XY) {
    e->n= gistD.n;
    e->x= gistD.x;
    e->y= gistD.y;
    if (e->xlog) { GmFree(e->xlog); e->xlog= 0; }
    if (e->ylog) { GmFree(e->ylog); e->ylog= 0; }
  }
  return 0;
}

void Gd_LinesSubSet(void *el)
{
  GeLines *e= el;
  e->el.hidden= gistD.hidden;
  e->l= gistA.l;
  e->dl= gistA.dl;
  e->m= gistA.m;
}

static int DisjointSet(void *el, int xyzChanged)
{
  GeDisjoint *e= el;
  e->el.hidden= gistD.hidden;
  e->el.legend= gistD.legend;
  e->n= gistD.n;
  e->x= gistD.x;
  e->y= gistD.y;
  e->xq= gistD.xq;
  e->yq= gistD.yq;
  e->l= gistA.l;
  if (xyzChanged & CHANGE_XY) {
    if (e->xlog) { GmFree(e->xlog); e->xlog= 0; }
    if (e->ylog) { GmFree(e->ylog); e->ylog= 0; }
    if (e->xqlog) { GmFree(e->xqlog); e->xqlog= 0; }
    if (e->yqlog) { GmFree(e->yqlog); e->yqlog= 0; }
  }
  return 0;
}

/* ARGSUSED */
static int TextSet(void *el, int xyzChanged)
{
  GeText *e= el;
  e->el.hidden= gistD.hidden;
  e->el.legend= gistD.legend;
  e->x0= gistD.x0;
  e->y0= gistD.y0;
  e->text= gistD.text;
  e->t= gistA.t;
  return 0;
}

/* ARGSUSED */
static int CellsSet(void *el, int xyzChanged)
{
  GeCells *e= el;
  e->el.hidden= gistD.hidden;
  e->el.legend= gistD.legend;
  e->px= gistD.px;
  e->py= gistD.py;
  e->qx= gistD.qx;
  e->qy= gistD.qy;
  e->width= gistD.width;
  e->height= gistD.height;
  e->colors= gistD.colors;
  return 0;
}

static int PolysSet(void *el, int xyzChanged)
{
  GePolys *e= el;
  Gd_LinesSubSet(el);
  e->el.legend= gistD.legend;
  if (xyzChanged & CHANGE_XY) {
    e->n= gistD.n;
    e->x= gistD.x;
    e->y= gistD.y;
    if (e->xlog) { GmFree(e->xlog); e->xlog= 0; }
    if (e->ylog) { GmFree(e->ylog); e->ylog= 0; }
  }
  e->pn= gistD.pn;
  e->colors= gistD.colors;
  return 0;
}

static void MeshXYSet(void *vMeshEl, int xyzChanged)
{
  GeMesh *meshEl= vMeshEl;
  GaQuadMesh *mesh= &meshEl->mesh;
  meshEl->el.legend= gistD.legend;
  meshEl->el.hidden= gistD.hidden;
  meshEl->noCopy= gistD.noCopy;
  mesh->iMax= gistD.mesh.iMax;
  mesh->jMax= gistD.mesh.jMax;
  mesh->x= gistD.mesh.x;
  mesh->y= gistD.mesh.y;
  mesh->reg= gistD.mesh.reg;
  mesh->triangle= gistD.mesh.triangle;
  if (xyzChanged & CHANGE_XY) {
    if (meshEl->xlog) { GmFree(meshEl->xlog); meshEl->xlog= 0; }
    if (meshEl->ylog) { GmFree(meshEl->ylog); meshEl->ylog= 0; }
  }
  meshEl->region= gistD.region;
}

static int MeshSet(void *el, int xyzChanged)
{
  GeMesh *e= el;
  MeshXYSet(el, xyzChanged);
  e->boundary= gistD.boundary;
  e->inhibit= gistD.inhibit;
  e->l= gistA.l;
  return 0;
}

static int FilledSet(void *el, int xyzChanged)
{
  GeFill *e= el;
  MeshXYSet(el, xyzChanged);
  e->nColumns= gistD.nColumns;
  e->colors= gistD.colors;
  e->e= gistA.e;
  return 0;
}

static int VectorsSet(void *el, int xyzChanged)
{
  GeVectors *e= el;
  MeshXYSet(el, xyzChanged);
  e->u= gistD.u;
  e->v= gistD.v;
  e->scale= gistD.scale;
  e->l= gistA.l;
  e->f= gistA.f;
  e->vect= gistA.vect;
  return 0;
}

static int ContoursSet(void *el, int xyzChanged)
{
  GeContours *e= el;
  int oldN= e->nLevels;
  MeshXYSet(el, xyzChanged);
  e->z= gistD.z;
  e->nLevels= gistD.nLevels;
  e->levels= gistD.levels;
  e->l= gistA.l;
  e->dl= gistA.dl;
  e->m= gistA.m;
  if (xyzChanged & CHANGE_Z) {
    if (e->groups) {
      KillConGrps(e->groups, oldN);
      if (oldN!=gistD.nLevels) {
	GmFree(e->groups);
	e->groups= 0;
      }
    }
    if (gistD.nLevels>0) {
      if (!e->groups)
	e->groups= (GeLines **)GmMalloc(sizeof(GeLines *)*gistD.nLevels);
      if (!e->groups || Gd_MakeContours(e)) return 1;
    }
  }
  return 0;
}

/* ARGSUSED */
static int SystemSet(void *el, int xyzChanged)
{
  GeSystem *e= el;
  e->el.hidden= gistD.hidden;
  e->el.legend= gistD.legend;
  return 0;
}

/* ------------------------------------------------------------------------ */
/* Draw virtual function calls Gp and Ga level rendering routines */

static int LinesDraw(void *el, int xIsLog, int yIsLog)
{
  GeLines *e= el;
  GpReal *px= xIsLog? e->xlog : e->x;
  GpReal *py= yIsLog? e->ylog : e->y;
  long n= e->n;
  if (e->el.hidden || n<=0) return 0;
  gistA.l= e->l;
  gistA.dl= e->dl;
  gistA.m= e->m;
  return GaLines(n, px, py);
}

static int DisjointDraw(void *el, int xIsLog, int yIsLog)
{
  GeDisjoint *e= el;
  GpReal *px= xIsLog? e->xlog : e->x;
  GpReal *py= yIsLog? e->ylog : e->y;
  GpReal *qx= xIsLog? e->xqlog : e->xq;
  GpReal *qy= yIsLog? e->yqlog : e->yq;
  long n= e->n;
  if (e->el.hidden || n<=0) return 0;
  gistA.l= e->l;
  return GpDisjoint(n, px, py, qx, qy);
}

static int TextDraw(void *el, int xIsLog, int yIsLog)
{
  GeText *e= el;
  GpReal x0, y0;
  if (e->el.hidden || !e->text) return 0;
  x0= xIsLog? SAFELOG(e->x0) : e->x0;
  y0= yIsLog? SAFELOG(e->y0) : e->y0;
  gistA.t= e->t;
  return GpText(x0, y0, e->text);
}

static int CellsDraw(void *el, int xIsLog, int yIsLog)
{
  GeCells *e= el;
  GpReal px, py, qx, qy;
  if (e->el.hidden) return 0;
  px= xIsLog? SAFELOG(e->px) : e->px;
  py= yIsLog? SAFELOG(e->py) : e->py;
  qx= xIsLog? SAFELOG(e->qx) : e->qx;
  qy= yIsLog? SAFELOG(e->qy) : e->qy;
  return GpCells(px, py, qx, qy, e->width, e->height, e->width, e->colors);
}

static int PolysDraw(void *el, int xIsLog, int yIsLog)
{
  int result= 0;
  GePolys *e= el;
  GpReal *px= xIsLog? e->xlog : e->x;
  GpReal *py= yIsLog? e->ylog : e->y;
  GpColor *colors= e->colors;
  long n= e->n;
  long *pn= e->pn;
  long i;
  if (e->el.hidden || n<=0) return 0;
  gistA.e= e->e;
  if (n<2 || pn[1]>1) {
    for (i=0 ; i<n ; i++) {
      gistA.f.color= colors? colors[i] : BG_COLOR;
      result|= GpFill(pn[i], px, py);
      px+= pn[i];
      py+= pn[i];
    }
  } else {
    long i0= pn[0]-1;
    for (i=1 ; i<n ; i++) {
      gistA.f.color= colors? colors[i] : BG_COLOR;
      if (gistA.f.color > 245) gistA.f.color-= 256;
      result|= GaFillMarker(pn[0], px, py, px[i0+i], py[i0+i]);
    }
  }
  return result;
}

static void MeshDrawSet(GaQuadMesh *mesh, void *vMeshEl,
			int xIsLog, int yIsLog)
{
  GeMesh *meshEl= vMeshEl;
  GaQuadMesh *msh= &meshEl->mesh;
  GpReal *x= xIsLog? meshEl->xlog : msh->x;
  GpReal *y= yIsLog? meshEl->ylog : msh->y;
  mesh->iMax= msh->iMax;
  mesh->jMax= msh->jMax;
  mesh->x= x;
  mesh->y= y;
  mesh->reg= msh->reg;
  mesh->triangle= msh->triangle;
}

static int MeshDraw(void *el, int xIsLog, int yIsLog)
{
  GeMesh *e= el;
  GaQuadMesh mesh;
  if (e->el.hidden) return 0;
  MeshDrawSet(&mesh, el, xIsLog, yIsLog);
  gistA.l= e->l;
  return GaMesh(&mesh, e->region, e->boundary, e->inhibit);
}

static int FilledDraw(void *el, int xIsLog, int yIsLog)
{
  GeFill *e= el;
  GaQuadMesh mesh;
  if (e->el.hidden) return 0;
  MeshDrawSet(&mesh, el, xIsLog, yIsLog);
  gistA.e= e->e;
  return GaFillMesh(&mesh, e->region, e->colors, e->nColumns);
}

static int VectorsDraw(void *el, int xIsLog, int yIsLog)
{
  GeVectors *e= el;
  GaQuadMesh mesh;
  if (e->el.hidden) return 0;
  MeshDrawSet(&mesh, el, xIsLog, yIsLog);
  gistA.l= e->l;
  gistA.f= e->f;
  gistA.vect= e->vect;
  return GaVectors(&mesh, e->region, e->u, e->v, e->scale);
}

static int ContoursDraw(void *el, int xIsLog, int yIsLog)
{
  GeContours *e= el;
  int nLevels= e->nLevels;
  GeLines **groups= e->groups;
  int value= 0;
  if (e->el.hidden || nLevels<=0) return 0;
  if (!groups) return 1;
  while (nLevels--) value|= Gd_DrawRing(*groups++, xIsLog, yIsLog, 0, 1);
  return value;
}

static int SystemDraw(void *el, int xIsLog, int yIsLog)
     /* NOTE: xIsLog input is used in a very non-standard way as the
	index of this system.  This is possible since xIsLog and yIsLog
	are otherwise meaningless in SystemDraw.  */
{
  GeSystem *e= el;
  int vflags, hflags= e->flags;
  int systemCounter= xIsLog;   /* Yes, this is non-standard usage */
  GpBox port, *tickIn;
  if (e->el.hidden || !e->elements) return 0;

  xIsLog= hflags & D_LOGX;
  yIsLog= hflags & D_LOGY;
  GpSetTrans(&e->trans);

  /* In order to prevent needless GaTick calls, a special feature
     is built into GdBeginSy.  */
  hflags= e->ticks.horiz.flags;
  vflags= e->ticks.vert.flags;
  if (vflags & TICK_C || hflags & TICK_C) tickIn= 0;
  else {
    GpReal tlen= e->ticks.vert.tickLen[0];
    GpReal twid= 0.5*e->ticks.vert.tickStyle.width*DEFAULT_LINE_WIDTH;
    tlen= (vflags&TICK_IN)? ((vflags&TICK_OUT)? 0.5 : 1.0)*tlen : 0.0;

    tickIn= &port;
    port= e->trans.viewport;
    if (vflags & TICK_L) port.xmin-= e->ticks.vert.tickOff + tlen + twid;
    if (vflags & TICK_U) port.xmax+= e->ticks.vert.tickOff - tlen - twid;

    tlen= e->ticks.horiz.tickLen[0];
    twid= 0.5*e->ticks.horiz.tickStyle.width*DEFAULT_LINE_WIDTH;
    tlen= (hflags&TICK_IN)? ((hflags&TICK_OUT)? 0.5 : 1.0)*tlen : 0.0;
    if (hflags & TICK_L) port.ymin-= e->ticks.horiz.tickOff + tlen + twid;
    if (hflags & TICK_U) port.ymax+= e->ticks.horiz.tickOff - tlen - twid;
  }

  hflags= GdBeginSy(&e->el.box, tickIn, &e->trans.viewport,
		    e->el.number, systemCounter);

  /* Draw the elements for this coordinate system before the ticks.  */
  gistClip= 1;   /* turn on clipping for elements */
  if (hflags & 1) Gd_DrawRing(e->elements, xIsLog, yIsLog, e, 0);

  /* Draw tick marks on top of elements.  If the user has chosen a style
     where the ticks overlap the viewport, he probably wants the ticks
     to obstruct his data anyway. */
  gistClip= 0;   /* turn off clipping for ticks */
  if (hflags & 2) GaAltTick(&e->ticks, xIsLog, yIsLog,
			    e->xtick, e->xlabel, e->ytick, e->ylabel);
  return 0;
}

/* ------------------------------------------------------------------------ */
/* Scan virtual function gets logs if reqd, sets box, scans xy  */

static int ScanMn(long n, GpReal *x, GpReal *y, GpReal ymin,
		  GpReal *xmin, GpReal *xmax)
{
  GpReal xn, xx;
  long i;
  for (i=0 ; i<n ; i++) if (y[i]>=ymin) break;
  if (i>=n) return 0;
  xn= xx= x[i++];
  for ( ; i<n ; i++) if (y[i]>=ymin) {
    if (x[i]<xn) xn= x[i];
    else if (x[i]>xx) xx= x[i];
  }
  *xmin= xn;
  *xmax= xx;
  return 1;
}

static int ScanMx(long n, GpReal *x, GpReal *y, GpReal ymax,
		  GpReal *xmin, GpReal *xmax)
{
  GpReal xn, xx;
  long i;
  for (i=0 ; i<n ; i++) if (y[i]<=ymax) break;
  if (i>=n) return 0;
  xn= xx= x[i++];
  for ( ; i<n ; i++) if (y[i]<=ymax) {
    if (x[i]<xn) xn= x[i];
    else if (x[i]>xx) xx= x[i];
  }
  *xmin= xn;
  *xmax= xx;
  return 1;
}

static int ScanMnMx(long n, GpReal *x, GpReal *y, GpReal ymin, GpReal ymax,
		    GpReal *xmin, GpReal *xmax)
{
  GpReal xn, xx;
  long i;
  for (i=0 ; i<n ; i++) if (y[i]>=ymin && y[i]<=ymax) break;
  if (i>=n) return 0;
  xn= xx= x[i++];
  for ( ; i<n ; i++) if (y[i]>=ymin && y[i]<=ymax) {
    if (x[i]<xn) xn= x[i];
    else if (x[i]>xx) xx= x[i];
  }
  *xmin= xn;
  *xmax= xx;
  return 1;
}

static void ScanRXY(long n, GpReal *x, GpReal *y,
		    int flags, GpBox *limits, GpBox *box)
{
  int dxmin= flags & D_XMIN,  dxmax= flags & D_XMAX;
  int dymin= flags & D_YMIN,  dymax= flags & D_YMAX;
  int any;

  if (dxmin || dxmax) {
    GpReal xmin, xmax;
    if (dymin) {
      if (dymax) { xmin= box->xmin;  xmax= box->xmax; any= 1; }
      else if (box->ymin>limits->ymax) any= 0;
      else any= ScanMx(n, x, y, limits->ymax, &xmin, &xmax);
    } else if (box->ymax<limits->ymin) {
      any= 0;
    } else {
      if (dymax) any= ScanMn(n, x, y, limits->ymin, &xmin, &xmax);
      else if (box->ymin>limits->ymax) any= 0;
      else any= ScanMnMx(n, x, y, limits->ymin, limits->ymax, &xmin, &xmax);
    }
    if (any) {
      if (dxmin) limits->xmin= xmin;
      if (dxmax) limits->xmax= xmax;
    } else {  /* GdScan requires min>max if no curves visible */
      if (dxmin) {
	if (dxmax) limits->xmax= 0.0;
	if (limits->xmax>0.0) limits->xmin= 1.1*limits->xmax;
	else limits->xmin= 0.9*limits->xmax+1.0;
      } else { /* dxmax is set */
	if (limits->xmin>0.0) limits->xmax= 0.9*limits->xmin;
	else limits->xmax= 1.1*limits->xmin-1.0;
      }
    }
  }
  if (dymin || dymax) {
    GpReal ymin, ymax;
    if (dxmin) {
      if (dxmax) { ymin= box->ymin;  ymax= box->ymax; any= 1; }
      else if (box->xmin>limits->xmax) any= 0;
      else any= ScanMx(n, y, x, limits->xmax, &ymin, &ymax);
    } else if (box->xmax<limits->xmin) {
      any= 0;
    } else {
      if (dxmax) any= ScanMn(n, y, x, limits->xmin, &ymin, &ymax);
      else if (box->xmin>limits->xmax) any= 0;
      else any= ScanMnMx(n, y, x, limits->xmin, limits->xmax, &ymin, &ymax);
    }
    if (any) {
      if (dymin) limits->ymin= ymin;
      if (dymax) limits->ymax= ymax;
    } else {  /* GdScan requires min>max if no curves visible */
      if (dymin) {
	if (dymax) limits->ymax= 0.0;
	if (limits->ymax>0.0) limits->ymin= 1.1*limits->ymax;
	else limits->ymin= 0.9*limits->ymax+1.0;
      } else { /* dymax is set */
	if (limits->ymin>0.0) limits->ymax= 0.9*limits->ymin;
	else limits->ymax= 1.1*limits->ymin-1.0;
      }
    }
  }
}

static int GetLogZ(long n, GpReal *z, GpReal **zlog,
		   GpReal *zmin, GpReal *zmax)
{
  GpReal *zl= (GpReal *)GmMalloc(sizeof(GpReal)*n);
  *zlog= zl;
  if (zl) {
    long i;
    for (i=0 ; i<n ; i++) zl[i]= SAFELOG(z[i]);
    if (zmin) Gd_ScanZ(n, zl, zmin, zmax);
  } else {
    strcpy(gistError, "memory manager failed in Gd log function");
    return -1;
  }
  return 0;
}

static int Get_LogZ(long n, long nndc, GpReal *z, GpReal **zlog,
		    GpReal *zmin, GpReal *zmax)
{
  GpReal *zl= (GpReal *)GmMalloc(sizeof(GpReal)*n);
  *zlog= zl;
  if (zl) {
    long i;
    for (i=0 ; i<nndc ; i++) zl[i]= z[i];
    for ( ; i<n ; i++) zl[i]= SAFELOG(z[i]);
    if (zmin) Gd_ScanZ(n-nndc, zl+nndc, zmin, zmax);
  } else {
    strcpy(gistError, "memory manager failed in Gd_log function");
    return -1;
  }
  return 0;
}

static int LinesScan(void *el, int flags, GpBox *limits)
{
  GeLines *e= el;
  GpReal *x, *y;

  /* First, get log values if necessary, and set box */
  if (flags & D_LOGX) {
    if (!e->xlog && GetLogZ(e->n, e->x, &e->xlog,
			    &e->logBox.xmin, &e->logBox.xmax)) return 1;
    e->el.box.xmin= e->logBox.xmin;
    e->el.box.xmax= e->logBox.xmax;
    x= e->xlog;
  } else {
    e->el.box.xmin= e->linBox.xmin;
    e->el.box.xmax= e->linBox.xmax;
    x= e->x;
  }
  if (flags & D_LOGY) {
    if (!e->ylog && GetLogZ(e->n, e->y, &e->ylog,
			    &e->logBox.ymin, &e->logBox.ymax)) return 1;
    e->el.box.ymin= e->logBox.ymin;
    e->el.box.ymax= e->logBox.ymax;
    y= e->ylog;
  } else {
    e->el.box.ymin= e->linBox.ymin;
    e->el.box.ymax= e->linBox.ymax;
    y= e->y;
  }

  if (flags & D_RESTRICT) {
    /* Scan points, restricting x limits to lie within fixed y limits
       and vice-versa.  Assume that limits.min<limits.max.  */
    ScanRXY(e->n, x, y, flags, limits, &e->el.box);
  } else {
    /* Unrestricted limits are either fixed or same as bounding box.  */
    if (flags & D_XMIN) limits->xmin= e->el.box.xmin;
    if (flags & D_XMAX) limits->xmax= e->el.box.xmax;
    if (flags & D_YMIN) limits->ymin= e->el.box.ymin;
    if (flags & D_YMAX) limits->ymax= e->el.box.ymax;
  }

  return 0;
}

static int DisjointScan(void *el, int flags, GpBox *limits)
{
  GeDisjoint *e= el;
  GpReal *x, *y, *xq, *yq;
  GpReal xymin, xymax;

  /* First, get log values if necessary, and set box */
  if (flags & D_LOGX) {
    if (!e->xlog && GetLogZ(e->n, e->x, &e->xlog,
			    &e->logBox.xmin, &e->logBox.xmax)) return 1;
    e->el.box.xmin= e->logBox.xmin;
    e->el.box.xmax= e->logBox.xmax;
    x= e->xlog;
    if (!e->xqlog && GetLogZ(e->n, e->xq, &e->xqlog,
			     &xymin, &xymax)) return 1;
    if (xymin<e->el.box.xmin) e->el.box.xmin= e->logBox.xmin;
    if (xymax>e->el.box.xmax) e->el.box.xmax= e->logBox.xmax;
    xq= e->xqlog;
  } else {
    e->el.box.xmin= e->linBox.xmin;
    e->el.box.xmax= e->linBox.xmax;
    x= e->x;
    xq= e->xq;
  }
  if (flags & D_LOGY) {
    if (!e->ylog && GetLogZ(e->n, e->y, &e->ylog,
			    &e->logBox.ymin, &e->logBox.ymax)) return 1;
    e->el.box.ymin= e->logBox.ymin;
    e->el.box.ymax= e->logBox.ymax;
    y= e->ylog;
    if (!e->yqlog && GetLogZ(e->n, e->yq, &e->yqlog,
			     &xymin, &xymax)) return 1;
    if (xymin<e->el.box.ymin) e->el.box.ymin= e->logBox.ymin;
    if (xymax>e->el.box.ymax) e->el.box.ymax= e->logBox.ymax;
    yq= e->yqlog;
  } else {
    e->el.box.ymin= e->linBox.ymin;
    e->el.box.ymax= e->linBox.ymax;
    y= e->y;
    yq= e->yq;
  }

  if (flags & D_RESTRICT) {
    /* Scan points, restricting x limits to lie within fixed y limits
       and vice-versa.  Assume that limits.min<limits.max.  */
    GpBox box;
    ScanRXY(e->n, x, y, flags, limits, &e->el.box);
    ScanRXY(e->n, xq, yq, flags, &box, &e->el.box);
    GpSwallow(limits, &box);
  } else {
    /* Unrestricted limits are either fixed or same as bounding box.  */
    if (flags & D_XMIN) limits->xmin= e->el.box.xmin;
    if (flags & D_XMAX) limits->xmax= e->el.box.xmax;
    if (flags & D_YMIN) limits->ymin= e->el.box.ymin;
    if (flags & D_YMAX) limits->ymax= e->el.box.ymax;
  }

  return 0;
}

static int TextScan(void *el, int flags, GpBox *limits)
{
  GeText *e= el;
  GpReal x0= e->x0;
  GpReal y0= e->y0;

  if (flags & D_LOGX) x0= SAFELOG(x0);
  if (flags & D_LOGY) y0= SAFELOG(y0);

  if (flags & D_XMIN) limits->xmin= x0;
  if (flags & D_XMAX) limits->xmax= x0;
  if (flags & D_YMIN) limits->ymin= y0;
  if (flags & D_YMAX) limits->ymax= y0;
  return 0;
}

static int CellsScan(void *el, int flags, GpBox *limits)
{
  GeCells *e= el;
  GpReal x[2], y[2];

  if (e->px<e->qx) { x[0]= e->px;  x[1]= e->qx; }
  else { x[0]= e->qx;  x[1]= e->px; }
  if (e->py<e->qy) { y[0]= e->py;  y[1]= e->qy; }
  else { y[0]= e->qy;  y[1]= e->py; }

  /* First, get log values if necessary, and set box */
  if (flags & D_LOGX) {
    e->el.box.xmin= SAFELOG(x[0]);
    e->el.box.xmax= SAFELOG(x[1]);
  } else {
    e->el.box.xmin= x[0];
    e->el.box.xmax= x[1];
  }
  if (flags & D_LOGY) {
    e->el.box.ymin= SAFELOG(y[0]);
    e->el.box.ymax= SAFELOG(y[1]);
  } else {
    e->el.box.ymin= y[0];
    e->el.box.ymax= y[1];
  }

  if (flags & D_XMIN) limits->xmin= e->el.box.xmin;
  if (flags & D_XMAX) limits->xmax= e->el.box.xmax;
  if (flags & D_YMIN) limits->ymin= e->el.box.ymin;
  if (flags & D_YMAX) limits->ymax= e->el.box.ymax;

  return 0;
}

static int PolysScan(void *el, int flags, GpBox *limits)
{
  GePolys *e= el;
  GpReal *x, *y;
  long i, ntot= 0;
  long nndc= (e->n<2 || e->pn[1]>1)? 0 : e->pn[0];

  /* compute total number of points */
  for (i=0 ; i<e->n ; i++) ntot+= e->pn[i];

  /* First, get log values if necessary, and set box */
  if (flags & D_LOGX) {
    if (!e->xlog && Get_LogZ(ntot, nndc, e->x, &e->xlog,
			     &e->logBox.xmin, &e->logBox.xmax)) return 1;
    e->el.box.xmin= e->logBox.xmin;
    e->el.box.xmax= e->logBox.xmax;
    x= e->xlog;
  } else {
    e->el.box.xmin= e->linBox.xmin;
    e->el.box.xmax= e->linBox.xmax;
    x= e->x;
  }
  if (flags & D_LOGY) {
    if (!e->ylog && Get_LogZ(ntot, nndc, e->y, &e->ylog,
			     &e->logBox.ymin, &e->logBox.ymax)) return 1;
    e->el.box.ymin= e->logBox.ymin;
    e->el.box.ymax= e->logBox.ymax;
    y= e->ylog;
  } else {
    e->el.box.ymin= e->linBox.ymin;
    e->el.box.ymax= e->linBox.ymax;
    y= e->y;
  }

  if (flags & D_RESTRICT) {
    /* Scan points, restricting x limits to lie within fixed y limits
       and vice-versa.  Assume that limits.min<limits.max.  */
    ScanRXY(ntot-nndc, x+nndc, y+nndc, flags, limits, &e->el.box);
  } else {
    /* Unrestricted limits are either fixed or same as bounding box.  */
    if (flags & D_XMIN) limits->xmin= e->el.box.xmin;
    if (flags & D_XMAX) limits->xmax= e->el.box.xmax;
    if (flags & D_YMIN) limits->ymin= e->el.box.ymin;
    if (flags & D_YMAX) limits->ymax= e->el.box.ymax;
  }

  return 0;
}

static int MeshXYScan(void *vMeshEl, int flags, GpBox *limits, GpBox *box)
{
  GeMesh *meshEl= vMeshEl;
  GaQuadMesh *mesh= &meshEl->mesh;
  GpReal *x, *y;

  /* First, get log values if necessary, and set box */
  if (flags & D_LOGX) {
    long len= mesh->iMax*mesh->jMax;
    int region= meshEl->region;
    GpReal xmin, xmax;
    long i, j, iMax= mesh->iMax;
    int *reg= mesh->reg, first= 1;

    if (!meshEl->xlog && GetLogZ(len, mesh->x, &meshEl->xlog, 0, 0))
      return 1;
    for (i=0 ; i<len ; ) {
      Gd_NextMeshBlock(&i, &j, len, iMax, reg, region);
      if (i>=len) break;
      Gd_ScanZ(j-i, meshEl->xlog+i, &xmin, &xmax);
      if (first) {
	meshEl->logBox.xmin= xmin;
	meshEl->logBox.xmax= xmax;
	first= 0;
      } else {
	if (xmin<meshEl->logBox.xmin) meshEl->logBox.xmin= xmin;
	if (xmax>meshEl->logBox.xmax) meshEl->logBox.xmax= xmax;
      }
      i= j+1;
    }
    box->xmin= meshEl->logBox.xmin;
    box->xmax= meshEl->logBox.xmax;
    x= meshEl->xlog;
  } else {
    box->xmin= meshEl->linBox.xmin;
    box->xmax= meshEl->linBox.xmax;
    x= mesh->x;
  }
  if (flags & D_LOGY) {
    long len= mesh->iMax*mesh->jMax;
    int region= meshEl->region;
    GpReal ymin, ymax;
    long i, j, iMax= mesh->iMax;
    int *reg= mesh->reg, first= 1;

    if (!meshEl->ylog && GetLogZ(len, mesh->y, &meshEl->ylog, 0, 0))
      return 1;
    for (i=0 ; i<len ; ) {
      Gd_NextMeshBlock(&i, &j, len, iMax, reg, region);
      if (i>=len) break;
      Gd_ScanZ(j-i, meshEl->ylog+i, &ymin, &ymax);
      if (first) {
	meshEl->logBox.ymin= ymin;
	meshEl->logBox.ymax= ymax;
	first= 0;
      } else {
	if (ymin<meshEl->logBox.ymin) meshEl->logBox.ymin= ymin;
	if (ymax>meshEl->logBox.ymax) meshEl->logBox.ymax= ymax;
      }
      i= j+1;
    }
    box->ymin= meshEl->logBox.ymin;
    box->ymax= meshEl->logBox.ymax;
    y= meshEl->ylog;
  } else {
    box->ymin= meshEl->linBox.ymin;
    box->ymax= meshEl->linBox.ymax;
    y= mesh->y;
  }

  if (flags & D_RESTRICT) {
    /* Scan points, restricting x limits to lie within fixed y limits
       and vice-versa.  Assume that limits.min<limits.max.  */
    long len= mesh->iMax*mesh->jMax;
    int region= meshEl->region;
    GpBox tbox;
    long i, j, iMax= mesh->iMax;
    int *reg= mesh->reg, first= 1;
    tbox= *limits;
    for (i=0 ; i<len ; ) {
      Gd_NextMeshBlock(&i, &j, len, iMax, reg, region);
      if (i>=len) break;
      ScanRXY(j-i, x+i, y+i, flags, limits, &tbox);
      if (first) { *box= tbox;  first= 0; }
      else GpSwallow(box, &tbox);
      i= j+1;
    }
  } else {
    /* Unrestricted limits are either fixed or same as bounding box.  */
    if (flags & D_XMIN) limits->xmin= box->xmin;
    if (flags & D_XMAX) limits->xmax= box->xmax;
    if (flags & D_YMIN) limits->ymin= box->ymin;
    if (flags & D_YMAX) limits->ymax= box->ymax;
  }

  return 0;
}

static int MeshScan(void *el, int flags, GpBox *limits)
{
  GeMesh *e= el;
  return MeshXYScan(el, flags, limits, &e->el.box);
}

static int FilledScan(void *el, int flags, GpBox *limits)
{
  GeFill *e= el;
  return MeshXYScan(el, flags, limits, &e->el.box);
}

static int VectorsScan(void *el, int flags, GpBox *limits)
{
  GeVectors *e= el;
  return MeshXYScan(el, flags, limits, &e->el.box);
}

static int ContoursScan(void *el, int flags, GpBox *limits)
{
  GeContours *e= el;
  GpBox lims= *limits;
  GeLines *elx, *el0, **groups= e->groups;
  int i, value= 0, none= 1;
  for (i=0 ; i<e->nLevels ; i++) {
    el0= *groups++;
    if ((elx= el0)) do {
      value|= LinesScan(elx, flags, &lims);
      if (none) { *limits= lims;   e->el.box= lims; }
      else { GpSwallow(limits, &lims);   GpSwallow(&e->el.box, &lims); }
      none= 0;
      elx= (GeLines *)elx->el.next;
    } while (elx != el0);
  }
  if (none) value= MeshXYScan(el, flags, limits, &e->el.box);
  return value;
}

/* ARGSUSED */
static int SystemScan(void *el, int flags, GpBox *limits)
{
  return 0;   /* cannot ever happen... */
}

/* ------------------------------------------------------------------------ */
/* Margin virtual function returns margin box */

/* ARGSUSED */
static void NoMargin(void *el, GpBox *margin)
{
  margin->xmin= margin->xmax= margin->ymin= margin->ymax= 0.0;
}

static void LinesMargin(void *el, GpBox *margin)
{
  /* This only accounts for line width, ignoring other decorations--
     other decorations seem unlikely outside coordinate systems */
  GeLines *lines= el;
  margin->xmin= margin->xmax= margin->ymin= margin->ymax=
    0.5*lines->l.width*DEFAULT_LINE_WIDTH;
}

static void DisjointMargin(void *el, GpBox *margin)
{
  GeDisjoint *lines= el;
  margin->xmin= margin->xmax= margin->ymin= margin->ymax=
    0.5*lines->l.width*DEFAULT_LINE_WIDTH;
}

static void TextMargin(void *el, GpBox *margin)
{
  /* The actual text box cannot be computed without text metric data--
     the following is a crude guess based on character counts and an
     assumed width/height ratio of 0.6 (as in 9x15) and descent/height
     ratio of 0.2.  This should be close for Courier, but it probably
     way off in width for the proportional fonts.  */
  GeText *text= el;
  GpReal width, x0, y0, dx, dy;
  int alignH, alignV;
  int nLines= GtTextShape(text->text, &text->t, (WidthFunction)0, &width);

  dx= text->t.height*width*0.6;
  dy= text->t.height*((GpReal)nLines);

  GtGetAlignment(&text->t, &alignH, &alignV);
  if (alignH==TH_LEFT) {
    x0= 0.0;
  } else if (alignH==TH_CENTER) {
    x0= -0.5*dx;
  } else {
    x0= -dx;
  }
  if (alignV==TV_TOP || alignV==TV_CAP) {
    y0= -dy;
  } else if (alignH==TV_HALF) {
    y0= -0.1*text->t.height - 0.5*dy;
  } else if (alignH==TV_BASE) {
    y0= -0.2*text->t.height;
  } else {
    y0= 0.0;
  }

  margin->xmin= x0;
  margin->xmax= x0 + dx;
  margin->ymin= y0;
  margin->ymax= y0 + dy;
}

static void MeshMargin(void *el, GpBox *margin)
{
  GeMesh *mesh= el;
  margin->xmin= margin->xmax= margin->ymin= margin->ymax=
    0.5*mesh->l.width*DEFAULT_LINE_WIDTH;
}

/* ARGSUSED */
static void VectorsMargin(void *el, GpBox *margin)
{
  /* This is a wild guess-- otherwise must scan (u, v) --
     should never arise in practice */
  /* GeVectors *vec= el; */
  margin->xmin= margin->xmax= margin->ymin= margin->ymax= 0.05;
}

static void ContoursMargin(void *el, GpBox *margin)
{
  /* Should never actually happen */
  GeContours *con= el;
  margin->xmin= margin->xmax= margin->ymin= margin->ymax=
    0.5*con->l.width*DEFAULT_LINE_WIDTH;
}

/* ------------------------------------------------------------------------ */
