/*
 * PS.C
 *
 * $Id$
 *
 * Implement the PostScript engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "pstdio.h"
#include "ps.h"
#include "gtext.h"

extern p_file *GistOpen(const char *name);  /* from gread.c */

#include "play.h"
#include <string.h>
#include <time.h>

#ifndef CREATE_PS
#define CREATE_PS(name) p_fopen(name, "w")
#endif

static char *psType= "PostScript";

static char line[80];  /* no lines longer than 78 characters! */

#define PS_PER_POINT 20
#define NDC_TO_PS (20.0/ONE_POINT)
#define DEFAULT_PS_WIDTH (DEFAULT_LINE_WIDTH*NDC_TO_PS)
#define N_PSFONTS 20
#define N_PSDASHES 5

static int PutPrologLine(p_file *file);
static p_file *CopyProlog(const char *name, const char *title);
static void InitBB(GpsBBox *bb);
static void SetPageDefaults(PSEngine *psEngine);
static void SetPSTransform(GpTransform *toPixels, int landscape);
static void GetEPSFBox(int landscape, GpsBBox *bb,
                  int *xll, int *yll, int *xur, int *yur);
static int PutLine(PSEngine *psEngine);
static int Append(PSEngine *psEngine, const char *s);
static int BeginPage(PSEngine *psEngine);
static int EndClip(PSEngine *psEngine);
static int BeginClip(PSEngine *psEngine, GpTransform *trans);
static int EndPage(PSEngine *psEngine);
static int PutPoints(PSEngine *psEngine, GpPoint *points, long nPoints,
                     int margin);
static int ChangePalette(Engine *engine);
static void Kill(Engine *engine);
static int Clear(Engine *engine, int always);
static int Flush(Engine *engine);
static int SetupColor(PSEngine *psEngine, unsigned long color);
static int SetupLine(PSEngine *psEngine, GpLineAttribs *gistAl);
static int CheckClip(PSEngine *psEngine);
static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth);
static int SetupFont(PSEngine *psEngine, GpReal height);
static int DrawMarkers(Engine *engine, long n, const GpReal *px,
                       const GpReal *py);
static int SetupText(PSEngine *psEngine);
static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text);
static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py);
static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors);
static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy);

static int ps_fputs(p_file *file, char *buf);

/* ------------------------------------------------------------------------ */

static const char *titleIs;
static int needUser, needDate;
static void *pf_stdout = &titleIs;  /* any non-zero address will do */

static int
ps_fputs(p_file *file, char *buf)
{
  if (file != pf_stdout) return p_fputs(file, buf);
  if (g_stdout) g_stdout(buf);
  return 0;
}

static int PutPrologLine(p_file *file)
{
  if (titleIs && strncmp(line, "%%Title:", 8L)==0) {
    line[8]= ' ';
    line[9]= '\0';
    strncat(line, titleIs, 60L);
    strcat(line, "\n");
    titleIs= 0;
  } else if (needUser && strncmp(line, "%%For:", 6L)==0) {
    line[6]= ' ';
    line[7]= '\0';
    strncat(line, p_getuser(), 60L);
    strcat(line, "\n");
    needUser= 0;
  } else if (needDate && strncmp(line, "%%CreationDate:", 15L)==0) {
    time_t t= time((time_t *)0);
    /* ctime returns 26 chars, e.g.- "Sun Jan  3 15:14:13 1988\n\0" */
    char *st= (t==-1)? "\n" : ctime(&t);
    line[15]= ' ';
    line[16]= '\0';
    strcat(line, st? st : "\n");
    needDate= 0;
  }
  return ps_fputs(file, line);
}

static p_file *CopyProlog(const char *name, const char *title)
{
  p_file *psps= GistOpen("ps.ps");
  p_file *file= strcmp(name, "*stdout*")? CREATE_PS(name) : pf_stdout;
  if (!psps) strcpy(gistError, "unable to open PostScript prolog ps.ps");
  if (!file) strcpy(gistError, "unable to create PostScript output file");

  if (psps && file) {
    titleIs= title;
    needUser= needDate= 1;
    for (;;) {
      if (!p_fgets(psps, line, 79) || PutPrologLine(file)<0) {
        if (file!=pf_stdout) p_fclose(file);
        file= 0;
        strcpy(gistError, "bad PostScript prolog format in ps.ps??");
        break;
      }
      if (strncmp(line, "%%EndSetup", 10L)==0) break;
    }

  } else if (file) {
    if (file!=pf_stdout) p_fclose(file);
    file= 0;
  }

  if (psps) p_fclose(psps);
  return file;
}

static void InitBB(GpsBBox *bb)
{
  bb->xll= bb->yll= 0x7ff0;
  bb->xur= bb->yur= 0;
}

static void SetPageDefaults(PSEngine *psEngine)
{
  /* Set current state to match state set by GI procedure in ps.ps */
  psEngine->curClip= 0;
  psEngine->curColor= P_FG;
  psEngine->curType= L_SOLID;
  psEngine->curWidth= 1.0;
  psEngine->curFont= T_COURIER;
  psEngine->curHeight= 12.0*ONE_POINT;
  psEngine->curAlignH= TH_LEFT;
  psEngine->curAlignV= TV_BASE;
  psEngine->curOpaque= 0;
  InitBB(&psEngine->pageBB);
}

static void SetPSTransform(GpTransform *toPixels, int landscape)
{
  /* PostScript thinks an 8.5-by-11 inch page is 612-by-792 points */
  toPixels->window.xmin= toPixels->window.ymin= 0.0;
  if (landscape) {
    toPixels->window.xmax= 15840.0;
    toPixels->window.ymax= 12240.0;
  } else {
    toPixels->window.xmax= 12240.0;
    toPixels->window.ymax= 15840.0;
  }
  toPixels->viewport.xmin= toPixels->viewport.ymin= 0.0;
  toPixels->viewport.xmax= toPixels->window.xmax*(1.0/NDC_TO_PS);
  toPixels->viewport.ymax= toPixels->window.ymax*(1.0/NDC_TO_PS);
}

static void GetEPSFBox(int landscape, GpsBBox *bb,
                       int *xll, int *yll, int *xur, int *yur)
{
  /* Transform bounding box to
     Document Structure Comment coordinates for use as EPSF file */
  int xl, yl, xu, yu;

  if (bb->xll < bb->xur) {
    /* this is a valid bounding box */
    xl= bb->xll/PS_PER_POINT;
    yl= bb->yll/PS_PER_POINT;
    xu= 1+(bb->xur-1)/PS_PER_POINT;
    yu= 1+(bb->yur-1)/PS_PER_POINT;

  } else {
    /* this is not a valid bounding box, return whole page */
    xl= yl= 0;
    xu= 612;    /* In PostScript, 1 inch is exactly 72 points */
    yu= 792;
  }

  if (landscape) {
    *xll= 612-yu;  *yll= xl;
    *xur= 612-yl;  *yur= xu;
  } else {
    *xll= xl;  *yll= yl;
    *xur= xu;  *yur= yu;
  }
}

static int PutLine(PSEngine *psEngine)
{
  p_file *file= psEngine->file;
  char *line= psEngine->line;
  long nchars= psEngine->nchars;

  if (!file) {
    if (psEngine->closed) return 1;
    file= psEngine->file= CopyProlog(psEngine->filename, psEngine->e.name);
    if (!file) { psEngine->closed= 1;  return 1; }
    psEngine->currentPage= 1;
    SetPageDefaults(psEngine);
    InitBB(&psEngine->docBB);
  }

  line[nchars++]= '\n';
  line[nchars]= '\0';
  if (ps_fputs(file, line)<0) {
    if (file!=pf_stdout) p_fclose(file);
    psEngine->file= 0;
    psEngine->closed= 1;
    strcpy(gistError, "p_fputs failed writing PostScript file");
    return 1;
  }
  line[0]= '\0';
  psEngine->nchars= 0;
  return 0;
}

static int Append(PSEngine *psEngine, const char *s)
{
  long len= s? strlen(s) : 0;
  long nchars= psEngine->nchars;
  char *line= psEngine->line;

  if (nchars+1+len>78) {
    if (PutLine(psEngine)) return 1;
    nchars= 0;
  } else if (nchars>0) {
    line[nchars++]= ' ';
  }
  strcpy(line+nchars, s);
  psEngine->nchars= nchars+len;

  return 0;
}

static int BeginPage(PSEngine *psEngine)
{
  int currentPage= psEngine->currentPage;

  psEngine->e.marked= 1;

  /* A change in color table can take place only at the beginning
     of a page.  ChangePalette strobes the palette in the Engine base
     class part into the PSEngine palette.  */
  psEngine->nColors= 0;      /* reset to mono mode */
  ChangePalette((Engine *)psEngine);

  if (psEngine->nchars && PutLine(psEngine)) return 1;
  if (PutLine(psEngine)) return 1;

  sprintf(line, "%%%%Page: %d %d", currentPage, currentPage);
  if (Append(psEngine, line) || PutLine(psEngine)) return 1;

  if (Append(psEngine, "%%PageBoundingBox: (atend)") ||
      PutLine(psEngine)) return 1;

  if (Append(psEngine, "GistPrimitives begin /PG save def GI") ||
      PutLine(psEngine)) return 1;

  /* Set transform viewport to reflect current page orientation */
  if (psEngine->landscape != psEngine->e.landscape) {
    SetPSTransform(&psEngine->e.transform, psEngine->e.landscape);
    psEngine->landscape= psEngine->e.landscape;
  }
  if (psEngine->landscape) {
    if (Append(psEngine, "LAND") ||
        PutLine(psEngine)) return 1;
  }

  if (psEngine->e.colorMode && psEngine->e.palette &&
      psEngine->e.nColors>0) {
    int i, nColors= psEngine->e.nColors;
    GpColorCell *palette= psEngine->e.palette;
    long color;
    sprintf(line, "%d CT", nColors);
    if (Append(psEngine, line) || PutLine(psEngine)) return 1;
    for (i=0 ; i<nColors ; i++) {
      color= P_R(palette[i])<<16 | P_G(palette[i])<<8 | P_B(palette[i]);
      sprintf(line, "%06lx", color);
      if (Append(psEngine, line)) return 1;
    }
    if (psEngine->nchars && PutLine(psEngine)) return 1;
    psEngine->colorMode= 1;  /* color table has been written */
    psEngine->nColors= nColors;
  } else {
    psEngine->colorMode= 0;  /* NO color table exists on this page */
    /* But, if there is a palette, still want to use grays */
    if (psEngine->e.palette && psEngine->e.nColors>0)
      psEngine->nColors= psEngine->e.nColors;
  }

  if (Append(psEngine, "%%EndPageSetup") || PutLine(psEngine)) return 1;
  return 0;
}

static int EndClip(PSEngine *psEngine)
{
  if (psEngine->curClip) {
    if ((psEngine->nchars && PutLine(psEngine)) ||
        Append(psEngine, "CLOF")) return 1;
    psEngine->curClip= 0;
    /* Restore state at time of CLON */
    psEngine->curColor= psEngine->clipColor;
    psEngine->curType= psEngine->clipType;
    psEngine->curWidth= psEngine->clipWidth;
    psEngine->curFont= psEngine->clipFont;
    psEngine->curHeight= psEngine->clipHeight;
  }
  return 0;
}

static int BeginClip(PSEngine *psEngine, GpTransform *trans)
{
  GpReal x[2], y[2];
  GpPoint *points;
  int xll, yll, xur, yur;
  GpBox *port= &trans->viewport;
  GpBox *box= &psEngine->clipBox;

  if (!psEngine->e.marked && BeginPage(psEngine)) return 1;

  if (psEngine->curClip) {
    if (port->xmin==box->xmin && port->xmax==box->xmax &&
        port->ymin==box->ymin && port->ymax==box->ymax) return 0;
    if (EndClip(psEngine)) return 1;
  }

  x[0]= trans->window.xmin;  x[1]= trans->window.xmax;
  y[0]= trans->window.ymin;  y[1]= trans->window.ymax;

  GpIntPoints(&psEngine->e.map, 3, 2, x, y, &points);
  if (points[0].x > points[1].x) {
    xll= points[1].x;  xur= points[0].x;
  } else {
    xll= points[0].x;  xur= points[1].x;
  }
  if (points[0].y > points[1].y) {
    yll= points[1].y;  yur= points[0].y;
  } else {
    yll= points[0].y;  yur= points[1].y;
  }

  sprintf(line, "%d %d %d %d CLON", xur-xll, yur-yll, xll, yll);
  if (Append(psEngine, line)) return 1;
  psEngine->curClip= 1;
  *box= *port;

  /* Must save state at time of CLON, since CLOF does grestore */
  psEngine->clipColor= psEngine->curColor;
  psEngine->clipType= psEngine->curType;
  psEngine->clipWidth= psEngine->curWidth;
  psEngine->clipFont= psEngine->curFont;
  psEngine->clipHeight= psEngine->curHeight;
  /* Note that text alignment/opacity is not affected by grestore */

  /* Expand page boundary to include clip boundary */
  if (xll < psEngine->pageBB.xll) psEngine->pageBB.xll= xll;
  if (yll < psEngine->pageBB.yll) psEngine->pageBB.yll= yll;
  if (xur > psEngine->pageBB.xur) psEngine->pageBB.xur= xur;
  if (yur > psEngine->pageBB.yur) psEngine->pageBB.yur= yur;
  return 0;
}

static int EndPage(PSEngine *psEngine)
{
  char *line= psEngine->line;
  int xll, yll, xur, yur;

  if (EndClip(psEngine)) return 1;

  if ((psEngine->nchars && PutLine(psEngine)) ||
      Append(psEngine, "PG restore") || PutLine(psEngine)) return 1;
  if (Append(psEngine, "showpage") || PutLine(psEngine)) return 1;
  if (Append(psEngine, "end") || PutLine(psEngine)) return 1;
  if (Append(psEngine, "%%PageTrailer") || PutLine(psEngine)) return 1;

  GetEPSFBox(psEngine->e.landscape, &psEngine->pageBB,
             &xll, &yll, &xur, &yur);
  if (xll < psEngine->docBB.xll) psEngine->docBB.xll= xll;
  if (yll < psEngine->docBB.yll) psEngine->docBB.yll= yll;
  if (xur > psEngine->docBB.xur) psEngine->docBB.xur= xur;
  if (yur > psEngine->docBB.yur) psEngine->docBB.yur= yur;
  sprintf(line, "%%%%PageBoundingBox: %d %d %d %d", xll, yll, xur, yur);
  if (Append(psEngine, line) || PutLine(psEngine)) return 1;

  psEngine->currentPage++;
  psEngine->e.marked= 0;
  SetPageDefaults(psEngine);
  if (psEngine->file!=pf_stdout) p_fflush(psEngine->file);
  return 0;
}

static char hexChar[17]= "0123456789abcdef";

static int PutPoints(PSEngine *psEngine, GpPoint *points, long nPoints,
                     int margin)
{
  int ix, iy, i;
  int xll= 0x7ff0, yll= 0x7ff0, xur= 0, yur= 0;
  char *now, *line= psEngine->line;
  int perLine, nchars= psEngine->nchars;

  if (!psEngine->e.marked && BeginPage(psEngine)) return 1;

  if (nchars<=0) {
    now= line;
    perLine= 9;
  } else if (nchars<70) {
    now= line+nchars;
    perLine= (78-nchars)/8;
  } else {
    if (PutLine(psEngine)) return 1;
    nchars= 0;
    now= line;
    perLine= 9;
  }

  /* Dump the points in hex 9 (72 chars) to a line */
  while (nPoints>0) {
    for (i=0 ; i<perLine ; i++) {
      ix= points->x;
      iy= points->y;
      points++;
      if (ix<0) ix= xll;           /* be sure x in range, set xll, xur */
      else if (ix>0x7ff0) ix= xur;
      else if (ix<xll) xll= ix;
      else if (ix>xur) xur= ix;
      if (iy<0) iy= yll;           /* be sure y in range, set yll, yur */
      else if (iy>0x7ff0) iy= yur;
      else if (iy<yll) yll= iy;
      else if (iy>yur) yur= iy;
      *now++= hexChar[ix>>12];
      *now++= hexChar[(ix>>8)&0xf];
      *now++= hexChar[(ix>>4)&0xf];
      *now++= hexChar[ix&0xf];
      *now++= hexChar[iy>>12];
      *now++= hexChar[(iy>>8)&0xf];
      *now++= hexChar[(iy>>4)&0xf];
      *now++= hexChar[iy&0xf];
      nchars+= 8;
      nPoints--;
      if (nPoints==0) break;
    }
    psEngine->nchars= nchars;
    if (PutLine(psEngine)) return 1;
    nchars= 0;
    now= line;
    perLine= 9;
  }

  /* Adjust the bounding box for the current page */
  if (xll<xur && !psEngine->curClip) {
    xll-= margin;
    xur+= margin;
    yll-= margin;
    yur+= margin;
    if (xll<psEngine->pageBB.xll) psEngine->pageBB.xll= xll;
    if (xur>psEngine->pageBB.xur) psEngine->pageBB.xur= xur;
    if (yll<psEngine->pageBB.yll) psEngine->pageBB.yll= yll;
    if (yur>psEngine->pageBB.yur) psEngine->pageBB.yur= yur;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int ChangePalette(Engine *engine)
{
  PSEngine *psEngine= (PSEngine *)engine;
  int nColors= engine->nColors;

  if (nColors<=0 || !engine->palette) {
    /* new color table is void-- don't allow indexed colors */
    psEngine->colorMode= 0;
    psEngine->nColors= 0;

  } else {
    /* remember current color mode, then set to mono mode until
       new color table can be written in BeginPage */
    psEngine->colorMode= 0;
    psEngine->nColors= 0;    /* don't index into table before BeginPage */
  }

  engine->colorChange= 0;
  return 256;
}

/* ------------------------------------------------------------------------ */

static char *psFontNames[N_PSFONTS]= {
  "Courier", "Courier-Bold", "Courier-Oblique", "Courier-BoldOblique",
  "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic",
  "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique",
  "Symbol", "Symbol", "Symbol", "Symbol",
  "NewCenturySchlbk-Roman", "NewCenturySchlbk-Bold",
  "NewCenturySchlbk-Italic", "NewCenturySchlbk-BoldItalic" };

static void Kill(Engine *engine)
{
  PSEngine *psEngine= (PSEngine *)engine;
  long fonts= psEngine->fonts;
  int i, xll, yll, xur, yur;
  int bad= 0;

  if (psEngine->e.marked) bad= EndPage(psEngine);

  if (psEngine->file) {
    if (!bad) bad= PutLine(psEngine);
    if (!bad) bad= Append(psEngine, "%%Trailer");
    if (!bad) bad= PutLine(psEngine);

    sprintf(line, "%%%%Pages: %d", psEngine->currentPage-1);
    if (!bad) bad= Append(psEngine, line);
    if (!bad) bad= PutLine(psEngine);

    xll= psEngine->docBB.xll;
    xur= psEngine->docBB.xur;
    if (xll < xur) {
      yll= psEngine->docBB.yll;
      yur= psEngine->docBB.yur;
    } else {
      xll= yll= 0;
      xur= 612;
      yur= 792;
    }
    sprintf(line, "%%%%BoundingBox: %d %d %d %d", xll, yll, xur, yur);
    if (!bad) bad= Append(psEngine, line);
    if (!bad) bad= PutLine(psEngine);

    strcpy(line, "%%DocumentFonts: ");
    for (i=0 ; i<N_PSFONTS ; i++) {
      if ((1<<i) & fonts) {
        strcat(line, psFontNames[i]);
        if (!bad) bad= Append(psEngine, line);
        if (!bad) bad= PutLine(psEngine);
        strcpy(line, "%%+ ");
      }
    }

    if (psEngine->file!=pf_stdout) p_fclose(psEngine->file);
  }
  GpDelEngine(engine);
}

static int Clear(Engine *engine, int always)
{
  PSEngine *psEngine= (PSEngine *)engine;
  if (always || engine->marked) EndPage(psEngine);
  engine->marked= 0;
  return 0;
}

static int Flush(Engine *engine)
{
  PSEngine *psEngine= (PSEngine *)engine;
  if (psEngine->file && psEngine->file!=pf_stdout)
    p_fflush(psEngine->file);
  return 0;
}

/* ------------------------------------------------------------------------ */

static char *colorCommands[14]= {
  "BG", "FG", "BLK", "WHT", "RED", "GRN", "BLU", "CYA", "MAG", "YEL",
  "GYD", "GYC", "GYB", "GYA" };

static int SetupColor(PSEngine *psEngine, unsigned long color)
{
  int nColors= psEngine->nColors;

  if (!psEngine->e.marked && BeginPage(psEngine)) return 1;
  if (color==psEngine->curColor) return 0;

  if (color<240UL) {
    /* "CI index C"-- "CI" omitted if current color is indexed */
    unsigned long c;
    GpColorCell *palette= psEngine->e.palette;
    if (nColors>0) {
      if (color>=(unsigned long)nColors) color= nColors-1;
      if (psEngine->colorMode) c= color;  /* this page has a color table */
      else c= (P_R(palette[color])+
               P_G(palette[color])+P_B(palette[color]))/3;  /* mono page */
    } else {
      if (color>255) color= 255;
      c= color;         /* No palette ==> no color table on page */
    }
    if (psEngine->curColor>=240UL) {
      if (Append(psEngine, "CI")) return 1;
    }
    sprintf(line, "%ld C", c);
    if (Append(psEngine, line)) return 1;
  } else if (color<256UL) {
    /* standard color command FG, RED, etc. */
    if (color<P_GRAYA) color= P_FG;
    if (Append(psEngine, colorCommands[255UL-color])) return 1;
  } else {
    sprintf(line, "16#%lx C", color);
    if (Append(psEngine, line)) return 1;
  }
  psEngine->curColor= color;

  return 0;
}

static int SetupLine(PSEngine *psEngine, GpLineAttribs *gistAl)
{
  int changeLW= (psEngine->curWidth!=gistAl->width);
  if (SetupColor(psEngine, gistAl->color)) return 1;

  if (changeLW) {
    int lwidth;

    lwidth= (int)(DEFAULT_PS_WIDTH*gistAl->width);
    sprintf(line, "%d LW", lwidth);
    if (Append(psEngine, line)) return 1;

    psEngine->curWidth= gistAl->width;
  }

  /* WARNING--
     the dash pattern is a function of the line width (see pscom.ps) */
  if (psEngine->curType!=gistAl->type ||
      (changeLW && gistAl->type!=L_SOLID)) {
    int ltype= gistAl->type;
    if (ltype==L_NONE) return 1;

    if (ltype<0 || ltype>N_PSDASHES) ltype= L_SOLID;
    sprintf(line, "%d DSH", ltype-1);
    if (Append(psEngine, line)) return 1;

    psEngine->curType= gistAl->type;
  }

  return 0;
}

static int CheckClip(PSEngine *psEngine)
{
  if (gistClip) return BeginClip(psEngine, &gistT);
  else if (psEngine->curClip) return EndClip(psEngine);
  return 0;
}

static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth)
{
  PSEngine *psEngine= (PSEngine *)engine;
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints;
  long np= n + (closed?1:0);
  int firstPass= 1, markEnd= 0;
  GpPoint firstPoint, *points;
  int size;

  if (CheckClip(psEngine)) return 1;
  if (n<1) return 0;
  if (SetupLine(psEngine, &gistA.l)) return 1;
  if (psEngine->curClip) size= 0;
  else size= (int)(psEngine->curWidth*0.5*DEFAULT_PS_WIDTH);

  if (np>90) {
    long nLines= (np-1)/9 + 1;  /* 9 points is 72 characters */
    if (psEngine->nchars && PutLine(psEngine)) return 1;
    sprintf(line, "%%%%BeginData: %ld ASCII Lines", nLines+1);
    if (Append(psEngine, line) || PutLine(psEngine)) return 1;
    markEnd= 1;
  }
  sprintf(line, smooth? "%ld LS" : "%ld L", np);
  if (Append(psEngine, line) || PutLine(psEngine)) return 1;

  while ((nPoints=
          GpIntPoints(map, maxPoints, n, px, py, &points))) {
    if (closed) {
      if (firstPass) {
        firstPoint= points[0];
        firstPass= 0;
      }
      if (n==nPoints) {
        n++;
        points[nPoints++]= firstPoint;
      }
    }
    if (PutPoints(psEngine, points, nPoints, size)) return 1;
    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
  }

  if (markEnd) {
    if (Append(psEngine, "%%EndData") || PutLine(psEngine)) return 1;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static char *psFontCommands[N_PSFONTS]= {
  "0 Cour", "1 Cour", "2 Cour", "3 Cour",
  "0 Tims", "1 Tims", "2 Tims", "3 Tims",
  "0 Helv", "1 Helv", "2 Helv", "3 Helv",
  "0 Symb", "1 Symb", "2 Symb", "3 Symb",
  "0 NCen", "1 NCen", "2 NCen", "3 NCen" };

static int SetupFont(PSEngine *psEngine, GpReal height)
{
  int font= gistA.t.font;
  if (font<0 || font>=N_PSFONTS) font= 0;
  if (psEngine->curFont!=font ||
      psEngine->curHeight!=height) {
    int ptSz= (int)(gistA.t.height*NDC_TO_PS+0.5);
    int lnSp= ptSz;   /* LnSp same as PtSz for now */
    if (Append(psEngine, psFontCommands[font])) return 1;
    sprintf(line, "%d %d FNT", ptSz, lnSp);
    if (Append(psEngine, line)) return 1;
    psEngine->curFont= font;
    psEngine->curHeight= height;
    psEngine->fonts|= (1<<font);
  }
  return 0;
}

static int DrawMarkers(Engine *engine, long n, const GpReal *px,
                       const GpReal *py)
{
  PSEngine *psEngine= (PSEngine *)engine;
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints;
  int type, markEnd= 0;
  GpPoint *points;
  int size;
  char typeString[8];

  if (n<1 || gistA.m.type<=0) return 0;
  if (CheckClip(psEngine)) return 1;
  if (SetupColor(psEngine, gistA.m.color) ||
      SetupFont(psEngine, gistA.m.size*DEFAULT_MARKER_SIZE)) return 1;
  if (psEngine->curClip) size= 0;
  else size= (int)(psEngine->curHeight*NDC_TO_PS);

  if (gistA.m.type>32 && gistA.m.type<127) {
    char *now= typeString;
    *now++= '(';
    if (gistA.m.type=='(' || gistA.m.type==')' || gistA.m.type=='\\')
      *now++= '\\';
    *now++= type= gistA.m.type;
    *now++= ')';
    *now++= '\0';
  } else {
    if (gistA.m.type<0 || gistA.m.type>M_CROSS) type= M_ASTERISK;
    else type= gistA.m.type;
    sprintf(typeString, "%d", type-1);
  }

  if (n>90) {
    long nLines= (n-1)/9 + 1;  /* 9 points is 72 characters */
    if (psEngine->nchars && PutLine(psEngine)) return 1;
    sprintf(line, "%%%%BeginData: %ld ASCII Lines", nLines+1);
    if (Append(psEngine, line) || PutLine(psEngine)) return 1;
    markEnd= 1;
  }
  if (Append(psEngine, typeString)) return 1;  /* "(A)" or "1", e.g. */
  sprintf(line, type<32? "%ld MS" : "%ld M", n);
  if (Append(psEngine, line) || PutLine(psEngine)) return 1;

  while ((nPoints=
          GpIntPoints(map, maxPoints, n, px, py, &points))) {
    if (PutPoints(psEngine, points, nPoints, size)) return 1;
    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
  }

  if (markEnd) {
    if (Append(psEngine, "%%EndData") || PutLine(psEngine)) return 1;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static char *psHCommands[3]= { "/LF", "/CN", "/RT" };
static char *psVCommands[5]= { "/TP", "/CP", "/HF", "/BA", "/BT" };

static int SetupText(PSEngine *psEngine)
{
  int opq;
  int h= gistA.t.alignH, v= gistA.t.alignV;
  GtGetAlignment(&gistA.t, &h, &v);

  if (SetupColor(psEngine, gistA.t.color)) return 1;

  if (psEngine->curAlignH!=h || psEngine->curAlignV!=v) {
    sprintf(line, "%s %s JUS", psHCommands[h-1], psVCommands[v-1]);
    if (Append(psEngine, line)) return 1;
    psEngine->curAlignH= h;
    psEngine->curAlignV= v;
  }

  opq= gistA.t.opaque;
  /* if (gistA.t.orient!=TX_RIGHT) opq= 1; let this be X only limitation */
  if (psEngine->curOpaque != (opq!=0)) {
    if (opq) Append(psEngine, "1 OPQ");
    else if (Append(psEngine, "0 OPQ")) return 1;
    psEngine->curOpaque= (opq!=0);
  }

  if (psEngine->curFont!=gistA.t.font ||
      psEngine->curHeight!=gistA.t.height) {
    if (SetupFont(psEngine, gistA.t.height)) return 1;
  }
  return 0;
}

static int break_line_now(PSEngine *psEngine, char **nw, int *nchs);

static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text)
{
  PSEngine *psEngine= (PSEngine *)engine;
  int ix, iy, nlines;
  GpReal width, height, lineHeight;
  GpXYMap *map= &engine->map;
  GpBox *wind= &engine->transform.window;
  GpReal xmin, xmax, ymin, ymax;
  int xll, yll, xur, yur, alignH, alignV;
  int count, nchars, firstPass;
  char *now, c;
  const char *t;
  int state= 0;

  if (CheckClip(psEngine) || SetupText(psEngine)) return 1;
  alignH= psEngine->curAlignH;
  alignV= psEngine->curAlignV;
  /* Compute text location in PostScript coordinates.  */
  x0= map->x.scale*x0 + map->x.offset;
  y0= map->y.scale*y0 + map->y.offset;

  /* handle multi-line strings */
  nlines= GtTextShape(text, &gistA.t, (WidthFunction)0, &width);
  lineHeight= gistA.t.height * NDC_TO_PS;
  width*= 0.6*lineHeight;
  height= lineHeight*(GpReal)nlines;

  /* Reject if and only if the specified point is off of the current
     page by more than the approximate size of the text.  Note that
     the entire block of text is either accepted or rejected --
     PostScript does the clipping.  */
  if (wind->xmax>wind->xmin) { xmin= wind->xmin; xmax= wind->xmax; }
  else { xmin= wind->xmax; xmax= wind->xmin; }
  if (wind->ymax>wind->ymin) { ymin= wind->ymin; ymax= wind->ymax; }
  else { ymin= wind->ymax; ymax= wind->ymin; }
  if (gistA.t.orient==TX_RIGHT || gistA.t.orient==TX_LEFT) {
    if (x0<xmin-width || x0>xmax+width ||
        y0<ymin-height || y0>ymax+height) return 0;
  } else {
    if (x0<xmin-height || x0>xmax+height ||
        y0<ymin-width || y0>ymax+width) return 0;
  }

  /* Adjust y0 (or x0) to represent topmost line */
  if (nlines > 1) {
    if (gistA.t.orient==TX_RIGHT) {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) y0+= height-lineHeight;
      if (alignV==TV_HALF) y0+= 0.5*(height-lineHeight);
    } else if (gistA.t.orient==TX_LEFT) {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) y0-= height-lineHeight;
      if (alignV==TV_HALF) y0-= 0.5*(height-lineHeight);
    } else if (gistA.t.orient==TX_UP) {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) x0-= height-lineHeight;
      if (alignV==TV_HALF) x0-= 0.5*(height-lineHeight);
    } else {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) x0+= height-lineHeight;
      if (alignV==TV_HALF) x0+= 0.5*(height-lineHeight);
    }
  }

  ix= (int)x0;
  iy= (int)y0;

  /* Guess at the bounding box for the text */
  if (!psEngine->curClip) {
    if (gistA.t.orient==TX_RIGHT) {

      if (alignH==TH_CENTER) x0 -= 0.5*width;
      else if (alignH==TH_RIGHT) x0 -= width;

      if (alignV==TV_TOP || alignV==TV_CAP) y0 -= height;
      else if (alignV==TV_HALF) y0 -= height-0.4*lineHeight;
      else if (alignV==TV_BASE) y0 -= height-0.8*lineHeight;
      else if (alignV==TV_BOTTOM) y0 -= height-lineHeight;

    } else if (gistA.t.orient==TX_LEFT) {

      if (alignH==TH_CENTER) x0 -= 0.5*width;
      else if (alignH==TH_LEFT) x0 -= width;

      if (alignV==TV_HALF) y0 -= .4*lineHeight;
      else if (alignV==TV_BASE) y0 -= 0.8*lineHeight;
      else if (alignV==TV_BOTTOM) y0 -= lineHeight;

    } else if (gistA.t.orient==TX_UP) {
      
      if (alignH==TH_CENTER) y0 -= 0.5*width;
      else if (alignH==TH_RIGHT) y0 -= width;

      if (alignV==TV_HALF) x0 -= 0.4*lineHeight;
      else if (alignV==TV_BASE) x0 -= 0.8*lineHeight;
      else if (alignV==TV_BOTTOM) x0 -= lineHeight;

    } else { /* TX_DOWN */

      if (alignH==TH_CENTER) y0 -= 0.5*width;
      else if (alignH==TH_LEFT) y0 -= width;

      if (alignV==TV_TOP || alignV==TV_CAP) x0 -= height;
      else if (alignV==TV_HALF) x0 -= height-0.4*lineHeight;
      else if (alignV==TV_BASE) x0 -= height-0.8*lineHeight;
      else if (alignV==TV_BOTTOM) x0 -= height-lineHeight;

    }

    if (x0>xmin) xmin= x0;
    if (x0+width<xmax) xmax= x0+width;
    if (y0>ymin) ymin= y0;
    if (y0+height<ymax) ymax= y0+height;    

    xll= (int)xmin;
    xur= (int)xmax;
    yll= (int)ymin;
    yur= (int)ymax;
    if (xll<psEngine->pageBB.xll) psEngine->pageBB.xll= xll;
    if (xur>psEngine->pageBB.xur) psEngine->pageBB.xur= xur;
    if (yll<psEngine->pageBB.yll) psEngine->pageBB.yll= yll;
    if (yur>psEngine->pageBB.yur) psEngine->pageBB.yur= yur;
  }

  if (nlines>1 || gistA.t.orient!=TX_RIGHT) {
    if (Append(psEngine, "[")) return 1;
  }

  nchars= psEngine->nchars;
  now= psEngine->line+nchars;
  firstPass= 1;
  while ((t= GtNextLine(text, &count, gistA.t.orient)) || firstPass) {
    text= t;
    firstPass= 0;
    state= 0;

    if (nchars > 70) {
      psEngine->nchars= nchars;
      if (PutLine(psEngine)) return 1;
      nchars= 0;
      now= psEngine->line;
    }
    *now++= ' ';
    *now++= '(';
    nchars+= 2;

    while (count--) {
      if (nchars>72 && break_line_now(psEngine, &now, &nchars)) return 1;

      c= *text;
      if (c>=32 && c<127) {
        if (gtDoEscapes) {
          if (c=='!' && count) {
            c= text[1];  /* peek at next char */
            if (c=='!' || c=='^' || c=='_') {
              count--;
              text++;
              goto normal;
            }
            if (nchars>68 &&
                break_line_now(psEngine, &now, &nchars)) return 1;
            strcpy(now, "\\024");  /* DC4 is ps.ps escape char */
            now+= 4;
            nchars+= 4;
            psEngine->fonts|= (1<<T_SYMBOL);  /* symbol font used */
            if (c==']') {
              *now++= '^';         /* !] means ^ (perp) in symbol font */
              nchars++;
              text+= 2;
              count--;
              continue;
            }
            goto escape;
          } else if (c=='^') {
            if (state!=0) {
              /* terminate current super or subscript escape */
              strcpy(now, "\\021");  /* DC1 is ps.ps escape char */
              now+= 4;
              nchars+= 4;
            }
            if (state!=1) {
              /* intialize superscript escape */
              if (nchars>64 &&
                  break_line_now(psEngine, &now, &nchars)) return 1;
              strcpy(now, "\\021\\022");  /* DC1DC2 is ps.ps escape seq */
              now+= 8;
              nchars+= 8;
              state= 1;
            } else {
              state= 0;
            }
            goto escape;
          } else if (c=='_') {
            if (state!=0) {
              /* terminate current super or subscript escape */
              if (nchars>68 &&
                  break_line_now(psEngine, &now, &nchars)) return 1;
              strcpy(now, "\\021");  /* DC1 is ps.ps escape char */
              now+= 4;
              nchars+= 4;
            }
            if (state!=2) {
              /* intialize subscript escape */
              if (nchars>64 &&
                  break_line_now(psEngine, &now, &nchars)) return 1;
              strcpy(now, "\\021\\023");  /* DC1DC3 is ps.ps escape seq */
              now+= 8;
              nchars+= 8;
              state= 2;
            } else {
              state= 0;
            }
            goto escape;
          }
        }
        /* ordinary printing character */
        if (c=='(' || c==')' || c=='\\') {
          *now++= '\\';
          nchars++;
        }
      normal:
        *now++= c;
        nchars++;
      } else {
        /* non-printing characters rendered as \ooo */
        if (c=='\t') {
          *now++= '\\';
          *now++= 't';
          nchars+= 2;
        } else if (c<'\021' || c>'\024') {
          /* DC1 through DC4 have special meaning in ps.ps, skip them */
          sprintf(now, "\\%03o", (int)((unsigned char)c));
          now+= 4;
          nchars+= 4;
        }
      }
    escape:
      text++;
    }
    if (state!=0) {
      /* terminate current super or subscript escape */
      strcpy(now, "\\021");  /* DC1 is ps.ps escape char */
      now+= 4;
      nchars+= 4;
    }

    *now++= ')';
    nchars++;
    *now= '\0';
  }
  psEngine->nchars= nchars;

  if (gistA.t.orient==TX_RIGHT) {
    sprintf(line, nlines>1? "] %d %d TA" : "%d %d T", ix, iy);
  } else {
    int angle;
    if (gistA.t.orient==TX_LEFT) angle= 180;
    else if (gistA.t.orient==TX_UP) angle= 90;
    else angle= 270;
    sprintf(line, "] %d %d %d TR", angle, ix, iy);
  }
  if (Append(psEngine, line)) return 1;

  return 0;
}

static int break_line_now(PSEngine *psEngine, char **nw, int *nchs)
{
  char *now= *nw;
  *now++= '\\';
  *now= '\0';
  psEngine->nchars= *nchs+1;
  if (PutLine(psEngine)) return 1;
  *nchs= 0;
  *nw= psEngine->line;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py)
{
  PSEngine *psEngine= (PSEngine *)engine;
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints;
  int markEnd= 0;
  GpPoint *points;
  int value= 0;

  /* For now, only FillSolid style supported */

  if (n<1) return 0;
  if (CheckClip(psEngine) || SetupColor(psEngine, gistA.f.color)) return 1;

  if (n>90) {
    long nLines= (n-1)/9 + 1;  /* 9 points is 72 characters */
    if (psEngine->nchars && PutLine(psEngine)) return 1;
    sprintf(line, "%%%%BeginData: %ld ASCII Lines", nLines+1);
    if (Append(psEngine, line) || PutLine(psEngine)) return 1;
    markEnd= 1;
  }
  if (gistA.e.type==L_NONE) sprintf(line, "%ld F", n);
  else sprintf(line, "%ld E", n);
  if (Append(psEngine, line) || PutLine(psEngine)) return 1;

  while ((nPoints=
          GpIntPoints(map, maxPoints, n, px, py, &points))) {
    if (PutPoints(psEngine, points, nPoints, 0)) return 1;
    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
    value= 1;   /* Polygons with >4050 sides won't be filled correctly */
  }
  if (gistA.e.type!=L_NONE) {
    /* setup for edge (usually different color than fill), then draw it */
    if (SetupLine(psEngine, &gistA.e)) return 1;
    if (Append(psEngine, "0 E") || PutLine(psEngine)) return 1;
  }

  if (markEnd) {
    if (Append(psEngine, "%%EndData") || PutLine(psEngine)) return 1;
  }

  return value;
}

/* ------------------------------------------------------------------------ */

static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors)
{
  PSEngine *psEngine= (PSEngine *)engine;
  GpXYMap *map= &psEngine->e.map;
  int nColors= psEngine->nColors;
  GpColorCell *palette;
  int ix, iy, idx, idy, depth;
  long i, j, off;
  int markEnd= 0;
  long nLines;
  char *now= psEngine->line;
  int nc, ncmax, color, colorMode;

  if (!psEngine->e.marked && BeginPage(psEngine)) return 1;
  if (CheckClip(psEngine)) return 1;

  /* Transform corner coordinates, clipping and adjusting width,
     height, nColumns, and colors as necessary.  */
  width = GpClipCells(&map->x, &px, &qx,
                      gistT.window.xmin, gistT.window.xmax, width, &off);
  colors += gistA.rgb? 3*off : off;
  height = GpClipCells(&map->y, &py, &qy,
                       gistT.window.ymin, gistT.window.ymax, height, &off);
  colors += gistA.rgb? 3*nColumns*off : nColumns*off;

  if (width<=0 || height<=0) return 0;
  ix= (int)px;
  iy= (int)py;
  idx= (int)(qx-px);
  idy= (int)(qy-py);

  /* Set bounding box for image if necessary */
  if (!psEngine->curClip) {
    GpBox *wind= &engine->transform.window;
    GpReal xmin, xmax, ymin, ymax;
    int xll, yll, xur, yur;
    if (wind->xmax>wind->xmin) { xmin= wind->xmin; xmax= wind->xmax; }
    else { xmin= wind->xmax; xmax= wind->xmin; }
    if (wind->ymax>wind->ymin) { ymin= wind->ymin; ymax= wind->ymax; }
    else { ymin= wind->ymax; ymax= wind->ymin; }

    if (px<qx) {
      if (px>xmin) xmin= px;
      if (qx<xmax) xmax= qx;
    } else {
      if (qx>xmin) xmin= qx;
      if (px<xmax) xmax= px;
    }
    if (py<qy) {
      if (py>ymin) ymin= py;
      if (qy<ymax) ymax= qy;
    } else {
      if (qy>ymin) ymin= qy;
      if (py<ymax) ymax= py;
    }

    xll= (int)xmin;
    xur= (int)xmax;
    yll= (int)ymin;
    yur= (int)ymax;
    if (xll<psEngine->pageBB.xll) psEngine->pageBB.xll= xll;
    if (xur>psEngine->pageBB.xur) psEngine->pageBB.xur= xur;
    if (yll<psEngine->pageBB.yll) psEngine->pageBB.yll= yll;
    if (yur>psEngine->pageBB.yur) psEngine->pageBB.yur= yur;
  }

  /* Use 4 bit depth if the color table is small, otherwise use
     8 bit depth.  */
  if (nColors>0) {
    /* Image value will be either  */
    colorMode= psEngine->colorMode;
    if (colorMode) {
      palette= 0;                     /* palette already written */
      if (nColors>16) depth= 8;
      else depth= 4;
    } else {
      palette= psEngine->e.palette;   /* must lookup gray level now */
      depth= 8;
    }
  } else {
    /* Must assume image varies over maximum possible range */
    colorMode= 1;  /* That is, use index without trying palette lookup */
    depth= 8;
    palette= 0;
  }

  if (gistA.rgb) {
    /* Write the 6 arguments to the J procedure */
    sprintf(line, "%d %d %d %d %d %d",
            (int)width, (int)height, idx, idy, ix, iy);
  } else {
    /* Write the 7 arguments to the I procedure */
    sprintf(line, "%d %d %d %d %d %d %d",
            (int)width, (int)height, depth, idx, idy, ix, iy);
  }
  if (Append(psEngine, line)) return 1;

  if (gistA.rgb) {
    nLines= 6*width*height;
    ncmax = 72;   /* 12 cells (72 chars) per line */
  } else {
    nLines= width*height;
    depth= (depth==8);
    if (depth) nLines*= 2;
    else if (nLines&1L) nLines++;
    ncmax = 76;   /* Will put 76 or 38 cells per line */
  }
  nLines= 1+(nLines-1)/ncmax;
  if (nLines>10) {
    if (psEngine->nchars && PutLine(psEngine)) return 1;
    sprintf(line, "%%%%BeginData: %ld ASCII Lines", nLines+1);
    if (Append(psEngine, line) || PutLine(psEngine)) return 1;
    markEnd= 1;
  }
  if (Append(psEngine, gistA.rgb?"J":"I") || PutLine(psEngine)) return 1;

  i= j= 0;
  while (nLines--) {
    for (nc=0 ; nc<ncmax ; ) {
      if (i>=width) {
        height--;
        if (height<=0) break;
        i= 0;
        j+= nColumns;
      }
      if (gistA.rgb) {
        const GpColor *ccell = &colors[3*(i+j)];
        now[nc++] = hexChar[ccell[0]>>4];
        now[nc++] = hexChar[ccell[0]&0xf];
        now[nc++] = hexChar[ccell[1]>>4];
        now[nc++] = hexChar[ccell[1]&0xf];
        now[nc++] = hexChar[ccell[2]>>4];
        now[nc++] = hexChar[ccell[2]&0xf];
      } else {
        color= colors[i+j];
        i++;
        if (color>=nColors && nColors>0) color= nColors-1;
        if (!colorMode) color= (P_R(palette[color])+
                                P_G(palette[color])+P_B(palette[color]))/3;
        if (depth) {  /* 2 hex chars per cell */
          now[nc++]= hexChar[color>>4];
          now[nc++]= hexChar[color&0xf];
        } else {      /* 1 hex char per cell */
          now[nc++]= hexChar[color];
        }
      }
    }
    now[nc]= '\0';
    psEngine->nchars= nc;
    if (PutLine(psEngine)) return 1;
  }

  if (markEnd) {
    if (Append(psEngine, "%%EndData") || PutLine(psEngine)) return 1;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy)
{
  PSEngine *psEngine= (PSEngine *)engine;
  GpXYMap *map= &engine->map;
  long maxSegs= 2025, nSegs;
  int markEnd= 0;
  GpSegment *segs;
  int size;

  if (CheckClip(psEngine)) return 1;
  if (n<1 || SetupLine(psEngine, &gistA.l)) return 0;
  if (psEngine->curClip) size= 0;
  else size= (int)(psEngine->curWidth*0.5*DEFAULT_PS_WIDTH);

  if (n>45) {
    long nLines= (2*n-1)/9 + 1;  /* 9 points is 72 characters */
    if (psEngine->nchars && PutLine(psEngine)) return 1;
    sprintf(line, "%%%%BeginData: %ld ASCII Lines", nLines+1);
    if (Append(psEngine, line) || PutLine(psEngine)) return 1;
    markEnd= 1;
  }
  sprintf(line, "%ld D", n);
  if (Append(psEngine, line) || PutLine(psEngine)) return 1;

  while ((nSegs=
          GpIntSegs(map, maxSegs, n, px, py, qx, qy, &segs))) {
    /* This cast from GpPoint to GpSegment is lazy programming,
       but I don't know of any platforms it will not work on... */
    if (PutPoints(psEngine, (GpPoint *)segs, 2*nSegs, size)) return 1;
    if (n==nSegs) break;
    n-= nSegs;
    px+= nSegs;
    py+= nSegs;
    qx+= nSegs;
    qy+= nSegs;
  }

  if (markEnd) {
    if (Append(psEngine, "%%EndData") || PutLine(psEngine)) return 1;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

Engine *GpPSEngine(char *name, int landscape, int mode, char *file)
{
  PSEngine *psEngine;
  long flen= file? strlen(file) : 0;
  long engineSize= sizeof(PSEngine)+flen+1;
  GpTransform toPixels;

  if (flen<=0) return 0;

  SetPSTransform(&toPixels, landscape);

  psEngine=
    (PSEngine *)GpNewEngine(engineSize, name, psType, &toPixels, landscape,
                            &Kill, &Clear, &Flush, &GpComposeMap,
                            &ChangePalette, &DrawLines, &DrawMarkers,
                            &DrwText, &DrawFill, &DrawCells,
                            &DrawDisjoint);

  if (!psEngine) {
    strcpy(gistError, "memory manager failed in GpPSEngine");
    return 0;
  }

  psEngine->filename= (char *)(psEngine+1);
  strcpy(psEngine->filename, file);
  psEngine->file= 0;
  psEngine->closed= 0;

  SetPageDefaults(psEngine);
  psEngine->e.colorMode= mode;
  psEngine->colorMode= 0;
  psEngine->nColors= 0;

  psEngine->landscape= landscape;
  psEngine->docBB.xll= psEngine->docBB.xur=
    psEngine->docBB.yll= psEngine->docBB.yur= 0;
  psEngine->currentPage= 1;
  psEngine->fonts|= (1<<T_COURIER);

  psEngine->clipColor= psEngine->curColor;
  psEngine->clipType= psEngine->curType;
  psEngine->clipWidth= psEngine->curWidth;
  psEngine->clipFont= psEngine->curFont;
  psEngine->clipHeight= psEngine->curHeight;

  psEngine->clipBox.xmin= psEngine->clipBox.xmax=
    psEngine->clipBox.ymin= psEngine->clipBox.ymax= 0.0;

  psEngine->line[0]= '\0';
  psEngine->nchars= 0;

  return (Engine *)psEngine;
}

PSEngine *GisPSEngine(Engine *engine)
{
  return (engine && engine->type==psType)? (PSEngine *)engine : 0;
}

/* ------------------------------------------------------------------------ */
