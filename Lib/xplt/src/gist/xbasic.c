/*
 * XBASIC.C
 *
 * $Id$
 *
 * Implement the basic X windows engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "xbasic.h"
#include "gtext.h"
#include "pstdlib.h"

#include <string.h>

static char *xType= "PLAY Window";

static int ChangePalette(Engine *engine);
static GpReal TextWidth(const char *text, int nc, const GpTextAttribs *t);
static void Kill(Engine *engine);
static void ShutDown(XEngine *xEngine);
static int Clear(Engine *engine, int always);
static int Flush(Engine *engine);
static void GetXRectangle(GpXYMap *map, GpBox *box,
                          int *x0, int *y0, int *x1, int *y1);
static void ClearArea(Engine *engine, GpBox *box);
static void SetXTransform(GpTransform *trans, int landscape, int dpi);
static void gx_translate(GpBox *trans_window, int x, int y);
static GpBox *DamageClip(GpBox *damage);
static void ChangeMap(Engine *engine);
static void chk_clipping(XEngine *xeng);
static int SetupLine(XEngine *xeng, GpLineAttribs *gistAl, int join);
static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth);
static int DrawMarkers(Engine *engine, long n, const GpReal *px,
                       const GpReal *py);
static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text);
static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py);
static int GetCells(GpMap *map, GpReal xmin, GpReal xmax,
                    GpReal px, GpReal qx, long width,
                    int *i0, int *di, int *ncols, int *x0, int *x1);
static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors);
static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy);
static void GetVisibleNDC(XEngine *xeng,
                          GpReal *xn, GpReal *xx, GpReal *yn, GpReal *yx);

/* Hook for hlevel.c error handling.  */
extern void (*HLevelHook)(Engine *engine);


extern int GxJustifyText(GpXYMap *map, GpReal x0, GpReal y0, const char *text,
                         int *ix, int *iy, int xbox[], int ybox[]);
extern int GxJustifyNext(const char **text, int *ix, int *iy);

static int gxErrorFlag= 0;
static void GxErrorHandler(void);

/* ------------------------------------------------------------------------ */

static int
ChangePalette(Engine *engine)
{
  XEngine *xeng= (XEngine *)engine;
  p_scr *s = xeng->s;
  p_win *win = xeng->win;
  GpColorCell *palette= engine->palette;  /* requested palette */
  int nColors= engine->nColors;
  int width, height, depth;

  if (!s) return 0;
  depth = p_sshape(s, &width, &height);
  if (depth > 8) depth = 8;
  if (nColors > 256) nColors= 256;

  p_palette(win, palette, nColors);

  xeng->e.colorChange= 0;

  return (1<<depth);
}

/* ------------------------------------------------------------------------ */

static int nChars, lineHeight, textHeight, prevWidth, maxWidth, alignH;
static int textAscent;
static int firstTextLine= 0;

static p_scr *current_scr;
static p_win *current_win;
static int current_state, current_fsym, current_fsize;
static int chunkWidth, nChunk, supersub, dy_super, dy_sub, x_chunks;

/* ARGSUSED */
static GpReal
TextWidth(const char *text, int nc, const GpTextAttribs *t)
{
  int width= 0;
  if (firstTextLine) nChars= nc;
  if (!gtDoEscapes) {
    width = p_txwidth(current_scr, text, nc, gistA.t.font, current_fsize);
    if (firstTextLine) {
      /* record width of first line */
      prevWidth= chunkWidth= width;
      nChunk= nc;
      firstTextLine= 0;
    }
  } else if (nc>0) {
    int firstChunk= firstTextLine;
    const char *txt= text;
    char c;
    for ( ; nc-- ; text++) {
      c= text[0];
      if ((nc && c=='!') || c=='^' || c=='_') {
        if (txt<text)
          width += p_txwidth(current_scr, txt, (int)(text-txt),
                             gistA.t.font, current_fsize);
        if (firstChunk) {
          nChunk= (int)(text-txt);
          chunkWidth= width;
          firstChunk= 0;
        }
        txt= text+1;
        if (c=='!') {
          /* process single character escapes immediately */
          text= txt;
          nc--;
          c= txt[0];
          if (c!='!' && c!='^' && c!='_') {
            if (c==']') c= '^';
            width += p_txwidth(current_scr, &c, 1,
                               current_fsym, current_fsize);
            txt++;
          }
        } else if (c=='^') {
          supersub|= 1;
        } else {
          supersub|= 2;
        }
      }
    }
    if (txt<text)
      width += p_txwidth(current_scr, txt, (int)(text-txt),
                         gistA.t.font, current_fsize);
    if (firstChunk) {
      nChunk= (int)(text-txt);
      chunkWidth= width;
    }
    if (firstTextLine) {
      /* record width of first line */
      prevWidth= width;  /* this is whole line */
      firstTextLine= 0;
    }
  }
  return (GpReal)width;
}

int
GxJustifyText(GpXYMap *map, GpReal x0, GpReal y0, const char *text,
              int *ix, int *iy, int xbox[], int ybox[])
{
  int nLines, alignV, dx, dy, xmin, xmax, ymin, ymax, ix0, iy0;
  GpReal rwidest;

  /* abort if string screen position is ridiculous */
  x0 = map->x.scale*x0 + map->x.offset;
  y0 = map->y.scale*y0 + map->y.offset;
  if (x0<-16000. || x0>16000. || y0<-16000. || y0>16000.) return -1;
  current_state = 0;

  /* call p_font before any calls to TextWidth
   * - may improve efficiency of p_txwidth, p_txheight */
  p_font(current_win, gistA.t.font, current_fsize, gistA.t.orient);

  /* Set nLines, maxWidth, nChars, prevWidth */
  firstTextLine = 1;
  nChars = prevWidth = chunkWidth = nChunk = supersub = 0;
  nLines = GtTextShape(text, &gistA.t, &TextWidth, &rwidest);
  maxWidth = (int)rwidest;

  /* state bits:
     1 - this chunk is superscript
     2 - this chunk is subscript (1 and 2 together meaningless)
     4 - this chunk is single symbol char
   */
  x_chunks= 0;

  /* Compute height of one line */
  textHeight = p_txheight(current_scr, gistA.t.font, current_fsize,
                          &textAscent);
  dy_super = (supersub&1)? textAscent/3 : 0;
  dy_sub = (supersub&2)? textAscent/3 : 0;
  lineHeight = textHeight+dy_sub+dy_super;

  /* Compute displacement and bounding box of entire string,
     relative to specified location */
  GtGetAlignment(&gistA.t, &alignH, &alignV);
  /* Note: xmax= xmin+maxWidth  */
  if (alignH==TH_LEFT) {
    dx= xmin= 0;
    xmax= maxWidth;
  } else if (alignH==TH_CENTER) {
    xmax= maxWidth/2;
    xmin= -xmax;
    dx= -prevWidth/2;
  } else {
    xmax= 0;
    xmin= -maxWidth;
    dx= -prevWidth;
  }

  /* Note: ymax= ymin+nLines*lineHeight  and  ymin= dy-textAscent  */
  if (alignV<=TV_CAP) {
    dy= textAscent + dy_super;
    ymin= 0;
    ymax= nLines*lineHeight;
  } else if (alignV==TV_HALF) {
    ymin= textAscent/2;
    ymax= nLines*lineHeight;
    dy= ymin-(ymax-lineHeight)/2;
    ymin= dy-textAscent-dy_super;
    ymax+= ymin;
  } else if (alignV==TV_BASE) {
    dy= -(nLines-1)*lineHeight;
    ymin= dy-textAscent-dy_super;
    ymax= textHeight-textAscent + dy_sub;
  } else {
    ymin= dy_sub-nLines*lineHeight;
    dy= ymin+textAscent+dy_super;
    ymax= dy_sub;
  }

  /* handle orientation (path) of text */
  if (gistA.t.orient==TX_LEFT) {         /* upside down */
    int tmp;
    dx= -dx;  dy= -dy;
    tmp= xmin;  xmin= -xmax;  xmax= -tmp;
    tmp= ymin;  ymin= -ymax;  ymax= -tmp;
  } else if (gistA.t.orient==TX_UP) {    /* reading upwards */
    int tmp;
    tmp= dx;  dx= dy;  dy= -tmp;
    tmp= xmin;  xmin= ymin;  ymin= -xmax;
                xmax= ymax;  ymax= -tmp;
  } else if (gistA.t.orient==TX_DOWN) {  /* reading downwards */
    int tmp;
    tmp= dx;  dx= -dy;  dy= tmp;
    tmp= xmin;  xmin= -ymax;  ymax= xmax;
                xmax= -ymin;  ymin= tmp;
  }

  /* get bounding box and adjusted reference point */
  ix0 = (int)x0;
  iy0 = (int)y0;

  xbox[0] = ix0+xmin;   ybox[0] = iy0+ymin;
  xbox[1] = ix0+xmax;   ybox[1] = iy0+ymax;

  *ix= dx + ix0;
  *iy= dy + iy0;

  if (nChunk) p_font(current_win,
                     gistA.t.font, current_fsize, gistA.t.orient);
  return nChunk;
}

int
GxJustifyNext(const char **text, int *ix, int *iy)
{
  const char *txt= *text+nChunk;
  int xadj= 0, yadj= 0;
  char c;

  nChars-= nChunk;
  if (!nChars) {
    /* last chunk was the last one on a line */
    txt= GtNextLine(txt, &nChars, 0);
    if (!txt) return -1;

    *text= txt;

    /* scan for end of first chunk */
    if (gtDoEscapes) {
      for (nChunk=0 ; nChunk<nChars ; nChunk++) {
        c= txt[nChunk];
        if ((nChunk+1<nChars && c=='!') || c=='^' || c=='_') break;
      }
    } else {
      nChunk= nChars;
    }

    /* compute width of this chunk if necessary, compute width
       of whole line if necessary for justification */
    if (alignH!=TH_LEFT || gistA.t.orient!=TX_RIGHT) {
      /* need to compute width of entire line */
      int width= prevWidth;
      firstTextLine= 1;
      prevWidth= (int)TextWidth(txt, nChars, &gistA.t);
      if (alignH==TH_CENTER) xadj= (width-prevWidth)/2;
      else if (alignH==TH_RIGHT) xadj= width-prevWidth;
      /* TextWidth sets chunkWidth */
    } else if (nChunk<nChars) {
      /* just need width of this chunk */
      if (nChunk) chunkWidth = p_txwidth(current_scr, txt, nChunk,
                                         gistA.t.font, current_fsize);
      else chunkWidth= 0;  /* unused */
    }

    /* reset state, adjusting (possibly rotated) x and y as well */
    xadj-= x_chunks;
    yadj= lineHeight;
    if (current_state&1) yadj+= dy_super;
    else if (current_state&2) yadj-= dy_sub;
    if (nChunk && (current_state&4))
      p_font(current_win, gistA.t.font, current_fsize, gistA.t.orient);
    current_state= 0;

  } else {
    /* previous chunk ended with an escape character, or was single
       escaped symbol character -- can't get here unles gtDoEscapes */
    char c1= '\0';
    xadj= chunkWidth;      /* width of previous chunk */
    x_chunks+= chunkWidth; /* accumulate all chunks except last */
    yadj= 0;
    if (!(current_state&4)) {
      c1= *txt++;
      nChars--;
    }
    *text= txt;

    if (c1=='!') {
      /* this chunk begins with escaped character */
      nChunk= 1;
      c= txt[0];
      if (c=='!' || c=='^' || c=='_') {
        /* chunk is just ordinary text */
        for ( ; nChunk<nChars ; nChunk++) {
          c= txt[nChunk];
          if ((nChunk+1<nChars && c=='!') || c=='^' || c=='_') break;
        }
        p_font(current_win, gistA.t.font, current_fsize, gistA.t.orient);
        current_state&= 3;
      } else {
        /* chunk is single symbol char */
        p_font(current_win, current_fsym, current_fsize, gistA.t.orient);
        current_state|= 4;
      }

    } else {
      for (nChunk=0 ; nChunk<nChars ; nChunk++) {
        c= txt[nChunk];
        if ((nChunk+1<nChars && c=='!') || c=='^' || c=='_') break;
      }
      if (nChunk)
        p_font(current_win, gistA.t.font, current_fsize, gistA.t.orient);
      if (c1=='^') {
        if (current_state&1) {
          yadj+= dy_super;  /* return from super to normal */
          current_state= 0;
        } else {
          if (current_state&2) yadj-= dy_sub;
          yadj-= dy_super;  /* move to superscript */
          current_state= 1;
        }
      } else if (c1=='_') {
        if (current_state&2) {
          yadj-= dy_sub;  /* return from sub to normal */
          current_state= 0;
        } else {
          if (current_state&1) yadj+= dy_super;
          yadj+= dy_sub;  /* move to subscript */
          current_state= 2;
        }
      } else {
        /* just finished a symbol char */
        current_state&= 3;
      }
    }

    if (nChunk &&
        (nChunk<nChars || alignH!=TH_LEFT || gistA.t.orient!=TX_RIGHT)) {
      char caret= '^';
      if (nChunk==1 && (current_state&4) && txt[0]==']') txt= &caret;
      chunkWidth = p_txwidth(current_scr, txt, nChunk,
                             (current_state&4)? current_fsym : gistA.t.font,
                             current_fsize);
    } else {
      chunkWidth= 0;  /* unused */
    }
  }

  if (gistA.t.orient==TX_RIGHT) {
    *iy+= yadj;
    *ix+= xadj;
  } else if (gistA.t.orient==TX_LEFT) {
    *iy-= yadj;
    *ix-= xadj;
  } else if (gistA.t.orient==TX_UP) {
    *ix+= yadj;
    *iy-= xadj;
  } else {
    *ix-= yadj;
    *iy+= xadj;
  }

  return nChunk;
}

/* ------------------------------------------------------------------------ */

/* notes for Kill() and g_on_destroy
 * (1) Kill() is program-driven (e.g.- winkill)
 *     g_on_destroy() is event-driven (e.g.- mouse click)
 * (2) however, the p_destroy() function MIGHT or MIGHT NOT
 *     result in a call to g_on_destroy()
 *     under Windows, g_on_destroy() is naturally called
 *       by the call p_destroy() must make to close the window
 *     under UNIX/X, the play event handler calls p_destroy()
 *       after g_on_destroy() returns
 * (3) therefore g_on_destroy() must not call p_destroy(), since
 *       this would cause infinite recursion on Windows, and a
 *       double call on X.
 * (4) worse, if this is the final window on an X display, gist
 *       should disconnect from the display, but that can only
 *       be done after the final p_destroy(), which is after
 *       the g_on_destroy() in the event-driven case
 *       -- thus, the p_disconnect must take place on the next
 *          event, which is during the GhBeforeWait hlevel.c function
 *       -- unfortunately, in the meantime, a new window could
 *          have been created on that screen,
 *          hence the g_test_pending() routine
 *       -- the mess is further compounded by the fact that there
 *          could be an unlimited number of different screens with
 *          pending disconnections, after the 4th one, the current
 *          code simply gives up and forgets about disconnecting
 */
extern void (*g_pending_task)(void);
void (*g_pending_task)(void) = 0;
static void g_test_pending(p_scr *s);

static void
Kill(Engine *engine)
{
  XEngine *xeng= (XEngine *)engine;
  p_win *w = xeng->win;
  ShutDown(xeng);
  if (w) {
    xeng->w = 0;
    p_destroy(w);
  }
  /* for program-driven Kill(), can take care of p_disconnect immediately */
  if (g_pending_task) g_pending_task();
}

static int
Clear(Engine *engine, int always)
{
  XEngine *xeng = (XEngine *)engine;
  if (!xeng->w) return 1;
  if ((always || xeng->e.marked) && xeng->w==xeng->win) {
    int tm = xeng->topMargin;
    int lm = xeng->leftMargin;
    if (tm || lm) {
      int xmax = (int)xeng->swapped.window.xmax;
      int ymax = (int)xeng->swapped.window.ymin;
      if (xmax > lm+xeng->wtop) xmax = lm+xeng->wtop;
      if (ymax > tm+xeng->htop) ymax = tm+xeng->htop;
      if (xeng->clipping) {
        p_clip(xeng->w, 0,0,0,0);
        xeng->clipping = 0;
      }
      p_color(xeng->w, P_BG);
      p_rect(xeng->w, lm, tm, xmax, ymax, 0);
    } else {
      p_clear(xeng->w);
    }
  }
  if (xeng->e.colorChange) ChangePalette(engine);
  xeng->e.marked = 0;
  return 0;
}

static int
Flush(Engine *engine)
{
  XEngine *xeng = (XEngine *)engine;
  if (!xeng->w) return 1;
  p_flush(xeng->w);
  /* test whether an X error has been reported */
  if (gxErrorFlag) GxErrorHandler();
  return 0;
}

static void
GetXRectangle(GpXYMap *map, GpBox *box,
              int *x0, int *y0, int *x1, int *y1)
{
  /* get corners of clip rectangle in pixel coordinates */
  int wmin= (int)(map->x.scale*box->xmin+map->x.offset);
  int wmax= (int)(map->x.scale*box->xmax+map->x.offset);
  if (wmax>=wmin) {
    *x0 = wmin;
    *x1 = wmax+1;
  } else {
    *x0 = wmax;
    *x1 = wmin+1;
  }
  wmin= (int)(map->y.scale*box->ymin+map->y.offset);
  wmax= (int)(map->y.scale*box->ymax+map->y.offset);
  if (wmax>=wmin) {
    *y0 = wmin;
    *y1 = wmax+1;
  } else {
    *y0 = wmax;
    *y1 = wmin+1;
  }
}

static void
ClearArea(Engine *engine, GpBox *box)
{
  XEngine *xeng= (XEngine *)engine;
  p_win *w = xeng->w;
  int x0, y0, x1, y1;
  if (!w) return;
  /* if this is animation mode, do not try to clear window */
  if (w==xeng->win) {
    int lm = xeng->leftMargin;
    int tm = xeng->topMargin;
    GetXRectangle(&engine->devMap, box, &x0, &y0, &x1, &y1);
    if (x0 < lm) x0 = lm;
    if (x1 > lm+xeng->wtop) x1 = lm+xeng->wtop;
    if (y0 < tm) y0 = tm;
    if (y1 > tm+xeng->htop) y1 = tm+xeng->htop;
    p_color(w, P_BG);
    p_rect(w, x0, y0, x1, y1, 0);
  }
}

static void
SetXTransform(GpTransform *trans, int landscape, int dpi)
{
  trans->viewport= landscape? gLandscape : gPortrait;
  trans->window.xmin= 0.0;
  trans->window.xmax= PixelsPerNDC(dpi)*trans->viewport.xmax;
  trans->window.ymin= PixelsPerNDC(dpi)*trans->viewport.ymax;
  trans->window.ymax= 0.0;
}

static void
gx_translate(GpBox *trans_window, int x, int y)
{
  trans_window->xmax += x - trans_window->xmin;
  trans_window->xmin = x;
  trans_window->ymin += y - trans_window->ymax;
  trans_window->ymax = y;
}

static GpBox cPort;

static GpBox *
DamageClip(GpBox *damage)
{
  cPort= gistT.viewport;
  if (cPort.xmin>cPort.xmax)
    { GpReal tmp= cPort.xmin; cPort.xmin= cPort.xmax; cPort.xmax= tmp; }
  if (cPort.ymin>cPort.ymax)
    { GpReal tmp= cPort.ymin; cPort.ymin= cPort.ymax; cPort.ymax= tmp; }
  /* (assume damage box is properly ordered) */
  if (damage->xmin>cPort.xmin) cPort.xmin= damage->xmin;
  if (damage->xmax<cPort.xmax) cPort.xmax= damage->xmax;
  if (damage->ymin>cPort.ymin) cPort.ymin= damage->ymin;
  if (damage->ymax<cPort.ymax) cPort.ymax= damage->ymax;
  if (cPort.xmin>cPort.xmax || cPort.ymin>cPort.ymax) return 0;
  else return &cPort;
}

static void
ChangeMap(Engine *engine)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  int landscape = xeng->width > xeng->height;
  int x0, y0, x1, y1;
  GpBox *clipport;
  if (!w) return;

  /* check to be sure that landscape/portrait mode hasn't changed */
  if (landscape!=xeng->e.landscape) {
    /* this is probably insane if in animation mode... */
    SetXTransform(&xeng->e.transform, xeng->e.landscape, xeng->dpi);
    xeng->width = (int)xeng->e.transform.window.xmax;
    xeng->height = (int)xeng->e.transform.window.ymin;
    xeng->swapped = xeng->e.transform;
    /* make adjustments to allow for SetXTransform, then recenter */
    if (xeng->w != xeng->win) {
      xeng->a_x += xeng->x+1;
      xeng->a_y += xeng->y+1;
    }
    xeng->x = xeng->y = -1;
    GxRecenter(xeng, xeng->wtop+xeng->leftMargin, xeng->htop+xeng->topMargin);
  }

  /* do generic change map */
  GpComposeMap(engine);

  /* get current clip window */
  if (xeng->e.damaged) clipport = DamageClip(&xeng->e.damage);
  else clipport = &gistT.viewport;
  if (clipport) {
    /* set clipping rectangle for this XEngine */
    GetXRectangle(&engine->devMap, clipport, &x0, &y0, &x1, &y1);
    if (xeng->w == xeng->win) {
      /* additional restriction for vignetting by window borders */
      int lm = xeng->leftMargin;
      int tm = xeng->topMargin;
      if (x0 < lm) x0 = lm;
      if (x1 > lm+xeng->wtop) x1 = lm+xeng->wtop;
      if (y0 < tm) y0 = tm;
      if (y1 > tm+xeng->htop) y1 = tm+xeng->htop;
      xeng->clipping = 1;
    }
    if (x1<=x0) x1 = x0+1;
    if (y1<=y0) y1 = y0+1;
    p_clip(xeng->w, x0, y0, x1, y1);
  }
}

static void
chk_clipping(XEngine *xeng)
{
  p_win *w = xeng->win;
  if (!xeng->clipping) {
    int x0, y0, x1, y1;
    int lm = xeng->leftMargin;
    int tm = xeng->topMargin;
    if (xeng->e.damaged) {
      GpBox *box = DamageClip(&xeng->e.damage);
      GpXYMap map;
      if (xeng->w != w)
        GpSetMap(&xeng->swapped.viewport, &xeng->swapped.window, &map);
      else
        map = xeng->e.devMap;
      GetXRectangle(&map, box, &x0, &y0, &x1, &y1);
      /* additional restriction for vignetting by window borders */
      if (x0 < lm) x0 = lm;
      if (x1 > lm+xeng->wtop) x1 = lm+xeng->wtop;
      if (y0 < tm) y0 = tm;
      if (y1 > tm+xeng->htop) y1 = tm+xeng->htop;
    } else {
      x0 = lm;
      x1 = lm+xeng->wtop;
      y0 = tm;
      y1 = tm+xeng->htop;
    }
    xeng->clipping = 1;
    if (x1<=x0) x1 = x0+1;
    if (y1<=y0) y1 = y0+1;
    p_clip(xeng->w, x0, y0, x1, y1);
  }
}

/* ------------------------------------------------------------------------ */

static int
SetupLine(XEngine *xeng, GpLineAttribs *gistAl, int join)
{
  GpXYMap *map= &xeng->e.map;
  double xt[2], yt[2];
  xt[0] = map->x.scale;
  xt[1] = map->x.offset;
  yt[0] = map->y.scale;
  yt[1] = map->y.offset;
  p_d_map(xeng->w, xt, yt, 1);
  chk_clipping(xeng);
  if (gistAl->type != L_NONE) {
    int type = gistAl->type-1;
    int width = (unsigned int)(DEFAULT_LINE_INCHES*xeng->dpi*gistAl->width);
    if (join) type |= P_SQUARE;
    p_pen(xeng->w, width, type);
    p_color(xeng->w, gistAl->color);
    return 0;
  } else {
    return 1;
  }
}

static int
DrawLines(Engine *engine, long n, const GpReal *px,
          const GpReal *py, int closed, int smooth)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  long i, imax;
  int npts;

  if (!w) return 1;
  if (n<=0 || SetupLine(xeng, &gistA.l, 0)) return 0;

  closed = (closed && n>1 && (px[0]!=px[n-1] || py[0]!=py[n-1]));
  for (i=0 ; i<n ; i=imax) {
    imax = i+2047;
    npts = (imax<=n)? 2047 : (int)(n-i);
    p_d_pnts(w, px+i, py+i, npts);
    if (closed && imax>=n)
      p_d_pnts(w, px, py, -1);
    p_lines(w);
  }

  xeng->e.marked = 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

/* play has no polymarker primitive, use GpPseudoMark-- unless
   user requested tiny points, in which case we use p_dots */
static int
DrawMarkers(Engine *engine, long n, const GpReal *px, const GpReal *py)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  if (!w || !xeng->mapped) return 1;

  xeng->e.marked = 1;
  if (gistA.m.type!=M_POINT || gistA.m.size>1.5) {
    return GpPseudoMark(engine, n, px, py);

  } else {
    long i, imax;
    int npts;
    GpXYMap *map= &xeng->e.map;
    double xt[2], yt[2];
    xt[0] = map->x.scale;
    xt[1] = map->x.offset;
    yt[0] = map->y.scale;
    yt[1] = map->y.offset;
    p_d_map(w, xt, yt, 1);
    chk_clipping(xeng);
    p_color(w, gistA.m.color);

    for (i=0 ; i<n ; i=imax) {
      imax = i+2048;
      npts = (imax<=n)? 2048 : (int)(n-i);
      p_d_pnts(w, px+i, py+i, npts);
      p_dots(w);
    }

    return 0;
  }
}

/* ------------------------------------------------------------------------ */

static int
DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  GpXYMap *map = &xeng->e.map;
  int ix, iy, len, xbox[2], ybox[2], xn, xx, yn, yx;
  const char *txt;
  char *caret= "^";

  if (!w || !xeng->mapped) return 1;
  chk_clipping(xeng);

  current_fsize = (int)((xeng->dpi/ONE_INCH)*gistA.t.height);
  if (current_fsize<4) current_fsize = 4;     /* totally illegible */
  if (current_fsize>180) current_fsize = 180; /* way too big */
  current_fsym = T_SYMBOL | (gistA.t.font&3);
  current_scr = xeng->s;
  current_win = w;

  /* get current window */
  xn = (int)(gistT.window.xmin*map->x.scale + map->x.offset);
  xx = (int)(gistT.window.xmax*map->x.scale + map->x.offset);
  yn = (int)(gistT.window.ymin*map->y.scale + map->y.offset);
  yx = (int)(gistT.window.ymax*map->y.scale + map->y.offset);
  if (yn > yx) { int tmp=yn ; yn=yx; yx=tmp ; }

  /* handle multi-line strings */
  len = GxJustifyText(map, x0, y0, text, &ix, &iy, xbox, ybox);
  if (len < 0) return 0;

  /* consider whether string is completely clipped */
  if (ybox[0]>yx || ybox[1]<yn || xbox[0]>xx || xbox[1]<xn) return 0;

  /* erase background if string is opaque */
  if (gistA.t.opaque) {
    p_color(w, P_BG);
    p_rect(w, xbox[0], ybox[0], xbox[1], ybox[1], 0);
  }
  p_color(w, gistA.t.color);

  do {
    if (len>0) {
      if (len==1 && (current_state&4) && text[0]==']') txt = caret;
      else txt = text;
      p_text(w, ix, iy, txt, len);
    }
    len = GxJustifyNext(&text, &ix, &iy);
  } while (len>=0);

  xeng->e.marked = 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int
DrawFill(Engine *engine, long n, const GpReal *px,
         const GpReal *py)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  long i, imax;
  int npts, has_edge;

  if (!w || !xeng->mapped) return 1;
  has_edge = !SetupLine(xeng, &gistA.e, 0); /* but shouldn't set color */
  p_color(w, gistA.f.color);

  /* This gives incorrect results if more than one pass through the loop,
   * but there is no way to give the correct result, so may as well...  */
  for (i=0 ; i<n ; i=imax) {
    imax = i+2048;
    npts = (imax<=n)? 2048 : (int)(n-i);
    p_d_pnts(w, px+i, py+i, npts);
    /* Can Nonconvex or Convex be detected? */
    p_fill(w, 0);
  }

  xeng->e.marked = 1;
  if (has_edge) {
    p_color(w, gistA.e.color);
    for (i=0 ; i<n ; i=imax) {
      imax = i+2047;
      npts = (imax<=n)? 2047 : (int)(n-i);
      p_d_pnts(w, px+i, py+i, npts);
      if (imax>=n)
        p_d_pnts(w, px, py, -1);
      p_lines(w);
    }
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int
GetCells(GpMap *map, GpReal xmin, GpReal xmax,
         GpReal px, GpReal qx, long width,
         int *i0, int *di, int *ncols, int *x0, int *x1)
{
  GpReal scale = map->scale;
  GpReal offset = map->offset;
  GpReal x, dx = (qx-px)/width;
  int imin, imax;

  if (xmin>xmax) x=xmin, xmin=xmax, xmax=x;

  if (dx>=0.0) {
    if (qx<xmin || px>xmax) return 0;
    if (dx > 0.0) {
      x = (xmin - px)/dx;
      if (x < 1.) imin = 0;
      else imin=(int)x, px+=imin*dx;
      x = (qx - xmax)/dx;
      if (x < 1.) imax = 0;
      else imax=(int)x, qx-=imax*dx;
      width -= imin + imax;
      if (width < 3) { /* see comment below */
        if (qx <= xmax) {
          if (imin) px = xmin;
        } else if (px >= xmin) {
          if (imax) qx = xmax;
        } else if (width < 2) {
          px = xmin;
          qx = xmax;
        } else if (qx-xmax <= xmin-px) {
          px += qx - xmax;
          qx = xmax;
        } else {
          qx += qx - xmin;
          px = xmin;
        }
      }
    } else {
      imin = width/2;
      imax = width - imin - 1;
      width = 1;
    }
  } else {
    if (px<xmin || qx>xmax) return 0;
    dx = -dx;
    x = (px - xmax)/dx;
    if (x < 1.) imin = 0;
    else imin=(int)x, px-=imin*dx;
    x = (xmin - qx)/dx;
    if (x < 1.) imax = 0;
    else imax=(int)x, qx+=imax*dx;
    width -= imin + imax;
    if (width < 3) { /* see comment below */
      if (px <= xmax) {
        if (imax) qx = xmin;
      } else if (qx >= xmin) {
        if (imin) px = xmax;
      } else if (width < 2) {
        qx = xmin;
        px = xmax;
      } else if (px-xmax <= xmin-qx) {
        qx += px - xmax;
        px = xmax;
      } else {
        px += qx - xmin;
        qx = xmin;
      }
    }
  }
  /* width<3 logic above guarantees these will not overflow as int */
  px = px*scale + offset;
  qx = qx*scale + offset;
  if (qx >= px) {
    *i0 = imin;
    *di = 1;
    *ncols = width;
    *x0 = (int)px;
    *x1 = (int)qx;
  } else {
    *i0 = imin + width - 1;
    *di = -1;
    *ncols = width;
    *x0 = (int)qx;
    *x1 = (int)px;
  }
  /* cell array always at least 1 pixel */
  if (*x1 == *x0) *x1 += 1;

  return 1;
}

static int
DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
          GpReal qy, long width, long height, long nColumns,
          const GpColor *colors)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  GpXYMap *map= &xeng->e.map;
  int i0, j0, di, dj, ncols, nrows, x0, y0, x1, y1;

  if (!w || !xeng->mapped) return 1;
  chk_clipping(xeng);

  if (GetCells(&map->x, gistT.window.xmin, gistT.window.xmax,
               px, qx, width, &i0, &di, &ncols, &x0, &x1) &&
      GetCells(&map->y, gistT.window.ymin, gistT.window.ymax,
               py, qy, height, &j0, &dj, &nrows, &y0, &y1)) {
    unsigned char *ndxs = (unsigned char *)colors;
    if (di<0 || dj<0 || ncols!=width || nrows!=height || nColumns!=width) {
      int r, c, nr, nc, i, j;
      ndxs = p_malloc(gistA.rgb?3*ncols*nrows:ncols*nrows);
      j0 *= nColumns;
      dj *= nColumns;
      if (gistA.rgb) {
        for (j=j0,c=0,nr=nrows ; nr-- ; j+=dj,c+=ncols) {
          nc = ncols;
          for (i=i0,r=0 ; nc-- ; i+=di,r++) {
            ndxs[3*(c+r)] = colors[3*(j+i)];
            ndxs[3*(c+r)+1] = colors[3*(j+i)+1];
            ndxs[3*(c+r)+2] = colors[3*(j+i)+2];
          }
        }
      } else {
        for (j=j0,c=0,nr=nrows ; nr-- ; j+=dj,c+=ncols) {
          nc = ncols;
          for (i=i0,r=0 ; nc-- ; i+=di,r++)
            ndxs[c+r] = colors[j+i];
        }
      }
    }
    if (ncols && nrows) {
      if (gistA.rgb)
        p_rgb_cell(w, ndxs, ncols, nrows, x0, y0, x1, y1);
      else
        p_ndx_cell(w, ndxs, ncols, nrows, x0, y0, x1, y1);
    }
    if (ndxs!=colors) p_free(ndxs);
  }

  xeng->e.marked = 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int
DrawDisjoint(Engine *engine, long n, const GpReal *px,
             const GpReal *py, const GpReal *qx, const GpReal *qy)
{
  XEngine *xeng = (XEngine *)engine;
  p_win *w = xeng->w;
  long i, imax;
  int nseg;

  if (!w || !xeng->mapped) return 1;
  if (SetupLine(xeng, &gistA.l, 1)) return 0;

  p_d_pnts(w, px, py, 0);
  for (i=0 ; i<n ;) {
    imax = i+1024;
    nseg = (imax<=n)? 1024 : (int)(n-i);
    while (nseg--) {
      p_d_pnts(w, px+i, py+i, -1);
      p_d_pnts(w, qx+i, qy+i, -1);
      i++;
    }
    p_segments(w);
  }

  xeng->e.marked = 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

extern int (*gg_on_expose)(void *c, int *xy);
extern int (*gg_on_destroy)(void *c);
extern int (*gg_on_resize)(void *c,int w,int h);
extern int (*gg_on_focus)(void *c,int in);
extern int (*gg_on_key)(void *c,int k,int md);
extern int (*gg_on_click)(void *c,int b,int md,int x,int y, unsigned long ms);
extern int (*gg_on_motion)(void *c,int md,int x,int y);
extern int (*gg_on_deselect)(void *c);
int (*gg_on_expose)(void *c, int *xy) = 0;
int (*gg_on_destroy)(void *c) = 0;
int (*gg_on_resize)(void *c,int w,int h) = 0;
int (*gg_on_focus)(void *c,int in) = 0;
int (*gg_on_key)(void *c,int k,int md);
int (*gg_on_click)(void *c,int b,int md,int x,int y, unsigned long ms) = 0;
int (*gg_on_motion)(void *c,int md,int x,int y) = 0;
int (*gg_on_deselect)(void *c) = 0;

static void g_on_expose(void *c, int *xy);
static void g_on_destroy(void *c);
static void g_on_resize(void *c,int w,int h);
static void g_on_focus(void *c,int in);
static void g_on_key(void *c,int k,int md);
static void g_on_click(void *c,int b,int md,int x,int y, unsigned long ms);
static void g_on_motion(void *c,int md,int x,int y);
static void g_on_deselect(void *c);
static void g_on_panic(p_scr *screen);

static Engine *waiting_for = 0;
static void (*wait_callback)(void) = 0;

extern int gist_expose_wait(Engine *eng, void (*e_callback)(void));

int
gist_expose_wait(Engine *eng, void (*e_callback)(void))
{
  if (waiting_for) {
    waiting_for = 0;
    wait_callback = 0;
    return 1;
  } else {
    XEngine *xeng = GisXEngine(eng);
    if (!xeng || !xeng->w) return 1;
    if (xeng->mapped) return 2;
    waiting_for = eng;
    wait_callback = e_callback;
    return 0;
  }
}

static void
g_on_expose(void *c, int *xy)
{
  if (!gg_on_expose || gg_on_expose(c,xy)) {
    XEngine *xeng = c;
    if (c && c == waiting_for) {
      waiting_for = 0;
      if (wait_callback) wait_callback();
      wait_callback = 0;
    }
    if (!xeng->w) return;
    xeng->mapped = 1;
    if (xeng->HandleExpose)
      /* the alternate handler should probably call GxExpose */
      xeng->HandleExpose(&xeng->e, xeng->e.drawing, xy);
    else
      GxExpose(&xeng->e, xeng->e.drawing, xy);
  }
}

static void
g_on_click(void *c,int b,int md,int x,int y, unsigned long ms)
{
  if (!gg_on_click || gg_on_click(c,b,md,x,y,ms)) {
    XEngine *xeng = c;
    if (!xeng->w) return;
    if (xeng->HandleClick)
      xeng->HandleClick(&xeng->e, b, md, x, y, ms);
  }
}

static void
g_on_motion(void *c,int md,int x,int y)
{
  if (!gg_on_motion || gg_on_motion(c,md,x,y)) {
    XEngine *xeng = c;
    if (!xeng->w) return;
    if (xeng->HandleMotion)
      xeng->HandleMotion(&xeng->e, md, x, y);
  }
}

static void
g_on_destroy(void *c)
{
  if (!gg_on_destroy || gg_on_destroy(c)) {
    XEngine *xeng = c;
    /* if xeng->win==0, assume ShutDown already called */
    if (xeng->win) ShutDown(xeng);
  }
}

static void
g_on_resize(void *c,int w,int h)
{
  if (!gg_on_resize || gg_on_resize(c,w,h)) {
    XEngine *xeng = c;
    if (xeng->w) GxRecenter(xeng, w, h);
  }
}

static void
g_on_focus(void *c,int in)
{
  if (!gg_on_focus || gg_on_focus(c,in)) {
    XEngine *xeng = c;
    if (xeng->w && xeng->HandleMotion && in==2)
      xeng->HandleMotion(&xeng->e, 0, -1, -1);
  }
}

static void
g_on_key(void *c,int k,int md)
{
  if (!gg_on_key || gg_on_key(c,k,md)) {
    XEngine *xeng = c;
    if (xeng->w && xeng->HandleKey)
      xeng->HandleKey(&xeng->e, k, md);
  }
}

static void g_on_deselect(void *c)
{
  if (!gg_on_deselect || gg_on_deselect(c)) {
  }
}

void
GxExpose(Engine *engine, Drauing *drawing, int *xy)
{
  XEngine *xeng = (XEngine *)engine;
  GpBox damage;
  if (!drawing || !xeng->w) return;
  /* xy=0 to redraw all, otherwise x0,y0,x1,y1 */
  if (xy) {
    GpXYMap *map = &engine->devMap;
    damage.xmin= (xy[0]-map->x.offset)/map->x.scale;
    damage.xmax= (xy[2]-map->x.offset)/map->x.scale;
    damage.ymax= (xy[1]-map->y.offset)/map->y.scale;
    damage.ymin= (xy[3]-map->y.offset)/map->y.scale;
  } else {
    damage.xmin = xeng->swapped.viewport.xmin;
    damage.xmax = xeng->swapped.viewport.xmax;
    damage.ymin = xeng->swapped.viewport.ymin;
    damage.ymax = xeng->swapped.viewport.ymax;
  }
  if (engine->damaged) {
    GpSwallow(&engine->damage, &damage);
  } else {
    engine->damage = damage;
    engine->damaged = 1;
  }
  GdSetDrawing(drawing);
  GpPreempt(engine);
  GdDraw(1);
  GpPreempt(0);        /* not correct if damaged during a preempt... */
  GdSetDrawing(0);
}

void
GxRecenter(XEngine *xeng, int width, int height)
{
  int x, y;
  int eWidth = xeng->width;
  int eHeight = xeng->height;
  width -= xeng->leftMargin;
  height -= xeng->topMargin;
  xeng->wtop = width;
  xeng->htop = height;
  x = (eWidth-width)/2;
  /* put center of page at center of landscape window */
  if (eWidth>eHeight) y = (eHeight-height)/2;
  /* put center of upper square of page at center of portrait window */
  else y = (eWidth-height)/2;
  /* once either dimension is big enough for whole picture, stop moving it */
  if (y<0) y = 0;
  if (x<0) x = 0;
  if (x!=xeng->x || y!=xeng->y) {
    int tmargin = xeng->topMargin;
    int lmargin = xeng->leftMargin;
    gx_translate(&xeng->swapped.window, -x+lmargin, -y+tmargin);
    if (xeng->w == xeng->win) {
      gx_translate(&xeng->e.transform.window, -x+lmargin, -y+tmargin);
      GpDeviceMap(&xeng->e);
    } else {
      xeng->a_x -= x - xeng->x;
      xeng->a_y -= y - xeng->y;
      lmargin = tmargin = 0;
    }
    xeng->x = x;
    xeng->y = y;
    if (xeng->wtop>0) x = xeng->wtop+lmargin;
    else x = lmargin+1;
    if (xeng->htop>0) y = xeng->htop+tmargin;
    else y = tmargin+1;
    xeng->clipping = 1;
    p_clip(xeng->win, lmargin, tmargin, x, y);
  }
}

/* ------------------------------------------------------------------------ */

typedef struct g_scr g_scr;
struct g_scr {
  char *name;
  int number;
  p_scr *s;
};
static g_scr *g_screens = 0;
static int n_screens = 0;

/* ARGSUSED */
void
g_initializer(int *pargc, char *argv[])
{
  extern char *g_argv0;
  g_argv0 = argv? argv[0] : 0;
  p_gui(&g_on_expose, &g_on_destroy, &g_on_resize, &g_on_focus,
        &g_on_key, &g_on_click, &g_on_motion, &g_on_deselect,
        &g_on_panic);
}

p_scr *
g_connect(char *displayName)
{
  p_scr *s = 0;
  int i, j, i0=-1, len=0, number=0;

  /* split display into base name and screen number (separated by dot) */
  if (displayName) while (displayName[len]) len++;
  if (len) {
    for (i=len-1 ; i>=0 ; i--) if (displayName[i]=='.') break;
    if (i>=0) {
      int i0 = i;
      for (i++ ; i<len && displayName[i]<='9' && displayName[i]>='0' ; i++)
        number = 10*number + (displayName[i]-'0');
      if (i == len) len = i0;
      else number = 0;
    }
  }
  if (!len) displayName = 0;
  if (g_screens) {
    for (i=0 ; i<n_screens ; i++) {
      j = 0;
      if (g_screens[i].name) {
        for ( ; j<len ; j++)
          if (g_screens[i].s && g_screens[i].name[j]!=displayName[j]) break;
      }
      if (j==len && (len? (!g_screens[i].name[j]) : !g_screens[i].name)) {
        if (number == g_screens[i].number) break;
        else if (i0<0) i0 = i;
      }
    }
    if (i<n_screens) s = g_screens[i].s;
  }
  if (!s) {
    if (i0<0) s = p_connect(displayName);
    else s = p_multihead(g_screens[i0].s, number);
    if (!s) return s;
    g_test_pending(s);
    for (i=0 ; i<n_screens ; i++) if (!g_screens[i].s) break;
    if (i==n_screens && !(i & (i-1))) {
      int n = i? 2*i : 1;
      g_screens = p_realloc(g_screens, sizeof(g_scr)*n);
    }
    g_screens[i].number = number;
    g_screens[i].name = displayName? p_strncat(0, displayName, len) : 0;
    g_screens[i].s = s;
    if (i==n_screens) n_screens++;
  }

  return s;
}

void
g_disconnect(p_scr *s)
{
  if (s) {
    int i;
    char *name;
    for (i=0 ; i<n_screens ; i++) {
      if (g_screens[i].s == s) {
        name = g_screens[i].name;
        g_screens[i].name = 0;
        g_screens[i].s = 0;
        p_free(name);
      }
    }
    p_disconnect(s);
  }
}

XEngine *
GxEngine(p_scr *s, char *name, GpTransform *toPixels,
         int x, int y, int topMargin, int leftMargin, long engineSize)
{
  XEngine *xEngine;
  unsigned int width, height;
  GpReal pixels_per_page;
  int dpi;

  if (!s) return 0;

  /* Graphics window will have dimensions of toPixels transform window */
  if (toPixels->window.xmin<toPixels->window.xmax)
    width = (unsigned int)(toPixels->window.xmax - toPixels->window.xmin);
  else
    width = (unsigned int)(toPixels->window.xmin - toPixels->window.xmax);
  if (toPixels->window.ymin<toPixels->window.ymax)
    height = (unsigned int)(toPixels->window.ymax - toPixels->window.ymin);
  else
    height = (unsigned int)(toPixels->window.ymin - toPixels->window.ymax);

  /* Reconstruct dpi (dots per inch) from toPixels transform */
  pixels_per_page = toPixels->window.ymin;
  if (pixels_per_page < toPixels->window.xmax)
    pixels_per_page = toPixels->window.xmax;
  dpi = (int)(0.01 + pixels_per_page*ONE_INCH/gPortrait.ymax);

  /* adjust VDC window so GpDeviceMap in GpNewEngine sets proper
   * transform, which will have xmin=x<0, ymax=y<0 */
  gx_translate(&toPixels->window, x+leftMargin, y+topMargin);

  xEngine =
    (XEngine *)GpNewEngine(engineSize, name, xType, toPixels, width>height,
                           &Kill, &Clear, &Flush, &ChangeMap,
                           &ChangePalette, &DrawLines, &DrawMarkers,
                           &DrwText, &DrawFill, &DrawCells,
                           &DrawDisjoint);
  if (!xEngine) {
    strcpy(gistError, "memory manager failed in GxEngine");
    return 0;
  }

  /* XEngines can repair damage */
  xEngine->e.ClearArea = &ClearArea;

  /* Fill in Engine properties specific to XEngine */
  xEngine->s = s;
  xEngine->win = 0;
  xEngine->width = width;
  xEngine->height = height;
  xEngine->topMargin = topMargin;
  xEngine->leftMargin = leftMargin;
  xEngine->x = -x;
  xEngine->y = -y;
  xEngine->mapped = xEngine->clipping = 0;
  xEngine->dpi = dpi;

  xEngine->e.colorMode = 0;

  xEngine->w = 0;
  xEngine->a_width = xEngine->a_height= 0;
  xEngine->a_x = xEngine->a_y= 0;
  xEngine->swapped = xEngine->e.transform;

  xEngine->HandleExpose = 0;
  xEngine->HandleClick = 0;
  xEngine->HandleMotion = 0;
  xEngine->HandleKey = 0;

  return xEngine;
}

/* default top window represents 6 inch square */
int gx75width = 450;
int gx100width = 600;
int gx75height = 450;
int gx100height = 600;

int gist_private_map = 0;
int gist_input_hint = 0;
int gist_rgb_hint = 0;

Engine *
GpBXEngine(char *name, int landscape, int dpi, char *displayName)
{
  p_scr *s = g_connect(displayName);
  int topWidth = DefaultTopWidth(dpi);
  int topHeight = DefaultTopHeight(dpi);
  GpTransform toPixels;
  int x, y, width, height, hints;
  XEngine *xeng;

  if (!s) return 0;

  SetXTransform(&toPixels, landscape, dpi);
  width = (int)toPixels.window.xmax;
  height = (int)toPixels.window.ymin;
  x = (width-topWidth)/2;
  if (landscape) y = (height-topHeight)/2;
  else y = (width-topHeight)/2;
  if (y<0) y = 0;
  if (x<0) x = 0;
  xeng = GxEngine(s, name, &toPixels, -x,-y,0,0, sizeof(XEngine));

  xeng->wtop = topWidth;
  xeng->htop = topHeight;
  /* possibly want optional P_RGBMODEL as well */
  hints = (gist_private_map?P_PRIVMAP:0) | (gist_input_hint?0:P_NOKEY) |
    (gist_rgb_hint?P_RGBMODEL:0);
  xeng->win = xeng->w =
    p_window(s, topWidth, topHeight, name, P_BG, hints, xeng);
  if (!xeng->win) {
    GpDelEngine(&xeng->e);
    return 0;
  }

  return &xeng->e;
}

int
GxInput(Engine *engine,
        void (*HandleExpose)(Engine *, Drauing *, int *),
        void (*HandleClick)(Engine *,int,int,int,int,unsigned long),
        void (*HandleMotion)(Engine *,int,int,int),
        void (*HandleKey)(Engine *,int,int))
{
  XEngine *xeng = GisXEngine(engine);
  if (!xeng) return 1;
  xeng->HandleExpose = HandleExpose;
  xeng->HandleClick = HandleClick;
  xeng->HandleMotion = HandleMotion;
  xeng->HandleKey = HandleKey;
  return 0;
}

XEngine *
GisXEngine(Engine *engine)
{
  return (engine && engine->type==xType)? (XEngine *)engine : 0;
}

/* ------------------------------------------------------------------------ */

int
GxAnimate(Engine *engine, GpBox *viewport)
{
  XEngine *xeng = GisXEngine(engine);
  int x, y, x0, y0, x1, y1;
  GpBox *v, *w;
  GpReal xmin, xmax, ymin, ymax;
  GpReal scalx, offx, scaly, offy;

  if (!xeng || !xeng->w) return 1;
  if (xeng->w!=xeng->win) GxDirect(engine);

  v = &xeng->e.transform.viewport;  /* NDC */
  w = &xeng->e.transform.window;    /* pixels */

  /* get NDC-->pixel mapping coefficients */
  scalx = xeng->e.devMap.x.scale;
  offx = xeng->e.devMap.x.offset;
  scaly = xeng->e.devMap.y.scale;
  offy = xeng->e.devMap.y.offset;

  /* clip given viewport to portion of NDC space which is actually
   * visible now -- note that v is either gLandscape or gPortrait,
   * so that min<max for v; must also be true for input viewport */
  GetVisibleNDC(xeng, &xmin, &xmax, &ymin, &ymax);
  if (viewport->xmin>xmin) xmin = viewport->xmin;
  if (viewport->xmax<xmax) xmax = viewport->xmax;
  if (viewport->ymin>ymin) ymin = viewport->ymin;
  if (viewport->ymax<ymax) ymax = viewport->ymax;

  /* install NDC-->pixel transform for animation pixmap */
  v->xmin = xmin;
  v->xmax = xmax;
  v->ymin = ymin;
  v->ymax = ymax;

  /* set the engine transform to map the specified viewport into
   * offscreen pixels, and get (x,y) offset from full window pixels
   * to offscreen pixels */
  w->xmin = scalx*xmin+offx;
  w->xmax = scalx*xmax+offx;
  if (w->xmax > w->xmin) {
    x = (int)w->xmin;
    w->xmax -= w->xmin;
    w->xmin = 0.0;
  } else {
    x = (int)w->xmax;
    w->xmin -= w->xmax;
    w->xmax = 0.0;
  }
  w->ymin = scaly*ymin+offy;
  w->ymax = scaly*ymax+offy;
  if (w->ymax > w->ymin) {
    y = (int)w->ymin;
    w->ymax -= w->ymin;
    w->ymin = 0.0;
  } else {
    y = (int)w->ymax;
    w->ymin -= w->ymax;
    w->ymax = 0.0;
  }
  GpDeviceMap((Engine *)xeng);
  GetXRectangle(&xeng->e.devMap, v, &x0, &y0, &x1, &y1);
  x1 -= x0;
  y1 -= y0;

  /* create the offscreen pixmap */
  xeng->w = p_offscreen(xeng->win, x1, y1);
  if (!xeng->w) {
    xeng->w = xeng->win;
    xeng->e.transform = xeng->swapped;
    GpDeviceMap((Engine *)xeng);
    return 2;
  }
  xeng->a_width = x1;
  xeng->a_height = y1;
  xeng->a_x = x;
  xeng->a_y = y;

  /* set coordinate mapping for offscreen */
  ChangeMap((Engine *)xeng);

  /* reset mapping clip to whole visible window */
  if (xeng->wtop>0) x1 = xeng->wtop+xeng->leftMargin;
  else x1 = xeng->leftMargin+1;
  if (xeng->htop>0) y1 = xeng->htop+xeng->topMargin;
  else y1 = xeng->topMargin+1;
  xeng->clipping = 1;
  p_clip(xeng->win, xeng->leftMargin, xeng->topMargin, x1, y1);

  p_clear(xeng->w);
  return 0;
}

static void
GetVisibleNDC(XEngine *xeng,
              GpReal *xn, GpReal *xx, GpReal *yn, GpReal *yx)
{
  GpReal scalx = xeng->e.devMap.x.scale;
  GpReal offx = xeng->e.devMap.x.offset;
  GpReal scaly = xeng->e.devMap.y.scale;
  GpReal offy = xeng->e.devMap.y.offset;
  int xmin, xmax, ymin, ymax;

  xmin = xeng->leftMargin;
  xmax = xmin+xeng->wtop;
  ymax = xeng->topMargin;
  ymin = ymax+xeng->htop;

  /* Convert pixels to NDC coordinates */
  *xn = (xmin-offx)/scalx;
  *xx = (xmax-offx)/scalx;
  *yn = (ymin-offy)/scaly;
  *yx = (ymax-offy)/scaly;
}

int
GxStrobe(Engine *engine, int clear)
{
  XEngine *xeng = GisXEngine(engine);

  if (!xeng || !xeng->w || xeng->w==xeng->win) return 1;

  p_bitblt(xeng->win, xeng->a_x, xeng->a_y, xeng->w,
           0, 0, xeng->a_width, xeng->a_height);
  if (clear) p_clear(xeng->w);

  return 0;
}

int
GxDirect(Engine *engine)
{
  XEngine *xeng = GisXEngine(engine);

  if (!xeng || !xeng->w || xeng->w==xeng->win) return 1;

  p_destroy(xeng->w);
  xeng->w = xeng->win;

  /* set coordinate map and clipping to values for graphics window */
  xeng->e.transform = xeng->swapped;
  GpDeviceMap((Engine *)xeng);
  ChangeMap((Engine *)xeng);

  return 0;
}

/* ------------------------------------------------------------------------ */

void (*HLevelHook)(Engine *engine)= 0;

/* hack to disconnect if last engine destroyed (see GhBeforeWait) */
#define G_N_PENDING 4
static p_scr *g_pending_scr[G_N_PENDING];
static void g_do_disconnect(void);
static void
g_do_disconnect(void)
{
  p_scr *s;
  int i;
  for (i=0 ; i<=G_N_PENDING ; i++) {
    s = g_pending_scr[i];
    g_pending_scr[i] = 0;
    if (s) g_disconnect(s);
  }
  g_pending_task = 0;
}

static void
g_test_pending(p_scr *s)
{
  int i;
  if (g_pending_task == g_do_disconnect) {
    for (i=0 ; i<=G_N_PENDING ; i++)
      if (g_pending_scr[i] == s) {
        g_pending_scr[i] = 0;
        break;
      }
  } else {
    g_pending_task = 0;
    for (i=0 ; i<=G_N_PENDING ; i++) g_pending_scr[i] = 0;
  }
}

static void
ShutDown(XEngine *xeng)
{
  p_scr *s = xeng->s;
  p_win *w = xeng->w;
  p_win *win = xeng->win;
  xeng->mapped= 0;
  /* destroy any hlevel references without further ado */
  if (HLevelHook) HLevelHook((Engine *)xeng);
  xeng->w = xeng->win = 0;
  xeng->s = 0;
  if (w && w!=win) p_destroy(w);
  GpDelEngine(&xeng->e);
  if (s) {
    Engine *eng;
    XEngine *xe2;
    for (eng=GpNextEngine(0) ; eng ; eng=GpNextEngine(eng)) {
      xe2 = GisXEngine(eng);
      if (xe2 && xe2->s==s) break;
    }
    if (!eng) {
      int i;
      if (g_pending_task == g_do_disconnect) {
        for (i=0 ; i<=G_N_PENDING ; i++)
          if (g_pending_scr[i] == s) break;
        if (i >= G_N_PENDING) {
          for (i=0 ; i<=G_N_PENDING ; i++)
            if (!g_pending_scr[i]) break;
          if (i < G_N_PENDING) g_pending_scr[i] = s;
        }
      } else {
        g_pending_scr[0] = s;
        for (i=1 ; i<=G_N_PENDING ; i++) g_pending_scr[i] = 0;
        g_pending_task = g_do_disconnect;
      }
    }
  }
}

static void (*XErrHandler)(char *errMsg)= 0;

static void
g_on_panic(p_scr *screen)
{
  Engine *eng = 0;
  XEngine *xeng = 0;
  do {
    for (eng=GpNextEngine(eng) ; eng ; eng=GpNextEngine(eng)) {
      xeng= GisXEngine(eng);
      if (xeng && xeng->s==screen) break;
    }
    if (eng) {
      xeng->s = 0;  /* be sure not to call p_disconnect */
      Kill(eng);
    }
  } while (eng);
  XErrHandler("play on_panic called (screen graphics engines killed)");
}

/* this routine actually calls the XErrHandler, which may not return
   and/or which may trigger additional X protocol requests */
static void GxErrorHandler(void)
{
  char msg[80];
  gxErrorFlag= 0;
  XErrHandler(msg);
}

int GpSetXHandler(void (*ErrHandler)(char *errMsg))
{
  /* install X error handlers which don't call exit */
  XErrHandler= ErrHandler;
  return 0;
}

/* ------------------------------------------------------------------------ */

int
g_rgb_read(Engine *eng, GpColor *rgb, long *nx, long *ny)
{
  XEngine *xeng = GisXEngine(eng);
  if (!xeng || !xeng->w || !xeng->win) return 1;
  if (!rgb) {
    *nx = xeng->wtop;
    *ny = xeng->htop;
  } else {
    p_rgb_read(xeng->win, rgb, xeng->leftMargin, xeng->topMargin,
               xeng->leftMargin+xeng->wtop, xeng->topMargin+xeng->htop);
  }
  return 0;
}

/* ------------------------------------------------------------------------ */
