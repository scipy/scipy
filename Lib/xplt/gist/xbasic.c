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
#include "xfont.h"
#include "gtext.h"
#include "dispax.h"

/* XUniqueContext is front for XrmUniqueQuark */
#include <X11/Xresource.h>

/* Memory manager functions (use long instead of indeterminate size_t) */
extern void *malloc(long);  /* required for image->data below */
extern char *strcpy(char *, const char *);

static char *xType= "X Window";

static int ChangePalette(Engine *engine);
static GpReal TextWidth(const char *text, int nc, const GpTextAttribs *t);
static void Kill(Engine *engine);
static void ShutDown(XEngine *xEngine);
static int Clear(Engine *engine, int always);
static int Flush(Engine *engine);
static void GetXRectangle(GpXYMap *map, GpBox *box, XRectangle *xrect);
static void ClearArea(Engine *engine, GpBox *box);
static void SetXTransform(GpTransform *trans, int landscape, int dpi);
static GpBox *DamageClip(GpBox *damage);
static void ChangeMap(Engine *engine);
static int SetupLine(XEngine *xEngine, Display *display, GC gc, int join,
		     GpLineAttribs *gistAl);
static int DrawLines(Engine *engine, long n, const GpReal *px,
		     const GpReal *py, int closed, int smooth);
static int DrawMarkers(Engine *engine, long n, const GpReal *px,
		       const GpReal *py);
static int DrawText(Engine *engine, GpReal x0, GpReal y0, const char *text);
static int DrawFill(Engine *engine, long n, const GpReal *px,
		    const GpReal *py);
static long GetCells(GpMap *map, GpReal px, GpReal dx,
		     GpReal xmin, GpReal xmax, long nmax, short *x, long *o);
static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
		     GpReal qy, long width, long height, long nColumns,
		     const GpColor *colors);
static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
			const GpReal *py, const GpReal *qx, const GpReal *qy);
extern int GxBasicXHandler(XEvent *event);
static void ClearPixmap(XEngine *xeng);
static void GetVisibleNDC(XEngine *xeng,
			  GpReal *xn, GpReal *xx, GpReal *yn, GpReal *yx);
static int YXError(Display *dpy, XErrorEvent *xev);
static int YXIOError(Display *dpy);

/* Hook for hlevel.c error handling.  */
extern void (*HLevelHook)(Engine *engine);

extern int GxJustifyText(GpXYMap *map, XFontStruct *font,
			 GpReal x0, GpReal y0, const char *text,
			 int *ix, int *iy);
extern int GxJustifyNext(XFontStruct *font, const char **text,
			 int *ix, int *iy);

static int gxErrorFlag= 0;
static void GxErrorHandler(void);

/* ------------------------------------------------------------------------ */

static int ChangePalette(Engine *engine)
{
  XEngine *xeng= (XEngine *)engine;
  GxScreen *xscr= xeng->xscr;
  GpColorCell *palette= engine->palette;  /* requested palette */
  int nColors= engine->nColors;
  int class, mapSize;

  if (!xscr) return 0;
  class= xscr->v->class;
  mapSize= class==DirectColor? xscr->dcinfo.rwMapSize :
    xscr->v->colormap_size;

  if (mapSize>256) mapSize= 256;
  if (nColors>256) nColors= 256;

  /* free any shared colors currently used by this XEngine */
  if (xeng->pixelMap && !xeng->private)
    XFreeColors(xscr->display, xscr->colormap,
		xeng->pixelMap, xeng->nColors, 0L);

  if (engine->colorMode) {
    /* get exact colors requested in palette */
    int wasPublic= (xeng->private==0);
    xeng->pixelMap= GxExactColors(xscr, palette, nColors, xeng->pixelMap,
				  &xeng->private, xeng->nColors);
    if (xeng->private && wasPublic)
      XSetWindowColormap(xscr->display, xeng->top, xeng->private);
  } else {
    /* get rid of any private colormap, and get shared colors */
    if (xeng->private) XFreeColormap(xscr->display, xeng->private);
    xeng->private= 0;
    xeng->pixelMap= GxShareColors(xscr, palette, nColors, xeng->pixelMap);
  }
  xeng->e.colorChange= 0;

  if (!xeng->pixelMap) {
    xeng->nColors= 0;
    return 0;
  }

  xeng->nColors= nColors;
  return mapSize;
}

/* ------------------------------------------------------------------------ */

/* These might be useful in future extensions.  */
extern void GxSetColor(XEngine *xEngine, int color);
extern XFontStruct *GxSetFont(GxScreen *xscr, int dpi, GC gc, int id);

static XFontStruct *currentFont, *symbol_font;
static int current_id, symbol_id;
static XEngine *font_engine;
static XFontStruct *GxToggleFont(int symbol, int set_gc);

void GxSetColor(XEngine *xEngine, int color)
{
  GxScreen *xscr= xEngine->xscr;
  GC gc= xEngine->gc;
  unsigned long pixel;
  int ncolors;
  if (!xscr) return;
  if (xEngine->e.colorChange) ChangePalette((Engine *)xEngine);
  ncolors= xEngine->nColors;

  if (color>=ncolors) pixel= xscr->stdColors[-1-FG_COLOR].pixel;
  else if (color>=0) {
    unsigned long *pixelMap= xEngine->pixelMap;
    pixel= pixelMap[color];
  } else if (color>=YELLOW_COLOR) pixel= xscr->stdColors[-1-color].pixel;
  else pixel= xscr->stdColors[-1-FG_COLOR].pixel;

  XSetForeground(xscr->display, gc, pixel);
}

XFontStruct *GxSetFont(GxScreen *xscr, int dpi, GC gc, int id)
{
  GxDisplay *xdpy= xscr->owner;
  XFontStruct *font, **loadedFonts= xdpy->loadedFonts;
  int *loadedID= xdpy->loadedID;
  int i, j;

  /* Check if this is permanently loaded font */
  if (id==PERM_FONT_ID) {
    /* Check permanently loaded font */
    if (xdpy->permFont) {
      if (gc) XSetFont(xscr->display, gc, xdpy->permFont->fid);
      return xdpy->permFont;
    }
  }

  /* Check if this is ordinary loaded font */
  font= 0;
  for (i=0 ; i<5 ; i++) {
    font= loadedFonts[i];
    if (!font) break;
    if (loadedID[i]==id) break;
  }

  /* Free space in font cache list if necessary */
  if (i>=5) {
    XFreeFont(xscr->display, font);
    i= 4;
    font= loadedFonts[i]= 0;
  } else if (!font && i<4) {
    loadedFonts[i+1]= 0;
  }

  /* Place this font first in font cache list */
  for (j=i ; j>0 ; j--) {
    loadedFonts[j]= loadedFonts[j-1];
    loadedID[j]= loadedID[j-1];
  }

  /* Load this font if necessary */
  if (!font) {
    char *name= GxNameFont(id);
    if (name) font= XLoadQueryFont(xscr->display, name);
    else font= xdpy->defaultFont; /* but no way to be sure this is set?? */
  }

  /* Set this font in GC and cache list */
  loadedFonts[0]= font;
  loadedID[0]= id;
  if (gc) XSetFont(xscr->display, gc, font->fid);
  return font;
}

static XFontStruct *GxToggleFont(int symbol, int set_gc)
{
  XFontStruct *xfont= symbol? symbol_font : currentFont;
  int font= symbol? T_SYMBOL : gistA.t.font;
  int xfontid;
  if (!xfont) {
    if (font_engine->lastOp.tfont!=font ||
	font_engine->lastOp.tsize!=gistA.t.height) {
      xfontid= GxIDFont(font_engine->xdpy, font_engine->dpi,
			font, gistA.t.height);
    } else {
      xfontid= font_engine->lastOp.fontID;
    }
    if (symbol) symbol_id= xfontid;
    else current_id= xfontid;
    xfont= GxSetFont(font_engine->xscr, font_engine->dpi, (GC)0, xfontid);
    if (symbol) symbol_font= xfont;
    else currentFont= xfont;
  } else {
    xfontid= symbol? symbol_id : current_id;
  }
  if (set_gc) {
    if (font_engine->lastOp.tfont!=font ||
	font_engine->lastOp.tsize!=gistA.t.height) {
      font_engine->lastOp.tfont= font;
      font_engine->lastOp.tsize= gistA.t.height;
      font_engine->lastOp.fontID= xfontid;
    }
    XSetFont(font_engine->xscr->display, font_engine->gc, xfont->fid);
  }
  return xfont;
}

/* ------------------------------------------------------------------------ */

static int nChars, lineHeight, textHeight, prevWidth, maxWidth, alignH;
static int firstTextLine= 0;

static int current_state;
static int chunkWidth, nChunk, supersub, dy_super, dy_sub, x_chunks;

/* ARGSUSED */
static GpReal TextWidth(const char *text, int nc, const GpTextAttribs *t)
{
  int width= 0;
  if (firstTextLine) nChars= nc;
  if (!gtDoEscapes) {
    width= XTextWidth(currentFont, text, nc);
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
	  width+= XTextWidth(currentFont, txt, (int)(text-txt));
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
	    if (!symbol_font) GxToggleFont(1, 0);
	    width+= XTextWidth(symbol_font, &c, 1);
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
      width+= XTextWidth(currentFont, txt, (int)(text-txt));
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

int GxJustifyText(GpXYMap *map, XFontStruct *font,
		  GpReal x0, GpReal y0, const char *text, int *ix, int *iy)
{
  int nLines, alignV, dx, dy, xmin, xmax, ymin, ymax;
  GpReal rwidest;

  /* Set nLines, maxWidth, nChars, prevWidth */
  firstTextLine= 1;
  currentFont= font;
  nChars= prevWidth= chunkWidth= nChunk= supersub= 0;
  nLines= GtTextShape(text, &gistA.t, &TextWidth, &rwidest);
  maxWidth= (int)rwidest;

  /* state bits:
     1 - this chunk is superscript
     2 - this chunk is subscript (1 and 2 together meaningless)
     4 - this chunk is single symbol char
   */
  current_state= 0;
  x_chunks= 0;

  /* Compute height of one line */
  textHeight= font->ascent+font->descent;
  dy_super= (supersub&1)? font->ascent/3 : 0;
  dy_sub= (supersub&2)? font->ascent/3 : 0;
  lineHeight= textHeight+dy_sub+dy_super;

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

  /* Note: ymax= ymin+nLines*lineHeight  and  ymin= dy-font->ascent  */
  if (alignV<=TV_CAP) {
    dy= font->ascent + dy_super;
    ymin= 0;
    ymax= nLines*lineHeight;
  } else if (alignV==TV_HALF) {
    ymin= font->ascent/2;
    ymax= nLines*lineHeight;
    dy= ymin-(ymax-lineHeight)/2;
    ymin= dy-font->ascent-dy_super;
    ymax+= ymin;
  } else if (alignV==TV_BASE) {
    dy= -(nLines-1)*lineHeight;
    ymin= dy-font->ascent-dy_super;
    ymax= font->descent + dy_sub;
  } else {
    ymin= dy_sub-nLines*lineHeight;
    dy= ymin+font->ascent+dy_super;
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

  /* Check to see whether string is clipped */
  /* Note: if window.?max<window.?min, there may be a 1 pixel error
     in this clipping criterion... */
  if ( ((x0+xmax/map->x.scale<gistT.window.xmin) ^
	(x0+xmin/map->x.scale>gistT.window.xmax)) ||
       ((y0+ymax/map->y.scale<gistT.window.ymin) ^
	(y0+ymin/map->y.scale>gistT.window.ymax)) ) return -1;

  /* Individual lines may still be clipped, but let XDrawString
     take care of that.  */
  *ix= dx + (map->x.scale*x0+map->x.offset);
  *iy= dy + (map->y.scale*y0+map->y.offset);

  if (nChunk) GxToggleFont(0, 1);
  return nChunk;
}

int GxJustifyNext(XFontStruct *font, const char **text, int *ix, int *iy)
{
  const char *txt= *text+nChunk;
  int xadj= 0, yadj= 0;
  int count;
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
      if (nChunk) chunkWidth= XTextWidth(currentFont, txt, nChunk);
      else chunkWidth= 0;  /* unused */
    }

    /* reset state, adjusting (possibly rotated) x and y as well */
    xadj-= x_chunks;
    yadj= lineHeight;
    if (current_state&1) yadj+= dy_super;
    else if (current_state&2) yadj-= dy_sub;
    if (nChunk && (current_state&4)) GxToggleFont(0, 1);
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
	GxToggleFont(0, 1);
	current_state&= 3;
      } else {
	/* chunk is single symbol char */
	GxToggleFont(1, 1);
	current_state|= 4;
      }

    } else {
      for (nChunk=0 ; nChunk<nChars ; nChunk++) {
	c= txt[nChunk];
	if ((nChunk+1<nChars && c=='!') || c=='^' || c=='_') break;
      }
      if (nChunk) GxToggleFont(0, 1);
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
      chunkWidth= XTextWidth((current_state&4)? symbol_font : currentFont,
			     txt, nChunk);
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

static void Kill(Engine *engine)
{
  XEngine *xEngine= (XEngine *)engine;
  ShutDown(xEngine);
  GpDelEngine(engine);
}

static int shutDownPending= 0;
static Display *shutDownDisplay= 0;
static XEngine *shutDownEngine= 0;
void (*HLevelHook)(Engine *engine)= 0;

static void ShutDown(XEngine *xEngine)
{
  XEngine *xeng= xEngine;
  GxScreen *xscr;
  Display *display;
  xEngine->mapped= 0;   /* stop trying to draw to this turkey -- important
			   if this shutdown procedure fails */
  /* destroy any hlevel references without further ado */
  if (HLevelHook) HLevelHook((Engine *)xEngine);
  if (shutDownPending && shutDownEngine!=xEngine) xEngine= shutDownEngine;
  while (xEngine) {
    xscr= xEngine->xscr;
    if (!xscr) break;
    display= xscr->display;
    if (shutDownPending<1) {
      shutDownEngine= xEngine;
      shutDownDisplay= display;
      shutDownPending= 1;
    }
    if (xEngine->drawable!=xEngine->graphics && shutDownPending<3) {
      if (shutDownPending<2) {
	shutDownPending= 2;
	XFreePixmap(display, xEngine->drawable);
      }
      shutDownPending= 3;
      XFreeGC(display, xEngine->gca);
    }
    if (xEngine->pixelMap && shutDownPending<5) {
      if (!xEngine->private && shutDownPending<4) {
	shutDownPending= 4;
	XFreeColors(xscr->display, xscr->colormap,
		    xEngine->pixelMap, xEngine->nColors, 0L);
      }
      shutDownPending= 5;
      GmFree(xEngine->pixelMap);
      xEngine->pixelMap= 0;
    }
    if (xEngine->private && shutDownPending<6) {
      shutDownPending= 6;
      XFreeColormap(display, xEngine->private);
    }
    if (shutDownPending<7) {
      shutDownPending= 7;
      XDestroyWindow(display, xEngine->top);
    }
    if (shutDownPending<8) {
      shutDownPending= 8;
      XFreeGC(display, xEngine->gc);
    }
    if (shutDownPending<9) {
      shutDownPending= 9;
      xEngine->top= None;
      /* destroy any X resources associated with fancy X engines */
      if (xEngine->e.Kill!=&Kill) xEngine->e.Kill(&xEngine->e);
    }
    if (shutDownPending<10) {
      shutDownPending= 10;
      if (xscr->owner->references<=0) RemoveXDispatcher(display);
      GxDisconnect(xscr);
    }
    xEngine->xscr= 0;
    shutDownEngine= 0;
    shutDownDisplay= 0;
    shutDownPending= 0;
    if (xEngine!=xeng) xEngine= xeng;
    else xEngine= 0;
  }
}

static int Clear(Engine *engine, int always)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  if (!xscr) return 1;
  if ((always || xEngine->e.marked) && xEngine->graphics==xEngine->drawable)
    XClearWindow(xscr->display, xEngine->graphics);
  if (xEngine->e.colorChange) ChangePalette(engine);
  xEngine->e.marked= 0;
  return 0;
}

static int Flush(Engine *engine)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  if (!xscr) return 1;
  XFlush(xscr->display);
  /* test whether an X error has been reported */
  if (gxErrorFlag) GxErrorHandler();
  return 0;
}

static void GetXRectangle(GpXYMap *map, GpBox *box, XRectangle *xrect)
{
  short wmin, wmax;

  /* Get corners of clip rectangle in pixel coordinates */
  wmin= (short)(map->x.scale*box->xmin+map->x.offset);
  wmax= (short)(map->x.scale*box->xmax+map->x.offset);
  if (wmax>=wmin) {
    xrect->x= (short)wmin;
    xrect->width= ((short)wmax)-xrect->x+1;
  } else {
    xrect->x= (short)wmax;
    xrect->width= ((short)wmin)-xrect->x+1;
  }
  wmin= (short)(map->y.scale*box->ymin+map->y.offset);
  wmax= (short)(map->y.scale*box->ymax+map->y.offset);
  if (wmax>=wmin) {
    xrect->y= (short)wmin;
    xrect->height= ((short)wmax)-xrect->y+1;
  } else {
    xrect->y= (short)wmax;
    xrect->height= ((short)wmin)-xrect->y+1;
  }
}

static void ClearArea(Engine *engine, GpBox *box)
{
  XEngine *xEngine= (XEngine *)engine;
  XRectangle r;
  GxScreen *xscr= xEngine->xscr;
  if (!xscr) return;
  /* If this is animation mode, do not try to clear window */
  if (xEngine->drawable!=xEngine->graphics) return;
  GetXRectangle(&engine->devMap, box, &r);
  XClearArea(xEngine->xscr->display, xEngine->graphics,
	     r.x, r.y, r.width, r.height, False);
}

static void SetXTransform(GpTransform *trans, int landscape, int dpi)
{
  trans->viewport= landscape? gLandscape : gPortrait;
  trans->window.xmin= 0.0;
  trans->window.xmax= PixelsPerNDC(dpi)*trans->viewport.xmax;
  trans->window.ymin= PixelsPerNDC(dpi)*trans->viewport.ymax;
  trans->window.ymax= 0.0;
}

static GpBox cPort;

static GpBox *DamageClip(GpBox *damage)
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

static void ChangeMap(Engine *engine)
{
  XEngine *xEngine= (XEngine *)engine;
  XRectangle xrect;
  GpBox *clipport;
  int landscape= xEngine->width>xEngine->height;
  GxScreen *xscr= xEngine->xscr;
  if (!xscr) return;

  /* Check to be sure that landscape/portrait mode hasn't changed */
  if (landscape!=engine->landscape) {
    Window root;
    int x, y;
    unsigned int width, height, border, depth;
    XWindowChanges changes;
    SetXTransform(&engine->transform, engine->landscape, xEngine->dpi);
    GpDeviceMap(engine);
    xEngine->swapped= xEngine->e.transform;
    changes.width= xEngine->width= (int)engine->transform.window.xmax;
    changes.height= xEngine->height= (int)engine->transform.window.ymin;
    /* XConfigureWindow will NOT generate a ConfigureNotify event, because
       StructureNotifyMask is not set for the graphics window.  Therefore,
       call GxRecenter explicitly here.  */
    XConfigureWindow(xscr->display, xEngine->graphics,
		     CWWidth | CWHeight, &changes);
    if (!xscr ||
	!XGetGeometry(xscr->display, xEngine->top, &root, &x, &y,
		      &width, &height, &border, &depth)) {
      /* Use default window size on failure.  */
      width= DefaultTopWidth(xEngine->dpi);
      height= DefaultTopHeight(xEngine->dpi);
    }
    GxRecenter(xEngine, width, height);
  }

  /* Do generic change map */
  GpComposeMap(engine);

  /* Get XRectangle for current clip window.  */
  if (engine->damaged) clipport= DamageClip(&engine->damage);
  else clipport= &gistT.viewport;
  if (!clipport) return;

  /* Set clipping rectangle in gc for this XEngine */
  GetXRectangle(&engine->devMap, clipport, &xrect);
  XSetClipRectangles(xscr->display, xEngine->gc, 0, 0,
		     &xrect, 1, YXBanded);
}

/* ------------------------------------------------------------------------ */

static char dashed[]= { 5, 5 };
static char dotted[]= { 1, 3 };
static char dashdot[]= { 5, 2, 1, 2 };
static char dashdotdot[]= { 5, 2, 1, 2, 1, 2 };
char *gxDash[]= { dashed, dotted, dashdot, dashdotdot };
int gxDashN[]= { 2, 2, 4, 6 };

static int SetupLine(XEngine *xEngine, Display *display, GC gc, int join,
		     GpLineAttribs *gistAl)
{
  int widthChange= (xEngine->lastOp.lwidth!=gistAl->width);
  int typeChange= (xEngine->lastOp.ltype!=gistAl->type);
  int joinChange= (xEngine->lastOp.ljoin!=join);

  if (widthChange || typeChange || joinChange) {
    int lstyle;
    unsigned int lwidth;
    if (gistAl->type==L_NONE) return 1;
    lwidth= (unsigned int)(DEFAULT_LINE_WIDTH*
      gistAl->width*(xEngine->dpi==75?798.289:1064.385));
    lstyle= gistAl->type==L_SOLID? LineSolid : LineOnOffDash;
    xEngine->lastOp.lwidth= gistAl->width;
    xEngine->lastOp.ltype= gistAl->type;
    xEngine->lastOp.ljoin= 0;
    if (gistAl->type!=L_SOLID) {
      int iDash= gistAl->type-2;
      int dashN= gxDashN[iDash];
      char *dash, dashS[6];
      if (lwidth<2) {
	dash= gxDash[iDash];
      } else {  /* dash pattern must scale with line width */
	int i;
	for (i=0 ; i<dashN ; i++)
	  dashS[i]= gxDash[iDash][i]>1? lwidth*gxDash[iDash][i] : 1;
	dash= dashS;
      }
      XSetDashes(display, gc, 0, dash, dashN);
    }
    XSetLineAttributes(display, gc, lwidth, lstyle,
		       join? CapProjecting : CapRound,
		       join? JoinMiter : JoinRound);
  }
  if (xEngine->lastOp.color!=gistAl->color) {
    xEngine->lastOp.color= gistAl->color;
    GxSetColor(xEngine, gistAl->color);
  }

  return 0;
}

static int DrawLines(Engine *engine, long n, const GpReal *px,
		     const GpReal *py, int closed, int smooth)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  Drawable window= xEngine->drawable;
  Display *display= xscr? xscr->display : 0;
  GC gc= xEngine->gc;
  XPoint *points;
  GpXYMap *map= &xEngine->e.map;
  long maxPoints, nPoints;
  int firstPass= 1;
  XPoint firstPoint;

  if (!xscr || !xEngine->mapped) return 1;
  if (SetupLine(xEngine, display, gc, 0, &gistA.l)) return 0;

  maxPoints= xscr->maxRequest-4; /* allow for closure pt */
  while ((nPoints=
	  GpIntPoints(map, maxPoints, n, px, py, (GpPoint**)&points))) {
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
    XDrawLines(display, window, gc, points, (int)nPoints, CoordModeOrigin);
    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
  }

  xEngine->e.marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

/* X11 has no polymarker primitive, use GpPseudoMark-- unless
   user requested tiny points, in which case we use XDrawPoints */
static int DrawMarkers(Engine *engine, long n, const GpReal *px,
		       const GpReal *py)
{
  XEngine *xEngine= (XEngine *)engine;
  if (!xEngine->xscr || !xEngine->mapped) return 1;

  xEngine->e.marked= 1;
  if (gistA.m.type!=M_POINT || gistA.m.size>1.5) {
    return GpPseudoMark(engine, n, px, py);

  } else {
    GxScreen *xscr= xEngine->xscr;
    Drawable window= xEngine->drawable;
    Display *display= xscr->display;
    GC gc= xEngine->gc;
    XPoint *points;
    GpXYMap *map= &xEngine->e.map;
    long maxPoints, nPoints;

    if (xEngine->lastOp.color!=gistA.m.color) {
      xEngine->lastOp.color= gistA.m.color;
      GxSetColor(xEngine, gistA.m.color);
    }

    maxPoints= xscr->maxRequest-3;
    while ((nPoints=
	    GpIntPoints(map, maxPoints, n, px, py, (GpPoint**)&points))) {
      XDrawPoints(display, window, gc, points, (int)nPoints, CoordModeOrigin);
      if (n==nPoints) break;
      n-= nPoints;
      px+= nPoints;
      py+= nPoints;
    }

    return 0;
  }
}

/* ------------------------------------------------------------------------ */

static XImage *nimage=0, *rimage=0;

static int DrawText(Engine *engine, GpReal x0, GpReal y0, const char *text)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  Drawable window= xEngine->drawable;
  Display *display= xscr? xscr->display : 0;
  GC gc= xEngine->gc;
  int ix, iy, len;
  XFontStruct *font;
  Pixmap npixmap=None, rpixmap=None;
  const char *txt;
  char caret= '^';

  if (!xscr || !xEngine->mapped) return 1;

  /* There is no reasonable way to support text at angles under X.
     O'Reilly suggests that the "correct" way to get letters that
     read from the side is to design your own fonts with prone letters,
     then install them on all servers...
     Therefore, we adopt the unreasonable technique of transposing
     an image of a text string one pixel at a time.  The idea is to
     create a pixmap on the server, draw the string to it, read it back
     from the server as an image, destroy the pixmap on the server,
     create a second image, perform the rotation, destroy the first
     image, copy the rotated image to the final destination, destroy
     the rotated image.  One drawback is that only opaque=1 text can
     be correctly drawn by this algorithm -- it turns out that it is
     possible for XGetImage to be unable to return the current contents
     of the destination window (if it is obscured), so there is no way
     to get opaque=0 rotated text.  */
  /* cleanup after possible interrupt of a previous operation */
  if (nimage) {
    XImage *im= nimage;
    nimage= 0;
    XDestroyImage(im);
  }
  if (rimage) {
    XImage *im= rimage;
    rimage= 0;
    XDestroyImage(im);
  }
  if (xscr->rot_normal!=None) {
    npixmap= xscr->rot_normal;
    xscr->rot_normal= None;
    XFreePixmap(display, npixmap);
  }
  if (xscr->rot_rotated!=None) {
    rpixmap= xscr->rot_rotated;
    xscr->rot_rotated= None;
    XFreePixmap(display, rpixmap);
  }

  if (xEngine->lastOp.color!=gistA.t.color) {
    xEngine->lastOp.color= gistA.t.color;
    GxSetColor(xEngine, gistA.t.color);
  }

  /* handle multi-line strings */
  font_engine= xEngine;
  currentFont= symbol_font= 0;
  font= GxToggleFont(0, 0);
  len= GxJustifyText(&xEngine->e.map, font, x0, y0, text, &ix, &iy);

  if (gistA.t.orient!=TX_RIGHT) {
    npixmap= XCreatePixmap(display, xscr->root, maxWidth, textHeight,
			   DefaultDepth(display, xscr->screen));
    xscr->rot_normal= npixmap;
    if (gistA.t.orient==TX_LEFT) {
      rpixmap= XCreatePixmap(display, xscr->root, maxWidth, textHeight,
			     DefaultDepth(display, xscr->screen));
      xscr->rot_rotated= rpixmap;
      rimage= XGetImage(display, rpixmap, 0, 0, maxWidth, textHeight,
			AllPlanes, XYPixmap);
    } else {
      rpixmap= XCreatePixmap(display, xscr->root, textHeight, maxWidth,
			     DefaultDepth(display, xscr->screen));
      xscr->rot_rotated= rpixmap;
      rimage= XGetImage(display, rpixmap, 0, 0, textHeight, maxWidth,
			AllPlanes, XYPixmap);
    }
  }

  while (len>=0) {
    if (len>0) {
      if (len==1 && (current_state&4) && text[0]==']') txt= &caret;
      else txt= text;
      if (gistA.t.orient==TX_RIGHT) {
	if (gistA.t.opaque)
	  XDrawImageString(display, window, gc, ix, iy, txt, len);
	else
	  XDrawString(display, window, gc, ix, iy, txt, len);
      } else if (chunkWidth>0) {
	/* rotated text is a terrible mess
	   -- this kludge does not support non-opaque rotated text */
	XImage *im;
	int i, j;
	XDrawImageString(display, npixmap, gc, 0, font->ascent, txt, len);
	nimage= XGetImage(display, npixmap, 0, 0, chunkWidth, textHeight,
			  AllPlanes, XYPixmap);
	if (gistA.t.orient==TX_LEFT) {
	  for (j=0 ; j<textHeight ; ++j)
	    for (i=0 ; i<chunkWidth ; ++i)
	      XPutPixel(rimage, chunkWidth-1-i, textHeight-1-j,
			XGetPixel(nimage,i,j));
	  XPutImage(display, rpixmap, gc, rimage, 0,0,
		    0,0,chunkWidth,textHeight);
	  XCopyArea(display, rpixmap, window, gc, 0,0,chunkWidth,textHeight,
		    ix-chunkWidth,iy-font->descent);
	} else if (gistA.t.orient==TX_UP) {
	  for (j=0 ; j<textHeight ; ++j)
	    for (i=0 ; i<chunkWidth ; ++i)
	      XPutPixel(rimage, j, chunkWidth-1-i, XGetPixel(nimage,i,j));
	  XPutImage(display, rpixmap, gc, rimage, 0,0,
		    0,0,textHeight,chunkWidth);
	  XCopyArea(display, rpixmap, window, gc, 0,0,textHeight,chunkWidth,
		    ix-font->ascent,iy-chunkWidth);
	} else {
	  for (j=0 ; j<textHeight ; ++j)
	    for (i=0 ; i<chunkWidth ; ++i)
	      XPutPixel(rimage, textHeight-1-j, i, XGetPixel(nimage,i,j));
	  XPutImage(display, rpixmap, gc, rimage, 0,0,
		    0,0,textHeight,chunkWidth);
	  XCopyArea(display, rpixmap, window, gc, 0,0,textHeight,chunkWidth,
		    ix-font->descent,iy);
	}
	im= nimage;
	nimage= 0;
	XDestroyImage(im);
      }
    }
    len= GxJustifyNext(font, &text, &ix, &iy);
  }

  if (rimage) {
    XImage *im= rimage;
    xscr->rot_normal= None;
    XFreePixmap(display, npixmap);
    xscr->rot_rotated= None;
    XFreePixmap(display, rpixmap);
    rimage= 0;
    XDestroyImage(im);
  }

  xEngine->e.marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawFill(Engine *engine, long n, const GpReal *px,
		    const GpReal *py)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  Drawable window= xEngine->drawable;
  Display *display= xscr? xscr->display : 0;
  GC gc= xEngine->gc;
  XPoint *points;
  GpXYMap *map= &xEngine->e.map;
  long maxPoints, nPoints;
  int value= 0;

  if (!xscr || !xEngine->mapped) return 1;

  /* For now, only FillSolid style supported */

  if (xEngine->lastOp.color!=gistA.f.color) {
    xEngine->lastOp.color= gistA.f.color;
    GxSetColor(xEngine, gistA.f.color);
  }

  /* This give incorrect results if more than one pass through the loop,
     but there is no way to give the correct result, so may as well...  */
  maxPoints= (xscr->maxRequest-3)/2;
  while ((nPoints=
	  GpIntPoints(map, maxPoints, n, px, py, (GpPoint**)&points))) {
    /* Can Nonconvex or Convex be detected? */
    if (nPoints==4 &&
	((points[0].x==points[3].x && points[0].y==points[1].y &&
	  points[1].x==points[2].x && points[2].y==points[3].y) ||
	 (points[0].y==points[3].y && points[0].x==points[1].x &&
	  points[1].y==points[2].y && points[2].x==points[3].x))) {
      /* Detecting rectangles saves a huge amount of drawing time for
         the case of a rectangular mesh (GaFillMesh being used as a cell
	 array).  */
      int x, y, dx, dy;
      dx= points[1].x-points[0].x;
      if (dx!=0) {
	dy= points[3].y-points[0].y;
	if (dx>=0) x= points[0].x;
	else { dx= -dx;  x= points[1].x; }
	if (dy>=0) y= points[0].y;
	else { dy= -dy;  y= points[3].y; }
      } else {
	dx= points[3].x-points[0].x;
	dy= points[1].y-points[0].y;
	if (dx>=0) x= points[0].x;
	else { dx= -dx;  x= points[3].x; }
	if (dy>=0) y= points[0].y;
	else { dy= -dy;  y= points[1].y; }
      }
      XFillRectangle(display, window, gc, x, y,
		     (unsigned int)dx, (unsigned int)dy);
    } else {
      /* Do general case */
      XFillPolygon(display, window, gc, points, (int)nPoints,
		   Complex, CoordModeOrigin);
    }
    if (gistA.e.type!=L_NONE) {
      if (SetupLine(xEngine, display, gc, 0, &gistA.e)) return 1;
      n++;
      points[nPoints++]= points[0];
      XDrawLines(display, window, gc, points, (int)nPoints, CoordModeOrigin);
    }
    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
    value= 1;
  }

  xEngine->e.marked= 1;
  return value;
}

/* ------------------------------------------------------------------------ */

/* NB-
   The closest thing to a cell array primitive in X windows is XPutImage.
   However, it is not clear exactly how this is supposed to be used to
   do simple pseudo-color image rendering.  The problem is that, as far
   as I can tell, you need special coding for each possible screen
   depth.  The comments in the Xlib source indicate that depths of
   1, 4, 8, 12, 16, 20, 24, and 32 are known.

   Actually, if there are many pixels per cell, it will surely be faster
   to use XFillRectangle instead of XPutImage.  Assume each cell
   requires a different color from the preceding cell, then:
     1) PutImage request requires 6 + nPixels/(pixels per 4-byte word)
        words
     2) ChangeGC request to change color takes 4 words,
        PolyFillRectangle request takes 5 words, for a total of
	9 words per cell
   Hence, for 1-bit deep displays, the breakeven point is at 288 pixels
   per cell; for 8-bit deep displays it is at 36 pixels per cell; for
   32-bit deep displays it is at 9 pixels per cell.
 */

static long GetCells(GpMap *map, GpReal px, GpReal dx,
		     GpReal xmin, GpReal xmax, long nmax, short *x, long *o)
{
  long j, i= 0;
  GpReal scale= map->scale;
  GpReal offset= map->offset;

  xmin= xmin*scale+offset;
  xmax= xmax*scale+offset;
  if (xmin>xmax) {GpReal tmp=xmin; xmin=xmax; xmax=tmp;}
  px= px*scale+offset;
  dx*= scale;

  if (dx>0.0) {
    for (j=0 ; j<nmax && px<xmin ; j++) px+= dx;
    *o= j>0? j-1 : 0;
    if (j>0 && j<nmax) x[i++]= (short)(xmin);
    for ( ; j<nmax && px<=xmax ; j++) {
      x[i++]= (short)(px);
      px+= dx;
    }
    if (j>0 && j<nmax) x[i++]= (short)(xmax);
  } else if (dx<0.0) {
    for (j=0 ; j<nmax && px>xmax ; j++) px+= dx;
    *o= j>0? j-1 : 0;
    if (j>0 && j<nmax) x[i++]= (short)(xmax);
    for ( ; j<nmax && px>=xmin ; j++) {
      x[i++]= (short)(px);
      px+= dx;
    }
    if (j>0 && j<nmax) x[i++]= (short)(xmin);
  }
  return i;
}

static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
		     GpReal qy, long width, long height, long nColumns,
		     const GpColor *colors)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  Display *display= xscr? xscr->display : 0;
  Drawable window= xEngine->drawable;
  GC gc= xEngine->gc;
  GpXYMap *map= &xEngine->e.map;
  /* Create image now to get correct bits_per_pixel */
  XImage *image= xscr? XCreateImage(display, xscr->v->visual, xscr->v->depth,
				    ZPixmap, 0, 0, 1, 1,
				    BitmapPad(display), 0) : 0;
  short *x= (short *)GmMalloc(sizeof(short)*(width+height+2));
  short *y= x+width+1;
  long i, j, k, nx, ny, ox, oy;
  unsigned long *pixelMap= xEngine->pixelMap;
  int xpix, ypix;
  int ip, jp, xn, xx, yn, yx, dxn, dxx, dyn, dyx, x0, y0;

  if (!xEngine->mapped || !pixelMap || !image || !x) {
    if (image) XDestroyImage(image);
    if (x) GmFree(x);
    return 1;
  }

  /* Find cell boundaries (clip if necessary) */
  nx= GetCells(&map->x, px, (qx-px)/(GpReal)width,
	       gistT.window.xmin, gistT.window.xmax, width+1, x, &ox);
  ny= GetCells(&map->y, py, (qy-py)/(GpReal)height,
	       gistT.window.ymin, gistT.window.ymax, height+1, y, &oy);
  if (nx<2 || ny<2) {
    XDestroyImage(image);
    GmFree(x);
    return 0;
  }

  xpix= x[nx-1]-x[0];
  if (xpix<0) {
    xpix= -xpix;
    x0= x[nx-1];
    x[0]++;     /* largest value doesn't get filled */
    dxn= 1;
    dxx= 0;
  } else {
    x0= x[0];
    x[nx-1]++;  /* largest value doesn't get filled */
    dxn= 0;
    dxx= 1;
  }
  ypix= y[ny-1]-y[0];
  if (ypix<0) {
    ypix= -ypix;
    y0= y[ny-1];
    y[0]++;     /* largest value doesn't get filled */
    dyn= 1;
    dyx= 0;
  } else {
    y0= y[0];
    y[ny-1]++;  /* largest value doesn't get filled */
    dyn= 0;
    dyx= 1;
  }

  if ((long)xpix*(long)ypix*image->bits_per_pixel > 288*nx*ny) {
    /* Less X protocol traffic in series of filled rectangles */
    unsigned int w, h;
    int prev= xEngine->lastOp.color;
    k= ox+nColumns*oy;
    for (j=0 ; j<ny-1 ; j++) {
      yn= y[j+dyn];
      h= y[j+dyx]-yn;
      for (i=0 ; i<nx-1 ; i++) {
	xn= x[i+dxn];
	w= x[i+dxx]-xn;
	if (colors[k+i]!=prev) {
	  prev= colors[k+i];
	  XSetForeground(display, gc, pixelMap[prev]);
	}
	XFillRectangle(display, window, gc, xn, yn, w, h);
      }
      k+= nColumns;
    }
    xEngine->lastOp.color= prev;

  } else {
    /* Less X protocol traffic with single XPutImage */
    char *data;
    int bytesPerLine;
    
    bytesPerLine= ((xpix/image->bitmap_pad + 1)*
		   image->bitmap_pad*image->bits_per_pixel)/8;
    xpix++;
    ypix++;
    image->width= xpix;
    image->height= ypix;
    image->bytes_per_line= bytesPerLine;
    /* This is NOT GmMalloc, since free is called by XDestroyImage.  */
    image->data= data= (char *)malloc((long)(bytesPerLine*ypix));
    if (!data) {
      XDestroyImage(image);
      GmFree(x);
      return 1;
    }

    for (i=0 ; i<nx ; i++) x[i]-= x0;
    for (i=0 ; i<ny ; i++) y[i]-= y0;

    if (image->bits_per_pixel==8) {
      /* Easy to fill up pixels */
      char cval;
      k= ox+nColumns*oy;
      for (j=0 ; j<ny-1 ; j++) {
	yn= y[j+dyn];
	yx= y[j+dyx];
	for (i=0 ; i<nx-1 ; i++) {
	  cval= (char)pixelMap[colors[k+i]];
	  xn= x[i+dxn];
	  xx= x[i+dxx];
	  for (jp=yn ; jp<yx ; jp++)
	    for (ip=xn ; ip<xx ; ip++) data[jp*bytesPerLine+ip]= cval;
	}
	k+= nColumns;
      }
    } else {
      /* Use generic XPutPixel */
      unsigned long pixel;
      k= ox+nColumns*oy;
      for (j=0 ; j<ny-1 ; j++) {
	yn= y[j+dyn];
	yx= y[j+dyx];
	for (i=0 ; i<nx-1 ; i++) {
	  pixel= pixelMap[colors[k+i]];
	  xn= x[i+dxn];
	  xx= x[i+dxx];
	  for (jp=yn ; jp<yx ; jp++)
	    for (ip=xn ; ip<xx ; ip++) XPutPixel(image, ip, jp, pixel);
	}
	k+= nColumns;
      }
    }

    XPutImage(display, window, gc, image, 0, 0, x0, y0, xpix, ypix);
    /* free(image->data); done in XDestroyImage */
  }

  XDestroyImage(image);
  GmFree(x);

  xEngine->e.marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
			const GpReal *py, const GpReal *qx, const GpReal *qy)
{
  XEngine *xEngine= (XEngine *)engine;
  GxScreen *xscr= xEngine->xscr;
  Drawable window= xEngine->drawable;
  Display *display= xscr? xscr->display : 0;
  GC gc= xEngine->gc;
  XSegment *segs;
  GpXYMap *map= &xEngine->e.map;
  long maxSegs, nSegs;

  if (!xscr || !xEngine->mapped) return 1;
  if (SetupLine(xEngine, display, gc, 1, &gistA.l)) return 0;

  maxSegs= (xscr->maxRequest-3)/2;
  while ((nSegs=
	  GpIntSegs(map, maxSegs, n, px, py, qx, qy, (GpSegment**)&segs))) {
    XDrawSegments(display, window, gc, segs, (int)nSegs);
    if (n==nSegs) break;
    n-= nSegs;
    px+= nSegs;
    py+= nSegs;
    qx+= nSegs;
    qy+= nSegs;
  }

  xEngine->e.marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static XContext gxEngine;
static int gxEngineSet= 0;

XEngine *GxGetEngine(Display *display, Window window)
{
  Window orig= window;
  Window root, parent, *children;
  unsigned int nChildren;
  XEngine *xeng;
  while (XFindContext(display, window, gxEngine, (caddr_t *)&xeng)) {
    xeng= 0;
    if (!XQueryTree(display, window,
		    &root, &parent, &children, &nChildren)) break;
    XFree((caddr_t)children);
    if (root==parent || root==window) break;
    window= parent;
  }
  if (xeng && orig!=window)
    XSaveContext(display, orig, gxEngine, (void *)xeng);
  return xeng;
}

int GxBasicXHandler(XEvent *event)
{
  Display *dpy= event->xany.display;
  switch (event->type) {
  case Expose:
    /* redraw entire window */
    { XEngine *xEngine= GxGetEngine(dpy, event->xexpose.window);
      if (xEngine) {
	xEngine->mapped= 1;
	if (xEngine->HandleExpose)
	  /* the alternate handler should probably call GxExpose */
	  xEngine->HandleExpose(&xEngine->e, xEngine->e.drawing, event);
	else
	  GxExpose(&xEngine->e, xEngine->e.drawing, event);
      }
    }
    break;
  case UnmapNotify:
    /* reset mapped flag */
    { XEngine *xEngine= GxGetEngine(dpy, event->xunmap.window);
      if (xEngine) xEngine->mapped= 0;
    }
    break;
  case ConfigureNotify:
    /* Recenter graphics window within top level window */
    { XEngine *xEngine= GxGetEngine(dpy, event->xconfigure.window);
      if (xEngine) {
	if (xEngine->HandleResize)
	  xEngine->HandleResize(&xEngine->e, xEngine->e.drawing, event);
	else
	  GxRecenter(xEngine, event->xconfigure.width,
		     event->xconfigure.height);
      }
    }
    break;
  case ClientMessage:
    { XEngine *xEngine= GxGetEngine(dpy, event->xclient.window);
      /* the only message expected is in response to WM_DELETE_WINDOW
	 protocol -- to avoid having to retrieve the associated atom,
	 just shut down the engine without asking any questions

	 Sigh.  Some window managers apparently send unsolicited
	 ClientMessage events, so we have no choice but to check to
	 be sure this is really a WM_DELETE_WINDOW message.  */
      if (xEngine &&
	  event->xclient.format==32 &&
	  event->xclient.message_type==XInternAtom(dpy,"WM_PROTOCOLS",1) &&
	  event->xclient.data.l[0]==XInternAtom(dpy,"WM_DELETE_WINDOW",1)) {
	ShutDown(xEngine);
      }
    }
    break;
  default:
    /* ignore anything else unless there is a handler */
    { XEngine *xEngine= GxGetEngine(dpy, event->xany.window);
      if (xEngine) {
	int type= event->type;
	if (!xEngine->mapped) {
	  /* Sometimes a server or window manager bug (? ? ?) can prevent
	     the xEngine from ever seeing its expose event.  If it gets
	     a button press or keypress, we'll take that as evidence that
	     it is, in fact, exposed.  The idea is that if any events at
	     all are being delivered to the window, a click or a keypress
	     will "jolt" it back from zombiehood.  */
	  if (type==KeyPress || type==ButtonPress || type==MotionNotify)
	    xEngine->mapped= 1;
	}
	if (type==DestroyNotify)
	  XDeleteContext(dpy, event->xany.window, gxEngine);
	if (xEngine->HandleOther)
	  xEngine->HandleOther(&xEngine->e, xEngine->e.drawing, event);
      }
    }
    break;
  }

  /* test whether an X error has been reported */
  if (gxErrorFlag) GxErrorHandler();
  return 0;
}

void GxExpose(Engine *engine, Drawing *drawing, XEvent *event)
{
  GpBox damage;
  GpXYMap *map= &engine->devMap;
  XExposeEvent *xev= &event->xexpose;
  XEngine *xEngine= (XEngine *)engine;
  if (!drawing || event->xexpose.window!=xEngine->graphics) return;
  damage.xmin= (xev->x-map->x.offset)/map->x.scale;
  damage.xmax= (xev->x+xev->width-map->x.offset)/map->x.scale;
  damage.ymax= (xev->y-map->y.offset)/map->y.scale;
  damage.ymin= (xev->y+xev->height-map->y.offset)/map->y.scale;
  if (engine->damaged) {
    GpSwallow(&engine->damage, &damage);
  } else {
    engine->damage= damage;
    engine->damaged= 1;
  }
  if (event->xexpose.count != 0) return;
  GdSetDrawing(drawing);
  GpPreempt(engine);
  GdDraw(1);
  GpPreempt(0);        /* not correct if damaged during a preempt... */
  GdSetDrawing(0);
}

void GxRecenter(XEngine *xEngine, int width, int height)
{
  int x, y;
  long eWidth= (int)xEngine->width;    /* unsigned int deadly... */
  long eHeight= (int)xEngine->height;
  width-= xEngine->leftMargin;
  height-= xEngine->topMargin;
  x= (eWidth-width)/2;
  /* put center of page at center of landscape window */
  if (eWidth>eHeight) y= (eHeight-height)/2;
  /* put center of upper square of page at center of portrait window */
  else y= (eWidth-height)/2;
  /* once either dimension is big enough for whole picture, stop moving it */
  if (y<0) y= 0;
  if (x<0) x= 0;
  if (x!=xEngine->x || y!=xEngine->y) {
    XWindowChanges changes;
    changes.x= -x + xEngine->leftMargin -4/*border*/;
    changes.y= -y + xEngine->topMargin -4/*border*/;
    XConfigureWindow(xEngine->xscr->display, xEngine->graphics,
		     CWX | CWY, &changes);
    xEngine->x= x;
    xEngine->y= y;
  }
}

/* ------------------------------------------------------------------------ */

XGCValues gxXGCValues;

XEngine *GxEngine(char *name, GpTransform *toPixels, GxScreen *xscr, 
		  Window top, int x, int y,
		  int topMargin, int leftMargin, long engineSize)
{
  Display *display= xscr? xscr->display : 0;
  XEngine *xEngine;
  unsigned long bg= xscr? xscr->stdColors[0].pixel : 0;
  Window graphics;
  unsigned int width, height;
  GpReal rdpi;
  int dpi;
  XSetWindowAttributes cwa;

  if (!xscr) return 0;

  /* Graphics window will have dimensions of toPixels transform window */
  if (toPixels->window.xmin<toPixels->window.xmax)
    width= (unsigned int)(toPixels->window.xmax - toPixels->window.xmin);
  else
    width= (unsigned int)(toPixels->window.xmin - toPixels->window.xmax);
  if (toPixels->window.ymin<toPixels->window.ymax)
    height= (unsigned int)(toPixels->window.ymax - toPixels->window.ymin);
  else
    height= (unsigned int)(toPixels->window.ymin - toPixels->window.ymax);

  /* Reconstruct dpi (dots per inch) from toPixels transform */
  rdpi= (toPixels->window.ymax - toPixels->window.ymin) /
        (toPixels->viewport.ymin - toPixels->viewport.ymax);
  if (rdpi<0.0) rdpi= -rdpi;  /* almost certainly a bug, though... */
  rdpi*= ONE_INCH;
  dpi= rdpi<87.5? 75 : 100;

  xEngine=
    (XEngine *)GpNewEngine(engineSize, name, xType, toPixels, width>height,
			   &Kill, &Clear, &Flush, &ChangeMap,
			   &ChangePalette, &DrawLines, &DrawMarkers,
			   &DrawText, &DrawFill, &DrawCells,
			   &DrawDisjoint);
  if (!xEngine) {
    strcpy(gistError, "memory manager failed in GxEngine");
    return 0;
  }

  /* XEngines can repair damage */
  xEngine->e.ClearArea= &ClearArea;

  /* Window in which we draw is the full size of a sheet of paper */
  cwa.background_pixel= bg;
  cwa.border_pixel= xscr->stdColors[1].pixel; /* foreground */
#ifndef EXPOSE_DEBUG
  cwa.backing_store= WhenMapped;
#else
  cwa.backing_store= NotUseful;
#endif
  graphics= XCreateWindow(display, top, x+leftMargin-4, y+topMargin-4,
			  width, height, 4,
			  CopyFromParent, InputOutput, CopyFromParent,
			  CWBackPixel | CWBorderPixel | CWBackingStore, &cwa);

  /* Fill in Engine properties specific to XEngine */
  xEngine->xscr= xscr;
  xEngine->xdpy= xscr->owner;
  xEngine->top= top;
  xEngine->graphics= graphics;
  xEngine->width= width;
  xEngine->height= height;
  xEngine->topMargin= topMargin;
  xEngine->leftMargin= leftMargin;
  xEngine->x= -x;
  xEngine->y= -y;
  xEngine->mapped= 0;
  xEngine->dpi= dpi;

  gxXGCValues.foreground= xscr->stdColors[1].pixel;
  gxXGCValues.background= xscr->stdColors[0].pixel;
  gxXGCValues.line_width= (int)(dpi*10.64385*DEFAULT_LINE_WIDTH); /* 0 */
  gxXGCValues.line_style= LineSolid;
  gxXGCValues.cap_style= CapRound;
  gxXGCValues.join_style= JoinRound;
  gxXGCValues.fill_style= FillSolid;
  gxXGCValues.fill_rule= WindingRule;
  xEngine->gc= XCreateGC(display, graphics, GCForeground | GCBackground |
			 GCLineWidth | GCLineStyle | GCCapStyle |
			 GCJoinStyle | GCFillStyle | GCFillRule,
			 &gxXGCValues);
  xEngine->lastOp.color= FG_COLOR;
  xEngine->lastOp.lwidth= 1.0;
  xEngine->lastOp.ltype= L_SOLID;
  xEngine->lastOp.ljoin= 0;
  xEngine->lastOp.tfont= -1;  /* be sure not to match any font... */
  xEngine->lastOp.tsize= 0.0; /* ...or font size... */
  xEngine->lastOp.fontID= -1; /* ...or GIST fontID */

  xEngine->nColors= 0;
  xEngine->pixelMap= 0;
  xEngine->private= 0;  /* GxExactColors guarantees this is not legal */
  xEngine->e.colorMode= 0;

  xEngine->drawable= graphics;
  xEngine->aWidth= xEngine->aHeight= 0;
  xEngine->graphicsX= xEngine->graphicsY= 0;
  xEngine->swapped= xEngine->e.transform;

  xEngine->HandleExpose= xEngine->HandleResize= xEngine->HandleOther= 0;

  /* Note-- If top is completely occluded by its children, it will
     never actually receive any Expose events.  Therefore, ExposeMask
     should be selected for the graphics window as well to assure that
     an Expose event will actually arrive.  */
  XSelectInput(display, top, ExposureMask | StructureNotifyMask);
  XSelectInput(display, graphics, ExposureMask);

  /* Set context to associate XEngine* to top and graphics */
  if (!gxEngineSet) { gxEngine= XUniqueContext(); gxEngineSet= 1; }
  XSaveContext(display, top, gxEngine, (caddr_t)xEngine);
  XSaveContext(display, graphics, gxEngine, (caddr_t)xEngine);

  /* Copying the dispatcher onto the graphics window saves a round trip
     to the server to do an XQueryTree.  */
  CopyXDispatcher(display, top, graphics);

  return xEngine;
}

GxScreen *GxBasic(char *name, char *displayName, int width, int height,
		  Window *top)
{
  GxScreen *xscr= GxConnect(displayName);
  Display *display= xscr? xscr->display : 0;
  unsigned long bg= xscr? xscr->stdColors[0].pixel : 0;
  XSetWindowAttributes cwa;

  if (!xscr) return 0;

  /* Top level window vignettes actual drawing window */
  cwa.background_pixel= bg;
  cwa.border_pixel= bg;
#ifndef EXPOSE_DEBUG
  cwa.backing_store= WhenMapped;
#else
  cwa.backing_store= NotUseful;
#endif
  *top= XCreateWindow(display, xscr->root, 0, 0, width, height,
		      2, CopyFromParent, InputOutput, CopyFromParent,
		      CWBackPixel | CWBorderPixel | CWBackingStore, &cwa);

  /* Connect the basic dispatcher to the top level window */
  if (AddXDispatcher(display, *top, &GxBasicXHandler)) {
    XDestroyWindow(display, *top);
    GxDisconnect(xscr);
    return 0;
  }

  /* Set various properties for window manager */
  GxSetProperties(name, display, *top, width, height);

  return xscr;
}

/* default top window represents 6 inch square */
int gx75width= 450;
int gx100width= 600;
int gx75height= 450;
int gx100height= 600;

Engine *GpBXEngine(char *name, int landscape, int dpi, char *displayName)
{
  int topWidth= DefaultTopWidth(dpi);
  int topHeight= DefaultTopHeight(dpi);
  Window top;
  GxScreen *xscr= GxBasic(name, displayName, topWidth, topHeight, &top);
  Display *display= xscr? xscr->display : 0;
  GpTransform toPixels;
  int x, y;
  XEngine *xEngine;
  int width;
  int height;

  if (!xscr) return 0;

  SetXTransform(&toPixels, landscape, dpi);
  width= (int)toPixels.window.xmax;
  height= (int)toPixels.window.ymin;
  x= (width-topWidth)/2;
  if (landscape) y= (height-topHeight)/2;
  else y= (width-topHeight)/2;
  if (y<0) y= 0;
  if (x<0) x= 0;
  xEngine= GxEngine(name, &toPixels, xscr, top, -x,-y,0,0, sizeof(XEngine));

  /* Only maps children of top, no grandchildren */
  XMapSubwindows(display, top);

  /* Map top level window, then wait for resulting Expose event(?).  */
  XMapWindow(display, top);
  XSync(display, False);

  return (Engine *)xEngine;
}

int GxInput(Engine *engine, GxHandler HandleExpose,
		   GxHandler HandleResize, GxHandler HandleOther,
		   long eventMask)
{
  XEngine *xEngine= GisXEngine(engine);
  if (!xEngine) return 1;
  xEngine->HandleExpose= HandleExpose;
  xEngine->HandleResize= HandleResize;
  xEngine->HandleOther= HandleOther;
  if (HandleOther) XSelectInput(xEngine->xscr->display,
				xEngine->graphics, ExposureMask | eventMask);
  return 0;
}

XEngine *GisXEngine(Engine *engine)
{
  return (engine && engine->type==xType)? (XEngine *)engine : 0;
}

/* ------------------------------------------------------------------------ */

static void ClearPixmap(XEngine *xeng)
{
  /* Incredibly, (1) XClearWindow doesn't work for pixmaps, and
     (2) there is no way to inquire the current fg and bg for a
     GC which is not R4 specific.  Sigh.  */
  GxScreen *xscr= xeng->xscr;
  GC gc= xeng->gc;
  Display *display= xscr? xscr->display : 0;
  if (!xscr) return;

  /* Set foreground to background, fill pixmap, restore last color */
  XSetForeground(display, gc, xscr->stdColors[0].pixel);
  /* Note that XFillRectangle does not fill the last row and column of
     the specified rectangle... */
  XFillRectangle(display, xeng->drawable, gc,
		 0, 0, xeng->aWidth+1, xeng->aHeight+1);
  GxSetColor(xeng, xeng->lastOp.color);
}

int GxAnimate(Engine *engine, GpBox *viewport)
{
  XEngine *xeng= GisXEngine(engine);
  int x, y;
  GpBox *v, *w;
  GpReal xmin, xmax, ymin, ymax;
  GpReal scalx, offx, scaly, offy;
  XRectangle xrect;

  if (!xeng || !xeng->xscr) return 1;
  if (xeng->drawable!=xeng->graphics) GxDirect(engine);

  v= &xeng->e.transform.viewport;  /* NDC */
  w= &xeng->e.transform.window;    /* pixels */

  /* get NDC-to-pixel mapping coefficients */
  scalx= xeng->e.devMap.x.scale;
  offx= xeng->e.devMap.x.offset;
  scaly= xeng->e.devMap.y.scale;
  offy= xeng->e.devMap.y.offset;

  /* Clip given viewport to portion of NDC space which is actually
     visible now -- note that v is either gLandscape or gPortrait,
     so that min<max for v; must also be true for input viewport */
  GetVisibleNDC(xeng, &xmin, &xmax, &ymin, &ymax);
  if (viewport->xmin>xmin) xmin= viewport->xmin;
  if (viewport->xmax<xmax) xmax= viewport->xmax;
  if (viewport->ymin>ymin) ymin= viewport->ymin;
  if (viewport->ymax<ymax) ymax= viewport->ymax;

  /* Install NDC->pixel transform for animation pixmap */
  v->xmin= xmin;
  v->xmax= xmax;
  v->ymin= ymin;
  v->ymax= ymax;

  /* Set the engine transform to map the specified viewport into
     Pixmap pixels, and get (x,y) offset from full graphics window
     pixels to Pixmap pixels */
  w->xmin= scalx*xmin+offx;
  w->xmax= scalx*xmax+offx;
  if (w->xmax > w->xmin) {
    x= (int)w->xmin;
    w->xmax-= w->xmin;
    w->xmin= 0.0;
  } else {
    x= (int)w->xmax;
    w->xmin-= w->xmax;
    w->xmax= 0.0;
  }
  w->ymin= scaly*ymin+offy;
  w->ymax= scaly*ymax+offy;
  if (w->ymax > w->ymin) {
    y= (int)w->ymin;
    w->ymax-= w->ymin;
    w->ymin= 0.0;
  } else {
    y= (int)w->ymax;
    w->ymin-= w->ymax;
    w->ymax= 0.0;
  }
  GpDeviceMap((Engine *)xeng);
  GetXRectangle(&xeng->e.devMap, v, &xrect);

  if (xeng->drawable!=xeng->graphics) {
    XFreePixmap(xeng->xscr->display, xeng->drawable);
  } else {
    xeng->gca=
      XCreateGC(xeng->xscr->display, xeng->graphics, 0, &gxXGCValues);
  }

  xeng->drawable= XCreatePixmap(xeng->xscr->display, xeng->graphics,
				xrect.width, xrect.height,
				xeng->xscr->v->depth);
  xeng->aWidth= xrect.width;
  xeng->aHeight= xrect.height;
  xeng->graphicsX= x;
  xeng->graphicsY= y;

  /* Set coordinate mapping and clipping for pixmap--
     clipping includes both clipping for drawing to pixmap (gc),
     and clipping for copying pixmap to graphics window (gca).  */
  ChangeMap((Engine *)xeng);
  xrect.x= x;
  xrect.y= y;
  XSetClipRectangles(xeng->xscr->display, xeng->gca, 0, 0,
		     &xrect, 1, YXBanded);

  ClearPixmap(xeng);
  return 0;
}

static void GetVisibleNDC(XEngine *xeng,
			  GpReal *xn, GpReal *xx, GpReal *yn, GpReal *yx)
{
  GpReal scalx= xeng->e.devMap.x.scale;
  GpReal offx= xeng->e.devMap.x.offset;
  GpReal scaly= xeng->e.devMap.y.scale;
  GpReal offy= xeng->e.devMap.y.offset;
  int xmin, xmax, ymin, ymax;

  Window root;
  int x, y;
  unsigned int width, height, border, depth;
  if (!xeng->xscr ||
      !XGetGeometry(xeng->xscr->display, xeng->top, &root, &x, &y,
		    &width, &height, &border, &depth)) {
    /* Use default window size on failure.  */
    width= DefaultTopWidth(xeng->dpi);
    height= DefaultTopHeight(xeng->dpi);
  }

  xmin= xeng->x;
  xmax= xmin + (int)width - xeng->leftMargin;
  ymax= xeng->y;
  ymin= ymax + (int)height - xeng->topMargin;

  /* Convert pixels to NDC coordinates */
  *xn= (xmin-offx)/scalx;
  *xx= (xmax-offx)/scalx;
  *yn= (ymin-offy)/scaly;
  *yx= (ymax-offy)/scaly;
}

int GxStrobe(Engine *engine, int clear)
{
  XEngine *xeng= GisXEngine(engine);

  if (!xeng || xeng->drawable==xeng->graphics || !xeng->xscr) return 1;

  XCopyArea(xeng->xscr->display, xeng->drawable, xeng->graphics, xeng->gca,
	    0, 0, xeng->aWidth, xeng->aHeight,
	    xeng->graphicsX, xeng->graphicsY);

  if (clear) ClearPixmap(xeng);
  return 0;
}

int GxDirect(Engine *engine)
{
  XEngine *xeng= GisXEngine(engine);

  if (!xeng || xeng->drawable==xeng->graphics || !xeng->xscr) return 1;

  XFreePixmap(xeng->xscr->display, xeng->drawable);
  XFreeGC(xeng->xscr->display, xeng->gca);
  xeng->drawable= xeng->graphics;

  /* Set coordinate map and clipping in GC back to value appropriate
     for graphics window.  */
  xeng->e.transform= xeng->swapped;
  GpDeviceMap((Engine *)xeng);
  ChangeMap((Engine *)xeng);

  return 0;
}

/* ------------------------------------------------------------------------ */

static void (*XErrHandler)(char *errMsg)= 0;
static int gxErrorCode= 0;
static Display *gxErrorDisplay= 0;

/* this routine actually calls the XErrHandler, which may not return
   and/or which may trigger additional X protocol requests */
static void GxErrorHandler(void)
{
  char msg[80];
  gxErrorFlag= 0;
  XGetErrorText(gxErrorDisplay, gxErrorCode, msg, 80);
  XErrHandler(msg);
}

/* this routine is called from Xlib and must return without triggering
   any X protocol requests or attempting to read X input events, hence
   all it can do is to set a flag which will trigger the actual high
   level error handling routine at GpFlush time (fma or wait for input)
   or event processing (GxBasicXHandler) -- hopefully, anything which
   triggered the error will be followed by a fairly prompt call to
   the GxErrorHandler that alerts the caller that something is wrong */
static int YXError(Display *dpy, XErrorEvent *xev)
{
  if (shutDownDisplay!=dpy) {
    if (!gxErrorFlag) {
      gxErrorDisplay= dpy;
      gxErrorCode= xev->error_code;
    }
    gxErrorFlag++;  /* may as well keep count */
  } else {
    gxErrorFlag= 0;
  }
  return 1;
}

static int YXIOError(Display *dpy)
{
  /* panic time -- need to kill every X engine on this display */
  Engine *eng= 0;
  XEngine *xeng;
  GxDisplay *xdpy= 0;
  RemoveXDispatcher(dpy);
  do {
    for (eng=GpNextEngine(eng) ; eng ; eng=GpNextEngine(eng)) {
      xeng= GisXEngine(eng);
      if (xeng && xeng->xdpy->display==dpy) break;
    }
    if (eng) {
      if (HLevelHook) HLevelHook((Engine *)xeng);
      GmFree(xeng->pixelMap);
      xeng->pixelMap= 0;
      xeng->xscr= 0;
      xdpy= xeng->xdpy;
      xeng->xdpy= 0;
      GpDelEngine(eng);
    }
  } while (eng);
  if (xdpy) {
    extern void GxUnlink(GxDisplay *owner);
    GmFree(xdpy->normalizedName);
    GmFree(xdpy->screens);
    XFree((char *)xdpy->v);
    GxUnlink(xdpy);
    /* do not attempt to call XCloseDisplay */
    GmFree(xdpy);
  }
  XErrHandler("Xlib I/O error (X engines killed)");
  return 1;
}

int GpSetXHandler(void (*ErrHandler)(char *errMsg))
{
  /* install X error handlers which don't call exit */
  XErrHandler= ErrHandler;
  XSetErrorHandler(&YXError);
  XSetIOErrorHandler(&YXIOError);
  return 0;
}

int GxWaitForExpose(Engine *engine)
{
  XEngine *xeng= GisXEngine(engine);
  if (!xeng || !xeng->xscr) return 1;
  if (!xeng->mapped) {
    Display *dpy= xeng->xscr->display;
    XEvent event;
    /* NB--
       xeng->top is completely occluded by its children and never
       receives an expose event.  */
    XWindowEvent(dpy, xeng->graphics, ExposureMask, &event);
    xeng->mapped= 1;
    if (xeng->HandleExpose)
      /* the alternate handler should probably call GxExpose */
      xeng->HandleExpose(&xeng->e, xeng->e.drawing, &event);
    else
      GxExpose(&xeng->e, xeng->e.drawing, &event);
  }
  return 0;
}

/* ------------------------------------------------------------------------ */
