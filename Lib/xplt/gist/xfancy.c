/*
 * XFANCY.C
 *
 * $Id$
 *
 * Implement the basic X windows engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "xfancy.h"
#include "xfont.h"
#include "draw.h"
#include "dispat.h"

#include <X11/cursorfont.h>

#ifdef STDC_HEADERS
#include <string.h>
#else
#ifndef SIZE_T_TYPE
#define SIZE_T_TYPE unsigned long
#endif
extern SIZE_T_TYPE strlen(const char *);
extern char *strcpy(char *, const char *);
#endif

extern int sprintf(char *s, const char *format, ...);

extern double log10(double);

#ifdef HAS_EXP10
extern double exp10(double);
#else
extern double pow(double, double);
#define exp10(x) pow(10.,x)
#endif

/* ------------------------------------------------------------------------ */

static void FXKill(Engine *engine);

static void HandleExpose(Engine *engine, Drawing *drawing, XEvent *event);
static void HandleResize(Engine *engine, Drawing *drawing, XEvent *event);
static void HandleOther(Engine *engine, Drawing *drawing, XEvent *event);

static void RedrawMessage(FXEngine *fxEngine);
static void MovePointer(FXEngine *fxEngine, Drawing *drawing,
			XMotionEvent *event);
static void PressZoom(FXEngine *fxEngine, Drawing *drawing,
		      XButtonEvent *event);
static void ReleaseZoom(FXEngine *fxEngine, Drawing *drawing,
			XButtonEvent *event);

static void RedrawButton(FXEngine *fxEngine);
static void EnterButton(FXEngine *fxEngine, XCrossingEvent *event);
static void LeaveButton(FXEngine *fxEngine);
static void PressButton(FXEngine *fxEngine, XButtonEvent *event);
static void ReleaseButton(FXEngine *fxEngine, Drawing *drawing,
			  XButtonEvent *event);

static void ButtonAction(FXEngine *fxEngine, Drawing *drawing);
static void HighlightButton(FXEngine *fxEngine);
static void UnHighlightButton(FXEngine *fxEngine);
static void InvertButton(FXEngine *fxEngine);

static void ResetZoom(FXEngine *fxEngine);
static void DoZoom(GpReal factor, GpReal w0, GpReal w1,
		   GpReal *wmin, GpReal *wmax);
static void AltZoom(GpReal w0, GpReal w1, GpReal *wmin, GpReal *wmax);
static int FindSystem(FXEngine *fxEngine, Drawing *drawing, int x, int y,
		      GeSystem **system, GpReal *xr, GpReal *yr);
static void Find1System(FXEngine *fxEngine, Drawing *drawing, int iSystem,
			int x, int y, GeSystem **system,
			GpReal *xr, GpReal *yr);
static GeSystem *GetSystemN(Drawing *drawing, int n);
static int FindAxis(GeSystem *system, GpReal x, GpReal y);
static void FindCoordinates(GeSystem *system, GpReal xNDC, GpReal yNDC,
			    GpReal *xWC, GpReal *yWC);
static GpReal GetFormat(char *format, GpReal w,
			GpReal wmin, GpReal wmax, int isLog);
static void DrawRubber(FXEngine *fxEngine, int x, int y);
static void SetRubberPixel(FXEngine *fxEngine, int which);

/* ------------------------------------------------------------------------ */

static void (*XKill)(Engine *engine)= 0;

static void FXKill(Engine *engine)
{
  FXEngine *fxEngine= (FXEngine *)engine;

  if (fxEngine->xe.top!=None) {
    /* On first pass, call the default XKill routine (xbasic.c).
       After destroying all of the basic X resources, this will set
       the top window to None, and recursively call this routine to
       delete the additional resources used by the fancy engine.  */
    XKill(engine);

  } else {
    /* This is reached on the recursive pass of an ordinary Kill, but
       can be reached directly on a window-manager initiated window
       kill -- this baroque coding would require additional kill hooks
       in the XEngine data structure to be avoided...  Someday.  */
    GxScreen *xscr= fxEngine->xe.xscr;
    if (!xscr) return;
    XFreeGC(xscr->display, fxEngine->gc);
    XFreeCursor(xscr->display, fxEngine->cursor);
  }
}

Engine *GpFXEngine(char *name, int landscape, int dpi, char *displayName)
{
  int topWidth= DefaultTopWidth(dpi);   /* not including button, message */
  int topHeight= DefaultTopHeight(dpi);
  int heightButton, widthButton, baseline;
  Window top, button, message;
  GxScreen *tmpscr= GxConnect(displayName);  /* GxBasic needs window height */
  Display *display= tmpscr? tmpscr->display : 0;
  GxDisplay *xdpy= tmpscr? tmpscr->owner : 0;
  GxScreen *xscr;
  GpTransform toPixels;
  int x, y;
  FXEngine *fxEngine;
  int width, height, direction, ascent, descent;
  XFontStruct *font;
  XCharStruct overall;
  XSetWindowAttributes cwa;

  if (!tmpscr) return 0;
  font= xdpy->permFont;
  if (!font) font= xdpy->defaultFont;
  if (!font) {  /* forget it, no fonts... */
    Engine *engine= GpBXEngine(name, landscape, dpi, displayName);
    GxDisconnect(tmpscr);
    return engine;
  }
  baseline= font->ascent+2;
  heightButton= baseline+font->descent+4;  /* leave at least 2 lines
					      above and below text */
  xscr= GxBasic(name, displayName, topWidth, topHeight+heightButton+2, &top);
  GxDisconnect(tmpscr);   /* now that xscr set, don't need this */
  XTextExtents(font, "System", 6, &direction, &descent, &ascent, &overall);
  widthButton= overall.width+8;  /* leave at least 4 pixels before and after */

  /* set toPixels as by SetXTransform(&toPixels, landscape, dpi) */
  toPixels.viewport= landscape? gLandscape : gPortrait;
  toPixels.window.xmin= 0.0;
  toPixels.window.xmax= PixelsPerNDC(dpi)*toPixels.viewport.xmax;
  toPixels.window.ymin= PixelsPerNDC(dpi)*toPixels.viewport.ymax;
  toPixels.window.ymax= 0.0;
  width= (int)toPixels.window.xmax;
  height= (int)toPixels.window.ymin;
  x= (width-topWidth)/2;
  if (landscape) y= (height-topHeight)/2;
  else y= (width-topHeight)/2;
  if (y<0) y= 0;
  if (x<0) x= 0;
  fxEngine= (FXEngine *)GxEngine(name, &toPixels, xscr, top,
				 -x,-y,heightButton+2,0,
				 sizeof(FXEngine));
  XKill= fxEngine->xe.e.Kill;
  fxEngine->xe.e.Kill= &FXKill;

  cwa.background_pixel= xscr->stdColors[0].pixel;  /* background */
  cwa.border_pixel= xscr->stdColors[1].pixel;      /* foreground */
  /* cwa.backing_store= WhenMapped; */
  cwa.backing_store= NotUseful;
  button= XCreateWindow(display, top, 0, 0,
			widthButton, heightButton, 1,
			CopyFromParent, InputOutput, CopyFromParent,
			CWBackPixel | CWBorderPixel | CWBackingStore, &cwa);

  message= XCreateWindow(display, top, widthButton, 0,
			 width-widthButton, heightButton, 1,
			 CopyFromParent, InputOutput, CopyFromParent,
			 CWBackPixel | CWBorderPixel | CWBackingStore, &cwa);

  /* set the extra information in the fxEngine structure */
  fxEngine->button= button;
  fxEngine->message= message;
  fxEngine->baseline= baseline;
  fxEngine->heightButton= heightButton;
  fxEngine->widthButton= widthButton;
  gxXGCValues.foreground= xscr->stdColors[1].pixel;
  gxXGCValues.background= xscr->stdColors[0].pixel;
  gxXGCValues.font= font->fid;
  gxXGCValues.line_width= 3;
  fxEngine->gc= XCreateGC(display, message, GCForeground | GCBackground |
			  GCFont | GCLineWidth, &gxXGCValues);
  fxEngine->font= font;
  fxEngine->cursor= XCreateFontCursor(display, XC_crosshair);
  fxEngine->buttonState= 0;
  fxEngine->iSystem= -1;
  strcpy(fxEngine->msgText, "Press 1, 2, 3 to zoom in, pan, zoom out");
  XTextExtents(font, "Press 1, 2, 3 to zoom in, pan, zoom out", 39,
	       &direction, &descent, &ascent, &overall);
  fxEngine->msgWidth= overall.width;
  fxEngine->zoomState= fxEngine->zoomSystem= fxEngine->zoomAxis= 0;
  fxEngine->zoomX= fxEngine->zoomY= 0.0;

  /* change cursor to crosshair in graphics window */
  XDefineCursor(display, fxEngine->xe.graphics, fxEngine->cursor);

  /* Add button and motion events to graphics window,
     and set the fancy event handlers.  */
  GxInput((Engine *)fxEngine, &HandleExpose, &HandleResize, &HandleOther,
	  ButtonPressMask | ButtonReleaseMask | PointerMotionMask);

  XSelectInput(display, button, ExposureMask | EnterWindowMask |
	       LeaveWindowMask | ButtonPressMask | ButtonReleaseMask);

  XSelectInput(display, message, ExposureMask);

  /* Only maps children of top, no grandchildren */
  XMapSubwindows(display, top);

  /* Map top level window, then wait for resulting Expose event(?).  */
  XMapWindow(display, top);
  XSync(display, False);

  return (Engine *)fxEngine;
}

/* ------------------------------------------------------------------------ */

static void HandleExpose(Engine *engine, Drawing *drawing, XEvent *event)
{
  /* Already know that event->type is Expose */
  FXEngine *fxEngine= (FXEngine *)engine;
  Window w= event->xexpose.window;
  if (w==fxEngine->xe.graphics) {
    GxExpose(engine, drawing, event);
  } else if (w==fxEngine->button) {
    RedrawButton(fxEngine);
  } else if (w==fxEngine->message) {
    RedrawMessage(fxEngine);
  }
}

static void HandleResize(Engine *engine, Drawing *drawing, XEvent *event)
{
  /* Already know that event->type is ConfigureNotify */
  XEngine *xEngine= (XEngine *)engine;
  if (event->xconfigure.window==xEngine->top)
    GxRecenter(xEngine, event->xconfigure.width, event->xconfigure.height);
}

static void HandleOther(Engine *engine, Drawing *drawing, XEvent *event)
{
  FXEngine *fxEngine= (FXEngine *)engine;
  Window w= event->xany.window;
  if (w==fxEngine->xe.graphics) {
    if (event->type==MotionNotify)
      MovePointer(fxEngine, drawing, &event->xmotion);
    else if (event->type==ButtonPress)
      PressZoom(fxEngine, drawing, &event->xbutton);
    else if (event->type==ButtonRelease)
      ReleaseZoom(fxEngine, drawing, &event->xbutton);

  } else if (w==fxEngine->button) {
    if (event->type==EnterNotify)
      EnterButton(fxEngine, &event->xcrossing);
    else if (event->type==LeaveNotify)
      LeaveButton(fxEngine);
    else if (event->type==ButtonPress)
      PressButton(fxEngine, &event->xbutton);
    else if (event->type==ButtonRelease)
      ReleaseButton(fxEngine, drawing, &event->xbutton);
  }
}

/* ------------------------------------------------------------------------ */

static void RedrawButton(FXEngine *fxEngine)
{
  if (fxEngine->xe.xscr) {
    Display *dpy= fxEngine->xe.xscr->display;
    Window button= fxEngine->button;
    XClearWindow(dpy, button);
    XDrawString(dpy, button,
		fxEngine->gc, 3, fxEngine->baseline,
		"System", 6);
    if (fxEngine->buttonState) {
      HighlightButton(fxEngine);
      if (fxEngine->buttonState==2) InvertButton(fxEngine);
    }
  }
}

static void HighlightButton(FXEngine *fxEngine)
{
  if (fxEngine->xe.xscr) {
    Display *dpy= fxEngine->xe.xscr->display;
    Window button= fxEngine->button;
    XDrawRectangle(dpy, button, fxEngine->gc, 1, 1,
		   fxEngine->widthButton-4, fxEngine->heightButton-3);
  }
}

static void UnHighlightButton(FXEngine *fxEngine)
{
  GxScreen *xscr= fxEngine->xe.xscr;
  if (xscr) {
    unsigned long fg= xscr->stdColors[1].pixel;
    Display *dpy= xscr->display;
    Window button= fxEngine->button;
    GC gc= fxEngine->gc;
    XSetFunction(dpy, gc, GXxor);
    XSetForeground(dpy, gc, fg^xscr->stdColors[0].pixel);
    XDrawRectangle(dpy, button, gc, 1, 1,
		   fxEngine->widthButton-4, fxEngine->heightButton-3);
    XSetForeground(dpy, gc, fg);
    XSetFunction(dpy, gc, GXcopy);
  }
}

static void InvertButton(FXEngine *fxEngine)
{
  GxScreen *xscr= fxEngine->xe.xscr;
  if (xscr) {
    unsigned long fg= xscr->stdColors[1].pixel;
    Display *dpy= xscr->display;
    Window button= fxEngine->button;
    GC gc= fxEngine->gc;
    XSetFunction(dpy, gc, GXxor);
    XSetForeground(dpy, gc, fg^xscr->stdColors[0].pixel);
    XFillRectangle(dpy, button, gc, 0, 0,
		   fxEngine->widthButton, fxEngine->heightButton);
    XSetForeground(dpy, gc, fg);
    XSetFunction(dpy, gc, GXcopy);
  }
}

static void EnterButton(FXEngine *fxEngine, XCrossingEvent *event)
{
  if ((event->state&(Button1Mask|Button2Mask|Button3Mask))==0 &&
      fxEngine->buttonState==0) {
    fxEngine->buttonState= 1;
    HighlightButton(fxEngine);
  } else if (fxEngine->buttonState!=0) {
    fxEngine->buttonState= 0;
    RedrawButton(fxEngine);
  }
}

static void LeaveButton(FXEngine *fxEngine)
{
  int state= fxEngine->buttonState;
  fxEngine->buttonState= 0;
  if (state==1) UnHighlightButton(fxEngine);
  else if (state) RedrawButton(fxEngine);
}

static void PressButton(FXEngine *fxEngine, XButtonEvent *event)
{
  if ((event->state&(Button1Mask|Button2Mask|Button3Mask))==0 &&
      fxEngine->buttonState==1) {
    InvertButton(fxEngine);
    fxEngine->buttonState= 2;
  } else {
    LeaveButton(fxEngine);
  }
}

static void ReleaseButton(FXEngine *fxEngine, Drawing *drawing,
			  XButtonEvent *event)
{
  if (fxEngine->buttonState==2) {
    ButtonAction(fxEngine, drawing);
    InvertButton(fxEngine);
    fxEngine->buttonState= 1;
  } else {
    LeaveButton(fxEngine);
  }
}

/* ------------------------------------------------------------------------ */

static void ButtonAction(FXEngine *fxEngine, Drawing *drawing)
{
  int nSystems= drawing? drawing->nSystems : 0;
  int iSystem= fxEngine->iSystem;
  iSystem++;
  if (iSystem<=nSystems) fxEngine->iSystem= iSystem;
  else fxEngine->iSystem= iSystem= -1;
  sprintf(fxEngine->msgText, "%s%2d", iSystem>=0?"=":":",
	  iSystem>=0 ? iSystem : 0);
  RedrawMessage(fxEngine);
}

/* ------------------------------------------------------------------------ */

static void RedrawMessage(FXEngine *fxEngine)
{
  if (fxEngine->xe.xscr) {
    Display *dpy= fxEngine->xe.xscr->display;
    Window message= fxEngine->message;
    char *msgText= fxEngine->msgText;
    int len= (int)strlen(fxEngine->msgText);
    int direction, ascent, descent;
    XCharStruct overall;

    XTextExtents(fxEngine->font, msgText, len,
		 &direction, &descent, &ascent, &overall);
    if (fxEngine->msgWidth!=overall.width) {
      XClearWindow(dpy, message);
      fxEngine->msgWidth= overall.width;
    }

    XDrawImageString(dpy, message, fxEngine->gc,
		     4, fxEngine->baseline, msgText, len);
  }
}

static char stdFormat[]= "%7.4f";
static int rubberBanding= 0;
static unsigned long rubberPixel;
static int anchorX, anchorY, oldX, oldY;

static void MovePointer(FXEngine *fxEngine, Drawing *drawing,
			XMotionEvent *event)
{
  int iSystem= fxEngine->iSystem;
  int locked, logX= 0, logY= 0;
  char format[24];  /* e.g.- "%s%2d (%11.3e, %11.3e)" */
  char xFormat[16], yFormat[16], *f1, *f2;
  GpReal xWC, yWC;
  GeSystem *system;
  XEvent later;

  /* Allow for slow response by cleaning off any pending motion
     events -- no sense in displaying any but the latest.  */
  while (XCheckTypedEvent(event->display, MotionNotify, &later))
    event= &later.xmotion;

  if (!drawing || fxEngine->buttonState) return;

  /* find the system number and world coordinates */
  if (iSystem>=0) {
    /* coordinate system is locked */
    Find1System(fxEngine, drawing, iSystem, event->x, event->y,
		&system, &xWC, &yWC);
    if (!system) iSystem= 0;
    locked= 1;
  } else {
    /* select coordinate system under pointer */
    iSystem= FindSystem(fxEngine, drawing, event->x, event->y,
			&system, &xWC, &yWC);
    locked= 0;
  }

  if (system) {
    logX= system->flags&D_LOGX;
    logY= system->flags&D_LOGY;
    xWC= GetFormat(xFormat, xWC, system->trans.window.xmin,
		   system->trans.window.xmax, logX);
    f1= xFormat;
    yWC= GetFormat(yFormat, yWC, system->trans.window.ymin,
		   system->trans.window.ymax, logY);
    f2= yFormat;
  } else {
    f1= f2= stdFormat;
  }
  sprintf(format, "%%s%%2d (%s, %s)", f1, f2);
  sprintf(fxEngine->msgText, format,
	  locked?"=":":", iSystem, xWC, yWC);

  RedrawMessage(fxEngine);
  if (rubberBanding) DrawRubber(fxEngine, event->x, event->y);
}

GpReal gxZoomFactor= 1.5;

static int (*PtClCallBack)(Engine *engine, int system,
			   int release, GpReal x, GpReal y,
			   int butmod, GpReal xn, GpReal yn)= 0;
static int ptClStyle= 0, ptClSystem= 0, ptClCount= 0;

static unsigned int cShapes[3]= {
  XC_sb_h_double_arrow, XC_sb_v_double_arrow, XC_fleur };

static void PressZoom(FXEngine *fxEngine, Drawing *drawing,
		      XButtonEvent *event)
{
  int but= event->button;
  if (!drawing || !fxEngine->xe.xscr) return;
  if (fxEngine->zoomState==0) {
    /* record button number as zoomState */
    if (but==Button1) fxEngine->zoomState= 1;
    else if (but==Button2) fxEngine->zoomState= 2;
    else if (but==Button3) fxEngine->zoomState= 3;
    if (fxEngine->zoomState) {
      if (PtClCallBack)
	fxEngine->zoomState|= 8;
      else if (event->state&(ControlMask|ShiftMask|Mod1Mask))
	fxEngine->zoomState|= 4;
    }
    if (fxEngine->zoomState) {
      /* record system and axis, x and y, change cursor */
      Display *dpy= fxEngine->xe.xscr->display;
      GeSystem *system;
      int iSystem, axis;
      if (!PtClCallBack || ptClSystem<0) {
	iSystem= FindSystem(fxEngine, drawing, event->x, event->y,
			    &system, &fxEngine->zoomX, &fxEngine->zoomY);
	axis= FindAxis(system, fxEngine->zoomX, fxEngine->zoomY);
      } else {
	iSystem= ptClSystem;
	Find1System(fxEngine, drawing, iSystem, event->x, event->y,
		    &system, &fxEngine->zoomX, &fxEngine->zoomY);
	if (!system) iSystem= ptClSystem= 0;
	axis= 3;
      }
      fxEngine->zoomSystem= iSystem;

      if (PtClCallBack) {
	GpXYMap *map= &fxEngine->xe.e.map;   /* NDC->VDC (x,y) mapping */
	GpReal xNDC= ((GpReal)event->x - map->x.offset)/map->x.scale;
	GpReal yNDC= ((GpReal)event->y - map->y.offset)/map->y.scale;
	int button= fxEngine->zoomState&3;
	ptClCount--;
	if (PtClCallBack(&fxEngine->xe.e, iSystem, 0,
			 fxEngine->zoomX, fxEngine->zoomY,
			 button, xNDC, yNDC)) {
	  /* callback routine signals abort */
	  ResetZoom(fxEngine);
	  return;
	}
      }

      if ((fxEngine->zoomAxis= axis)) {
	if (fxEngine->zoomState<4) {
	  fxEngine->cursor= XCreateFontCursor(dpy, cShapes[axis-1]);
	} else if (!PtClCallBack || ptClStyle) {
	  fxEngine->zoomAxis= 3;
	  fxEngine->cursor= XCreateFontCursor(dpy, XC_left_ptr);
	  SetRubberPixel(fxEngine, (fxEngine->zoomState&3)-1);
	  anchorX= oldX= event->x;
	  anchorY= oldY= event->y;
	  rubberBanding= PtClCallBack? ptClStyle : 1;
	  DrawRubber(fxEngine, anchorX, anchorY);
	}
	XDefineCursor(dpy, fxEngine->xe.graphics, fxEngine->cursor);
      } else if (!PtClCallBack) {
	/* no such thing as zoom or pan outside of all systems */
	fxEngine->zoomState= 0;
      }
    }

  } else {
    /* abort if 2nd button pressed */
    ResetZoom(fxEngine);
  }
}

static void ReleaseZoom(FXEngine *fxEngine, Drawing *drawing,
			XButtonEvent *event)
{
  int zoomState= fxEngine->zoomState;
  if (zoomState) {
    /* perform the indicated zoom operation */
    int iSystem= fxEngine->zoomSystem;
    GeSystem *system= GetSystemN(drawing, iSystem);
    GpXYMap *map= &fxEngine->xe.e.map;   /* NDC->VDC (x,y) mapping */
    GpReal x, xNDC= ((GpReal)event->x - map->x.offset)/map->x.scale;
    GpReal y, yNDC= ((GpReal)event->y - map->y.offset)/map->y.scale;

    /* get current pointer position (in world coordinates) */
    if (system &&
	/* be sure system has been scanned if limits extreme */
	(!(system->rescan || system->unscanned>=0) || !GdScan(system))) {
      FindCoordinates(system, xNDC, yNDC, &x, &y);
    } else {
      x= xNDC;
      y= yNDC;
      iSystem= ptClSystem= 0;
    }

    if (!PtClCallBack) {
      int axis= fxEngine->zoomAxis;
      GpReal factor= 1.0;

      if (zoomState==1) factor= 1.0/gxZoomFactor;
      else if (zoomState==3) factor= gxZoomFactor;

      /* the redraw triggered here can mess up a rubber band box/line */
      if (rubberBanding) {
	DrawRubber(fxEngine, anchorX, anchorY);
	rubberBanding= 0;
      }

      /* switch to current drawing and engine temporarily, then save
	 limits if this is first mouse-driven zoom */
      GdSetDrawing(drawing);
      GpPreempt((Engine *)fxEngine);
      if (!(system->flags&D_ZOOMED)) GdSaveLimits(0);

      if (zoomState<4) {
	if (axis&1) DoZoom(factor, fxEngine->zoomX, x,
			   &system->trans.window.xmin,
			   &system->trans.window.xmax);
	if (axis&2) DoZoom(factor, fxEngine->zoomY, y,
			   &system->trans.window.ymin,
			   &system->trans.window.ymax);
      } else {
	if (axis&1) AltZoom(fxEngine->zoomX, x, &system->trans.window.xmin,
			    &system->trans.window.xmax);
	if (axis&2) AltZoom(fxEngine->zoomY, y, &system->trans.window.ymin,
			    &system->trans.window.ymax);
      }
      system->flags&= ~(D_XMIN | D_XMAX | D_YMIN | D_YMAX);
      system->flags|= D_ZOOMED;
      system->rescan= 1;

      /* The redraw must be done immediately -- nothing else triggers it.  */
      GdDraw(1);
      GpPreempt(0);     /* not correct if zoomed during a preempt... */
      GdSetDrawing(0);

    } else {
      int state= event->state;
      int modifier= 0;
      if (state&ShiftMask) modifier|= 1;
      else if (state&LockMask) modifier|= 2;
      else if (state&ControlMask) modifier|= 4;
      else if (state&Mod1Mask) modifier|= 8;
      else if (state&Mod2Mask) modifier|= 16;
      else if (state&Mod3Mask) modifier|= 32;
      else if (state&Mod4Mask) modifier|= 64;
      else if (state&Mod5Mask) modifier|= 128;
      ptClCount--;
      PtClCallBack(&fxEngine->xe.e, iSystem, 1, x, y, modifier, xNDC, yNDC);
    }

    /* free the zoom/pan cursor and reset zoomState */
    ResetZoom(fxEngine);
  }
}

static void ResetZoom(FXEngine *fxEngine)
{
  /* free the zoom cursor and reset the zoom state */
  if (rubberBanding) {
    DrawRubber(fxEngine, anchorX, anchorY);
    rubberBanding= 0;
  }
  if (fxEngine->zoomState && fxEngine->xe.xscr) {
    Display *dpy= fxEngine->xe.xscr->display;
    XFreeCursor(dpy, fxEngine->cursor);
    fxEngine->cursor= XCreateFontCursor(dpy, XC_crosshair);
    XDefineCursor(dpy, fxEngine->xe.graphics, fxEngine->cursor);
  }
  fxEngine->zoomState= 0;
  PtClCallBack= 0;
}

static void DoZoom(GpReal factor, GpReal w0, GpReal w1,
		   GpReal *wmin, GpReal *wmax)
{
  GpReal wn= *wmin;
  GpReal wx= *wmax;
  /* With origin at the release point, expand the scale, then put this origin
     at the press point.  This is a translation (pan) if factor==1.0.  */
  *wmin= w0-factor*(w1-wn);
  *wmax= w0+factor*(wx-w1);
}

static void AltZoom(GpReal w0, GpReal w1, GpReal *wmin, GpReal *wmax)
{
  GpReal wn= *wmin;
  GpReal wx= *wmax;
  /* Expand marked interval to full scale.  */
  if (wn<wx) {
    if (w0<w1) { *wmin= w0; *wmax= w1; }
    else if (w0>w1) { *wmin= w1; *wmax= w0; }
    else {
      wx= wx-wn;
      if (!wx) { if (wn) wx=0.0001*wn; else wx= 1.e-6; }
      else wx*= 0.01;
      *wmin= w0-wx; *wmax= w1+wx;
    }
  } else {
    if (w0<w1) { *wmin= w1; *wmax= w0; }
    else if (w0>w1) { *wmin= w0; *wmax= w1; }
    else {
      wx= wn-wx;
      if (!wx) { if (wn) wx=0.0001*wn; else wx= 1.e-6; }
      else wx*= 0.01;
      *wmin= w0+wx; *wmax= w1-wx;
    }
  }
}

static int FindSystem(FXEngine *fxEngine, Drawing *drawing, int x, int y,
		      GeSystem **system, GpReal *xr, GpReal *yr)
{
  GeSystem *sys= drawing->systems, *thesys=sys;
  int nSystems= drawing->nSystems;
  GpXYMap *map= &fxEngine->xe.e.map;  /* NDC->VDC (x,y) mapping */
  GpReal xn= ((GpReal)x - map->x.offset)/map->x.scale;
  GpReal yn= ((GpReal)y - map->y.offset)/map->y.scale;
  GpBox *box;
  int i, iSystem=0;
  GpReal min=9., tmp; /* assume all viewports have area<9 */
  for (i=nSystems ; i>0 ; i--) {
    sys= (GeSystem *)sys->el.prev;
    if (!sys->elements ||
	/* be sure system has been scanned if limits extreme */
	((sys->rescan || sys->unscanned>=0) &&
	 GdScan(sys))) continue;
    box= &sys->trans.viewport; 
    if (xn>=box->xmin && xn<=box->xmax && yn>=box->ymin && yn<=box->ymax)
      { 
	tmp= (box->xmax-box->xmin)*(box->ymax-box->ymin);
	if(tmp<0) tmp= -tmp;
	if(tmp<min) { 
	  min= tmp;
	  iSystem= i;
	  thesys= sys;
	}
      }
  }
  if (!iSystem) { /* look for nearest axis */
    min= 9.; /* now assume mouse is closer to an axis than 9 units */
    for (i=nSystems ; i>0 ; i--) {
      sys= (GeSystem *)sys->el.prev;
      box= &sys->trans.viewport; 
      if (yn>=box->ymin && yn<=box->ymax) {
	tmp= xn-box->xmax;
	if(tmp<min && tmp>0) { min= tmp; iSystem= i; thesys= sys; }
	tmp= box->xmin-xn;
	if(tmp<min && tmp>0) { min= tmp; iSystem= i; thesys= sys; }
      }
      if (xn>=box->xmin && xn<=box->xmax) {
	tmp= yn-box->ymax;
	if(tmp<min && tmp>0) { min= tmp; iSystem= i; thesys= sys; }
	tmp= box->ymin-yn;
	if(tmp<min && tmp>0) { min= tmp; iSystem= i; thesys= sys; }
      }
    }
  }
  if (iSystem) {
    sys= thesys;
    *system= sys;
    FindCoordinates(sys, xn, yn, xr, yr);
  } else {
    *system= 0;
    *xr= xn;  *yr= yn;
  }
  return iSystem;
}

static void Find1System(FXEngine *fxEngine, Drawing *drawing, int iSystem,
			int x, int y, GeSystem **system,
			GpReal *xr, GpReal *yr)
{
  GpXYMap *map= &fxEngine->xe.e.map; /* NDC->VDC (x,y) mapping */
  GpReal xn= ((GpReal)x - map->x.offset)/map->x.scale;
  GpReal yn= ((GpReal)y - map->y.offset)/map->y.scale;
  GeSystem *sys= GetSystemN(drawing, iSystem);
  if (sys && (!(sys->rescan || sys->unscanned>=0) ||
	      !GdScan(sys))) {
    FindCoordinates(sys, xn, yn, xr, yr);
    *system= sys;
  } else {
    *xr= xn;
    *yr= yn;
    *system= 0;
  }
}

static GeSystem *GetSystemN(Drawing *drawing, int n)
{
  if (n<=0 || n>drawing->nSystems) return 0;
  else {
    GeSystem *sys= drawing->systems;
    while (--n) sys= (GeSystem *)sys->el.next;
    return (sys && sys->elements)? sys : 0;
  }
}

static int FindAxis(GeSystem *system, GpReal x, GpReal y)
{
  if (system) {
    int xout= (x==system->trans.window.xmin || x==system->trans.window.xmax);
    int yout= (y==system->trans.window.ymin || y==system->trans.window.ymax);
    if (xout) return yout? 3 : 2;
    else return yout? 1 : 3;
  }
  return 0;
}

static void FindCoordinates(GeSystem *system, GpReal xNDC, GpReal yNDC,
			    GpReal *xWC, GpReal *yWC)
{
  GpReal x, y;
  GpXYMap map;
  GpSetMap(&system->trans.viewport, &system->trans.window, &map);
  x= map.x.scale*xNDC + map.x.offset;
  y= map.y.scale*yNDC + map.y.offset;
  if (system->trans.window.xmin<system->trans.window.xmax) {
    if (x<system->trans.window.xmin) x= system->trans.window.xmin;
    else if (x>system->trans.window.xmax) x= system->trans.window.xmax;
  } else {
    if (x<system->trans.window.xmax) x= system->trans.window.xmax;
    else if (x>system->trans.window.xmin) x= system->trans.window.xmin;
  }
  if (system->trans.window.ymin<system->trans.window.ymax) {
    if (y<system->trans.window.ymin) y= system->trans.window.ymin;
    else if (y>system->trans.window.ymax) y= system->trans.window.ymax;
  } else {
    if (y<system->trans.window.ymax) y= system->trans.window.ymax;
    else if (y>system->trans.window.ymin) y= system->trans.window.ymin;
  }
  *xWC= x;
  *yWC= y;
}

static GpReal GetFormat(char *format, GpReal w,
			GpReal wmin, GpReal wmax, int isLog)
{
  GpReal delta;

  if (isLog) {
    w= exp10(w);
    wmin= exp10(wmin);
    wmax= exp10(wmax);
  }
  delta= wmax-wmin;
  if (delta<0.0) delta= -delta;
  if (wmin<0.0) wmin= -wmin;
  if (wmax<0.0) wmax= -wmax;
  if (wmax<wmin) {
    wmax= wmin;
    wmin-= delta;
  }
  if (!isLog) wmin= wmax; /* only for deciding between e and f format */

  if (wmax>=10000.0 || wmin<0.001) {
    /* use e format */
    int n;
    if (delta>0.0) n= 3+(int)log10(wmax/delta + 0.50001);
    else n= 3;
    sprintf(format, "%%%d.%de", 8+n, n);
  } else {
    /* use f format */
    int pre, post;
    if (wmax>1000.0) pre= 4;
    else if (wmax>100.0) pre= 3;
    else if (wmax>10.0) pre= 2;
    else pre= 1;
    if (delta<0.1 && delta>0.0) post= 3-(int)log10(delta);
    else if (isLog && wmin<0.01) post= 5;
    else post= 4;
    sprintf(format, "%%%d.%df", pre+post+3, post);
  }

  return w;
}

/* ------------------------------------------------------------------------ */

static void DrawRubber(FXEngine *fxEngine, int x, int y)
{
  GC gc= fxEngine->gc;
  Display *dpy;
  int iPass= 2;

  if (!fxEngine->xe.xscr) return;

  dpy= fxEngine->xe.xscr->display;
  XSetFunction(dpy, gc, GXxor);
  XSetForeground(dpy, gc, rubberPixel);
  XSetLineAttributes(dpy, gc, 0, LineSolid, CapButt, JoinMiter);

  /* first undraw previous box or line, then draw new box or line */
  while (iPass--) {
    if (anchorX!=oldX || anchorY!=oldY) {
      if (rubberBanding==1) {
	/* this is a rubber box */
	int xul= oldX>anchorX? anchorX : oldX;
	unsigned int width= oldX>anchorX? oldX-anchorX : anchorX-oldX;
	int yul= oldY>anchorY? anchorY : oldY;
	unsigned int height= oldY>anchorY? oldY-anchorY : anchorY-oldY;
	XDrawRectangle(dpy, fxEngine->xe.graphics, gc,
		       xul, yul, width, height);
      } else {
	/* this is a rubber line */
	XDrawLine(dpy, fxEngine->xe.graphics, gc,
		  anchorX, anchorY, oldX, oldY);
      }
    }
    oldX= x;
    oldY= y;
  }

  XSetLineAttributes(dpy, gc, 3, LineSolid, CapButt, JoinMiter);
  XSetForeground(dpy, gc, fxEngine->xe.xscr->stdColors[1].pixel);
  XSetFunction(dpy, gc, GXcopy);
}

static void SetRubberPixel(FXEngine *fxEngine, int which)
{
  GxScreen *xscr= fxEngine->xe.xscr;
  if (xscr) {
    unsigned long bg2fg= xscr->stdColors[1].pixel^xscr->stdColors[0].pixel;
    if (which==0) rubberPixel= bg2fg;
    else if (which==1) rubberPixel= bg2fg==0x81L? 1L : 0x81L;
    else rubberPixel= bg2fg==0xff? 1L : 0xffL;
  }
}

/* ------------------------------------------------------------------------ */

extern int GxBasicXHandler(XEvent *event);

int GxPointClick(Engine *engine, int style, int system,
		 int (*CallBack)(Engine *engine, int system,
				 int release, GpReal x, GpReal y,
				 int butmod, GpReal xn, GpReal yn))
{
  XEngine *xeng= GisXEngine(engine);
  Display *dpy;
  int fd;
  if (!xeng || !xeng->xscr) return 1;
  dpy= xeng->xscr->display;
  fd= ConnectionNumber(dpy);  /* to check against DispatchNext value */

  /* set up state variables for point-and-click sequence */
  if (!(PtClCallBack= CallBack)) return 1;
  if (style==1 || style==2) ptClStyle= style;
  else ptClStyle= 0;
  if (system<0) ptClSystem= -1;
  else ptClSystem= system;
  ptClCount= 2;

  while (DispatchNext()==fd && PtClCallBack && ptClCount);
  PtClCallBack= 0;

  return ptClCount>0;
}

/* ------------------------------------------------------------------------ */
