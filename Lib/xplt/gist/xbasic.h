/*
 * XBASIC.H
 *
 * $Id$
 *
 * Declare the basic X windows engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef XBASIC_H
#define XBASIC_H

#include "gist.h"
#include "engine.h"
#include "xicky.h"

typedef struct GxLastOp GxLastOp;
struct GxLastOp {
  int color;

  GpReal lwidth;
  int ltype;
  int ljoin;

  int tfont;
  GpReal tsize;
  int fontID;    /* X font id code for best match to tfont, tsize */

  /* tiles and stipples may eventually be mapped to patterns and hatches */
};

typedef void (*GxHandler)(Engine *, Drawing *, XEvent *);

typedef struct XEngine XEngine;
struct XEngine {
  Engine e;

  /* --------------- Specific to XEngine ------------------- */

  GxScreen *xscr;
  GxDisplay *xdpy;
  Window top, graphics;
  unsigned int width, height;  /* of graphics window */
  int topMargin;   /* height of top menu bar, if any */
  int leftMargin;  /* width of left menu bar, if any */
  int x, y;        /* position of graphics relative to top (<0 usually) */
  int dpi;         /* resolution of X window (dots per inch, 75 or 100) */
  int mapped;

  GC gc;
  GxLastOp lastOp;  /* reflects current gc settings */

  /* If drawable!=graphics, this is animation mode */
  Drawable drawable;
  unsigned int aWidth, aHeight; /* of animation Pixmap */
  int graphicsX, graphicsY;     /* where it goes on graphics window */
  GC gca;                       /* GC for XCopyArea (for clipping only) */
  GpTransform swapped;          /* transform for graphics window while
				   in animation mode */

  /* Current mapping from GpColor index to X pixels, if any */
  int nColors;
  unsigned long *pixelMap;

  /* Private colormap, if any.  GxExactColors guaratees that this will
     be non-zero if it exists, so zero means there is no private map.  */
  Colormap private;

  /* If non-zero, these handlers can deal with X input.  Set using
     GxInput.  */
  GxHandler HandleExpose, HandleResize, HandleOther;
};

/* GxBasic creates the basic top level window for a GpBXEngine.
   The window manager properties are filled in.
   Use the DefaultTopSize macro to set the width and
   height appropriate for a particular resolution (dpi).
   DefaultTopSize is 6 inches (450 pixels at 75 dpi, 600 at 100 dpi).  */
extern GxScreen *GxBasic(char *name, char *display, int width, int height,
			 Window *top);

extern int gx75width, gx100width;    /* defaults are 450 and 600 pixels */
#define DefaultTopWidth(dpi) ((dpi)<88? gx75width : gx100width)
extern int gx75height, gx100height;  /* defaults are 450 and 600 pixels */
#define DefaultTopHeight(dpi) ((dpi)<88? gx75height : gx100height)
#define PixelsPerNDC(dpi) ((dpi)<88? 75.0/ONE_INCH : 100.0/ONE_INCH)

/* GxEngine creates an XEngine and adds it to the list of GIST engines.
   The top window will generally be smaller than the graphics
   window created by GxEngine; specific engines are responsible for
   scrolling of the graphics window relative to the top window, although
   the initial location is passed in via (x, y).  The size argument is
   sizeof(XEngine), or sizeof some derived engine class.  */
extern XEngine *GxEngine(char *name, GpTransform *toPixels, GxScreen *xscr, 
			 Window top, int x, int y,
			 int topMargin, int leftMargin, long size);

/* GxInput sets optional event handlers, and calls XSelectInput with
   the given eventMask.  HandleExpose, if non-zero, will be called
   instead of redrawing the Drawing associated with the Engine, which
   is the default action.  HandleResize, if non-zero, will be called
   instead of the default resize action (which is to recenter the
   graphics window).  HandleOther, if non-zero, will be called for
   keyboard, button, or other events not recogized by the default
   handler.  */
extern int GxInput(Engine *engine, GxHandler HandleExpose,
		   GxHandler HandleResize, GxHandler HandleOther,
		   long eventMask);

extern XEngine *GxGetEngine(Display *display, Window window);
extern XEngine *GisXEngine(Engine *engine);

/* If HandleExpose or HandleResize is zero, the basic event handling
   routines are called.  Any replacements should almost certainly call
   these default handlers early in the course of their own processing.
   Note that GxRecenter is not a *GxHandler; the width and height are
   event->xconfigure.width, event->xconfigure.height.  */
extern void GxExpose(Engine *engine, Drawing *drawing, XEvent *event);
extern void GxRecenter(XEngine *xEngine, int width, int height);

/* GxAnimate creates an offscreen pixmap for the specified region of
   the window.  Subsequent drawing takes place on the pixmap, not
   on the graphics window.  GxStrobe copies the pixmap to the screen,
   optionally clearing it for the next frame of the animation.
   The viewport should be large enough to cover everything that will
   change as the animation proceeds, but no larger to get peak speed.
   GxDirect restores the usual direct-to-screen drawing mode.  */
extern int GxAnimate(Engine *engine, GpBox *viewport);
extern int GxStrobe(Engine *engine, int clear);
extern int GxDirect(Engine *engine);

/* As a convenience for creating GCs, here is the required struct */
extern XGCValues gxXGCValues;

/* Unless you wait for the top level window to be exposed, drawing into
   an X engine will silently abort.  */
extern int GxWaitForExpose(Engine *engine);

#endif
