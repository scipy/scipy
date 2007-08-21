/*
 * XBASIC.H
 *
 * $Id$
 *
 * Declare the basic play engine for GIST.
 *
 */
/*    Copyright (c) 2000.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef XBASIC_H
#define XBASIC_H

#include "gist.h"
#include "engine.h"
#include "play.h"

typedef struct XEngine XEngine;
struct XEngine {
  Engine e;

  /* --------------- Specific to XEngine ------------------- */

  p_scr *s;
  p_win *win;
  int width, height;  /* of (virtual page) graphics window */
  int wtop, htop;     /* of actual top-level window */
  int topMargin;   /* height of top menu bar, if any */
  int leftMargin;  /* width of left menu bar, if any */
  int x, y;        /* position of graphics relative to win */
  int dpi;         /* resolution of X window (dots per inch, 75 or 100) */
  int mapped, clipping;

  /* if w!=win, this is animation mode */
  p_win *w;
  int a_width, a_height;        /* of animation Pixmap */
  int a_x, a_y;                 /* where it goes on graphics window */
  GpTransform swapped;          /* transform for graphics window while
                                 * in animation mode */

  /* if non-zero, these handlers can deal with input events */
  void (*HandleExpose)(Engine *engine, Drauing *drawing, int *xy);
  void (*HandleClick)(Engine *e,int b,int md,int x,int y, unsigned long ms);
  void (*HandleMotion)(Engine *e,int md,int x,int y);
  void (*HandleKey)(Engine *e,int k,int md);
};

/* GxBasic creates the basic top level window for a GpBXEngine.
   The window manager properties are filled in.
   Use the DefaultTopSize macro to set the width and
   height appropriate for a particular resolution (dpi).
   DefaultTopSize is 6 inches (450 pixels at 75 dpi, 600 at 100 dpi).  */
extern p_scr *GxBasic(char *name, char *display, int width, int height,
                      p_win **win);

extern p_scr *g_connect(char *displayName);
extern void g_disconnect(p_scr *s);

extern int gx75width, gx100width;    /* defaults are 450 and 600 pixels */
#define DefaultTopWidth(dpi) \
  (gx75width<gx100width?((dpi)*gx100width)/100:gx100width)
extern int gx75height, gx100height;  /* defaults are 450 and 600 pixels */
#define DefaultTopHeight(dpi) \
  (gx75width<gx100width?((dpi)*gx100height)/100:gx100height)
#define PixelsPerNDC(dpi) ((dpi)/ONE_INCH)

/* GxEngine creates an XEngine and adds it to the list of GIST engines.
   The top window will generally be smaller than the graphics
   window created by GxEngine; specific engines are responsible for
   scrolling of the graphics window relative to the top window, although
   the initial location is passed in via (x, y).  The size argument is
   sizeof(XEngine), or sizeof some derived engine class.  */
extern XEngine *GxEngine(p_scr *s, char *name, GpTransform *toPixels,
                         int x, int y,
                         int topMargin, int leftMargin, long size);

/* GxInput sets optional event handlers, and calls XSelectInput with
   the given eventMask.  HandleExpose, if non-zero, will be called
   instead of redrawing the Drauing associated with the Engine, which
   is the default action.  HandleResize, if non-zero, will be called
   instead of the default resize action (which is to recenter the
   graphics window).  HandleOther, if non-zero, will be called for
   keyboard, button, or other events not recogized by the default
   handler.  */
extern int GxInput(Engine *engine,
                   void (*HandleExpose)(Engine *, Drauing *, int *),
                   void (*HandleClick)(Engine *,int,int,int,int,unsigned long),
                   void (*HandleMotion)(Engine *,int,int,int),
                   void (*HandleKey)(Engine *,int,int));

extern XEngine *GisXEngine(Engine *engine);

extern void GxExpose(Engine *engine, Drauing *drawing, int *xy);
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

extern int g_rgb_read(Engine *eng, GpColor *rgb, long *nx, long *ny);

#endif
