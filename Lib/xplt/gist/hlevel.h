/*
 * HLEVEL.H
 *
 * $Id$
 *
 * Declare routines for recommended GIST interactive interface
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef HLEVEL_H
#define HLEVEL_H

#include "gist.h"

/* See README for description of these control functions */

extern void GhBeforeWait(void);
extern void GhFMA(void);
extern void GhRedraw(void);
extern void GhHCP(void);
extern void GhFMAMode(int hcp, int animate); /* 0 off, 1 on, 2 nc, 3 toggle */

/* Ensure that display engine is ready to draw -- waiting for expose
   event if necessary.  */
extern void GhWaitDisplay(void);

/* The pldevice call should create the necessary engines, and set their
   Engine pointers in ghDevices, then call GhSetPlotter to set the
   current device, deactivating the old device.  GhSetPlotter returns
   0 if successful, 1 if neither display nor hcp engine has been defined
   for the requested device.  GhGetPlotter returns the number of the
   current plot device, or -1 if none has been set.  */

typedef struct GhDevice GhDevice;
struct GhDevice {
  Drawing *drawing;
  Engine *display, *hcp;
  int doLegends;
  int fmaCount;
  void *hook;
};

/* Allow up to 8 windows per application */
extern GhDevice ghDevices[8];

extern int GhSetPlotter(int number);
extern int GhGetPlotter(void);

/* The default hardcopy device is used for hcp commands whenever the
   current device has no hcp engine of its own.  */
extern Engine *hcpDefault;

extern void GhDumpColors(int n, int hcp, int private);
extern int GhGetColorMode(Engine *engine);
extern void GhSetPalette(int n, GpColorCell *palette, int nColors);
extern int GhReadPalette(int n, const char *gpFile,
			 GpColorCell **palette, int maxColors);
extern int GhGetPalette(int n, GpColorCell **palette);
extern int GhGetColorMode(Engine *engine);
extern void SetHCPPalette(void);
extern void GhDeletePalette(int n);

/* A high-level error handler takes down an X-window before calling
   the user-installed error handler.  This prevents a huge blast of
   errors when a window is detroyed bby a window manager (for example),
   but is obviously a litle more fragile than a smart error handler
   could be.  */
extern int GhSetXHandler(void (*XHandler)(char *msg));

/* For each of the D level drawing primitives, a set of
   default parameter settings is maintained, and can be installed
   with the appropriate GhGet routine.  The GhSet routines set the
   defaults themselves.  GdCells does not use any attributes,
   and GdContours uses the same attributes as GdLines.
   GdFillMesh uses line attributes for edges (if any).  */
extern void GhGetLines(void);
extern void GhGetText(void);
extern void GhGetMesh(void);
extern void GhGetVectors(void);
extern void GhGetFill(void);

extern void GhSetLines(void);
extern void GhSetText(void);
extern void GhSetMesh(void);
extern void GhSetVectors(void);
extern void GhSetFill(void);

/* The GpFXEngine (fancy X engine) has controls for the zoom factor
   and a function for initiating a point-and-click sequence.  */

extern GpReal gxZoomFactor;   /* should be >1.0, default is 1.5 */

/* The GxPointClick function initiates an interactive point-and-click
   session with the window -- it will not return until a button has
   been pressed, then released.  It returns non-zero if the operation
   was aborted by pressing a second button before releasing the first.
     engine --   an X engine whose display is to be used
     style --    1 to draw a rubber box, 2 to draw a rubber line,
                 otherwise, no visible indication of operation
     system --   system number to which the world coordinates should
                 be transformed, or -1 to use the system under the
		 pointer -- the release coordinates are always in the
		 same system as the press coordinates
     CallBack -- function to be called twice, first when the button is
                 pressed, next when it is released -- operation will
		 be aborted if CallBack returns non-zero
		 Arguments passed to CallBack:
		   engine  -- in which press/release occurred
		   system  -- system under pointer, if above system -1
		   release -- 0 on press, 1 on release
		   x, y    -- coordinates of pointer relative to system
		   butmod  -- 1 - 5 on press to tell which button
		              mask to tell which modifiers on release:
			      1 shift, 2 lock, 4 control, 8 - 128 mod1-5
		   xn, yn  -- NDC coordinates of pointer
 */
extern int GxPointClick(Engine *engine, int style, int system,
			int (*CallBack)(Engine *engine, int system,
					int release, GpReal x, GpReal y,
					int butmod, GpReal xn, GpReal yn));

#endif
