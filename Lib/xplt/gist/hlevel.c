/*
 * HLEVEL.C
 *
 * $Id$
 *
 * Define routines for recommended GIST interactive interface
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "hlevel.h"

#ifndef NO_XLIB
#ifndef ANIMATE_H
#include "xbasic.h"
#else
#include ANIMATE_H
#endif
#else
#include "engine.h"
static int GxAnimate(Engine *engine, GpBox *viewport);
static int GxStrobe(Engine *engine, int clear);
static int GxDirect(Engine *engine);
static int GxAnimate(Engine *engine, GpBox *viewport)
{ return 0; }
static int GxStrobe(Engine *engine, int clear)
{ return 0; }
static int GxDirect(Engine *engine)
{ return 0; }
#endif

static void UpdateOrRedraw(int changesOnly);

/* Here is a kludge to detect when GdDraw is about to be called to
   walk a Drawing.  The how argument has these meanings:
     even: called just before GdDraw
     odd:  called after following flush (or other related actions)
     0/1:  called from UpdateOrRedraw (just before prompt)
     2/3:  called from GhFMA (explicit frame advance)
     4/5:  called from GhHCP (with hardcopy Engine, not display engine)
 */
extern void (*gdraw_hook)(Engine *display, int how);
void (*gdraw_hook)(Engine *display, int how)= 0;

/* ------------------------------------------------------------------------ */
/* See README for description of these control functions */

GhDevice ghDevices[8];

Engine *hcpDefault= 0;

static int currentDevice= -1;

static int hcpOn= 0;
static int animateOn= 0;

static int fmaCount= 0;

static void UpdateOrRedraw(int changesOnly)
{
  Engine *display= currentDevice<0? 0 : ghDevices[currentDevice].display;
  if (!display) return;
  GpPreempt(display);
  if (gdraw_hook) gdraw_hook(display, 0);
  GdDraw(changesOnly);
  GpFlush(0);
  if (gdraw_hook) gdraw_hook(display, 1);
  GpPreempt(0);
}

void GhBeforeWait(void)
{
  if (currentDevice<0 || !ghDevices[currentDevice].display ||
      animateOn) return;  /* nothing happens in animate mode until
			     explicit call to GhFMA */
  UpdateOrRedraw(1);
}

void GhFMA(void)
{
  Engine *display;
  Engine *hcp= 0;

  if (currentDevice<0) return;
  display= ghDevices[currentDevice].display;
  if (animateOn && !display) animateOn= 0;

  if (hcpOn) {
    hcp= ghDevices[currentDevice].hcp;
    if (!hcp) hcp= hcpDefault;
    if (hcp) GpActivate(hcp);
  }

  if (gdraw_hook) gdraw_hook(display, 2);
  GdDraw(1);
  if (hcpOn && hcp && ghDevices[currentDevice].doLegends)
    GdDrawLegends(hcp);
  if (animateOn) GxStrobe(display, 1);
  GpFlush(0);
  if (animateOn!=1) GdClear(0);
  else GdClearSystem();
  if (gdraw_hook) gdraw_hook(display, 3);

  if (hcpOn && hcp) {
    GpClear(hcp, CONDITIONALLY);
    GpDeactivate(hcp);
  }

  ghDevices[currentDevice].fmaCount++;
  if (++fmaCount > 100) {  /* clean house once in a while */
    fmaCount= 0;
    GaFreeScratch();
  }
}

void GhRedraw(void)
{
  UpdateOrRedraw(-1);
}

void GhHCP(void)
{
  Engine *hcp= currentDevice<0? 0 : ghDevices[currentDevice].hcp;
  if (!hcp) hcp= hcpDefault;
  if (!hcp) return;
  GpPreempt(hcp);
  if (gdraw_hook) gdraw_hook(hcp, 4);
  GdDraw(0);
  /* NB- must be very careful not to Preempt twice with GdDrawLegends */
  if (ghDevices[currentDevice].doLegends) GdDrawLegends(0);
  GpClear(0, ALWAYS);
  GpFlush(0);
  if (gdraw_hook) gdraw_hook(hcp, 5);
  GpPreempt(0);
}

void GhFMAMode(int hcp, int animate)
{
  /* 0 off, 1 on, 2 no change, 3 toggle */
  if (hcp&2) hcpOn^= (hcp&1);   /* if 2 bit, XOR operation */
  else hcpOn= (hcp&1);          /* else, COPY operation */

  if ((animate&3)!=2) {
    Engine *display= currentDevice<0? 0 : ghDevices[currentDevice].display;
    if (!display) return;

    if ((animate&2) || (!animateOn)!=(!(animate&1))) {
      /* animation mode will actually change */
      animateOn= !animateOn;
      if (animateOn) {
	GpBox aport;
	GpBox *port= GdClearSystem();
	/* Can only animate using GdClearSystem if there is a current
	   system; GdClearSystem returns appropriate box to animate,
	   or 0 if can't.
	   If no current system, then animate entire picture using
	   ordinary GpClear.  */
	aport.xmin= 0.0;
	aport.xmax= 2.0;
	aport.ymin= 0.0;
	aport.ymax= 2.0;
	if (!port) {
	  port= &aport;
	  animateOn= 2;
	}
	if (GxAnimate(display, port)) animateOn= 0;
      } else {
	GxDirect(display);
      }
    }
  }
}

int GhSetPlotter(int number)
{
  if (number<0 || number>7) return 1;

  if (currentDevice>=0) {
    if (ghDevices[currentDevice].display) {
      GdSetDrawing(ghDevices[currentDevice].drawing);
      GhBeforeWait();
      GpDeactivate(ghDevices[currentDevice].display);
    }
    if (ghDevices[currentDevice].hcp)
      GpDeactivate(ghDevices[currentDevice].hcp);
  }
  if (hcpDefault) GpDeactivate(hcpDefault);

  currentDevice= number;
  if (ghDevices[number].display) GpActivate(ghDevices[number].display);
  return GdSetDrawing(ghDevices[number].drawing);
}

int GhGetPlotter(void)
{
  return currentDevice;
}

#ifndef NO_XLIB
/* xbasic.c supplies a hook in its error handlers to allow the hlevel
   to clear its display devices */
static void ShutDownDev(Engine *engine);
extern void (*HLevelHook)(Engine *engine);

static void ShutDownDev(Engine *engine)
{
  int i;

  if (hcpDefault==engine) hcpDefault= 0;
  for (i=0 ; i<8 ; i++) {
    if (ghDevices[i].display==engine) {
      if (i==currentDevice) currentDevice= -1;
      ghDevices[i].display= 0;
    }
    if (ghDevices[i].hcp==engine) {
      if (!ghDevices[i].display && i==currentDevice) currentDevice= -1;
      ghDevices[i].hcp= 0;
    }
  }
}

int GhSetXHandler(void (*ErrHandler)(char *errMsg))
{
  GpSetXHandler(ErrHandler);
  HLevelHook= &ShutDownDev;
  return 0;
}

void GhWaitDisplay(void)
{
  if (currentDevice<0) return;
  GxWaitForExpose(ghDevices[currentDevice].display);
}

#else
/* ARGSUSED */
int GhSetXHandler(void (*ErrHandler)(char *errMsg))
{
  /* no-op */
  return 0;
}

void GhWaitDisplay(void)
{
  /* no-op */
  return;
}
#endif

/* ------------------------------------------------------------------------ */
/* Default management */

static GpLineAttribs lDefault= { FG_COLOR, L_SOLID, 1.0 };
static GpMarkerAttribs mDefault= { FG_COLOR, 0, 1.0 };
static GpFillAttribs fDefault= { FG_COLOR, F_SOLID, 0, 0.01, 0.01, 0.0, 0.0 };
static GpTextAttribs tDefault= { FG_COLOR, 0, 0.0156,
				   TX_RIGHT, TH_NORMAL, TV_NORMAL };
static GaLineAttribs dlDefault= { 0, 0, 0, 0.16, 0.14, 0,
				    0.13, 0.11375, 1.0, 1.0 };
static GaVectAttribs vectDefault= { 0, 0.125 };
static GpLineAttribs edgeDefault= { FG_COLOR, L_NONE, 1.0 };

void GhGetLines(void)
{
  gistA.l= lDefault;
  gistA.m= mDefault;
  gistA.dl= dlDefault;
}

void GhGetText(void)
{
  gistA.t= tDefault;
}

void GhGetMesh(void)
{
  gistA.l= lDefault;
}

void GhGetVectors(void)
{
  gistA.l= lDefault;
  gistA.f= fDefault;
  gistA.vect= vectDefault;
}

void GhGetFill(void)
{
  gistA.e= edgeDefault;
}

void GhSetLines(void)
{
  lDefault= gistA.l;
  mDefault= gistA.m;
  mDefault.type= 0;    /* never a default marker */
  dlDefault= gistA.dl;
}

void GhSetText(void)
{
  tDefault= gistA.t;
}

void GhSetMesh(void)
{
  lDefault= gistA.l;
}

void GhSetVectors(void)
{
  lDefault= gistA.l;
  fDefault= gistA.f;
  vectDefault= gistA.vect;
}

void GhSetFill(void)
{
  edgeDefault= gistA.e;
}

/* ------------------------------------------------------------------------ */

void GhDumpColors(int n, int hcp, int private)
{
  Engine *engine;
  if (n>=0 && n<8) {
    if (hcp) engine= ghDevices[n].hcp;
    else engine= ghDevices[n].display;
  } else {
    engine= hcpDefault;
  }
  if (engine) GpDumpColors(engine, private);
}

int GhGetColorMode(Engine *engine)
{
  return engine->colorMode;
}

void GhSetPalette(int n, GpColorCell *palette, int nColors)
{
  if (ghDevices[n].display && ghDevices[n].display->palette!=palette) {
    GpSetPalette(ghDevices[n].display, palette, nColors);
    if (!ghDevices[n].display->colorMode) GhRedraw();
  }
  if (ghDevices[n].hcp && ghDevices[n].hcp->palette!=palette)
    GpSetPalette(ghDevices[n].hcp, palette, nColors);
}

int GhReadPalette(int n, const char *gpFile,
		  GpColorCell **palette, int maxColors)
{
  int paletteSize= 0;
  if (ghDevices[n].display) {
    paletteSize= GpReadPalette(ghDevices[n].display, gpFile,
			       &ghDevices[n].display->palette, maxColors);
    if (ghDevices[n].hcp)
      GpSetPalette(ghDevices[n].hcp,
		   ghDevices[n].display->palette,
		   ghDevices[n].display->nColors);
    if (palette) *palette= ghDevices[n].display->palette;
    /* override the (possibly truncated) value returned by GpReadPalette */
    paletteSize= ghDevices[n].display->nColors;
    if (!ghDevices[n].display->colorMode) GhRedraw();
  } else if (ghDevices[n].hcp) {
    paletteSize= GpReadPalette(ghDevices[n].hcp, gpFile,
			       &ghDevices[n].hcp->palette, maxColors);
    if (palette) *palette= ghDevices[n].hcp->palette;
    /* override the (possibly truncated) value returned by GpReadPalette */
    paletteSize= ghDevices[n].hcp->nColors;
  }
  return paletteSize;
}

int GhGetPalette(int n, GpColorCell **palette)
{
  *palette= 0;
  if (n<0 || n>7) return 0;
  if (ghDevices[n].display)
    return GpGetPalette(ghDevices[n].display, palette);
  else if (ghDevices[n].hcp)
    return GpGetPalette(ghDevices[n].hcp, palette);
  else
    return 0;
}

void GhDeletePalette(int n)
{
  GpColorCell *palette= 0;
  if (n<0 || n>7) return;
  if (ghDevices[n].display) palette= ghDevices[n].display->palette;
  else if (ghDevices[n].hcp) palette= ghDevices[n].hcp->palette;
  if (palette) {
    int i;

    /* clear palette for this device */
    if (ghDevices[n].display)
      GpSetPalette(ghDevices[n].display, (GpColorCell *)0, 0);
    if (ghDevices[n].hcp)
      GpSetPalette(ghDevices[n].hcp, (GpColorCell *)0, 0);

    /* free the palette if there are no other references to it */
    for (i=0 ; i<8 ; i++)
      if ((ghDevices[i].display &&
	   ghDevices[i].display->palette==palette) ||
	  (ghDevices[i].hcp &&
	   ghDevices[i].hcp->palette==palette)) break;
    if (i>=8) {
      if (hcpDefault && palette==hcpDefault->palette)
	GpSetPalette(hcpDefault, (GpColorCell *)0, 0);
      GmFree(palette);
    }
  }
}

void SetHCPPalette(void)
{
  if (hcpDefault && currentDevice>=0) {
    GpColorCell *palette= 0;
    int nColors= 0;
    if (ghDevices[currentDevice].display) {
      palette= ghDevices[currentDevice].display->palette;
      nColors= ghDevices[currentDevice].display->nColors;
    } else if (ghDevices[currentDevice].hcp) {
      palette= ghDevices[currentDevice].hcp->palette;
      nColors= ghDevices[currentDevice].hcp->nColors;
    }
    GpSetPalette(hcpDefault, palette, nColors);
  }
}

/* ------------------------------------------------------------------------ */
