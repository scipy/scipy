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
#include "draw.h"

#include <string.h>

/* possibly should just include stdio.h, math.h */
extern int sprintf(char *s, const char *format, ...);
extern double log10(double);

#ifndef NO_EXP10
  extern double exp10(double);
#else
# define exp10(x) pow(10.,x)
  extern double pow(double,double);
#endif

/* ------------------------------------------------------------------------ */

static void HandleExpose(Engine *engine, Drauing *drawing, int *xy);
static void HandleClick(Engine *e,int b,int md,int x,int y, unsigned long ms);
static void HandleMotion(Engine *e,int md,int x,int y);
static void HandleKey(Engine *e,int k,int md);

static void RedrawMessage(FXEngine *fxe);
static void MovePointer(FXEngine *fxe, Drauing *drawing,
                        int md,int x,int y);
static void PressZoom(FXEngine *fxe, Drauing *drawing,
                      int b,int md,int x,int y, unsigned long ms);
static void ReleaseZoom(FXEngine *fxe, Drauing *drawing,
                        int b,int md,int x,int y, unsigned long ms);

static void redraw_seps(FXEngine *fxe);
static void check_clipping(FXEngine *fxe);
static void RedrawButton(FXEngine *fxe);
static void EnterButton(FXEngine *fxe);
static void LeaveButton(FXEngine *fxe);
static void PressButton(FXEngine *fxe,
                        int b,int md,int x,int y, unsigned long ms);
static void ReleaseButton(FXEngine *fxe, Drauing *drawing,
                          int b,int md,int x,int y, unsigned long ms);

static void ButtonAction(FXEngine *fxe, Drauing *drawing);
static void HighlightButton(FXEngine *fxe);
static void UnHighlightButton(FXEngine *fxe);

static void ResetZoom(FXEngine *fxe);
static void DoZoom(GpReal factor, GpReal w0, GpReal w1,
                   GpReal *wmin, GpReal *wmax);
static void AltZoom(GpReal w0, GpReal w1, GpReal *wmin, GpReal *wmax);
static int FindSystem(FXEngine *fxe, Drauing *drawing, int x, int y,
                      GeSystem **system, GpReal *xr, GpReal *yr);
static void Find1System(FXEngine *fxe, Drauing *drawing, int iSystem,
                        int x, int y, GeSystem **system,
                        GpReal *xr, GpReal *yr);
static GeSystem *GetSystemN(Drauing *drawing, int n);
static int FindAxis(GeSystem *system, GpReal x, GpReal y);
static void FindCoordinates(GeSystem *system, GpReal xNDC, GpReal yNDC,
                            GpReal *xWC, GpReal *yWC);
static GpReal GetFormat(char *format, GpReal w,
                        GpReal wmin, GpReal wmax, int isLog);
static void DrawRubber(FXEngine *fxe, int x, int y);

/* ------------------------------------------------------------------------ */

Engine *
GpFXEngine(char *name, int landscape, int dpi, char *displayName)
{
  p_scr *s = g_connect(displayName);
  int topWidth= DefaultTopWidth(dpi);   /* not including button, message */
  int topHeight= DefaultTopHeight(dpi);
  GpTransform toPixels;
  int x, y;
  FXEngine *fxe;
  int heightButton, widthButton, baseline;
  int width, height, ascent, descent, hints;

  if (!s) return 0;
  descent = p_txheight(s, P_GUI_FONT, 15, &ascent);
  descent -= ascent;
  width = p_txwidth(s, "System", 6, P_GUI_FONT, 15);
  baseline = ascent+2;
  heightButton = baseline+descent+4;  /* leave at least 2 lines
                                       * above and below text */
  widthButton = width+8;  /* leave at least 4 pixels before and after */

  /* set toPixels as by SetXTransform(&toPixels, landscape, dpi) */
  toPixels.viewport = landscape? gLandscape : gPortrait;
  toPixels.window.xmin = 0.0;
  toPixels.window.xmax = PixelsPerNDC(dpi)*toPixels.viewport.xmax;
  toPixels.window.ymin = PixelsPerNDC(dpi)*toPixels.viewport.ymax;
  toPixels.window.ymax = 0.0;
  width = (int)toPixels.window.xmax;
  height = (int)toPixels.window.ymin;
  x = (width-topWidth)/2;
  if (landscape) y = (height-topHeight)/2;
  else y = (width-topHeight)/2;
  if (y<0) y = 0;
  if (x<0) x = 0;
  fxe = (FXEngine *)GxEngine(s, name, &toPixels, -x,-y,heightButton+2,0,
                             sizeof(FXEngine));

  fxe->xe.wtop = topWidth;
  fxe->xe.htop = topHeight;
  /* possibly want optional P_RGBMODEL as well */
  hints = (gist_private_map?P_PRIVMAP:0) | (gist_input_hint?0:P_NOKEY) |
    (gist_rgb_hint?P_RGBMODEL:0);
  fxe->xe.win = fxe->xe.w =
    p_window(s, topWidth, topHeight+heightButton+2, name, P_BG, hints, fxe);
  if (!fxe->xe.win) {
    GpDelEngine(&fxe->xe.e);
    return 0;
  }

  /* set the extra information in the fxe structure */
  fxe->baseline = baseline;
  fxe->heightButton = heightButton;
  fxe->widthButton = widthButton;

  fxe->xmv = fxe->ymv = -1;
  fxe->pressed = 0;
  fxe->buttonState = 0;
  fxe->iSystem = -1;
  strcpy(fxe->msgText, "Press 1, 2, 3 to zoom in, pan, zoom out");
  fxe->msglen = 0;
  fxe->zoomState = fxe->zoomSystem = fxe->zoomAxis = 0;
  fxe->zoomX = fxe->zoomY = 0.0;

  /* set the fancy event handlers */
  GxInput((Engine *)fxe, &HandleExpose, &HandleClick, &HandleMotion,
          &HandleKey);

  return (Engine *)fxe;
}

/* ------------------------------------------------------------------------ */

static void
HandleExpose(Engine *engine, Drauing *drawing, int *xy)
{
  FXEngine *fxe = (FXEngine *)engine;
  redraw_seps(fxe);
  RedrawButton(fxe);
  RedrawMessage(fxe);
  if (!xy || xy[3]>=fxe->xe.topMargin)
    GxExpose(engine, drawing, xy);
}

static int xf_region(FXEngine *fxe, int x, int y);

static void
HandleClick(Engine *e,int b,int md,int x,int y, unsigned long ms)
{
  FXEngine *fxe = (FXEngine *)e;
  if (x!=fxe->xmv || y!=fxe->ymv)
    HandleMotion(e, md, x, y);
  if (md & (1<<(b+2))) {  /* this is button release */
    if (fxe->pressed==1)
      ReleaseButton(fxe, fxe->xe.e.drawing, b, md, x, y, ms);
    else if (fxe->pressed==2)
      ReleaseZoom(fxe, fxe->xe.e.drawing, b, md, x, y, ms);

  } else {                /* this is button press */
    int region = xf_region(fxe, x, y);
    if (region == 3)
      PressZoom(fxe, fxe->xe.e.drawing, b, md, x, y, ms);
    else if (region == 1)
      PressButton(fxe, b, md, x, y, ms);
  }
}

static void
HandleMotion(Engine *e,int md,int x,int y)
{
  FXEngine *fxe = (FXEngine *)e;
  int region = xf_region(fxe, x, y);
  int old_region = xf_region(fxe, fxe->xmv, fxe->ymv);
  if (!fxe->pressed) {
    if (region != old_region) {
      if (old_region==3) p_cursor(fxe->xe.win, P_SELECT);
      else if (region==3) p_cursor(fxe->xe.win, P_CROSSHAIR);
    }
    if (region==1) {
      if (old_region!=1) EnterButton(fxe);
    } else {
      if (old_region==1) LeaveButton(fxe);
    }
  }
  if (fxe->pressed==1) {
    if (region!=1 && old_region==1)
      LeaveButton(fxe);
  } else if (region==3 || fxe->pressed==2) {
    MovePointer(fxe, fxe->xe.e.drawing, md, x, y);
  }
  fxe->xmv = x;
  fxe->ymv = y;
}

static int
xf_region(FXEngine *fxe, int x, int y)
{
  if (y>=0 && x>=0 && x<fxe->xe.leftMargin+fxe->xe.wtop) {
    y -= fxe->xe.topMargin;
    if (y < 0)
      return (x < fxe->widthButton)? 1 : 2;
    else if (y < fxe->xe.htop)
      return 3;
  }
  return 0;
}

static void
HandleKey(Engine *e,int k,int md)
{
  FXEngine *fxe = (FXEngine *)e;
  if (!fxe->pressed && !fxe->buttonState) {
    if (!fxe->msglen) fxe->msgText[0] = '\0';
    if (k>=' ' && k<'\177') {
      /* append printing characters */
      if (fxe->msglen>=94) fxe->msglen = 0;  /* crude overflow handling */
      fxe->msgText[fxe->msglen++] = k;
      fxe->msgText[fxe->msglen] = '\0';
    } else if (k=='\177' || k=='\010') {
      /* delete or backspace kills char */
      if (fxe->msglen)
        fxe->msgText[--fxe->msglen] = '\0';
    } else if (k=='\027') {
      /* C-w kills word */
      int n = fxe->msglen;
      char c = n? fxe->msgText[n-1] : '\0';
      if (c<'0' || (c>'9' && c<'A') || (c>'Z' && c<'a' && c!='_') || c>'z') {
        if (fxe->msglen)
          fxe->msgText[--fxe->msglen] = '\0';
      } else {
        while (--n) {
          c = fxe->msgText[n-1];
          if (c<'0' || (c>'9' && c<'A') || (c>'Z' && c<'a' && c!='_') ||
              c>'z') break;
        }
        fxe->msgText[n] = '\0';
        fxe->msglen = n;
      }
    } else if (k=='\025') {
      /* C-u kills line */
      fxe->msgText[0] = '\0';
      fxe->msglen = 0;
    } else if (k=='\012' || k=='\015') {
      /* linefeed or carriage return sends line to interpreter */
      int n = fxe->msglen;
      fxe->msgText[n] = '\0';
      if (g_on_keyline) g_on_keyline(fxe->msgText);
      fxe->msgText[0] = '\0';
      fxe->msglen = 0;
    }
    RedrawMessage(fxe);
  }
}

/* ------------------------------------------------------------------------ */

static void
check_clipping(FXEngine *fxe)
{
  p_win *w = fxe->xe.win;
  if (w) {
    if (fxe->xe.clipping) {
      p_clip(w, 0, 0, 0, 0);
      fxe->xe.clipping = 0;
    }
  }
}

static void
redraw_seps(FXEngine *fxe)
{
  p_win *w = fxe->xe.win;
  if (w) {
    int xpage = (int)fxe->xe.swapped.window.xmax;
    int ypage = (int)fxe->xe.swapped.window.ymin;
    int x[4], y[4];
    check_clipping(fxe);
    p_color(w, P_FG);
    if (fxe->xe.leftMargin+fxe->xe.wtop > xpage) {
      p_pen(w, 4, P_SOLID);
      x[0] = x[1] = xpage+2;
      y[0] = 0;  y[1] = ypage+2;
      p_i_pnts(w, x, y, 2);
      p_segments(w);
    }
    if (fxe->xe.topMargin+fxe->xe.htop > ypage) {
      p_pen(w, 4, P_SOLID);
      x[0] = 0;  x[1] = xpage+2;
      y[0] = y[1] = ypage+2;
      p_i_pnts(w, x, y, 2);
      p_segments(w);
    }
    p_pen(w, 1, P_SOLID);
    x[0] = 0;  x[1] = fxe->xe.wtop;
    y[0] = y[1] = fxe->xe.topMargin-1;
    x[2] = x[3] = fxe->widthButton;
    y[2] = 0;  y[3] = fxe->xe.topMargin-1;
    p_i_pnts(w, x, y, 4);
    p_segments(w);
  }
}

static void
RedrawButton(FXEngine *fxe)
{
  p_win *w = fxe->xe.win;
  if (w) {
    int fg = (fxe->buttonState==2)? P_BG : P_FG;
    int bg = (fxe->buttonState==2)? P_FG : P_BG;
    check_clipping(fxe);
    p_color(w, bg);
    p_rect(w, 0, 0, fxe->widthButton, fxe->xe.topMargin-1, 0);
    if (fxe->buttonState) HighlightButton(fxe);
    else p_color(w, fg);
    p_font(w, P_GUI_FONT, 15, 0);
    p_text(w, 3, fxe->baseline, "System", 6);
  }
}

static void
HighlightButton(FXEngine *fxe)
{
  p_win *w = fxe->xe.win;
  if (w && fxe->buttonState) {
    check_clipping(fxe);
    p_color(w, (fxe->buttonState==2)? P_BG : P_FG);
    p_pen(w, 3, P_SOLID);
    p_rect(w, 1, 1, fxe->widthButton-2, fxe->xe.topMargin-3, 1);
  }
}

static void
UnHighlightButton(FXEngine *fxe)
{
  p_win *w = fxe->xe.win;
  if (w) {
    check_clipping(fxe);
    p_color(w, (fxe->buttonState==2)? P_FG : P_BG);
    p_pen(w, 3, P_SOLID);
    p_rect(w, 1, 1, fxe->widthButton-2, fxe->xe.topMargin-3, 1);
  }
}

static void
EnterButton(FXEngine *fxe)
{
  if (fxe->buttonState == 0) {
    fxe->buttonState = 1;
    HighlightButton(fxe);
  } else if (fxe->buttonState != 0) {
    fxe->buttonState = 0;
    RedrawButton(fxe);
  }
}

static void
LeaveButton(FXEngine *fxe)
{
  int state = fxe->buttonState;
  fxe->buttonState = fxe->pressed = 0;
  if (state==1) UnHighlightButton(fxe);
  else if (state) RedrawButton(fxe);
}

static void
PressButton(FXEngine *fxe,int b,int md,int x,int y, unsigned long ms)
{
  if (!(md&0370) && fxe->buttonState==1) {
    fxe->buttonState = 2;
    fxe->pressed = 1;
    RedrawButton(fxe);
  } else {
    LeaveButton(fxe);
  }
}

static void
ReleaseButton(FXEngine *fxe, Drauing *drawing,
              int b,int md,int x,int y, unsigned long ms)
{
  if (fxe->buttonState==2) {
    fxe->buttonState = 1;
    fxe->pressed = 0;
    RedrawButton(fxe);
    ButtonAction(fxe, drawing);
  } else {
    LeaveButton(fxe);
  }
}

/* ------------------------------------------------------------------------ */

static void ButtonAction(FXEngine *fxe, Drauing *drawing)
{
  int nSystems= drawing? drawing->nSystems : 0;
  int iSystem= fxe->iSystem;
  iSystem++;
  if (iSystem<=nSystems) fxe->iSystem= iSystem;
  else fxe->iSystem= iSystem= -1;
  sprintf(fxe->msgText, "%s%2d", iSystem>=0?"=":":",
          iSystem>=0 ? iSystem : 0);
  RedrawMessage(fxe);
}

/* ------------------------------------------------------------------------ */

static void RedrawMessage(FXEngine *fxe)
{
  p_win *w = fxe->xe.win;
  if (w) {
    char *msg= fxe->msgText;
    int len= (int)strlen(msg);
    check_clipping(fxe);
    /* NOTE: may need to animate this if flickering is annoying */
    p_color(w, P_BG);
    p_rect(w, fxe->widthButton+1, 0, fxe->xe.wtop, fxe->xe.topMargin-2, 0);
    p_color(w, P_FG);
    p_font(w, P_GUI_FONT, 15, 0);
    p_text(w, fxe->widthButton+4, fxe->baseline, msg, len);
  }
}

static char stdFormat[] = "%7.4f";
static int rubberBanding = 0;
static int anchorX, anchorY, oldX, oldY;

static void
MovePointer(FXEngine *fxe, Drauing *drawing,
            int md,int x,int y)
{
  int iSystem = fxe->iSystem;
  int locked, logX = 0, logY = 0;
  char format[24];  /* e.g.- "%s%2d (%11.3e, %11.3e)" */
  char xFormat[16], yFormat[16], *f1, *f2;
  GpReal xWC, yWC;
  GeSystem *system;

  if (!drawing || fxe->buttonState) return;

  /* find the system number and world coordinates */
  if (iSystem >= 0) {
    /* coordinate system is locked */
    Find1System(fxe, drawing, iSystem, x, y, &system, &xWC, &yWC);
    if (!system) iSystem = 0;
    locked = 1;
  } else {
    /* select coordinate system under pointer */
    iSystem = FindSystem(fxe, drawing, x, y, &system, &xWC, &yWC);
    locked = 0;
  }

  if (!fxe->msglen) {
    if (system) {
      logX = system->flags&D_LOGX;
      logY = system->flags&D_LOGY;
      xWC = GetFormat(xFormat, xWC, system->trans.window.xmin,
                      system->trans.window.xmax, logX);
      f1 = xFormat;
      yWC = GetFormat(yFormat, yWC, system->trans.window.ymin,
                      system->trans.window.ymax, logY);
      f2 = yFormat;
    } else {
      f1 = f2 = stdFormat;
    }
    sprintf(format, "%%s%%2d (%s, %s)", f1, f2);
    sprintf(fxe->msgText, format, locked? "=" : ":", iSystem, xWC, yWC);

    RedrawMessage(fxe);
  }
  if (rubberBanding) DrawRubber(fxe, x, y);
}

GpReal gxZoomFactor= 1.5;

static int (*PtClCallBack)(Engine *engine, int system,
                           int release, GpReal x, GpReal y,
                           int butmod, GpReal xn, GpReal yn) = 0;
static int ptClStyle = 0, ptClSystem = 0, ptClCount = 0;

static unsigned int cShapes[3] = { P_EW, P_NS, P_NSEW };

static void
PressZoom(FXEngine *fxe, Drauing *drawing,
          int b,int md,int x,int y, unsigned long ms)
{
  if (!drawing || !fxe->xe.win) return;
  if (fxe->zoomState == 0 && fxe->pressed == 0) {
    /* record button number as zoomState */
    if ((b==1)? (md&P_SHIFT) : (b==2)) fxe->zoomState = 2;
    else if ((b==1)? (md&P_META) : (b==3)) fxe->zoomState = 3;
    else if (b==1) fxe->zoomState = 1;
    if (fxe->zoomState) {
      if (PtClCallBack)
        fxe->zoomState |= 8;
      else if (md&P_CONTROL)
        fxe->zoomState |= 4;
    }
    if (fxe->zoomState) {
      /* record system and axis, x and y, change cursor */
      GeSystem *system;
      int iSystem, axis;
      if (!PtClCallBack || ptClSystem<0) {
        iSystem = FindSystem(fxe, drawing, x, y,
                             &system, &fxe->zoomX, &fxe->zoomY);
        axis = FindAxis(system, fxe->zoomX, fxe->zoomY);
      } else {
        iSystem = ptClSystem;
        Find1System(fxe, drawing, iSystem, x, y,
                    &system, &fxe->zoomX, &fxe->zoomY);
        if (!system) iSystem = ptClSystem = 0;
        axis = 3;
      }
      fxe->zoomSystem = iSystem;

      fxe->pressed = 2;

      if (PtClCallBack) {
        GpXYMap *map = &fxe->xe.e.map;   /* NDC->VDC (x,y) mapping */
        GpReal xNDC = ((GpReal)x - map->x.offset)/map->x.scale;
        GpReal yNDC = ((GpReal)y - map->y.offset)/map->y.scale;
        int button = fxe->zoomState&3;
        ptClCount--;
        if (PtClCallBack(&fxe->xe.e, iSystem, 0,
                         fxe->zoomX, fxe->zoomY,
                         button, xNDC, yNDC)) {
          /* callback routine signals abort */
          ResetZoom(fxe);
          return;
        }
      }

      if ((fxe->zoomAxis = axis)) {
        if (fxe->zoomState<4) {
          p_cursor(fxe->xe.win, cShapes[axis-1]);
        } else if (!PtClCallBack || ptClStyle) {
          fxe->zoomAxis = 3;
          p_cursor(fxe->xe.win, P_SELECT);
          anchorX = oldX = x;
          anchorY = oldY = y;
          rubberBanding = PtClCallBack? ptClStyle : 1;
          DrawRubber(fxe, anchorX, anchorY);
        }
      } else if (!PtClCallBack) {
        /* no such thing as zoom or pan outside of all systems */
        fxe->zoomState = 0;
      }
    }

  } else {
    /* abort if 2nd button pressed */
    fxe->pressed = 0;
    ResetZoom(fxe);
  }
}

static void
ReleaseZoom(FXEngine *fxe, Drauing *drawing,
            int b,int md,int ix,int iy, unsigned long ms)
{
  int zoomState = fxe->zoomState;
  if (zoomState) {
    /* perform the indicated zoom operation */
    int iSystem = fxe->zoomSystem;
    GeSystem *system = GetSystemN(drawing, iSystem);
    GpXYMap *map = &fxe->xe.e.map;   /* NDC->VDC (x,y) mapping */
    GpReal x, xNDC = ((GpReal)ix - map->x.offset)/map->x.scale;
    GpReal y, yNDC = ((GpReal)iy - map->y.offset)/map->y.scale;

    /* get current pointer position (in world coordinates) */
    if (system &&
        /* be sure system has been scanned if limits extreme */
        (!(system->rescan || system->unscanned>=0) || !GdScan(system))) {
      FindCoordinates(system, xNDC, yNDC, &x, &y);
    } else {
      x = xNDC;
      y = yNDC;
      iSystem = ptClSystem = 0;
    }

    if (!PtClCallBack) {
      int axis = fxe->zoomAxis;
      GpReal factor = 1.0;

      if (!iSystem) {
        fxe->pressed = 0;
        ResetZoom(fxe);
        return;
      }

      if (zoomState==1) factor = 1.0/gxZoomFactor;
      else if (zoomState==3) factor = gxZoomFactor;

      /* the redraw triggered here can mess up a rubber band box/line */
      if (rubberBanding) {
        DrawRubber(fxe, anchorX, anchorY);
        rubberBanding = 0;
      }

      /* switch to current drawing and engine temporarily, then save
         limits if this is first mouse-driven zoom */
      GdSetDrawing(drawing);
      GpPreempt((Engine *)fxe);
      if (!(system->flags&D_ZOOMED)) GdSaveLimits(0);

      if (zoomState<4) {
        if (axis&1) DoZoom(factor, fxe->zoomX, x,
                           &system->trans.window.xmin,
                           &system->trans.window.xmax);
        if (axis&2) DoZoom(factor, fxe->zoomY, y,
                           &system->trans.window.ymin,
                           &system->trans.window.ymax);
      } else {
        if (axis&1) AltZoom(fxe->zoomX, x, &system->trans.window.xmin,
                            &system->trans.window.xmax);
        if (axis&2) AltZoom(fxe->zoomY, y, &system->trans.window.ymin,
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
      int modifier = 0;
      if (md&P_SHIFT) modifier|= 1;
      /* else if (state&LockMask) modifier|= 2; */
      else if (md&P_CONTROL) modifier|= 4;
      else if (md&P_META) modifier|= 8;
      else if (md&P_ALT) modifier|= 16;
      else if (md&P_COMPOSE) modifier|= 32;
      else if (md&P_KEYPAD) modifier|= 64;
      /* else if (state&Mod5Mask) modifier|= 128; */
      ptClCount--;
      PtClCallBack(&fxe->xe.e, iSystem, 1, x, y, modifier, xNDC, yNDC);
    }

    /* free the zoom/pan cursor and reset zoomState */
    fxe->pressed = 0;
    ResetZoom(fxe);
  } else if (fxe->pressed==2) {
    fxe->pressed = 0;
  }
}

static void
ResetZoom(FXEngine *fxe)
{
  int (*cback)(Engine *engine, int system,
               int release, GpReal x, GpReal y,
               int butmod, GpReal xn, GpReal yn) = PtClCallBack;
  /* free the zoom cursor and reset the zoom state */
  if (rubberBanding) {
    DrawRubber(fxe, anchorX, anchorY);
    rubberBanding = 0;
  }
  if (fxe->zoomState && fxe->xe.win)
    p_cursor(fxe->xe.win, P_CROSSHAIR);
  fxe->zoomState = 0;
  PtClCallBack = 0;
  if (cback) cback(0, -1, -1, 0., 0., -1, 0., 0.);
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

static int FindSystem(FXEngine *fxe, Drauing *drawing, int x, int y,
                      GeSystem **system, GpReal *xr, GpReal *yr)
{
  GeSystem *sys= drawing->systems, *thesys=sys;
  int nSystems= drawing->nSystems;
  GpXYMap *map= &fxe->xe.e.map;  /* NDC->VDC (x,y) mapping */
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
      if (!sys->elements) continue;
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

static void Find1System(FXEngine *fxe, Drauing *drawing, int iSystem,
                        int x, int y, GeSystem **system,
                        GpReal *xr, GpReal *yr)
{
  GpXYMap *map= &fxe->xe.e.map; /* NDC->VDC (x,y) mapping */
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

static GeSystem *GetSystemN(Drauing *drawing, int n)
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

static void
DrawRubber(FXEngine *fxe, int x, int y)
{
  p_win *w = fxe->xe.win;
  int iPass= 2;

  if (!w) return;

  p_color(w, P_XOR);
  p_pen(w, 1, P_SOLID);

  /* first undraw previous box or line, then draw new box or line */
  while (iPass--) {
    if (anchorX!=oldX || anchorY!=oldY) {
      if (rubberBanding==1) {
        /* this is a rubber box */
        p_rect(w, anchorX, anchorY, oldX, oldY, 1);
      } else {
        /* this is a rubber line */
        int x[2], y[2];
        x[0] = anchorX;  x[1] = oldX;
        y[0] = anchorY;  y[1] = oldY;
        p_i_pnts(w, x, y, 2);
        p_lines(w);
      }
    }
    oldX = x;
    oldY = y;
  }
}

/* ------------------------------------------------------------------------ */

int GxPointClick(Engine *engine, int style, int system,
                 int (*CallBack)(Engine *engine, int system,
                                 int release, GpReal x, GpReal y,
                                 int butmod, GpReal xn, GpReal yn))
{
  XEngine *xeng= GisXEngine(engine);
  if (!xeng || !xeng->w) return 1;

  /* set up state variables for point-and-click sequence */
  if (!(PtClCallBack= CallBack)) return 1;
  if (style==1 || style==2) ptClStyle= style;
  else ptClStyle= 0;
  if (system<0) ptClSystem= -1;
  else ptClSystem= system;
  ptClCount= 2;

  return 0;
}

/* ------------------------------------------------------------------------ */
