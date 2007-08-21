/*
 * pwin.c -- $Id$
 * X11 window management procedures
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"

#include "pstdlib.h"

#include <X11/Xutil.h>

/* during debug, generate as many Expose events as possible */
#ifndef P_DEBUG
# define BACKING_STORE WhenMapped
#else
# define BACKING_STORE NotUseful
#endif

#define PWIN_PLAIN  0
#define PWIN_MENU   1
#define PWIN_PIXMAP 2

static p_win *x_create(p_scr *s, Window parent, int hints, void *ctx,
                       int x, int y, int width, int height, int border,
                       p_col_t bg, int pwin_type);

static p_win *
x_create(p_scr *s, Window parent, int hints, void *ctx,
         int x, int y, int width, int height, int border,
         p_col_t bg, int pwin_type)
{
  /* these would be different for OpenGL windows -- not p_win* */
  int depth = CopyFromParent;
  Visual *v = CopyFromParent; /* (sic) X type consistency */

  x_display *xdpy = s->xdpy;
  Display *dpy = xdpy->dpy;
  XSetWindowAttributes cwa;
  unsigned long mask;
  int i;
  p_win *w = p_malloc(sizeof(p_win));
  if (!w) return 0;

  if (pwin_type==PWIN_MENU && !x_wire_events) goto fail0;

  w->s = s;
  w->context = ctx;

  if (pwin_type!=PWIN_PIXMAP) {
    w->parent = 0;
    cwa.background_pixel = s->colors[(bg>=255 || bg<242)? 0 : 255-bg].pixel;
    cwa.border_pixel = s->colors[1].pixel;
    cwa.backing_store = BACKING_STORE;
    mask = CWBackPixel | CWBorderPixel | CWBackingStore | CWCursor;
    if (pwin_type==PWIN_MENU) {
      cwa.override_redirect = cwa.save_under = True;
      mask |= CWOverrideRedirect | CWSaveUnder;
      cwa.cursor = x_cursor(s, P_W);
    } else {
      cwa.cursor = x_cursor(s, P_SELECT);
    }
    w->d = XCreateWindow(dpy, parent, x, y, width, height, border,
                         depth, InputOutput, v, mask, &cwa);
    if (w->d==None) goto fail0;

    if (p_hinsert(xdpy->id2pwin, P_IHASH(w->d), w)) goto fail1;

    if (x_wire_events) {
      /* note: if parent is completely occluded by its children, it will
       *       never actually receive any Expose events */
      mask = ExposureMask | StructureNotifyMask | FocusChangeMask |
        ButtonPressMask | ButtonReleaseMask | OwnerGrabButtonMask |
        EnterWindowMask | LeaveWindowMask;
      if (!(hints&P_NOKEY)) mask |= KeyPressMask;
      if (!(hints&P_NOMOTION)) mask |= PointerMotionMask;
      XSelectInput(dpy, w->d, mask);
    }

    /* load up empty palette colors */
    w->pixels = p_malloc(sizeof(p_col_t)*256);
    if (!w->pixels) goto fail1;
    for (i=0 ; i<242 ; i++) w->pixels[i] = s->colors[1].pixel;  /* fg */
    for (; i<256 ; i++) w->pixels[i] = s->colors[255-i].pixel;
    w->pixels[P_XOR] ^= s->colors[0].pixel;  /* bg */

  } else { /* offscreen pixmap */
    w->d = XCreatePixmap(dpy, parent, width, height, s->depth);
    w->pixels = 0;
  }
  w->n_palette = 0;
  w->rgb_pixels = 0;
  w->cmap = None;

  /* initialize window geometry (maintained by StructureNotify events) */
  w->x = x;
  w->y = y;
  w->width = width;
  w->height = height;
  w->xyclip[0] = w->xyclip[1] = w->xyclip[2] = w->xyclip[3] = 0;

  if (p_signalling) p_abort();
  return w;

  {
  fail1: XDestroyWindow(dpy, w->d);
  fail0: p_free(w);
  }
  if (p_signalling) p_abort();
  return 0;
}

p_win *
x_pwin(x_display *xdpy, Drawable d)
{
  return p_hfind(xdpy->id2pwin, P_IHASH(d));
}

void
p_winloc(p_win *w, int *x, int *y)
{
  *x = w->x;
  *y = w->y;
}

static XSizeHints *size_hints = 0;
static XWMHints *wm_hints = 0;
static XClassHint *class_hint = 0;

p_win *
p_window(p_scr *s, int width, int height, char *title,
         p_col_t bg, int hints, void *ctx)
{
  p_win *w = x_create(s, s->root, hints, ctx,
                      0, 0, width, height, 2, bg, PWIN_PLAIN);
  if (w) {
    XTextProperty x_title, *px_title;
    Display *dpy = s->xdpy->dpy;
    Window xwin = w->d;

    /* create and initialize private colormap if requested
     * - do it here instead of on demand (by p_palette), so that
     *   OpenGL child window can find the colormap if necessary
     * - GLX creation should call p_palette to install 5x9x5 RGB */
    if ((hints&P_PRIVMAP) && s->vclass==PseudoColor) {
      Visual *visual = DefaultVisual(dpy,s->scr_num);
      w->cmap = XCreateColormap(dpy, s->root, visual, AllocAll);

      if (w->cmap!=None) {
        XSetWindowAttributes attr;
        int i;

        /* copy entire default colormap so no flashing initially */
        Colormap cmap = DefaultColormap(dpy,s->scr_num);
        XColor map[256];
        int map_size = visual->map_entries;
        if (map_size>256) map_size = 256;
        for (i=0 ; i<map_size ; i++) map[i].pixel = i;
        XQueryColors(dpy, cmap, map, map_size);
        for (i=0 ; i<map_size ; i++) {
          map[i].flags = DoRed | DoGreen | DoBlue;
          XStoreColor(dpy, w->cmap, &map[i]);
        }

        /* attach colormap to window (window manager will install it) */
        attr.colormap = w->cmap;
        XChangeWindowAttributes(dpy, w->d, CWColormap, &attr);
      }
    }
    if (hints&P_RGBMODEL) x_rgb_palette(w);

    /* windows need to interact properly with window manager */

    if (s->xdpy->wm_delete!=None && x_wire_events)
      XSetWMProtocols(dpy, xwin, &s->xdpy->wm_delete, 1);

    if (!size_hints) size_hints = XAllocSizeHints();
    if (size_hints) {
      size_hints->x = size_hints->y = 0;
      size_hints->width = width;
      size_hints->height = height;
      size_hints->flags = PSize;
      if (hints&P_NORESIZE) {
        size_hints->min_width = size_hints->max_width = width;
        size_hints->min_height = size_hints->max_height = height;
        size_hints->flags |= PMinSize | PMaxSize;
      }
    }
    if (!wm_hints) wm_hints = XAllocWMHints();
    if (wm_hints) {
      wm_hints->initial_state = NormalState;
      wm_hints->input = (hints&P_NOKEY)? False : True;
      wm_hints->flags = StateHint | InputHint;
    }
    if (!class_hint) class_hint = XAllocClassHint();
    if (class_hint) {
      class_hint->res_name = 0;
      class_hint->res_class = (hints&P_DIALOG)? "GistDialog" : "Gist";
    }

    if (!title || !title[0] ||
        XStringListToTextProperty((char **)&title, 1, &x_title)==0)
      px_title = 0;
    else
      px_title = &x_title;

    XSetWMProperties(dpy, xwin, px_title, px_title, (char **)0, 0,
                     size_hints, wm_hints, class_hint);
    if (px_title) XFree((char *)x_title.value);

    /* no obvious "main window" to set transient property for dialogs
     * - use class to distinguish above */
    /* XSetTransientForHint(dpy, w->d, main_window); */

    if (!x_wire_events) XSelectInput(dpy, w->d, ExposureMask);

    XMapWindow(dpy, w->d);

    if (!x_wire_events) {
      /* support event-free operation when p_gui never called */
      XEvent event;
      XWindowEvent(dpy, w->d, ExposureMask, &event);
      XSelectInput(dpy, w->d, 0L);
      XSync(dpy, True);
    }
    if (p_signalling) p_abort();
  }

  return w;
}

p_win *
p_menu(p_scr *s, int width, int height, int x, int y,
       p_col_t bg, void *ctx)
{
  p_win *w = x_create(s, s->root, 0, ctx,
                      x, y, width, height, 0, bg, PWIN_MENU);
  if (w) {
    x_display *xdpy = s->xdpy;
    w->is_menu = 1;
    XMapWindow(xdpy->dpy, w->d);
    if (!xdpy->n_menus++ &&
        XGrabPointer(xdpy->dpy, w->d, True, ButtonPressMask |
                     ButtonReleaseMask | PointerMotionMask |
                     EnterWindowMask | LeaveWindowMask, GrabModeAsync,
                     GrabModeAsync, None, None, CurrentTime)!=GrabSuccess) {
      xdpy->n_menus = 0;
      w->is_menu = 0;  /* make this look like it is not a menu */
      p_destroy(w);
      w = 0;
    }
    if (p_signalling) p_abort();
  }
  return w;
}

p_win *
p_offscreen(p_win *parent, int width, int height)
{
  p_win *w = x_create(parent->s, parent->d, 0, 0,
                      0, 0, width, height, 0, P_BG, PWIN_PIXMAP);
  if (w) {
    w->is_menu = 0;
    w->parent = parent;
    p_clear(w);  /* otherwise initial contents undefined */
  }
  return w;
}

void
p_bitblt(p_win *w, int x, int y, p_win *offscreen,
         int x0, int y0, int x1, int y1)
{
  if (w && w==offscreen->parent) {
    p_scr *s = w->s;
    GC gc = x_getgc(s, w, FillSolid);
    XCopyArea(s->xdpy->dpy, offscreen->d, w->d, gc,
              x0, y0, x1-x0, y1-y0, x, y);
    if (p_signalling) p_abort();
  }
}

p_win *
p_metafile(p_win *parent, char *filename,
           int x0, int y0, int width, int height, int hints)
{
  return 0;
}

void
p_destroy(p_win *w)
{
  x_display *xdpy = w->s->xdpy;
  Display *dpy = xdpy->dpy;
  int connected = dpy && !xdpy->panic;

  if (w->is_menu && connected && !(--xdpy->n_menus))
    XUngrabPointer(dpy, CurrentTime);

  if (!w->parent && connected) {
    p_col_t *rgb_pixels = w->rgb_pixels;
    if (w->cmap) x_cmzap(dpy, &w->cmap);
    else p_palette(w, (p_col_t *)0, 0);
    if (rgb_pixels) {
      x_tmpzap(&w->pixels);
      w->rgb_pixels = 0;
      w->pixels = rgb_pixels;
      w->n_palette = 225;
      p_palette(w, (p_col_t *)0, 0);
    }
  }
  x_tmpzap(&w->pixels);
  x_tmpzap(&w->rgb_pixels);

  /* XDestroyWindow will clear the selection */
  if (xdpy->sel_owner==w) xdpy->sel_owner = 0;

  if (connected) {
    Window xwin = w->d;
    if (!w->parent) {
      /* make sure no more events will be delivered */
      p_hinsert(xdpy->id2pwin, P_IHASH(xwin), (void *)0);
      w->d = None;
      XDestroyWindow(dpy, xwin);
    } else {
      w->d = None;
      XFreePixmap(dpy, xwin);
    }
  }

  p_free(w);
}

void
p_resize(p_win *w, int width, int height)
{
  if (!w->parent) XResizeWindow(w->s->xdpy->dpy, w->d, width, height);
  if (p_signalling) p_abort();
}

void
p_raise(p_win *w)
{
  if (!w->parent) {
    Display *dpy = w->s->xdpy->dpy;
    Window xwin = w->d;
    XMapWindow(dpy, xwin);   /* user may have iconified it */
    XRaiseWindow(dpy, xwin);
    if (p_signalling) p_abort();
  }
}

void
p_clear(p_win *w)
{
  Display *dpy = w->s->xdpy->dpy;
  if (!w->parent) {
    XClearWindow(dpy, w->d);
  } else {
    GC gc = x_getgc(w->s, w, FillSolid);
    p_color(w, P_BG);
    XFillRectangle(dpy, w->d, gc, 0,0, w->width+1, w->height+1);
  }
  if (p_signalling) p_abort();
}

void
p_flush(p_win *w)
{
  XFlush(w->s->xdpy->dpy);
  if (p_signalling) p_abort();
}
