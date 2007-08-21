/*
 * playx.h -- $Id$
 * declare routines and structs to use an X11 server
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

/* X11 implementation files include this instead of play.h */
#include "play.h"
#include "phash.h"
#include <X11/Xlib.h>

#define N_FONT_CACHE 6

/* the Display struct may be shared among several root windows,
 * in the unusual case that a single server drives several screens
 * - fonts, however, are shared among all screens of a server
 *   (I ignore a possible exception of the default font, which could
 *    in principle differ from screen to screen on a server -- this
 *    may be a severe problem if the two screens differ greatly in
 *    resolution, so that what is legible on one is not on the other
 *    -- hopefully this is even rarer than multiple screens) */

typedef struct x_display x_display;
struct x_display {
  int panic;
  p_scr *screens;       /* list of screens on this server (for panic) */
  x_display *next;      /* list of all servers */
  Display *dpy;

  Atom wm_protocols, wm_delete;
  p_hashtab *id2pwin;   /* use phash instead of XContext */

  XFontStruct *font;    /* default font to use on this server */
  int unload_font;      /* non-0 if font must be unloaded */

  struct {
    XFontStruct *f;
    int font, pixsize, next;
  } cached[N_FONT_CACHE];
  int most_recent;

  struct {
    int nsizes, *sizes;
    char **names;
  } available[20];

  Cursor cursors[14];

  /* number of motion events queued when previous motion callback
   * completed -- all but the last will be skipped */
  int motion_q;

  unsigned int meta_state, alt_state;  /* masks for modifier keys */

  /* selection data */
  p_win *sel_owner;
  char *sel_string;

  /* menu count for pointer grabs */
  int n_menus;
};

extern x_display *x_displays;
typedef struct x_cshared x_cshared;

struct p_scr {
  x_display *xdpy;  /* may be accessing multiple screens on server */
  p_scr *next;       /* keep list of all screens on this server */

  int scr_num;          /* screen number on this server */
  Window root;          /* root window on this screen */
  int width, height, depth;  /* of root window */

  int vclass;           /* visual class for this screen */
  /* pixels==0 (part of p_win) for PseudoColor visual
   * for all others, points to pixels[256] */
  p_col_t *pixels;
  /* red, green, and blue masks for TrueColor and DirectColor */
  p_col_t rmask, gmask, bmask;
  Colormap cmap;        /* None except for GrayScale and DirectColor */

  XColor colors[14];    /* standard colors */
  int free_colors;      /* bit flags for which ones need to be freed */
  Pixmap gray;          /* in case GRAYA-GRAYD are stipples */
  int gui_flags;        /* marker for stippled grays */
  x_cshared *shared;    /* tables for shared PseudoColor colors */

  /* generic graphics context and its current state */
  GC gc;
  p_col_t gc_color;
  int gc_fillstyle;     /* invalid for stippled colors, see colors.c */
  p_win *gc_w_clip;     /* gc contains clipping for this window */
  int gc_width, gc_type;
  int gc_font, gc_pixsize;

  /* temporaries required for rotating fonts (see textout.c) */
  void *tmp;
  XImage *image;
  int own_image_data;
  Pixmap pixmap;
  GC rotgc;
  int rotgc_font, rotgc_pixsize, rotgc_orient;
};

struct p_win {
  void *context;     /* application context for event callbacks */
  p_scr *s;

  Drawable d;
  p_win *parent;     /* non-0 only for offscreen pixmaps */
  int is_menu;       /* non-0 only for menus */

  Colormap cmap;
  p_col_t *pixels, *rgb_pixels;
  int n_palette;  /* number of pixels[] belonging to palette */
  int x, y, width, height, xyclip[4];
};

/* point list for polylines, fills, dots, segments */
extern XPoint x_pt_list[2050];
extern int x_pt_count;

/* retrieve p_win* given Window id number, Display* (for event handling)
 * - the p_win context can be used to back up the hierarchy further
 * - this could be implemented using XContext mechanism */
extern x_display *x_dpy(Display *dpy);
extern p_win *x_pwin(x_display *xdpy, Drawable d);

/* ordinary and I/O X error handlers */
extern int x_err_handler(Display *dpy, XErrorEvent *event);
extern int x_panic(Display *dpy);
extern void (*x_on_panic)(p_scr *s);

/* arrange to deliver X events for an X window to the event
 * handler for the corresponding p_win
 * this is virtual to allow a simpler mode in case p_gui is never called */
extern void (*x_wire_events)(x_display *xdpy, int disconnect);

/* routines to convert Gist colors, linestyles, and fonts to X11 */
extern int x_rgb_palette(p_win *w);
extern XFontStruct *x_font(x_display *xdpy, int font, int pixsize);
extern void x_clip(Display *dpy, GC gc, int x0, int y0, int x1, int y1);
extern GC x_getgc(p_scr *s, p_win *w, int fillstyle);
extern p_col_t x_getpixel(p_win *w, p_col_t color);
extern void x_nuke_shared(p_scr *s);

/* optional X resource values (class Gist) */
extern char *x_xfont;       /* boldfont, font, Font */
extern char *x_foreground;  /* foreground, Foreground */
extern char *x_background;  /* background, Background */
extern char *x_guibg;       /* guibg */
extern char *x_guifg;       /* guifg */
extern char *x_guihi;       /* guihi */
extern char *x_guilo;       /* guilo */

extern Cursor x_cursor(p_scr *s, int cursor);

/* simple destructors that zero input pointer */
extern void x_rotzap(p_scr *s);
extern void x_tmpzap(void *ptmp);   /* x_tmpzap(anytype **ptmp) */
extern void x_gczap(Display *dpy, GC *pgc);
extern void x_imzap(p_scr *s);
extern void x_pxzap(Display *dpy, Pixmap *ppx);
extern void x_cmzap(Display *dpy, Colormap *pcm);
