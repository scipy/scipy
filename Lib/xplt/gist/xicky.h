/*
 * XICKY.H
 *
 * $Id$
 *
 * Declare icky X windows utilities for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef XICKY_H
#define XICKY_H

#include <X11/Xlib.h>
/* XVisualInfo declared in Xutil.h */
#include <X11/Xutil.h>

typedef struct GxFontProps GxFontProps;
struct GxFontProps {
#define FONT_FAMILIES 5
#define FONT_SIZES 6
#define FONT_FALLBACKS 3
  char *foundry;    /* "Adobe" or 0 if none exist */
  char *facename;   /* courier, times, helvetica, symbol,
		       or new century schoolbook */
  char *slantname;  /* italic or oblique (or 0) */
  int masks[FONT_SIZES];     /* for point sizes 8, 10, 12, 14, 18, 24
				according to following-- */

#define HAS_NORMAL 1
#define HAS_SLANT 2
#define HAS_BOLD 4
#define HAS_BOLD_SLANT 8
};

typedef struct GxDisplay GxDisplay;
typedef struct GxScreen GxScreen;

struct GxDisplay {
  GxDisplay *next;  /* All X server connections kept in a list */
  int references;    /* number of GxScreen pointers passed out - 1 */

  Display *display;
  char *normalizedName;

  /* List of screens and their defaults.  */
  int nScreens;
  GxScreen *screens;

  /* List returned by XGetVisualInfo.  */
  int nVisuals;
  XVisualInfo *v;

  /* fonts for this display */
  GxFontProps fontProps75[FONT_FAMILIES], fontProps100[FONT_FAMILIES];
  int fontFallbacks;  /* FONT_FALLBACKS bit masks, set if present */
  XFontStruct *loadedFonts[5];
  int loadedID[5];    /* font id is gist font<<(X_SIZE_BITS+2)
			 ored with masks below */
  XFontStruct *defaultFont, *permFont;

#define X_SIZE_BITS 4
#define X_SIZE_MASK 0xf
#define IS_FALLBACK 1
#define IS_100DPI 2
#define PERM_FONT_ID -1
};

/* DirectColor visual class requires more information than the
   X data structure provides.  In order to construct colors, the location
   of the red_mask, green_mask, and blue_mask are required.  Furthermore,
   colormap_size reflects the largest of the three colormaps, while
   the maximum number of colors available for XAllocColorCells is the
   length of the smallest of the three color tables.
   Unless the server enforces some sort of policy (unstated?) for how
   it allocates new colors -- such as "all read/write cells shall have
   equal red, green, and blue indices" -- DirectColor could get very
   messy.  */
typedef struct DirectColorInfo DirectColorInfo;
struct DirectColorInfo {
  int rwMapSize;      /* 2^min(red bits, green bits, blue bits) */
  int red_shift, green_shift, blue_shift;  /* location of masks */
  int red_size, green_size, blue_size;    /* 2^(red_bits), etc. */
  /* red_mask, green_mask, and blue_mask are in XVisualInfo structure */
};

struct GxScreen {
  GxDisplay *owner;

  /* Display, maxRequest are actually server-wide properties belonging
     to the owner-- put here for convenience.  */
  Display *display;
  long maxRequest;

  /* Any given server may support several screens.  */
  int screen;
  Window root;
  int width, height;

  /* Each screen may support several visuals; v points to the default
     visual on this screen.  You can find others by examining the other
     elements in the owner->v list.  */
  XVisualInfo *v;  /* members are:
		      Visual *visual;
		      VisualID visualid;
		      int screen;
		      unsigned int depth;
		      int class;
		      unsigned long red_mask, green_mask, blue_mask;
		      int colormap_size;   (same as visual->map_entries)
		      int bits_per_rgb;
		    */
  DirectColorInfo dcinfo;

  /* default colormap and pixels for this screen  */
  XColor stdColors[10];
    /* the 10 are:
       bg, fg, black, white, red, green, blue, cyan, magenta, yellow */
  Colormap colormap;

  /* rotated text requires two pixmaps and a GC */
  Pixmap rot_normal, rot_rotated;
  GC rot_gc;
};

extern GxDisplay *gistX;

/* Note-- If the display is already open, it is safe to call GxConnect
   with DisplayString(display).  */
extern GxScreen *GxConnect(char *displayName);
/* GxDisconnect returns 1 if the server connection was actually closed,
   0 if there are still open connections */
extern int GxDisconnect(GxScreen *xscr);

extern void GxSetProperties(char *name, Display *display, Window window,
			    unsigned int winx, unsigned int winy);
extern void GxInitialize(int *argc, char **argv);


#ifndef GIST_H
/* definition of GpColorCell from gist.h */
typedef struct GpColorCell GpColorCell;
struct GpColorCell {
  unsigned char red, green, blue, gray;
};
#endif

/* Return a pixelMap[256] (or use the space supplied).
   The first nColors entries will be filled with the "best" shared
   colors equal to or approximating the given palette; the remaining
   entries will be filled with the pixel corresponding to FG_COLOR.  */
extern unsigned long *GxShareColors(GxScreen *xscr, GpColorCell *palette,
				    int nColors, unsigned long *pixelMap);

/* Return a pixelMap[256] (or use the space supplied).
   The first nColors entries will be filled with read/write pixel
   values corresponding to the specified palette; the remaining
   entries will be filled with the pixel corresponding to FG_COLOR.
   If fewer than nColors cells were available in the default colormap,
   a private colormap is returned.  Conversely, if *private!=0 on
   input, but nColors were available in the default map, then
   *private is freed and returned as 0.
   (Note-- Special coding prevents *private==0 for a valid colormap.)  */
extern unsigned long *GxExactColors(GxScreen *xscr, GpColorCell *palette,
				    int nColors, unsigned long *pixelMap,
				    Colormap *private, int oldNColors);

/* Return a list of all colors in the default colormap which can be
   shared (worker for GxShareColors).  */
extern int GxGetSharable(GxScreen *xscr, XColor **shared, int *nShared);

/* Free the XColor* returned by GxGetSharable.  */
extern void GxFreeSharable(GxScreen *xscr, XColor *shared, int nShared);

#endif
