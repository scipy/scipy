/*
 * XICKY.C
 *
 * $Id$
 *
 * Implement icky X windows utilities for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "xicky.h"
#include "xfont.h"

/* Get BSD-style (index instead of strchr) string functions via Xos.h... */
#include <X11/Xos.h>

#ifdef STDC_HEADERS
#include <string.h>
#else
#ifndef SIZE_T_TYPE
#define SIZE_T_TYPE unsigned long
#endif
extern int strcmp(const char *, const char *);
extern int strncmp(const char *, const char *, SIZE_T_TYPE);
extern SIZE_T_TYPE strlen(const char *);
extern char *strcpy(char *, const char *);
extern char *strncpy(char *, const char *, SIZE_T_TYPE);
extern char *strncat(char *, const char *, SIZE_T_TYPE);
#endif

extern void *(*GmMalloc)(long);
extern void (*GmFree)(void *);
extern long strtol(const char *, char **, int);

extern char gistError[128];  /* most recent error message */

/* _XGetHostname is in MIT internal library.  If this external is
   unsatisfied, the code is reproduced in host.c.  */
extern int _XGetHostname(char *, int);

#define SQ(x) ((x)*(x))

static char *NormalizeDisplay(char *displayName);

extern void GxDirectColor(const XVisualInfo *vinfo, DirectColorInfo *dcinfo);

/* ------------------------------------------------------------------------ */

GxDisplay *gistX= 0;

/* ------------------------------------------------------------------------ */

#define MAX_DPY_NAME 256

static char *NormalizeDisplay(char *displayName)
{
  /* displayName may be NULL or host:server.screen, where

     host is either the hostname OR
     (1) local ("best" connection to local host)
     (2) unix  (unix socket connection to local host)
     (3) (nil) (same as unix)

     NormalizeDisplay returns a copy of displayName up to but not
     including the period before the screen number, UNLESS displayName
     is NULL or one of the 3 special cases.  If displayName is NULL,
     NormalizeDisplay calls XDisplayName (which retrieves the DISPLAY
     environment variable), then procedes as before.  If one of the 3
     special cases is encountered, NormalizeDisplay calls the
     internal library routine _XGetHostname (supplied in host.c if
     unavailable on your X implementation) and substitutes that for the
     host field in the display name.

     The resulting string should be freed using GmFree.

     The point of all this is to be able to recognize by name a server
     that we already have made a connection to.
   */

  char *dName= XDisplayName(displayName);
  char *colon= dName? index(dName, ':') : dName;
  char *dot= colon? index(colon, '.') : 0;
  char dpyName[MAX_DPY_NAME+1], *normal;
  int hlen, slen;

  if (colon==dName || strncmp(dName, "unix:", 5)==0 ||
                      strncmp(dName, "local:", 6)==0) {
    dName= dpyName;
    hlen= _XGetHostname(dName, MAX_DPY_NAME);
  } else {
    hlen= colon? colon-dName : 0;
  }
  slen= dot? dot-colon : (colon? strlen(colon) : 0);

  normal= (char *)GmMalloc(hlen+slen+1);
  if (!normal) return 0;
  if (hlen) strncpy(normal, dName, hlen);
  if (slen) strncpy(normal+hlen, colon, slen);
  normal[hlen+slen]= '\0';

  return normal;
}

/* ------------------------------------------------------------------------ */

unsigned long *GxShareColors(GxScreen *xscr, GpColorCell *palette,
			     int nColors, unsigned long *pixelMap)
{
  int class= xscr->v->class;
  int gray= (class==StaticGray || class==GrayScale);
  int i, irev, iinc, nGot, nShared;
  XColor color, *shared= 0;

  /* get a pixelMap if none supplied */
  if (!pixelMap) {
    /* always grab enough space for a maximum length 256 element map */
    pixelMap= (unsigned long *)GmMalloc(sizeof(unsigned long)*256);
    if (!pixelMap) {
      strcpy(gistError, "memory manager failed in GxShareColors");
      return 0;
    }
  }

  /* Allocate colors in bit reversed order to give approximately
     equal coverage of the palette if the allocation process
     fails to get every requested color.  */
  if (nColors>1) {
    iinc= 1;
    for (i=nColors-1 ; i<128 ; i<<= 1) iinc<<= 1;
  } else {
    iinc= 256;
  }
  nGot= -1;
  for (i=0 ; i<256 ; i+=iinc) {
    /* bit reverse i */
    irev= (i>>7) | ((i>>5)&2) | ((i>>3)&4) | ((i>>1)&8) |
      ((i<<1)&16) | ((i<<3)&32) | ((i<<5)&64) | ((i<<7)&128);
    if (irev>=nColors) continue;

    /* If nGot has not been set, XAllocColor has not yet failed, so
       go ahead and keep trying.  */
    if (nGot<0) {
      if (gray) {
	color.red= color.green= color.blue= palette[irev].gray<<8;
      } else {
	color.red= palette[irev].red<<8;
	color.green= palette[irev].green<<8;
	color.blue= palette[irev].blue<<8;
      }
      /* Note: each entry in pixelMap always corresponds to one successful
	 call to XAllocColor */
      if (XAllocColor(xscr->display, xscr->colormap, &color))
	pixelMap[irev]= color.pixel;
      else
	nGot= i;
    }

    /* If more than half of the requested colors have been exactly
       allocated, just fill in the remaining ones by their nearest
       neighbor.  This will give poor results if the palette contains
       wildly skipping colors (as, for example, a way to get contour
       lines on an image), but it is much faster than the "closest
       fit" algorithm.  Also, by placing outlier colors at even
       palette addresses, you can guarantee they will already have
       been satisfied.  */
    if (2*nGot>=nColors) {
      /* The color immediately below this color has already been
	 allocated.  Use it again.  The repeated call to XAllocColor
	 allows simpler logic for the eventual call to XFreeColors--
	 O'Reilly Volume 0 protocol request FreeColors
         states that the pixel must be specified to XFreeColors
         the same number of times it was specified to XAllocColor.  */
      pixelMap[irev]= pixelMap[irev-1];
      if (gray) {
	color.red= color.green= color.blue= palette[irev-1].gray<<8;
      } else {
	color.red= palette[irev-1].red<<8;
	color.green= palette[irev-1].green<<8;
	color.blue= palette[irev-1].blue<<8;
      }
      XAllocColor(xscr->display, xscr->colormap, &color);

    /* If less than half of the requested table could be allocated
       as shared colors, try to get the sharable color cells closest
       to the requested colors.  This is a time consuming operation,
       and can quite easily give terrible results, but we're desperate
       to do something sensible.  */
    } else if (nGot>=0) {
      long dc, mindc;
      int j, jmin;
      if (!shared) {
	if (GxGetSharable(xscr, &shared, &nShared)) {
	  GmFree(pixelMap);
	  return 0;   /* memory manager failure */
	}
      }
      if (gray) {
	color.red= color.green= color.blue= palette[irev].gray<<8;
      } else {
	color.red= palette[irev].red<<8;
	color.green= palette[irev].green<<8;
	color.blue= palette[irev].blue<<8;
      }
      mindc= SQ(((long)shared[0].red-(long)color.red)>>8)+
	SQ(((long)shared[0].green-(long)color.green)>>8)+
	  SQ(((long)shared[0].blue-(long)color.blue)>>8);
      jmin= 0;
      for (j=1 ; j<nShared ; j++) {
	dc= SQ(((long)shared[j].red-(long)color.red)>>8)+
	  SQ(((long)shared[j].green-(long)color.green)>>8)+
	    SQ(((long)shared[j].blue-(long)color.blue)>>8);
	if (dc<mindc) {
	  mindc= dc;
	  jmin= j;
	}
      }
      pixelMap[irev]= shared[jmin].pixel;
      XAllocColor(xscr->display, xscr->colormap, &shared[jmin]);
    }
  }

  /* clean up after closest sharable color calculation */
  if (shared) GxFreeSharable(xscr, shared, nShared);

  /* fill remainder of pixelMap with foreground color (for DrawImage) */
  for (i=nColors ; i<256 ; i++)
    pixelMap[i]= xscr->stdColors[1].pixel;

  return pixelMap;
}

int GxGetSharable(GxScreen *xscr, XColor **shared, int *nShared)
{
  int class= xscr->v->class;
  int mapSize= class==DirectColor? xscr->dcinfo.rwMapSize :
                                   xscr->v->colormap_size;
  XColor *s;
  int i, n;
  unsigned long rmask, gmask, bmask;
  int rshft, gshft, bshft;

  if (mapSize>256) mapSize= 256;
  s= (XColor *)GmMalloc(sizeof(XColor)*mapSize);
  if (!s) {
    strcpy(gistError, "memory manager failed in GxGetSharable");
    *shared= 0;
    *nShared= 0;
    return 1;
  }

  /* query to find every color in default colormap */
  if (class==DirectColor) {
    rmask= xscr->v->red_mask;
    gmask= xscr->v->green_mask;
    bmask= xscr->v->blue_mask;
    rshft= xscr->dcinfo.red_shift;
    gshft= xscr->dcinfo.green_shift;
    bshft= xscr->dcinfo.blue_shift;
  } else {
    /* prevent compiler from complaining about uninitialized variables */
    rmask= gmask= bmask= 0;
    rshft= gshft= bshft= 0;
  }
  for (i=0 ; i<mapSize ; i++) {
    if (class!=DirectColor) {
      s[i].pixel= i;
    } else {
      /* NB: After 20 readings of all of the X documentation I can get
	 my hands on, I'm still not sure this is guaranteed to produce
	 a list of all possible pixel values with XQueryColors...  */
      s[i].pixel= ((i<<rshft)&rmask) | ((i<<gshft)&gmask) |
	((i<<bshft)&bmask);
    }
  }
  XQueryColors(xscr->display, xscr->colormap, s, mapSize);

  /* Second pass tries to allocate each of these colors to see
     whether it is sharable -- this may produce duplicates if some
     privately owned color happens to coincide with a shared color.
     Hopefully, including such a duplicate multiple times in the
     list passed to XFreeColors will cancel the effect of the
     multiple calls to XAllocColor...  The same thing could happen
     with unallocated colors, but GxGetSharable is intended to be
     called after all unallocated colors have been exhausted.  */
  n= 0;
  for (i=0 ; i<mapSize ; i++) {
    if (XAllocColor(xscr->display, xscr->colormap, &s[i])) {
      /* this color is sharable */
      if (n<i) s[n]= s[i];   /* compress list */
      n++;
    } /* else this color must be privately owned */
  }

  *shared= s;
  *nShared= n;
  return 0;
}

void GxFreeSharable(GxScreen *xscr, XColor *shared, int nShared)
{
  unsigned long *pixels= &shared->pixel;
  int i;
  if (!shared) return;
  for (i=0 ; i<nShared ; i++) pixels[i]= shared[i].pixel;
  XFreeColors(xscr->display, xscr->colormap, pixels, nShared, 0L);
  GmFree(shared);
}

unsigned long *GxExactColors(GxScreen *xscr, GpColorCell *palette,
			     int nColors, unsigned long *pixelMap,
			     Colormap *private, int oldNColors)
{
  int class= xscr->v->class;
  int i;
  XColor transfer;

  /* GxShareColors should be perfectly suitable for immutable colormaps */
  if (class==StaticGray || class==StaticColor || class==TrueColor)
    return GxShareColors(xscr, palette, nColors, pixelMap);

  /* get a pixelMap if none supplied */
  if (!pixelMap) {
    /* always grab enough space for a maximum length 256 element map */
    pixelMap= (unsigned long *)GmMalloc(sizeof(unsigned long)*256);
    if (!pixelMap) {
      strcpy(gistError, "memory manager failed in GxExactColors");
      return 0;
    }
  }

  if (nColors<=0) {
    if (*private) {
      XFreeColormap(xscr->display, *private);
      *private= 0;
    }
    for (i=nColors ; i<256 ; i++)
      pixelMap[i]= xscr->stdColors[1].pixel;
    return pixelMap;
  }

  if ((*private && oldNColors==nColors) ||
      XAllocColorCells(xscr->display, xscr->colormap, False, (void *)0,
		       0, pixelMap, nColors)) {
    /* the default colormap is big enough to cover everything, or the
       current private colormap is just right */
    Colormap xpalette;
    if (*private) {
      if (oldNColors!=nColors) {
	XFreeColormap(xscr->display, *private);
	*private= 0;
	xpalette= xscr->colormap;
      } else {
	xpalette= *private;
      }
    } else {
      xpalette= xscr->colormap;
    }
    for (i=0 ; i<nColors ; i++) {
      transfer.pixel= pixelMap[i];
      transfer.red= palette[i].red<<8;
      transfer.green= palette[i].green<<8;
      transfer.blue= palette[i].blue<<8;
      transfer.flags= DoRed | DoGreen | DoBlue;
      XStoreColor(xscr->display, xpalette, &transfer);
    }

  } else {
    /* must create (or reuse) a private colormap */
    int mapSize= class==DirectColor? xscr->dcinfo.rwMapSize :
                                     xscr->v->colormap_size;
    int j, k, nstd, skip;
    unsigned long pixel;
    unsigned long rmask, gmask, bmask;
    int rshft, gshft, bshft;
    Colormap xpalette;

    if (!*private) {
      *private= xpalette= XCreateColormap(xscr->display, xscr->root,
					  xscr->v->visual, AllocAll);
      /* I can't find any value of a Colormap (an X resource ID) which
	 is guaranteed invalid, so here I guarantee that xpalette is
	 non-zero.  In the unlikely event that xpalette==0, just ask
	 for another colormap.  */
      if (!xpalette) {
	 xpalette= XCreateColormap(xscr->display, xscr->root,
				   xscr->v->visual, AllocAll);
	 XFreeColormap(xscr->display, *private);
	 *private= xpalette;
       }
    } else {
      xpalette= *private;
    }

    /* Allocate black pixel, white pixel, and as many stdColors as
       possible to prevent flickering when this colormap is installed.  */
    nstd= mapSize-nColors<10? mapSize-nColors : 10;
    if (class==DirectColor) {
      rmask= xscr->v->red_mask;
      gmask= xscr->v->green_mask;
      bmask= xscr->v->blue_mask;
      rshft= xscr->dcinfo.red_shift;
      gshft= xscr->dcinfo.green_shift;
      bshft= xscr->dcinfo.blue_shift;
    } else {
      /* prevent compiler from complaining about uninitialized variables */
      rmask= gmask= bmask= 0;
      rshft= gshft= bshft= 0;
    }
    /* Loop downwards on the assumption that the server will pass out
       colors to other clients by looping upwards.  */
    k= nColors-1;
    for (i=mapSize-1 ; i>=0 ; i--) {
      if (class!=DirectColor) {
	pixel= i;
      } else {
	/* This is still not quite correct, since the standard colors
	   may overlap allocated colors for one of r g or b...
	   Must rely on server to not give away read/write colorcells
	   with different r, g, b indices, and on standard colors
	   having a lot of overlap with one another (in principle,
	   all 8 standard colors require only values found in BlackPixel
	   and WhitePixel).  */
	pixel= ((i<<rshft)&rmask) | ((i<<gshft)&gmask) |
	  ((i<<bshft)&bmask);
      }
      if (k>=0) {
	if (i > k) {
	  /* skip this if it matches a stdColor */
	  for (j=0 ; j<nstd ; j++)
	    if (pixel==xscr->stdColors[j].pixel) break;
	  skip= (j<nstd);
	} else skip= 0;
      } else {
	skip= 1;
      }
      /* Either set the pixel value to the next requested color, or
         to the matching standard color.  */
      if (!skip) {
	pixelMap[k]= transfer.pixel= pixel;
	transfer.red= palette[k].red<<8;
	transfer.green= palette[k].green<<8;
	transfer.blue= palette[k].blue<<8;
	transfer.flags= DoRed | DoGreen | DoBlue;
	k--;
      } else {
	transfer.pixel= pixel;
	XQueryColor(xscr->display, xscr->colormap, &transfer);
      }
      XStoreColor(xscr->display, xpalette, &transfer);
    }
  }

  /* fill remainder of pixelMap with foreground color (for DrawImage) */
  for (i=nColors ; i<256 ; i++)
    pixelMap[i]= xscr->stdColors[1].pixel;

  return pixelMap;
}

void GxDirectColor(const XVisualInfo *vinfo, DirectColorInfo *dcinfo)
{
  int i, size;

  for (i=0 ; i<sizeof(unsigned long) ; i++)
    if (vinfo->red_mask & (1<<i)) break;
  dcinfo->red_shift= i;
  size= 2;
  for (i++ ; i<sizeof(unsigned long) ; i++) {
    if (vinfo->red_mask & (1<<i)) size*= 2;
    else break;
  }
  dcinfo->red_size= size;
  dcinfo->rwMapSize= size;

  for (i=0 ; i<sizeof(unsigned long) ; i++)
    if (vinfo->green_mask & (1<<i)) break;
  dcinfo->green_shift= i;
  size= 2;
  for (i++ ; i<sizeof(unsigned long) ; i++) {
    if (vinfo->green_mask & (1<<i)) size*= 2;
    else break;
  }
  dcinfo->green_size= size;
  if (size<dcinfo->rwMapSize) dcinfo->rwMapSize= size;

  for (i=0 ; i<sizeof(unsigned long) ; i++)
    if (vinfo->blue_mask & (1<<i)) break;
  dcinfo->blue_shift= i;
  size= 2;
  for (i++ ; i<sizeof(unsigned long) ; i++) {
    if (vinfo->blue_mask & (1<<i)) size*= 2;
    else break;
  }
  dcinfo->blue_size= size;
  if (size<dcinfo->rwMapSize) dcinfo->rwMapSize= size;
}

/* ------------------------------------------------------------------------ */

extern void GxUnlink(GxDisplay *owner);
void GxUnlink(GxDisplay *owner)
{
  if (gistX==owner) {
    gistX= owner->next;
  } else {
    GxDisplay *xdpy, *prev= gistX;
    for (xdpy=prev->next ; xdpy ; xdpy=xdpy->next) {
      if (xdpy==owner) {
	prev->next= xdpy->next;
	break;
      }
      prev= xdpy;
    }
  }
}

GxScreen *GxConnect(char *displayName)
{
  char *normalizedName= NormalizeDisplay(displayName);
  GxDisplay *xdpy= gistX;
  Display *display;
  int nVisuals, nScreens, i, j;
  GxScreen *screens= 0;
  XVisualInfo template, *vinfo= 0;
  char *fgOption, *bgOption, *fontOption;
  Colormap colormap;
  XColor color;

  /* Is there already a connection to this server?  */
  while (xdpy) {
    if (strcmp(xdpy->normalizedName, normalizedName)==0) break;
    xdpy= xdpy->next;
  }

  if (xdpy) {
    /* A connection to this server has already been made.  */
    int screen= 0;    /* default screen number is 0 */
    GmFree(normalizedName);
    if (displayName) {
      char *scr= index(displayName, ':');
      if (scr) scr= index(scr, '.');
      if (scr) {
	char *p= scr;
	for (p++ ; (*p)>='0' && (*p)<='9' ; p++);
	if (*p) screen= xdpy->nScreens; /* will cause error below */
	else if (p>scr) screen= (int)strtol(scr+1, (char **)0, 10);
      }
    }
    if (screen<xdpy->nScreens) {
      /* Increment reference counter and return screen.  */
      xdpy->references++;
      return xdpy->screens+screen;
    } else {
      return 0;
    }
  } /* (never drops out of this if) */

  /* No connection to this server yet exists.  */
  if (!(display= XOpenDisplay(displayName)) ||
      (nScreens= ScreenCount(display))<=0 ||
      !(screens= (GxScreen *)GmMalloc(sizeof(GxScreen)*nScreens)) ||
      !(vinfo= XGetVisualInfo(display, VisualNoMask, &template, &nVisuals)) ||
      !(xdpy= (GxDisplay *)GmMalloc(sizeof(GxDisplay)))) {
    if (vinfo) XFree((char *)vinfo);
    if (screens) GmFree(screens);
    if (display) XCloseDisplay(display);
    GmFree(normalizedName);
    strcpy(gistError, "failed to connect to X server ");
    strncat(gistError, XDisplayName(displayName), 40);
    return 0;
  }

  /* Fill in GxDisplay structure.  */
  xdpy->references= 0;
  xdpy->display= display;
  xdpy->normalizedName= normalizedName;
  xdpy->nScreens= nScreens;
  xdpy->screens= screens;
  xdpy->nVisuals= nVisuals;
  xdpy->v= vinfo;

  /* Check for user preferred foreground and background colors and font
     on this display.
     NB- XGetDefault uses the class name "Name" for all components
     except the program (first), which gets "Program".  Since these
     are unlikely to match anything interesting, both the specific and
     class names for foreground and background must be tried.  */
  fgOption= XGetDefault(display, "Gist", "foreground");
  if (!fgOption) fgOption= XGetDefault(display, "Gist", "Foreground");
  bgOption= XGetDefault(display, "Gist", "background");
  if (!bgOption) bgOption= XGetDefault(display, "Gist", "Background");
  fontOption= XGetDefault(display, "Gist", "font");
  if (!fontOption) fontOption= XGetDefault(display, "Gist", "Font");

  /* Fill in fonts available to this server.  */
  GxGrabFonts(xdpy, fontOption);

  /* Fill in screen specific data for each screen of this display.  */
  for (i=0 ; i<nScreens ; i++) {
    screens[i].owner= xdpy;
    screens[i].display= display;
    screens[i].maxRequest= XMaxRequestSize(display);
    screens[i].screen= i;
    screens[i].root= RootWindow(display, i);
    screens[i].width= DisplayWidth(display, i);
    screens[i].height= DisplayHeight(display, i);
    screens[i].colormap= colormap= DefaultColormap(display, i);

    /* Scan through list of visuals to find this one.  */
    for (j=0 ; j<nVisuals ; j++)
      if (vinfo[j].visual==DefaultVisual(display,i)) break;
    if (j<nVisuals) screens[i].v= vinfo+j;
    else screens[i].v= vinfo; /* this is actually a fairly serious bug... */

    /* If this is DirectColor, fill out DirectColorInfo.  This is also
       helpful for TrueColor.  */
    if (screens[i].v->class==DirectColor || screens[i].v->class==TrueColor)
      GxDirectColor(screens[i].v, &screens[i].dcinfo);

    /* Get named foreground and background preferences. If none, default
       background to white, foreground to black.  */
    if (!fgOption ||
	!XAllocNamedColor(display, colormap, fgOption,
			  &screens[i].stdColors[1], &color)) {
      screens[i].stdColors[1].pixel= BlackPixel(display, i);
      XQueryColor(display, colormap, &screens[i].stdColors[1]);
    }
    if (!bgOption ||
	!XAllocNamedColor(display, colormap, bgOption,
			  &screens[i].stdColors[0], &color)) {
      screens[i].stdColors[0].pixel= WhitePixel(display, i);
      XQueryColor(display, colormap, &screens[i].stdColors[0]);
    }
    screens[i].stdColors[2].pixel= BlackPixel(display, i);
    XQueryColor(display, colormap, &screens[i].stdColors[2]);
    screens[i].stdColors[3].pixel= WhitePixel(display, i);
    XQueryColor(display, colormap, &screens[i].stdColors[3]);

    /* Get primary colors.  If don't exist, default to foreground.  */
    for (j=4 ; j<10 ; j++) screens[i].stdColors[j]= screens[i].stdColors[1];
    if (!XAllocNamedColor(display, colormap, "red",
			  &screens[i].stdColors[4], &color))
      screens[i].stdColors[4]= screens[i].stdColors[1];
    if (!XAllocNamedColor(display, colormap, "green",
			  &screens[i].stdColors[5], &color))
      screens[i].stdColors[5]= screens[i].stdColors[1];
    if (!XAllocNamedColor(display, colormap, "blue",
			  &screens[i].stdColors[6], &color))
      screens[i].stdColors[6]= screens[i].stdColors[1];
    if (!XAllocNamedColor(display, colormap, "cyan",
			  &screens[i].stdColors[7], &color))
      screens[i].stdColors[7]= screens[i].stdColors[1];
    if (!XAllocNamedColor(display, colormap, "magenta",
			  &screens[i].stdColors[8], &color))
      screens[i].stdColors[8]= screens[i].stdColors[1];
    if (!XAllocNamedColor(display, colormap, "yellow",
			  &screens[i].stdColors[9], &color))
      screens[i].stdColors[9]= screens[i].stdColors[1];

    screens[i].rot_normal= screens[i].rot_rotated= None;
    screens[i].rot_gc= 0;
  }

  /* Add this display to the list of all displays.  */
  xdpy->next= gistX;
  gistX= xdpy;

  /* Return default (or requested) screen.  */
  return xdpy->screens+DefaultScreen(display);
}

/* ------------------------------------------------------------------------ */

int GxDisconnect(GxScreen *xscr)
{
  int i, j, n;
  GxDisplay *owner= xscr? xscr->owner : 0;
  unsigned long *pixels, black, white;
  XColor *stdColors;

  /* Decrement reference counter, and quit if any still outstanding.  */
  if (!owner || (owner->references--)>0) return 0;

  for (i=0 ; i<owner->nScreens ; i++) {
    n= 0;
    stdColors= owner->screens[i].stdColors;
    black= stdColors[2].pixel;
    white= stdColors[3].pixel;
    pixels= &stdColors[0].pixel;
    if (stdColors[0].pixel!=black && stdColors[0].pixel!=white)
      pixels[n++]= stdColors[0].pixel;
    if (stdColors[1].pixel!=black && stdColors[1].pixel!=white)
      pixels[n++]= stdColors[1].pixel;
    for (j=4 ; j<10 ; j++) pixels[n++]= stdColors[j].pixel;
    XFreeColors(owner->display, owner->screens[i].colormap, pixels, n, 0L);

    if (owner->screens[i].rot_normal!=None)
      XFreePixmap(owner->display, owner->screens[i].rot_normal);
    if (owner->screens[i].rot_rotated!=None)
      XFreePixmap(owner->display, owner->screens[i].rot_rotated);
    if (owner->screens[i].rot_gc)
      XFreeGC(owner->display, owner->screens[i].rot_gc);
  }

  /* Free storage associated with GxDisplay structure.  */
  GmFree(owner->normalizedName);
  GmFree(owner->screens);
  XFree((char *)owner->v);
  for (i=0 ; i<5 ; i++) {
    if (!owner->loadedFonts[i]) break;
    XFreeFont(owner->display, owner->loadedFonts[i]);
  }
  if (owner->permFont) XFreeFont(owner->display, owner->permFont);

  /* Unlink owner from list of open displays.  */
  GxUnlink(owner);

  /* Close the connection to the server and free the GxDisplay itself.  */
  XCloseDisplay(owner->display);
  GmFree(owner);
  return 1;
}

/* ------------------------------------------------------------------------ */

static int argcGist= 0;
static char **argvGist= 0;

void GxInitialize(int *argc, char **argv)
{
  if (argcGist) return;
  argcGist= *argc;
  argvGist= argv;
}

/* ------------------------------------------------------------------------ */

#ifndef X11R3
/* R4 standard is supposedly designed to be binary compatible with
   future releases, despite anticipated increases in number of hints.  */
static XSizeHints *sizeHints= 0;
static XWMHints *wmHints= 0;
static XClassHint *classHint= 0;

#else
/* R3 is not compatible with future releases; must recompile as R4.  */
static XSizeHints sizeHints;
static XWMHints wmHints;
static XClassHint classHint;
#endif

static char defaultName[]= "Gist Graphics";

extern int gist_input_hint;
int gist_input_hint= 1;

void GxSetProperties(char *name, Display *display, Window window,
		     unsigned int winx, unsigned int winy)
{
  char *pname= name? name : defaultName;

#ifndef X11R3
  {   /* recommended sequence for R4 */
    Atom wmDeleteAtom;
    XTextProperty xName, *xxName;
    if (XStringListToTextProperty(&pname, 1, &xName)==0) xxName= 0;
    else xxName= &xName;

    /* if (!sizeHints) sizeHints= XAllocSizeHints();
     * position hints confuse fvwm into placing window at (0,0)
     * -- supposedly these are obsolete anyway */
    if (!wmHints) wmHints= XAllocWMHints();
    if (!classHint) classHint= XAllocClassHint();
    if (sizeHints) {
      sizeHints->x= sizeHints->y= 0;
      sizeHints->width= winx;
      sizeHints->height= winy;
      sizeHints->flags= PPosition | PSize;
    }
    if (wmHints) {
      wmHints->initial_state= NormalState;
      wmHints->input= gist_input_hint? True : False;
      wmHints->flags= StateHint | InputHint;
    }
    if (classHint) {
      classHint->res_name= 0;
      classHint->res_class= "Gist";
    }

    XSetWMProperties(display, window, xxName, xxName, argvGist, argcGist,
		     sizeHints, wmHints, classHint);

    if (xxName) XFree((char *)xName.value);

    wmDeleteAtom= XInternAtom(display, "WM_DELETE_WINDOW", False);
    if (wmDeleteAtom!=None)
      XSetWMProtocols(display, window, &wmDeleteAtom, 1);
  }

#else
  {   /* recommended sequence for R3 */
    sizeHints.x= sizeHints.y= 0;
    sizeHints.width= winx;
    sizeHints.height= winy;
    sizeHints.flags= PPosition | PSize;
    XSetStandardProperties(display, window, &pname, &pname, None,
			   argvGist, argcGist, &sizeHints);

    wmHints.initial_state= NormalState;
    wmHints.input= True;
    wmHints.flags= StateHint | InputHint;
    XSetWMHints(display, window, &wmHints);

    classHint.res_name= 0;
    classHint.res_class= "Gist";
    XSetClassHint(display, window, &classHint);
  }
#endif
}

/* ------------------------------------------------------------------------ */
