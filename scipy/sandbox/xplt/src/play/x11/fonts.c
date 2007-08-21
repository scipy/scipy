/*
 * fonts.c -- $Id$
 * font management for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

/* see README in this directory for explanatory notes */

#include "config.h"
#include "playx.h"

#include "pstdlib.h"

#include <stdio.h>
#include <string.h>

/* x_parse_fonts is semi-private, used only in connect.c */
extern void x_parse_fonts(x_display *xdpy);

static int x_substitute(x_display *xdpy, int *font, int *pixsize);
static void x_analyze(x_display *xdpy, int fam);
static void tmp_free(void);

static int x_closest(int pixsize, int *sizes, int nsizes, int *ndx);
static int x_lookup(int pixsize, int *sizes, int nsizes);
static char *x_face(char *name, int *face);
static int x_insert(int pixsize, int *sizes, char **names, int nsizes);

void
x_parse_fonts(x_display *xdpy)
{
  x_analyze(xdpy, P_COURIER);
  x_analyze(xdpy, P_TIMES);
  x_analyze(xdpy, P_HELVETICA);
  x_analyze(xdpy, P_SYMBOL);
  x_analyze(xdpy, P_NEWCENTURY);
}

static char **tmp_fonts = 0;
static void
tmp_free(void)
{
  char **tmp = tmp_fonts;
  if (tmp) {
    tmp_fonts = 0;
    XFreeFontNames(tmp);
  }
}

XFontStruct *
x_font(x_display *xdpy, int font, int pixsize)
{
  char nm[128], *name;
  XFontStruct *f;
  if (tmp_fonts) tmp_free();

  if (font>=P_GUI_FONT || font<0 || pixsize<=0 || pixsize>180) {
    f = xdpy->font;

  } else {
    int siz, i, j, pass, *ip, *jp;
    for (pass=siz=0 ;; pass++) {
      ip = jp = 0;
      j = -1;
      i = xdpy->most_recent;
      while (i>=0) {
        if (!xdpy->cached[i].f) {
          /* handle interrupted unload operation */
          if (ip) *ip = -1;
          break;
        }
        if (xdpy->cached[i].font==font &&
            xdpy->cached[i].pixsize==pixsize) {
          if (ip) {
            /* CRITICAL section */
            *ip = xdpy->cached[i].next;
            xdpy->cached[i].next = xdpy->most_recent;
            xdpy->most_recent = i;
          }
          return xdpy->cached[i].f;
        }
        jp = ip;
        ip = &xdpy->cached[i].next;
        j = i;
        i = *ip;
      }
      if (pass) break;
      siz = x_substitute(xdpy, &font, &pixsize);
      if (font==P_GUI_FONT) return xdpy->font;
    }

    /* construct font name */
    name = xdpy->available[font].names[siz];
    if (!xdpy->available[font].sizes[siz]) {
      /* scalable, need to find specific instance */
      char *pnm = nm;
      int n = 7;
      while (n--) while ((*(pnm++)= *(name++))!='-');
      sprintf(pnm, "%d%n", pixsize, &n);
      strcpy(pnm+n, name);
      tmp_fonts = XListFonts(xdpy->dpy, nm, 4, &n);
      if (n<=0) return xdpy->font;  /* should never happen (X server bug) */
      strcpy(nm, tmp_fonts[0]);
      XFreeFontNames(tmp_fonts);
      tmp_free();
      name = nm;
    }

    /* should be able to load it */
    f = XLoadQueryFont(xdpy->dpy, name);
    if (!f) return xdpy->font;      /* should never happen (X server bug) */

    if (!xdpy->cached[0].f) {
      /* cache not yet full */
      for (j=0 ; j<N_FONT_CACHE-1 ; j++)
        if (xdpy->cached[j+1].f) break;
    } else {
      /* cache is full, need to unload one, j is least recent */
      XFontStruct *fold = xdpy->cached[j].f;
      xdpy->cached[j].f = 0;
      if (jp) *jp = -1;   /* jp pointed to j, now least recent */
      XFreeFont(xdpy->dpy, fold);
    }

    xdpy->cached[j].font = font;
    xdpy->cached[j].pixsize = pixsize;
    xdpy->cached[j].f = f;
    xdpy->cached[j].next = xdpy->most_recent;
    xdpy->most_recent = j;
  }

  if (p_signalling) p_abort();

  return f;
}

static int
x_substitute(x_display *xdpy, int *font, int *pixsize)
{
  int fnt = *font;
  int siz = *pixsize;
  int face = fnt&(P_BOLD|P_ITALIC);

  if (!xdpy->available[fnt].nsizes) {
    fnt ^= face;
    if (!face || !xdpy->available[fnt].nsizes) {
      int i;
      for (i=1 ; i<4 ; i++) if (xdpy->available[fnt|i].nsizes) break;
      if (i<4) {
        fnt |= i;
      } else {
        if (xdpy->available[P_TIMES|face].nsizes)
          fnt = P_TIMES|face;
        else if (xdpy->available[P_NEWCENTURY|face].nsizes)
          fnt = P_NEWCENTURY|face;
        else if (xdpy->available[P_HELVETICA|face].nsizes)
          fnt = P_HELVETICA|face;
        else if (xdpy->available[P_COURIER|face].nsizes)
          fnt = P_COURIER|face;
        else if (face) {
          if (xdpy->available[P_TIMES].nsizes)
            fnt = P_TIMES;
          else if (xdpy->available[P_NEWCENTURY].nsizes)
            fnt = P_NEWCENTURY;
          else if (xdpy->available[P_HELVETICA].nsizes)
            fnt = P_HELVETICA;
          else if (xdpy->available[P_COURIER].nsizes)
            fnt = P_COURIER;
          else
            fnt = P_GUI_FONT;
        } else
          fnt = P_GUI_FONT;
      }
    }
  }

  if (fnt!=P_GUI_FONT)
    siz = x_closest(siz, xdpy->available[fnt].sizes,
                    xdpy->available[fnt].nsizes, &face);
  else
    face = -1;

  *font = fnt;
  *pixsize = siz;

  return face;
}

/*
 -foundry-family-wgt-slant-wid--pixels-pts-hres-vres-spacing-avgwid-char-set
 */

static char *pattern[5] = {
  "-*-courier-*-*-normal--*-*-*-*-m-*-iso8859-1",
  "-*-times-*-*-normal--*-*-*-*-p-*-iso8859-1",
  "-*-helvetica-*-*-normal--*-*-*-*-p-*-iso8859-1",
  "-*-symbol-*-*-normal--*-*-*-*-p-*-*-*",
  "-*-new century schoolbook-*-*-normal--*-*-*-*-p-*-iso8859-1" };

static void
x_analyze(x_display *xdpy, int fam)
{
  int i, j, n, face, pixsize, nsizes;
  char *name;
  if (tmp_fonts) tmp_free();
  tmp_fonts = XListFonts(xdpy->dpy, pattern[((unsigned int)fam)>>2], 1024, &n);

  for (i=0 ; i<n ; i++) {
    name = x_face(tmp_fonts[i], &face);
    if (!name) continue;

    /* extract pixels field */
    pixsize = 0;
    if (name[0]!='*') while (name[0] && name[0]>='0' && name[0]<='9')
      pixsize = 10*pixsize + *(name++) - '0';
    else
      name++;
    if (name[0]!='-') continue;

    /* protect against superlong font names */
    if (!pixsize && strlen(tmp_fonts[i])>120) continue;

    face += fam;

    nsizes = xdpy->available[face].nsizes;
    if (x_lookup(pixsize, xdpy->available[face].sizes, nsizes)) continue;
    if (nsizes%12==0) {
      int *sizes = xdpy->available[face].sizes;
      char **names = xdpy->available[face].names;
      xdpy->available[face].sizes = p_realloc(sizes, sizeof(int)*(nsizes+12));
      if (!xdpy->available[face].sizes) {
        xdpy->available[face].sizes = sizes;
        return;
      }
      xdpy->available[face].names = p_realloc(names,sizeof(char*)*(nsizes+13));
      if (!xdpy->available[face].names) {
        xdpy->available[face].names = names;
        return;
      }
    }
    j = x_insert(pixsize, xdpy->available[face].sizes,
                 xdpy->available[face].names, nsizes);
    xdpy->available[face].nsizes++;
    if (pixsize) {
      xdpy->available[face].names[j] = p_strcpy(tmp_fonts[i]);
    } else {
      /* scalable font needs wildcard name */
      char nm[128], *pnm = nm;
      int n = 7;
      name = tmp_fonts[i];
      while (n--) while ((*(pnm++)= *(name++))!='-');
      /* skip over pixels, points fields */
      *(pnm++)= '*';   *(pnm++)= '-';  *(pnm++)= '*';   *(pnm++)= '-';
      for (n=2 ; n-- ;) while (name[0] && *(name++)!='-');
      /* copy hres, vres, spacing fields */
      for (n=3 ; n-- ;) while (name[0] && (*(pnm++)= *(name++))!='-');
      /* skip over average width field */
      *(pnm++)= '*';   *(pnm++)= '-';
      while (name[0] && *(name++)!='-');
      /* copy remainder (character set fields) */
      while ((*(pnm++)= *(name++)));
      xdpy->available[face].names[j] = p_strcpy(nm);
    }
  }

  tmp_free();
}

static char *
x_face(char *name, int *face)
{
  char *pat = "bold";

  /* skip to wgt field */
  int n = 3;
  while (n--) while (name[0] && *(name++)!='-');
  *face = 0;
  if (!name[0]) return 0;

  /* compare with "bold" -- otherwise assume non-bold */
  while (name[0] && name[0]==pat[0]) name++, pat++;
  if (!pat[0] && name[0]=='-') *face |= P_BOLD;

  /* skip to slant field, compare with "r" -- otherwise assume italic */
  while (name[0] && *(name++)!='-');
  if (!name[0] || !name[1]) return 0;
  if (name[0]!='r' || name[1]!='-') *face |= P_ITALIC;

  /* skip to pixels field */
  n = 3;
  while (n--) while (name[0] &&  *(name++)!='-');
  return name[0]? name : 0;
}

static int
x_insert(int pixsize, int *sizes, char **names, int nsizes)
{
  unsigned int i, j;
  char *nm;
  if (!nsizes || sizes[0]>pixsize) {
    j = 0;
  } else if (sizes[nsizes-1]<pixsize) {
    j = nsizes;
  } else {
    int k;
    for (i=0,j=nsizes-1,k=j>>1 ; k!=i ; k=(i+j)>>1) {
      if (sizes[k]>pixsize) j = k;
      else i = k;
    }
  }
  names[nsizes+1] = 0; /* see x_disconnect */
  for (i=nsizes ; i>j ; i--) {
    sizes[i] = sizes[i-1];
    nm = names[i-1];
    names[i-1] = 0;    /* never want two copies of nm in names[] */
    names[i] = nm;
  }
  sizes[j] = pixsize;
  return j;
}

static int
x_lookup(int pixsize, int *sizes, int nsizes)
{
  unsigned int i, j=nsizes-1, k;
  if (nsizes<=0 || sizes[0]>pixsize || sizes[j]<pixsize) return 0;
  if (sizes[j]==pixsize) return nsizes;
  if (!j) return 0;
  if (sizes[0]==pixsize) return 1;
  i = 0;
  for (k=j>>1 ; k!=i ; k=(i+j)>>1) {
    if (sizes[k]==pixsize) return k+1;
    if (sizes[k]>pixsize) j = k;
    else i = k;
  }
  return 0;
}

static int
x_closest(int pixsize, int *sizes, int nsizes, int *ndx)
{
  unsigned int i, j=nsizes-1, k;
  if (nsizes<=0 || pixsize<=0) return (*ndx=-1, 0);
  if (sizes[j]<=pixsize) return (*ndx=j, sizes[j]);
  else if (sizes[0]>=pixsize) return (*ndx=0, sizes[0]);
  i = 0;
  for (k=j>>1 ; k!=i ; k=(i+j)>>1) {
    if (sizes[k]==pixsize) return (*ndx=k, pixsize);
    if (sizes[k]>pixsize) j = k;
    else i = k;
  }
  if (pixsize-sizes[i] < sizes[j]-pixsize) return (*ndx=i, sizes[i]);
  else return (*ndx=j, sizes[j]);
}
