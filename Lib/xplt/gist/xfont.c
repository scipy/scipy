/*
 * XFONT.C
 *
 * $Id$
 *
 * Implement X windows font utilities for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "xfont.h"

#ifdef STDC_HEADERS
#include "string.h"
#else
#ifndef SIZE_T_TYPE
#define SIZE_T_TYPE unsigned long
#endif
extern char *strcpy(char *, const char *);
extern char *strcat(char *, const char *);
extern int strncmp(const char *, const char *, SIZE_T_TYPE);
#endif

static char *fallbacks[FONT_FALLBACKS]= { "variable", "9x15", "fixed" };

static char *foundry= "Adobe";

static char *facenames[FONT_FAMILIES]= {
  "Courier", "Times", "Helvetica", "Symbol", "New Century" };

static char *sizenames[FONT_SIZES]= { "8", "10", "12", "14", "18", "24" };

static char *slantnames[3]= { "roman", "italic", "oblique" };
static int slants[FONT_FAMILIES]= { 2, 1, 2, 0, 1 };

static char *weightnames[2]= { "medium", "bold" };

static char *familyNames[FONT_FAMILIES]= {
  "-adobe-courier-",
  "-adobe-times-",
  "-adobe-helvetica-",
  "-adobe-symbol-",
  "-adobe-new century schoolbook-" };

static char *weightNames[2]= { "medium-", "bold-" };

static char *slantNames[3]= { "r-", "i-", "o-" };

static char *sizeNames[FONT_SIZES]= {
  "normal--*-80-", "normal--*-100-", "normal--*-120-",
  "normal--*-140-", "normal--*-180-", "normal--*-240-" };

static char *resNames[2]= { "75-75-*", "100-100-*" };

static char fullName[80];

char *GxNameFont(int id)
{
  if (!(id & IS_FALLBACK)) {
    int font= id>>(X_SIZE_BITS+4);
    int slant= (id & (T_ITALIC<<(X_SIZE_BITS+2)))? 1 : 0;
    int bold= (id & (T_BOLD<<(X_SIZE_BITS+2)))? 1 : 0;
    int size= (id>>2) & X_SIZE_MASK;
    int dpi= (id & IS_100DPI)? 1 : 0;
    strcpy(fullName, familyNames[font]);
    strcat(fullName, weightNames[bold]);
    strcat(fullName, slantNames[slant? slants[font] : 0]);
    strcat(fullName, sizeNames[size]);
    strcat(fullName, resNames[dpi]);
  } else {
    int i= id>>(X_SIZE_BITS+2);
    if (i<FONT_FALLBACKS) strcpy(fullName, fallbacks[i]);
    else return 0;
  }
  return fullName;
}

static GpReal cutoff[FONT_SIZES-1]= { /* halfway between font sizes */
  9.0*ONE_POINT, 11.0*ONE_POINT, 13.0*ONE_POINT,
  16.0*ONE_POINT, 21.0*ONE_POINT };

int GxIDFont(GxDisplay *xdpy, int dpi, int font, GpReal height)
{
  int i, j, k, up, is, bold, slant;
  GxFontProps *fontProps= dpi==75? xdpy->fontProps75 : xdpy->fontProps100;
  GxFontProps *altProps= dpi==75? xdpy->fontProps100 : xdpy->fontProps75;
  int firstTry, alt= dpi==75? -1 : +1;

  /* Get closest available font size.  */
  for (is=0 ; is<FONT_SIZES-1 ; is++) if (height<cutoff[is]) break;

  /* Break down font into bold flag, slant flag, and typeface.  */
  bold= (T_BOLD & font)? 1 : 0;
  slant= (T_ITALIC & font)? 1 : 0;
  font= font>>2;

  /* Do "spiral search" on sizes, beginning with requested size, and
     then one smaller, one larger, two smaller, and so on.  */
  up= k= is? 0 : 1;
  j= 1;
  firstTry= 1;
  for (;;) {
    /* Check for requested font at this size.  */
    if (fontProps[font].foundry && fontProps[font].masks[is]) break;
    if (firstTry) {
      int is0= is;
      is+= alt;
      if (is<0 || is>=FONT_SIZES) is= is0;
      if (altProps[font].foundry && altProps[font].masks[is]) {
	dpi= dpi==75? 100 : 75;
	break;
      }
      is= is0;
      firstTry= 0;
    }

    /* Try for a different typeface at this size.
       Backward search on faces keeps variable spacing if possible.  */
    for (i=FONT_FAMILIES-1 ; i>=0 ; i--) {
      if (!fontProps[i].foundry || i==font || i==3) continue;
      if (fontProps[i].masks[is]) break;
    }
    if (i>=0) { font= i;  break; }

    /* Implement spiral search pattern.  */
    if (up) {
      is+= j;
      if (is>=FONT_SIZES) break;
      if (!k) {
	if (is-j-1>=0) { up= 0; j++; }
	else { k= 1; j= 1; }
      }
    } else {
      is-= j;
      if (is<0) break;
      if (!k) {
	if (is+j+1<FONT_SIZES) { up= 1; j++; }
	else { k= 1; j= 1; }
      }
    }
  }

  if (is>=0 && is<FONT_SIZES) {
    /* Try to get slant, then bold */
    int masks= fontProps[font].masks[is];

    if (!slant && (masks & (HAS_NORMAL | HAS_BOLD))) {
      /* Can get unslanted as requested */
      if (bold) {
	if (!(masks & HAS_BOLD)) bold= 0;
      } else {
	if (!(masks & HAS_NORMAL)) bold= 1;
      }
    } else if (masks & (HAS_SLANT | HAS_BOLD_SLANT)) {
      /* Slanted text exists, and either slanted requested
	 or unslanted unavailable */
      slant= 1;
      if (bold) {
	if (!(masks & HAS_BOLD_SLANT)) bold= 0;
      } else {
	if (!(masks & HAS_SLANT)) bold= 1;
      }
    } else {
      /* Slanted text requested but unavailable */
      slant= 0;
      if (bold) {
	if (!(masks & HAS_BOLD)) bold= 0;
      } else {
	if (!(masks & HAS_NORMAL)) bold= 1;
      }
    }

    /* form loadedID for this font */
    font= font<<2;
    if (slant) font|= T_ITALIC;
    if (bold) font|= T_BOLD;
    font= font<<X_SIZE_BITS;
    font|= is;
    font= font<<2;
    if (dpi!=75) font|= IS_100DPI;

  } else {
    /* No Adobe fonts available, go with a fallback font.  */
    int fallbacks= xdpy->fontFallbacks;
    if (font==0) { /* Courier, try fixed spacing first */
      for (i=1 ; i<FONT_FALLBACKS ; i++) if ((1<<i)&fallbacks) break;
      if (i>=FONT_FALLBACKS && (1&fallbacks)) i= 0;
    } else {       /* try variable spacing first */
      for (i=0 ; i<FONT_FALLBACKS ; i++) if ((1<<i)&fallbacks) break;
    }
    /* Note that i>=FONT_FALLBACKS if no fallbacks */

    /* form loadedID for this font */
    font= (i<<(X_SIZE_BITS+2)) | IS_FALLBACK;
  }

  return font;
}

char **GxFontFaces(GxFontProps *fontProps, int isize,
		   int slant, int bold, int *nfonts, int *mask)
{
  int i, j, check, msk, mski;

  if (slant>=0) {
    slant= slant? 1:0;
    if (bold>=0) check= (bold?4:1)<<slant;
    else check= 5<<slant;
  } else {
    if (bold>=0) check= 3<<bold;
    else check= 15;
  }

  msk= 0;
  for (i=0 ; i<FONT_FAMILIES ; i++) {
    if (!fontProps[i].foundry) continue;
    mski= 1<<i;
    if (isize>=0) {
      if (fontProps[i].masks[isize]&check) msk|= mski;
    } else {
      for (j=0 ; j<FONT_SIZES ; j++) {
	if (fontProps[i].masks[j]&check) break;
      }
      if (j<FONT_SIZES) msk|= mski;
    }
  }

  *mask= msk;
  *nfonts= FONT_FAMILIES;
  return facenames;
}

char **GxFontSizes(GxFontProps *fontProps, int ifont,
		   int slant, int bold, int *nsizes, int *mask)
{
  int i, j, check, msk, mski;

  if (slant>=0) {
    slant= slant? 1:0;
    if (bold>=0) check= (bold?4:1)<<slant;
    else check= 5<<slant;
  } else {
    if (bold>=0) check= 3<<bold;
    else check= 15;
  }

  msk= 0;
  for (i=0 ; i<FONT_SIZES ; i++) {
    mski= 1<<i;
    if (ifont>=0) {
      if (fontProps[ifont].masks[i]&check) msk|= mski;
    } else {
      for (j=0 ; j<FONT_FAMILIES ; j++) {
	if (!fontProps[j].foundry) continue;
	if (fontProps[j].masks[i]&check) break;
      }
      if (j<FONT_SIZES) msk|= mski;
    }
  }

  *mask= msk;
  *nsizes= FONT_SIZES;
  return sizenames;
}

char **GxFontSlants(GxFontProps *fontProps, int ifont,
		    int isize, int bold, int *nslants, int *mask)
{
  int i, j, k, slant, check, msk, mski;

  if (bold>=0) bold= bold?4:1;
  else bold= 5;

  msk= 0;
  for (i=0 ; i<3 ; i++) {
    mski= 1<<i;
    slant= i? 1:0;
    check= bold<<slant;
    if (ifont>=0) {
      if (isize>=0) {
	if (fontProps[ifont].masks[isize]&check) msk|= mski;
      } else {
	for (j=0 ; j<FONT_SIZES ; j++) {
	  if (fontProps[ifont].masks[j]&check) break;
	}
	if (j<FONT_SIZES) msk|= mski;
      }
    } else {
      for (k=0 ; k<FONT_FAMILIES ; k++) {
	if (!fontProps[k].foundry) continue;
	if (isize>=0) {
	  if (fontProps[k].masks[isize]&check) break;
	} else {
	  for (j=0 ; j<FONT_SIZES ; j++) {
	    if (fontProps[k].masks[j]&check) break;
	  }
	  if (j<FONT_SIZES) break;
	}
      }
      if (k<FONT_FAMILIES) msk|= mski;
    }
  }

  *mask= msk;
  *nslants= 3;
  return slantnames;
}

char **GxFontWeights(GxFontProps *fontProps, int ifont,
		     int isize, int slant, int *nweights, int *mask)
{
  int i, j, k, bold, check, msk, mski;

  if (slant>=0) slant= slant?2:1;
  else slant= 3;

  msk= 0;
  for (i=0 ; i<2 ; i++) {
    mski= 1<<i;
    bold= i? 2:0;
    check= slant<<bold;
    if (ifont>=0) {
      if (isize>=0) {
	if (fontProps[ifont].masks[isize]&check) msk|= mski;
      } else {
	for (j=0 ; j<FONT_SIZES ; j++) {
	  if (fontProps[ifont].masks[j]&check) break;
	}
	if (j<FONT_SIZES) msk|= mski;
      }
    } else {
      for (k=0 ; k<FONT_FAMILIES ; k++) {
	if (!fontProps[k].foundry) continue;
	if (isize>=0) {
	  if (fontProps[k].masks[isize]&check) break;
	} else {
	  for (j=0 ; j<FONT_SIZES ; j++) {
	    if (fontProps[k].masks[j]&check) break;
	  }
	  if (j<FONT_SIZES) break;
	}
      }
      if (k<FONT_FAMILIES) msk|= mski;
    }
  }

  *mask= msk;
  *nweights= 2;
  return weightnames;
}

int gxFontSize[FONT_SIZES]= { 8, 10, 12, 14, 18, 24 };

/* ------------------------------------------------------------------------ */

static char *fontPats75[FONT_FAMILIES]= {
  "-adobe-courier-*-*-normal--*-*-75-75-*",
  "-adobe-times-*-*-normal--*-*-75-75-*",
  "-adobe-helvetica-*-*-normal--*-*-75-75-*",
  "-adobe-symbol-*-*-normal--*-*-75-75-*",
  "-adobe-new century schoolbook-*-*-normal--*-*-75-75-*" };

static char *fontPats100[FONT_FAMILIES]= {
  "-adobe-courier-*-*-normal--*-*-100-100-*",
  "-adobe-times-*-*-normal--*-*-100-100-*",
  "-adobe-helvetica-*-*-normal--*-*-100-100-*",
  "-adobe-symbol-*-*-normal--*-*-100-100-*",
  "-adobe-new century schoolbook-*-*-normal--*-*-100-100-*" };

static char *SkipHyphens(char *name, int nhyphens);
static int FindFontSize(char *font);
static int FindFontMask(char *font, int ifont);
static void DecodeFontList(char **fonts, int nfonts,
			   GxFontProps *fontProps, int ifont);

static char *SkipHyphens(char *name, int nhyphens)
{
  while (nhyphens>0) {
    while (*name && *name!='-') name++;
    if (!*name) break;
    name++;
    nhyphens--;
  }
  return name;
}

static int FindFontSize(char *font)
{
  int index;
  font= SkipHyphens(font, 8);
  if (*font=='8') {
    index= 0;
  } else if (*font=='1') {
    font++;
    if (*font=='0') index= 1;
    else if (*font=='2') index= 2;
    else if (*font=='4') index= 3;
    else if (*font=='8') index= 4;
    else index= FONT_SIZES;
  } else if (*font=='2') {
    font++;
    if (*font=='4') index= 5;
    else index= FONT_SIZES;
  } else if (*font=='0') {
    index= -1;
    font--;   /* incremented below */
  } else {
    index= FONT_SIZES;
  }
  if (index<FONT_SIZES) {
    font++;
    if (*font++!='0' || *font!='-') index= FONT_SIZES;
  }
  return index;
}

static int FindFontMask(char *font, int ifont)
{
  int bs;
  font= SkipHyphens(font, 3);

  if (strncmp(font, weightNames[0], 7)==0) bs= 0;
  else if (strncmp(font, weightNames[1], 5)==0) bs= 2;
  else return 0;

  font+= bs? 5 : 7;
  if (strncmp(font, slantNames[slants[ifont]], 2)==0) bs|= 1;
  else if (strncmp(font, slantNames[0], 2)!=0) return 0;

  return 1<<bs;
}

static void DecodeFontList(char **fonts, int nfonts,
			   GxFontProps *fontProps, int ifont)
{
  int i, j, mask, total= 0;

  fontProps[ifont].facename= facenames[ifont];
  fontProps[ifont].slantname= slantnames[slants[ifont]];

  for (i=0 ; i<FONT_SIZES ; i++) fontProps[ifont].masks[i]= 0;
  for (i=0 ; i<nfonts ; i++) {
    j= FindFontSize(fonts[i]);
    if (j>=FONT_SIZES) continue;
    if (j<0) {
      /* R5 scalable font available */
      for (j=0 ; j<FONT_SIZES ; j++) {
	mask= FindFontMask(fonts[i], ifont);
	fontProps[ifont].masks[j]|= mask;
	if (mask) total++;
      }
    } else {
      /* R4 fonts */
      mask= FindFontMask(fonts[i], ifont);
      fontProps[ifont].masks[j]|= mask;
      if (mask) total++;
    }
  }

  fontProps[ifont].foundry= total? foundry : 0;
}

void GxGrabFonts(GxDisplay *xdpy, char *permFont)
{
  Display *display= xdpy->display;
  char **fonts;
  int nfonts;
  int i;

  /* Get all 75 dpi fonts */
  for (i=0 ; i<FONT_FAMILIES ; i++) {
    fonts= XListFonts(display, fontPats75[i], 1024, &nfonts);
    DecodeFontList(fonts, nfonts, xdpy->fontProps75, i);
    XFreeFontNames(fonts);
  }

  /* Get all 100 dpi fonts */
  for (i=0 ; i<FONT_FAMILIES ; i++) {
    fonts= XListFonts(display, fontPats100[i], 1024, &nfonts);
    DecodeFontList(fonts, nfonts, xdpy->fontProps100, i);
    XFreeFontNames(fonts);
  }

  /* Get the fallback fonts */
  xdpy->fontFallbacks= 0;
  for (i=0 ; i<FONT_FALLBACKS ; i++) {
    if ((fonts= XListFonts(display, fallbacks[i], 1, &nfonts))) {
      XFreeFontNames(fonts);
      xdpy->fontFallbacks|= 1<<i;
    }
  }

  /* If no fallbacks, get the default font for screen 0 */
  if (xdpy->fontFallbacks==0) {
    xdpy->defaultFont= XQueryFont(xdpy->display,
				  XGContextFromGC(DefaultGC(xdpy->display,
							    0)));
  } else {
    xdpy->defaultFont= 0;
  }

  /* no fonts are loaded yet */
  xdpy->loadedFonts[0]= 0;

  /* get the permanent font (for use in menus, etc.) if possible */
  if (permFont)
    xdpy->permFont= XLoadQueryFont(xdpy->display, permFont);
  else if (xdpy->fontFallbacks & 4)
    xdpy->permFont= XLoadQueryFont(xdpy->display, fallbacks[2]);
  else if (xdpy->fontFallbacks & 2)
    xdpy->permFont= XLoadQueryFont(xdpy->display, fallbacks[1]);
  else if (xdpy->fontFallbacks & 1)
    xdpy->permFont= XLoadQueryFont(xdpy->display, fallbacks[0]);
  else
    xdpy->permFont= 0;  /* Use xdpy->defaultFont if required */
}
