/*
 * GREAD.C
 *
 * $Id$
 *
 * Define Drawing gread read routine for GIST
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "gist.h"
#include <stdio.h>

extern void GdKillSystems(void);  /* defined in draw.c */

#ifndef GISTPATH
#define GISTPATH "~/Gist:/usr/local/Gist"
#endif
char *gistPathDefault= GISTPATH;

/* ------------------------------------------------------------------------ */

#ifdef STDC_HEADERS
#include "string.h"
#else
#ifndef SIZE_T_TYPE
#define SIZE_T_TYPE unsigned long
#endif
extern SIZE_T_TYPE strlen(const char *);
extern char *strcpy(char *, const char *);
extern char *strncpy(char *, const char *, SIZE_T_TYPE);
extern char *strcat(char *, const char *);
extern char *strncat(char *, const char *, SIZE_T_TYPE);
extern int strcmp(const char *, const char *);
extern SIZE_T_TYPE strspn(const char *, const char *);
extern SIZE_T_TYPE strcspn(const char *, const char *);
extern char *strtok(char *, const char *);
#endif

extern double strtod(const char *, char **);
extern long strtol(const char *, char **, int);

extern char *getenv(const char *);

struct GsysRead {
  char *legend;
  GpBox viewport;
  GaTickStyle ticks;
} modelSystem, tempSystem;

struct GlegRead {
  GpReal x, y, dx, dy;
  GpTextAttribs textStyle;
  int nchars, nlines, nwrap;
} modelLegends;

static char *FormGistPath(void);
static char *FormFullName(char *gistPath, const char *name);
static void FormatError(FILE *fp, const char *name, const char *id);
static int SnarfColor(char *token);
static int SnarfRGB(char *token, GpColorCell *cell);
static int SnarfGray(GpColorCell *cell, int lookingFor4);
static char *WhiteSkip(char *input);
static char *DelimitRead(char *input, int *closed, int nlOK);
static char *IntRead(char *input, int *dest);
static char *RealRead(char *input, GpReal *dest);
static char *StringRead(char *input, char **dest);
static char *MemberRead(char *input, char **member);
static char *ArrayRead(char *input, GpReal *dest, int narray);
static char *LineRead(char *input, GpLineAttribs *dest);
static char *TextRead(char *input, GpTextAttribs *dest);
static char *AxisRead(char *input, GaAxisStyle *dest);
static char *TickRead(char *input, GaTickStyle *dest);
static char *SystemRead(char *input, struct GsysRead *dest);
static char *LegendsRead(char *input, struct GlegRead *dest);

/* ------------------------------------------------------------------------ */
/* A palette file (.gp) or style file (.gs) will be found if it:

   1. Is in the current working directory.
   2. Is the first file of its name encountered in any of the directories
      named in the GISTPATH environment variable.
   3. Ditto for the GISTPATH default string built in at compile time.
      Note that the environment variable is in addition to the compile
      time variable, not in place of it.

   The path name list should consist of directory names (with or without
   trailing '/'), separated by ':' with no intervening whitespace.  The
   symbol '~', if it is the first symbol of a directory name, will be
   expanded to the value of the HOME environment variable, but other
   environment variable expansions are not recognized.

   If the given filename begins with '/' the path search is
   not done.  '~' is NOT recognized in the given filename.
 */

static char *scratch;

static char *FormGistPath(void)
{
  char *gistPath= getenv("GISTPATH");
  int len= gistPath? strlen(gistPath) : 0;
  char *place;

  /* Get enough scratch space to hold 1023 characters of current path name,
     plus the concatenation of the GISTPATH environment variable and the
     GISTPATH compile-time option */
  scratch= GmMalloc(len+strlen(gistPathDefault)+1026);
  if (!scratch) return 0;

  place= scratch+1024;
  if (gistPath) {
    strcpy(place, gistPath);
    place+= len;
    *place= ':';
    place++;
  }
  strcpy(place, gistPathDefault);

  return scratch+1024;
}

static char *FormFullName(char *gistPath, const char *name)
{
  int nlen= strlen(name);
  int len, elen;
  char *now= scratch;

  for (;;) {
    /* Skip past any components of the GISTPATH which result in impossibly
       long path names */
    do len= strcspn(gistPath, ":"); while (!len);
    if (!len) break;
    elen= len;

    now= scratch;
    if (gistPath[0]=='~') {
      /* Get name of home directory from HOME environment variable */
      char *home= getenv("HOME");
      int hlen;
      if (home && (hlen= strlen(home))<1024) {
	strcpy(now, home);
	now+= hlen;
	gistPath++;
	len--;
	elen+= hlen-1;
      }
    }

    if (elen+nlen<1023) break;

    gistPath+= len+1;
  }

  if (len) {
    strncpy(now, gistPath, len);
    now+= len;
    if (now[-1]!='/') *now++= '/';
    strcpy(now, name);
  } else {
    scratch[0]= '\0';
  }

  return gistPath+len + strspn(gistPath+len, ":");
}

extern FILE *GistOpen(const char *name);
FILE *GistOpen(const char *name)
{
  FILE *f;
  if (!name) return 0;

  f= fopen(name, "r");

  if (!f && name[0]!='/') {
    /* Try to find relative file names somewhere on GISTPATH or, failing
       that, in the default directory specified at compile time.  */
    char *gistPath= FormGistPath();
    if (gistPath) {
      do {
	gistPath= FormFullName(gistPath, name);
	f= fopen(scratch, "r");
      } while (!f && gistPath[0]);
      GmFree(scratch);
    }
  }

  if (!f) {
    strcpy(gistError, "unable to open file ");
    strncat(gistError, name, 100);
  }
  return f;
}

static void FormatError(FILE *fp, const char *name, const char *id)
{
  fclose(fp);
  strcpy(gistError, id);
  strcat(gistError, " file format error in ");
  strncat(gistError, name, 127-strlen(gistError));
}

static char line[137];  /* longest allowed line is 136 characters */

/* ------------------------------------------------------------------------ */

static int SnarfColor(char *token)
     /* returns -1 if not unsigned char, -2 if missing */
{
  int color;
  char *suffix;

  if (!token) token= strtok(0, " \t\n");
  if (token) color= (int)strtol(token, &suffix, 0);
  else return -2;
  if (suffix==token || color<0 || color>255) return -1;
  else return color;
}

static int SnarfRGB(char *token, GpColorCell *cell)
{
  int red, blue, green;
  red= SnarfColor(token);
  if (red<0) return 1;
  green= SnarfColor(0);
  if (green<0) return 1;
  blue= SnarfColor(0);
  if (blue<0) return 1;
  cell->red= red;
  cell->green= green;
  cell->blue= blue;
  return 0;
}

static int SnarfGray(GpColorCell *cell, int lookingFor4)
{
  int gray= SnarfColor(0);
  if (gray==-2) return lookingFor4;
  else if (gray<0 || !lookingFor4) return 1;
  cell->gray= gray;
  return 0;
}

int GpReadPalette(Engine *engine, const char *gpFile,
		  GpColorCell **palette, int maxColors)
{
  char *token, *suffix;
  GpColorCell *pal= 0;
  int iColor= -1,  nColors= 0,  ntsc= 0,  lookingFor4= 0;
  FILE *gp= GistOpen(gpFile);

  *palette= 0;
  if (!gp) return 0;

  for (;;) {  /* loop on lines in file */
    token= fgets(line, 137, gp);
    if (!token) break;                      /* eof (or error) */

    token= strtok(token, " =\t\n");
    if (!token || token[0]=='#') continue;  /* blank or comment line */

    if (iColor<=0) {
      int *dest= 0;
      if (strcmp(token, "ncolors")==0) dest= &nColors;
      else if (strcmp(token, "ntsc")==0) dest= &ntsc;

      if (dest) {
	/* this is ncolors=... or ntsc=... line */
	token= strtok(0, " =\t\n");
	if (token) *dest= (int)strtol(token, &suffix, 0);
	else goto err;
	if (suffix==token || strtok(0, " \t\n")) goto err;

      } else {
	/* this must be the first rgb line */
	int gray;

	/* previous ncolors= is mandatory so palette can be allocated */
	if (nColors<=0) goto err;
	pal= GmMalloc(sizeof(GpColorCell)*nColors);
	if (!pal) goto merr;

	/* if first rgb line has 4 numbers, all must have 4, else 3 */
	if (SnarfRGB(token, pal)) goto err;
	gray= SnarfColor(0);
	if (gray==-1) goto err;
	if (gray>=0) {
	  lookingFor4= 1;
	  pal->gray= gray;
	  if (SnarfGray(pal, 0)) goto err;  /* just check for eol */
	} else {
	  lookingFor4= 0;
	  /* already got eol */
	}

	iColor= 1;
      }

    } else if (iColor<nColors) {
      /* read next rgb line */
      if (SnarfRGB(token, pal+iColor)) goto err;
      if (SnarfGray(pal+iColor, lookingFor4)) goto err;
      iColor++;

    } else {
      goto err;                  /* too many rgb for specified ncolors */
    }
  }
  if (iColor<nColors) goto err;  /* too few rgb for specified ncolors */

  fclose(gp);

  if (nColors>maxColors && maxColors>1) {
    /* attempt to rescale the palette to maxColors */
    int oldCell, newCell, nextCell;
    double ratio= ((double)(nColors-1))/((double)(maxColors-1));
    double frac, frac1, old= 0.0;
    for (newCell=0 ; newCell<maxColors ; newCell++) {
      oldCell= (int)old;
      nextCell= oldCell+1;
      if (nextCell>=nColors) nextCell= oldCell;
      frac= old-(double)oldCell;
      frac1= 1.0-frac;
      pal[newCell].red= frac1*pal[oldCell].red+frac*pal[nextCell].red;
      pal[newCell].green= frac1*pal[oldCell].green+frac*pal[nextCell].green;
      pal[newCell].blue= frac1*pal[oldCell].blue+frac*pal[nextCell].blue;
      if (!lookingFor4)
	pal[newCell].gray= frac1*pal[oldCell].gray+frac*pal[nextCell].gray;
      old+= ratio;
    }
    nColors= maxColors;
  }

  if (!lookingFor4) {
    /* gray values were not explicitly specified */
    if (ntsc) GpPutNTSC(nColors, pal);
    else GpPutGray(nColors, pal);
  }

  *palette= pal;
  iColor= GpSetPalette(engine, pal, nColors);
  return iColor>nColors? nColors : iColor;

 err:
  FormatError(gp, gpFile, "palette");
  if (pal) GmFree(pal);
  return 0;

 merr:
  strcpy(gistError, "memory manager failed to get space for palette");
  fclose(gp);
  return 0;
}

/* ------------------------------------------------------------------------ */

#define OPEN_BRACE '{'
#define CLOSE_BRACE '}'

static FILE *gs= 0;

static char *WhiteSkip(char *input)
{
  input+= strspn(input, " \t\n");

  while (!input[0] || input[0]=='#') { /* rest of line missing or comment */
    input= fgets(line, 137, gs);
    if (input) input= line + strspn(line, " \t\n");
    else break;
  }

  return input;
}

static char *DelimitRead(char *input, int *closed, int nlOK)
{
  int nlFound= 0;

  if (nlOK) {
    input+= strspn(input, " \t");
    if (*input=='\n' || *input=='\0') nlFound= 1;
  }

  input= WhiteSkip(input);
  if (input) {
    if (*input == CLOSE_BRACE) {
      *closed= 1;
      input++;
    } else {
      *closed= 0;
      if (*input == ',') {
	input++;
      } else {
	if (!nlOK || !nlFound) input= 0;
      }
    }

  } else {
    /* distinguish end-of-file from comma not found */
    *closed= 1;
  }

  return input;
}

static char *IntRead(char *input, int *dest)
{
  int value;
  char *suffix;

  input= WhiteSkip(input);  /* may be on a new line */
  value= (int)strtol(input, &suffix, 0);
  if (suffix==input) return 0;

  *dest= value;
  return suffix;
}

static char *RealRead(char *input, GpReal *dest)
{
  GpReal value;
  char *suffix;

  input= WhiteSkip(input);  /* may be on a new line */
  value= (GpReal)strtod(input, &suffix);
  if (suffix==input) return 0;

  *dest= value;
  return suffix;
}

char legendString[41];

static char *StringRead(char *input, char **dest)
{
  input= WhiteSkip(input);
  if (input) {
    if (*input=='0') {
      *dest= 0;
      input++;
    } else if (*input=='\"') {
      long len= strcspn(++input, "\"");
      int nc= len>40? 40 : len;
      strncpy(legendString, input, nc);
      input+= len;
      if (*input=='\"') { *dest= legendString;  input++; }
      else input= 0;
    } else {
      input= 0;
    }
  }
  return input;
}

static char *MemberRead(char *input, char **member)
{
  input= WhiteSkip(input);
  *member= input;
  if (input) {
    int gotEqual= 0;
    input+= strcspn(input, "= \t\n");
    if (*input == '=') gotEqual= 1;
    if (*input) *input++= '\0';
    if (!gotEqual) {
      input= WhiteSkip(input);
      if (input && *input++!='=') input= 0;
    }
  }
  return input;
}

static char *ArrayRead(char *input, GpReal *dest, int narray)
{
  int foundClose;

  input= WhiteSkip(input);
  if (!input) return 0;

  if (*input++ != OPEN_BRACE) return 0;  /* no open brace */
  input= WhiteSkip(input);
  if (!input) return 0;                  /* eof after open brace */

  for (narray-- ; ; narray--) {
    if (narray<0) return 0;           /* too many numbers in aggregate */

    input= RealRead(input, dest++);
    if (!input) return 0;             /* token was not a number */

    input= DelimitRead(input, &foundClose, 0);
    if (!input) return 0;             /* neither comma nor close brace */
    if (foundClose) break;
  }

  return input;
}

static char *LineRead(char *input, GpLineAttribs *dest)
{
  int foundClose;
  char *member;

  input= WhiteSkip(input);
  if (!input || *input++!=OPEN_BRACE) return 0;

  for (;;) {
    input= MemberRead(input, &member);
    if (!input) return 0;             /* couldn't find member = */

    if (strcmp(member, "color")==0) {
      input= IntRead(input, &dest->color);
    } else if (strcmp(member, "type")==0) {
      input= IntRead(input, &dest->type);
    } else if (strcmp(member, "width")==0) {
      input= RealRead(input, &dest->width);
    } else {
      return 0;                       /* unknown member */
    }
    if (!input) return 0;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) return 0;             /* not comma, nl, or close brace */
    if (foundClose) break;
  }

  return input;
}

static char *TextRead(char *input, GpTextAttribs *dest)
{
  int foundClose;
  char *member;
  int ijunk;
  GpReal rjunk;

  input= WhiteSkip(input);
  if (!input || *input++!=OPEN_BRACE) return 0;

  for (;;) {
    input= MemberRead(input, &member);
    if (!input) return 0;             /* couldn't find member = */

    if (strcmp(member, "color")==0) {
      input= IntRead(input, &dest->color);
    } else if (strcmp(member, "font")==0) {
      input= IntRead(input, &dest->font);
    } else if (strcmp(member, "prec")==0) {
      input= IntRead(input, &ijunk);
    } else if (strcmp(member, "height")==0) {
      input= RealRead(input, &dest->height);
    } else if (strcmp(member, "expand")==0) {
      input= RealRead(input, &rjunk);
    } else if (strcmp(member, "spacing")==0) {
      input= RealRead(input, &rjunk);
    } else if (strcmp(member, "upX")==0) {
      input= RealRead(input, &rjunk);
    } else if (strcmp(member, "upY")==0) {
      input= RealRead(input, &rjunk);
    } else if (strcmp(member, "path")==0 || strcmp(member, "orient")==0) {
      input= IntRead(input, &dest->orient);
    } else if (strcmp(member, "alignH")==0) {
      input= IntRead(input, &dest->alignH);
    } else if (strcmp(member, "alignV")==0) {
      input= IntRead(input, &dest->alignV);
    } else if (strcmp(member, "opaque")==0) {
      input= IntRead(input, &dest->opaque);
    } else {
      return 0;                       /* unknown member */
    }
    if (!input) return 0;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) return 0;             /* not comma, nl, or close brace */
    if (foundClose) break;
  }

  return input;
}

static char *AxisRead(char *input, GaAxisStyle *dest)
{
  int foundClose;
  char *member;

  input= WhiteSkip(input);
  if (!input || *input++!=OPEN_BRACE) return 0;

  for (;;) {
    input= MemberRead(input, &member);
    if (!input) return 0;             /* couldn't find member = */

    if (strcmp(member, "nMajor")==0) {
      input= RealRead(input, &dest->nMajor);
    } else if (strcmp(member, "nMinor")==0) {
      input= RealRead(input, &dest->nMinor);
    } else if (strcmp(member, "logAdjMajor")==0) {
      input= RealRead(input, &dest->logAdjMajor);
    } else if (strcmp(member, "logAdjMinor")==0) {
      input= RealRead(input, &dest->logAdjMinor);
    } else if (strcmp(member, "nDigits")==0) {
      input= IntRead(input, &dest->nDigits);
    } else if (strcmp(member, "gridLevel")==0) {
      input= IntRead(input, &dest->gridLevel);
    } else if (strcmp(member, "flags")==0) {
      input= IntRead(input, &dest->flags);
    } else if (strcmp(member, "tickOff")==0) {
      input= RealRead(input, &dest->tickOff);
    } else if (strcmp(member, "labelOff")==0) {
      input= RealRead(input, &dest->labelOff);
    } else if (strcmp(member, "tickLen")==0) {
      input= ArrayRead(input, dest->tickLen, 5);
    } else if (strcmp(member, "tickStyle")==0) {
      input= LineRead(input, &dest->tickStyle);
    } else if (strcmp(member, "gridStyle")==0) {
      input= LineRead(input, &dest->gridStyle);
    } else if (strcmp(member, "textStyle")==0) {
      input= TextRead(input, &dest->textStyle);
    } else if (strcmp(member, "xOver")==0) {
      input= RealRead(input, &dest->xOver);
    } else if (strcmp(member, "yOver")==0) {
      input= RealRead(input, &dest->yOver);
    } else {
      return 0;                       /* unknown member */
    }
    if (!input) return 0;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) return 0;             /* not comma, nl, or close brace */
    if (foundClose) break;
  }

  return input;
}

static char *TickRead(char *input, GaTickStyle *dest)
{
  int foundClose;
  char *member;

  input= WhiteSkip(input);
  if (!input || *input++!=OPEN_BRACE) return 0;

  for (;;) {
    input= MemberRead(input, &member);
    if (!input) return 0;             /* couldn't find member = */

    if (strcmp(member, "horiz")==0) {
      input= AxisRead(input, &dest->horiz);
    } else if (strcmp(member, "vert")==0) {
      input= AxisRead(input, &dest->vert);
    } else if (strcmp(member, "frame")==0) {
      input= IntRead(input, &dest->frame);
    } else if (strcmp(member, "frameStyle")==0) {
      input= LineRead(input, &dest->frameStyle);
    } else {
      return 0;                       /* unknown member */
    }
    if (!input) return 0;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) return 0;             /* not comma, nl, or close brace */
    if (foundClose) break;
  }

  return input;
}

/* defaultSystem is initialized to reasonable value for portrait mode */
#define DEF_XMIN 0.25
#define DEF_XMAX 0.60
#define DEF_YMIN 0.50
#define DEF_YMAX 0.85
struct GsysRead defaultSystem= {
  0, { DEF_XMIN, DEF_XMAX, DEF_YMIN, DEF_YMAX },
  {

    {7.5, 50., 1.2, 1.2, 3, 1, TICK_L|TICK_U|TICK_OUT|LABEL_L,
     0.0, 14.0*ONE_POINT,
     {12.*ONE_POINT, 8.*ONE_POINT, 5.*ONE_POINT, 3.*ONE_POINT, 2.*ONE_POINT},
     {FG_COLOR, L_SOLID, DEFAULT_LINE_WIDTH},
     {FG_COLOR, L_DOT, DEFAULT_LINE_WIDTH},
     {FG_COLOR, T_HELVETICA, 14.*ONE_POINT, TX_RIGHT, TH_NORMAL, TV_NORMAL, 1},
     0.5*(DEF_XMIN+DEF_XMAX), DEF_YMIN-52.*ONE_POINT},

    {7.5, 50., 1.2, 1.2, 4, 1, TICK_L|TICK_U|TICK_OUT|LABEL_L,
     0.0, 14.0*ONE_POINT,
     {12.*ONE_POINT, 8.*ONE_POINT, 5.*ONE_POINT, 3.*ONE_POINT, 2.*ONE_POINT},
     {FG_COLOR, L_SOLID, DEFAULT_LINE_WIDTH},
     {FG_COLOR, L_DOT, DEFAULT_LINE_WIDTH},
     {FG_COLOR, T_HELVETICA, 14.*ONE_POINT, TX_RIGHT, TH_NORMAL, TV_NORMAL, 1},
     DEF_XMIN, DEF_YMIN-52.*ONE_POINT},

  0, {FG_COLOR, L_SOLID, DEFAULT_LINE_WIDTH}
  }
};

struct GlegRead defaultLegends[2]= {
  /* Ordinary legends form two 36x22 character columns below viewport */
  { 0.5*ONE_INCH, DEF_YMIN-64.*ONE_POINT, 3.875*ONE_INCH, 0.0,
    {FG_COLOR, T_COURIER, 12.*ONE_POINT, TX_RIGHT, TH_LEFT, TV_TOP, 1},
    36, 22, 2 },
  /* Contour legends get a single 14x28 column to left of viewport */
  { DEF_XMAX+14.*ONE_POINT, DEF_YMAX+12.*ONE_POINT, 0.0, 0.0,
    {FG_COLOR, T_COURIER, 12.*ONE_POINT, TX_RIGHT, TH_LEFT, TV_TOP, 1},
    14, 28, 1 },
};

static char *SystemRead(char *input, struct GsysRead *dest)
{
  int foundClose;
  char *member;

  input= WhiteSkip(input);
  if (!input || *input++!=OPEN_BRACE) return 0;

  for (;;) {
    input= MemberRead(input, &member);
    if (!input) return 0;             /* couldn't find member = */

    if (strcmp(member, "viewport")==0) {
      GpReal box[4];
      box[0]= box[1]= box[2]= box[3]= -1.0;
      input= ArrayRead(input, box, 4);
      if (box[3]<0.0) input= 0;       /* all four required */
      else {
	dest->viewport.xmin= box[0];
	dest->viewport.xmax= box[1];
	dest->viewport.ymin= box[2];
	dest->viewport.ymax= box[3];
      }
    } else if (strcmp(member, "ticks")==0) {
      input= TickRead(input, &dest->ticks);
    } else if (strcmp(member, "legend")==0) {
      input= StringRead(input, &dest->legend);
    } else {
      return 0;                       /* unknown member */
    }
    if (!input) return 0;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) return 0;             /* not comma, nl, or close brace */
    if (foundClose) break;
  }

  return input;
}

static char *LegendsRead(char *input, struct GlegRead *dest)
{
  int foundClose;
  char *member;

  input= WhiteSkip(input);
  if (!input || *input++!=OPEN_BRACE) return 0;

  for (;;) {
    input= MemberRead(input, &member);
    if (!input) return 0;             /* couldn't find member = */

    if (strcmp(member, "x")==0) {
      input= RealRead(input, &dest->x);
    } else if (strcmp(member, "y")==0) {
      input= RealRead(input, &dest->y);
    } else if (strcmp(member, "dx")==0) {
      input= RealRead(input, &dest->dx);
    } else if (strcmp(member, "dy")==0) {
      input= RealRead(input, &dest->dy);
    } else if (strcmp(member, "textStyle")==0) {
      input= TextRead(input, &dest->textStyle);
    } else if (strcmp(member, "nchars")==0) {
      input= IntRead(input, &dest->nchars);
    } else if (strcmp(member, "nlines")==0) {
      input= IntRead(input, &dest->nlines);
    } else if (strcmp(member, "nwrap")==0) {
      input= IntRead(input, &dest->nwrap);
    } else {
      return 0;                       /* unknown member */
    }
    if (!input) return 0;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) return 0;             /* not comma, nl, or close brace */
    if (foundClose) break;
  }

  return input;
}

int GdReadStyle(Drawing *drawing, const char *gsFile)
{
  int foundClose, sysIndex, landscape;
  char *input, *keyword;

  if (!gsFile) return 0;

  gs= GistOpen(gsFile);
  if (!gs) return 1;

  tempSystem= defaultSystem;
  landscape= 0;

  input= fgets(line, 137, gs);
  if (!input) goto err;                      /* eof (or error) */

  GdKillSystems();

  for (;;) {
    input= WhiteSkip(input);
    if (!input) break;

    input= MemberRead(input, &keyword);
    if (!input) goto err;             /* couldn't find keyword = */

    if (strcmp(keyword, "default")==0) {
      input= SystemRead(input, &tempSystem);
    } else if (strcmp(keyword, "system")==0) {
      modelSystem= tempSystem;
      input= SystemRead(input, &modelSystem);
      gistD.hidden= 0;
      gistD.legend= modelSystem.legend;
      sysIndex= GdNewSystem(&modelSystem.viewport, &modelSystem.ticks);
      if (sysIndex<0) return 1;
    } else if (strcmp(keyword, "landscape")==0) {
      input= IntRead(input, &landscape);
    } else if (strcmp(keyword, "legends")==0) {
      modelLegends= defaultLegends[0];
      input= LegendsRead(input, &modelLegends);
      if (input) GdLegendBox(0, modelLegends.x, modelLegends.y,
			     modelLegends.dx, modelLegends.dy,
			     &modelLegends.textStyle, modelLegends.nchars,
			     modelLegends.nlines, modelLegends.nwrap);
    } else if (strcmp(keyword, "clegends")==0) {
      modelLegends= defaultLegends[1];
      input= LegendsRead(input, &modelLegends);
      if (input) GdLegendBox(1, modelLegends.x, modelLegends.y,
			     modelLegends.dx, modelLegends.dy,
			     &modelLegends.textStyle, modelLegends.nchars,
			     modelLegends.nlines, modelLegends.nwrap);
    } else {
      goto err;                       /* unknown keyword */
    }
    if (!input) goto err;             /* illegal format */

    input= DelimitRead(input, &foundClose, 1);
    if (!input) {
      if (foundClose) break;
      goto err;                       /* not comma, nl, or eof */
    }
    if (foundClose) goto err;         /* close brace not legal here */
  }

  if (landscape) GdLandscape(1);
  fclose(gs);
  return 0;

 err:
  FormatError(gs, gsFile, "drawing style");
  return 1;
}

/* ------------------------------------------------------------------------ */
