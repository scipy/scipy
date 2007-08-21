/*
 * CGM.C
 *
 * $Id$
 *
 * Implement the CGM binary metafile engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "play.h"
#include "pstdlib.h"
#include "pstdio.h"
#include "cgm.h"
#include "gtext.h"

#include <string.h>
#include <time.h>

static char *cgmType= "binary CGM";

#define N_CGMFONTS 20

extern double floor(double);

#ifndef CREATE_CGM
#define CREATE_CGM(name) p_fopen(name, "wb")
#endif

/* ------------------------------------------------------------------------ */

static Octet *FormCommand(Octet *buffer, int klass, int id, long nbytes,
                          long *lpart);
static Octet *NextPartition(Octet *buffer, long nbytes, long *lpart);
static Octet *Pascalify(Octet *buffer, const char *text, long len, int pad);
static Octet *FormWords(Octet *buffer, const short *words, long nwords);
static void FormReal(short *fixed, GpReal x);

static int WriteB(p_file *file, void *buffer, long nbytes);
static void WriteError(CGMEngine *cgm, const char *msg);

static int BeginMetafile(CGMEngine *cgm);
static void SetPageDefaults(CGMEngine *cgmEngine);
static void SetCGMTransform(GpTransform *toPixels, int landscape,
                            GpReal scale);
static int BeginPage(CGMEngine *cgmEngine);
static void EndClip(CGMEngine *cgmEngine);
static void BeginClip(CGMEngine *cgmEngine, GpTransform *trans);
static int EndPage(CGMEngine *cgmEngine);

static int ChangePalette(Engine *engine);
static void Kill(Engine *engine);
static int Clear(Engine *engine, int always);
static int Flush(Engine *engine);

static int SetupColor(CGMEngine *cgmEngine, unsigned long color, int which);
static int SetupLine(CGMEngine *cgmEngine);
static int SetupEdge(CGMEngine *cgmEngine);
static void CheckClip(CGMEngine *cgmEngine);
static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth);
static int SetupMarker(CGMEngine *cgmEngine);
static int DrawMarkers(Engine *engine, long n, const GpReal *px,
                       const GpReal *py);
static int SetupText(CGMEngine *cgmEngine);
static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text);
static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py);
static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors);
static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy);

static void IncrementName(char *filename);

/* ------------------------------------------------------------------------ */
/* Allow for architectures in which a short does not have the same
   size or order as presumed in the CGM binary encoding standard
   (a short presumed to be 2 bytes in big-endian order, i.e.- MSB first).  */

#define CGM_WORD_ORDER(words, nwords) if(cgm_not)cgm_swap(words,nwords)

static void cgm_swap(short *words, long nwords);
static int cgm_not = -1;

static void
cgm_swap(short *words, long nwords)
{
  short word;
  char *pos;
  if (cgm_not<0) {
    word = 1;
    pos = (char *)&word;
    cgm_not = (sizeof(short)!=2 || pos[0]);
    if (!cgm_not) return;
  }
  pos = (char *)words;   /* reordering is done in place */
  while (nwords--) {
    word = *(words++);
    *(pos++) = word >> 8;
    *(pos++) = word & 0xff;
  }
}

/* ------------------------------------------------------------------------ */

static Octet *FormCommand(Octet *buffer, int klass, int id, long nbytes,
                          long *lpart)
{
  long nhere;
  *buffer++= (Octet)((klass<<4) | (id>>3));
  if (nbytes < 0x1f) {
    *buffer++= (Octet)((id<<5) | nbytes);
    nhere= nbytes;
  } else {
    *buffer++= (Octet)((id<<5) | (nbytes<0x1f? nbytes : 0x1f));
    if (nbytes < MAX_PARTITION) {
      nhere= nbytes;
      *buffer++= (Octet)(nhere>>8);
      *buffer++= (Octet)(nhere&0xff);
    } else {
      nhere= MAX_PARTITION;
      *buffer++= (Octet)((nhere>>8) | 0x80);
      *buffer++= (Octet)(nhere&0xff);
    }
  }
  *lpart= nhere;
  return buffer;
}

static Octet *NextPartition(Octet *buffer, long nbytes, long *lpart)
{
  long nhere;
  if (nbytes < MAX_PARTITION) {
    nhere= nbytes;
    *buffer++= (Octet)(nhere>>8);
    *buffer++= (Octet)(nhere&0xff);
  } else {
    nhere= MAX_PARTITION;
    *buffer++= (Octet)((nhere>>8) | 0x80);
    *buffer++= (Octet)(nhere&0xff);
  }
  *lpart= nhere;
  return buffer;
}

static Octet *Pascalify(Octet *buffer, const char *text, long len, int pad)
{
  /* Note-- strlen(text)<255 for this to work properly; writes
     strlen(text)+1 bytes into buffer (+2 if pad set and len even) */
  len= len>=0? len : (text? strlen(text) : 0);
  if (len>254) len= 254;     /* ignore text beyond 254 characters */
  pad= (pad && !(len&1));
  *buffer++= (Octet)len;
  while (len--) *buffer++= *text++;
  if (pad) *buffer++= '\0';  /* force total number of Octets to be even */
  return buffer;
}

static Octet *FormWords(Octet *buffer, const short *words, long nwords)
{
  while (nwords--) {
    *buffer++= (Octet)(*words >> 8);
    *buffer++= (Octet)(*words++ & 0xff);
  }
  return buffer;
}

static void FormReal(short *fixed, GpReal x)
{
  unsigned short *fx= (unsigned short *)fixed;
  GpReal ix= floor(x);
  fx[0]= (unsigned short)((int)ix);
  fx[1]= (unsigned short)(65536.0*(x-ix));
}

static int WriteB(p_file *file, void *buffer, long nbytes)
{
  if (!file) return 1;
  return p_fwrite(file, buffer, nbytes) != (unsigned long)nbytes;
}

static void WriteError(CGMEngine *cgm, const char *msg)
{
  if (!cgm->file) return;
  strcpy(gistError, msg);
  p_fclose(cgm->file);
  cgm->file= 0;
  cgm->state= 1;
}

/* ------------------------------------------------------------------------ */

/* Gist font numbers match those used by the Apple LaserWriter (for no
   particular reason).  The GPLOT program (from Pittsburgh supercomputer
   center) does a poor job with fonts, as does the CGTODL program on
   LLNL Crays.  The most important thing here is to make the Gist fonts
   used by the standard Gist style sheets (e.g. work.gs), Helvetica and
   Courier, be legible across these dysfunctional CGM interpreters.  */
extern int cgmFontNumbers[N_CGMFONTS];
int cgmFontNumbers[N_CGMFONTS]= {
  1,5,9,13,  3,7,11,15,  2,6,10,14,  4,8,12,16,  17,18,19,20 };

/* Use standard PostScript font names for now.  */
extern char *cgmFontNames[N_CGMFONTS];
char *cgmFontNames[N_CGMFONTS]= {
  "Courier", "Helvetica", "Times-Roman", "Symbol",
  "Courier-Bold", "Helvetica-Bold", "Times-Bold", "Symbol",
  "Courier-Oblique", "Helvetica-Oblique", "Times-Italic", "Symbol",
  "Courier-BoldOblique", "Helvetica-BoldOblique",
                                      "Times-BoldItalic", "Symbol",
  "NewCenturySchlbk-Roman", "NewCenturySchlbk-Bold",
  "NewCenturySchlbk-Italic", "NewCenturySchlbk-BoldItalic" };

static int BeginMetafile(CGMEngine *cgm)
{
  p_file *file;
  int i;
  long lcmnd, lpart, len, lfont[N_CGMFONTS], lfonts;
  char description[88];
  Octet *descriptor, *now, *font;
  short value[3];
  time_t t= time((time_t *)0);
  /* ctime returns 26 chars, e.g.- "Sun Jan  3 15:14:13  1988\n" */
  char *st= (t==-1)? "\n" : ctime(&t);

  if (!cgm || cgm->state) return 1;

  lcmnd= cgm->e.name? strlen(cgm->e.name) : 0;
  if (lcmnd>254) lcmnd= 254;
  cgm->e.name[lcmnd]= '\0';
  strcpy(description, "Gist;  ");
  strcpy(description+7, st? st : "\n");
  strcpy(description+31, ";  For: ");
  strncat(description, p_getuser(), 50L);
  lpart= strlen(description);
  lfonts= 0;
  for (i=0 ; i<N_CGMFONTS ; i++)
    lfonts+= lfont[i]= strlen(cgmFontNames[i]);
  lfonts+= N_CGMFONTS;  /* allow for Pascalification... */

  len= (lcmnd<31? 2 : 4) + lcmnd+1;
  if (len&1) len++;
  len+= 4;
  len+= 4 + lpart+1;
  if (len&1) len++;
  len+= 6*4;
  len+= 8;
  len+= 784;
  len+= 4 + lfonts;
  if (len&1) len++;
  descriptor= (Octet *)p_malloc(len+2);
  if (!descriptor) {
    strcpy(gistError, "memory manager failed in BeginMetafile");
    cgm->state= 1;     /* Metafile closed */
    return 1;
  }

  file= CREATE_CGM(cgm->filename);
  if (!file) {
    strcpy(gistError, "unable to create CGM output");
    cgm->state= 1;     /* Metafile closed */
    p_free(descriptor);
    return 1;
  }

  now= descriptor;

  now= FormCommand(now, 0, 1, lcmnd+1, &len);          /* BEGIN METAFILE */
  now= Pascalify(now, cgm->e.name, lcmnd, 1);

  now= FormCommand(now, 1, 1, 2L, &len);             /* METAFILE VERSION */
  value[0]= 1;
  now= FormWords(now, value, 1L);

  now= FormCommand(now, 1, 2, lpart+1, &len);    /* METAFILE DESCRIPTION */
  now= Pascalify(now, description, lpart, 1);

  now= FormCommand(now, 1, 3, 2L, &len);                     /* VDC TYPE */
  value[0]= 0;
  now= FormWords(now, value, 1L);

  now= FormCommand(now, 1, 4, 2L, &len);            /* INTEGER PRECISION */
  value[0]= 16;
  now= FormWords(now, value, 1L);

  now= FormCommand(now, 1, 6, 2L, &len);              /* INDEX PRECISION */
  value[0]= 16;
  now= FormWords(now, value, 1L);

  now= FormCommand(now, 1, 7, 2L, &len);              /* COLOR PRECISION */
  value[0]= 8;
  now= FormWords(now, value, 1L);

  now= FormCommand(now, 1, 8, 2L, &len);        /* COLOR INDEX PRECISION */
  value[0]= 8;
  now= FormWords(now, value, 1L);

  now= FormCommand(now, 1, 9, 1L, &len);          /* MAXIMUM COLOR INDEX */
  *now++= 255;
  *now++= 0;

  now= FormCommand(now, 1, 11, 6L, &len);       /* METAFILE ELEMENT LIST */
  value[0]= 1;
  value[1]= -1;
  value[2]= 1;    /* drawing plus control set */
  now= FormWords(now, value, 3L);

  now= FormCommand(now, 1, 12, 780L, &len);
                                        /* METAFILE DEFAULTS REPLACEMENT */
  now= FormCommand(now, 3, 6, 2L, &len);          /* CLIP INDICATOR */
  value[0]= 0;  /* clip off by default */
  now= FormWords(now, value, 1L);
  now= FormCommand(now, 5, 22, 2L, &len);         /* INTERIOR STYLE */
  value[0]= 1;  /* solid */
  now= FormWords(now, value, 1L);
  now= FormCommand(now, 5, 11, 2L, &len);         /* TEXT PRECISION */
  value[0]= 2;  /* stroke */
  now= FormWords(now, value, 1L);
  /* CGM color 0 is (by default) the background color, while color 1 is
     the foreground color.  Define the 8 standard Gist colors as
     colors 2-9.  This leaves colors 10-255 up for grabs.  The index
     into the Gist palette will therefore be incremented by 10 before
     being written into a CGM file.  Initialize these colors to a
     240 color gray scale, so that even without dumping the color table,
     we get sensible results.  */
  now= FormCommand(now, 5, 34, 763L, &len);          /* COLOR TABLE */
  *now++= 2;   /* begin at color index 2 (BLACK_COLOR, which is -3) */
  now[0]= now[1]= now[2]= 0;       /* black */
  now+= 3;
  now[0]= now[1]= now[2]= 255;     /* white */
  now+= 3;
  now[0]= 255; now[1]= now[2]= 0;  /* red */
  now+= 3;
  now[1]= 255; now[0]= now[2]= 0;  /* green */
  now+= 3;
  now[2]= 255; now[0]= now[1]= 0;  /* blue */
  now+= 3;
  now[0]= 0; now[1]= now[2]= 255;  /* cyan */
  now+= 3;
  now[1]= 0; now[0]= now[2]= 255;  /* magenta */
  now+= 3;
  now[2]= 0; now[0]= now[1]= 255;  /* yellow */
  now+= 3;
  for (i=0 ; i<240 ; i++) {
    now[0]= now[1]= now[2]= i + ((i+8)/15);  /* 240 grays */
    now+= 3;
  }
  for (i=0 ; i<6 ; i++) {
    now[0]= now[1]= now[2]= 255;   /* extra whites at top */
    now+= 3;
  }
  *now++= 0;  /* Fill up to even number of bytes */

  now= FormCommand(now, 1, 13, lfonts, &len);             /* FONT LIST */
  font= now;
  for (i=0 ; i<N_CGMFONTS ; i++)
    now= Pascalify(now, cgmFontNames[i], lfont[i], 0);
  if ((now-font)&1) *now++= 0;  /* even number of bytes */

  cgm->file= file;
  if (WriteB(file, descriptor, now-descriptor)) {
    WriteError(cgm, "write to CGM failed in BeginMetafile");
    p_free(descriptor);
    return 1;
  }
  p_free(descriptor);

  cgm->state= 2;  /* Metafile Description (actually finished) */
  return 0;
}

/* ------------------------------------------------------------------------ */

static void SetPageDefaults(CGMEngine *cgmEngine)
{
  int i;
  /* Set current state to match state set by GI procedure in ps.ps */
  cgmEngine->curClip= 0;
  cgmEngine->clipBox.xmin= cgmEngine->clipBox.xmax=
    cgmEngine->clipBox.ymin= cgmEngine->clipBox.ymax= 0.0;
  for (i=0 ; i<5 ; i++)
    cgmEngine->curColor[i]= FG_COLOR; /* CGM color index #1 */
  cgmEngine->curType= L_SOLID;
  cgmEngine->curWidth= 1.0;
  cgmEngine->curMark= M_ASTERISK;
  cgmEngine->curSize= 1.0;
  cgmEngine->curFont= T_TIMES;   /* CGM font #1 */
  cgmEngine->curHeight= 0.0;     /* Default is unusable (1 VDC) */
  cgmEngine->curAlignH= TH_NORMAL;
  cgmEngine->curAlignV= TV_NORMAL;
  cgmEngine->curPath= TX_RIGHT;
  cgmEngine->curOpaque= 0;
  cgmEngine->curEtype= L_NONE;
  cgmEngine->curEwidth= 1.0;
}

static void SetCGMTransform(GpTransform *toPixels, int landscape,
                            GpReal scale)
{
  toPixels->viewport.xmin= toPixels->viewport.ymin= 0.0;
  if (landscape) {
    toPixels->viewport.xmax= 1.033461;   /* 11 in */
    toPixels->viewport.ymax= 0.798584;   /* 8.5 in */
  } else {
    toPixels->viewport.xmax= 0.798584;   /* 8.5 in */
    toPixels->viewport.ymax= 1.033461;   /* 11 in */
  }
  toPixels->window.xmin= toPixels->window.ymin= 0.0;
  toPixels->window.xmax= toPixels->viewport.xmax*scale;
  toPixels->window.ymax= toPixels->viewport.ymax*scale;
}

static int BeginPage(CGMEngine *cgmEngine)
{
  p_file *file= cgmEngine->file;
  long len, lpage;
  Octet command[44], *now;
  char page[30];
  short xy[4];

  if (!cgmEngine) return 1;
  if (!file) {
    BeginMetafile(cgmEngine);
    file= cgmEngine->file;
    if (!file) return 1;
  }

  if (cgmEngine->state!=2 && cgmEngine->state!=5) {
    WriteError(cgmEngine, "CGM driver bug found in BeginPage");
    return 1;
  }

  /* A change in color table can take place only at the beginning
     of a page.  ChangePalette strobes the palette in the Engine base
     class part into the cgmEngine palette.  */
  cgmEngine->nColors= 0;      /* reset to mono mode */
  ChangePalette((Engine *)cgmEngine);

  /* Set transform viewport to reflect current page orientation */
  if (cgmEngine->landscape != cgmEngine->e.landscape) {
    SetCGMTransform(&cgmEngine->e.transform, cgmEngine->e.landscape,
                    cgmEngine->scale);
    cgmEngine->landscape= cgmEngine->e.landscape;
  }

  now= command;

  sprintf(page, "Page %d", cgmEngine->currentPage);
  lpage= strlen(page);
  now= FormCommand(now, 0, 3, lpage+1, &len);         /* BEGIN PICTURE */
  now= Pascalify(now, page, lpage, 1);

  now= FormCommand(now, 2, 6, 8L, &len);                 /* VDC EXTENT */
  xy[0]= (short)cgmEngine->e.transform.window.xmin;
  xy[1]= (short)cgmEngine->e.transform.window.ymin;
  xy[2]= (short)cgmEngine->e.transform.window.xmax;
  xy[3]= (short)cgmEngine->e.transform.window.ymax;
  now= FormWords(now, xy, 4L);

  now= FormCommand(now, 0, 4, 0L, &len);         /* BEGIN PICTURE BODY */
  if (WriteB(file, command, now-command)) {
    WriteError(cgmEngine, "write to CGM failed in BeginPage");
    return 1;
  }

  if (cgmEngine->e.colorMode && cgmEngine->e.palette &&
      cgmEngine->e.nColors>0) {
    int i, nColors= cgmEngine->e.nColors;
    GpColorCell *palette= cgmEngine->e.palette;
    Octet *col;
    if (nColors>246) nColors= 246;

    col= (Octet *)p_malloc(3*nColors+6);
    if (!col) {
      WriteError(cgmEngine, "memory manager failed in CGM BeginPage");
      return 1;
    }
    now= col;

    now= FormCommand(now, 5, 34, 3*nColors+1, &len);    /* COLOR TABLE */
    *now++= 10;   /* begin at color index 10 */
    for (i=0 ; i<nColors ; i++) {
      *now++= (Octet)(P_R(palette[i]));
      *now++= (Octet)(P_G(palette[i]));
      *now++= (Octet)(P_B(palette[i]));
    }
    if (!(nColors&1)) *now++= 0;

    if (WriteB(file, col, now-col)) {
      p_free(col);
      WriteError(cgmEngine, "write CT to CGM failed in BeginPage");
      return 1;
    }
    p_free(col);
    now= command;

    cgmEngine->colorMode= 1;  /* color table has been written */
    cgmEngine->nColors= nColors;

  } else {
    cgmEngine->colorMode= 0;  /* NO color table exists on this page */
    /* But, if there is a palette, still want to use grays */
    if (cgmEngine->e.palette && cgmEngine->e.nColors>0)
      cgmEngine->nColors= cgmEngine->e.nColors;
  }

  cgmEngine->state= 4;   /* Picture open */
  cgmEngine->e.marked= 1;

  return 0;
}

static void EndClip(CGMEngine *cgmEngine)
{
  if (!cgmEngine || cgmEngine->state!=4) return;
  if (cgmEngine->curClip) {
    p_file *file= cgmEngine->file;
    long len;
    Octet command[4];
    short clip= 0;
    FormCommand(command, 3, 6, 2L, &len);            /* CLIP INDICATOR */
    FormWords(command+2, &clip, 1L);
    if (WriteB(file, command, 4L))
      WriteError(cgmEngine, "write to CGM failed in EndClip");
    cgmEngine->curClip= 0;
  }
}

static void BeginClip(CGMEngine *cgmEngine, GpTransform *trans)
{
  GpPoint *points;
  GpBox *port= &trans->viewport;
  GpBox *box= &cgmEngine->clipBox;
  Octet command[16], *now;
  p_file *file;
  long len;

  if (!cgmEngine || cgmEngine->state!=4) return;
  file= cgmEngine->file;

  now= command;

  if (port->xmin!=box->xmin || port->xmax!=box->xmax ||
      port->ymin!=box->ymin || port->ymax!=box->ymax) {
    GpReal x[2], y[2];
    x[0]= trans->window.xmin;  x[1]= trans->window.xmax;
    y[0]= trans->window.ymin;  y[1]= trans->window.ymax;

    GpIntPoints(&cgmEngine->e.map, 3, 2, x, y, &points);
    if (points[0].x > points[1].x) {
      GpReal xt= points[0].x;
      points[0].x= points[1].x;
      points[1].x= (short)xt;
    }
    if (points[0].y > points[1].y) {
      GpReal yt= points[0].y;
      points[0].y= points[1].y;
      points[1].y= (short)yt;
    }

    now= FormCommand(now, 3, 5, 8L, &len);           /* CLIP RECTANGLE */
    /* Note-- this cast from GpPoint* to short* works on all machines
       that I know about, but may fail someday... */
    now= FormWords(now, (short *)points, 4L);
    *box= *port;
  }

  if (!cgmEngine->curClip) {
    short clip= 1;
    now= FormCommand(now, 3, 6, 2L, &len);           /* CLIP INDICATOR */
    now= FormWords(now, &clip, 1L);
    cgmEngine->curClip= 1;
  } else if (now==command) {
    return;
  }

  if (WriteB(file, command, now-command))
    WriteError(cgmEngine, "write to CGM failed in BeginClip");
}

static int EndPage(CGMEngine *cgmEngine)
{
  p_file *file;
  long len;
  Octet command[4];
  unsigned long pos = 0;
  int oops = 0;

  if (!cgmEngine) return 1;
  if (cgmEngine->state!=4) {
    if (cgmEngine->state==2 || cgmEngine->state==5) return 0;
    WriteError(cgmEngine, "CGM driver bug found in EndPage");
    return 1;
  }

  file= cgmEngine->file;
  FormCommand(command, 0, 5, 0L, &len);                 /* END PICTURE */
  FormCommand(command+2, 0, 2, 0L, &len);              /* END METAFILE */

  /* Write END PICTURE END METAFILE, then backspace over
     END METAFILE.  Combined with the p_fflush, this assures
     that a legal metafile is present on disk after EndPage.  */
  oops = WriteB(file, command, 4L);
  if (!oops) pos = p_ftell(file), oops = (pos==-1L || p_fseek(file, pos-2));
  if (oops) {
    WriteError(cgmEngine, "write to CGM failed in EndPage");
    return 1;
  }

  if (p_ftell(file) < (unsigned long)cgmEngine->fileSize) {
    /* Just flush buffers if file still reasonable size */
    p_fflush(file);
    cgmEngine->state= 5;   /* Picture closed */
  } else {
    /* File is too big, start a new one */
    p_fclose(file);
    cgmEngine->file= 0;
    cgmEngine->IncrementName(cgmEngine->filename);
    cgmEngine->state= 0;   /* No state yet */
  }

  cgmEngine->e.marked= 0;
  cgmEngine->currentPage++;
  SetPageDefaults(cgmEngine);
  return 0;
}

/* ------------------------------------------------------------------------ */

static int ChangePalette(Engine *engine)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  int nColors= engine->nColors;

  if (nColors<=0 || !engine->palette) {
    /* new color table is void-- don't allow indexed colors */
    cgmEngine->colorMode= 0;
    cgmEngine->nColors= 0;
    if (nColors>246) engine->nColors= 246;

  } else {
    /* remember current color mode, then set to mono mode until
       new color table can be written in BeginPage */
    cgmEngine->colorMode= 0;
    cgmEngine->nColors= 0;    /* don't index into table before BeginPage */
  }

  engine->colorChange= 0;
  return 246;
}

/* ------------------------------------------------------------------------ */

static void Kill(Engine *engine)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;

  if (engine->marked) EndPage(cgmEngine);

  if (cgmEngine->file) {
    /* Write directory of page locations */

    p_fclose(cgmEngine->file);
  }
  GpDelEngine(engine);
}

static int Clear(Engine *engine, int always)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  if (always || engine->marked) EndPage(cgmEngine);
  engine->marked= 0;
  return 0;
}

static int Flush(Engine *engine)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  if (cgmEngine->file) p_fflush(cgmEngine->file);
  return 0;
}

/* ------------------------------------------------------------------------ */

static int cgmColorID[]= { 4, 8, 14, 23, 29 };

static int SetupColor(CGMEngine *cgmEngine, unsigned long color, int which)
{
  int nColors= cgmEngine->nColors;
  p_file *file= cgmEngine->file;
  long len;
  short c = 0;
  Octet command[4];

  if (!cgmEngine->e.marked && BeginPage(cgmEngine)) return 1;
  if (color==cgmEngine->curColor[which]) return 0;

  if (color<240UL) {
    /* "CI index C"-- "CI" omitted if current color is indexed */
    GpColorCell *palette= cgmEngine->e.palette;
    if (nColors>0) {
      if (color>=(unsigned long)nColors) color= nColors-1;
      if (cgmEngine->colorMode) c= (short)color;  /* page has color table */
      else {                              /* this is a 240 gray mono page */
        c= (short)((P_R(palette[color]) +
                    P_G(palette[color]) + P_B(palette[color]))/3);
        c-= (c+8)/16;
      }
    } else {
      if (color>255) color= 255;
      c= (short)color;         /* No palette ==> no color table on page */
    }
    c+= 10;

  } else if (color<256UL) {
    /* standard color command FG, RED, etc. */
    if (color<YELLOW_COLOR) color= FG_COLOR;
    c= (short)(255UL-color);  /* standard colors are 0 thru 9 */
  } else {
    c = (short)(255UL-FG_COLOR);
  }

  FormCommand(command, 5, cgmColorID[which], 1L, &len);
        /* LINE COLOR (0) MARKER COLOR (1) TEXT COLOR (2) FILL COLOR (3) */
  command[2]= (Octet)c;
  command[3]= 0;
  if (WriteB(file, command, 4L)) {
    WriteError(cgmEngine, "write to CGM failed in SetupColor");
    return 1;
  }

  cgmEngine->curColor[which]= color;
  return 0;
}

static int SetupLine(CGMEngine *cgmEngine)
{
  p_file *file= cgmEngine->file;
  long len;
  Octet command[14], *now;

  if (cgmEngine->state!=4) return 1;
  if (SetupColor(cgmEngine, gistA.l.color, 0)) return 1;
  now= command;

  if (cgmEngine->curType!=gistA.l.type) {
    short ltype= (short)gistA.l.type;
    if (ltype==L_NONE) return 1;

    if (ltype<=0 || ltype>L_DASHDOTDOT) ltype= L_SOLID;

    now= FormCommand(now, 5, 2, 2L, &len);                  /* LINE TYPE */
    now= FormWords(now, &ltype, 1L);

    cgmEngine->curType= gistA.l.type;
  }

  if (cgmEngine->curWidth!=gistA.l.width) {
    short lwidth[2];

    now= FormCommand(now, 5, 3, 4L, &len);                 /* LINE WIDTH */
    FormReal(lwidth, gistA.l.width);
    now= FormWords(now, lwidth, 2L);

    cgmEngine->curWidth= gistA.l.width;
  }

  if (cgmEngine->curOpaque!=0 && cgmEngine->curType!=L_SOLID) {
    short trans= 0;

    now= FormCommand(now, 3, 4, 2L, &len);               /* TRANSPARENCY */
    now= FormWords(now, &trans, 1L);

    cgmEngine->curOpaque= trans;
  }

  if (now!=command && WriteB(file, command, now-command)) {
    WriteError(cgmEngine, "write to CGM failed in SetupLine");
    return 1;
  }

  return 0;
}

static int SetupEdge(CGMEngine *cgmEngine)
{
  p_file *file= cgmEngine->file;
  long len;
  Octet command[14], *now;

  if (cgmEngine->state!=4) return 1;
  if (SetupColor(cgmEngine, gistA.e.color, 4)) return 1;
  now= command;

  if (cgmEngine->curEtype!=gistA.e.type) {
    short ltype= (short)gistA.e.type;
    short visible= (ltype!=L_NONE);
    int change= visible^(cgmEngine->curEtype!=L_NONE);

    if (visible) {
      if (ltype<=0 || ltype>L_DASHDOTDOT) ltype= L_SOLID;
      now= FormCommand(now, 5, 27, 2L, &len);               /* EDGE TYPE */
      now= FormWords(now, &ltype, 1L);
    }
    if (change) {
      now= FormCommand(now, 5, 30, 2L, &len);         /* EDGE VISIBILITY */
      now= FormWords(now, &visible, 1L);
    }

    cgmEngine->curEtype= gistA.e.type;
  }

  if (cgmEngine->curEwidth!=gistA.e.width) {
    short lwidth[2];

    now= FormCommand(now, 5, 28, 4L, &len);                /* EDGE WIDTH */
    FormReal(lwidth, gistA.e.width);
    now= FormWords(now, lwidth, 2L);

    cgmEngine->curEwidth= gistA.e.width;
  }

  if (cgmEngine->curOpaque!=0 && cgmEngine->curEtype!=L_SOLID) {
    short trans= 0;

    now= FormCommand(now, 3, 4, 2L, &len);               /* TRANSPARENCY */
    now= FormWords(now, &trans, 1L);

    cgmEngine->curOpaque= trans;
  }

  if (now!=command && WriteB(file, command, now-command)) {
    WriteError(cgmEngine, "write to CGM failed in SetupEdge");
    return 1;
  }

  return 0;
}

static void CheckClip(CGMEngine *cgmEngine)
{
  if (!cgmEngine->e.marked) BeginPage(cgmEngine);
  if (gistClip) BeginClip(cgmEngine, &gistT);
  else if (cgmEngine->curClip) EndClip(cgmEngine);
}

static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  p_file *file;
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints;
  int firstPass= 1;
  GpPoint firstPoint, *points;
  long len;
  Octet command[10], *now;

  CheckClip(cgmEngine);
  if (n<1) return 0;
  if (SetupLine(cgmEngine)) return 1;
  file= cgmEngine->file;

  while ((nPoints=
          GpIntPoints(map, maxPoints, n, px, py, &points))) {
    if (closed) {
      if (firstPass) {
        firstPoint= points[0];
        firstPass= 0;
      }
      if (n==nPoints) {
        n++;
        points[nPoints++]= firstPoint;
      }
    }

    if (smooth) {
      /* There may be a GDP for Bezier curves, but I don't know how
         standard it is.  This is sure to be non-portable, but also
         sure to produces SOME acceptable output (identical, in fact,
         to what the Gist X engine displays).  */
      short value= smooth;  /* identifier */
      now= FormCommand(command, 7, 2, 3L, &len);     /* APPLICATION DATA */
      now= FormWords(now, &value, 1L);
      now= Pascalify(now, "", 0, 1);
    } else {
      now= command;
    }
    now= FormCommand(now, 4, 1, 4*nPoints, &len);            /* POLYLINE */
    if (WriteB(file, command, now-command)) {
      WriteError(cgmEngine, "write to CGM failed in DrawLines");
      return 1;
    }
    /* The cast from GpPoint* to short* works on all machines I know of,
       but may fail someday...  Note that maxPoints is set low enough
       that multiple partitions will never be required. */
    CGM_WORD_ORDER((short *)points, nPoints<<1);
    if (WriteB(file, (Octet *)points, nPoints<<2)) {
      WriteError(cgmEngine, "write to CGM failed in DrawLines");
      return 1;
    }

    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int SetupMarker(CGMEngine *cgmEngine)
{
  p_file *file= cgmEngine->file;
  long len;
  Octet command[14], *now;

  if (cgmEngine->state!=4) return 1;
  if (SetupColor(cgmEngine, gistA.m.color, 1)) return 1;
  now= command;

  if (cgmEngine->curMark!=gistA.m.type) {
    short mtype= (short)gistA.m.type;

    now= FormCommand(now, 5, 6, 2L, &len);                /* MARKER TYPE */
    now= FormWords(now, &mtype, 1L);

    cgmEngine->curMark= gistA.m.type;
  }

  if (cgmEngine->curSize!=gistA.m.size) {
    short msize[2];

    now= FormCommand(now, 5, 7, 4L, &len);                /* MARKER SIZE */
    FormReal(msize, gistA.m.size);
    now= FormWords(now, msize, 2L);

    cgmEngine->curSize= gistA.m.size;
  }

  if (cgmEngine->curOpaque!=0) {
    short trans= 0;

    now= FormCommand(now, 3, 4, 2L, &len);               /* TRANSPARENCY */
    now= FormWords(now, &trans, 1L);

    cgmEngine->curOpaque= trans;
  }

  if (now!=command && WriteB(file, command, now-command)) {
    WriteError(cgmEngine, "write to CGM failed in SetupMarker");
    return 1;
  }

  return 0;
}

static int DrawMarkers(Engine *engine, long n, const GpReal *px,
                       const GpReal *py)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  p_file *file;
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints;
  GpPoint *points;
  long len;
  Octet command[4], *now;

  if (n<1 || gistA.m.type<=0) return 0;
  CheckClip(cgmEngine);
  if (SetupMarker(cgmEngine)) return 1;
  file= cgmEngine->file;

  while ((nPoints=
          GpIntPoints(map, maxPoints, n, px, py, &points))) {

    now= FormCommand(command, 4, 3, 4*nPoints, &len);      /* POLYMARKER */
    if (WriteB(file, command, now-command)) {
      WriteError(cgmEngine, "write to CGM failed in DrawMarkers");
      return 1;
    }
    /* The cast from GpPoint* to short* works on all machines I know of,
       but may fail someday...  Note that maxPoints is set low enough
       that multiple partitions will never be required. */
    CGM_WORD_ORDER((short *)points, nPoints<<1);
    if (WriteB(file, (Octet *)points, nPoints<<2)) {
      WriteError(cgmEngine, "write to CGM failed in DrawMarkers");
      return 1;
    }

    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int SetupText(CGMEngine *cgmEngine)
{
  p_file *file= cgmEngine->file;
  int opq;
  long len;
  Octet command[36], *now;
  int h= gistA.t.alignH, v= gistA.t.alignV;
  GtGetAlignment(&gistA.t, &h, &v);

  if (cgmEngine->state!=4) return 1;
  if (SetupColor(cgmEngine, gistA.t.color, 2)) return 1;
  now= command;

  if (cgmEngine->curFont!=gistA.t.font) {
    int fn= (gistA.t.font>=0 && gistA.t.font<N_CGMFONTS)?
      gistA.t.font : 0;
    short font= (short)cgmFontNumbers[fn];

    now= FormCommand(now, 5, 10, 2L, &len);           /* TEXT FONT INDEX */
    now= FormWords(now, &font, 1L);

    cgmEngine->curFont= gistA.t.font;
  }

  if (cgmEngine->curHeight!=gistA.t.height) {
    /* -----------------------WARNING----------------------------- */
    /* The CGM standard specifies the character height, that is, the
       distance from the baseline to the capline of a line of text.
       This is different from the standard typestters' point size,
       which is the minimum vertical distance between lines of text
       required for good legibility.  The Gist "height" text attribute
       follows the ordinary definition of point size, as in PostScript,
       so here it must be shrunk by an arbitrary factor (0.8) to
       allow for the messed up CGM definition.  */
    /* -----------------------WARNING----------------------------- */
    short ch= (short)(0.8*gistA.t.height*cgmEngine->scale+0.5);

    now= FormCommand(now, 5, 15, 2L, &len);          /* CHARACTER HEIGHT */
    now= FormWords(now, &ch, 1L);

    cgmEngine->curHeight= gistA.t.height;
  }

  if (cgmEngine->curAlignH!=h || cgmEngine->curAlignV!=v) {
    short align[6];

    now= FormCommand(now, 5, 18, 12L, &len);           /* TEXT ALIGNMENT */
    align[0]= (short)h;
    align[1]= (short)v;
    align[2]= align[3]= align[4]= align[5]= 0;
    now= FormWords(now, align, 6L);

    cgmEngine->curAlignH= h;
    cgmEngine->curAlignV= v;
  }

  if (cgmEngine->curPath!=gistA.t.orient) {
    short orient[4];

    now= FormCommand(now, 5, 16, 8L, &len);     /* CHARACTER ORIENTATION */
    if (gistA.t.orient==TX_RIGHT) {
      orient[0]= 0;
      orient[1]= 1;
      orient[2]= 1;
      orient[3]= 0;
    } else if (gistA.t.orient==TX_UP) {
      orient[0]= -1;
      orient[1]= 0;
      orient[2]= 0;
      orient[3]= 1;
    } else if (gistA.t.orient==TX_LEFT) {
      orient[0]= 0;
      orient[1]= -1;
      orient[2]= -1;
      orient[3]= 0;
    } else {
      orient[0]= 1;
      orient[1]= 0;
      orient[2]= 0;
      orient[3]= -1;
    }
    now= FormWords(now, orient, 4L);

    cgmEngine->curPath= gistA.t.orient;
  }

  opq= gistA.t.opaque;
  /* if (gistA.t.orient!=TX_RIGHT) opq= 1; let this be X only limitation */
  if (cgmEngine->curOpaque != (opq!=0)) {
    short trans= (short)(opq==0);

    now= FormCommand(now, 3, 4, 2L, &len);               /* TRANSPARENCY */
    now= FormWords(now, &trans, 1L);

    cgmEngine->curOpaque= !trans;
  }

  if (now!=command && WriteB(file, command, now-command)) {
    WriteError(cgmEngine, "write to CGM failed in SetupText");
    return 1;
  }

  return 0;
}

static Octet cgmText[264];  /* header + point + flag + 254 pascal chars */

static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  p_file *file;
  int nlines;
  GpReal width, height, lineHeight;
  GpXYMap *map= &engine->map;
  GpBox *wind= &engine->transform.window;
  GpReal xmin, xmax, ymin, ymax;
  int count, alignV /*, alignH*/;
  short xyf[3];
  Octet *now;
  long len;
  const char *t;

  CheckClip(cgmEngine);
  if (SetupText(cgmEngine)) return 1;
  file= cgmEngine->file;
  /*alignH= cgmEngine->curAlignH;*/
  alignV= cgmEngine->curAlignV;

  /* handle multi-line strings */
  nlines= GtTextShape(text, &gistA.t, (WidthFunction)0, &width);
  lineHeight= gistA.t.height * cgmEngine->scale;
  width*= 0.6*lineHeight;
  height= lineHeight*(GpReal)nlines;

  /* Compute text position in VDC coordinates */
  x0= map->x.scale*x0+map->x.offset;
  y0= map->y.scale*y0+map->y.offset;

  /* Reject if and only if the specified point is off of the current
     page by more than the approximate size of the text.  Note that
     the entire block of text is either accepted or rejected --
     the CGM interpreter does the clipping.  */
  if (wind->xmax>wind->xmin) { xmin= wind->xmin; xmax= wind->xmax; }
  else { xmin= wind->xmax; xmax= wind->xmin; }
  if (wind->ymax>wind->ymin) { ymin= wind->ymin; ymax= wind->ymax; }
  else { ymin= wind->ymax; ymax= wind->ymin; }
  if (gistA.t.orient==TX_RIGHT || gistA.t.orient==TX_LEFT) {
    if (x0<xmin-width || x0>xmax+width ||
        y0<ymin-height || y0>ymax+height) return 0;
  } else {
    if (x0<xmin-height || x0>xmax+height ||
        y0<ymin-width || y0>ymax+width) return 0;
  }

  /* Adjust y0 to represent topmost line */
  if (nlines > 1) {
    if (gistA.t.orient==TX_RIGHT) {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) y0+= height-lineHeight;
      if (alignV==TV_HALF) y0+= 0.5*(height-lineHeight);
    } else if (gistA.t.orient==TX_LEFT) {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) y0-= height-lineHeight;
      if (alignV==TV_HALF) y0-= 0.5*(height-lineHeight);
    } else if (gistA.t.orient==TX_UP) {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) x0-= height-lineHeight;
      if (alignV==TV_HALF) x0-= 0.5*(height-lineHeight);
    } else {
      if (alignV==TV_BASE || alignV==TV_BOTTOM) x0+= height-lineHeight;
      if (alignV==TV_HALF) x0+= 0.5*(height-lineHeight);
    }
  }

  xyf[2]= 1;  /* never use append text, all strings final */
  while ((t= GtNextLine(text, &count, TX_RIGHT))) {
    if (count>0) {
      now= FormCommand(cgmText, 4, 4, (long)count+7, &len);      /* TEXT */
      xyf[0]= (short)x0;
      xyf[1]= (short)y0;
      now= FormWords(now, xyf, 3L);
      now= Pascalify(now, t, (long)count, 1);
      if (WriteB(file, cgmText, now-cgmText)) {
        WriteError(cgmEngine, "write to CGM failed in DrwText");
        return 1;
      }
    }

    text= t+count;
    if (gistA.t.orient==TX_RIGHT) y0-= lineHeight;
    else if (gistA.t.orient==TX_UP) x0+= lineHeight;
    else if (gistA.t.orient==TX_LEFT) y0+= lineHeight;
    else x0-= lineHeight;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  p_file *file;
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints;
  GpPoint *points;
  int value= 0;
  long len;
  Octet command[4], *now;

  /* For now, only solid interior style supported */

  if (n<3) return 0;
  CheckClip(cgmEngine);
  if (SetupColor(cgmEngine, gistA.f.color, 3) ||
      SetupEdge(cgmEngine)) return 1;
  file= cgmEngine->file;

  while ((nPoints=
          GpIntPoints(map, maxPoints, n, px, py, &points))) {

    now= FormCommand(command, 4, 7, 4*nPoints, &len);         /* POLYGON */
    if (WriteB(file, command, now-command)) {
      WriteError(cgmEngine, "write to CGM failed in DrawFill");
      return 1;
    }
    /* The cast from GpPoint* to short* works on all machines I know of,
       but may fail someday...  Note that maxPoints is set low enough
       that multiple partitions will never be required. */
    CGM_WORD_ORDER((short *)points, nPoints<<1);
    if (WriteB(file, (Octet *)points, nPoints<<2)) {
      WriteError(cgmEngine, "write to CGM failed in DrawFill");
      return 1;
    }

    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
    value= 1;   /* Polygons with >4050 sides won't be filled correctly */
  }

  return value;
}

/* ------------------------------------------------------------------------ */

static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  p_file *file;
  GpXYMap *map= &cgmEngine->e.map;
  int nColors= cgmEngine->nColors;
  GpColorCell *palette;
  short pqrwhd[10];
  long i, j, off, nCells, lPart;
  int pad, color, colorMode;
  Octet *buffer, *now;

  if (!cgmEngine->e.marked && BeginPage(cgmEngine)) return 1;
  CheckClip(cgmEngine);
  file= cgmEngine->file;

  /* Transform corner coordinates, clipping and adjusting width,
     height, nColumns, and colors as necessary.  */
  width= GpClipCells(&map->x, &px, &qx,
                     gistT.window.xmin, gistT.window.xmax, width, &off);
  colors+= off;
  height= GpClipCells(&map->y, &py, &qy,
                      gistT.window.ymin, gistT.window.ymax, height, &off);
  colors+= nColumns*off;
  if (width>0x7ffe) width= 0x7ffe;
  if (height>0x7ffe) height= 0x7ffe;
  pad= ((width & 1)!=0);

  if (width<=0 || height<=0) return 0;

  nCells= (pad? width+1 : width)*height;
  lPart= nCells+20;
  if (lPart>MAX_PARTITION) lPart= MAX_PARTITION;
  if (lPart & 1) lPart++;
  buffer= (Octet *)p_malloc(lPart+4);
  if (!buffer) {
    WriteError(cgmEngine, "memory manager failed in CGM DrawCells");
    return 1;
  }

  now= buffer;
  now= FormCommand(now, 4, 9, lPart, &lPart);              /* CELL ARRAY */
  pqrwhd[0]= (short)px;
  pqrwhd[1]= (short)py;
  pqrwhd[2]= (short)qx;
  pqrwhd[3]= (short)qy;
  pqrwhd[4]= (short)qx;
  pqrwhd[5]= (short)py;
  pqrwhd[6]= (short)width;  /* Note limitation to < 32767-by-32767 */
  pqrwhd[7]= (short)height;
  pqrwhd[8]= 0;  /* default depth is color index precision */
  pqrwhd[9]= 1;  /* NOT run length encoded */
  now= FormWords(now, pqrwhd, 10L);

  if (nColors>0) {
    /* Image value will be either index or palette gray */
    colorMode= cgmEngine->colorMode;
    if (colorMode) {
      palette= 0;                      /* palette already written */
    } else {
      palette= cgmEngine->e.palette;   /* must lookup gray level now */
    }
  } else {
    /* Must assume image varies over maximum possible range */
    colorMode= 1;  /* That is, use index without trying palette lookup */
    palette= 0;
  }

  i= j= 0;
  lPart-= 20;  /* first partition has 20 octets of pqrwhd */
  for (;;) {
    while (lPart--) {
      if (i>=width) {
        if (pad) {
          *now++= 0;
          nCells--;
          if (!lPart--) break;
        }
        i= 0;
        j+= nColumns;
      }
      if (gistA.rgb) {
        color = (colors[3*(i+j)]+colors[3*(i+j)+1]+colors[3*(i+j)+2])/3;
        color = (color*nColors)>>8;
      } else {
        color= colors[i+j];
      }
      i++;
      if (color>=nColors && nColors>0) color= nColors-1;
      if (!colorMode) {  /* this is default 240 gray mono page */
        color= (P_R(palette[color]) +
                P_G(palette[color]) + P_B(palette[color]))/3;
        color-= (color+8)/16;
      }

      *now++= 10+color;  /* skip the 10 standard colors */
      nCells--;
    }

    if (WriteB(file, buffer, now-buffer)) {
      p_free(buffer);
      WriteError(cgmEngine, "write to CGM failed in DrawCells");
      return 1;
    }

    if (nCells<=0) break;

    now= NextPartition(buffer, nCells, &lPart);
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy)
{
  CGMEngine *cgmEngine= (CGMEngine *)engine;
  p_file *file;
  GpXYMap *map= &engine->map;
  long maxSegs= 2025, nSegs;
  GpSegment *segs;
  long len;
  Octet command[4], *now;

  CheckClip(cgmEngine);
  if (n<1) return 0;
  if (SetupLine(cgmEngine)) return 1;
  file= cgmEngine->file;

  while ((nSegs=
          GpIntSegs(map, maxSegs, n, px, py, qx, qy, &segs))) {

    now= FormCommand(command, 4, 2, 8*nSegs, &len); /* DISJOINT POLYLINE */
    if (WriteB(file, command, now-command)) {
      WriteError(cgmEngine, "write to CGM failed in DrawDisjoint");
      return 1;
    }
    /* The cast from GpSegment* to short* works on all machines I know of,
       but may fail someday...  Note that maxSegs is set low enough
       that multiple partitions will never be required. */
    CGM_WORD_ORDER((short *)segs, nSegs<<2);
    if (WriteB(file, (Octet *)segs, nSegs<<3)) {
      WriteError(cgmEngine, "write to CGM failed in DrawDisjoint");
      return 1;
    }

    if (n==nSegs) break;
    n-= nSegs;
    px+= nSegs;
    py+= nSegs;
    qx+= nSegs;
    qy+= nSegs;
  }

  return 0;
}

/* ------------------------------------------------------------------------ */

static void IncrementName(char *filename)
{
  int i, len= filename? strlen(filename) : 0;
  if (len>4 && (strcmp(filename+len-4, ".cgm")==0 ||
                strcmp(filename+len-4, ".CGM")==0)) i= len-4;
  else i= len;
  while (i-- > 0) {
    if (filename[i]=='9') {
      filename[i]= '0';
    } else {
      if (filename[i]=='Z' || filename[i]>='z') filename[i]= '0';
      else filename[i]++;
      break;
    }
  }
}

/* ------------------------------------------------------------------------ */

GpReal gCGMScale= 25545.24;   /* default CGM scale is 2400 dpi (round up) */
long gCGMFileSize= 10000000;  /* default max file size is about 10 Meg */

Engine *GpCGMEngine(char *name, int landscape, int mode, char *file)
{
  CGMEngine *cgmEngine;
  long flen= file? strlen(file) : 0;
  long engineSize= sizeof(CGMEngine)+flen+1;
  GpTransform toPixels;

  if (flen<=0) return 0;

  SetCGMTransform(&toPixels, landscape, gCGMScale);

  cgmEngine=
    (CGMEngine *)GpNewEngine(engineSize, name, cgmType, &toPixels, landscape,
                             &Kill, &Clear, &Flush, &GpComposeMap,
                             &ChangePalette, &DrawLines, &DrawMarkers,
                             &DrwText, &DrawFill, &DrawCells,
                             &DrawDisjoint);

  if (!cgmEngine) {
    strcpy(gistError, "memory manager failed in GpCGMEngine");
    return 0;
  }

  cgmEngine->filename= (char *)(cgmEngine+1);
  strcpy(cgmEngine->filename, file);
  cgmEngine->scale= gCGMScale;
  cgmEngine->fileSize= gCGMFileSize;
  cgmEngine->IncrementName= &IncrementName;
  cgmEngine->file= 0;
  cgmEngine->state= 0;

  SetPageDefaults(cgmEngine);
  cgmEngine->e.colorMode= mode;
  cgmEngine->colorMode= 0;
  cgmEngine->nColors= 0;

  cgmEngine->landscape= landscape;
  cgmEngine->currentPage= 1;

  return (Engine *)cgmEngine;
}

CGMEngine *GisCGMEngine(Engine *engine)
{
  return (engine && engine->type==cgmType)? (CGMEngine *)engine : 0;
}

void GcgmSetScale(CGMEngine *cgmEngine, GpReal scale)
{
  GpTransform toPixels;

  if (!cgmEngine || cgmEngine->state!=0) return;

  cgmEngine->scale= scale;
  SetCGMTransform(&toPixels, cgmEngine->landscape, scale);
}

/* ------------------------------------------------------------------------ */
