/*
 * EPS.C
 *
 * $Id$
 *
 * Define the Encapsulated PostScript pseudo-engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include <string.h>
#include <stdio.h>
#include "pstdio.h"
#include "pstdlib.h"

#include "eps.h"

#include "ps.h"
#include "gtext.h"

#define PS_PER_POINT 20
#define PS_TO_NDC (ONE_POINT/20.0)

static char epsType[]= "EPS Pseudo-engine";

static p_file *psFile= 0, *epsFile= 0;
static unsigned char *epsPreview= 0;
static char line[256];
static int xll, xur, yll, yur, epsLandscape;
static GpReal epsXScale, epsYScale;

static void Kill(Engine *engine);
static int Clear(Engine *engine, int always);
static int Flush(Engine *engine);
static int ChangePalette(Engine *engine);
static void Rasterize(GpPoint *p0, GpPoint *p1);
static void ScaleFont(GpReal ndcHeight);
static void RenderChars(const char *text, int nChars, int orient);
static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth);
static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text);
static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py);
static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors);
static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy);

/* ------------------------------------------------------------------------ */

/* FrameMaker 3.0X has a bug in its EPS preview display to X Windows,
   which causes the preview to appear inverted (although the file prints
   correctly).  Set epsFMbug to compensate for this.  */
int epsFMbug= 0;

static char hexDigits[]= "0123456789abcdef";

static void Kill(Engine *engine)
{
  int i, j;

  /* Add EPSF to first comment line */
  if (!p_fgets(psFile, line, 256)) goto err;
  p_fputs(epsFile, "%!PS-Adobe-3.0 EPSF-3.0\n");

  /* Copy other header comment lines, moving bounding box to front */
  do {
    if (!p_fgets(psFile, line, 256)) goto err;
    if (!strncmp("%%BoundingBox:", line, 14L))
      sprintf(line, "%%%%BoundingBox: %d %d %d %d\n", xll, yll, xur, yur);
    p_fputs(epsFile, line);
  } while (strncmp("%%EndComments", line, 13L));

  /* Write preview comments */
  p_fputs(epsFile, "%%BeginPreview: 256 256 1 256\n");
  line[0]= '%';
  line[1]= ' ';
  line[66]= '\n';
  line[67]= '\0';
  if (epsFMbug) {
    for (j=32*255 ; j>=0 ; j-= 32) {
      for (i=0 ; i<32 ; i++) {
        line[2*i+2]= hexDigits[epsPreview[j+i]>>4];
        line[2*i+3]= hexDigits[epsPreview[j+i]&0xf];
      }
      p_fputs(epsFile, line);
    }
  } else {
    for (j=0 ; j<32*256 ; j+= 32) {
      for (i=0 ; i<32 ; i++) {
        line[2*i+2]= hexDigits[epsPreview[j+i]>>4];
        line[2*i+3]= hexDigits[epsPreview[j+i]&0xf];
      }
      p_fputs(epsFile, line);
    }
  }
  p_fputs(epsFile, "%%EndPreview\n");

  /* Copy rest of PS file, omitting bounding box at end */
  for (;;) {
    if (!p_fgets(psFile, line, 256)) break;
    if (strncmp("%%BoundingBox:", line, 14L)) p_fputs(epsFile, line);
  }

 err:
  p_fclose(psFile);
  psFile= 0;
  p_remove("_tmp.eps");
  p_fclose(epsFile);
  epsFile= 0;
  p_free(epsPreview);
  epsPreview= 0;
  GpDelEngine(engine);
}

static int Clear(Engine *engine, int always)
{
  engine->marked= 0;
  return 0;
}

static int Flush(Engine *engine)
{
  return 0;
}

static int ChangePalette(Engine *engine)
{
  engine->colorChange= 0;
  return 256;
}

/* ------------------------------------------------------------------------ */

#define X_MAX 255
#define Y_MAX (255<<5)

static void Rasterize(GpPoint *p0, GpPoint *p1)
{
  int lr, bt;
  short p0x, p0y, p1x, p1y, dx, dy, x, y, last, test, dt;

  if (epsLandscape) {
    p0x= p0->y;   p0y= p0->x;
    p1x= p1->y;   p1y= p1->x;
  } else {
    p0x= p0->x;   p0y= p0->y;
    p1x= p1->x;   p1y= p1->y;
  }

  lr= p1x>=p0x;
  dx= lr? p1x-p0x : p0x-p1x;
  bt= p1y>=p0y;
  dy= bt? p1y-p0y : p0y-p1y;

  if (dx < dy) {
    /* Turn on 1 pixel for each y coordinate value (slope>1) */
    if (bt) {
      y= p0y<<5;
      last= p1y<<5;
      if (p1x<0 || p1x>X_MAX) return;
      x= p0x;
    } else {
      y= p1y<<5;
      last= p0y<<5;
      if (p0x<0 || p0x>X_MAX) return;
      x= p1x;
      lr= !lr;
    }
    if (x<0 || x>X_MAX || y<0 || last>Y_MAX) return;
    test= dx - (dy>>1);
    dt= dx - dy;
    epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
    if (lr) {
      for (y+=32 ; y<=last ; y+=32) {
        test+= dt;
        if (test>0) {
          x++;
        } else {
          test+= dy;
        }
        epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
      }
    } else {
      for (y+=32 ; y<=last ; y+=32) {
        test+= dt;
        if (test>0) {
          x--;
        } else {
          test+= dy;
        }
        epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
      }
    }

  } else {
    /* Turn on 1 pixel for each x coordinate value (slope<=1) */
    if (lr) {
      x= p0x;
      last= p1x;
      if (p1y<0 || (p1y<<5)>Y_MAX) return;
      y= p0y<<5;
    } else {
      x= p1x;
      last= p0x;
      if (p0y<0 || (p0y<<5)>Y_MAX) return;
      y= p1y<<5;
      bt= !bt;
    }
    if (x<0 || last>X_MAX || y<0 || y>Y_MAX) return;
    test= dy - (dx>>1);
    dt= dy - dx;
    epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
    if (bt) {
      for (x++ ; x<=last ; x++) {
        test+= dt;
        if (test>0) {
          y+= 32;
        } else {
          test+= dx;
        }
        epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
      }
    } else {
      for (x++ ; x<=last ; x++) {
        test+= dt;
        if (test>0) {
          y-= 32;
        } else {
          test+= dx;
        }
        epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
      }
    }
  }
}

/* ------------------------------------------------------------------------ */

/* Stroked characters are conform exactly to the "guesses" about
   text bounding boxes made elsewhere in Gist-- namely, if the
   point size is 10, the character width is 6, the bottom line
   is 2 below the baseline, and the top line is 8 above the
   baseline.  (The cap and half lines differ slightly from the
   assumptions of 8 above and 4 above the baseline...)  */

static GpReal xAlpha[4]= { 0.05, 0.30, 0.55, 0.60 };
static GpReal yAlpha[7]= { -0.18, 0.00, 0.18, 0.36, 0.54, 0.72, 0.80 };
static GpReal xChar[4], yChar[7], xCurrent, yCurrent;
static GpXYMap *charMap;

static void ScaleFont(GpReal ndcHeight)
{
  GpReal yScale= ndcHeight*epsYScale;
  GpReal xScale= ndcHeight*epsXScale;
  int i;
  for (i=0 ; i<4 ; i++) xChar[i]= xScale*xAlpha[i];
  for (i=0 ; i<7 ; i++) yChar[i]= yScale*yAlpha[i];
}

/* Characters are represented as sequences of 4+n short integers,
   where n is the number of points required to draw the character.
   The first 4 numbers are counts of the number of points to
   connect in a polyline, while the subsequent n numbers are
   indices into the xAlpha and yAlpha arrays as follows:

    capline   20   21   22
              16   17   18
              12   13   14
               8    9   10
   baseline    4    5    6
               0    1    2

   Only the 94 printing (non-blank) ASCII characters are represented,
   character codes 0x21 through 0x7e (33 through 126) inclusive.
 */

static short c01[]= { 2,2,0,0,  21,9, 5,5 };                 /* ! ASCII 33 */ 
static short c02[]= { 2,2,0,0,  21,17, 22,18 };
static short c03[]= { 2,2,2,2,  21,4, 22,5, 16,18, 8,10 };
static short c04[]= { 2,6,0,0,  21,5, 18,16,12,14,10,8 };
static short c05[]= { 2,4,4,0,  22,4, 20,12,17,20, 14,9,6,14 };
static short c06[]= { 8,0,0,0,  6,16,21,18,8,4,5,10 };
static short c07[]= { 2,0,0,0,  21,18 };
static short c08[]= { 4,0,0,0,  22,17,9,6 };
static short c09[]= { 4,0,0,0,  20,17,9,4 };
static short c10[]= { 2,2,2,0,  17,9, 18,8, 16,10 };
static short c11[]= { 2,2,0,0,  17,9, 12,14 };
static short c12[]= { 2,0,0,0,  5,0 };
static short c13[]= { 2,0,0,0,  12,14 };
static short c14[]= { 2,0,0,0,  5,5 };
static short c15[]= { 2,0,0,0,  22,4 };
static short c16[]= { 7,2,0,0,  21,16,8,5,10,18,21, 13,13 }; /* 0 ASCII 48 */
static short c17[]= { 3,2,0,0,  16,21,5, 4,6 };
static short c18[]= { 6,0,0,0,  16,21,18,8,4,6 };
static short c19[]= { 6,2,0,0,  20,21,18,10,5,4, 13,14 };
static short c20[]= { 4,0,0,0,  6,22,12,14 };
static short c21[]= { 7,0,0,0,  22,20,12,13,10,5,4 };
static short c22[]= { 8,0,0,0,  22,21,16,8,5,10,13,12 };
static short c23[]= { 3,0,0,0,  20,22,5 };
static short c24[]= { 7,0,0,0,  21,16,10,5,8,18,21 };
static short c25[]= { 8,0,0,0,  4,5,10,18,21,16,13,14 };
static short c26[]= { 2,2,0,0,  17,17, 9,9 };
static short c27[]= { 2,2,0,0,  9,9, 5,0 };
static short c28[]= { 3,0,0,0,  22,12,6 };
static short c29[]= { 2,2,0,0,  16,18, 8,10 };
static short c30[]= { 3,0,0,0,  20,14,4 };
static short c31[]= { 5,2,0,0,  16,21,18,13,9, 5,5 };
static short c32[]= { 9,0,0,0,  18,13,14,22,21,16,8,5,6 };
static short c33[]= { 5,2,0,0,  4,16,21,18,6, 12,14 };       /* A ASCII 65 */
static short c34[]= { 9,0,0,0,  13,18,21,20,4,5,10,13,12 };
static short c35[]= { 6,0,0,0,  22,21,16,8,5,6 };
static short c36[]= { 7,0,0,0,  20,4,5,10,18,21,20 };
static short c37[]= { 4,2,0,0,  22,20,4,6, 12,14 };
static short c38[]= { 3,2,0,0,  22,20,4, 12,14 };
static short c39[]= { 8,0,0,0,  22,21,16,8,5,6,14,13 };
static short c40[]= { 2,2,2,0,  20,4, 22,6, 12,14 };
static short c41[]= { 2,2,2,0,  21,5, 20,22, 4,6 };
static short c42[]= { 4,0,0,0,  8,5,10,22 };
static short c43[]= { 3,3,0,0,  20,12,6, 22,12,4 };
static short c44[]= { 3,0,0,0,  20,4,6 };
static short c45[]= { 5,0,0,0,  4,20,13,22,6 };
static short c46[]= { 4,0,0,0,  4,20,6,22 };
static short c47[]= { 7,0,0,0,  21,16,8,5,10,18,21 };
static short c48[]= { 6,0,0,0,  12,13,18,21,20,4 };
static short c49[]= { 7,2,0,0,  21,16,8,5,10,18,21, 9,6 };
static short c50[]= { 7,2,0,0,  12,13,18,21,20,4, 13,6 };
static short c51[]= { 6,0,0,0,  22,21,16,10,5,4 };
static short c52[]= { 2,2,0,0,  20,22, 21,5 };
static short c53[]= { 5,0,0,0,  20,8,5,10,22 };
static short c54[]= { 3,0,0,0,  20,5,22 };
static short c55[]= { 5,0,0,0,  20,4,13,6,22 };
static short c56[]= { 2,2,0,0,  20,6, 4,22 };
static short c57[]= { 3,2,0,0,  20,13,22, 13,5 };
static short c58[]= { 2,0,0,0,  20,22,4,6 };                 /* Z ASCII 90 */
static short c59[]= { 4,0,0,0,  22,21,5,6 };
static short c60[]= { 2,0,0,0,  20,6 };
static short c61[]= { 4,0,0,0,  20,21,5,4 };
static short c62[]= { 3,0,0,0,  12,21,14 };
static short c63[]= { 2,0,0,0,  4,6 };
static short c64[]= { 2,0,0,0,  21,16 };
static short c65[]= { 6,0,0,0,  14,6,5,8,13,10 };            /* a ASCII 97 */
static short c66[]= { 6,0,0,0,  20,4,5,10,13,12 };
static short c67[]= { 5,0,0,0,  14,13,8,5,6 };
static short c68[]= { 6,0,0,0,  14,13,8,5,6,22 };
static short c69[]= { 7,0,0,0,  8,10,14,13,8,5,6 };
static short c70[]= { 3,2,0,0,  18,21,5, 13,14 };
static short c71[]= { 7,0,0,0,  1,2,14,13,8,5,6 };
static short c72[]= { 2,4,0,0,  20,4, 12,13,10,6 };
static short c73[]= { 2,2,0,0,  17,17, 13,5 };
static short c74[]= { 2,3,0,0,  17,17, 13,1,0 };
static short c75[]= { 3,3,0,0,  20,8,6, 4,8,14 };
static short c76[]= { 2,0,0,0,  21,5 };
static short c77[]= { 4,3,0,0,  4,12,13,5, 9,10,2 };
static short c78[]= { 5,0,0,0,  4,12,13,10,6 };
static short c79[]= { 5,0,0,0,  13,8,5,10,13 };
static short c80[]= { 6,0,0,0,  0,12,13,10,5,4 };
static short c81[]= { 6,0,0,0,  2,14,13,8,5,6 };
static short c82[]= { 2,3,0,0,  4,12, 8,13,14 };
static short c83[]= { 4,0,0,0,  14,12,6,4 };
static short c84[]= { 2,2,0,0,  17,5, 12,14 };
static short c85[]= { 5,0,0,0,  12,8,5,10,14 };
static short c86[]= { 3,0,0,0,  12,5,14 };
static short c87[]= { 5,0,0,0,  12,4,9,6,14 };
static short c88[]= { 2,2,0,0,  12,6, 4,14 };
static short c89[]= { 2,3,0,0,  12,5, 0,5,14 };
static short c90[]= { 4,0,0,0,  12,14,4,6 };                /* z ASCII 122 */
static short c91[]= { 4,2,0,0,  22,21,5,6, 12,13 };
static short c92[]= { 2,2,2,0,  21,17, 13,13, 9,5 };
static short c93[]= { 4,2,0,0,  20,21,5,4, 13,14 };
static short c94[]= { 4,0,0,0,  16,21,13,18 };

static short *cN[94]= {
  c01, c02, c03, c04, c05, c06, c07, c08, c09, c10,
  c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
  c21, c22, c23, c24, c25, c26, c27, c28, c29, c30,
  c31, c32, c33, c34, c35, c36, c37, c38, c39, c40,
  c41, c42, c43, c44, c45, c46, c47, c48, c49, c50,
  c51, c52, c53, c54, c55, c56, c57, c58, c59, c60,
  c61, c62, c63, c64, c65, c66, c67, c68, c69, c70,
  c71, c72, c73, c74, c75, c76, c77, c78, c79, c80,
  c81, c82, c83, c84, c85, c86, c87, c88, c89, c90,
  c91, c92, c93, c94 };

static void RenderChars(const char *text, int nChars, int orient)
{
  GpPoint xy[10];
  unsigned char c;
  short i, j, n, *nstroke, *stroke;

  while (nChars-- > 0) {
    c= (unsigned char)(*text++);
    if (c>32 && c<127) {
      nstroke= cN[c-33];
      stroke= nstroke+4;
      for (i=0 ; (n=nstroke[i]) && i<4 ; i++) {
        for (j=0 ; j<n ; j++) {
          if (orient==TX_RIGHT) {
            xy[j].x= (short)(xCurrent+xChar[stroke[j]&3]);
            xy[j].y= (short)(yCurrent+yChar[stroke[j]>>2]);
          } else if (orient==TX_UP) {
            xy[j].y= (short)(yCurrent+xChar[stroke[j]&3]);
            xy[j].x= (short)(xCurrent-yChar[stroke[j]>>2]);
          } else if (orient==TX_LEFT) {
            xy[j].x= (short)(xCurrent-xChar[stroke[j]&3]);
            xy[j].y= (short)(yCurrent-yChar[stroke[j]>>2]);
          } else {
            xy[j].y= (short)(yCurrent-xChar[stroke[j]&3]);
            xy[j].x= (short)(xCurrent+yChar[stroke[j]>>2]);
          }
        }
        for (j=0 ; j<n-1 ; j++) Rasterize(xy+j, xy+j+1);
        stroke+= n;
      }
    }
    xCurrent+= xChar[3];
  }
}

/* ------------------------------------------------------------------------ */

static int DrawLines(Engine *engine, long n, const GpReal *px,
                     const GpReal *py, int closed, int smooth)
{
  GpXYMap *map= &engine->map;
  long maxPoints= 4050, nPoints, i;
  int firstPass= 1;
  GpPoint firstPoint, *points;

  if (n<1 || gistA.l.type==L_NONE) return 0;

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
    for (i=0 ; i<nPoints-1 ; i++) Rasterize(points+i, points+i+1);
    if (n==nPoints) break;
    n-= nPoints;
    px+= nPoints;
    py+= nPoints;
  }

  engine->marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrwText(Engine *engine, GpReal x0, GpReal y0, const char *text)
{
  int nlines;
  GpReal width, height, lineHeight;
  GpBox *wind= &engine->transform.window;
  GpReal xmin, xmax, ymin, ymax;
  int alignH= gistA.t.alignH, alignV= gistA.t.alignV;
  int count;
  const char *t;

  charMap= &engine->map;
  ScaleFont(gistA.t.height);

  /* Compute text location in EPS Preview coordinates.  */
  x0= charMap->x.scale*x0 + charMap->x.offset;
  y0= charMap->y.scale*y0 + charMap->y.offset;

  /* handle multi-line strings */
  GtGetAlignment(&gistA.t, &alignH, &alignV);
  nlines= GtTextShape(text, &gistA.t, (WidthFunction)0, &width);
  lineHeight= yChar[6];
  if (lineHeight<0.0) lineHeight= -lineHeight;
  width*= xChar[3];
  height= lineHeight*(GpReal)nlines;

  /* Reject if and only if the specified point is off of the current
     page by more than the approximate size of the text.  Note that
     the entire block of text is either accepted or rejected.  */
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

  /* Adjust y0 (or x0) to represent topmost line */
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

  /* Adjust y0 to represent baseline of topmost line */
  if (gistA.t.orient==TX_RIGHT) {
    if (alignV==TV_TOP) y0-= yChar[6];
    else if (alignV==TV_CAP) y0-= yChar[5];
    else if (alignV==TV_HALF) y0-= yChar[3];
    else if (alignV==TV_BOTTOM) y0-= yChar[0];
  } else if (gistA.t.orient==TX_LEFT) {
    if (alignV==TV_TOP) y0+= yChar[6];
    else if (alignV==TV_CAP) y0+= yChar[5];
    else if (alignV==TV_HALF) y0+= yChar[3];
    else if (alignV==TV_BOTTOM) y0+= yChar[0];
  } else if (gistA.t.orient==TX_UP) {
    if (alignV==TV_TOP) x0+= yChar[6];
    else if (alignV==TV_CAP) x0+= yChar[5];
    else if (alignV==TV_HALF) x0+= yChar[3];
    else if (alignV==TV_BOTTOM) x0+= yChar[0];
  } else {
    if (alignV==TV_TOP) x0-= yChar[6];
    else if (alignV==TV_CAP) x0-= yChar[5];
    else if (alignV==TV_HALF) x0-= yChar[3];
    else if (alignV==TV_BOTTOM) x0-= yChar[0];
  }

  /* Set initial position, set up for horiz centering */
  xCurrent= x0;
  yCurrent= y0;
  width= xChar[3];

  while ((t= GtNextLine(text, &count, gistA.t.orient))) {
    /* Set xCurrent to left edge of current line */
    if (gistA.t.orient==TX_RIGHT) {
      if (alignH==TH_CENTER) xCurrent= x0 - 0.5*width*(GpReal)count;
      else if (alignH==TH_RIGHT) xCurrent= x0 - width*(GpReal)count;
      else xCurrent= x0;
    } else if (gistA.t.orient==TX_UP) {
      if (alignH==TH_CENTER) yCurrent= y0 + 0.5*width*(GpReal)count;
      else if (alignH==TH_RIGHT) yCurrent= y0 + width*(GpReal)count;
      else yCurrent= y0;
    } else if (gistA.t.orient==TX_LEFT) {
      if (alignH==TH_CENTER) xCurrent= x0 + 0.5*width*(GpReal)count;
      else if (alignH==TH_RIGHT) xCurrent= x0 + width*(GpReal)count;
      else xCurrent= x0;
    } else {
      if (alignH==TH_CENTER) yCurrent= y0 - 0.5*width*(GpReal)count;
      else if (alignH==TH_RIGHT) yCurrent= y0 - width*(GpReal)count;
      else yCurrent= y0;
    }

    RenderChars(t, count, gistA.t.orient);
    text= t+count;

    if (gistA.t.orient==TX_RIGHT) yCurrent-= lineHeight;
    else if (gistA.t.orient==TX_UP) xCurrent+= lineHeight;
    else if (gistA.t.orient==TX_LEFT) yCurrent+= lineHeight;
    else xCurrent-= lineHeight;
  }

  engine->marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawFill(Engine *engine, long n, const GpReal *px,
                    const GpReal *py)
{
  return DrawLines(engine, n, px, py, 1, 0);
}

/* ------------------------------------------------------------------------ */

/* ARGSUSED */
static int DrawCells(Engine *engine, GpReal px, GpReal py, GpReal qx,
                     GpReal qy, long width, long height, long nColumns,
                     const GpColor *colors)
{
  GpXYMap *map= &engine->map;
  short ipx, ipy, iqx, iqy, x, y;

  /* Compute image location in EPS Preview coordinates.  */
  ipx= (short)(map->x.scale*px + map->x.offset);
  iqx= (short)(map->x.scale*qx + map->x.offset);
  ipy= (short)(map->y.scale*py + map->y.offset);
  iqy= (short)(map->y.scale*qy + map->y.offset);

  if (ipx>iqx) {short tmp=ipx;  ipx=iqx;  iqx=tmp;}
  if (ipy>iqy) {short tmp=ipy;  ipy=iqy;  iqy=tmp;}

  if (ipx<0) {
    if (iqx<0) return 0;
    ipx= 0;
  }
  if (iqx>255) {
    if (ipx>255) return 0;
    iqx= 255;
  }
  if (ipy<0) {
    if (iqy<0) return 0;
    ipy= 0;
  }
  if (iqy>255) {
    if (ipy>255) return 0;
    iqy= 255;
  }

  /* Represent image as light gray stipple (one cell every 3) */
  ipy= ipy<<5;
  iqy= iqy<<5;
  for (x=ipx ; x<=iqx ; x+=3) {
    for (y=ipy ; y<=iqy ; y+=96)
      epsPreview[y+(x>>3)]|= 0x80 >> (x&0x7);
  }

  engine->marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

static int DrawDisjoint(Engine *engine, long n, const GpReal *px,
                        const GpReal *py, const GpReal *qx, const GpReal *qy)
{
  GpXYMap *map= &engine->map;
  long maxSegs= 2025, nSegs, i;
  GpSegment *segs;
  GpPoint *points;

  if (n<1 || gistA.l.type==L_NONE) return 0;

  while ((nSegs=
          GpIntSegs(map, maxSegs, n, px, py, qx, qy, &segs))) {
    /* This cast from GpSegment to GpPoint is lazy programming,
       but I don't know of any platforms it will not work on... */
    points= (GpPoint *)segs;
    for (i=0 ; i<nSegs ; i++) Rasterize(points+2*i, points+2*i+1);
    if (n==nSegs) break;
    n-= nSegs;
    px+= nSegs;
    py+= nSegs;
    qx+= nSegs;
    qy+= nSegs;
  }

  engine->marked= 1;
  return 0;
}

/* ------------------------------------------------------------------------ */

Engine *EPSPreview(Engine *engine, char *file)
{
  PSEngine *psEngine= (PSEngine *)engine;
  GpTransform toPixels;
  int i;

  /* Grab bounding box from the PS Engine */
  epsLandscape= engine->landscape;
  xll= psEngine->pageBB.xll;
  xur= psEngine->pageBB.xur;
  if (xll < xur) {
    yll= psEngine->pageBB.yll;
    yur= psEngine->pageBB.yur;
  } else {
    xll= yll= 0;
    xur= 612*PS_PER_POINT;
    yur= 792*PS_PER_POINT;
  }
  /* (xll,yll), (xur,yur) are in points, convert to NDC */
  toPixels.viewport.xmin= (GpReal)xll*PS_TO_NDC;
  toPixels.viewport.xmax= (GpReal)xur*PS_TO_NDC;
  toPixels.viewport.ymin= (GpReal)yll*PS_TO_NDC;
  toPixels.viewport.ymax= (GpReal)yur*PS_TO_NDC;
  toPixels.window.xmin= 0.0;
  toPixels.window.xmax= 255.99;
  toPixels.window.ymin= 0.0;
  toPixels.window.ymax= 255.99;

  /* Duplicate code in ps.c:GetEPSFBox */
  xll/= PS_PER_POINT;
  yll/= PS_PER_POINT;
  xur/= PS_PER_POINT;
  yur/= PS_PER_POINT;
  if (epsLandscape) {
    int xl= xll;
    int xu= xur;
    xll= 612-yur;
    xur= 612-yll;
    yll= xl;
    yur= xu;
  }

  /* Compute EPS units per NDC unit for DrwText */
  epsXScale= toPixels.window.xmax /
    (toPixels.viewport.xmax-toPixels.viewport.xmin);
  epsYScale= toPixels.window.ymax /
    (toPixels.viewport.ymax-toPixels.viewport.ymin);

  /* Finish and close the PS file "_tmp.eps" */
  GpKillEngine(engine);

  epsFile= p_fopen(file, "w");
  if (!epsFile) return 0;
  epsPreview= (unsigned char *)p_malloc(sizeof(char)*32*256);
  for (i=0 ; i<32*256 ; i++) epsPreview[i]= 0;
  if (!epsPreview) {
    p_fclose(epsFile);
    return 0;
  }

  /* Create the EPS Engine */
  engine= GpNewEngine(sizeof(Engine), "Gist EPS", epsType, &toPixels, 0,
                      &Kill, &Clear, &Flush, &GpComposeMap,
                      &ChangePalette, &DrawLines, &GpPseudoMark,
                      &DrwText, &DrawFill, &DrawCells,
                      &DrawDisjoint);

  if (!engine) {
    p_free(epsPreview);
    p_fclose(epsFile);
    return 0;
  }

  psFile= p_fopen("_tmp.eps", "r");
  if (!psFile) {
    GpDelEngine(engine);
    p_free(epsPreview);
    p_fclose(epsFile);
    return 0;
  }

  return engine;
}

/* ------------------------------------------------------------------------ */
