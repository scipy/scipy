/*
 * GTEXT.C
 *
 * $Id$
 *
 * Define GIST text utilities
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "gtext.h"

extern long strcspn(const char *, const char *);

int gtDoEscapes= 1;

/* Return t->alignH, t->alignV, guaranteed not TH_NORMAL or TV_NORMAL */
void GtGetAlignment(const GpTextAttribs *t, int *alignH, int *alignV)
{
  *alignH= t->alignH;
  *alignV= t->alignV;
  if (*alignH==TH_NORMAL) *alignH= TH_LEFT;
  if (*alignV==TV_NORMAL) *alignV= TV_BASE;
}

/* Get shape of text input to GdText, given a function Width which can
   compute the width of a simple text string (no imbedded \n).  Returns
   largest value of Width for any line, and a line count.  */
int GtTextShape(const char *text, const GpTextAttribs *t,
                WidthFunction Width, GpReal *widest)
{
  int path= t->orient;
  GpReal wdest, thisWid;
  int nLines, nChars;

  /* Count lines in this text, find widest line */
  nLines= 0;
  wdest= 0.0;
  while ((text= GtNextLine(text, &nChars, path))) {
    nLines++;
    if (Width) thisWid= Width(text, nChars, t);
    else thisWid= (GpReal)nChars;
    if (thisWid>wdest) wdest= thisWid;
    text+= nChars;
  }

  *widest= wdest;
  return nLines;
}

/* Return the next line of text--
   returns text and a count of the characters in the
   line, nChars.  If text is '\0', or '\n', returns text+1 and
   a count of the number of characters to the next '\n' or '\0'.  */
const char *GtNextLine(const char *text, int *nChars, int path)
{
  char first= text[0];

  if (!first) {
    *nChars= 0;
    return 0;
  }

  if (first=='\n') text+= 1;
  *nChars= strcspn(text, "\n");

  return text;
}
