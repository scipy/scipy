/*
 * GTEXT.H
 *
 * $Id$
 *
 * Declare GIST text utilities
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef GTEXT_H
#define GTEXT_H

#include "gist.h"

/* Return t->alignH, t->alignV, guaranteed not TH_NORMAL or TV_NORMAL */
extern void GtGetAlignment(const GpTextAttribs *t,
			   int *alignH, int *alignV);

/* Get shape of text input to GdText, given a function Width which can
   compute the width of a simple text string (no imbedded \n).  Returns
   largest value of Width for any line, and a line count.
   If gtDoEscapes==1 (default), and Width!=0, the Width function must
   handle the !, ^, and _ escape sequences.
   If Width==0, the default Width function is the number of
   non-escaped characters.  */
typedef GpReal (*WidthFunction)(const char *text, int nChars,
				const GpTextAttribs *t);
extern int GtTextShape(const char *text, const GpTextAttribs *t,
		       WidthFunction Width, GpReal *widest);

/* Return the next line of text-- if text[0] is not '\n' and path not
   T_UP or T_DOWN, returns text and a count of the characters in the
   line, nChars (always 1 if T_UP or T_DOWN).  If text is '\0', or '\n'
   with path T_UP or T_DOWN, returns 0.  Otherwise, returns text+1 and
   a count of the number of characters to the next '\n' or '\0'.  */
extern const char *GtNextLine(const char *text, int *nChars, int path);

/* Use ^ (superscript), _ (subscript), and ! (symbol) escape sequences
   in GpText.  These are on (gtDoEscapes==1) by default.  */
extern int gtDoEscapes;

#endif
