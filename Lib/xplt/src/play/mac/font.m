/*
 * font.m
 * Font handling routines for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

static char* families[5][4] =
  {{"CourierNewPSMT",
    "CourierNewPS-BoldMT",
    "CourierNewPS-ItalicMT",
    "CourierNewPS-BoldItalicMT"},
   {"Times-Roman",
    "Times-Bold",
    "Times-Italic",
    "Times-BoldItalic"},
   {"ArialMS",
    "Arial-BoldMS",
    "Arial-ItalicMS",
    "Arial-BoldItalicMS"},
   {"Symbol",
    "Symbol",
    "Symbol",
    "Symbol"},
   {"CenturySchoolbook",
    "CenturySchoolbook-Bold",
    "CenturySchoolbook-Italic",
    "CenturySchoolbook-BoldItalic"}};

static void
m_font(p_scr* s, int fontid, CGFontRef* fontref, ATSFontMetrics** metric);
    
void
p_font(p_win *w, int fontid, int pixsize, int orient)
{ p_scr* s = w->s;
  CGContextRef cr = w->cr;
  if (p_signalling) {
    p_abort();
    return;
  }
  if (fontid != w->fontid)
  { m_font(s, fontid, &(w->fontref), NULL);
    w->fontid = fontid;
  }
  w->orient = orient;
  CGContextSetFont(cr, w->fontref);
  CGContextSetFontSize(cr, (float) pixsize);
}

int
p_txheight(p_scr *s, int fontid, int pixsize, int *baseline)
{ ATSFontMetrics* metric;
  int height = pixsize;
  int base = pixsize;
  if (!p_signalling) {
    float ascent, descent;
    m_font(s, fontid, NULL, &metric);
    ascent  = pixsize * (metric->ascent);
    descent = pixsize * (metric->descent);
    base = (int) ascent;
    height = (int) (ascent - descent);
  } else {
    p_abort();
  }
  if (baseline) *baseline = base;
  return height;
}

int
p_txwidth(p_scr *s, const char *text, int n, int fontid, int pixsize)
{ int width = pixsize;
  if (!p_signalling) {
    int i;
    CGPoint point;
    CGContextRef cr = s->scratch;
    CGFontRef font;
    m_font(s, fontid, &font, NULL);
    CGContextSetFont(cr, font);
    CGContextSetFontSize(cr, (float) pixsize);
    CGContextSetTextPosition(cr, 0.0, 0.0);
    for (i = 0; i < n; i++)
    { CGGlyph c = (unsigned short)(text[i]) - 29;
      CGContextShowGlyphs(cr, &c, 1);
    }
    point = CGContextGetTextPosition(cr);
    width = (int) (point.x);
  } else {
    p_abort();
  }
  return width;
}

static void 
m_font(p_scr* s, int fontid, CGFontRef* fontref, ATSFontMetrics** metric)
{ int index;
  int face = ((unsigned int)(fontid&0x1c))>>2;
  if (face >= 5) {
    if(fontref) *fontref = s->font;
    if(metric) *metric = &(s->metric);
    return;
  }
  else index = (fontid&P_ITALIC) + (fontid&P_BOLD);
  if (!(s->fonts[face][index]))
  { CFStringRef family = CFStringCreateWithCString(kCFAllocatorDefault,
                                                   families[face][index],
                                                   kCFStringEncodingMacRoman);
    ATSFontRef atsfont = ATSFontFindFromPostScriptName(family,
                                                       kATSOptionFlagsDefault);
    /* Some fonts are not found, although the family name is correct and
     * the font can be found with [NSFont fontWithName: size:] */
    ATSFontMetrics* newmetric = &(s->metrics[face][index]);
    ATSFontGetHorizontalMetrics(atsfont, kATSOptionFlagsDefault, newmetric);
    s->fonts[face][index] = CGFontCreateWithPlatformFont((void*)&atsfont);
    CFRelease(family);
    /* atsfont is an unsigned integer (UInt32); I assume it should not be
     * released. */
  }
  if(fontref) *fontref = s->fonts[face][index];
  if(metric) *metric = &(s->metrics[face][index]);
}
