/*
 * getdc.c -- $Id$
 * get a DC to draw into a window
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

/* WARNING
 * implementations of p_txwidth and p_txheight assume
 * that p_font will be called afterwards
 */
static HDC wp_font(p_win *w, int font, int pixsize, int orient);
static void w_font(p_scr *s, HDC dc, int font, int pixsize, int orient);

void
p_font(p_win *w, int font, int pixsize, int orient)
{
  wp_font(w, font, pixsize, orient);
}

static HDC
wp_font(p_win *w, int font, int pixsize, int orient)
{
  HDC dc = w_getdc(w, 0);
  if (dc) {
    if (!w->s->does_rotext) {
      w->orient = orient;
      orient = 0;
    }
    if (w->font!=font || w->pixsize!=pixsize) {
      w_font(w->s, dc, font, pixsize, orient);
      w->font = font;
      w->pixsize = pixsize;
      w->s->font_win = w;
    }
  }
  return dc;
}

static void
w_font(p_scr *s, HDC dc, int font, int pixsize, int orient)
{
  int face = ((unsigned int)(font&0x1c))>>2;

  if (face >= 5) {
    SelectObject(dc, s->gui_font);

  } else {
    int i, j;

    for (i=j=0 ; i<W_FONTS_CACHED && s->font_cache[i].hfont ; i++) {
      j = s->font_order[i];
      if (s->font_cache[j].font==font &&
          s->font_cache[j].pixsize==pixsize &&
          s->font_cache[j].orient==orient) {
        /* install cached font and make it most recently used */
        for (; i>0 ; i--) s->font_order[i] = s->font_order[i-1];
        s->font_order[0] = j;
        SelectObject(dc, s->font_cache[j].hfont);
        return;
      }
    }
    if (i<W_FONTS_CACHED) j = i++;

    {
      static DWORD families[5] = {
        FF_MODERN|FIXED_PITCH, FF_ROMAN|VARIABLE_PITCH,
        FF_SWISS|VARIABLE_PITCH, VARIABLE_PITCH, FF_ROMAN|VARIABLE_PITCH };
      /* note: my MS Windows box does not have Helvetica -- is Arial
       *       equivalent? -- hopefully substitution will be reasonable */
      static char *names[5] = { "Courier", "Times New Roman", "Helvetica",
                               "Symbol", "Century Schoolbook" };
      static int angles[4] = { 0, 900, 1800, 2700 };
      int ang = angles[orient];
      int weight = (font&P_BOLD)? FW_BOLD : FW_NORMAL;
      DWORD charset = ((font&P_SYMBOL)!=P_SYMBOL)? ANSI_CHARSET:SYMBOL_CHARSET;
      HFONT hfont = CreateFont(pixsize, 0, ang, ang, weight,
                              (font&P_ITALIC)!=0, 0, 0, charset,
                              OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                              PROOF_QUALITY, families[face], names[face]);
      HFONT fold = s->font_cache[j].hfont;
      for (i-- ; i>0 ; i--)
        s->font_order[i] = s->font_order[i-1];
      s->font_order[0] = j;
      s->font_cache[j].hfont = hfont;
      s->font_cache[j].font = font;
      s->font_cache[j].pixsize = pixsize;
      s->font_cache[j].orient = orient;
      if (fold) DeleteObject(fold);
      SelectObject(dc, hfont);
    }
  }
}

int
p_txheight(p_scr *s, int font, int pixsize, int *baseline)
{
  int height = pixsize;
  int base = pixsize;
  if (!p_signalling) {
    HFONT orig = 0;
    HDC dc;
    TEXTMETRIC tm;
    if (s->font_win) {
      dc = wp_font(s->font_win, font, pixsize, s->font_win->orient);
    } else {
      dc = GetDC(0);
      orig = dc? SelectObject(dc, GetStockObject(ANSI_FIXED_FONT)) : 0;
      w_font(s, dc, font, pixsize, 0);
    }
    if (dc && GetTextMetrics(dc, &tm)) {
      base = tm.tmAscent;
      height = tm.tmHeight + 1;
    }
    if (orig) SelectObject(dc, orig);
  } else {
    p_abort();
  }
  if (baseline) *baseline = base;
  return height;
}

int
p_txwidth(p_scr *s, const char *text, int n, int font, int pixsize)
{
  int width = pixsize;
  if (!p_signalling) {
    HFONT orig = 0;
    HDC dc;
    SIZE sz;
    if (s->font_win) {
      dc = wp_font(s->font_win, font, pixsize, s->font_win->orient);
    } else {
      dc = GetDC(0);
      orig = dc? SelectObject(dc, GetStockObject(ANSI_FIXED_FONT)) : 0;
      w_font(s, dc, font, pixsize, 0);
    }
    if (dc && GetTextExtentPoint32(dc, text, n, &sz))
      width = sz.cx;
    else
      width = n*(pixsize*3/5);
    if (orig) SelectObject(dc, orig);
  } else {
    p_abort();
  }
  return width;
}

void
p_pen(p_win *w, int width, int type)
{
  if (p_signalling) {
    p_abort();
  } else {
    if (width<1) width = 1;
    else if (width>100) width = 100;
    if (w->pen_width!=width || w->pen_type!=type) {
      w->pen_width = width;
      w->pen_type = type;
      w->pen = 0;
    }
  }
}

void
p_color(p_win *w, unsigned long color)
{
  if (p_signalling) p_abort(), color = w->color;
  if (w->color != color) {
    int de_xor = (w->color == P_XOR);
    w->color = color;
    if (de_xor || color==P_XOR) {
      if (w->dc)
        SetROP2(w->dc, de_xor? R2_COPYPEN : R2_NOT);  /* or R2_XORPEN */
    }
  }
}

COLORREF
w_color(p_win *w, unsigned long color)
{
  if (w->parent) w = w->parent;
  if (!P_IS_RGB(color)) {
    return w->pixels[color];
  } else {
    unsigned int r = P_R(color);
    unsigned int g = P_G(color);
    unsigned int b = P_B(color);
    if (w->w && !w->menu && w->s->sys_pal) {
      if (!w->rgb_mode) {
        p_palette(w, p_595, 225);
        w->rgb_mode = 1;
      }
      /* must be consistent with p_palette in pals.c */
      r = (r+32)>>6;
      g = (g+16)>>5;
      b = (b+32)>>6;
      g += b+(b<<3);
      return PALETTEINDEX(w->s->sys_offset+r+g+(g<<2));  /* r+5*g+45*b */
    } else {
      return RGB(r, g, b);
    }
  }
}

HDC
w_getdc(p_win *w, int flags)
{
  HDC dc = w->dc;
  if (p_signalling) p_abort(), dc = 0;
  if (dc) {
    if (flags) {
      COLORREF color = (flags&7)? w_color(w, w->color) : 0;

      if ((flags & 1) && (w->font_color!=w->color)) {
        /* text */
        w->font_color = w->color;
        SetTextColor(dc, color);
      }

      if (flags & 2) {
        /* pen */
        HPEN pen = 0;
        if (!w->pen || (w->pen_color != w->color)) {
          int width = w->pen_width;
#ifdef USE_GEOMETRIC_PEN
          /* NOTE: geometric pen is way to slow for practical use */
          int type = w->pen_type;
          static DWORD styles[8] = {
            PS_SOLID, PS_DASH, PS_DOT, PS_DASHDOT, PS_DASHDOTDOT, 0, 0, 0 };
          DWORD ltype = w->s->does_linetypes? styles[type&7] : PS_SOLID;
          DWORD cap = (type&P_SQUARE)? PS_ENDCAP_SQUARE : PS_ENDCAP_ROUND;
          DWORD join = (type&P_SQUARE)? PS_JOIN_MITER : PS_JOIN_ROUND;
          LOGBRUSH brush;
          brush.lbStyle = BS_SOLID;
          brush.lbColor = color;
          brush.lbHatch = 0;
          pen = ExtCreatePen(PS_GEOMETRIC | ltype | cap | join,
                             width, &brush, 0, 0);
#else
          /* downside of cosmetic pen is always round endcaps and joins,
           * so P_SQUARE type modifier is ignored
           * -- go for high performance instead */
          pen = CreatePen(PS_SOLID, width, color);
#endif
          if (!pen) return 0;
          w->pen = pen;
          w->pen_color = w->color;
        } else if (w->pbflag&1) {
          /* null pen has been installed */
          pen = w->pen;
        }
        if (pen) {
          pen = SelectObject(dc, pen);
          if (pen) DeleteObject(pen);
          w->pbflag &= ~1;
        }
      }

      if (flags & 4) {
        /* brush */
        HBRUSH b = 0;
        if (!w->brush || (w->brush_color != w->color)) {
          b = CreateSolidBrush(color);
          if (!b) return 0;
          w->brush = b;
          w->brush_color = w->color;
        } else if (w->pbflag&2) {
          b = w->brush;
        }
        if (b) {
          b = SelectObject(dc, b);
          if (b) DeleteObject(b);
          w->pbflag &= ~2;
        }
      }

      if (flags & 8) {
        /* install null pen */
        HPEN pen = CreatePen(PS_NULL, 0, 0);
        if (pen) {
          pen = SelectObject(dc, pen);
          if (pen && pen==w->pen) w->pbflag |= 1;
        }
      }

      if (flags & 16) {
        /* install null brush */
        HBRUSH b = SelectObject(dc, GetStockObject(NULL_BRUSH));
        if (b && b==w->brush) w->pbflag |= 2;
      }
    }
  }
  return dc;
}
