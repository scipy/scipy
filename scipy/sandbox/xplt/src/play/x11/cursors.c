/*
 * cursors.c -- $Id$
 * p_cursor for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"
#include <X11/cursorfont.h>

static unsigned int cursor_shape[12]= {
  XC_left_ptr, XC_crosshair, XC_xterm,
  XC_sb_up_arrow, XC_sb_down_arrow, XC_sb_right_arrow, XC_sb_left_arrow,
  XC_sb_v_double_arrow, XC_sb_h_double_arrow, XC_fleur,
  XC_exchange, XC_pirate };

#define hand_width 16
#define hand_height 16
#define hand_x_hot 3
#define hand_y_hot 2
static unsigned char hand_bits[] = {
   0x80, 0x01, 0x58, 0x0e, 0x64, 0x12, 0x64, 0x52, 0x48, 0xb2, 0x48, 0x92,
   0x16, 0x90, 0x19, 0x80, 0x11, 0x40, 0x02, 0x40, 0x04, 0x40, 0x04, 0x20,
   0x08, 0x20, 0x10, 0x10, 0x20, 0x10, 0x20, 0x10};
#define HAND_BITS (char *)hand_bits
#define hand_mask_width 16
#define hand_mask_height 16
static unsigned char hand_mask_bits[] = {
   0x80, 0x01, 0xd8, 0x0f, 0xfc, 0x1f, 0xfc, 0x5f, 0xf8, 0xff, 0xf8, 0xff,
   0xfe, 0xff, 0xff, 0xff, 0xff, 0x7f, 0xfe, 0x7f, 0xfc, 0x7f, 0xfc, 0x3f,
   0xf8, 0x3f, 0xf0, 0x1f, 0xe0, 0x1f, 0xe0, 0x1f};
#define HAND_MASK_BITS (char *)hand_mask_bits

void
p_cursor(p_win *w, int cursor)
{
  p_scr *s = w->s;
  Display *dpy = s->xdpy->dpy;
  if (dpy) XDefineCursor(dpy, w->d, x_cursor(s, cursor));
  if (p_signalling) p_abort();
}

Cursor
x_cursor(p_scr *s, int cursor)
{
  x_display *xdpy = s->xdpy;
  Display *dpy = xdpy->dpy;
  if (dpy && cursor>=0 && cursor<=13) {
    if (xdpy->cursors[cursor]==None) {
      if (cursor==P_HAND || cursor==P_NONE) {
        /* it is disturbing that you can call XCreateFontCursor without
         * knowing what screen the cursor will be on, but not
         * XCreatePixmapCursor -- I assume the cursor that is finally
         * created can be used on any screen of the display... */
        Window root = RootWindow(dpy,s->scr_num);
        Pixmap hand, mask;
        XColor *fg = &s->colors[1];
        XColor *bg = &s->colors[0];
        char hbits[32], mbits[32], *hptr, *mptr;
        int i;
        for (i=0 ; i<32 ; i++) hbits[i]= mbits[i]= '\0';
        if (cursor==P_HAND) {
          hptr = HAND_BITS;
          mptr = HAND_MASK_BITS;
        } else {
          hptr = hbits;
          mptr = mbits;
        }
        hand = XCreatePixmapFromBitmapData(dpy, root, hptr,
                       hand_width, hand_height, 1, 0, 1);
        mask = XCreatePixmapFromBitmapData(dpy, root, mptr,
                       hand_mask_width, hand_mask_height, 1, 0, 1);
        xdpy->cursors[cursor]= XCreatePixmapCursor(dpy, hand, mask, fg, bg,
                                                   hand_x_hot, hand_y_hot);
        XFreePixmap(dpy, hand);
        XFreePixmap(dpy, mask);
      } else {
        xdpy->cursors[cursor]= XCreateFontCursor(dpy, cursor_shape[cursor]);
      }
      if (p_signalling) p_abort();
    }
    return xdpy->cursors[cursor];
  } else {
    return None;
  }
}
