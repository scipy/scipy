/*
 * XFANCY.H
 *
 * $Id$
 *
 * Declare the fancy X windows engine for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef XFANCY_H
#define XFANCY_H

#include "xbasic.h"

typedef struct FXEngine FXEngine;
struct FXEngine {
  XEngine xe;

  /* --------------- Specific to FXEngine ------------------- */

  Window button;   /* to unlock/cycle coordinate display */
  Window message;  /* place to display coordinates, other messages */
  int baseline;    /* y coordinate for text in button and message windows */
  int heightButton, widthButton;  /* shape of button */
  /* height of both button and message windows is xe.topMargin */
  GC gc;           /* for button and message windows */
  XFontStruct *font;  /* font in gc -- DO NOT FREE */
  Cursor cursor;   /* cursor for xe.graphics window */

  int buttonState;  /* 0 -  button inactive
		       1 -  pointer in button, button ready
		       2 -  button pressed, button armed  */
  int iSystem;      /* <0 for unlocked, else locked system number */
  char msgText[64]; /* current text displayed in message window */
  int msgWidth;     /* width in pixels of msgText when last displayed */
  int zoomState;    /* button number if zoom or pan in progress, else 0 */
  int zoomSystem;   /* system number in xe.drawing */
  int zoomAxis;     /* 1 for x-axis, 2 for y-axis, 3 for both */
  GpReal zoomX, zoomY; /* initial coordinates for zoom/pan operation */
};

/* zoom factor for point-and-click zooming */
extern GpReal gxZoomFactor;

/* The GxPointClick function initiates an interactive point-and-click
   session with the window -- it will not return until a button has
   been pressed, then released.  It returns non-zero if the operation
   was aborted by pressing a second button before releasing the first.
     engine --   an X engine whose display is to be used
     style --    1 to draw a rubber box, 2 to draw a rubber line,
                 otherwise, no visible indication of operation
     system --   system number to which the world coordinates should
                 be transformed, or -1 to use the system under the
		 pointer -- the release coordinates are always in the
		 same system as the press coordinates
     CallBack -- function to be called twice, first when the button is
                 pressed, next when it is released -- operation will
		 be aborted if CallBack returns non-zero
		 Arguments passed to CallBack:
		   engine  -- in which press/release occurred
		   system  -- system under pointer, if above system -1
		   release -- 0 on press, 1 on release
		   x, y    -- coordinates of pointer relative to system
		   butmod  -- 1 - 5 on press to tell which button
		              mask to tell which modifiers on release:
			      1 shift, 2 lock, 4 control, 8 - 128 mod1-5
		   xn, yn  -- NDC coordinates of pointer
 */
extern int GxPointClick(Engine *engine, int style, int system,
			int (*CallBack)(Engine *engine, int system,
					int release, GpReal x, GpReal y,
					int butmod, GpReal xn, GpReal yn));

#endif
