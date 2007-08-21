/*
 * errors.c -- $Id$
 * X11 error handling
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"
#include "playu.h"
#include "pstdlib.h"

#include <string.h>

static char x11_errmsg[90];
static int x11_nerrs = 0;

int
x_err_handler(Display *dpy, XErrorEvent *event)
{
  if (!p_signalling) {
    strcpy(x11_errmsg, "Xlib: ");
    XGetErrorText(dpy, event->error_code, x11_errmsg+6, 83);
    x11_errmsg[89] = '\0';
    u_errmsg = x11_errmsg;
    p_signalling = PSIG_SOFT;
    x11_nerrs = 1;
  } else {
    x11_nerrs++;
  }
  /* Xlib apparently ignores return value
   * - must return here, or Xlib internal data structures trashed
   * - therefore actual error processing delayed (u_error must return) */
  return 1;
}

int
x_panic(Display *dpy)
{
  x_display *xdpy = x_dpy(dpy);

  if (xdpy) {
    /* attempt to close display to free all associated Xlib structs
     * -- study of xfree86 source indicates that first call to
     *    XCloseDisplay might well end up calling x_panic
     *    recursively, but a second call might succeed */
    xdpy->panic++;
    while (xdpy->screens) p_disconnect(xdpy->screens);
    if (xdpy->panic<3) XCloseDisplay(dpy);
    xdpy->dpy = 0;
    p_free(xdpy);  /* x_disconnect will not have done this */
  }

  p_abort();
  return 1;  /* Xlib will call exit! */
}
