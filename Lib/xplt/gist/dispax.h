/*
 * DISPAX.H
 *
 * $Id$
 *
 * Declare dispatcher routines for X window servers
 *
 * These routines allow a program to use more than one X toolkit
 * simultaneously.  This is tricky, because input events must be
 * directed to the event handler for the toolkit which owns the
 * window.
 *
 * Notes--
 *
 * 1) You need to write a special event dispatching routine for each
 *    toolkit you want to use, then call AddXDispatcher for each top
 *    level window.  Call DispatchEvents to run your main event loop.
 *    Note that the event has already been removed from the event queue
 *    when your dispatch method is called; you may need to put it back
 *    using XPutBackEvent (e.g.- for InterViews).
 *    (Use AddFDispatcher to add stream input instead of a toolkit
 *     specific routine such as X toolkit intrinsic's XtAppAddInput.)
 *
 * 2) WARNING--
 *    If you use XCopyArea or XCopyPlanes with a Pixmap as a destination
 *    (hopefully this is rare to begin with), you must treat the Pixmap
 *    as a top level window, calling AddXDispatcher with it.  If you
 *    don't do this, at least turn off the graphics exposure bit in
 *    the graphics context for the operation, since otherwise a
 *    BadWindow error will be generated when this package tries to
 *    use XQueryTree to find the parent of the Pixmap.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef DISPAX_H
#define DISPAX_H

#include <X11/Xlib.h>

/* When an event arrives for window on display, pass it to the
   given Dispatch method.  If no Dispatch method has been registered
   for this window, its ancestors are searched until a Dispatch method
   is found; if found, the window and its unregistered ancestors are
   registered with the Dispatch routine that was found.  Hence, you
   need only register the highest level window for which Dispatch
   is to be called.  */
extern int AddXDispatcher(Display *dpy, Window win,
			  int (*Dispatch)(XEvent *event));

extern void CopyXDispatcher(Display *dpy, Window src, Window dst);

/* You should call RemoveXDispatcher only if you are about to close
   the display.  Individual windows used in AddXDispatcher may be
   safely destroyed without informing the dispatcher.  */
extern void RemoveXDispatcher(Display *dpy);

/* need dispat.h for DispatchEvents declaration */
#include "dispat.h"

#endif
