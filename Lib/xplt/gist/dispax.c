/*
 * DISPAX.C
 *
 * $Id$
 *
 * Implement dispatcher routines for X window servers
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "dispax.h"

#include <X11/Xutil.h>

/* XUniqueContext is front for XrmUniqueQuark */
#include <X11/Xresource.h>

static int xContextSet= 0;
static XContext xDispatchContext;

typedef int (*XCallback)(XEvent *);

static XCallback GetXDispatcher(Display *dpy, Window win, Window orig);
static int DispatchX(void *context);
static void FlushX(void *context);
static int PendingX(void *context);

static XCallback GetXDispatcher(Display *dpy, Window win, Window orig)
{
  XCallback callback;
  if (XFindContext(dpy, win, xDispatchContext, (caddr_t *)&callback)) {
    /* No xDispatchContext associated with this window, recurse to parent */
    Window root, parent, *children;
    unsigned int nChildren;
    if (!XQueryTree(dpy, win, &root, &parent, &children, &nChildren))
      return 0;
    XFree((caddr_t)children);
    if (root==parent || root==win) return 0;
    callback= GetXDispatcher(dpy, parent, orig);
    if (callback) XSaveContext(dpy, orig, xDispatchContext, (void *)callback);
  }
  return callback;
}

static int DispatchX(void *context)
{
  Display *display= context;
  XEvent event;
  XCallback Dispatch;

  do {
    XNextEvent(display, &event);

    /* Note-- Although SelectionNotify and SelectionRequest events do
       not have a member named "window", their requestor and owner
       members serve the same purpose in this context.
       A more serious problem is presented by the GraphicsExpose
       and NoExpose events, whose drawable member might be a Pixmap.
       If it is, the XQueryTree call in GetXDispatcher will produce
       a BadWindow error.  This can be avoided by treating Pixmaps as
       top level windows, or by turning off graphics exposure events
       in the GC, or, most simply, by not using Pixmaps as the destination
       of XCopyArea or XCopyPlanes calls.  */

    Dispatch= GetXDispatcher(display, event.xany.window, event.xany.window);
    if (Dispatch && Dispatch(&event)) return 1;
  } while (QLength(display));
  return 0;
}

static void FlushX(void *context)
{
  Display *display= context;
  XFlush(display);
}

static int PendingX(void *context)
{
  Display *display= context;
  if (QLength(display)) return DispatchX(context);
  return 0;
}

int AddXDispatcher(Display *dpy, Window win, int (*Dispatch)(XEvent *event))
{
  int fd= ConnectionNumber(dpy);
  void *callback= (void *)Dispatch;

  /* First, create a Dispatcher for this display, if there isn't one */
  if (!HasDispatcher(fd)) {
    int val= AddDispatcher(fd, dpy, &PendingX, &FlushX, &DispatchX);
    if (val) return val;
  }

  /* Next, make sure that xDispatchContext is set */
  if (!xContextSet) {
    xDispatchContext= XUniqueContext();
    xContextSet= 1;
  }

  /* Finally, associate the given Dispatch method with the window */
  XSaveContext(dpy, win, xDispatchContext, callback);
  return 0;
}

void CopyXDispatcher(Display *dpy, Window src, Window dst)
{
  XCallback callback= dpy? GetXDispatcher(dpy, src, src) : 0;
  if (callback) XSaveContext(dpy, dst, xDispatchContext, (void *)callback);
}

void RemoveXDispatcher(Display *dpy)
{
  if (dpy) RemoveDispatcher(ConnectionNumber(dpy));
}
