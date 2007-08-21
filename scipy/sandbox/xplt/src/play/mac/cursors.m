/*
 * cursors.m
 * p_cursor for Mac OS X.
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playm.h"

static NSCursor* m_cursor(int cursor);

void
p_cursor(p_win *w, int cursor)
{
  p_scr* s = w->s;
  if (cursor==P_NONE) {
    [NSCursor hide];
    return;
  }
  if (cursor<0 || cursor>P_NONE) cursor = P_SELECT;
  if (!s->cursors[cursor]) s->cursors[cursor] = m_cursor(cursor);
  [s->cursors[cursor] set];
}

static NSCursor* m_cursor(int cursor)
{
  NSCursor* c = 0;
  SEL selector;
  switch (cursor)
  { case P_SELECT: selector = @selector(arrowCursor); break;
    case P_CROSSHAIR: selector = @selector(crosshairCursor); break;
    case P_TEXT: selector = @selector(IBeamCursor); break;
    case P_N: selector = @selector(resizeUpCursor); break;
    case P_S: selector = @selector(resizeDownCursor); break;
    case P_E: selector = @selector(resizeRightCursor); break;
    case P_W: selector = @selector(resizeLeftCursor); break;
    case P_NS: selector = @selector(resizeUpDownCursor); break;
    case P_EW: selector = @selector(resizeLeftRightCursor); break;
    case P_NSEW: selector = @selector(closedHandCursor); break;
    case P_ROTATE: selector = @selector(arrowCursor); break;
    case P_DEATH: selector = @selector(arrowCursor); break;
    default: selector = @selector(arrowCursor);
  }
  /* For P_NSEW, closedHandCursor seems suitable for usage in Pygist
   * (zooming operations). P_ROTATE and P_DEATH have arrowCursor as
   * place holders. */
  if ([NSCursor respondsToSelector:selector])
    c = [NSCursor performSelector: selector];
  else
    c = [NSCursor arrowCursor];
  /* Mac OS X older than 10.3 don't have all required cursors;
   * the arrowCursor is always present. */
  [c retain];
  return c;
}
