/*
 * pscr.m
 * routines to initialize graphics for Mac OS X.
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playm.h"
#include "pstdlib.h"

static NSAutoreleasePool *pool = 0;
static p_scr m_screen;

static int m_checksig(void);
static void (*m_exception)(int signal, char *errmsg) = 0;
static void (*m_abort_hook)(void) = 0;

void (*p_on_connect)(int dis, int fd) = 0;

p_scr *
p_connect(char *server_name)
{
  int i, j;
 
  if (server_name) return 0;

  NSApp = [NSApplication sharedApplication];
  pool = [[NSAutoreleasePool alloc] init];

  NSScreen* screen = [NSScreen mainScreen];
  NSRect r = [screen visibleFrame];

  m_screen.x0 = r.origin.x;
  m_screen.y0 = r.origin.y;
  m_screen.width = r.size.width;
  m_screen.height = r.size.height;
  m_screen.depth = NSBitsPerPixelFromDepth([screen depth]);

  m_screen.colorspace = CGColorSpaceCreateDeviceRGB();
  m_screen.sys_colors[255-P_BLACK][0] = 0.0;
  m_screen.sys_colors[255-P_BLACK][1] = 0.0;
  m_screen.sys_colors[255-P_BLACK][2] = 0.0;
  m_screen.sys_colors[255-P_BLACK][3] = 1.0;
  m_screen.sys_colors[255-P_WHITE][0] = 1.0;
  m_screen.sys_colors[255-P_WHITE][1] = 1.0;
  m_screen.sys_colors[255-P_WHITE][2] = 1.0;
  m_screen.sys_colors[255-P_WHITE][3] = 1.0;
  m_screen.sys_colors[255-P_RED][0]   = 1.0;
  m_screen.sys_colors[255-P_RED][1]   = 0.0;
  m_screen.sys_colors[255-P_RED][2]   = 0.0;
  m_screen.sys_colors[255-P_RED][3]   = 1.0;
  m_screen.sys_colors[255-P_GREEN][0] = 0.0;
  m_screen.sys_colors[255-P_GREEN][1] = 1.0;
  m_screen.sys_colors[255-P_GREEN][2] = 0.0;
  m_screen.sys_colors[255-P_GREEN][3] = 1.0;
  m_screen.sys_colors[255-P_BLUE][0]  = 0.0;
  m_screen.sys_colors[255-P_BLUE][1]  = 0.0;
  m_screen.sys_colors[255-P_BLUE][2]  = 1.0;
  m_screen.sys_colors[255-P_BLUE][3]  = 1.0;
  m_screen.sys_colors[255-P_CYAN][0]  = 0.0;
  m_screen.sys_colors[255-P_CYAN][1]  = 1.0;
  m_screen.sys_colors[255-P_CYAN][2]  = 1.0;
  m_screen.sys_colors[255-P_CYAN][3]  = 1.0;
  m_screen.sys_colors[255-P_MAGENTA][0] = 1.0;
  m_screen.sys_colors[255-P_MAGENTA][1] = 0.0;
  m_screen.sys_colors[255-P_MAGENTA][2] = 1.0;
  m_screen.sys_colors[255-P_MAGENTA][3] = 1.0;
  m_screen.sys_colors[255-P_YELLOW][0] = 1.0;
  m_screen.sys_colors[255-P_YELLOW][1] = 1.0;
  m_screen.sys_colors[255-P_YELLOW][2] = 0.0;
  m_screen.sys_colors[255-P_YELLOW][3] = 1.0;
  m_screen.sys_colors[255-P_BG][0] = 1.0;
  m_screen.sys_colors[255-P_BG][1] = 1.0;
  m_screen.sys_colors[255-P_BG][2] = 1.0;
  m_screen.sys_colors[255-P_BG][3] = 1.0;
  m_screen.sys_colors[255-P_FG][0] = 0.0;
  m_screen.sys_colors[255-P_FG][1] = 0.0;
  m_screen.sys_colors[255-P_FG][2] = 0.0;
  m_screen.sys_colors[255-P_FG][3] = 1.0;
  m_screen.sys_colors[255-P_GRAYD][0] = 100.0/255.0;
  m_screen.sys_colors[255-P_GRAYD][1] = 100.0/255.0;
  m_screen.sys_colors[255-P_GRAYD][2] = 100.0/255.0;
  m_screen.sys_colors[255-P_GRAYD][3] = 1.0;
  m_screen.sys_colors[255-P_GRAYC][0] = 150.0/255.0;
  m_screen.sys_colors[255-P_GRAYC][1] = 150.0/255.0;
  m_screen.sys_colors[255-P_GRAYC][2] = 150.0/255.0;
  m_screen.sys_colors[255-P_GRAYC][3] = 1.0;
  m_screen.sys_colors[255-P_GRAYB][0] = 190.0/255.0;
  m_screen.sys_colors[255-P_GRAYB][1] = 190.0/255.0;
  m_screen.sys_colors[255-P_GRAYB][2] = 190.0/255.0;
  m_screen.sys_colors[255-P_GRAYB][3] = 1.0;
  m_screen.sys_colors[255-P_GRAYA][0] = 214.0/255.0;
  m_screen.sys_colors[255-P_GRAYA][1] = 214.0/255.0;
  m_screen.sys_colors[255-P_GRAYA][2] = 214.0/255.0;
  m_screen.sys_colors[255-P_GRAYA][3] = 1.0;

  /* It seems that Cocoa doesn't support XOR drawing at this point.
  m_screen.sys_colors[255-P_XOR] = 0xffffff;
  */

  /* Initially no fonts are installed, except for the standard GUI font.
   * Fonts are installed as needed during the run. */
  for (i = 0; i < 5; i++)
    for (j = 0; j < 4; j++)
      m_screen.fonts[i][j] = NULL;
  /* The scratch space is used for drawing strings to determine their width. */
  CFURLRef url = CFURLCreateWithFileSystemPath(kCFAllocatorDefault, 
                                        CFSTR("/dev/null"),
                                        kCFURLPOSIXPathStyle,
                                        false);
  CGDataConsumerRef trashcan = CGDataConsumerCreateWithURL(url);
  m_screen.scratch = CGPDFContextCreate(trashcan, NULL, NULL);
  CGDataConsumerRelease(trashcan);
  CFRelease(url);
  const char* fontname = [[[NSFont userFixedPitchFontOfSize: 0.0] fontName] cString];
  CFStringRef family = CFStringCreateWithCString(kCFAllocatorDefault,
                                                 fontname,
                                                 kCFStringEncodingMacRoman);
  ATSFontRef atsfont = ATSFontFindFromPostScriptName(family,
                                                     kATSOptionFlagsDefault);
  ATSFontGetHorizontalMetrics(atsfont, kATSOptionFlagsDefault, &(m_screen.metric));
  m_screen.font = CGFontCreateWithPlatformFont((void*)&atsfont);
  CFRelease(family);

  if (p_on_connect) p_on_connect(0, -1);

  return &m_screen;
}

void
p_pending_events(void)
{
  NSEvent* event;
  View* view = m_screen.lockedView;
  m_checksig();
  if (view) 
  { [view unlockFocus];
    m_screen.lockedView = NULL;
    [[view window] flushWindow];
  }
  p_on_idle(0);
  while ((event = [NSApp nextEventMatchingMask: NSAnyEventMask
                                     untilDate: [NSDate distantPast]
                                        inMode: NSDefaultRunLoopMode
                                       dequeue: YES]))
  { [NSApp sendEvent: event];
    m_checksig();
  }
}

int
p_sshape(p_scr *s, int *width, int *height)
{
  *width = s->width;
  *height = s->height;
  return s->depth;
}

p_scr *
p_multihead(p_scr *other, int number)
{
  return 0;
}

void
p_disconnect(p_scr *s)
{
  int i, face, index;
  if (p_on_connect) p_on_connect(1, -1);
  if (s->lockedView)
  { [s->lockedView unlockFocus];
    s->lockedView = NULL;
  }
  CGColorSpaceRelease(s->colorspace);
  for (face = 0; face < 5; face++)
    for (index = 0; index < 4; index++)
      if (s->fonts[face][index]) CGFontRelease(s->fonts[face][index]);
  CGFontRelease(s->font);
  CGContextRelease(s->scratch);
  for (i=0 ; i< P_NONE; i++) {
    if (s->cursors[i]) [s->cursors[i] release];
  }
}

void
p_flush(p_win *pw)
{
  p_scr* s = pw->s;
  if (s->lockedView)
  { [s->lockedView unlockFocus];
    s->lockedView = NULL;
    CGContextFlush(pw->cr);
  }
}

void
p_clear(p_win *pw)
{
  View* view = pw->view;
  if (!p_signalling)
  { if (view) {
      NSRect r = pw->w ? [pw->view bounds] : [pw->parent->view bounds];
      CGContextRef cr = pw->cr;
      CGContextSetFillColor(cr, pw->s->sys_colors[255-P_BG]);
      CGContextFillRect(cr, *(CGRect*)(&r));
    }
  }
}

void
p_abort(void)
{
  if (!p_signalling) p_signalling = PSIG_SOFT;
  if (m_abort_hook) m_abort_hook();
}

void
p_xhandler(void (*abort_hook)(void),
           void (*on_exception)(int signal, char *errmsg))
{
  m_abort_hook = abort_hook;   /* replaces p_abort */
  m_exception = on_exception;  /* when p_signalling detected */
}

void
p_wait_while(int *flag)
{
  NSEvent* event;
  if (!m_checksig())
  { while (*flag)
    { event = [NSApp nextEventMatchingMask: NSAnyEventMask
                                 untilDate: [NSDate distantPast]
                                    inMode: NSDefaultRunLoopMode
                                   dequeue: YES];
      [NSApp sendEvent: event];
      if (m_checksig()) break;
    }
  }
}

static int
m_checksig(void)
{
  int sig = p_signalling;
  if (sig) {
    p_signalling = 0;
    if (m_exception) m_exception(sig, (char *)0);
  }
  return sig;
}
