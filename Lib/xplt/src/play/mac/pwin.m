/*
 * pwin.m
 * routines to create graphics devices for Mac OS X.
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playm.h"
#include "pstdlib.h"

volatile int p_signalling = 0;

static void (*mon_expose)(void *c,int *xy)= 0;
static void (*mon_destroy)(void *c)= 0;
static void (*mon_resize)(void *c,int w,int h)= 0;
static void (*mon_focus)(void *c,int in)= 0;
static void (*mon_key)(void *c,int k,int md)= 0;
static void (*mon_click)(void *c,int b,int md,int x,int y,unsigned long ms)=0;
static void (*mon_motion)(void *c,int md,int x,int y)= 0;
static void (*mon_deselect)(void *c)= 0;

static p_win *
m_pwin(void *ctx, p_scr *s, unsigned long bg);

void
p_gui(void (*on_expose)(void *c, int *xy),
      void (*on_destroy)(void *c),
      void (*on_resize)(void *c,int w,int h),
      void (*on_focus)(void *c,int in),
      void (*on_key)(void *c,int k,int md),
      void (*on_click)(void *c,int b,int md,int x,int y,
                       unsigned long ms),
      void (*on_motion)(void *c,int md,int x,int y),
      void (*on_deselect)(void *c),
      void (*on_panic)(p_scr *s))
{
  mon_expose = on_expose;
  mon_destroy = on_destroy;
  mon_resize = on_resize;
  mon_focus = on_focus;
  mon_key = on_key;
  mon_click = on_click;
  mon_motion = on_motion;
  mon_deselect = on_deselect;
}

@implementation View
- (View*) initWithFrame:(NSRect)rect window:(p_win*)w
{ mousebutton = 0;
  self->pw = w;
  return [super initWithFrame: rect];
}

- (void)expose:(id)dummy
{ void* ctx = pw->ctx;
  if (mon_expose && ctx) mon_expose(ctx, 0);
  CGContextFlush(pw->cr);
}

- (void)drawRect:(NSRect)rect
{ if (!pw) return; // offscreen window
  NSArray* mode = [NSArray arrayWithObject: NSDefaultRunLoopMode];
  [[NSRunLoop currentRunLoop] performSelector:@selector(expose:)
                                       target:(id)self
                                     argument:nil
                                       order:0
                                       modes:mode];
}

- (NSSize)windowWillResize: (NSWindow*) window toSize:(NSSize) size
{ if (mon_resize && pw->ctx)
  { int width = (int) size.width;
    int height = (int) size.height;
    mon_resize(pw->ctx, width, height);
    size.width = width;
    size.height = height;
  }
  return size;
}

- (void)windowDidResignKey:(NSNotification *)aNotification
{ if (mon_focus && pw->ctx) mon_focus(pw->ctx, 0);
}

- (void)mouseDown:(NSEvent*)event
{ int button = 1;
  if (mon_click && pw->ctx) {
    int state;
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    unsigned long time = (unsigned long) (1000.0 * [event timestamp]);
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask) /* Emulate a right-button click */
    { state |= P_BTN3;
      button = 3;
    }
    else if (modifier&NSAlternateKeyMask) /* Emulate a middle-button click */
    { state |= P_BTN2;
      button = 2;
    }
    else state |= P_BTN1;
    mousebutton |= P_BTN1;
    state ^= (1<<(button+2));  /* make consistent with X11 */
    mon_click(pw->ctx, button, state, x, y, time);
  }
}

- (void)mouseUp:(NSEvent*)event
{ int button = 1;
  if (mon_click && pw->ctx) {
    int state;
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    unsigned long time = (unsigned long) (1000.0 * [event timestamp]);
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask)
    { state &= ~(P_BTN3); /* Emulate a right-button release */
      button = 3;
    }
    else if (modifier&NSAlternateKeyMask) /* Emulate a middle-button release */
    { state &= ~(P_BTN2);
      button = 2;
    }
    else state &= ~(P_BTN1);
    mousebutton &= ~(P_BTN1);
    state ^= (1<<(button+2));  /* make consistent with X11 */
    mon_click(pw->ctx, button, state, x, y, time);
  }
}

- (void)otherMouseDown:(NSEvent*)event
{ const int button = 2;
  if (mon_click && pw->ctx) {
    int state;
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    unsigned long time = (unsigned long) (1000.0 * [event timestamp]);
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    mousebutton |= P_BTN2;
    state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSAlternateKeyMask) state |= P_ALT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask) state |= P_CONTROL;
    state ^= (1<<(button+2));  /* make consistent with X11 */
    mon_click(pw->ctx, button, state, x, y, time);
  }
}

- (void)otherMouseUp:(NSEvent*)event
{ const int button = 2;
  if (mon_click && pw->ctx) {
    int state;
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    unsigned long time = (unsigned long) (1000.0 * [event timestamp]);
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    mousebutton &= ~(P_BTN2);
    state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSAlternateKeyMask) state |= P_ALT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask) state |= P_CONTROL;
    state ^= (1<<(button+2));  /* make consistent with X11 */
    mon_click(pw->ctx, button, state, x, y, time);
  }
}

- (void)rightMouseDown:(NSEvent*)event
{ const int button = 3;
  if (mon_click && pw->ctx) {
    int state;
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    unsigned long time = (unsigned long) (1000.0 * [event timestamp]);
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    mousebutton |= P_BTN3;
    state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSAlternateKeyMask) state |= P_ALT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask) state |= P_CONTROL;
    state ^= (1<<(button+2));  /* make consistent with X11 */
    mon_click(pw->ctx, button, state, x, y, time);
  }
}

- (void)rightMouseUp:(NSEvent*)event
{ const int button = 3;
  if (mon_click && pw->ctx) {
    int state;
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    unsigned long time = (unsigned long) (1000.0 * [event timestamp]);
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    mousebutton &= ~(P_BTN3);
    state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSAlternateKeyMask) state |= P_ALT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask) state |= P_CONTROL;
    state ^= (1<<(button+2));  /* make consistent with X11 */
    mon_click(pw->ctx, button, state, x, y, time);
  }
}

- (void)mouseMoved:(NSEvent*)event
{ if (mon_motion && pw->ctx) {
    NSPoint location = [self convertPoint: [event locationInWindow] fromView: nil];
    int x = (int) location.x;
    int y = (int) location.y;
    unsigned int modifier = [event modifierFlags];
    int state = mousebutton;
    if (modifier&NSShiftKeyMask) state |= P_SHIFT;
    if (modifier&NSAlternateKeyMask) state |= P_ALT;
    if (modifier&NSCommandKeyMask) state |= P_META;
    if (modifier&NSControlKeyMask) state |= P_CONTROL;
    mon_motion(pw->ctx, state, x, y);
  }
}

- (void)keyDown:(NSEvent*)event
{ if (mon_key)
  { int key = 0;
    int mods = 0;
    NSString* s = [event characters];
    int i;
    const int n = [s length];
    unsigned int modifier = [event modifierFlags];
    if (modifier&NSShiftKeyMask) mods |= P_SHIFT;
    if (modifier&NSAlternateKeyMask) mods |= P_ALT;
    if (modifier&NSCommandKeyMask) mods |= P_META;
    if (modifier&NSControlKeyMask) mods |= P_CONTROL;
    if (modifier&NSNumericPadKeyMask) mods |= P_KEYPAD;
    for (i = 0; i < n; i++)  
    { unichar c = [s characterAtIndex:i];
      switch(c)
      { case NSUpArrowFunctionKey: key = P_UP; break;
        case NSDownArrowFunctionKey: key = P_DOWN; break;
        case NSLeftArrowFunctionKey: key = P_LEFT; break;
        case NSRightArrowFunctionKey: key = P_RIGHT; break;
        case NSF1FunctionKey: key = P_F1; break;
        case NSF2FunctionKey: key = P_F2; break;
        case NSF3FunctionKey: key = P_F3; break;
        case NSF4FunctionKey: key = P_F4; break;
        case NSF5FunctionKey: key = P_F5; break;
        case NSF6FunctionKey: key = P_F6; break;
        case NSF7FunctionKey: key = P_F7; break;
        case NSF8FunctionKey: key = P_F8; break;
        case NSF9FunctionKey: key = P_F9; break;
        case NSF10FunctionKey: key = P_F10; break;
        case NSF11FunctionKey: key = P_F11; break;
        case NSF12FunctionKey: key = P_F12; break;
        case NSInsertFunctionKey: key = P_INSERT; break;
        case NSDeleteFunctionKey: key = '\177'; break;
        case NSHomeFunctionKey: key = P_HOME; break;
        case NSEndFunctionKey: key = P_END; break;
        case NSPageUpFunctionKey: key = P_PGUP; break;
        case NSPageDownFunctionKey: key = P_PGDN; break;
        default: if (n==1) key = [s cString][0]; break;
      }
      if (key) mon_key(pw->ctx, key, mods);
    }
  }
}

- (BOOL)acceptsFirstResponder {
  return YES;
}

- (BOOL)isFlipped {
  return YES;
}
@end

p_win *
p_window(p_scr *s, int width, int height, char *title,
         unsigned long bg, int hints, void *ctx)
{
  NSRect rect = NSMakeRect(20.0,20.0,width,height);
  NSWindow* window;
  unsigned int style = NSTitledWindowMask
                     | NSClosableWindowMask
                     | NSMiniaturizableWindowMask;
  if (!(hints & P_NORESIZE)) style |= NSResizableWindowMask;

  p_win *pw = m_pwin(ctx, s, bg);

  if (!pw) return 0;

  window = [ [NSWindow alloc] initWithContentRect: rect
                                        styleMask: style
                                          backing: NSBackingStoreBuffered
                                            defer: NO];
  if (!window)
  { p_destroy(pw);
    return NULL;
  }
  View* view = [[View alloc] initWithFrame: rect window: pw];
  [window setContentView: view];
  [window setDelegate: view];
  [window setTitle: [NSString stringWithCString: title]];
  [window setBackgroundColor: [NSColor whiteColor]];
  [window setAcceptsMouseMovedEvents: YES];
  [view allocateGState];

  pw->w = window;
  pw->view = view;
  [view lockFocus];
  NSGraphicsContext* gc = [NSGraphicsContext currentContext];
  pw->cr = (CGContextRef) [gc graphicsPort];
  [view unlockFocus];

  if (hints & P_RGBMODEL) {
    p_palette(pw, p_595, 225);
  }
    
  /* Calling orderFront or makeKeyAndOrderFront causes a
   * call to the drawRect method. In drawRect, we use
   * [NSRunLoop performSelector: target: argument: order: mode:]
   * to ensure that the call to mon_expose is postponed until
   * Python returns to the event loop. Instead, if we call
   * mon_expose here, the view would be drawn before the
   * initialization has completed. */

  if (hints&P_NOKEY) [window orderFront: nil];
  else [window makeKeyAndOrderFront: nil];

  return pw;
}

static p_win *
m_pwin(void *ctx, p_scr *s, unsigned long bg)
{
  p_win *pw = p_malloc(sizeof(p_win));
  if (pw) {
    int i;
    pw->ctx = ctx;
    pw->s = s;
    pw->w = 0;
    pw->view = NULL;

    pw->parent = NULL;

    pw->color = P_FG;
    pw->components[0] = s->sys_colors[255-P_FG][0];
    pw->components[1] = s->sys_colors[255-P_FG][1];
    pw->components[2] = s->sys_colors[255-P_FG][2];
    pw->components[3] = s->sys_colors[255-P_FG][3];
    for (i=0 ; i<=P_EXTRA ; i++)
    { pw->pixels[i][0] = s->sys_colors[255-P_FG][0];
      pw->pixels[i][1] = s->sys_colors[255-P_FG][1];
      pw->pixels[i][2] = s->sys_colors[255-P_FG][2];
      pw->pixels[i][3] = s->sys_colors[255-P_FG][3];
    }
    for (; i<256 ; i++)
    { pw->pixels[i][0] = s->sys_colors[255-i][0];
      pw->pixels[i][1] = s->sys_colors[255-i][1];
      pw->pixels[i][2] = s->sys_colors[255-i][2];
      pw->pixels[i][3] = s->sys_colors[255-i][3];
    }
    pw->n_pixels = 0;

    pw->color = P_FG;
    pw->bg = bg;

    if (p_signalling) {
      p_abort();
      p_destroy(pw);
    }
  }
  return pw;
}

void
p_destroy(p_win *pw)
{
  p_scr* s = pw->s;
  View* view = pw->view;
  NSWindow* window = pw->w;
  [view releaseGState];
  if (view)
  { if (s && s->lockedView==view)
    { [view unlockFocus];
      s->lockedView = NULL;
    }
    [view release];
  }
  if (window)
  { [window close];
    [window release];
  }
  pw->ctx = 0;           /* do not call on_destroy handler */
  pw->w = 0;
  pw->view = 0;
  p_free(pw);
}

p_win *
p_offscreen(p_win *parent, int width, int height)
{
  p_win *pw = NULL;
  if (p_signalling) {
    p_abort();
  }
  pw = m_pwin(0, parent->s, parent->bg);
  NSRect rect = NSMakeRect(0,0,width,height);
  pw->view = [[View alloc] initWithFrame: rect];
  [pw->view allocateGState];
  pw->w = [[NSWindow alloc] initWithContentRect: rect
                                      styleMask: NSBorderlessWindowMask
                                        backing: NSBackingStoreRetained
                                          defer: NO];
  [pw->w setContentView: pw->view];
  [pw->view lockFocus];
  NSGraphicsContext* gc = [NSGraphicsContext currentContext];
  pw->cr = (CGContextRef) [gc graphicsPort];
  [pw->view unlockFocus];
  pw->parent = parent;
  return pw;
}
