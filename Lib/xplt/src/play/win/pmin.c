/*
 * pmin.c -- $Id$
 * minimally intrusive play event handling (no stdio)
 */

#include "playw.h"
#include "pmin.h"

#ifndef CYGWIN
#include <conio.h>
#endif

static void (*w_abort_hook)(void) = 0;
static void (*w_exception)(int signal, char *errmsg) = 0;
static int w_checksig(void);

volatile int p_signalling = 0;

static int w_hook_recurse = 0;

void
p_abort(void)
{
  if (!p_signalling) p_signalling = PSIG_SOFT;
  w_hook_recurse = 0;
  w_abort_hook();   /* blow up if w_abort_hook == 0 */
}

#ifndef CYGWIN
/* apparently, it is impossible to wait for stdin in Windows
 * without the cooperation of the code which will actually read stdin
 * (despite the documentation of MsgWaitForMultipleObjects)
 * -- if you do control stdin reading code, then need to spawn
 *    a separate thread to actually read stdin as in conterm.c */
int
p_wait_stdin(void)
{
  for (;;) {
    p_pending_events();
    /* poll stdin and drop out if input has arrived
     * - this will freeze all graphics windows (and all windows events)
     *   if fgets is called after this, but should work if this routine
     *   is called back after each individual character arrives
     *   (as, for example, by readline) */
    if (_kbhit()) break;
    /* the point of this routine is to call the play on_idle routine */
    p_on_idle(0);
    /* pause 51 milliseconds or until an event arrives */
    MsgWaitForMultipleObjects(0, 0, 0, 51, QS_ALLINPUT);
  }
}
#endif

/* does not include stdin, but irrelevant for graphics events */
void
p_pending_events(void)
{
  MSG msg;
  w_checksig();
  while (PeekMessage(&msg, 0, 0,0, PM_REMOVE)) {
    if (msg.message == WM_QUIT) break;
    TranslateMessage(&msg);
    DispatchMessage(&msg);
    w_checksig();
  }
}

void
p_wait_while(int *flag)
{
  MSG msg;
  if (!w_checksig()) {
    while (*flag) {
      GetMessage(&msg, 0, 0,0);
      if (msg.message == WM_QUIT) break;
      TranslateMessage(&msg);
      DispatchMessage(&msg);
      if (w_checksig()) break;
    }
  }
}

static int
w_checksig(void)
{
  int sig = p_signalling;
  if (sig) {
    p_signalling = 0;
    if (w_exception) w_exception(sig, (char *)0);
  }
  return sig;
}

/* the WH_KEYBOARD hook sounds a little scary -- MSDN warns:
 *   Before terminating, an application must call the UnhookWindowsHookEx
 *   function to free system resources associated with the hook.
 * unclear whether this warning applies only to global hooks
 */
#if USE_KB_HOOK
static LRESULT CALLBACK w_hook(int n, WPARAM w, LPARAM l);
static HHOOK w_hhook = 0;
#endif

void
p_xhandler(void (*abort_hook)(void),
           void (*on_exception)(int signal, char *errmsg))
{
  w_abort_hook = abort_hook;   /* replaces p_abort */
  w_exception = on_exception;  /* when p_signalling detected */
#if USE_KB_HOOK
  if (!w_hhook)
    w_hhook = SetWindowsHookEx(WH_KEYBOARD, &w_hook, 0,
                               GetCurrentThreadId());
#endif
}

#if USE_KB_HOOK
static LRESULT CALLBACK
w_hook(int n, WPARAM w, LPARAM l)
{
  p_wait_stdin();
  return CallNextHookEx(w_hhook, n, w, l);
}
#endif
