/*
 * handler.c -- $Id$
 * MS Windows exception handling
 * - w_interrupt raises a signal in worker thread from boss thread,
 *   idea from MSDN "Win32 Q&A" KillThrd library
 *   by Jeffrey Richter in Microsoft Systems Journal
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

#include <process.h>
#include <float.h>

#define W_SIGINT_DELAY 1000

static UINT wm_exception = 0;
static void (*w_on_exception)(int signal, char *errmsg)= 0;

static void w_handle_exception(MSG *msg);

static int sigint_active = 0;
static HANDLE sigint_abort = 0;
static DWORD WINAPI sigint_main(LPVOID lp);
static void w_interrupt(void);

volatile int p_signalling = 0;

void
p_handler(void (*on_exception)(int signal, char *errmsg))
{
  w_on_exception = on_exception;
  wm_exception = w_add_msg(w_handle_exception);
  sigint_abort = CreateEvent(0,0,0,0);
  w_fpu_setup();
  w_siginit();
}

/* called by worker if an actual exception has been raised, or in
 * in response to wm_exception warning message */
void
w_caught(void)
{
  MSG msg;
  int sig = p_signalling;
  p_signalling = 0;
  if (sigint_active) PulseEvent(sigint_abort);
  w_on_exception(sig, (char *)0);  /* blows up if no handler */
  /* clear any pending exception warning messages off worker queue */
  while (PeekMessage(&msg, NULL, wm_exception, wm_exception,
                     PM_REMOVE));
}

/* called as response to wm_exception warning message */
/* ARGSUSED */
static void
w_handle_exception(MSG *msg)
{
  if (p_signalling) w_caught();
}

void
w_fpu_setup(void)
{
  _fpreset();
  _controlfp(_MCW_EM &
             ~(_EM_ZERODIVIDE | _EM_OVERFLOW | _EM_INVALID), _MCW_EM);
}

int
w_sigint(int delay)
{
  if (!w_on_exception) return 0;
  if (delay && !p_signalling && !sigint_active) {
    /* interrupt worker thread after W_SIGINT_DELAY */
    HANDLE h;
    UINT id;
    p_signalling = PSIG_INT;
    h = CreateThread(0,0, sigint_main, 0, 0, &id);
    if (h) {
      Sleep(0);
      CloseHandle(h);
      sigint_active = 1;
    } else {
      delay = 0;
    }
  }
  if (!delay) w_interrupt();
  return 1;
}

static DWORD WINAPI
sigint_main(LPVOID lp)
{
  PostThreadMessage(w_id_worker, wm_exception, 0, 0);
  if (WaitForSingleObject(sigint_abort, W_SIGINT_DELAY)==WAIT_TIMEOUT &&
      p_signalling==PSIG_INT)
    w_interrupt();
  sigint_active = 0;
  return 0;
}

#if defined(_X86_)
# define PC_NAME Eip
#elif defined(_MIPS_)
# define PC_NAME Fir
#elif defined(_ALPHA_)
# define PC_NAME Fir
#elif defined(_PPC_)
# define PC_NAME Iar
#else
# error need name of program counter in CONTEXT struct for this platform
#endif

static void
w_interrupt(void)
{
  /* stop the worker, change its program counter to p_abort, resume
   * -- claim is this always works under NT, but may sometimes fail
   *    under 95 if the thread is suspended in a bad place */
  SuspendThread(w_worker);
  if (WaitForSingleObject(w_worker, 0) == WAIT_TIMEOUT) {
    CONTEXT ctx;
    ctx.ContextFlags = CONTEXT_CONTROL;
    GetThreadContext(w_worker, &ctx);
    ctx.PC_NAME = (DWORD)p_abort;
    SetThreadContext(w_worker, &ctx);
    p_signalling = PSIG_INT;
    ResumeThread(w_worker);
  }
}
