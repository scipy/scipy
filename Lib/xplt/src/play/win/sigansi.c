/*
 * sigansi.c -- $Id$
 * signal handing using POSIX/ANSI signals (GNU/cygwin)
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

#include <signal.h>
#include <setjmp.h>

static int w_quitting = 0;

static jmp_buf w_mainloop;

static void sig_any(int sig);
static void sig_fpe(int sig);
static void sig_int(int sig);
static int w_catch_count = 0;

static int w_sigdbg = 0xffff;

int
w_protect(int (*run)(void))
{
  int result = 0;
  if (setjmp(w_mainloop)) {
    w_catch_count++;
    w_caught();
    w_catch_count = 0;
  }
  if (!w_quitting && w_catch_count<6)
    result = run();
  w_quitting = 1;
  return result;
}

void (*w_abort_hook)(void) = 0;

void
p_abort(void)
{
  if (!p_signalling) p_signalling = PSIG_SOFT;
  if (w_abort_hook) w_abort_hook();
  /* Microsoft documentation warns that msvcrt.dll
   * may not be POSIX compliant for using longjmp out of an
   * interrupt handling routine -- oh well */
  if (!w_quitting)
    longjmp(w_mainloop, 1);
}

void
w_siginit(void)
{
  if (w_sigdbg&1) {
    signal(SIGFPE, &sig_fpe);
    w_fpu_setup();
  }
  /* SIGINT handled by SetConsoleCtrlHandler */
  if (w_sigdbg&4) signal(SIGSEGV, &sig_any);
#ifdef SIGBUS
  if (w_sigdbg&8) signal(SIGBUS, &sig_any);    /* not POSIX */
# define MY_SIGBUS SIGBUS
#else
# define MY_SIGBUS 0
#endif
  if (w_sigdbg&16) signal(SIGILL, &sig_any);
#ifdef SIGBUS
  if (w_sigdbg&32) signal(SIGPIPE, &sig_any);  /* not ANSI C */
# define MY_SIGPIPE SIGPIPE
#else
# define MY_SIGPIPE 0
#endif
}

static int sig_table[] = {
  0, 0, SIGINT, SIGFPE, SIGSEGV, SIGILL, MY_SIGBUS, MY_SIGPIPE };

static void
sig_any(int sig)
{
  int i;
  for (i=1 ; i<PSIG_OTHER ; i++) if (sig==sig_table[i]) break;
  p_signalling = i;
  if (!w_quitting) signal(sig, &sig_any);
  p_abort();
}

static void
sig_fpe(int sig)
{
  p_signalling = PSIG_FPE;
  if (!w_quitting) signal(SIGFPE, &sig_fpe);
  w_fpu_setup();
  p_abort();
}
