/*
 * handler.c -- $Id$
 * exception handling for UNIX machines
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#ifndef _XOPEN_SOURCE
/* new Linux SIGINT broken -- use sysv_signal by defining _XOPEN_SOURCE */
#define _XOPEN_SOURCE 1
#endif

#include "config.h"
#include "pstdlib.h"
#include "play.h"
#include "playu.h"

#include <signal.h>

extern void u_sigfpe(int sig);  /* may be needed by u_fpu_setup */
static void u_sigint(int sig);
static void u_sigalrm(int sig);
static void u_sigany(int sig);

static int u_sigdbg = 0xffff;

void
p_handler(void (*on_exception)(int signal, char *errmsg))
{
  u_exception = on_exception;

  if (u_sigdbg&1) signal(SIGFPE, &u_sigfpe);
  if (u_sigdbg&2) {
    signal(SIGINT, &u_sigint);
    signal(SIGALRM, &u_sigalrm);
  }
  if (u_sigdbg&4) signal(SIGSEGV, &u_sigany);
#ifdef SIGBUS
  if (u_sigdbg&8) signal(SIGBUS, &u_sigany);    /* not POSIX */
# define MY_SIGBUS SIGBUS
#else
# define MY_SIGBUS 0
#endif
  if (u_sigdbg&16) signal(SIGILL, &u_sigany);
  if (u_sigdbg&32) signal(SIGPIPE, &u_sigany);  /* not ANSI C */
}

static int sig_table[] = {
  0, 0, SIGINT, SIGFPE, SIGSEGV, SIGILL, MY_SIGBUS, SIGPIPE };

static void
u_sigany(int sig)
{
  int i;
  for (i=1 ; i<PSIG_OTHER ; i++) if (sig==sig_table[i]) break;
  p_signalling = i;
  signal(sig, &u_sigany);
  p_abort();
}

void
u_sigfpe(int sig)
{
  p_signalling = PSIG_FPE;
  u_fpu_setup(1);
  signal(SIGFPE, &u_sigfpe);
  p_abort();
}

static void
u_sigint(int sig)
{
  if (!p_signalling) {
    /* not all systems have ualarm function for microsecond resolution */
    extern unsigned int alarm(unsigned int);       /* POSIX <unistd.h> */
    p_signalling = PSIG_INT;
    alarm(1);                        /* about 1/2 sec would be better? */
  }
  signal(SIGINT, &u_sigint);
}

static void
u_sigalrm(int sig)
{
  signal(sig, &u_sigalrm);
  if (p_signalling==PSIG_INT) p_abort();
}
