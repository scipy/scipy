/*
 * umain.c -- $Id$
 * UNIX objects referenced by main.c that goes with play model
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "pstdlib.h"
#include "play.h"
#include "playu.h"

#include <setjmp.h>

void (*u_abort_hook)(void)= 0;
void (*u_exception)(int, char *)= 0;
char *u_errmsg = 0;
volatile int p_signalling = 0;

static int (*u_quitter)(void)= 0;
static int u_quitting = 0;
static int u_launched = 0;
static jmp_buf u_mainloop;
static int fault_loop = 0;

int
u_main_loop(int (*on_launch)(int,char**), int argc, char **argv)
{
  u_fpu_setup(-1);
  if (setjmp(u_mainloop)) u_fpu_setup(0);
  if (!u_quitting && !u_launched) {
    int result;
    if (argc>0 && !u_launched)
      argv[0] = p_strcpy(u_track_link(u_find_exe(argv[0])));
    u_launched = 1;
    result = on_launch(argc, argv);
    if (result) return result;
  }
  while (!u_quitting) u_waiter(1);
  p_signalling = 0;  /* ignore signals after u_quitting flag is set */
  return u_quitter? u_quitter() : 0;
}

void
p_abort(void)
{
  if (!p_signalling) p_signalling = PSIG_SOFT;
  if (u_abort_hook) u_abort_hook();
  longjmp(u_mainloop, 1);
}

void
p_quitter(int (*on_quit)(void))
{
  u_quitter = on_quit;
}

void
p_quit(void)
{
  u_quitting = 1;
}

int
u_waiter(int wait)
{
  int serviced_event = 0;

  if (p_signalling) {
    /* first priority is to catch any pending signals */
    int i = p_signalling;
    p_signalling = 0;
    if (!fault_loop && u_exception) {
      fault_loop = 1;    /* don't trust u_exception not to fault */
      u_exception(i, u_errmsg);
      serviced_event = 1;
    }
    u_errmsg = 0;

  } else {
    int have_timeout = 0;
    serviced_event = u_poll(0);   /* anything pending? */
    if (!serviced_event) {        /* if not, wait for input */
      int timeout;
      double wait_secs = p_timeout();
      have_timeout = (wait_secs>0.);
      if (wait && wait_secs) {
        do { /* int timeout may not handle > 32.767 s at once */
          if (wait_secs < 0.0)          timeout = -1;
          else if (wait_secs < 32.767)  timeout = (int)(1000.*wait_secs);
          else                          timeout = 32767;
          serviced_event = u_poll(timeout);
          if (p_signalling) return 0;
          if (serviced_event) break;
          wait_secs -= 32.767;
        } while (wait_secs > 0.0);
      }
    }
    if (serviced_event) {
      if (serviced_event==-3) p_quit();
      p_on_idle(1);
    } else {
      p_on_idle(0);
      serviced_event = have_timeout;
    }
    fault_loop = 0;
  }

  return serviced_event;
}
