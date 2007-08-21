/*
 * ugetc.c -- $Id$
 * play interface for non-event-driven programs
 * -- incompatible with stdinit.c functions
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */
#ifndef _POSIX_SOURCE
/* to get fileno declared */
#define _POSIX_SOURCE 1
#endif

#include "config.h"

#include "playu.h"
#include <stdio.h>
#include "ugetc.h"
#include "pmin.h"

static void u_fd0_ready(void *c);
static void u_nowait(void);

static FILE *u_stream = 0;
static FILE *u_fd0_init = 0;
/* ARGSUSED */
static
void u_fd0_ready(void *c)
{
  u_stream = c;
}

int
u_getc(FILE *stream)
{
  u_waitfor(stream);
  return getc(stream);
}

char *
u_fgets(char *s, int size, FILE *stream)
{
  u_waitfor(stream);
  return fgets(s, size, stream);
}

int
u_waitfor(FILE *stream)
{
  if (stream != u_fd0_init) {
    u_nowait();
    u_event_src(fileno(stream), &u_fd0_ready, stream);
    u_fd0_init = stream;
  }

  u_stream = 0;
  while (!u_stream) u_waiter(1);
  stream = u_stream;
  u_stream = 0;
  return (stream != u_fd0_init);  /* 0 on success */
}

static void
u_nowait(void)
{
  if (u_fd0_init) {
    u_event_src(fileno(u_fd0_init), (void (*)(void *c))0, u_fd0_init);
    u_fd0_init = 0;
  }
}

int
p_wait_stdin(void)
{
  return u_waitfor(stdin);
}

void
p_pending_events(void)
{
  /* do not handle events on u_fd0_init -- just everything else */
  u_nowait();
  while (u_waiter(0));
}

void
p_wait_while(int *flag)
{
  p_pending_events();
  while (*flag) u_waiter(1);
}

void
p_xhandler(void (*abort_hook)(void),
           void (*on_exception)(int signal, char *errmsg))
{
  u_abort_hook = abort_hook;   /* replaces p_abort */
  u_exception = on_exception;  /* when u_waiter detects p_signalling */
}
