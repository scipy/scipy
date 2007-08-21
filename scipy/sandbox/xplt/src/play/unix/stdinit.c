/*
 * stdinit.c -- $Id$
 * UNIX version of play terminal I/O
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "play.h"
#include "pstdlib.h"
#include "pstdio.h"
#include "playu.h"

#include <stdio.h>
#include <string.h>

static void u_fd0_ready(void *);

static char *u_input_line = 0;
static int u_input_max = 0;
static int u_term_errs = 0;

static void (*u_stdin)(char *input_line);

void
p_stdinit(void (*on_stdin)(char *input_line))
{
  u_stdin = on_stdin;
  u_event_src(0, &u_fd0_ready, (void *)0);
}

void
p_stdout(char *output_line)
{
  if (!p_signalling) {
    fputs(output_line, stdout);
    fflush(stdout);
  }
}

void
p_stderr(char *output_line)
{
  if (!p_signalling) {
    fputs(output_line, stderr);
    fflush(stderr);
  }
}

/* ARGSUSED */
static void
u_fd0_ready(void *c)
{
  char *line = u_input_line;
  int n;

  /* before calling fgets, check to be sure that this process is
   * not in the background (via UNIX shell job control) */
  if (u_in_background()) return;

  do {
    if (u_input_max < (line-u_input_line)+1024) {
      if (u_input_max > 16000) break;
      n = line - u_input_line;
      u_input_line = p_realloc(u_input_line, u_input_max+1024);
      u_input_max += 1024;
      line = u_input_line + n;
    }
    if (!fgets(line, 1024, stdin)) {
      int at_eof = feof(stdin);
      clearerr(stdin);
      if (++u_term_errs>3 || at_eof) {
        /* cannot read stdin -- maybe serious error, remove it */
        u_event_src(0, (void (*)(void *))0, (void *)0);
        return;
      }
    } else {
      u_term_errs = 0;  /* reset error counter on each successful read */
    }
    n = strlen(line);
    line += n;
  } while (n==1023 && line[-1]!='\n');
  if (line[-1]=='\n') line[-1] = '\0';

  if (u_stdin) u_stdin(u_input_line);
  else p_stderr("\a");  /* beep to indicate rejection */

  if (u_input_max>1024) {
    u_input_line = p_realloc(u_input_line, 1024);
    u_input_max = 1024;
  }
}
