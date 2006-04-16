/*
 * ugetc.h -- $Id$
 * play interface for non-event-driven programs
 * -- incompatible with p_stdinit/p_stdout functions
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

/* assumes <stdio.h> has been included */

#include "extern_c.h"

/* non-play main programs may be able to use these as drop-in
 *   replacements for getc and/or fgets in order to use gist
 * - in particular, u_getc can be the readline rl_getc_function
 * - u_waitfor and u_wait_stdin return 0 if input available
 * - such programs must also define u_abort_hook as appropriate
 * - the code may call u_pending_events() to handle all pending
 *   events and expired alarms without blocking */
extern int u_getc(FILE *stream);
extern char *u_fgets(char *s, int size, FILE *stream);
extern int u_waitfor(FILE *stream);

/* call p_pending_events() to handle all pending
 *   events and expired alarms without blocking
 * call p_wait_while to set a flag and wait for play graphics events
 *   until flag is reset to zero
 * call p_xhandler to set abort_hook and on_exception routines
 *   a less obtrusive alternative to p_handler */
extern void p_pending_events(void);
extern void p_wait_while(int *flag);
extern void p_xhandler(void (*abort_hook)(void),
                       void (*on_exception)(int signal, char *errmsg));
extern int p_wait_stdin(void);

END_EXTERN_C
