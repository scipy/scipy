/*
 * playu.h -- $Id$
 * UNIX-private portability layer declarations
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "extern_c.h"

extern void (*u_on_idle)(void);

/* global objects used by main() for basic program flow control
 * and ANSI C main() semantics */
extern int u_main_loop(int (*on_launch)(int,char**), int, char **);
extern int u_waiter(int wait);
extern void (*u_exception)(int signal, char *errmsg);
extern char *u_errmsg;

/* if set non-zero, p_abort will call this */
extern void (*u_abort_hook)(void);

extern void u_fpu_setup(int when);

/* arrange to call callback(context) when input or error arrives on fd
 * - callback==0 cancels */
extern void u_event_src(int fd, void (*callback)(void *), void *context);
extern void u_prepoll(int (*conditional)(void *), void *context);

/* wait at most timeout milliseconds for input on any descriptor
 * for which u_event_src has been called; timeout==-1 means forever
 * - returns 0 on timeout or signal, 1 if callback made, <0 if error
 * - at most one callback made per call to u_poll */
extern int u_poll(int timeout);

/* wait at most timeout milliseconds for input to arrive on fd only
 * - returns 1 if input arrived, 0 if timeout or signal, <0 if error
 * - no callback is made, you have to read fd after this if >0 */
extern int u_poll1(int fd, int timeout);

/* check whether this program is running in background */
extern int u_in_background(void);

/* expand ~ and $ENV_VAR in pathnames
 * - return value is in p_wkspc, be careful not to clobber */
extern char *u_pathname(const char *pathname);

/* find argv[0] on PATH, track symbolic links to actual dir entry */
extern char *u_find_exe(const char *argv0);
extern char *u_track_link(const char *name);

END_EXTERN_C
