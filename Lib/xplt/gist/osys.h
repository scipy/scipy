/*
 * OSYS.H
 *
 * $Id$
 *
 * Declare operating system service routines for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef OSYS_H
#define OSYS_H

/* The returned pointers point to scratch space-- the strings should
   be copied immediately.  */
extern char *GetUserName(void);
extern char *GetCurrentDate(void);

/* The G_poll function tries to be slightly more general than Gist
   requires, with a timeout value instead of just a "no wait" flag.
   Use timeout -1 to block forever, 0 to not block.
   The maxfd and mask arguments are as for the select function,
   except that maxfd is the actual maximum file descriptor, rather
   than the first file desriptor NOT to check, as in select, and
   timeout is in milliseconds if it is >0 (Gist never does this).  */
extern int G_poll(long maxfd, unsigned long *mask, long timeout);

#endif
