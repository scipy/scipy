/*
 * upoll.h -- $Id$
 *
 * UNIX poll.h header, when poll() must be implemented using select
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

struct pollfd {
  int fd;         /* file descriptor */
  short events;   /* events desired */
  short revents;  /* events returned */
};

/* events or revents */
#define POLLIN          01              /* ready to read */
#define POLLPRI         02              /* urgently ready to read */
#define POLLOUT         04              /* fd is writable */

/* revents only */
#define POLLERR         010             /* error */
#define POLLHUP         020             /* hangup */
#define POLLNVAL        040             /* fd not open */

/* timeout is in milliseconds, -1 to wait forever
 * returns number of fds with non-zero revents,
 *   or -1 and sets errno to
 *     EAGAIN, EFAULT, EINTR (signal during poll), EINVAL */
extern int u__poll(struct pollfd *fds, unsigned long int nfds, int timeout);
