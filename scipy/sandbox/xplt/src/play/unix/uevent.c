/*
 * uevent.c -- $Id$
 * UNIX input event handling
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#ifndef _HPUX_SOURCE
/* HPUX turns off poll.h without this (_INCLUDE_AES_SOURCE) */
#define _HPUX_SOURCE 1
#endif
#ifndef _XOPEN_SOURCE_EXTENDED
/* Digital UNIX needs _XOPEN_SOURCE_EXTENDED for poll.h to work */
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include "config.h"
#include "play.h"
#include "playu.h"
#include "pstdlib.h"

#ifdef TEST_POLL
#define p_realloc (void *)realloc
#endif

/* errno determines if a signal caused premature return from poll */
#include <errno.h>

#ifndef USE_SELECT
# ifndef USE_SYS_POLL_H
#  include <poll.h>
# else
#  include <sys/poll.h>
# endif
# define u__poll poll
#else
# include "upoll.h"
# include "upoll.c"
#endif

typedef struct u_handler u_handler;
struct u_handler {
  void (*callback)(void *);
  void *context;
};

static u_handler *poll_handler = 0;
static struct pollfd *poll_fds = 0;
static unsigned long poll_nfds = 0;
static unsigned long poll_maxfd = 0;

typedef struct u_prehandler u_prehandler;
struct u_prehandler {
  int (*conditional)(void *);
  void *context;
};

static u_prehandler *prepoll = 0;
static int prepoll_i = 0;
static int prepoll_n = 0;
static int prepoll_mx = 0;

void
u_event_src(int fd, void (*callback)(void *), void *context)
{
  if (callback) {
    if (poll_nfds>=poll_maxfd) {
      unsigned long nn = poll_maxfd + 4;
      poll_fds = p_realloc(poll_fds, sizeof(struct pollfd)*nn);
      poll_handler = p_realloc(poll_handler, sizeof(u_handler)*nn);
      poll_maxfd = nn;
    }
    poll_fds[poll_nfds].fd = fd;
    poll_fds[poll_nfds].events = POLLIN | POLLPRI;
    poll_fds[poll_nfds].revents = 0;
    poll_handler[poll_nfds].callback = callback;
    poll_handler[poll_nfds].context = context;
    poll_nfds++;

  } else {
    int i;
    int n = poll_nfds-1;
    for (i=0 ; i<poll_nfds ; i++)
      if (poll_fds[i].fd==fd) {
        /* CRITICAL section */
        if (n) {
          poll_handler[i].callback = poll_handler[n].callback;
          poll_handler[i].context = poll_handler[n].context;
          poll_fds[i].fd = poll_fds[n].fd;
        }
        poll_nfds = n;
        break;
      }
  }
}

void
u_prepoll(int (*conditional)(void *), void *context)
{
  if (conditional) {
    if (prepoll_n>=prepoll_mx) {
      prepoll = p_realloc(prepoll, sizeof(u_prehandler)*(prepoll_mx+4));
      prepoll_mx += 4;
    }
    prepoll[prepoll_n].conditional = conditional;
    prepoll[prepoll_n].context = context;
    prepoll_n++;

  } else {
    int i;
    int n = prepoll_n-1;
    for (i=0 ; i<prepoll_n ; i++)
      if (prepoll[i].context==context) {
        /* CRITICAL section */
        if (n) {
          prepoll[i].conditional = prepoll[n].conditional;
          prepoll[i].context = prepoll[n].context;
        }
        prepoll_n = n;
        break;
      }
  }
}

int
u_poll(int timeout)
{
  int i, n = 0;

  /* check prepolling conditionals first
   *   prepoll_i counter makes this happen in round-robin fashion,
   *   so the first one won't get all the calls */
  for (i=prepoll_n ; i-- ; ) {
    if ((++prepoll_i) >= prepoll_n) prepoll_i = 0;
    if (prepoll[prepoll_i].conditional(prepoll[prepoll_i].context))
      return 1;
  }

  /* refuse to wait forever with no event sources */
  if (!poll_nfds && timeout<0) return -3;

  /* check for any events which arrived on previous poll before
   * calling poll again -- assures that first fd in poll_fds cannot
   * prevent any other fd from being serviced */
  do {
    for (i=0 ; i<poll_nfds ; i++)
      if (poll_fds[i].revents & (POLLIN|POLLPRI|POLLERR|POLLHUP)) {
        poll_fds[i].revents = 0;
        poll_handler[i].callback(poll_handler[i].context);
        return 1;
      }
    if (n) return -2;

    if (timeout<0) timeout = -1;
    n = u__poll(poll_fds, poll_nfds, timeout);
    if (n<0 && errno==EINTR) n = 0;
  } while (n>0);

  return n;
}

/* required by x11 to wait only finite time for X selection response */
int
u_poll1(int fd, int timeout)
{
  int n;
  struct pollfd pfd;
  pfd.fd = fd;
  pfd.events = POLLIN | POLLPRI;
  pfd.revents = 0;
  n = u__poll(&pfd, (unsigned long)1, timeout);
  if (n<0 && errno!=EINTR) return n;
  return n>0;
}
