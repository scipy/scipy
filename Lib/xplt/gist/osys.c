/*
 * OSYS.C
 *
 * $Id$
 *
 * Implement operating system service routines for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "osys.h"

/* cuserid is used to get the user name */
#include <stdio.h>

/* The time and ctime routines seem to be fairly universal...
   May need <sys/types.h>, <sys/time.h>, though ANSI is just <time.h>.  */
#include <time.h>

#ifdef HAS_CUSERID
#ifdef NEED_CUSERID
extern char *cuserid(char *);
#endif
#else
static char *cuserid(char *);
static char *cuserid(char *x) { return (char *)0; }
#endif

static char noUser[]= "";
static char noDate[]= "\n";

char *GetUserName(void)
{
  /* HP threatens to withdraw cuserid at some unspecified future date.
     The recommended replacement is much more complicated,
     invloving the getlogin function (in <unistd.h>), which is
     rumored to fail for processes (like xterm or emacs) running
     under a pty.  The fallback after failure is the complicated part.  */
  char *result= cuserid((char *)0);
  return result? result : noUser;
}

char *GetCurrentDate(void)
{
  char *result;
  time_t t= time((time_t *)0);
  if (t==-1) return noDate;
  result= ctime(&t);  /* e.g.- "Sun Jan  3 15:14:13  1988\n" */
  return result? result : noDate;
}

/* ------------------------------------------------------------------------ */
/* Gist needs to be able to block waiting for input on any of
   several I/O channels (e.g.- an X window socket or stdin).  It also
   needs to be able to check whether input is available.
   In BSD UNIX, this is the select function.
   POSIX allows select as an extension, but does not require that this
   functionality even exist.  Most SYSV operating systems seem to
   supply a working select function.
   As of March 1993, the select on SGI workstations frequently fails
   to detect input from an X socket.
   An alternative to select is the poll function, which is part of the
   streamio(4) package.
 */

/* The G_poll function tries to be slightly more general than Gist
   requires, with a timeout value instead of just a "no wait" flag.
   Use timeout -1 to block forever, 0 to not block.  */

/* errno determines if a signal caused premature return from select */
#include <errno.h>

/* ------------------------------------------------------------------------ */

#ifndef USE_POLL

/* Hopefully, the true timeval struct (in sys/time.h on Suns) will
   not be any longer than twice this guess, hence the bogus two element
   array.  Its value just needs to represent zero time.  */
struct timevalFake {
  long	tv_sec;		/* seconds */
  long	tv_usec;	/* and microseconds */
};

int G_poll(long maxfd, unsigned long *mask, long timeout)
{
  int nReady;
  long i;
  struct timevalFake zeroDelay[2];
  if (timeout>0) {
    zeroDelay[0].tv_sec= timeout/1000;
    zeroDelay[0].tv_usec= 1000*(timeout%1000);
  } else {
    zeroDelay[0].tv_sec= 0;
    zeroDelay[0].tv_usec= 0;
  }
  zeroDelay[1].tv_sec= 0;
  zeroDelay[1].tv_usec= 0;

  nReady= select(maxfd+1, mask, (void *)0, (void *)0,
		 timeout<0? (void *)0 : zeroDelay);

  if (nReady<0) {
    /* An error has occurred.  This is serious unless the "error"
       is merely that a signal has arrived before any of the selected
       events.  */
    maxfd= maxfd/sizeof(long) + 1;
    for (i=0 ; i<maxfd ; i++) mask[i]= 0;
    if (errno==EINTR) nReady= 0;
  }
  return nReady;
}

#else

#include <poll.h>

/* copy mask macros from Xlibos.h (MIT R4 distribution for Suns) */

/* May need to define "MIT_R4_SYSV" (was "att") appropriately... */
#ifdef MIT_R4_SYSV
/*
 * UNIX System V Release 3.2
 */
#include <sys/param.h>
#define MSKBITCNT NOFILES_MAX

#else
/*
 * 4.2BSD-based systems
 */
#include <sys/param.h> /* needed for XConnDis.c */
#define MSKBITCNT NOFILE
#endif

/* Does this require sizeof(unsigned long)==32?  Probably... */
#define MSKCNT ((MSKBITCNT+31) / 32)

#if (MSKCNT==1)
#define BITMASK(i) (1 << (i))
#define MASKIDX(i) 0
#endif
#if (MSKCNT>1)
/* (ditto above comment)... */
#define BITMASK(i) (1 << ((i) & 31))
#define MASKIDX(i) ((i) >> 5)
#endif

#define MASKWORD(buf, i) buf[MASKIDX(i)]
#define BITSET(buf, i) MASKWORD(buf, i) |= BITMASK(i)
#define BITCLEAR(buf, i) MASKWORD(buf, i) &= ~BITMASK(i)
#define GETBIT(buf, i) (MASKWORD(buf, i) & BITMASK(i))

int G_poll(long maxfd, unsigned long *mask, long timeout)
{
  int nReady;
  struct pollfd fds[16];  /* limited to 16 file descriptors */
  long nfds, nwords;
  int i;
  int wait= timeout;

  nfds= 0;
  for (i=0 ; i<=maxfd ; i++) {
    if (GETBIT(mask, i)) {
      if (nfds>=16) break;
      fds[nfds].fd= i;
      fds[nfds].events= POLLIN;
      fds[nfds].revents= 0;
      nfds++;
    } 
  }

  nReady= poll(fds, nfds, wait);
  if (i<=maxfd) nReady= -1;

  nwords= maxfd/sizeof(long) + 1;
  for (i=0 ; i<nwords ; i++) mask[i]= 0;
  if (nReady>=0) {
    for (i=0 ; i<nfds ; i++) {
      if (fds[i].revents&(POLLIN|POLLPRI|POLLERR))
	BITSET(mask, fds[i].fd);
    }
  } else {
    /* An error has occurred.  This is serious unless the "error"
       is merely that a signal has arrived before any of the selected
       events.  */
    if (errno==EINTR) nReady= 0;
  }
  return nReady;
}
#endif
