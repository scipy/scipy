/*
 * DISPAT.C
 *
 * $Id$
 *
 * Implement generic dispatcher routines for event driven programs
 *
 * The selection mask macros were taken from the MIT X11 R4 source
 * and bear the following copyright notice:
 *
 * Copyright 1989 Massachusetts Institute of Technology
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose and without fee is hereby granted, provided
 * that the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of M.I.T. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.  M.I.T. makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "dispat.h"

#include "osys.h"

extern void *(*GmMalloc)(long);
extern void (*GmFree)(void *);

/* ------------------------------------------------------------------------ */
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
#ifndef NO_SYS_PARAM_H
#include <sys/param.h> /* needed for XConnDis.c */
#endif
#ifdef NOFILE
#define MSKBITCNT NOFILE
#else
/* go ahead and try to get it with the other convention */
#define MSKBITCNT NOFILES_MAX
#endif
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

#undef NULL
#define NULL ((void *)0)

/* ------------------------------------------------------------------------ */

typedef struct Dispatcher Dispatcher;
struct Dispatcher {
  Dispatcher *next;
  int fd;                    /* UNIX file descriptor for this device */
  void *context;         /* argument passed to dispatching functions */
  int (*Pending)(void *context);      /* called before entering loop */
  void (*Flush)(void *context);   /* called before pausing for input */
  int (*Dispatch)(void *context);  /* called when input is available,
				      return non-zero to terminate
				      DispatchEvents loop */
  int cantRead;     /* set to indicate that read would cause SIGTTIN */
};

static int DoDispatch(int wait, int step);
static int DefaultError(void);

static Dispatcher *dispatcherList= 0;

int (*DispatchWorker)(void)= 0;

int (*DispatchError)(void)= 0;

/* ------------------------------------------------------------------------ */

int AddDispatcher(int fd, void *context, int (*Pending)(void *),
		  void (*Flush)(void *), int (*Dispatch)(void *))
{
  Dispatcher *dispat= dispatcherList;
  if (fd>=MSKBITCNT || fd<0) return 1;  /* impossible file descriptor */
  while (dispat) {
    if (dispat->fd==fd) break;
    dispat= dispat->next;
  }
  if (!dispat) {
    dispat= (Dispatcher *)GmMalloc(sizeof(Dispatcher));
    if (!dispat) return 2;                    /* memory manger failed */
  }
  dispat->next= dispatcherList;
  dispat->fd= fd;
  dispat->context= context;
  dispat->Pending= Pending;
  dispat->Flush= Flush;
  dispat->Dispatch= Dispatch;
  dispat->cantRead= 0;
  dispatcherList= dispat;
  return 0;                                    /* return successfully */
}

static Dispatcher *nextDispatcher;   /* see DispatchEvents */

void *RemoveDispatcher(int fd)
{
  void *context= 0;
  Dispatcher *prev= 0;
  Dispatcher *dispat= dispatcherList;
  while (dispat) {
    if (dispat->fd==fd) break;
    prev= dispat;
    dispat= dispat->next;
  }
  if (dispat) {
    if (prev) prev->next= dispat->next;
    else dispatcherList= dispat->next;
    if (dispat==nextDispatcher) nextDispatcher= dispat->next;
    context= dispat->context;
    GmFree(dispat);
  }
  return context;
}

void DispatchEvents(void) { DoDispatch(1, 0); }

void MaybeDispatch(void) { DoDispatch(0, 0); }

int DispatchNext(void) { return DoDispatch(1, 1); }

static int DoDispatch(int wait, int step)
{
  Dispatcher *dispat;
  int i, nReady, nReadable, lastPass, maxfd;
  unsigned long mask[MSKCNT];

  /* Handle pending events before entering perpetual loop */
  dispat= dispatcherList;
  while (dispat) {
    if (dispat->Pending && dispat->Pending(dispat->context)>0) return -3;
    dispat= dispat->next;
  }

  lastPass= 0;
  while (dispatcherList) {

    /* Flush output */
    for (dispat=dispatcherList ; dispat ; dispat=dispat->next) {
      if (dispat->Flush) dispat->Flush(dispat->context);
    }
    for (i=0 ; i<MSKCNT ; i++) mask[i]= 0;

    /* If there is a worker routine, invoke it after polling for input */
    if (DispatchWorker) {
      do {
	maxfd= -1;
	for (dispat=dispatcherList ; dispat ; dispat=dispat->next) {
	  BITSET(mask, dispat->fd);
	  if (dispat->fd>maxfd) maxfd= dispat->fd;
	}
	nReady= G_poll(maxfd, mask, 0L);
      } while (!nReady && DispatchWorker());
    } else nReady= 0;

    /* Pause for new input */
    if (!nReady) {
	maxfd= -1;
	for (dispat=dispatcherList ; dispat ; dispat=dispat->next) {
	  BITSET(mask, dispat->fd);
	  if (dispat->fd>maxfd) maxfd= dispat->fd;
	}
	nReady= G_poll(maxfd, mask, wait? -1L : 0L);
    }
    if (nReady<0) {
      /* A serious error has occurred.  */
      if (!DispatchError) DefaultError();
      else if (DispatchError()) break;
    }

    /* Dispatch events to devices for which input has arrived */
    nReadable= 0;             /* count number of readable devices */
    dispat= dispatcherList;
    while (dispat) {
      /* allow RemoveDispatcher to function properly if called by a
	 Dispatch routine */
      nextDispatcher= dispat->next;
      /* check bit mask returned by select */
      if (GETBIT(mask, dispat->fd) && dispat->Dispatch) {
	int status= dispat->Dispatch(dispat->context);
	if (status>0) return -2;
	if (status<0) dispat->cantRead= 1;
	else if (step) return dispat->fd;
	else nReadable++;
      } else if (!dispat->cantRead) nReadable++;
      dispat= nextDispatcher;
    }

    if (!nReadable) {
      /* cantRead is set on all devices, try one more pass... */
      if (lastPass) break;
      lastPass= 1;
    } else {
      lastPass= 0;
      if (!wait && !nReady) break;
    }
  }

  return -1;
}

int HasDispatcher(int fd)
{
  Dispatcher *dispat= dispatcherList;
  while (dispat) {
    if (dispat->fd==fd) return 1;
    dispat= dispat->next;
  }
  return 0;
}

#undef NULL
#include <stdio.h>

static int DefaultError(void)
{
  fprintf(stderr, "\n\n**WARNING** select failed in DispatchEvents\n\n");
  return 0;
}
