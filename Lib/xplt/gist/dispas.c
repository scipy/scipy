/*
 * DISPAS.C
 *
 * $Id$
 *
 * Implement dispatcher routines for ordinary file i/o streams
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "dispas.h"

/* fileno is a macro which should be defined in <stdio.h>, but use
   a second macro here just in case something else is needed on some
   other platform.  */
#ifndef GetFD
#define GetFD(file) fileno(file)
#endif

#undef NEVER_USE
#ifdef NEVER_USE
/* Use getpgrp to determine this process's process group, and
   tcgetpgrp to determine whether this file represents the controlling
   terminal of a different process group (see below for more).  */
extern int getpgrp(int);
extern int tcgetpgrp(int);
#undef NEED_TCGETPGRP
#endif

extern void *(*GmMalloc)(long);
extern void (*GmFree)(void *);

typedef struct StreamContext StreamContext;
struct StreamContext {
  FILE *file;
  int (*Dispatch)(FILE *file, void *context);
  void *context;
};

static int DispatchS(void *context);

static int DispatchS(void *context)
{
  StreamContext *scon= context;
  FILE *file= scon->file;

#ifndef NO_TERMINALS
  int fd= GetFD(file);

  /* Before calling the user supplied dispatching routine, check
     to be sure that if this file descriptor is the controlling terminal
     for this process, the process is in the foreground.  If this is
     the controlling terminal and the process is in the background,
     any attempt to read will generate a SIGTTIN signal.  Read the
     intro(2) and termio(4) man pages for a more complete discussion.  */

  int tgid= tcgetpgrp(fd);     /* process group for this terminal--
				  could also use TIOCGPGRP ioctl */
  if (tgid>0 && !feof(file) && !ferror(file)) {
    /* This might be the tty of a background process, check */
    int pgid= getpgrp(0);      /* the System V version of getpgrp() has
				  no parameters, but this should be OK  */
    if (pgid!=tgid) return -1; /* attempting to read would cause SIGTTIN */
  }
#endif

  return scon->Dispatch(file, scon->context);
}

int AddFDispatcher(FILE *file, int (*Dispatch)(FILE *file, void *context),
		   void *context)
{
  int fd, value;
  StreamContext *scon;
  if (!file) return 1;
  scon= (StreamContext *)GmMalloc(sizeof(StreamContext));
  if (!scon) return 2;
  scon->file= file;
  scon->Dispatch= Dispatch;
  scon->context= context;
  fd= GetFD(file);
  value= AddDispatcher(fd, scon, 0, 0, &DispatchS);
  if (value) GmFree(scon);
  return value;
}

void RemoveFDispatcher(FILE *file)
{
  int fd= GetFD(file);
  StreamContext *scon= RemoveDispatcher(fd);
  if (scon) GmFree(scon);
}

#ifdef NEED_TCGETPGRP
#ifdef HAVE_SYS_TERMIOS_H
#include <sys/termios.h>
#else
#include <sgtty.h>
#endif
int tcgetpgrp(int fd)
{
  int tgid;
  if (ioctl(fd, TIOCGPGRP, &tgid)<0) return -1;
  else return tgid;
}
#endif
