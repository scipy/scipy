/* Copyright (C) 1991-1994,1996-2001,2003,2004,2005,2007,2009,2010
   Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	POSIX Standard: 3.2.1 Wait for Process Termination	<sys/wait.h>
 */

#ifndef	_SYS_WAIT_H
#define	_SYS_WAIT_H	1

#include <features.h>

__BEGIN_DECLS

#include <signal.h>
#include <sys/resource.h>

/* These macros could also be defined in <stdlib.h>.  */
#if !defined _STDLIB_H || (!defined __USE_XOPEN && !defined __USE_XOPEN2K8)
/* This will define the `W*' macros for the flag
   bits to `waitpid', `wait3', and `wait4'.  */
# include <bits/waitflags.h>

# ifdef	__USE_BSD

/* Lots of hair to allow traditional BSD use of `union wait'
   as well as POSIX.1 use of `int' for the status word.  */

#  if defined __GNUC__ && !defined __cplusplus
#   define __WAIT_INT(status) \
  (__extension__ (((union { __typeof(status) __in; int __i; }) \
		   { .__in = (status) }).__i))
#  else
#   define __WAIT_INT(status)	(*(__const int *) &(status))
#  endif

/* This is the type of the argument to `wait'.  The funky union
   causes redeclarations with either `int *' or `union wait *' to be
   allowed without complaint.  __WAIT_STATUS_DEFN is the type used in
   the actual function definitions.  */

#  if !defined __GNUC__ || __GNUC__ < 2 || defined __cplusplus
#   define __WAIT_STATUS	void *
#   define __WAIT_STATUS_DEFN	void *
#  else
/* This works in GCC 2.6.1 and later.  */
typedef union
  {
    union wait *__uptr;
    int *__iptr;
  } __WAIT_STATUS __attribute__ ((__transparent_union__));
#   define __WAIT_STATUS_DEFN	int *
#  endif

# else /* Don't use BSD.  */

#  define __WAIT_INT(status)	(status)
#  define __WAIT_STATUS		int *
#  define __WAIT_STATUS_DEFN	int *

# endif /* Use BSD.  */

/* This will define all the `__W*' macros.  */
# include <bits/waitstatus.h>

# define WEXITSTATUS(status)	__WEXITSTATUS (__WAIT_INT (status))
# define WTERMSIG(status)	__WTERMSIG (__WAIT_INT (status))
# define WSTOPSIG(status)	__WSTOPSIG (__WAIT_INT (status))
# define WIFEXITED(status)	__WIFEXITED (__WAIT_INT (status))
# define WIFSIGNALED(status)	__WIFSIGNALED (__WAIT_INT (status))
# define WIFSTOPPED(status)	__WIFSTOPPED (__WAIT_INT (status))
# ifdef __WIFCONTINUED
#  define WIFCONTINUED(status)	__WIFCONTINUED (__WAIT_INT (status))
# endif
#endif	/* <stdlib.h> not included.  */

#ifdef	__USE_BSD
# define WCOREFLAG		__WCOREFLAG
# define WCOREDUMP(status)	__WCOREDUMP (__WAIT_INT (status))
# define W_EXITCODE(ret, sig)	__W_EXITCODE (ret, sig)
# define W_STOPCODE(sig)	__W_STOPCODE (sig)
#endif

/* The following values are used by the `waitid' function.  */
#if defined __USE_SVID || defined __USE_XOPEN || defined __USE_XOPEN2K8
typedef enum
{
  P_ALL,		/* Wait for any child.  */
  P_PID,		/* Wait for specified process.  */
  P_PGID		/* Wait for members of process group.  */
} idtype_t;
#endif


/* Wait for a child to die.  When one does, put its status in *STAT_LOC
   and return its process ID.  For errors, return (pid_t) -1.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern __pid_t wait (__WAIT_STATUS __stat_loc);

#ifdef	__USE_BSD
/* Special values for the PID argument to `waitpid' and `wait4'.  */
# define WAIT_ANY	(-1)	/* Any process.  */
# define WAIT_MYPGRP	0	/* Any process in my process group.  */
#endif

/* Wait for a child matching PID to die.
   If PID is greater than 0, match any process whose process ID is PID.
   If PID is (pid_t) -1, match any process.
   If PID is (pid_t) 0, match any process with the
   same process group as the current process.
   If PID is less than -1, match any process whose
   process group is the absolute value of PID.
   If the WNOHANG bit is set in OPTIONS, and that child
   is not already dead, return (pid_t) 0.  If successful,
   return PID and store the dead child's status in STAT_LOC.
   Return (pid_t) -1 for errors.  If the WUNTRACED bit is
   set in OPTIONS, return status for stopped children; otherwise don't.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern __pid_t waitpid (__pid_t __pid, int *__stat_loc, int __options);

#if defined __USE_SVID || defined __USE_XOPEN || defined __USE_XOPEN2K8
# define __need_siginfo_t
# include <bits/siginfo.h>
/* Wait for a childing matching IDTYPE and ID to change the status and
   place appropriate information in *INFOP.
   If IDTYPE is P_PID, match any process whose process ID is ID.
   If IDTYPE is P_PGID, match any process whose process group is ID.
   If IDTYPE is P_ALL, match any process.
   If the WNOHANG bit is set in OPTIONS, and that child
   is not already dead, clear *INFOP and return 0.  If successful, store
   exit code and status in *INFOP.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern int waitid (idtype_t __idtype, __id_t __id, siginfo_t *__infop,
		   int __options);
#endif

#if defined __USE_BSD || defined __USE_XOPEN_EXTENDED
/* This being here makes the prototypes valid whether or not
   we have already included <sys/resource.h> to define `struct rusage'.  */
struct rusage;

/* Wait for a child to exit.  When one does, put its status in *STAT_LOC and
   return its process ID.  For errors return (pid_t) -1.  If USAGE is not
   nil, store information about the child's resource usage there.  If the
   WUNTRACED bit is set in OPTIONS, return status for stopped children;
   otherwise don't.  */
extern __pid_t wait3 (__WAIT_STATUS __stat_loc, int __options,
		      struct rusage * __usage) __THROW;
#endif

#ifdef __USE_BSD
/* PID is like waitpid.  Other args are like wait3.  */
extern __pid_t wait4 (__pid_t __pid, __WAIT_STATUS __stat_loc, int __options,
		      struct rusage *__usage) __THROW;
#endif /* Use BSD.  */


__END_DECLS

#endif /* sys/wait.h  */
