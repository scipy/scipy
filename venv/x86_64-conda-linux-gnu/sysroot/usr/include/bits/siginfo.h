/* siginfo_t, sigevent and constants.  Linux version.
   Copyright (C) 1997-2002, 2003 Free Software Foundation, Inc.
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

#if !defined _SIGNAL_H && !defined __need_siginfo_t \
    && !defined __need_sigevent_t
# error "Never include this file directly.  Use <signal.h> instead"
#endif

#include <bits/wordsize.h>

#if (!defined __have_sigval_t \
     && (defined _SIGNAL_H || defined __need_siginfo_t \
	 || defined __need_sigevent_t))
# define __have_sigval_t	1

/* Type for data associated with a signal.  */
typedef union sigval
  {
    int sival_int;
    void *sival_ptr;
  } sigval_t;
#endif

#if (!defined __have_siginfo_t \
     && (defined _SIGNAL_H || defined __need_siginfo_t))
# define __have_siginfo_t	1

# define __SI_MAX_SIZE     128
# if __WORDSIZE == 64
#  define __SI_PAD_SIZE     ((__SI_MAX_SIZE / sizeof (int)) - 4)
# else
#  define __SI_PAD_SIZE     ((__SI_MAX_SIZE / sizeof (int)) - 3)
# endif

typedef struct siginfo
  {
    int si_signo;		/* Signal number.  */
    int si_errno;		/* If non-zero, an errno value associated with
				   this signal, as defined in <errno.h>.  */
    int si_code;		/* Signal code.  */

    union
      {
	int _pad[__SI_PAD_SIZE];

	 /* kill().  */
	struct
	  {
	    __pid_t si_pid;	/* Sending process ID.  */
	    __uid_t si_uid;	/* Real user ID of sending process.  */
	  } _kill;

	/* POSIX.1b timers.  */
	struct
	  {
	    int si_tid;		/* Timer ID.  */
	    int si_overrun;	/* Overrun count.  */
	    sigval_t si_sigval;	/* Signal value.  */
	  } _timer;

	/* POSIX.1b signals.  */
	struct
	  {
	    __pid_t si_pid;	/* Sending process ID.  */
	    __uid_t si_uid;	/* Real user ID of sending process.  */
	    sigval_t si_sigval;	/* Signal value.  */
	  } _rt;

	/* SIGCHLD.  */
	struct
	  {
	    __pid_t si_pid;	/* Which child.  */
	    __uid_t si_uid;	/* Real user ID of sending process.  */
	    int si_status;	/* Exit value or signal.  */
	    __clock_t si_utime;
	    __clock_t si_stime;
	  } _sigchld;

	/* SIGILL, SIGFPE, SIGSEGV, SIGBUS.  */
	struct
	  {
	    void *si_addr;	/* Faulting insn/memory ref.  */
	  } _sigfault;

	/* SIGPOLL.  */
	struct
	  {
	    long int si_band;	/* Band event for SIGPOLL.  */
	    int si_fd;
	  } _sigpoll;
      } _sifields;
  } siginfo_t;


/* X/Open requires some more fields with fixed names.  */
# define si_pid		_sifields._kill.si_pid
# define si_uid		_sifields._kill.si_uid
# define si_timerid	_sifields._timer.si_tid
# define si_overrun	_sifields._timer.si_overrun
# define si_status	_sifields._sigchld.si_status
# define si_utime	_sifields._sigchld.si_utime
# define si_stime	_sifields._sigchld.si_stime
# define si_value	_sifields._rt.si_sigval
# define si_int		_sifields._rt.si_sigval.sival_int
# define si_ptr		_sifields._rt.si_sigval.sival_ptr
# define si_addr	_sifields._sigfault.si_addr
# define si_band	_sifields._sigpoll.si_band
# define si_fd		_sifields._sigpoll.si_fd


/* Values for `si_code'.  Positive values are reserved for kernel-generated
   signals.  */
enum
{
  SI_ASYNCNL = -60,		/* Sent by asynch name lookup completion.  */
# define SI_ASYNCNL	SI_ASYNCNL
  SI_TKILL = -6,		/* Sent by tkill.  */
# define SI_TKILL	SI_TKILL
  SI_SIGIO,			/* Sent by queued SIGIO. */
# define SI_SIGIO	SI_SIGIO
  SI_ASYNCIO,			/* Sent by AIO completion.  */
# define SI_ASYNCIO	SI_ASYNCIO
  SI_MESGQ,			/* Sent by real time mesq state change.  */
# define SI_MESGQ	SI_MESGQ
  SI_TIMER,			/* Sent by timer expiration.  */
# define SI_TIMER	SI_TIMER
  SI_QUEUE,			/* Sent by sigqueue.  */
# define SI_QUEUE	SI_QUEUE
  SI_USER,			/* Sent by kill, sigsend, raise.  */
# define SI_USER	SI_USER
  SI_KERNEL = 0x80		/* Send by kernel.  */
#define SI_KERNEL	SI_KERNEL
};


/* `si_code' values for SIGILL signal.  */
enum
{
  ILL_ILLOPC = 1,		/* Illegal opcode.  */
# define ILL_ILLOPC	ILL_ILLOPC
  ILL_ILLOPN,			/* Illegal operand.  */
# define ILL_ILLOPN	ILL_ILLOPN
  ILL_ILLADR,			/* Illegal addressing mode.  */
# define ILL_ILLADR	ILL_ILLADR
  ILL_ILLTRP,			/* Illegal trap. */
# define ILL_ILLTRP	ILL_ILLTRP
  ILL_PRVOPC,			/* Privileged opcode.  */
# define ILL_PRVOPC	ILL_PRVOPC
  ILL_PRVREG,			/* Privileged register.  */
# define ILL_PRVREG	ILL_PRVREG
  ILL_COPROC,			/* Coprocessor error.  */
# define ILL_COPROC	ILL_COPROC
  ILL_BADSTK			/* Internal stack error.  */
# define ILL_BADSTK	ILL_BADSTK
};

/* `si_code' values for SIGFPE signal.  */
enum
{
  FPE_INTDIV = 1,		/* Integer divide by zero.  */
# define FPE_INTDIV	FPE_INTDIV
  FPE_INTOVF,			/* Integer overflow.  */
# define FPE_INTOVF	FPE_INTOVF
  FPE_FLTDIV,			/* Floating point divide by zero.  */
# define FPE_FLTDIV	FPE_FLTDIV
  FPE_FLTOVF,			/* Floating point overflow.  */
# define FPE_FLTOVF	FPE_FLTOVF
  FPE_FLTUND,			/* Floating point underflow.  */
# define FPE_FLTUND	FPE_FLTUND
  FPE_FLTRES,			/* Floating point inexact result.  */
# define FPE_FLTRES	FPE_FLTRES
  FPE_FLTINV,			/* Floating point invalid operation.  */
# define FPE_FLTINV	FPE_FLTINV
  FPE_FLTSUB			/* Subscript out of range.  */
# define FPE_FLTSUB	FPE_FLTSUB
};

/* `si_code' values for SIGSEGV signal.  */
enum
{
  SEGV_MAPERR = 1,		/* Address not mapped to object.  */
# define SEGV_MAPERR	SEGV_MAPERR
  SEGV_ACCERR			/* Invalid permissions for mapped object.  */
# define SEGV_ACCERR	SEGV_ACCERR
};

/* `si_code' values for SIGBUS signal.  */
enum
{
  BUS_ADRALN = 1,		/* Invalid address alignment.  */
# define BUS_ADRALN	BUS_ADRALN
  BUS_ADRERR,			/* Non-existant physical address.  */
# define BUS_ADRERR	BUS_ADRERR
  BUS_OBJERR			/* Object specific hardware error.  */
# define BUS_OBJERR	BUS_OBJERR
};

/* `si_code' values for SIGTRAP signal.  */
enum
{
  TRAP_BRKPT = 1,		/* Process breakpoint.  */
# define TRAP_BRKPT	TRAP_BRKPT
  TRAP_TRACE			/* Process trace trap.  */
# define TRAP_TRACE	TRAP_TRACE
};

/* `si_code' values for SIGCHLD signal.  */
enum
{
  CLD_EXITED = 1,		/* Child has exited.  */
# define CLD_EXITED	CLD_EXITED
  CLD_KILLED,			/* Child was killed.  */
# define CLD_KILLED	CLD_KILLED
  CLD_DUMPED,			/* Child terminated abnormally.  */
# define CLD_DUMPED	CLD_DUMPED
  CLD_TRAPPED,			/* Traced child has trapped.  */
# define CLD_TRAPPED	CLD_TRAPPED
  CLD_STOPPED,			/* Child has stopped.  */
# define CLD_STOPPED	CLD_STOPPED
  CLD_CONTINUED			/* Stopped child has continued.  */
# define CLD_CONTINUED	CLD_CONTINUED
};

/* `si_code' values for SIGPOLL signal.  */
enum
{
  POLL_IN = 1,			/* Data input available.  */
# define POLL_IN	POLL_IN
  POLL_OUT,			/* Output buffers available.  */
# define POLL_OUT	POLL_OUT
  POLL_MSG,			/* Input message available.   */
# define POLL_MSG	POLL_MSG
  POLL_ERR,			/* I/O error.  */
# define POLL_ERR	POLL_ERR
  POLL_PRI,			/* High priority input available.  */
# define POLL_PRI	POLL_PRI
  POLL_HUP			/* Device disconnected.  */
# define POLL_HUP	POLL_HUP
};

# undef __need_siginfo_t
#endif	/* !have siginfo_t && (have _SIGNAL_H || need siginfo_t).  */


#if (defined _SIGNAL_H || defined __need_sigevent_t) \
    && !defined __have_sigevent_t
# define __have_sigevent_t	1

/* Structure to transport application-defined values with signals.  */
# define __SIGEV_MAX_SIZE	64
# if __WORDSIZE == 64
#  define __SIGEV_PAD_SIZE	((__SIGEV_MAX_SIZE / sizeof (int)) - 4)
# else
#  define __SIGEV_PAD_SIZE	((__SIGEV_MAX_SIZE / sizeof (int)) - 3)
# endif

typedef struct sigevent
  {
    sigval_t sigev_value;
    int sigev_signo;
    int sigev_notify;

    union
      {
	int _pad[__SIGEV_PAD_SIZE];

	/* When SIGEV_SIGNAL and SIGEV_THREAD_ID set, LWP ID of the
	   thread to receive the signal.  */
	__pid_t _tid;

	struct
	  {
	    void (*_function) (sigval_t);	/* Function to start.  */
	    void *_attribute;			/* Really pthread_attr_t.  */
	  } _sigev_thread;
      } _sigev_un;
  } sigevent_t;

/* POSIX names to access some of the members.  */
# define sigev_notify_function   _sigev_un._sigev_thread._function
# define sigev_notify_attributes _sigev_un._sigev_thread._attribute

/* `sigev_notify' values.  */
enum
{
  SIGEV_SIGNAL = 0,		/* Notify via signal.  */
# define SIGEV_SIGNAL	SIGEV_SIGNAL
  SIGEV_NONE,			/* Other notification: meaningless.  */
# define SIGEV_NONE	SIGEV_NONE
  SIGEV_THREAD,			/* Deliver via thread creation.  */
# define SIGEV_THREAD	SIGEV_THREAD

  SIGEV_THREAD_ID = 4		/* Send signal to specific thread.  */
#define SIGEV_THREAD_ID	SIGEV_THREAD_ID
};

#endif	/* have _SIGNAL_H.  */
