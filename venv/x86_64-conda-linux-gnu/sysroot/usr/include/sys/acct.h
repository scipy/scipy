/* Copyright (C) 1996, 1997, 1998, 1999, 2007 Free Software Foundation, Inc.
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

#ifndef _SYS_ACCT_H
#define _SYS_ACCT_H	1

#include <features.h>

#include <endian.h>
#define	__need_time_t
#include <time.h>
#include <sys/types.h>

__BEGIN_DECLS

#define ACCT_COMM 16

/*
  comp_t is a 16-bit "floating" point number with a 3-bit base 8
  exponent and a 13-bit fraction. See linux/kernel/acct.c for the
  specific encoding system used.
*/

typedef u_int16_t comp_t;

struct acct
{
  char ac_flag;			/* Flags.  */
  u_int16_t ac_uid;		/* Real user ID.  */
  u_int16_t ac_gid;		/* Real group ID.  */
  u_int16_t ac_tty;		/* Controlling terminal.  */
  u_int32_t ac_btime;		/* Beginning time.  */
  comp_t ac_utime;		/* User time.  */
  comp_t ac_stime;		/* System time.  */
  comp_t ac_etime;		/* Elapsed time.  */
  comp_t ac_mem;		/* Average memory usage.  */
  comp_t ac_io;			/* Chars transferred.  */
  comp_t ac_rw;			/* Blocks read or written.  */
  comp_t ac_minflt;		/* Minor pagefaults.  */
  comp_t ac_majflt;		/* Major pagefaults.  */
  comp_t ac_swaps;		/* Number of swaps.  */
  u_int32_t ac_exitcode;	/* Process exitcode.  */
  char ac_comm[ACCT_COMM+1];	/* Command name.  */
  char ac_pad[10];		/* Padding bytes.  */
};


struct acct_v3
{
  char ac_flag;			/* Flags */
  char ac_version;		/* Always set to ACCT_VERSION */
  u_int16_t ac_tty;		/* Control Terminal */
  u_int32_t ac_exitcode;	/* Exitcode */
  u_int32_t ac_uid;		/* Real User ID */
  u_int32_t ac_gid;		/* Real Group ID */
  u_int32_t ac_pid;		/* Process ID */
  u_int32_t ac_ppid;		/* Parent Process ID */
  u_int32_t ac_btime;		/* Process Creation Time */
  float ac_etime;		/* Elapsed Time */
  comp_t ac_utime;		/* User Time */
  comp_t ac_stime;		/* System Time */
  comp_t ac_mem;		/* Average Memory Usage */
  comp_t ac_io;			/* Chars Transferred */
  comp_t ac_rw;			/* Blocks Read or Written */
  comp_t ac_minflt;		/* Minor Pagefaults */
  comp_t ac_majflt;		/* Major Pagefaults */
  comp_t ac_swaps;		/* Number of Swaps */
  char ac_comm[ACCT_COMM];	/* Command Name */
};


enum
  {
    AFORK = 0x01,		/* Has executed fork, but no exec.  */
    ASU = 0x02,			/* Used super-user privileges.  */
    ACORE = 0x08,		/* Dumped core.  */
    AXSIG = 0x10		/* Killed by a signal.  */
  };

#if __BYTE_ORDER == __BIG_ENDIAN
# define ACCT_BYTEORDER 0x80	/* Accounting file is big endian.  */
#else
# define ACCT_BYTEORDER 0x00	/* Accounting file is little endian.  */
#endif

#define AHZ     100


/* Switch process accounting on and off.  */
extern int acct (__const char *__filename) __THROW;

__END_DECLS

#endif	/* sys/acct.h */
