/* Copyright (C) 2007, 2008 Free Software Foundation, Inc.
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

#ifndef	_SYS_SIGNALFD_H
#define	_SYS_SIGNALFD_H	1

#define __need_sigset_t
#include <signal.h>
#include <stdint.h>


struct signalfd_siginfo
{
  uint32_t ssi_signo;
  int32_t ssi_errno;
  int32_t ssi_code;
  uint32_t ssi_pid;
  uint32_t ssi_uid;
  int32_t ssi_fd;
  uint32_t ssi_tid;
  uint32_t ssi_band;
  uint32_t ssi_overrun;
  uint32_t ssi_trapno;
  int32_t ssi_status;
  int32_t ssi_int;
  uint64_t ssi_ptr;
  uint64_t ssi_utime;
  uint64_t ssi_stime;
  uint64_t ssi_addr;
  uint8_t __pad[48];
};

/* Flags for signalfd.  */
enum
  {
    SFD_CLOEXEC = 02000000,
#define SFD_CLOEXEC SFD_CLOEXEC
    SFD_NONBLOCK = 04000
#define SFD_NONBLOCK SFD_NONBLOCK
  };

__BEGIN_DECLS

/* Request notification for delivery of signals in MASK to be
   performed using descriptor FD.*/
extern int signalfd (int __fd, const sigset_t *__mask, int __flags)
  __THROW __nonnull ((2));

__END_DECLS

#endif /* sys/signalfd.h */
