/* Copyright (C) 2007, 2008, 2009 Free Software Foundation, Inc.
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

#ifndef	_SYS_EVENTFD_H
#define	_SYS_EVENTFD_H	1

#include <stdint.h>


/* Type for event counter.  */
typedef uint64_t eventfd_t;

/* Flags for signalfd.  */
enum
  {
    EFD_SEMAPHORE = 1,
#define EFD_SEMAPHORE EFD_SEMAPHORE
    EFD_CLOEXEC = 02000000,
#define EFD_CLOEXEC EFD_CLOEXEC
    EFD_NONBLOCK = 04000
#define EFD_NONBLOCK EFD_NONBLOCK
  };


__BEGIN_DECLS

/* Return file descriptor for generic event channel.  Set initial
   value to COUNT.  */
extern int eventfd (int __count, int __flags) __THROW;

/* Read event counter and possibly wait for events.  */
extern int eventfd_read (int __fd, eventfd_t *__value);

/* Increment event counter.  */
extern int eventfd_write (int __fd, eventfd_t __value);

__END_DECLS

#endif /* sys/eventfd.h */
