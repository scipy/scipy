/* Copyright (C) 1995,1996,1997,1999,2002 Free Software Foundation, Inc.
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

#ifndef _SYS_IPC_H
#define _SYS_IPC_H	1

#include <features.h>

#if !defined __USE_SVID && !defined __USE_XOPEN && __GNUC__ >= 2
# warning "Files using this header must be compiled with _SVID_SOURCE or _XOPEN_SOURCE"
#endif

/* Get system dependent definition of `struct ipc_perm' and more.  */
#include <bits/ipctypes.h>
#include <bits/ipc.h>

#ifndef __uid_t_defined
typedef __uid_t uid_t;
# define __uid_t_defined
#endif

#ifndef __gid_t_defined
typedef __gid_t gid_t;
# define __gid_t_defined
#endif

#ifndef __mode_t_defined
typedef __mode_t mode_t;
# define __mode_t_defined
#endif

#ifndef __key_t_defined
typedef __key_t key_t;
# define __key_t_defined
#endif

__BEGIN_DECLS

/* Generates key for System V style IPC.  */
extern key_t ftok (__const char *__pathname, int __proj_id) __THROW;

__END_DECLS

#endif /* sys/ipc.h */
