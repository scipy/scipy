/* Copyright (C) 2002 Free Software Foundation, Inc.
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

#ifndef _SYS_XATTR_H
#define _SYS_XATTR_H	1

#include <features.h>
#include <sys/types.h>


__BEGIN_DECLS

/* The following constants should be used for the fifth parameter of
   `*setxattr'.  */
enum
{
  XATTR_CREATE = 1,	/* set value, fail if attr already exists.  */
#define XATTR_CREATE	XATTR_CREATE
  XATTR_REPLACE = 2	/* set value, fail if attr does not exist.  */
#define XATTR_REPLACE	XATTR_REPLACE
};

/* Set the attribute NAME of the file pointed to by PATH to VALUE (which
   is SIZE bytes long).  Return 0 on success, -1 for errors.  */
extern int setxattr (__const char *__path, __const char *__name,
		     __const void *__value, size_t __size, int __flags)
	__THROW;

/* Set the attribute NAME of the file pointed to by PATH to VALUE (which is
   SIZE bytes long), not following symlinks for the last pathname component.
   Return 0 on success, -1 for errors.  */
extern int lsetxattr (__const char *__path, __const char *__name,
		      __const void *__value, size_t __size, int __flags)
	__THROW;

/* Set the attribute NAME of the file descriptor FD to VALUE (which is SIZE
   bytes long).  Return 0 on success, -1 for errors.  */
extern int fsetxattr (int __fd, __const char *__name, __const void *__value,
		      size_t __size, int __flags) __THROW;

/* Get the attribute NAME of the file pointed to by PATH to VALUE (which is
   SIZE bytes long).  Return 0 on success, -1 for errors.  */
extern ssize_t getxattr (__const char *__path, __const char *__name,
			 void *__value, size_t __size) __THROW;

/* Get the attribute NAME of the file pointed to by PATH to VALUE (which is
   SIZE bytes long), not following symlinks for the last pathname component.
   Return 0 on success, -1 for errors.  */
extern ssize_t lgetxattr (__const char *__path, __const char *__name,
			  void *__value, size_t __size) __THROW;

/* Get the attribute NAME of the file descriptor FD to VALUE (which is SIZE
   bytes long).  Return 0 on success, -1 for errors.  */
extern ssize_t fgetxattr (int __fd, __const char *__name, void *__value,
			  size_t __size) __THROW;

/* List attributes of the file pointed to by PATH into the user-supplied
   buffer LIST (which is SIZE bytes big).  Return 0 on success, -1 for
   errors.  */
extern ssize_t listxattr (__const char *__path, char *__list, size_t __size)
	__THROW;

/* List attributes of the file pointed to by PATH into the user-supplied
   buffer LIST (which is SIZE bytes big), not following symlinks for the
   last pathname component.  Return 0 on success, -1 for errors.  */
extern ssize_t llistxattr (__const char *__path, char *__list, size_t __size)
	__THROW;

/* List attributes of the file descriptor FD into the user-supplied buffer
   LIST (which is SIZE bytes big).  Return 0 on success, -1 for errors.  */
extern ssize_t flistxattr (int __fd, char *__list, size_t __size)
	__THROW;

/* Remove the attribute NAME from the file pointed to by PATH.  Return 0
   on success, -1 for errors.  */
extern int removexattr (__const char *__path, __const char *__name) __THROW;

/* Remove the attribute NAME from the file pointed to by PATH, not
   following symlinks for the last pathname component.  Return 0 on
   success, -1 for errors.  */
extern int lremovexattr (__const char *__path, __const char *__name) __THROW;

/* Remove the attribute NAME from the file descriptor FD.  Return 0 on
   success, -1 for errors.  */
extern int fremovexattr (int __fd, __const char *__name) __THROW;

__END_DECLS

#endif	/* sys/xattr.h  */
