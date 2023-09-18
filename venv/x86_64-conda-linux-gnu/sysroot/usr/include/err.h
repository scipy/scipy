/* 4.4BSD utility functions for error messages.
   Copyright (C) 1995,1996,1997,1998,1999,2003 Free Software Foundation, Inc.
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

#ifndef	_ERR_H
#define	_ERR_H	1

#include <features.h>

#define	__need___va_list
#include <stdarg.h>
#ifndef	__GNUC_VA_LIST
# define __gnuc_va_list	__ptr_t
#endif

__BEGIN_DECLS

/* Print "program: ", FORMAT, ": ", the standard error string for errno,
   and a newline, on stderr.  */
extern void warn (__const char *__format, ...)
     __attribute__ ((__format__ (__printf__, 1, 2)));
extern void vwarn (__const char *__format, __gnuc_va_list)
     __attribute__ ((__format__ (__printf__, 1, 0)));

/* Likewise, but without ": " and the standard error string.  */
extern void warnx (__const char *__format, ...)
     __attribute__ ((__format__ (__printf__, 1, 2)));
extern void vwarnx (__const char *__format, __gnuc_va_list)
     __attribute__ ((__format__ (__printf__, 1, 0)));

/* Likewise, and then exit with STATUS.  */
extern void err (int __status, __const char *__format, ...)
     __attribute__ ((__noreturn__, __format__ (__printf__, 2, 3)));
extern void verr (int __status, __const char *__format, __gnuc_va_list)
     __attribute__ ((__noreturn__, __format__ (__printf__, 2, 0)));
extern void errx (int __status, __const char *__format, ...)
     __attribute__ ((__noreturn__, __format__ (__printf__, 2, 3)));
extern void verrx (int __status, __const char *, __gnuc_va_list)
     __attribute__ ((__noreturn__, __format__ (__printf__, 2, 0)));

__END_DECLS

#endif	/* err.h */
