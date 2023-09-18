/* Header file for monetary value formatting functions.
   Copyright (C) 1996,1997,1998,1999,2000,2002,2006,2009
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

#ifndef	_MONETARY_H
#define	_MONETARY_H	1

#include <features.h>

/* Get needed types.  */
#define __need_size_t
#include <stddef.h>
#include <bits/types.h>

#ifndef	__ssize_t_defined
typedef __ssize_t ssize_t;
# define __ssize_t_defined
#endif


__BEGIN_DECLS

/* Formatting a monetary value according to the current locale.  */
extern ssize_t strfmon (char *__restrict __s, size_t __maxsize,
			__const char *__restrict __format, ...)
     __THROW __attribute_format_strfmon__ (3, 4);

#ifdef __USE_XOPEN2K8
# include <xlocale.h>

/* Formatting a monetary value according to the current locale.  */
extern ssize_t strfmon_l (char *__restrict __s, size_t __maxsize,
			  __locale_t __loc,
			  __const char *__restrict __format, ...)
     __THROW __attribute_format_strfmon__ (4, 5);
#endif

#ifdef __LDBL_COMPAT
# include <bits/monetary-ldbl.h>
#endif

__END_DECLS

#endif	/* monetary.h */
