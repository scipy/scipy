/* Copyright (C) 1996, 1997, 1999, 2001 Free Software Foundation, Inc.
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
 * Never include this file directly; use <limits.h> instead.
 */

/* Additional definitions from X/Open Portability Guide, Issue 4, Version 2
   System Interfaces and Headers, 4.16 <limits.h>

   Please note only the values which are not greater than the minimum
   stated in the standard document are listed.  The `sysconf' functions
   should be used to obtain the actual value.  */

#ifndef _XOPEN_LIM_H
#define _XOPEN_LIM_H	1

#define __need_IOV_MAX
#include <bits/stdio_lim.h>

/* We do not provide fixed values for

   ARG_MAX	Maximum length of argument to the `exec' function
		including environment data.

   ATEXIT_MAX	Maximum number of functions that may be registered
		with `atexit'.

   CHILD_MAX	Maximum number of simultaneous processes per real
		user ID.

   OPEN_MAX	Maximum number of files that one process can have open
		at anyone time.

   PAGESIZE
   PAGE_SIZE	Size of bytes of a page.

   PASS_MAX	Maximum number of significant bytes in a password.

   We only provide a fixed limit for

   IOV_MAX	Maximum number of `iovec' structures that one process has
		available for use with `readv' or writev'.

   if this is indeed fixed by the underlying system.
*/


/* Maximum number of `iovec' structures that one process has available
   for use with `readv' or writev'.  */
#define	_XOPEN_IOV_MAX	_POSIX_UIO_MAXIOV


/* Maximum value of `digit' in calls to the `printf' and `scanf'
   functions.  We have no limit, so return a reasonable value.  */
#define NL_ARGMAX	_POSIX_ARG_MAX

/* Maximum number of bytes in a `LANG' name.  We have no limit.  */
#define NL_LANGMAX	_POSIX2_LINE_MAX

/* Maximum message number.  We have no limit.  */
#define NL_MSGMAX	INT_MAX

/* Maximum number of bytes in N-to-1 collation mapping.  We have no
   limit.  */
#define NL_NMAX		INT_MAX

/* Maximum set number.  We have no limit.  */
#define NL_SETMAX	INT_MAX

/* Maximum number of bytes in a message.  We have no limit.  */
#define NL_TEXTMAX	INT_MAX

/* Default process priority.  */
#define NZERO		20


/* Number of bits in a word of type `int'.  */
#ifdef INT_MAX
# if INT_MAX == 32767
#  define WORD_BIT	16
# else
#  if INT_MAX == 2147483647
#   define WORD_BIT	32
#  else
/* Safe assumption.  */
#   define WORD_BIT	64
#  endif
# endif
#elif defined __INT_MAX__
# if __INT_MAX__ == 32767
#  define WORD_BIT	16
# else
#  if __INT_MAX__ == 2147483647
#   define WORD_BIT	32
#  else
/* Safe assumption.  */
#   define WORD_BIT	64
#  endif
# endif
#else
# define WORD_BIT	32
#endif

/* Number of bits in a word of type `long int'.  */
#ifdef LONG_MAX
# if LONG_MAX == 2147483647
#  define LONG_BIT	32
# else
/* Safe assumption.  */
#  define LONG_BIT	64
# endif
#elif defined __LONG_MAX__
# if __LONG_MAX__ == 2147483647
#  define LONG_BIT	32
# else
/* Safe assumption.  */
#  define LONG_BIT	64
# endif
#else
# include <bits/wordsize.h>
# if __WORDSIZE == 64
#  define LONG_BIT	64
# else
#  define LONG_BIT	32
# endif
#endif

#endif /* bits/xopen_lim.h */
