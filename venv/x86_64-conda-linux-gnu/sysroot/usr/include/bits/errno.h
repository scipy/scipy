/* Error constants.  Linux specific version.
   Copyright (C) 1996-1999, 2005, 2009 Free Software Foundation, Inc.
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

#ifdef _ERRNO_H

# undef EDOM
# undef EILSEQ
# undef ERANGE
# include <linux/errno.h>

/* Linux has no ENOTSUP error code.  */
# define ENOTSUP EOPNOTSUPP

/* Older Linux versions also had no ECANCELED error code.  */
# ifndef ECANCELED
#  define ECANCELED	125
# endif

/* Support for error codes to support robust mutexes was added later, too.  */
# ifndef EOWNERDEAD
#  define EOWNERDEAD		130
#  define ENOTRECOVERABLE	131
# endif

# ifndef ERFKILL
#  define ERFKILL		132
# endif

# ifndef __ASSEMBLER__
/* Function to get address of global `errno' variable.  */
extern int *__errno_location (void) __THROW __attribute__ ((__const__));

#  if !defined _LIBC || defined _LIBC_REENTRANT
/* When using threads, errno is a per-thread value.  */
#   define errno (*__errno_location ())
#  endif
# endif /* !__ASSEMBLER__ */
#endif /* _ERRNO_H */

#if !defined _ERRNO_H && defined __need_Emath
/* This is ugly but the kernel header is not clean enough.  We must
   define only the values EDOM, EILSEQ and ERANGE in case __need_Emath is
   defined.  */
# define EDOM	33	/* Math argument out of domain of function.  */
# define EILSEQ	84	/* Illegal byte sequence.  */
# define ERANGE	34	/* Math result not representable.  */
#endif /* !_ERRNO_H && __need_Emath */
