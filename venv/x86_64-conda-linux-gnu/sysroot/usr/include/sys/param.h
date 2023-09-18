/* Copyright (C) 1995-1997,2000,2001,2003,2008 Free Software Foundation, Inc.
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

#ifndef _SYS_PARAM_H
#define _SYS_PARAM_H	1

#ifndef ARG_MAX
# define __undef_ARG_MAX
#endif

#include <limits.h>
#include <linux/limits.h>
#include <linux/param.h>

/* The kernel headers defines ARG_MAX.  The value is wrong, though.  */
#ifndef __undef_ARG_MAX
# undef ARG_MAX
# undef __undef_ARG_MAX
#endif

/* BSD names for some <limits.h> values.  */

#define	NBBY		CHAR_BIT
#ifndef	NGROUPS
# define NGROUPS	NGROUPS_MAX
#endif
#define	MAXSYMLINKS	20
#define	CANBSIZ		MAX_CANON
#define MAXPATHLEN	PATH_MAX
/* The following are not really correct but it is a value we used for a
   long time and which seems to be usable.  People should not use NOFILE
   and NCARGS anyway.  */
#define NOFILE		256
#define	NCARGS		131072


#include <sys/types.h>

/* Bit map related macros.  */
#define	setbit(a,i)	((a)[(i)/NBBY] |= 1<<((i)%NBBY))
#define	clrbit(a,i)	((a)[(i)/NBBY] &= ~(1<<((i)%NBBY)))
#define	isset(a,i)	((a)[(i)/NBBY] & (1<<((i)%NBBY)))
#define	isclr(a,i)	(((a)[(i)/NBBY] & (1<<((i)%NBBY))) == 0)

/* Macros for counting and rounding.  */
#ifndef howmany
# define howmany(x, y)	(((x) + ((y) - 1)) / (y))
#endif
#ifdef __GNUC__
# define roundup(x, y)	(__builtin_constant_p (y) && powerof2 (y)	      \
			 ? (((x) + (y) - 1) & ~((y) - 1))		      \
			 : ((((x) + ((y) - 1)) / (y)) * (y)))
#else
# define roundup(x, y)	((((x) + ((y) - 1)) / (y)) * (y))
#endif
#define powerof2(x)	((((x) - 1) & (x)) == 0)

/* Macros for min/max.  */
#define	MIN(a,b) (((a)<(b))?(a):(b))
#define	MAX(a,b) (((a)>(b))?(a):(b))


/* Unit of `st_blocks'.  */
#define DEV_BSIZE       512


#endif	/* sys/param.h */
