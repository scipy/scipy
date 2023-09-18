/* `HUGE_VAL' constant for IEEE 754 machines (where it is infinity).
   Used by <stdlib.h> and <math.h> functions for overflow.
   Copyright (C) 1992, 1995, 1996, 1997, 1999, 2000, 2004
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

#ifndef _MATH_H
# error "Never use <bits/huge_val.h> directly; include <math.h> instead."
#endif

/* IEEE positive infinity (-HUGE_VAL is negative infinity).  */

#if __GNUC_PREREQ(3,3)
# define HUGE_VAL	(__builtin_huge_val())
#elif __GNUC_PREREQ(2,96)
# define HUGE_VAL	(__extension__ 0x1.0p2047)
#elif defined __GNUC__

# define HUGE_VAL \
  (__extension__							      \
   ((union { unsigned __l __attribute__((__mode__(__DI__))); double __d; })   \
    { __l: 0x7ff0000000000000ULL }).__d)

#else /* not GCC */

# include <endian.h>

typedef union { unsigned char __c[8]; double __d; } __huge_val_t;

# if __BYTE_ORDER == __BIG_ENDIAN
#  define __HUGE_VAL_bytes	{ 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
# endif
# if __BYTE_ORDER == __LITTLE_ENDIAN
#  define __HUGE_VAL_bytes	{ 0, 0, 0, 0, 0, 0, 0xf0, 0x7f }
# endif

static __huge_val_t __huge_val = { __HUGE_VAL_bytes };
# define HUGE_VAL	(__huge_val.__d)

#endif	/* GCC.  */
