/* Default `HUGE_VALL' constant.
   Used by <stdlib.h> and <math.h> functions for overflow.
   Copyright (C) 1992, 1996, 1997, 2004 Free Software Foundation, Inc.
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
# error "Never use <bits/huge_vall.h> directly; include <math.h> instead."
#endif

#if __GNUC_PREREQ(3,3)
# define HUGE_VALL	(__builtin_huge_vall())
#else
# define HUGE_VALL	((long double) HUGE_VAL)
#endif
