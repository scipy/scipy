/* -mlong-double-64 compatibility mode for <stdlib.h> functions.
   Copyright (C) 2006, 2007 Free Software Foundation, Inc.
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

#ifndef _STDLIB_H
# error "Never include <bits/stdlib-ldbl.h> directly; use <stdlib.h> instead."
#endif

#ifdef	__USE_ISOC99
__BEGIN_NAMESPACE_C99
__LDBL_REDIR1_DECL (strtold, strtod)
__END_NAMESPACE_C99
#endif

#ifdef __USE_GNU
__LDBL_REDIR1_DECL (strtold_l, strtod_l)
#endif

#if defined __USE_SVID || defined __USE_XOPEN_EXTENDED
# ifdef __USE_MISC
__LDBL_REDIR1_DECL (qecvt, ecvt)
__LDBL_REDIR1_DECL (qfcvt, fcvt)
__LDBL_REDIR1_DECL (qgcvt, gcvt)
__LDBL_REDIR1_DECL (qecvt_r, ecvt_r)
__LDBL_REDIR1_DECL (qfcvt_r, fcvt_r)
# endif
#endif
