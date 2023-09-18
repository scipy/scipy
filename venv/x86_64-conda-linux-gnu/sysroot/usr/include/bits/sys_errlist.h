/* Declare sys_errlist and sys_nerr, or don't.  Compatibility (do) version.
   Copyright (C) 2002 Free Software Foundation, Inc.
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

#ifndef _STDIO_H
# error "Never include <bits/sys_errlist.h> directly; use <stdio.h> instead."
#endif

/* sys_errlist and sys_nerr are deprecated.  Use strerror instead.  */

#ifdef  __USE_BSD
extern int sys_nerr;
extern __const char *__const sys_errlist[];
#endif
#ifdef  __USE_GNU
extern int _sys_nerr;
extern __const char *__const _sys_errlist[];
#endif
