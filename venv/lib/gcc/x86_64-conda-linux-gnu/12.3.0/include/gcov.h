/* GCOV interface routines.
   Copyright (C) 2017-2022 Free Software Foundation, Inc.

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 3, or (at your option) any later
   version.

   GCC is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   for more details.

   Under Section 7 of GPL version 3, you are granted additional
   permissions described in the GCC Runtime Library Exception, version
   3.1, as published by the Free Software Foundation.

   You should have received a copy of the GNU General Public License and
   a copy of the GCC Runtime Library Exception along with this program;
   see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
   <http://www.gnu.org/licenses/>.  */

#ifndef GCC_GCOV_H
#define GCC_GCOV_H

struct gcov_info;

/* Set all counters to zero.  */

extern void __gcov_reset (void);

/* Write profile information to a file.  */

extern void __gcov_dump (void);

/* Convert the gcov information referenced by INFO to a gcda data stream.
   The FILENAME_FN callback is called exactly once with the filename associated
   with the gcov information.  The filename may be NULL.  Afterwards, the
   DUMP_FN callback is subsequently called with chunks (the begin and length of
   the chunk are passed as the first two callback parameters) of the gcda data
   stream.  The ALLOCATE_FN callback shall allocate memory with a size in
   characters specified by the first callback parameter.  The ARG parameter is
   a user-provided argument passed as the last argument to the callback
   functions.  */

extern void
__gcov_info_to_gcda (const struct gcov_info *__info,
		     void (*__filename_fn) (const char *, void *),
		     void (*__dump_fn) (const void *, unsigned, void *),
		     void *(*__allocate_fn) (unsigned, void *),
		     void *__arg);

#endif /* GCC_GCOV_H */
