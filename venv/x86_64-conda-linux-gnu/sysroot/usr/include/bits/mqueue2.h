/* Checking macros for mq functions.
   Copyright (C) 2007 Free Software Foundation, Inc.
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

#ifndef	_FCNTL_H
# error "Never include <bits/mqueue2.h> directly; use <mqueue.h> instead."
#endif

/* Check that calls to mq_open with O_CREAT set have an appropriate third and fourth
   parameter.  */
extern mqd_t mq_open (__const char *__name, int __oflag, ...)
     __THROW __nonnull ((1));
extern mqd_t __mq_open_2 (__const char *__name, int __oflag) __nonnull ((1));
extern mqd_t __REDIRECT (__mq_open_alias, (__const char *__name, int __oflag, ...),
			 mq_open) __nonnull ((1));
__errordecl (__mq_open_wrong_number_of_args,
	     "mq_open can be called either with 2 or 4 arguments");
__errordecl (__mq_open_missing_mode_and_attr,
	     "mq_open with O_CREAT in second argument needs 4 arguments");

__extern_always_inline mqd_t
mq_open (__const char *__name, int __oflag, ...)
{
  if (__va_arg_pack_len () != 0 && __va_arg_pack_len () != 2)
    __mq_open_wrong_number_of_args ();

  if (__builtin_constant_p (__oflag))
    {
      if ((__oflag & O_CREAT) != 0 && __va_arg_pack_len () == 0)
	{
	  __mq_open_missing_mode_and_attr ();
	  return __mq_open_2 (__name, __oflag);
	}
      return __mq_open_alias (__name, __oflag, __va_arg_pack ());
    }

  if (__va_arg_pack_len () == 0)
    return __mq_open_2 (__name, __oflag);

  return __mq_open_alias (__name, __oflag, __va_arg_pack ());
}
