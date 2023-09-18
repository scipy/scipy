/* Copyright (C) 2004, 2005, 2007, 2009 Free Software Foundation, Inc.
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

#ifndef _STRING_H
# error "Never use <bits/string3.h> directly; include <string.h> instead."
#endif

__warndecl (__warn_memset_zero_len,
	    "memset used with constant zero length parameter; this could be due to transposed parameters");

#ifndef __cplusplus
/* XXX This is temporarily.  We should not redefine any of the symbols
   and instead integrate the error checking into the original
   definitions.  */
# undef memcpy
# undef memmove
# undef memset
# undef strcat
# undef strcpy
# undef strncat
# undef strncpy
# ifdef __USE_GNU
#  undef mempcpy
#  undef stpcpy
# endif
# ifdef __USE_BSD
#  undef bcopy
#  undef bzero
# endif
#endif


__extern_always_inline void *
__NTH (memcpy (void *__restrict __dest, __const void *__restrict __src,
	       size_t __len))
{
  return __builtin___memcpy_chk (__dest, __src, __len, __bos0 (__dest));
}

__extern_always_inline void *
__NTH (memmove (void *__dest, __const void *__src, size_t __len))
{
  return __builtin___memmove_chk (__dest, __src, __len, __bos0 (__dest));
}

#ifdef __USE_GNU
__extern_always_inline void *
__NTH (mempcpy (void *__restrict __dest, __const void *__restrict __src,
		size_t __len))
{
  return __builtin___mempcpy_chk (__dest, __src, __len, __bos0 (__dest));
}
#endif


/* The first two tests here help to catch a somewhat common problem
   where the second and third parameter are transposed.  This is
   especially problematic if the intended fill value is zero.  In this
   case no work is done at all.  We detect these problems by referring
   non-existing functions.  */
__extern_always_inline void *
__NTH (memset (void *__dest, int __ch, size_t __len))
{
  if (__builtin_constant_p (__len) && __len == 0
      && (!__builtin_constant_p (__ch) || __ch != 0))
    {
      __warn_memset_zero_len ();
      return __dest;
    }
  return __builtin___memset_chk (__dest, __ch, __len, __bos0 (__dest));
}

#ifdef __USE_BSD
__extern_always_inline void
__NTH (bcopy (__const void *__src, void *__dest, size_t __len))
{
  (void) __builtin___memmove_chk (__dest, __src, __len, __bos0 (__dest));
}

__extern_always_inline void
__NTH (bzero (void *__dest, size_t __len))
{
  (void) __builtin___memset_chk (__dest, '\0', __len, __bos0 (__dest));
}
#endif

__extern_always_inline char *
__NTH (strcpy (char *__restrict __dest, __const char *__restrict __src))
{
  return __builtin___strcpy_chk (__dest, __src, __bos (__dest));
}

#ifdef __USE_GNU
__extern_always_inline char *
__NTH (stpcpy (char *__restrict __dest, __const char *__restrict __src))
{
  return __builtin___stpcpy_chk (__dest, __src, __bos (__dest));
}
#endif


__extern_always_inline char *
__NTH (strncpy (char *__restrict __dest, __const char *__restrict __src,
		size_t __len))
{
  return __builtin___strncpy_chk (__dest, __src, __len, __bos (__dest));
}

// XXX We have no corresponding builtin yet.
extern char *__stpncpy_chk (char *__dest, __const char *__src, size_t __n,
			    size_t __destlen) __THROW;
extern char *__REDIRECT_NTH (__stpncpy_alias, (char *__dest,
					       __const char *__src,
					       size_t __n), stpncpy);

__extern_always_inline char *
__NTH (stpncpy (char *__dest, __const char *__src, size_t __n))
{
  if (__bos (__dest) != (size_t) -1
      && (!__builtin_constant_p (__n) || __n <= __bos (__dest)))
    return __stpncpy_chk (__dest, __src, __n, __bos (__dest));
  return __stpncpy_alias (__dest, __src, __n);
}


__extern_always_inline char *
__NTH (strcat (char *__restrict __dest, __const char *__restrict __src))
{
  return __builtin___strcat_chk (__dest, __src, __bos (__dest));
}


__extern_always_inline char *
__NTH (strncat (char *__restrict __dest, __const char *__restrict __src,
		size_t __len))
{
  return __builtin___strncat_chk (__dest, __src, __len, __bos (__dest));
}
