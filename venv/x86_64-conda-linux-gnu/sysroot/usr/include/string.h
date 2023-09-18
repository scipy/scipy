/* Copyright (C) 1991-1993,1995-2004,2007,2009,2010
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

/*
 *	ISO C99 Standard: 7.21 String handling	<string.h>
 */

#ifndef	_STRING_H
#define	_STRING_H	1

#include <features.h>

__BEGIN_DECLS

/* Get size_t and NULL from <stddef.h>.  */
#define	__need_size_t
#define	__need_NULL
#include <stddef.h>

/* Tell the caller that we provide correct C++ prototypes.  */
#if defined __cplusplus && __GNUC_PREREQ (4, 4)
# define __CORRECT_ISO_CPP_STRING_H_PROTO
#endif


__BEGIN_NAMESPACE_STD
/* Copy N bytes of SRC to DEST.  */
extern void *memcpy (void *__restrict __dest,
		     __const void *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));
/* Copy N bytes of SRC to DEST, guaranteeing
   correct behavior for overlapping strings.  */
extern void *memmove (void *__dest, __const void *__src, size_t __n)
     __THROW __nonnull ((1, 2));
__END_NAMESPACE_STD

/* Copy no more than N bytes of SRC to DEST, stopping when C is found.
   Return the position in DEST one byte past where C was copied,
   or NULL if C was not found in the first N bytes of SRC.  */
#if defined __USE_SVID || defined __USE_BSD || defined __USE_XOPEN
extern void *memccpy (void *__restrict __dest, __const void *__restrict __src,
		      int __c, size_t __n)
     __THROW __nonnull ((1, 2));
#endif /* SVID.  */


__BEGIN_NAMESPACE_STD
/* Set N bytes of S to C.  */
extern void *memset (void *__s, int __c, size_t __n) __THROW __nonnull ((1));

/* Compare N bytes of S1 and S2.  */
extern int memcmp (__const void *__s1, __const void *__s2, size_t __n)
     __THROW __attribute_pure__ __nonnull ((1, 2));

/* Search N bytes of S for C.  */
#ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern void *memchr (void *__s, int __c, size_t __n)
      __THROW __asm ("memchr") __attribute_pure__ __nonnull ((1));
extern __const void *memchr (__const void *__s, int __c, size_t __n)
      __THROW __asm ("memchr") __attribute_pure__ __nonnull ((1));

# ifdef __OPTIMIZE__
__extern_always_inline void *
memchr (void *__s, int __c, size_t __n) __THROW
{
  return __builtin_memchr (__s, __c, __n);
}

__extern_always_inline __const void *
memchr (__const void *__s, int __c, size_t __n) __THROW
{
  return __builtin_memchr (__s, __c, __n);
}
# endif
}
#else
extern void *memchr (__const void *__s, int __c, size_t __n)
      __THROW __attribute_pure__ __nonnull ((1));
#endif
__END_NAMESPACE_STD

#ifdef __USE_GNU
/* Search in S for C.  This is similar to `memchr' but there is no
   length limit.  */
# ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++" void *rawmemchr (void *__s, int __c)
     __THROW __asm ("rawmemchr") __attribute_pure__ __nonnull ((1));
extern "C++" __const void *rawmemchr (__const void *__s, int __c)
     __THROW __asm ("rawmemchr") __attribute_pure__ __nonnull ((1));
# else
extern void *rawmemchr (__const void *__s, int __c)
     __THROW __attribute_pure__ __nonnull ((1));
# endif

/* Search N bytes of S for the final occurrence of C.  */
# ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++" void *memrchr (void *__s, int __c, size_t __n)
      __THROW __asm ("memrchr") __attribute_pure__ __nonnull ((1));
extern "C++" __const void *memrchr (__const void *__s, int __c, size_t __n)
      __THROW __asm ("memrchr") __attribute_pure__ __nonnull ((1));
# else
extern void *memrchr (__const void *__s, int __c, size_t __n)
      __THROW __attribute_pure__ __nonnull ((1));
# endif
#endif


__BEGIN_NAMESPACE_STD
/* Copy SRC to DEST.  */
extern char *strcpy (char *__restrict __dest, __const char *__restrict __src)
     __THROW __nonnull ((1, 2));
/* Copy no more than N characters of SRC to DEST.  */
extern char *strncpy (char *__restrict __dest,
		      __const char *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));

/* Append SRC onto DEST.  */
extern char *strcat (char *__restrict __dest, __const char *__restrict __src)
     __THROW __nonnull ((1, 2));
/* Append no more than N characters from SRC onto DEST.  */
extern char *strncat (char *__restrict __dest, __const char *__restrict __src,
		      size_t __n) __THROW __nonnull ((1, 2));

/* Compare S1 and S2.  */
extern int strcmp (__const char *__s1, __const char *__s2)
     __THROW __attribute_pure__ __nonnull ((1, 2));
/* Compare N characters of S1 and S2.  */
extern int strncmp (__const char *__s1, __const char *__s2, size_t __n)
     __THROW __attribute_pure__ __nonnull ((1, 2));

/* Compare the collated forms of S1 and S2.  */
extern int strcoll (__const char *__s1, __const char *__s2)
     __THROW __attribute_pure__ __nonnull ((1, 2));
/* Put a transformation of SRC into no more than N bytes of DEST.  */
extern size_t strxfrm (char *__restrict __dest,
		       __const char *__restrict __src, size_t __n)
     __THROW __nonnull ((2));
__END_NAMESPACE_STD

#ifdef __USE_XOPEN2K8
/* The following functions are equivalent to the both above but they
   take the locale they use for the collation as an extra argument.
   This is not standardsized but something like will come.  */
# include <xlocale.h>

/* Compare the collated forms of S1 and S2 using rules from L.  */
extern int strcoll_l (__const char *__s1, __const char *__s2, __locale_t __l)
     __THROW __attribute_pure__ __nonnull ((1, 2, 3));
/* Put a transformation of SRC into no more than N bytes of DEST.  */
extern size_t strxfrm_l (char *__dest, __const char *__src, size_t __n,
			 __locale_t __l) __THROW __nonnull ((2, 4));
#endif

#if defined __USE_SVID || defined __USE_BSD || defined __USE_XOPEN_EXTENDED \
    || defined __USE_XOPEN2K8
/* Duplicate S, returning an identical malloc'd string.  */
extern char *strdup (__const char *__s)
     __THROW __attribute_malloc__ __nonnull ((1));
#endif

/* Return a malloc'd copy of at most N bytes of STRING.  The
   resultant string is terminated even if no null terminator
   appears before STRING[N].  */
#if defined __USE_XOPEN2K8
extern char *strndup (__const char *__string, size_t __n)
     __THROW __attribute_malloc__ __nonnull ((1));
#endif

#if defined __USE_GNU && defined __GNUC__
/* Duplicate S, returning an identical alloca'd string.  */
# define strdupa(s)							      \
  (__extension__							      \
    ({									      \
      __const char *__old = (s);					      \
      size_t __len = strlen (__old) + 1;				      \
      char *__new = (char *) __builtin_alloca (__len);			      \
      (char *) memcpy (__new, __old, __len);				      \
    }))

/* Return an alloca'd copy of at most N bytes of string.  */
# define strndupa(s, n)							      \
  (__extension__							      \
    ({									      \
      __const char *__old = (s);					      \
      size_t __len = strnlen (__old, (n));				      \
      char *__new = (char *) __builtin_alloca (__len + 1);		      \
      __new[__len] = '\0';						      \
      (char *) memcpy (__new, __old, __len);				      \
    }))
#endif

__BEGIN_NAMESPACE_STD
/* Find the first occurrence of C in S.  */
#ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern char *strchr (char *__s, int __c)
     __THROW __asm ("strchr") __attribute_pure__ __nonnull ((1));
extern __const char *strchr (__const char *__s, int __c)
     __THROW __asm ("strchr") __attribute_pure__ __nonnull ((1));

# ifdef __OPTIMIZE__
__extern_always_inline char *
strchr (char *__s, int __c) __THROW
{
  return __builtin_strchr (__s, __c);
}

__extern_always_inline __const char *
strchr (__const char *__s, int __c) __THROW
{
  return __builtin_strchr (__s, __c);
}
# endif
}
#else
extern char *strchr (__const char *__s, int __c)
     __THROW __attribute_pure__ __nonnull ((1));
#endif
/* Find the last occurrence of C in S.  */
#ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern char *strrchr (char *__s, int __c)
     __THROW __asm ("strrchr") __attribute_pure__ __nonnull ((1));
extern __const char *strrchr (__const char *__s, int __c)
     __THROW __asm ("strrchr") __attribute_pure__ __nonnull ((1));

# ifdef __OPTIMIZE__
__extern_always_inline char *
strrchr (char *__s, int __c) __THROW
{
  return __builtin_strrchr (__s, __c);
}

__extern_always_inline __const char *
strrchr (__const char *__s, int __c) __THROW
{
  return __builtin_strrchr (__s, __c);
}
# endif
}
#else
extern char *strrchr (__const char *__s, int __c)
     __THROW __attribute_pure__ __nonnull ((1));
#endif
__END_NAMESPACE_STD

#ifdef __USE_GNU
/* This function is similar to `strchr'.  But it returns a pointer to
   the closing NUL byte in case C is not found in S.  */
# ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++" char *strchrnul (char *__s, int __c)
     __THROW __asm ("strchrnul") __attribute_pure__ __nonnull ((1));
extern "C++" __const char *strchrnul (__const char *__s, int __c)
     __THROW __asm ("strchrnul") __attribute_pure__ __nonnull ((1));
# else
extern char *strchrnul (__const char *__s, int __c)
     __THROW __attribute_pure__ __nonnull ((1));
# endif
#endif

__BEGIN_NAMESPACE_STD
/* Return the length of the initial segment of S which
   consists entirely of characters not in REJECT.  */
extern size_t strcspn (__const char *__s, __const char *__reject)
     __THROW __attribute_pure__ __nonnull ((1, 2));
/* Return the length of the initial segment of S which
   consists entirely of characters in ACCEPT.  */
extern size_t strspn (__const char *__s, __const char *__accept)
     __THROW __attribute_pure__ __nonnull ((1, 2));
/* Find the first occurrence in S of any character in ACCEPT.  */
#ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern char *strpbrk (char *__s, __const char *__accept)
     __THROW __asm ("strpbrk") __attribute_pure__ __nonnull ((1, 2));
extern __const char *strpbrk (__const char *__s, __const char *__accept)
     __THROW __asm ("strpbrk") __attribute_pure__ __nonnull ((1, 2));

# ifdef __OPTIMIZE__
__extern_always_inline char *
strpbrk (char *__s, __const char *__accept) __THROW
{
  return __builtin_strpbrk (__s, __accept);
}

__extern_always_inline __const char *
strpbrk (__const char *__s, __const char *__accept) __THROW
{
  return __builtin_strpbrk (__s, __accept);
}
# endif
}
#else
extern char *strpbrk (__const char *__s, __const char *__accept)
     __THROW __attribute_pure__ __nonnull ((1, 2));
#endif
/* Find the first occurrence of NEEDLE in HAYSTACK.  */
#ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern char *strstr (char *__haystack, __const char *__needle)
     __THROW __asm ("strstr") __attribute_pure__ __nonnull ((1, 2));
extern __const char *strstr (__const char *__haystack,
			     __const char *__needle)
     __THROW __asm ("strstr") __attribute_pure__ __nonnull ((1, 2));

# ifdef __OPTIMIZE__
__extern_always_inline char *
strstr (char *__haystack, __const char *__needle) __THROW
{
  return __builtin_strstr (__haystack, __needle);
}

__extern_always_inline __const char *
strstr (__const char *__haystack, __const char *__needle) __THROW
{
  return __builtin_strstr (__haystack, __needle);
}
# endif
}
#else
extern char *strstr (__const char *__haystack, __const char *__needle)
     __THROW __attribute_pure__ __nonnull ((1, 2));
#endif


/* Divide S into tokens separated by characters in DELIM.  */
extern char *strtok (char *__restrict __s, __const char *__restrict __delim)
     __THROW __nonnull ((2));
__END_NAMESPACE_STD

/* Divide S into tokens separated by characters in DELIM.  Information
   passed between calls are stored in SAVE_PTR.  */
extern char *__strtok_r (char *__restrict __s,
			 __const char *__restrict __delim,
			 char **__restrict __save_ptr)
     __THROW __nonnull ((2, 3));
#if defined __USE_POSIX || defined __USE_MISC
extern char *strtok_r (char *__restrict __s, __const char *__restrict __delim,
		       char **__restrict __save_ptr)
     __THROW __nonnull ((2, 3));
#endif

#ifdef __USE_GNU
/* Similar to `strstr' but this function ignores the case of both strings.  */
# ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++" char *strcasestr (char *__haystack, __const char *__needle)
     __THROW __asm ("strcasestr") __attribute_pure__ __nonnull ((1, 2));
extern "C++" __const char *strcasestr (__const char *__haystack,
				       __const char *__needle)
     __THROW __asm ("strcasestr") __attribute_pure__ __nonnull ((1, 2));
# else
extern char *strcasestr (__const char *__haystack, __const char *__needle)
     __THROW __attribute_pure__ __nonnull ((1, 2));
# endif
#endif

#ifdef __USE_GNU
/* Find the first occurrence of NEEDLE in HAYSTACK.
   NEEDLE is NEEDLELEN bytes long;
   HAYSTACK is HAYSTACKLEN bytes long.  */
extern void *memmem (__const void *__haystack, size_t __haystacklen,
		     __const void *__needle, size_t __needlelen)
     __THROW __attribute_pure__ __nonnull ((1, 3));

/* Copy N bytes of SRC to DEST, return pointer to bytes after the
   last written byte.  */
extern void *__mempcpy (void *__restrict __dest,
			__const void *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));
extern void *mempcpy (void *__restrict __dest,
		      __const void *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));
#endif


__BEGIN_NAMESPACE_STD
/* Return the length of S.  */
extern size_t strlen (__const char *__s)
     __THROW __attribute_pure__ __nonnull ((1));
__END_NAMESPACE_STD

#ifdef	__USE_XOPEN2K8
/* Find the length of STRING, but scan at most MAXLEN characters.
   If no '\0' terminator is found in that many characters, return MAXLEN.  */
extern size_t strnlen (__const char *__string, size_t __maxlen)
     __THROW __attribute_pure__ __nonnull ((1));
#endif


__BEGIN_NAMESPACE_STD
/* Return a string describing the meaning of the `errno' code in ERRNUM.  */
extern char *strerror (int __errnum) __THROW;
__END_NAMESPACE_STD
#if defined __USE_XOPEN2K || defined __USE_MISC
/* Reentrant version of `strerror'.
   There are 2 flavors of `strerror_r', GNU which returns the string
   and may or may not use the supplied temporary buffer and POSIX one
   which fills the string into the buffer.
   To use the POSIX version, -D_XOPEN_SOURCE=600 or -D_POSIX_C_SOURCE=200112L
   without -D_GNU_SOURCE is needed, otherwise the GNU version is
   preferred.  */
# if defined __USE_XOPEN2K && !defined __USE_GNU
/* Fill BUF with a string describing the meaning of the `errno' code in
   ERRNUM.  */
#  ifdef __REDIRECT_NTH
extern int __REDIRECT_NTH (strerror_r,
			   (int __errnum, char *__buf, size_t __buflen),
			   __xpg_strerror_r) __nonnull ((2));
#  else
extern int __xpg_strerror_r (int __errnum, char *__buf, size_t __buflen)
     __THROW __nonnull ((2));
#   define strerror_r __xpg_strerror_r
#  endif
# else
/* If a temporary buffer is required, at most BUFLEN bytes of BUF will be
   used.  */
extern char *strerror_r (int __errnum, char *__buf, size_t __buflen)
     __THROW __nonnull ((2));
# endif
#endif

#ifdef __USE_XOPEN2K8
/* Translate error number to string according to the locale L.  */
extern char *strerror_l (int __errnum, __locale_t __l) __THROW;
#endif


/* We define this function always since `bzero' is sometimes needed when
   the namespace rules does not allow this.  */
extern void __bzero (void *__s, size_t __n) __THROW __nonnull ((1));

#ifdef __USE_BSD
/* Copy N bytes of SRC to DEST (like memmove, but args reversed).  */
extern void bcopy (__const void *__src, void *__dest, size_t __n)
     __THROW __nonnull ((1, 2));

/* Set N bytes of S to 0.  */
extern void bzero (void *__s, size_t __n) __THROW __nonnull ((1));

/* Compare N bytes of S1 and S2 (same as memcmp).  */
extern int bcmp (__const void *__s1, __const void *__s2, size_t __n)
     __THROW __attribute_pure__ __nonnull ((1, 2));

/* Find the first occurrence of C in S (same as strchr).  */
# ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern char *index (char *__s, int __c)
     __THROW __asm ("index") __attribute_pure__ __nonnull ((1));
extern __const char *index (__const char *__s, int __c)
     __THROW __asm ("index") __attribute_pure__ __nonnull ((1));

#  if defined __OPTIMIZE__ && !defined __CORRECT_ISO_CPP_STRINGS_H_PROTO
__extern_always_inline char *
index (char *__s, int __c) __THROW
{
  return __builtin_index (__s, __c);
}

__extern_always_inline __const char *
index (__const char *__s, int __c) __THROW
{
  return __builtin_index (__s, __c);
}
#  endif
}
# else
extern char *index (__const char *__s, int __c)
     __THROW __attribute_pure__ __nonnull ((1));
# endif

/* Find the last occurrence of C in S (same as strrchr).  */
# ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++"
{
extern char *rindex (char *__s, int __c)
     __THROW __asm ("rindex") __attribute_pure__ __nonnull ((1));
extern __const char *rindex (__const char *__s, int __c)
     __THROW __asm ("rindex") __attribute_pure__ __nonnull ((1));

#  if defined __OPTIMIZE__ && !defined __CORRECT_ISO_CPP_STRINGS_H_PROTO
__extern_always_inline char *
rindex (char *__s, int __c) __THROW
{
  return __builtin_rindex (__s, __c);
}

__extern_always_inline __const char *
rindex (__const char *__s, int __c) __THROW
{
  return __builtin_rindex (__s, __c);
}
#endif
}
# else
extern char *rindex (__const char *__s, int __c)
     __THROW __attribute_pure__ __nonnull ((1));
# endif

/* Return the position of the first bit set in I, or 0 if none are set.
   The least-significant bit is position 1, the most-significant 32.  */
extern int ffs (int __i) __THROW __attribute__ ((__const__));

/* The following two functions are non-standard but necessary for non-32 bit
   platforms.  */
# ifdef	__USE_GNU
extern int ffsl (long int __l) __THROW __attribute__ ((__const__));
#  ifdef __GNUC__
__extension__ extern int ffsll (long long int __ll)
     __THROW __attribute__ ((__const__));
#  endif
# endif

/* Compare S1 and S2, ignoring case.  */
extern int strcasecmp (__const char *__s1, __const char *__s2)
     __THROW __attribute_pure__ __nonnull ((1, 2));

/* Compare no more than N chars of S1 and S2, ignoring case.  */
extern int strncasecmp (__const char *__s1, __const char *__s2, size_t __n)
     __THROW __attribute_pure__ __nonnull ((1, 2));
#endif /* Use BSD.  */

#ifdef	__USE_GNU
/* Again versions of a few functions which use the given locale instead
   of the global one.  */
extern int strcasecmp_l (__const char *__s1, __const char *__s2,
			 __locale_t __loc)
     __THROW __attribute_pure__ __nonnull ((1, 2, 3));

extern int strncasecmp_l (__const char *__s1, __const char *__s2,
			  size_t __n, __locale_t __loc)
     __THROW __attribute_pure__ __nonnull ((1, 2, 4));
#endif

#ifdef	__USE_BSD
/* Return the next DELIM-delimited token from *STRINGP,
   terminating it with a '\0', and update *STRINGP to point past it.  */
extern char *strsep (char **__restrict __stringp,
		     __const char *__restrict __delim)
     __THROW __nonnull ((1, 2));
#endif

#ifdef	__USE_XOPEN2K8
/* Return a string describing the meaning of the signal number in SIG.  */
extern char *strsignal (int __sig) __THROW;

/* Copy SRC to DEST, returning the address of the terminating '\0' in DEST.  */
extern char *__stpcpy (char *__restrict __dest, __const char *__restrict __src)
     __THROW __nonnull ((1, 2));
extern char *stpcpy (char *__restrict __dest, __const char *__restrict __src)
     __THROW __nonnull ((1, 2));

/* Copy no more than N characters of SRC to DEST, returning the address of
   the last character written into DEST.  */
extern char *__stpncpy (char *__restrict __dest,
			__const char *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));
extern char *stpncpy (char *__restrict __dest,
		      __const char *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));
#endif

#ifdef	__USE_GNU
/* Compare S1 and S2 as strings holding name & indices/version numbers.  */
extern int strverscmp (__const char *__s1, __const char *__s2)
     __THROW __attribute_pure__ __nonnull ((1, 2));

/* Sautee STRING briskly.  */
extern char *strfry (char *__string) __THROW __nonnull ((1));

/* Frobnicate N bytes of S.  */
extern void *memfrob (void *__s, size_t __n) __THROW __nonnull ((1));

# ifndef basename
/* Return the file name within directory of FILENAME.  We don't
   declare the function if the `basename' macro is available (defined
   in <libgen.h>) which makes the XPG version of this function
   available.  */
#  ifdef __CORRECT_ISO_CPP_STRING_H_PROTO
extern "C++" char *basename (char *__filename)
     __THROW __asm ("basename") __nonnull ((1));
extern "C++" __const char *basename (__const char *__filename)
     __THROW __asm ("basename") __nonnull ((1));
#  else
extern char *basename (__const char *__filename) __THROW __nonnull ((1));
#  endif
# endif
#endif


#if defined __GNUC__ && __GNUC__ >= 2
# if defined __OPTIMIZE__ && !defined __OPTIMIZE_SIZE__ \
     && !defined __NO_INLINE__ && !defined __cplusplus
/* When using GNU CC we provide some optimized versions of selected
   functions from this header.  There are two kinds of optimizations:

   - machine-dependent optimizations, most probably using inline
     assembler code; these might be quite expensive since the code
     size can increase significantly.
     These optimizations are not used unless the symbol
	__USE_STRING_INLINES
     is defined before including this header.

   - machine-independent optimizations which do not increase the
     code size significantly and which optimize mainly situations
     where one or more arguments are compile-time constants.
     These optimizations are used always when the compiler is
     taught to optimize.

   One can inhibit all optimizations by defining __NO_STRING_INLINES.  */

/* Get the machine-dependent optimizations (if any).  */
#  include <bits/string.h>

/* These are generic optimizations which do not add too much inline code.  */
#  include <bits/string2.h>
# endif

# if __USE_FORTIFY_LEVEL > 0 && defined __extern_always_inline
/* Functions with security checks.  */
#  include <bits/string3.h>
# endif
#endif

__END_DECLS

#endif /* string.h  */
