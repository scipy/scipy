/* Copyright (C) 1995-2008, 2009, 2010 Free Software Foundation, Inc.
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
 *      ISO C99 Standard: 7.24
 *	Extended multibyte and wide character utilities	<wchar.h>
 */

#ifndef _WCHAR_H

#if !defined __need_mbstate_t && !defined __need_wint_t
# define _WCHAR_H 1
# include <features.h>
#endif

#ifdef _WCHAR_H
/* Get FILE definition.  */
# define __need___FILE
# if defined __USE_UNIX98 || defined __USE_XOPEN2K
#  define __need_FILE
# endif
# include <stdio.h>
/* Get va_list definition.  */
# define __need___va_list
# include <stdarg.h>

# include <bits/wchar.h>

/* Get size_t, wchar_t, wint_t and NULL from <stddef.h>.  */
# define __need_size_t
# define __need_wchar_t
# define __need_NULL
#endif
#if defined _WCHAR_H || defined __need_wint_t || !defined __WINT_TYPE__
# undef __need_wint_t
# define __need_wint_t
# include <stddef.h>

/* We try to get wint_t from <stddef.h>, but not all GCC versions define it
   there.  So define it ourselves if it remains undefined.  */
# ifndef _WINT_T
/* Integral type unchanged by default argument promotions that can
   hold any value corresponding to members of the extended character
   set, as well as at least one value that does not correspond to any
   member of the extended character set.  */
#  define _WINT_T
typedef unsigned int wint_t;
# else
/* Work around problems with the <stddef.h> file which doesn't put
   wint_t in the std namespace.  */
#  if defined __cplusplus && defined _GLIBCPP_USE_NAMESPACES \
      && defined __WINT_TYPE__
__BEGIN_NAMESPACE_STD
typedef __WINT_TYPE__ wint_t;
__END_NAMESPACE_STD
#  endif
# endif

/* Tell the caller that we provide correct C++ prototypes.  */
# if defined __cplusplus && __GNUC_PREREQ (4, 4)
#  define __CORRECT_ISO_CPP_WCHAR_H_PROTO
# endif
#endif

#if (defined _WCHAR_H || defined __need_mbstate_t) && !defined __mbstate_t_defined
# define __mbstate_t_defined	1
/* Conversion state information.  */
typedef struct
{
  int __count;
  union
  {
# ifdef __WINT_TYPE__
    __WINT_TYPE__ __wch;
# else
    wint_t __wch;
# endif
    char __wchb[4];
  } __value;		/* Value so far.  */
} __mbstate_t;
#endif
#undef __need_mbstate_t


/* The rest of the file is only used if used if __need_mbstate_t is not
   defined.  */
#ifdef _WCHAR_H

__BEGIN_NAMESPACE_C99
/* Public type.  */
typedef __mbstate_t mbstate_t;
__END_NAMESPACE_C99
#ifdef __USE_GNU
__USING_NAMESPACE_C99(mbstate_t)
#endif

#ifndef WCHAR_MIN
/* These constants might also be defined in <inttypes.h>.  */
# define WCHAR_MIN __WCHAR_MIN
# define WCHAR_MAX __WCHAR_MAX
#endif

#ifndef WEOF
# define WEOF (0xffffffffu)
#endif

/* For XPG4 compliance we have to define the stuff from <wctype.h> here
   as well.  */
#if defined __USE_XOPEN && !defined __USE_UNIX98
# include <wctype.h>
#endif


__BEGIN_DECLS

__BEGIN_NAMESPACE_STD
/* This incomplete type is defined in <time.h> but needed here because
   of `wcsftime'.  */
struct tm;
__END_NAMESPACE_STD
/* XXX We have to clean this up at some point.  Since tm is in the std
   namespace but wcsftime is in __c99 the type wouldn't be found
   without inserting it in the global namespace.  */
__USING_NAMESPACE_STD(tm)


__BEGIN_NAMESPACE_STD
/* Copy SRC to DEST.  */
extern wchar_t *wcscpy (wchar_t *__restrict __dest,
			__const wchar_t *__restrict __src) __THROW;
/* Copy no more than N wide-characters of SRC to DEST.  */
extern wchar_t *wcsncpy (wchar_t *__restrict __dest,
			 __const wchar_t *__restrict __src, size_t __n)
     __THROW;

/* Append SRC onto DEST.  */
extern wchar_t *wcscat (wchar_t *__restrict __dest,
			__const wchar_t *__restrict __src) __THROW;
/* Append no more than N wide-characters of SRC onto DEST.  */
extern wchar_t *wcsncat (wchar_t *__restrict __dest,
			 __const wchar_t *__restrict __src, size_t __n)
     __THROW;

/* Compare S1 and S2.  */
extern int wcscmp (__const wchar_t *__s1, __const wchar_t *__s2)
     __THROW __attribute_pure__;
/* Compare N wide-characters of S1 and S2.  */
extern int wcsncmp (__const wchar_t *__s1, __const wchar_t *__s2, size_t __n)
     __THROW __attribute_pure__;
__END_NAMESPACE_STD

#ifdef __USE_XOPEN2K8
/* Compare S1 and S2, ignoring case.  */
extern int wcscasecmp (__const wchar_t *__s1, __const wchar_t *__s2) __THROW;

/* Compare no more than N chars of S1 and S2, ignoring case.  */
extern int wcsncasecmp (__const wchar_t *__s1, __const wchar_t *__s2,
			size_t __n) __THROW;

/* Similar to the two functions above but take the information from
   the provided locale and not the global locale.  */
# include <xlocale.h>

extern int wcscasecmp_l (__const wchar_t *__s1, __const wchar_t *__s2,
			 __locale_t __loc) __THROW;

extern int wcsncasecmp_l (__const wchar_t *__s1, __const wchar_t *__s2,
			  size_t __n, __locale_t __loc) __THROW;
#endif

__BEGIN_NAMESPACE_STD
/* Compare S1 and S2, both interpreted as appropriate to the
   LC_COLLATE category of the current locale.  */
extern int wcscoll (__const wchar_t *__s1, __const wchar_t *__s2) __THROW;
/* Transform S2 into array pointed to by S1 such that if wcscmp is
   applied to two transformed strings the result is the as applying
   `wcscoll' to the original strings.  */
extern size_t wcsxfrm (wchar_t *__restrict __s1,
		       __const wchar_t *__restrict __s2, size_t __n) __THROW;
__END_NAMESPACE_STD

#ifdef __USE_XOPEN2K8
/* Similar to the two functions above but take the information from
   the provided locale and not the global locale.  */

/* Compare S1 and S2, both interpreted as appropriate to the
   LC_COLLATE category of the given locale.  */
extern int wcscoll_l (__const wchar_t *__s1, __const wchar_t *__s2,
		      __locale_t __loc) __THROW;

/* Transform S2 into array pointed to by S1 such that if wcscmp is
   applied to two transformed strings the result is the as applying
   `wcscoll' to the original strings.  */
extern size_t wcsxfrm_l (wchar_t *__s1, __const wchar_t *__s2,
			 size_t __n, __locale_t __loc) __THROW;

/* Duplicate S, returning an identical malloc'd string.  */
extern wchar_t *wcsdup (__const wchar_t *__s) __THROW __attribute_malloc__;
#endif

__BEGIN_NAMESPACE_STD
/* Find the first occurrence of WC in WCS.  */
#ifdef __CORRECT_ISO_CPP_WCHAR_H_PROTO
extern "C++" wchar_t *wcschr (wchar_t *__wcs, wchar_t __wc)
     __THROW __asm ("wcschr") __attribute_pure__;
extern "C++" __const wchar_t *wcschr (__const wchar_t *__wcs, wchar_t __wc)
     __THROW __asm ("wcschr") __attribute_pure__;
#else
extern wchar_t *wcschr (__const wchar_t *__wcs, wchar_t __wc)
     __THROW __attribute_pure__;
#endif
/* Find the last occurrence of WC in WCS.  */
#ifdef __CORRECT_ISO_CPP_WCHAR_H_PROTO
extern "C++" wchar_t *wcsrchr (wchar_t *__wcs, wchar_t __wc)
     __THROW __asm ("wcsrchr") __attribute_pure__;
extern "C++" __const wchar_t *wcsrchr (__const wchar_t *__wcs, wchar_t __wc)
     __THROW __asm ("wcsrchr") __attribute_pure__;
#else
extern wchar_t *wcsrchr (__const wchar_t *__wcs, wchar_t __wc)
     __THROW __attribute_pure__;
#endif
__END_NAMESPACE_STD

#ifdef __USE_GNU
/* This function is similar to `wcschr'.  But it returns a pointer to
   the closing NUL wide character in case C is not found in S.  */
extern wchar_t *wcschrnul (__const wchar_t *__s, wchar_t __wc)
     __THROW __attribute_pure__;
#endif

__BEGIN_NAMESPACE_STD
/* Return the length of the initial segmet of WCS which
   consists entirely of wide characters not in REJECT.  */
extern size_t wcscspn (__const wchar_t *__wcs, __const wchar_t *__reject)
     __THROW __attribute_pure__;
/* Return the length of the initial segmet of WCS which
   consists entirely of wide characters in  ACCEPT.  */
extern size_t wcsspn (__const wchar_t *__wcs, __const wchar_t *__accept)
     __THROW __attribute_pure__;
/* Find the first occurrence in WCS of any character in ACCEPT.  */
#ifdef __CORRECT_ISO_CPP_WCHAR_H_PROTO
extern "C++" wchar_t *wcspbrk (wchar_t *__wcs, __const wchar_t *__accept)
     __THROW __asm ("wcspbrk") __attribute_pure__;
extern "C++" __const wchar_t *wcspbrk (__const wchar_t *__wcs,
				       __const wchar_t *__accept)
     __THROW __asm ("wcspbrk") __attribute_pure__;
#else
extern wchar_t *wcspbrk (__const wchar_t *__wcs, __const wchar_t *__accept)
     __THROW __attribute_pure__;
#endif
/* Find the first occurrence of NEEDLE in HAYSTACK.  */
#ifdef __CORRECT_ISO_CPP_WCHAR_H_PROTO
extern "C++" wchar_t *wcsstr (wchar_t *__haystack, __const wchar_t *__needle)
     __THROW __asm ("wcsstr") __attribute_pure__;
extern "C++" __const wchar_t *wcsstr (__const wchar_t *__haystack,
				      __const wchar_t *__needle)
     __THROW __asm ("wcsstr") __attribute_pure__;
#else
extern wchar_t *wcsstr (__const wchar_t *__haystack, __const wchar_t *__needle)
     __THROW __attribute_pure__;
#endif

/* Divide WCS into tokens separated by characters in DELIM.  */
extern wchar_t *wcstok (wchar_t *__restrict __s,
			__const wchar_t *__restrict __delim,
			wchar_t **__restrict __ptr) __THROW;

/* Return the number of wide characters in S.  */
extern size_t wcslen (__const wchar_t *__s) __THROW __attribute_pure__;
__END_NAMESPACE_STD

#ifdef __USE_XOPEN
/* Another name for `wcsstr' from XPG4.  */
# ifdef __CORRECT_ISO_CPP_WCHAR_H_PROTO
extern "C++" wchar_t *wcswcs (wchar_t *__haystack, __const wchar_t *__needle)
     __THROW __asm ("wcswcs") __attribute_pure__;
extern "C++" __const wchar_t *wcswcs (__const wchar_t *__haystack,
				      __const wchar_t *__needle)
     __THROW __asm ("wcswcs") __attribute_pure__;
# else
extern wchar_t *wcswcs (__const wchar_t *__haystack, __const wchar_t *__needle)
     __THROW __attribute_pure__;
# endif
#endif

#ifdef __USE_XOPEN2K8
/* Return the number of wide characters in S, but at most MAXLEN.  */
extern size_t wcsnlen (__const wchar_t *__s, size_t __maxlen)
     __THROW __attribute_pure__;
#endif


__BEGIN_NAMESPACE_STD
/* Search N wide characters of S for C.  */
#ifdef __CORRECT_ISO_CPP_WCHAR_H_PROTO
extern "C++" wchar_t *wmemchr (wchar_t *__s, wchar_t __c, size_t __n)
     __THROW __asm ("wmemchr") __attribute_pure__;
extern "C++" __const wchar_t *wmemchr (__const wchar_t *__s, wchar_t __c,
				       size_t __n)
     __THROW __asm ("wmemchr") __attribute_pure__;
#else
extern wchar_t *wmemchr (__const wchar_t *__s, wchar_t __c, size_t __n)
     __THROW __attribute_pure__;
#endif

/* Compare N wide characters of S1 and S2.  */
extern int wmemcmp (__const wchar_t *__restrict __s1,
		    __const wchar_t *__restrict __s2, size_t __n)
     __THROW __attribute_pure__;

/* Copy N wide characters of SRC to DEST.  */
extern wchar_t *wmemcpy (wchar_t *__restrict __s1,
			 __const wchar_t *__restrict __s2, size_t __n) __THROW;

/* Copy N wide characters of SRC to DEST, guaranteeing
   correct behavior for overlapping strings.  */
extern wchar_t *wmemmove (wchar_t *__s1, __const wchar_t *__s2, size_t __n)
     __THROW;

/* Set N wide characters of S to C.  */
extern wchar_t *wmemset (wchar_t *__s, wchar_t __c, size_t __n) __THROW;
__END_NAMESPACE_STD

#ifdef __USE_GNU
/* Copy N wide characters of SRC to DEST and return pointer to following
   wide character.  */
extern wchar_t *wmempcpy (wchar_t *__restrict __s1,
			  __const wchar_t *__restrict __s2, size_t __n)
     __THROW;
#endif


__BEGIN_NAMESPACE_STD
/* Determine whether C constitutes a valid (one-byte) multibyte
   character.  */
extern wint_t btowc (int __c) __THROW;

/* Determine whether C corresponds to a member of the extended
   character set whose multibyte representation is a single byte.  */
extern int wctob (wint_t __c) __THROW;

/* Determine whether PS points to an object representing the initial
   state.  */
extern int mbsinit (__const mbstate_t *__ps) __THROW __attribute_pure__;

/* Write wide character representation of multibyte character pointed
   to by S to PWC.  */
extern size_t mbrtowc (wchar_t *__restrict __pwc,
		       __const char *__restrict __s, size_t __n,
		       mbstate_t *__p) __THROW;

/* Write multibyte representation of wide character WC to S.  */
extern size_t wcrtomb (char *__restrict __s, wchar_t __wc,
		       mbstate_t *__restrict __ps) __THROW;

/* Return number of bytes in multibyte character pointed to by S.  */
extern size_t __mbrlen (__const char *__restrict __s, size_t __n,
			mbstate_t *__restrict __ps) __THROW;
extern size_t mbrlen (__const char *__restrict __s, size_t __n,
		      mbstate_t *__restrict __ps) __THROW;
__END_NAMESPACE_STD

#ifdef __USE_EXTERN_INLINES
/* Define inline function as optimization.  */

/* We can use the BTOWC and WCTOB optimizations since we know that all
   locales must use ASCII encoding for the values in the ASCII range
   and because the wchar_t encoding is always ISO 10646.  */
extern wint_t __btowc_alias (int __c) __asm ("btowc");
__extern_inline wint_t
__NTH (btowc (int __c))
{ return (__builtin_constant_p (__c) && __c >= '\0' && __c <= '\x7f'
	  ? (wint_t) __c : __btowc_alias (__c)); }

extern int __wctob_alias (wint_t __c) __asm ("wctob");
__extern_inline int
__NTH (wctob (wint_t __wc))
{ return (__builtin_constant_p (__wc) && __wc >= L'\0' && __wc <= L'\x7f'
	  ? (int) __wc : __wctob_alias (__wc)); }

__extern_inline size_t
__NTH (mbrlen (__const char *__restrict __s, size_t __n,
	       mbstate_t *__restrict __ps))
{ return (__ps != NULL
	  ? mbrtowc (NULL, __s, __n, __ps) : __mbrlen (__s, __n, NULL)); }
#endif

__BEGIN_NAMESPACE_STD
/* Write wide character representation of multibyte character string
   SRC to DST.  */
extern size_t mbsrtowcs (wchar_t *__restrict __dst,
			 __const char **__restrict __src, size_t __len,
			 mbstate_t *__restrict __ps) __THROW;

/* Write multibyte character representation of wide character string
   SRC to DST.  */
extern size_t wcsrtombs (char *__restrict __dst,
			 __const wchar_t **__restrict __src, size_t __len,
			 mbstate_t *__restrict __ps) __THROW;
__END_NAMESPACE_STD


#ifdef	__USE_XOPEN2K8
/* Write wide character representation of at most NMC bytes of the
   multibyte character string SRC to DST.  */
extern size_t mbsnrtowcs (wchar_t *__restrict __dst,
			  __const char **__restrict __src, size_t __nmc,
			  size_t __len, mbstate_t *__restrict __ps) __THROW;

/* Write multibyte character representation of at most NWC characters
   from the wide character string SRC to DST.  */
extern size_t wcsnrtombs (char *__restrict __dst,
			  __const wchar_t **__restrict __src,
			  size_t __nwc, size_t __len,
			  mbstate_t *__restrict __ps) __THROW;
#endif	/* use POSIX 2008 */


/* The following functions are extensions found in X/Open CAE.  */
#ifdef __USE_XOPEN
/* Determine number of column positions required for C.  */
extern int wcwidth (wchar_t __c) __THROW;

/* Determine number of column positions required for first N wide
   characters (or fewer if S ends before this) in S.  */
extern int wcswidth (__const wchar_t *__s, size_t __n) __THROW;
#endif	/* Use X/Open.  */


__BEGIN_NAMESPACE_STD
/* Convert initial portion of the wide string NPTR to `double'
   representation.  */
extern double wcstod (__const wchar_t *__restrict __nptr,
		      wchar_t **__restrict __endptr) __THROW;
__END_NAMESPACE_STD

#ifdef __USE_ISOC99
__BEGIN_NAMESPACE_C99
/* Likewise for `float' and `long double' sizes of floating-point numbers.  */
extern float wcstof (__const wchar_t *__restrict __nptr,
		     wchar_t **__restrict __endptr) __THROW;
extern long double wcstold (__const wchar_t *__restrict __nptr,
			    wchar_t **__restrict __endptr) __THROW;
__END_NAMESPACE_C99
#endif /* C99 */


__BEGIN_NAMESPACE_STD
/* Convert initial portion of wide string NPTR to `long int'
   representation.  */
extern long int wcstol (__const wchar_t *__restrict __nptr,
			wchar_t **__restrict __endptr, int __base) __THROW;

/* Convert initial portion of wide string NPTR to `unsigned long int'
   representation.  */
extern unsigned long int wcstoul (__const wchar_t *__restrict __nptr,
				  wchar_t **__restrict __endptr, int __base)
     __THROW;
__END_NAMESPACE_STD

#if defined __USE_ISOC99 || (defined __GNUC__ && defined __USE_GNU)
__BEGIN_NAMESPACE_C99
/* Convert initial portion of wide string NPTR to `long long int'
   representation.  */
__extension__
extern long long int wcstoll (__const wchar_t *__restrict __nptr,
			      wchar_t **__restrict __endptr, int __base)
     __THROW;

/* Convert initial portion of wide string NPTR to `unsigned long long int'
   representation.  */
__extension__
extern unsigned long long int wcstoull (__const wchar_t *__restrict __nptr,
					wchar_t **__restrict __endptr,
					int __base) __THROW;
__END_NAMESPACE_C99
#endif /* ISO C99 or GCC and GNU.  */

#if defined __GNUC__ && defined __USE_GNU
/* Convert initial portion of wide string NPTR to `long long int'
   representation.  */
__extension__
extern long long int wcstoq (__const wchar_t *__restrict __nptr,
			     wchar_t **__restrict __endptr, int __base)
     __THROW;

/* Convert initial portion of wide string NPTR to `unsigned long long int'
   representation.  */
__extension__
extern unsigned long long int wcstouq (__const wchar_t *__restrict __nptr,
				       wchar_t **__restrict __endptr,
				       int __base) __THROW;
#endif /* GCC and use GNU.  */

#ifdef __USE_GNU
/* The concept of one static locale per category is not very well
   thought out.  Many applications will need to process its data using
   information from several different locales.  Another application is
   the implementation of the internationalization handling in the
   upcoming ISO C++ standard library.  To support this another set of
   the functions using locale data exist which have an additional
   argument.

   Attention: all these functions are *not* standardized in any form.
   This is a proof-of-concept implementation.  */

/* Structure for reentrant locale using functions.  This is an
   (almost) opaque type for the user level programs.  */
# include <xlocale.h>

/* Special versions of the functions above which take the locale to
   use as an additional parameter.  */
extern long int wcstol_l (__const wchar_t *__restrict __nptr,
			  wchar_t **__restrict __endptr, int __base,
			  __locale_t __loc) __THROW;

extern unsigned long int wcstoul_l (__const wchar_t *__restrict __nptr,
				    wchar_t **__restrict __endptr,
				    int __base, __locale_t __loc) __THROW;

__extension__
extern long long int wcstoll_l (__const wchar_t *__restrict __nptr,
				wchar_t **__restrict __endptr,
				int __base, __locale_t __loc) __THROW;

__extension__
extern unsigned long long int wcstoull_l (__const wchar_t *__restrict __nptr,
					  wchar_t **__restrict __endptr,
					  int __base, __locale_t __loc)
     __THROW;

extern double wcstod_l (__const wchar_t *__restrict __nptr,
			wchar_t **__restrict __endptr, __locale_t __loc)
     __THROW;

extern float wcstof_l (__const wchar_t *__restrict __nptr,
		       wchar_t **__restrict __endptr, __locale_t __loc)
     __THROW;

extern long double wcstold_l (__const wchar_t *__restrict __nptr,
			      wchar_t **__restrict __endptr,
			      __locale_t __loc) __THROW;


/* Copy SRC to DEST, returning the address of the terminating L'\0' in
   DEST.  */
extern wchar_t *wcpcpy (wchar_t *__restrict __dest,
			__const wchar_t *__restrict __src) __THROW;

/* Copy no more than N characters of SRC to DEST, returning the address of
   the last character written into DEST.  */
extern wchar_t *wcpncpy (wchar_t *__restrict __dest,
			 __const wchar_t *__restrict __src, size_t __n)
     __THROW;
#endif	/* use GNU */


/* Wide character I/O functions.  */

#ifdef	__USE_XOPEN2K8
/* Like OPEN_MEMSTREAM, but the stream is wide oriented and produces
   a wide character string.  */
extern __FILE *open_wmemstream (wchar_t **__bufloc, size_t *__sizeloc) __THROW;
#endif

#if defined __USE_ISOC95 || defined __USE_UNIX98
__BEGIN_NAMESPACE_STD

/* Select orientation for stream.  */
extern int fwide (__FILE *__fp, int __mode) __THROW;


/* Write formatted output to STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fwprintf (__FILE *__restrict __stream,
		     __const wchar_t *__restrict __format, ...)
     /* __attribute__ ((__format__ (__wprintf__, 2, 3))) */;
/* Write formatted output to stdout.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int wprintf (__const wchar_t *__restrict __format, ...)
     /* __attribute__ ((__format__ (__wprintf__, 1, 2))) */;
/* Write formatted output of at most N characters to S.  */
extern int swprintf (wchar_t *__restrict __s, size_t __n,
		     __const wchar_t *__restrict __format, ...)
     __THROW /* __attribute__ ((__format__ (__wprintf__, 3, 4))) */;

/* Write formatted output to S from argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vfwprintf (__FILE *__restrict __s,
		      __const wchar_t *__restrict __format,
		      __gnuc_va_list __arg)
     /* __attribute__ ((__format__ (__wprintf__, 2, 0))) */;
/* Write formatted output to stdout from argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vwprintf (__const wchar_t *__restrict __format,
		     __gnuc_va_list __arg)
     /* __attribute__ ((__format__ (__wprintf__, 1, 0))) */;
/* Write formatted output of at most N character to S from argument
   list ARG.  */
extern int vswprintf (wchar_t *__restrict __s, size_t __n,
		      __const wchar_t *__restrict __format,
		      __gnuc_va_list __arg)
     __THROW /* __attribute__ ((__format__ (__wprintf__, 3, 0))) */;


/* Read formatted input from STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fwscanf (__FILE *__restrict __stream,
		    __const wchar_t *__restrict __format, ...)
     /* __attribute__ ((__format__ (__wscanf__, 2, 3))) */;
/* Read formatted input from stdin.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int wscanf (__const wchar_t *__restrict __format, ...)
     /* __attribute__ ((__format__ (__wscanf__, 1, 2))) */;
/* Read formatted input from S.  */
extern int swscanf (__const wchar_t *__restrict __s,
		    __const wchar_t *__restrict __format, ...)
     __THROW /* __attribute__ ((__format__ (__wscanf__, 2, 3))) */;

# if defined __USE_ISOC99 && !defined __USE_GNU \
     && (!defined __LDBL_COMPAT || !defined __REDIRECT) \
     && (defined __STRICT_ANSI__ || defined __USE_XOPEN2K)
#  ifdef __REDIRECT
/* For strict ISO C99 or POSIX compliance disallow %as, %aS and %a[
   GNU extension which conflicts with valid %a followed by letter
   s, S or [.  */
extern int __REDIRECT (fwscanf, (__FILE *__restrict __stream,
				 __const wchar_t *__restrict __format, ...),
		       __isoc99_fwscanf)
     /* __attribute__ ((__format__ (__wscanf__, 2, 3))) */;
extern int __REDIRECT (wscanf, (__const wchar_t *__restrict __format, ...),
		       __isoc99_wscanf)
     /* __attribute__ ((__format__ (__wscanf__, 1, 2))) */;
extern int __REDIRECT (swscanf, (__const wchar_t *__restrict __s,
				 __const wchar_t *__restrict __format, ...),
		       __isoc99_swscanf)
     __THROW /* __attribute__ ((__format__ (__wscanf__, 2, 3))) */;
#  else
extern int __isoc99_fwscanf (__FILE *__restrict __stream,
			     __const wchar_t *__restrict __format, ...);
extern int __isoc99_wscanf (__const wchar_t *__restrict __format, ...);
extern int __isoc99_swscanf (__const wchar_t *__restrict __s,
			     __const wchar_t *__restrict __format, ...)
     __THROW;
#   define fwscanf __isoc99_fwscanf
#   define wscanf __isoc99_wscanf
#   define swscanf __isoc99_swscanf
#  endif
# endif

__END_NAMESPACE_STD
#endif /* Use ISO C95, C99 and Unix98. */

#ifdef __USE_ISOC99
__BEGIN_NAMESPACE_C99
/* Read formatted input from S into argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vfwscanf (__FILE *__restrict __s,
		     __const wchar_t *__restrict __format,
		     __gnuc_va_list __arg)
     /* __attribute__ ((__format__ (__wscanf__, 2, 0))) */;
/* Read formatted input from stdin into argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vwscanf (__const wchar_t *__restrict __format,
		    __gnuc_va_list __arg)
     /* __attribute__ ((__format__ (__wscanf__, 1, 0))) */;
/* Read formatted input from S into argument list ARG.  */
extern int vswscanf (__const wchar_t *__restrict __s,
		     __const wchar_t *__restrict __format,
		     __gnuc_va_list __arg)
     __THROW /* __attribute__ ((__format__ (__wscanf__, 2, 0))) */;

# if !defined __USE_GNU \
     && (!defined __LDBL_COMPAT || !defined __REDIRECT) \
     && (defined __STRICT_ANSI__ || defined __USE_XOPEN2K)
#  ifdef __REDIRECT
extern int __REDIRECT (vfwscanf, (__FILE *__restrict __s,
				  __const wchar_t *__restrict __format,
				  __gnuc_va_list __arg), __isoc99_vfwscanf)
     /* __attribute__ ((__format__ (__wscanf__, 2, 0))) */;
extern int __REDIRECT (vwscanf, (__const wchar_t *__restrict __format,
				 __gnuc_va_list __arg), __isoc99_vwscanf)
     /* __attribute__ ((__format__ (__wscanf__, 1, 0))) */;
extern int __REDIRECT (vswscanf, (__const wchar_t *__restrict __s,
				  __const wchar_t *__restrict __format,
				  __gnuc_va_list __arg), __isoc99_vswscanf)
     __THROW /* __attribute__ ((__format__ (__wscanf__, 2, 0))) */;
#  else
extern int __isoc99_vfwscanf (__FILE *__restrict __s,
			      __const wchar_t *__restrict __format,
			      __gnuc_va_list __arg);
extern int __isoc99_vwscanf (__const wchar_t *__restrict __format,
			     __gnuc_va_list __arg);
extern int __isoc99_vswscanf (__const wchar_t *__restrict __s,
			      __const wchar_t *__restrict __format,
			      __gnuc_va_list __arg) __THROW;
#   define vfwscanf __isoc99_vfwscanf
#   define vwscanf __isoc99_vwscanf
#   define vswscanf __isoc99_vswscanf
#  endif
# endif

__END_NAMESPACE_C99
#endif /* Use ISO C99. */


__BEGIN_NAMESPACE_STD
/* Read a character from STREAM.

   These functions are possible cancellation points and therefore not
   marked with __THROW.  */
extern wint_t fgetwc (__FILE *__stream);
extern wint_t getwc (__FILE *__stream);

/* Read a character from stdin.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern wint_t getwchar (void);


/* Write a character to STREAM.

   These functions are possible cancellation points and therefore not
   marked with __THROW.  */
extern wint_t fputwc (wchar_t __wc, __FILE *__stream);
extern wint_t putwc (wchar_t __wc, __FILE *__stream);

/* Write a character to stdout.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern wint_t putwchar (wchar_t __wc);


/* Get a newline-terminated wide character string of finite length
   from STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern wchar_t *fgetws (wchar_t *__restrict __ws, int __n,
			__FILE *__restrict __stream);

/* Write a string to STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern int fputws (__const wchar_t *__restrict __ws,
		   __FILE *__restrict __stream);


/* Push a character back onto the input buffer of STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern wint_t ungetwc (wint_t __wc, __FILE *__stream);
__END_NAMESPACE_STD


#ifdef __USE_GNU
/* These are defined to be equivalent to the `char' functions defined
   in POSIX.1:1996.

   These functions are not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation they are cancellation points and
   therefore not marked with __THROW.  */
extern wint_t getwc_unlocked (__FILE *__stream);
extern wint_t getwchar_unlocked (void);

/* This is the wide character version of a GNU extension.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern wint_t fgetwc_unlocked (__FILE *__stream);

/* Faster version when locking is not necessary.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern wint_t fputwc_unlocked (wchar_t __wc, __FILE *__stream);

/* These are defined to be equivalent to the `char' functions defined
   in POSIX.1:1996.

   These functions are not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation they are cancellation points and
   therefore not marked with __THROW.  */
extern wint_t putwc_unlocked (wchar_t __wc, __FILE *__stream);
extern wint_t putwchar_unlocked (wchar_t __wc);


/* This function does the same as `fgetws' but does not lock the stream.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern wchar_t *fgetws_unlocked (wchar_t *__restrict __ws, int __n,
				 __FILE *__restrict __stream);

/* This function does the same as `fputws' but does not lock the stream.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern int fputws_unlocked (__const wchar_t *__restrict __ws,
			    __FILE *__restrict __stream);
#endif


__BEGIN_NAMESPACE_C99
/* Format TP into S according to FORMAT.
   Write no more than MAXSIZE wide characters and return the number
   of wide characters written, or 0 if it would exceed MAXSIZE.  */
extern size_t wcsftime (wchar_t *__restrict __s, size_t __maxsize,
			__const wchar_t *__restrict __format,
			__const struct tm *__restrict __tp) __THROW;
__END_NAMESPACE_C99

# ifdef __USE_GNU
# include <xlocale.h>

/* Similar to `wcsftime' but takes the information from
   the provided locale and not the global locale.  */
extern size_t wcsftime_l (wchar_t *__restrict __s, size_t __maxsize,
			  __const wchar_t *__restrict __format,
			  __const struct tm *__restrict __tp,
			  __locale_t __loc) __THROW;
# endif

/* The X/Open standard demands that most of the functions defined in
   the <wctype.h> header must also appear here.  This is probably
   because some X/Open members wrote their implementation before the
   ISO C standard was published and introduced the better solution.
   We have to provide these definitions for compliance reasons but we
   do this nonsense only if really necessary.  */
#if defined __USE_UNIX98 && !defined __USE_GNU
# define __need_iswxxx
# include <wctype.h>
#endif

/* Define some macros helping to catch buffer overflows.  */
#if __USE_FORTIFY_LEVEL > 0 && defined __extern_always_inline
# include <bits/wchar2.h>
#endif

#ifdef __LDBL_COMPAT
# include <bits/wchar-ldbl.h>
#endif

__END_DECLS

#endif	/* _WCHAR_H defined */

#endif /* wchar.h  */

/* Undefine all __need_* constants in case we are included to get those
   constants but the whole file was already read.  */
#undef __need_mbstate_t
#undef __need_wint_t
