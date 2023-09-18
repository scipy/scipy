/* Machine-independant string function optimizations.
   Copyright (C) 1997-2003, 2004, 2007, 2008 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@cygnus.com>, 1997.

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
# error "Never use <bits/string2.h> directly; include <string.h> instead."
#endif

#if !defined __NO_STRING_INLINES && !defined __BOUNDED_POINTERS__

/* Unlike the definitions in the header <bits/string.h> the
   definitions contained here are not optimized down to assembler
   level.  Those optimizations are not always a good idea since this
   means the code size increases a lot.  Instead the definitions here
   optimize some functions in a way which do not dramatically
   increase the code size and which do not use assembler.  The main
   trick is to use GCC's `__builtin_constant_p' function.

   Every function XXX which has a defined version in
   <bits/string.h> must be accompanied by a symbol _HAVE_STRING_ARCH_XXX
   to make sure we don't get redefinitions.

   We must use here macros instead of inline functions since the
   trick won't work with the latter.  */

#ifndef __STRING_INLINE
# ifdef __cplusplus
#  define __STRING_INLINE inline
# else
#  define __STRING_INLINE __extern_inline
# endif
#endif

#if _STRING_ARCH_unaligned
/* If we can do unaligned memory accesses we must know the endianess.  */
# include <endian.h>
# include <bits/types.h>

# if __BYTE_ORDER == __LITTLE_ENDIAN
#  define __STRING2_SMALL_GET16(src, idx) \
     (((__const unsigned char *) (__const char *) (src))[idx + 1] << 8	      \
      | ((__const unsigned char *) (__const char *) (src))[idx])
#  define __STRING2_SMALL_GET32(src, idx) \
     (((((__const unsigned char *) (__const char *) (src))[idx + 3] << 8      \
	| ((__const unsigned char *) (__const char *) (src))[idx + 2]) << 8   \
       | ((__const unsigned char *) (__const char *) (src))[idx + 1]) << 8    \
      | ((__const unsigned char *) (__const char *) (src))[idx])
# else
#  define __STRING2_SMALL_GET16(src, idx) \
     (((__const unsigned char *) (__const char *) (src))[idx] << 8	      \
      | ((__const unsigned char *) (__const char *) (src))[idx + 1])
#  define __STRING2_SMALL_GET32(src, idx) \
     (((((__const unsigned char *) (__const char *) (src))[idx] << 8	      \
	| ((__const unsigned char *) (__const char *) (src))[idx + 1]) << 8   \
       | ((__const unsigned char *) (__const char *) (src))[idx + 2]) << 8    \
      | ((__const unsigned char *) (__const char *) (src))[idx + 3])
# endif
#else
/* These are a few types we need for the optimizations if we cannot
   use unaligned memory accesses.  */
# define __STRING2_COPY_TYPE(N) \
  typedef struct { unsigned char __arr[N]; }				      \
    __attribute__ ((__packed__)) __STRING2_COPY_ARR##N
__STRING2_COPY_TYPE (2);
__STRING2_COPY_TYPE (3);
__STRING2_COPY_TYPE (4);
__STRING2_COPY_TYPE (5);
__STRING2_COPY_TYPE (6);
__STRING2_COPY_TYPE (7);
__STRING2_COPY_TYPE (8);
# undef __STRING2_COPY_TYPE
#endif

/* Dereferencing a pointer arg to run sizeof on it fails for the void
   pointer case, so we use this instead.
   Note that __x is evaluated twice. */
#define __string2_1bptr_p(__x) \
  ((size_t)(const void *)((__x) + 1) - (size_t)(const void *)(__x) == 1)

/* Set N bytes of S to C.  */
#if !defined _HAVE_STRING_ARCH_memset
# if !__GNUC_PREREQ (3, 0)
#  if _STRING_ARCH_unaligned
#   define memset(s, c, n) \
  (__extension__ (__builtin_constant_p (n) && (n) <= 16			      \
		  ? ((n) == 1						      \
		     ? __memset_1 (s, c)				      \
		     : __memset_gc (s, c, n))				      \
		  : (__builtin_constant_p (c) && (c) == '\0'		      \
		     ? ({ void *__s = (s); __bzero (__s, n); __s; })	      \
		     : memset (s, c, n))))

#   define __memset_1(s, c) ({ void *__s = (s);				      \
			    *((__uint8_t *) __s) = (__uint8_t) c; __s; })

#   define __memset_gc(s, c, n) \
  ({ void *__s = (s);							      \
     union {								      \
       unsigned int __ui;						      \
       unsigned short int __usi;					      \
       unsigned char __uc;						      \
     } *__u = __s;							      \
     __uint8_t __c = (__uint8_t) (c);					      \
									      \
     /* This `switch' statement will be removed at compile-time.  */	      \
     switch ((unsigned int) (n))					      \
       {								      \
       case 15:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 11:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 7:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 3:								      \
	 __u->__usi = (unsigned short int) __c * 0x0101;		      \
	 __u = __extension__ ((void *) __u + 2);			      \
	 __u->__uc = (unsigned char) __c;				      \
	 break;								      \
									      \
       case 14:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 10:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 6:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 2:								      \
	 __u->__usi = (unsigned short int) __c * 0x0101;		      \
	 break;								      \
									      \
       case 13:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 9:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 5:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 1:								      \
	 __u->__uc = (unsigned char) __c;				      \
	 break;								      \
									      \
       case 16:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 12:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 8:								      \
	 __u->__ui = __c * 0x01010101;					      \
	 __u = __extension__ ((void *) __u + 4);			      \
       case 4:								      \
	 __u->__ui = __c * 0x01010101;					      \
       case 0:								      \
	 break;								      \
       }								      \
									      \
     __s; })
#  else
#   define memset(s, c, n) \
  (__extension__ (__builtin_constant_p (c) && (c) == '\0'		      \
		  ? ({ void *__s = (s); __bzero (__s, n); __s; })	      \
		  : memset (s, c, n)))
#  endif
# endif

/* GCC < 3.0 optimizes memset(s, 0, n) but not bzero(s, n).
   The optimization is broken before EGCS 1.1.
   GCC 3.0+ has __builtin_bzero as well, but at least till GCC 3.4
   if it decides to call the library function, it calls memset
   and not bzero.  */
# if __GNUC_PREREQ (2, 91)
#  define __bzero(s, n) __builtin_memset (s, '\0', n)
# endif

#endif


/* Copy N bytes from SRC to DEST, returning pointer to byte following the
   last copied.  */
#ifdef __USE_GNU
# if !defined _HAVE_STRING_ARCH_mempcpy || defined _FORCE_INLINES
#  ifndef _HAVE_STRING_ARCH_mempcpy
#   if __GNUC_PREREQ (3, 4)
#    define __mempcpy(dest, src, n) __builtin_mempcpy (dest, src, n)
#   elif __GNUC_PREREQ (3, 0)
#    define __mempcpy(dest, src, n) \
  (__extension__ (__builtin_constant_p (src) && __builtin_constant_p (n)      \
		  && __string2_1bptr_p (src) && n <= 8			      \
		  ? __builtin_memcpy (dest, src, n) + (n)		      \
		  : __mempcpy (dest, src, n)))
#   else
#    define __mempcpy(dest, src, n) \
  (__extension__ (__builtin_constant_p (src) && __builtin_constant_p (n)      \
		  && __string2_1bptr_p (src) && n <= 8			      \
		  ? __mempcpy_small (dest, __mempcpy_args (src), n)	      \
		  : __mempcpy (dest, src, n)))
#   endif
/* In glibc we use this function frequently but for namespace reasons
   we have to use the name `__mempcpy'.  */
#   define mempcpy(dest, src, n) __mempcpy (dest, src, n)
#  endif

#  if !__GNUC_PREREQ (3, 0) || defined _FORCE_INLINES
#   if _STRING_ARCH_unaligned
#    ifndef _FORCE_INLINES
#     define __mempcpy_args(src) \
     ((__const char *) (src))[0], ((__const char *) (src))[2],		      \
     ((__const char *) (src))[4], ((__const char *) (src))[6],		      \
     __extension__ __STRING2_SMALL_GET16 (src, 0),			      \
     __extension__ __STRING2_SMALL_GET16 (src, 4),			      \
     __extension__ __STRING2_SMALL_GET32 (src, 0),			      \
     __extension__ __STRING2_SMALL_GET32 (src, 4)
#    endif
__STRING_INLINE void *__mempcpy_small (void *, char, char, char, char,
				       __uint16_t, __uint16_t, __uint32_t,
				       __uint32_t, size_t);
__STRING_INLINE void *
__mempcpy_small (void *__dest1,
		 char __src0_1, char __src2_1, char __src4_1, char __src6_1,
		 __uint16_t __src0_2, __uint16_t __src4_2,
		 __uint32_t __src0_4, __uint32_t __src4_4,
		 size_t __srclen)
{
  union {
    __uint32_t __ui;
    __uint16_t __usi;
    unsigned char __uc;
    unsigned char __c;
  } *__u = __dest1;
  switch ((unsigned int) __srclen)
    {
    case 1:
      __u->__c = __src0_1;
      __u = __extension__ ((void *) __u + 1);
      break;
    case 2:
      __u->__usi = __src0_2;
      __u = __extension__ ((void *) __u + 2);
      break;
    case 3:
      __u->__usi = __src0_2;
      __u = __extension__ ((void *) __u + 2);
      __u->__c = __src2_1;
      __u = __extension__ ((void *) __u + 1);
      break;
    case 4:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      break;
    case 5:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__c = __src4_1;
      __u = __extension__ ((void *) __u + 1);
      break;
    case 6:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__usi = __src4_2;
      __u = __extension__ ((void *) __u + 2);
      break;
    case 7:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__usi = __src4_2;
      __u = __extension__ ((void *) __u + 2);
      __u->__c = __src6_1;
      __u = __extension__ ((void *) __u + 1);
      break;
    case 8:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__ui = __src4_4;
      __u = __extension__ ((void *) __u + 4);
      break;
    }
  return (void *) __u;
}
#   else
#    ifndef _FORCE_INLINES
#     define __mempcpy_args(src) \
     ((__const char *) (src))[0],					      \
     __extension__ ((__STRING2_COPY_ARR2)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1] } }),      \
     __extension__ ((__STRING2_COPY_ARR3)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2] } }),				      \
     __extension__ ((__STRING2_COPY_ARR4)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3] } }),      \
     __extension__ ((__STRING2_COPY_ARR5)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4] } }),				      \
     __extension__ ((__STRING2_COPY_ARR6)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5] } }),      \
     __extension__ ((__STRING2_COPY_ARR7)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5],	      \
	  ((__const char *) (src))[6] } }),				      \
     __extension__ ((__STRING2_COPY_ARR8)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5],	      \
	  ((__const char *) (src))[6], ((__const char *) (src))[7] } })
#    endif
__STRING_INLINE void *__mempcpy_small (void *, char, __STRING2_COPY_ARR2,
				       __STRING2_COPY_ARR3,
				       __STRING2_COPY_ARR4,
				       __STRING2_COPY_ARR5,
				       __STRING2_COPY_ARR6,
				       __STRING2_COPY_ARR7,
				       __STRING2_COPY_ARR8, size_t);
__STRING_INLINE void *
__mempcpy_small (void *__dest, char __src1,
		 __STRING2_COPY_ARR2 __src2, __STRING2_COPY_ARR3 __src3,
		 __STRING2_COPY_ARR4 __src4, __STRING2_COPY_ARR5 __src5,
		 __STRING2_COPY_ARR6 __src6, __STRING2_COPY_ARR7 __src7,
		 __STRING2_COPY_ARR8 __src8, size_t __srclen)
{
  union {
    char __c;
    __STRING2_COPY_ARR2 __sca2;
    __STRING2_COPY_ARR3 __sca3;
    __STRING2_COPY_ARR4 __sca4;
    __STRING2_COPY_ARR5 __sca5;
    __STRING2_COPY_ARR6 __sca6;
    __STRING2_COPY_ARR7 __sca7;
    __STRING2_COPY_ARR8 __sca8;
  } *__u = __dest;
  switch ((unsigned int) __srclen)
    {
    case 1:
      __u->__c = __src1;
      break;
    case 2:
      __extension__ __u->__sca2 = __src2;
      break;
    case 3:
      __extension__ __u->__sca3 = __src3;
      break;
    case 4:
      __extension__ __u->__sca4 = __src4;
      break;
    case 5:
      __extension__ __u->__sca5 = __src5;
      break;
    case 6:
      __extension__ __u->__sca6 = __src6;
      break;
    case 7:
      __extension__ __u->__sca7 = __src7;
      break;
    case 8:
      __extension__ __u->__sca8 = __src8;
      break;
    }
  return __extension__ ((void *) __u + __srclen);
}
#   endif
#  endif
# endif
#endif


/* Return pointer to C in S.  */
#ifndef _HAVE_STRING_ARCH_strchr
extern void *__rawmemchr (const void *__s, int __c);
# if __GNUC_PREREQ (3, 2)
#  define strchr(s, c) \
  (__extension__ (__builtin_constant_p (c) && !__builtin_constant_p (s)	      \
		  && (c) == '\0'					      \
		  ? (char *) __rawmemchr (s, c)				      \
		  : __builtin_strchr (s, c)))
# else
#  define strchr(s, c) \
  (__extension__ (__builtin_constant_p (c) && (c) == '\0'		      \
		  ? (char *) __rawmemchr (s, c)				      \
		  : strchr (s, c)))
# endif
#endif


/* Copy SRC to DEST.  */
#if (!defined _HAVE_STRING_ARCH_strcpy && !__GNUC_PREREQ (3, 0)) \
    || defined _FORCE_INLINES
# if !defined _HAVE_STRING_ARCH_strcpy && !__GNUC_PREREQ (3, 0)
#  define strcpy(dest, src) \
  (__extension__ (__builtin_constant_p (src)				      \
		  ? (__string2_1bptr_p (src) && strlen (src) + 1 <= 8	      \
		     ? __strcpy_small (dest, __strcpy_args (src),	      \
				       strlen (src) + 1)		      \
		     : (char *) memcpy (dest, src, strlen (src) + 1))	      \
		  : strcpy (dest, src)))
# endif

# if _STRING_ARCH_unaligned
#  ifndef _FORCE_INLINES
#   define __strcpy_args(src) \
     __extension__ __STRING2_SMALL_GET16 (src, 0),			      \
     __extension__ __STRING2_SMALL_GET16 (src, 4),			      \
     __extension__ __STRING2_SMALL_GET32 (src, 0),			      \
     __extension__ __STRING2_SMALL_GET32 (src, 4)
#  endif
__STRING_INLINE char *__strcpy_small (char *, __uint16_t, __uint16_t,
				      __uint32_t, __uint32_t, size_t);
__STRING_INLINE char *
__strcpy_small (char *__dest,
		__uint16_t __src0_2, __uint16_t __src4_2,
		__uint32_t __src0_4, __uint32_t __src4_4,
		size_t __srclen)
{
  union {
    __uint32_t __ui;
    __uint16_t __usi;
    unsigned char __uc;
  } *__u = (void *) __dest;
  switch ((unsigned int) __srclen)
    {
    case 1:
      __u->__uc = '\0';
      break;
    case 2:
      __u->__usi = __src0_2;
      break;
    case 3:
      __u->__usi = __src0_2;
      __u = __extension__ ((void *) __u + 2);
      __u->__uc = '\0';
      break;
    case 4:
      __u->__ui = __src0_4;
      break;
    case 5:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__uc = '\0';
      break;
    case 6:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__usi = __src4_2;
      break;
    case 7:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__usi = __src4_2;
      __u = __extension__ ((void *) __u + 2);
      __u->__uc = '\0';
      break;
    case 8:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__ui = __src4_4;
      break;
    }
  return __dest;
}
# else
#  ifndef _FORCE_INLINES
#   define __strcpy_args(src) \
     __extension__ ((__STRING2_COPY_ARR2)				      \
      { { ((__const char *) (src))[0], '\0' } }),			      \
     __extension__ ((__STRING2_COPY_ARR3)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  '\0' } }),							      \
     __extension__ ((__STRING2_COPY_ARR4)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], '\0' } }),			      \
     __extension__ ((__STRING2_COPY_ARR5)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  '\0' } }),							      \
     __extension__ ((__STRING2_COPY_ARR6)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], '\0' } }),			      \
     __extension__ ((__STRING2_COPY_ARR7)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5],	      \
	  '\0' } }),							      \
     __extension__ ((__STRING2_COPY_ARR8)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5],	      \
	  ((__const char *) (src))[6], '\0' } })
#  endif
__STRING_INLINE char *__strcpy_small (char *, __STRING2_COPY_ARR2,
				      __STRING2_COPY_ARR3,
				      __STRING2_COPY_ARR4,
				      __STRING2_COPY_ARR5,
				      __STRING2_COPY_ARR6,
				      __STRING2_COPY_ARR7,
				      __STRING2_COPY_ARR8, size_t);
__STRING_INLINE char *
__strcpy_small (char *__dest,
		__STRING2_COPY_ARR2 __src2, __STRING2_COPY_ARR3 __src3,
		__STRING2_COPY_ARR4 __src4, __STRING2_COPY_ARR5 __src5,
		__STRING2_COPY_ARR6 __src6, __STRING2_COPY_ARR7 __src7,
		__STRING2_COPY_ARR8 __src8, size_t __srclen)
{
  union {
    char __c;
    __STRING2_COPY_ARR2 __sca2;
    __STRING2_COPY_ARR3 __sca3;
    __STRING2_COPY_ARR4 __sca4;
    __STRING2_COPY_ARR5 __sca5;
    __STRING2_COPY_ARR6 __sca6;
    __STRING2_COPY_ARR7 __sca7;
    __STRING2_COPY_ARR8 __sca8;
  } *__u = (void *) __dest;
  switch ((unsigned int) __srclen)
    {
    case 1:
      __u->__c = '\0';
      break;
    case 2:
      __extension__ __u->__sca2 = __src2;
      break;
    case 3:
      __extension__ __u->__sca3 = __src3;
      break;
    case 4:
      __extension__ __u->__sca4 = __src4;
      break;
    case 5:
      __extension__ __u->__sca5 = __src5;
      break;
    case 6:
      __extension__ __u->__sca6 = __src6;
      break;
    case 7:
      __extension__ __u->__sca7 = __src7;
      break;
    case 8:
      __extension__ __u->__sca8 = __src8;
      break;
  }
  return __dest;
}
# endif
#endif


/* Copy SRC to DEST, returning pointer to final NUL byte.  */
#ifdef __USE_GNU
# if !defined _HAVE_STRING_ARCH_stpcpy || defined _FORCE_INLINES
#  ifndef _HAVE_STRING_ARCH_stpcpy
#   if __GNUC_PREREQ (3, 4)
#    define __stpcpy(dest, src) __builtin_stpcpy (dest, src)
#   elif __GNUC_PREREQ (3, 0)
#    define __stpcpy(dest, src) \
  (__extension__ (__builtin_constant_p (src)				      \
		  ? (__string2_1bptr_p (src) && strlen (src) + 1 <= 8	      \
		     ? __builtin_strcpy (dest, src) + strlen (src)	      \
		     : ((char *) (__mempcpy) (dest, src, strlen (src) + 1)    \
			- 1))						      \
		  : __stpcpy (dest, src)))
#   else
#    define __stpcpy(dest, src) \
  (__extension__ (__builtin_constant_p (src)				      \
		  ? (__string2_1bptr_p (src) && strlen (src) + 1 <= 8	      \
		     ? __stpcpy_small (dest, __stpcpy_args (src),	      \
				       strlen (src) + 1)		      \
		     : ((char *) (__mempcpy) (dest, src, strlen (src) + 1)    \
			- 1))						      \
		  : __stpcpy (dest, src)))
#   endif
/* In glibc we use this function frequently but for namespace reasons
   we have to use the name `__stpcpy'.  */
#   define stpcpy(dest, src) __stpcpy (dest, src)
#  endif

#  if !__GNUC_PREREQ (3, 0) || defined _FORCE_INLINES
#   if _STRING_ARCH_unaligned
#    ifndef _FORCE_INLINES
#     define __stpcpy_args(src) \
     __extension__ __STRING2_SMALL_GET16 (src, 0),			      \
     __extension__ __STRING2_SMALL_GET16 (src, 4),			      \
     __extension__ __STRING2_SMALL_GET32 (src, 0),			      \
     __extension__ __STRING2_SMALL_GET32 (src, 4)
#    endif
__STRING_INLINE char *__stpcpy_small (char *, __uint16_t, __uint16_t,
				      __uint32_t, __uint32_t, size_t);
__STRING_INLINE char *
__stpcpy_small (char *__dest,
		__uint16_t __src0_2, __uint16_t __src4_2,
		__uint32_t __src0_4, __uint32_t __src4_4,
		size_t __srclen)
{
  union {
    unsigned int __ui;
    unsigned short int __usi;
    unsigned char __uc;
    char __c;
  } *__u = (void *) __dest;
  switch ((unsigned int) __srclen)
    {
    case 1:
      __u->__uc = '\0';
      break;
    case 2:
      __u->__usi = __src0_2;
      __u = __extension__ ((void *) __u + 1);
      break;
    case 3:
      __u->__usi = __src0_2;
      __u = __extension__ ((void *) __u + 2);
      __u->__uc = '\0';
      break;
    case 4:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 3);
      break;
    case 5:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__uc = '\0';
      break;
    case 6:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__usi = __src4_2;
      __u = __extension__ ((void *) __u + 1);
      break;
    case 7:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__usi = __src4_2;
      __u = __extension__ ((void *) __u + 2);
      __u->__uc = '\0';
      break;
    case 8:
      __u->__ui = __src0_4;
      __u = __extension__ ((void *) __u + 4);
      __u->__ui = __src4_4;
      __u = __extension__ ((void *) __u + 3);
      break;
    }
  return &__u->__c;
}
#   else
#    ifndef _FORCE_INLINES
#     define __stpcpy_args(src) \
     __extension__ ((__STRING2_COPY_ARR2)				      \
      { { ((__const char *) (src))[0], '\0' } }),			      \
     __extension__ ((__STRING2_COPY_ARR3)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  '\0' } }),							      \
     __extension__ ((__STRING2_COPY_ARR4)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], '\0' } }),			      \
     __extension__ ((__STRING2_COPY_ARR5)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  '\0' } }),							      \
     __extension__ ((__STRING2_COPY_ARR6)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], '\0' } }),			      \
     __extension__ ((__STRING2_COPY_ARR7)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5],	      \
	  '\0' } }),							      \
     __extension__ ((__STRING2_COPY_ARR8)				      \
      { { ((__const char *) (src))[0], ((__const char *) (src))[1],	      \
	  ((__const char *) (src))[2], ((__const char *) (src))[3],	      \
	  ((__const char *) (src))[4], ((__const char *) (src))[5],	      \
	  ((__const char *) (src))[6], '\0' } })
#    endif
__STRING_INLINE char *__stpcpy_small (char *, __STRING2_COPY_ARR2,
				      __STRING2_COPY_ARR3,
				      __STRING2_COPY_ARR4,
				      __STRING2_COPY_ARR5,
				      __STRING2_COPY_ARR6,
				      __STRING2_COPY_ARR7,
				      __STRING2_COPY_ARR8, size_t);
__STRING_INLINE char *
__stpcpy_small (char *__dest,
		__STRING2_COPY_ARR2 __src2, __STRING2_COPY_ARR3 __src3,
		__STRING2_COPY_ARR4 __src4, __STRING2_COPY_ARR5 __src5,
		__STRING2_COPY_ARR6 __src6, __STRING2_COPY_ARR7 __src7,
		__STRING2_COPY_ARR8 __src8, size_t __srclen)
{
  union {
    char __c;
    __STRING2_COPY_ARR2 __sca2;
    __STRING2_COPY_ARR3 __sca3;
    __STRING2_COPY_ARR4 __sca4;
    __STRING2_COPY_ARR5 __sca5;
    __STRING2_COPY_ARR6 __sca6;
    __STRING2_COPY_ARR7 __sca7;
    __STRING2_COPY_ARR8 __sca8;
  } *__u = (void *) __dest;
  switch ((unsigned int) __srclen)
    {
    case 1:
      __u->__c = '\0';
      break;
    case 2:
      __extension__ __u->__sca2 = __src2;
      break;
    case 3:
      __extension__ __u->__sca3 = __src3;
      break;
    case 4:
      __extension__ __u->__sca4 = __src4;
      break;
    case 5:
      __extension__ __u->__sca5 = __src5;
      break;
    case 6:
      __extension__ __u->__sca6 = __src6;
      break;
    case 7:
      __extension__ __u->__sca7 = __src7;
      break;
    case 8:
      __extension__ __u->__sca8 = __src8;
      break;
  }
  return __dest + __srclen - 1;
}
#   endif
#  endif
# endif
#endif


/* Copy no more than N characters of SRC to DEST.  */
#ifndef _HAVE_STRING_ARCH_strncpy
# if __GNUC_PREREQ (3, 2)
#  define strncpy(dest, src, n) __builtin_strncpy (dest, src, n)
# else
#  define strncpy(dest, src, n) \
  (__extension__ (__builtin_constant_p (src) && __builtin_constant_p (n)      \
		  ? (strlen (src) + 1 >= ((size_t) (n))			      \
		     ? (char *) memcpy (dest, src, n)			      \
		     : strncpy (dest, src, n))				      \
		  : strncpy (dest, src, n)))
# endif
#endif


/* Append no more than N characters from SRC onto DEST.  */
#ifndef _HAVE_STRING_ARCH_strncat
# ifdef _USE_STRING_ARCH_strchr
#  define strncat(dest, src, n) \
  (__extension__ ({ char *__dest = (dest);				      \
		    __builtin_constant_p (src) && __builtin_constant_p (n)    \
		    ? (strlen (src) < ((size_t) (n))			      \
		       ? strcat (__dest, src)				      \
		       : (*((char *) __mempcpy (strchr (__dest, '\0'),	      \
						src, n)) = '\0', __dest))     \
		    : strncat (dest, src, n); }))
# elif __GNUC_PREREQ (3, 2)
#  define strncat(dest, src, n) __builtin_strncat (dest, src, n)
# else
#  define strncat(dest, src, n) \
  (__extension__ (__builtin_constant_p (src) && __builtin_constant_p (n)      \
		  ? (strlen (src) < ((size_t) (n))			      \
		     ? strcat (dest, src)				      \
		     : strncat (dest, src, n))				      \
		  : strncat (dest, src, n)))
# endif
#endif


/* Compare characters of S1 and S2.  */
#ifndef _HAVE_STRING_ARCH_strcmp
# if __GNUC_PREREQ (3, 2)
#  define strcmp(s1, s2) \
  __extension__								      \
  ({ size_t __s1_len, __s2_len;						      \
     (__builtin_constant_p (s1) && __builtin_constant_p (s2)		      \
      && (__s1_len = strlen (s1), __s2_len = strlen (s2),		      \
	  (!__string2_1bptr_p (s1) || __s1_len >= 4)			      \
	  && (!__string2_1bptr_p (s2) || __s2_len >= 4))		      \
      ? __builtin_strcmp (s1, s2)					      \
      : (__builtin_constant_p (s1) && __string2_1bptr_p (s1)		      \
	 && (__s1_len = strlen (s1), __s1_len < 4)			      \
	 ? (__builtin_constant_p (s2) && __string2_1bptr_p (s2)		      \
	    ? __builtin_strcmp (s1, s2)					      \
	    : __strcmp_cg (s1, s2, __s1_len))				      \
	 : (__builtin_constant_p (s2) && __string2_1bptr_p (s2)		      \
	    && (__s2_len = strlen (s2), __s2_len < 4)			      \
	    ? (__builtin_constant_p (s1) && __string2_1bptr_p (s1)	      \
	       ? __builtin_strcmp (s1, s2)				      \
	       : __strcmp_gc (s1, s2, __s2_len))			      \
	    : __builtin_strcmp (s1, s2)))); })
# else
#  define strcmp(s1, s2) \
  __extension__								      \
  ({ size_t __s1_len, __s2_len;						      \
     (__builtin_constant_p (s1) && __builtin_constant_p (s2)		      \
      && (__s1_len = strlen (s1), __s2_len = strlen (s2),		      \
	  (!__string2_1bptr_p (s1) || __s1_len >= 4)			      \
	  && (!__string2_1bptr_p (s2) || __s2_len >= 4))		      \
      ? memcmp ((__const char *) (s1), (__const char *) (s2),		      \
		(__s1_len < __s2_len ? __s1_len : __s2_len) + 1)	      \
      : (__builtin_constant_p (s1) && __string2_1bptr_p (s1)		      \
	 && (__s1_len = strlen (s1), __s1_len < 4)			      \
	 ? (__builtin_constant_p (s2) && __string2_1bptr_p (s2)		      \
	    ? __strcmp_cc (s1, s2, __s1_len)				      \
	    : __strcmp_cg (s1, s2, __s1_len))				      \
	 : (__builtin_constant_p (s2) && __string2_1bptr_p (s2)		      \
	    && (__s2_len = strlen (s2), __s2_len < 4)			      \
	    ? (__builtin_constant_p (s1) && __string2_1bptr_p (s1)	      \
	       ? __strcmp_cc (s1, s2, __s2_len)				      \
	       : __strcmp_gc (s1, s2, __s2_len))			      \
	    : strcmp (s1, s2)))); })
# endif

# define __strcmp_cc(s1, s2, l) \
  (__extension__ ({ register int __result =				      \
		      (((__const unsigned char *) (__const char *) (s1))[0]   \
		       - ((__const unsigned char *) (__const char *)(s2))[0]);\
		    if (l > 0 && __result == 0)				      \
		      {							      \
			__result = (((__const unsigned char *)		      \
				     (__const char *) (s1))[1]		      \
				    - ((__const unsigned char *)	      \
				       (__const char *) (s2))[1]);	      \
			if (l > 1 && __result == 0)			      \
			  {						      \
			    __result =					      \
			      (((__const unsigned char *)		      \
				(__const char *) (s1))[2]		      \
			       - ((__const unsigned char *)		      \
				  (__const char *) (s2))[2]);		      \
			    if (l > 2 && __result == 0)			      \
			      __result =				      \
				(((__const unsigned char *)		      \
				  (__const char *) (s1))[3]		      \
				 - ((__const unsigned char *)		      \
				    (__const char *) (s2))[3]);		      \
			  }						      \
		      }							      \
		    __result; }))

# define __strcmp_cg(s1, s2, l1) \
  (__extension__ ({ __const unsigned char *__s2 =			      \
		      (__const unsigned char *) (__const char *) (s2);	      \
		    register int __result =				      \
		      (((__const unsigned char *) (__const char *) (s1))[0]   \
		       - __s2[0]);					      \
		    if (l1 > 0 && __result == 0)			      \
		      {							      \
			__result = (((__const unsigned char *)		      \
				     (__const char *) (s1))[1] - __s2[1]);    \
			if (l1 > 1 && __result == 0)			      \
			  {						      \
			    __result = (((__const unsigned char *)	      \
					 (__const char *) (s1))[2] - __s2[2]);\
			    if (l1 > 2 && __result == 0)		      \
			      __result = (((__const unsigned char *)	      \
					  (__const char *)  (s1))[3]	      \
					  - __s2[3]);			      \
			  }						      \
		      }							      \
		    __result; }))

# define __strcmp_gc(s1, s2, l2) \
  (__extension__ ({ __const unsigned char *__s1 =			      \
		      (__const unsigned char *) (__const char *) (s1);	      \
		    register int __result =				      \
		      __s1[0] - ((__const unsigned char *)		      \
				 (__const char *) (s2))[0];		      \
		    if (l2 > 0 && __result == 0)			      \
		      {							      \
			__result = (__s1[1]				      \
				    - ((__const unsigned char *)	      \
				       (__const char *) (s2))[1]);	      \
			if (l2 > 1 && __result == 0)			      \
			  {						      \
			    __result =					      \
			      (__s1[2] - ((__const unsigned char *)	      \
					  (__const char *) (s2))[2]);	      \
			    if (l2 > 2 && __result == 0)		      \
			      __result =				      \
				(__s1[3]				      \
				 - ((__const unsigned char *)		      \
				    (__const char *) (s2))[3]);		      \
			  }						      \
		      }							      \
		    __result; }))
#endif


/* Compare N characters of S1 and S2.  */
#ifndef _HAVE_STRING_ARCH_strncmp
# define strncmp(s1, s2, n)						      \
  (__extension__ (__builtin_constant_p (n)				      \
		  && ((__builtin_constant_p (s1)			      \
		       && strlen (s1) < ((size_t) (n)))			      \
		      || (__builtin_constant_p (s2)			      \
			  && strlen (s2) < ((size_t) (n))))		      \
		  ? strcmp (s1, s2) : strncmp (s1, s2, n)))
#endif


/* Return the length of the initial segment of S which
   consists entirely of characters not in REJECT.  */
#if !defined _HAVE_STRING_ARCH_strcspn || defined _FORCE_INLINES
# ifndef _HAVE_STRING_ARCH_strcspn
#  if __GNUC_PREREQ (3, 2)
#   define strcspn(s, reject) \
  __extension__								      \
  ({ char __r0, __r1, __r2;						      \
     (__builtin_constant_p (reject) && __string2_1bptr_p (reject)	      \
      ? ((__builtin_constant_p (s) && __string2_1bptr_p (s))		      \
	 ? __builtin_strcspn (s, reject)				      \
	 : ((__r0 = ((__const char *) (reject))[0], __r0 == '\0')	      \
	    ? strlen (s)						      \
	    : ((__r1 = ((__const char *) (reject))[1], __r1 == '\0')	      \
	       ? __strcspn_c1 (s, __r0)					      \
	       : ((__r2 = ((__const char *) (reject))[2], __r2 == '\0')	      \
		  ? __strcspn_c2 (s, __r0, __r1)			      \
		  : (((__const char *) (reject))[3] == '\0'		      \
		     ? __strcspn_c3 (s, __r0, __r1, __r2)		      \
		     : __builtin_strcspn (s, reject))))))		      \
      : __builtin_strcspn (s, reject)); })
#  else
#   define strcspn(s, reject) \
  __extension__								      \
  ({ char __r0, __r1, __r2;						      \
     (__builtin_constant_p (reject) && __string2_1bptr_p (reject)	      \
      ? ((__r0 = ((__const char *) (reject))[0], __r0 == '\0')		      \
	 ? strlen (s)							      \
	 : ((__r1 = ((__const char *) (reject))[1], __r1 == '\0')	      \
	    ? __strcspn_c1 (s, __r0)					      \
	    : ((__r2 = ((__const char *) (reject))[2], __r2 == '\0')	      \
	       ? __strcspn_c2 (s, __r0, __r1)				      \
	       : (((__const char *) (reject))[3] == '\0'		      \
		  ? __strcspn_c3 (s, __r0, __r1, __r2)			      \
		  : strcspn (s, reject)))))				      \
      : strcspn (s, reject)); })
#  endif
# endif

__STRING_INLINE size_t __strcspn_c1 (__const char *__s, int __reject);
__STRING_INLINE size_t
__strcspn_c1 (__const char *__s, int __reject)
{
  register size_t __result = 0;
  while (__s[__result] != '\0' && __s[__result] != __reject)
    ++__result;
  return __result;
}

__STRING_INLINE size_t __strcspn_c2 (__const char *__s, int __reject1,
				     int __reject2);
__STRING_INLINE size_t
__strcspn_c2 (__const char *__s, int __reject1, int __reject2)
{
  register size_t __result = 0;
  while (__s[__result] != '\0' && __s[__result] != __reject1
	 && __s[__result] != __reject2)
    ++__result;
  return __result;
}

__STRING_INLINE size_t __strcspn_c3 (__const char *__s, int __reject1,
				     int __reject2, int __reject3);
__STRING_INLINE size_t
__strcspn_c3 (__const char *__s, int __reject1, int __reject2,
	      int __reject3)
{
  register size_t __result = 0;
  while (__s[__result] != '\0' && __s[__result] != __reject1
	 && __s[__result] != __reject2 && __s[__result] != __reject3)
    ++__result;
  return __result;
}
#endif


/* Return the length of the initial segment of S which
   consists entirely of characters in ACCEPT.  */
#if !defined _HAVE_STRING_ARCH_strspn || defined _FORCE_INLINES
# ifndef _HAVE_STRING_ARCH_strspn
#  if __GNUC_PREREQ (3, 2)
#   define strspn(s, accept) \
  __extension__								      \
  ({ char __a0, __a1, __a2;						      \
     (__builtin_constant_p (accept) && __string2_1bptr_p (accept)	      \
      ? ((__builtin_constant_p (s) && __string2_1bptr_p (s))		      \
	 ? __builtin_strspn (s, accept)					      \
	 : ((__a0 = ((__const char *) (accept))[0], __a0 == '\0')	      \
	    ? ((void) (s), 0)						      \
	    : ((__a1 = ((__const char *) (accept))[1], __a1 == '\0')	      \
	       ? __strspn_c1 (s, __a0)					      \
	       : ((__a2 = ((__const char *) (accept))[2], __a2 == '\0')	      \
		  ? __strspn_c2 (s, __a0, __a1)				      \
		  : (((__const char *) (accept))[3] == '\0'		      \
		     ? __strspn_c3 (s, __a0, __a1, __a2)		      \
		     : __builtin_strspn (s, accept))))))		      \
      : __builtin_strspn (s, accept)); })
#  else
#   define strspn(s, accept) \
  __extension__								      \
  ({ char __a0, __a1, __a2;						      \
     (__builtin_constant_p (accept) && __string2_1bptr_p (accept)	      \
      ? ((__a0 = ((__const char *) (accept))[0], __a0 == '\0')		      \
	 ? ((void) (s), 0)						      \
	 : ((__a1 = ((__const char *) (accept))[1], __a1 == '\0')	      \
	    ? __strspn_c1 (s, __a0)					      \
	    : ((__a2 = ((__const char *) (accept))[2], __a2 == '\0')	      \
	       ? __strspn_c2 (s, __a0, __a1)				      \
	       : (((__const char *) (accept))[3] == '\0'		      \
		  ? __strspn_c3 (s, __a0, __a1, __a2)			      \
		  : strspn (s, accept)))))				      \
      : strspn (s, accept)); })
#  endif
# endif

__STRING_INLINE size_t __strspn_c1 (__const char *__s, int __accept);
__STRING_INLINE size_t
__strspn_c1 (__const char *__s, int __accept)
{
  register size_t __result = 0;
  /* Please note that __accept never can be '\0'.  */
  while (__s[__result] == __accept)
    ++__result;
  return __result;
}

__STRING_INLINE size_t __strspn_c2 (__const char *__s, int __accept1,
				    int __accept2);
__STRING_INLINE size_t
__strspn_c2 (__const char *__s, int __accept1, int __accept2)
{
  register size_t __result = 0;
  /* Please note that __accept1 and __accept2 never can be '\0'.  */
  while (__s[__result] == __accept1 || __s[__result] == __accept2)
    ++__result;
  return __result;
}

__STRING_INLINE size_t __strspn_c3 (__const char *__s, int __accept1,
				    int __accept2, int __accept3);
__STRING_INLINE size_t
__strspn_c3 (__const char *__s, int __accept1, int __accept2, int __accept3)
{
  register size_t __result = 0;
  /* Please note that __accept1 to __accept3 never can be '\0'.  */
  while (__s[__result] == __accept1 || __s[__result] == __accept2
	 || __s[__result] == __accept3)
    ++__result;
  return __result;
}
#endif


/* Find the first occurrence in S of any character in ACCEPT.  */
#if !defined _HAVE_STRING_ARCH_strpbrk || defined _FORCE_INLINES
# ifndef _HAVE_STRING_ARCH_strpbrk
#  if __GNUC_PREREQ (3, 2)
#   define strpbrk(s, accept) \
  __extension__								      \
  ({ char __a0, __a1, __a2;						      \
     (__builtin_constant_p (accept) && __string2_1bptr_p (accept)	      \
      ? ((__builtin_constant_p (s) && __string2_1bptr_p (s))		      \
	 ? __builtin_strpbrk (s, accept)				      \
	 : ((__a0 = ((__const char  *) (accept))[0], __a0 == '\0')	      \
	    ? ((void) (s), (char *) NULL)				      \
	    : ((__a1 = ((__const char *) (accept))[1], __a1 == '\0')	      \
	       ? __builtin_strchr (s, __a0)				      \
	       : ((__a2 = ((__const char *) (accept))[2], __a2 == '\0')	      \
		  ? __strpbrk_c2 (s, __a0, __a1)			      \
		  : (((__const char *) (accept))[3] == '\0'		      \
		     ? __strpbrk_c3 (s, __a0, __a1, __a2)		      \
		     : __builtin_strpbrk (s, accept))))))		      \
      : __builtin_strpbrk (s, accept)); })
#  else
#   define strpbrk(s, accept) \
  __extension__								      \
  ({ char __a0, __a1, __a2;						      \
     (__builtin_constant_p (accept) && __string2_1bptr_p (accept)	      \
      ? ((__a0 = ((__const char  *) (accept))[0], __a0 == '\0')		      \
	 ? ((void) (s), (char *) NULL)					      \
	 : ((__a1 = ((__const char *) (accept))[1], __a1 == '\0')	      \
	    ? strchr (s, __a0)						      \
	    : ((__a2 = ((__const char *) (accept))[2], __a2 == '\0')	      \
	       ? __strpbrk_c2 (s, __a0, __a1)				      \
	       : (((__const char *) (accept))[3] == '\0'		      \
		  ? __strpbrk_c3 (s, __a0, __a1, __a2)			      \
		  : strpbrk (s, accept)))))				      \
      : strpbrk (s, accept)); })
#  endif
# endif

__STRING_INLINE char *__strpbrk_c2 (__const char *__s, int __accept1,
				     int __accept2);
__STRING_INLINE char *
__strpbrk_c2 (__const char *__s, int __accept1, int __accept2)
{
  /* Please note that __accept1 and __accept2 never can be '\0'.  */
  while (*__s != '\0' && *__s != __accept1 && *__s != __accept2)
    ++__s;
  return *__s == '\0' ? NULL : (char *) (size_t) __s;
}

__STRING_INLINE char *__strpbrk_c3 (__const char *__s, int __accept1,
				     int __accept2, int __accept3);
__STRING_INLINE char *
__strpbrk_c3 (__const char *__s, int __accept1, int __accept2,
	      int __accept3)
{
  /* Please note that __accept1 to __accept3 never can be '\0'.  */
  while (*__s != '\0' && *__s != __accept1 && *__s != __accept2
	 && *__s != __accept3)
    ++__s;
  return *__s == '\0' ? NULL : (char *) (size_t) __s;
}
#endif


/* Find the first occurrence of NEEDLE in HAYSTACK.  Newer gcc versions
   do this itself.  */
#if !defined _HAVE_STRING_ARCH_strstr && !__GNUC_PREREQ (2, 97)
# define strstr(haystack, needle) \
  (__extension__ (__builtin_constant_p (needle) && __string2_1bptr_p (needle) \
		  ? (((__const char *) (needle))[0] == '\0'		      \
		     ? (char *) (size_t) (haystack)			      \
		     : (((__const char *) (needle))[1] == '\0'		      \
			? strchr (haystack,				      \
				  ((__const char *) (needle))[0]) 	      \
			: strstr (haystack, needle)))			      \
		  : strstr (haystack, needle)))
#endif


#if !defined _HAVE_STRING_ARCH_strtok_r || defined _FORCE_INLINES
# ifndef _HAVE_STRING_ARCH_strtok_r
#  define __strtok_r(s, sep, nextp) \
  (__extension__ (__builtin_constant_p (sep) && __string2_1bptr_p (sep)	      \
		  && ((__const char *) (sep))[0] != '\0'		      \
		  && ((__const char *) (sep))[1] == '\0'		      \
		  ? __strtok_r_1c (s, ((__const char *) (sep))[0], nextp)     \
		  : __strtok_r (s, sep, nextp)))
# endif

__STRING_INLINE char *__strtok_r_1c (char *__s, char __sep, char **__nextp);
__STRING_INLINE char *
__strtok_r_1c (char *__s, char __sep, char **__nextp)
{
  char *__result;
  if (__s == NULL)
    __s = *__nextp;
  while (*__s == __sep)
    ++__s;
  __result = NULL;
  if (*__s != '\0')
    {
      __result = __s++;
      while (*__s != '\0')
	if (*__s++ == __sep)
	  {
	    __s[-1] = '\0';
	    break;
	  }
    }
  *__nextp = __s;
  return __result;
}
# if defined __USE_POSIX || defined __USE_MISC
#  define strtok_r(s, sep, nextp) __strtok_r (s, sep, nextp)
# endif
#endif


#if !defined _HAVE_STRING_ARCH_strsep || defined _FORCE_INLINES
# ifndef _HAVE_STRING_ARCH_strsep

extern char *__strsep_g (char **__stringp, __const char *__delim);
#  define __strsep(s, reject) \
  __extension__								      \
  ({ char __r0, __r1, __r2;						      \
     (__builtin_constant_p (reject) && __string2_1bptr_p (reject)	      \
      && (__r0 = ((__const char *) (reject))[0],			      \
	  ((__const char *) (reject))[0] != '\0')			      \
      ? ((__r1 = ((__const char *) (reject))[1],			      \
	 ((__const char *) (reject))[1] == '\0')			      \
	 ? __strsep_1c (s, __r0)					      \
	 : ((__r2 = ((__const char *) (reject))[2], __r2 == '\0')	      \
	    ? __strsep_2c (s, __r0, __r1)				      \
	    : (((__const char *) (reject))[3] == '\0'			      \
	       ? __strsep_3c (s, __r0, __r1, __r2)			      \
	       : __strsep_g (s, reject))))				      \
      : __strsep_g (s, reject)); })
# endif

__STRING_INLINE char *__strsep_1c (char **__s, char __reject);
__STRING_INLINE char *
__strsep_1c (char **__s, char __reject)
{
  register char *__retval = *__s;
  if (__retval != NULL && (*__s = strchr (__retval, __reject)) != NULL)
    *(*__s)++ = '\0';
  return __retval;
}

__STRING_INLINE char *__strsep_2c (char **__s, char __reject1, char __reject2);
__STRING_INLINE char *
__strsep_2c (char **__s, char __reject1, char __reject2)
{
  register char *__retval = *__s;
  if (__retval != NULL)
    {
      register char *__cp = __retval;
      while (1)
	{
	  if (*__cp == '\0')
	    {
	      __cp = NULL;
	  break;
	    }
	  if (*__cp == __reject1 || *__cp == __reject2)
	    {
	      *__cp++ = '\0';
	      break;
	    }
	  ++__cp;
	}
      *__s = __cp;
    }
  return __retval;
}

__STRING_INLINE char *__strsep_3c (char **__s, char __reject1, char __reject2,
				   char __reject3);
__STRING_INLINE char *
__strsep_3c (char **__s, char __reject1, char __reject2, char __reject3)
{
  register char *__retval = *__s;
  if (__retval != NULL)
    {
      register char *__cp = __retval;
      while (1)
	{
	  if (*__cp == '\0')
	    {
	      __cp = NULL;
	  break;
	    }
	  if (*__cp == __reject1 || *__cp == __reject2 || *__cp == __reject3)
	    {
	      *__cp++ = '\0';
	      break;
	    }
	  ++__cp;
	}
      *__s = __cp;
    }
  return __retval;
}
# ifdef __USE_BSD
#  define strsep(s, reject) __strsep (s, reject)
# endif
#endif

/* We need the memory allocation functions for inline strdup().
   Referring to stdlib.h (even minimally) is not allowed
   in any of the tight standards compliant modes.  */
#ifdef __USE_MISC

# if !defined _HAVE_STRING_ARCH_strdup || !defined _HAVE_STRING_ARCH_strndup
#  define __need_malloc_and_calloc
#  include <stdlib.h>
# endif

# ifndef _HAVE_STRING_ARCH_strdup

extern char *__strdup (__const char *__string) __THROW __attribute_malloc__;
#  define __strdup(s) \
  (__extension__ (__builtin_constant_p (s) && __string2_1bptr_p (s)	      \
		  ? (((__const char *) (s))[0] == '\0'			      \
		     ? (char *) calloc ((size_t) 1, (size_t) 1)		      \
		     : ({ size_t __len = strlen (s) + 1;		      \
			  char *__retval = (char *) malloc (__len);	      \
			  if (__retval != NULL)				      \
			    __retval = (char *) memcpy (__retval, s, __len);  \
			  __retval; }))					      \
		  : __strdup (s)))

#  if defined __USE_SVID || defined __USE_BSD || defined __USE_XOPEN_EXTENDED
#   define strdup(s) __strdup (s)
#  endif
# endif

# ifndef _HAVE_STRING_ARCH_strndup

extern char *__strndup (__const char *__string, size_t __n)
     __THROW __attribute_malloc__;
#  define __strndup(s, n) \
  (__extension__ (__builtin_constant_p (s) && __string2_1bptr_p (s)	      \
		  ? (((__const char *) (s))[0] == '\0'			      \
		     ? (char *) calloc ((size_t) 1, (size_t) 1)		      \
		     : ({ size_t __len = strlen (s) + 1;		      \
			  size_t __n = (n);				      \
			  char *__retval;				      \
			  if (__n < __len)				      \
			    __len = __n + 1;				      \
			  __retval = (char *) malloc (__len);		      \
			  if (__retval != NULL)				      \
			    {						      \
			      __retval[__len - 1] = '\0';		      \
			      __retval = (char *) memcpy (__retval, s,	      \
							  __len - 1);	      \
			    }						      \
			  __retval; }))					      \
		  : __strndup (s, n)))

#  ifdef __USE_GNU
#   define strndup(s, n) __strndup (s, n)
#  endif
# endif

#endif /* Use misc. or use GNU.  */

#ifndef _FORCE_INLINES
# undef __STRING_INLINE
#endif

#endif /* No string inlines.  */
