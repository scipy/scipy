/* __sig_atomic_t, __sigset_t, and related definitions.  Linux version.
   Copyright (C) 1991, 1992, 1994, 1996, 1997, 2007
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

#ifndef	_SIGSET_H_types
# define _SIGSET_H_types	1

typedef int __sig_atomic_t;

/* A `sigset_t' has a bit for each signal.  */

# define _SIGSET_NWORDS	(1024 / (8 * sizeof (unsigned long int)))
typedef struct
  {
    unsigned long int __val[_SIGSET_NWORDS];
  } __sigset_t;

#endif


/* We only want to define these functions if <signal.h> was actually
   included; otherwise we were included just to define the types.  Since we
   are namespace-clean, it wouldn't hurt to define extra macros.  But
   trouble can be caused by functions being defined (e.g., any global
   register vars declared later will cause compilation errors).  */

#if !defined _SIGSET_H_fns && defined _SIGNAL_H
# define _SIGSET_H_fns 1

# ifndef _EXTERN_INLINE
#  define _EXTERN_INLINE __extern_inline
# endif

/* Return a mask that includes the bit for SIG only.  */
# define __sigmask(sig) \
  (((unsigned long int) 1) << (((sig) - 1) % (8 * sizeof (unsigned long int))))

/* Return the word index for SIG.  */
# define __sigword(sig)	(((sig) - 1) / (8 * sizeof (unsigned long int)))

# if defined __GNUC__ && __GNUC__ >= 2
#  define __sigemptyset(set) \
  (__extension__ ({ int __cnt = _SIGSET_NWORDS;				      \
		    sigset_t *__set = (set);				      \
		    while (--__cnt >= 0) __set->__val[__cnt] = 0;	      \
		    0; }))
#  define __sigfillset(set) \
  (__extension__ ({ int __cnt = _SIGSET_NWORDS;				      \
		    sigset_t *__set = (set);				      \
		    while (--__cnt >= 0) __set->__val[__cnt] = ~0UL;	      \
		    0; }))

#  ifdef __USE_GNU
/* The POSIX does not specify for handling the whole signal set in one
   command.  This is often wanted and so we define three more functions
   here.  */
#   define __sigisemptyset(set) \
  (__extension__ ({ int __cnt = _SIGSET_NWORDS;				      \
		    const sigset_t *__set = (set);			      \
		    int __ret = __set->__val[--__cnt];			      \
		    while (!__ret && --__cnt >= 0)			      \
			__ret = __set->__val[__cnt];			      \
		    __ret == 0; }))
#   define __sigandset(dest, left, right) \
  (__extension__ ({ int __cnt = _SIGSET_NWORDS;				      \
		    sigset_t *__dest = (dest);				      \
		    const sigset_t *__left = (left);			      \
		    const sigset_t *__right = (right);			      \
		    while (--__cnt >= 0)				      \
		      __dest->__val[__cnt] = (__left->__val[__cnt]	      \
					      & __right->__val[__cnt]);	      \
		    0; }))
#   define __sigorset(dest, left, right) \
  (__extension__ ({ int __cnt = _SIGSET_NWORDS;				      \
		    sigset_t *__dest = (dest);				      \
		    const sigset_t *__left = (left);			      \
		    const sigset_t *__right = (right);			      \
		    while (--__cnt >= 0)				      \
		      __dest->__val[__cnt] = (__left->__val[__cnt]	      \
					      | __right->__val[__cnt]);	      \
		    0; }))
#  endif
# endif

/* These functions needn't check for a bogus signal number -- error
   checking is done in the non __ versions.  */

extern int __sigismember (__const __sigset_t *, int);
extern int __sigaddset (__sigset_t *, int);
extern int __sigdelset (__sigset_t *, int);

# ifdef __USE_EXTERN_INLINES
#  define __SIGSETFN(NAME, BODY, CONST)					      \
  _EXTERN_INLINE int							      \
  NAME (CONST __sigset_t *__set, int __sig)				      \
  {									      \
    unsigned long int __mask = __sigmask (__sig);			      \
    unsigned long int __word = __sigword (__sig);			      \
    return BODY;							      \
  }

__SIGSETFN (__sigismember, (__set->__val[__word] & __mask) ? 1 : 0, __const)
__SIGSETFN (__sigaddset, ((__set->__val[__word] |= __mask), 0), )
__SIGSETFN (__sigdelset, ((__set->__val[__word] &= ~__mask), 0), )

#  undef __SIGSETFN
# endif


#endif /* ! _SIGSET_H_fns.  */
