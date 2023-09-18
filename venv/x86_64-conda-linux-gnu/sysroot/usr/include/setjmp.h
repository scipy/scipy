/* Copyright (C) 1991-1999,2001,2002,2007,2009 Free Software Foundation, Inc.
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
 *	ISO C99 Standard: 7.13 Nonlocal jumps	<setjmp.h>
 */

#ifndef	_SETJMP_H
#define	_SETJMP_H	1

#include <features.h>

__BEGIN_DECLS

#include <bits/setjmp.h>		/* Get `__jmp_buf'.  */
#include <bits/sigset.h>		/* Get `__sigset_t'.  */


/* Calling environment, plus possibly a saved signal mask.  */
struct __jmp_buf_tag
  {
    /* NOTE: The machine-dependent definitions of `__sigsetjmp'
       assume that a `jmp_buf' begins with a `__jmp_buf' and that
       `__mask_was_saved' follows it.  Do not move these members
       or add others before it.  */
    __jmp_buf __jmpbuf;		/* Calling environment.  */
    int __mask_was_saved;	/* Saved the signal mask?  */
    __sigset_t __saved_mask;	/* Saved signal mask.  */
  };


__BEGIN_NAMESPACE_STD

typedef struct __jmp_buf_tag jmp_buf[1];

/* Store the calling environment in ENV, also saving the signal mask.
   Return 0.  */
extern int setjmp (jmp_buf __env) __THROW;

__END_NAMESPACE_STD

/* Store the calling environment in ENV, also saving the
   signal mask if SAVEMASK is nonzero.  Return 0.
   This is the internal name for `sigsetjmp'.  */
extern int __sigsetjmp (struct __jmp_buf_tag __env[1], int __savemask) __THROW;

#ifndef	__FAVOR_BSD
/* Store the calling environment in ENV, not saving the signal mask.
   Return 0.  */
extern int _setjmp (struct __jmp_buf_tag __env[1]) __THROW;

/* Do not save the signal mask.  This is equivalent to the `_setjmp'
   BSD function.  */
# define setjmp(env)	_setjmp (env)
#else
/* We are in 4.3 BSD-compatibility mode in which `setjmp'
   saves the signal mask like `sigsetjmp (ENV, 1)'.  We have to
   define a macro since ISO C says `setjmp' is one.  */
# define setjmp(env)	setjmp (env)
#endif /* Favor BSD.  */


__BEGIN_NAMESPACE_STD

/* Jump to the environment saved in ENV, making the
   `setjmp' call there return VAL, or 1 if VAL is 0.  */
extern void longjmp (struct __jmp_buf_tag __env[1], int __val)
     __THROW __attribute__ ((__noreturn__));

__END_NAMESPACE_STD

#if defined __USE_BSD || defined __USE_XOPEN
/* Same.  Usually `_longjmp' is used with `_setjmp', which does not save
   the signal mask.  But it is how ENV was saved that determines whether
   `longjmp' restores the mask; `_longjmp' is just an alias.  */
extern void _longjmp (struct __jmp_buf_tag __env[1], int __val)
     __THROW __attribute__ ((__noreturn__));
#endif


#ifdef	__USE_POSIX
/* Use the same type for `jmp_buf' and `sigjmp_buf'.
   The `__mask_was_saved' flag determines whether
   or not `longjmp' will restore the signal mask.  */
typedef struct __jmp_buf_tag sigjmp_buf[1];

/* Store the calling environment in ENV, also saving the
   signal mask if SAVEMASK is nonzero.  Return 0.  */
# define sigsetjmp(env, savemask)	__sigsetjmp (env, savemask)

/* Jump to the environment saved in ENV, making the
   sigsetjmp call there return VAL, or 1 if VAL is 0.
   Restore the signal mask if that sigsetjmp call saved it.
   This is just an alias `longjmp'.  */
extern void siglongjmp (sigjmp_buf __env, int __val)
     __THROW __attribute__ ((__noreturn__));
#endif /* Use POSIX.  */


/* Define helper functions to catch unsafe code.  */
#if __USE_FORTIFY_LEVEL > 0
# include <bits/setjmp2.h>
#endif

__END_DECLS

#endif /* setjmp.h  */
