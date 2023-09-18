/* Prototypes and definition for malloc implementation.
   Copyright (C) 1996, 1997, 1999, 2000, 2002-2004, 2005, 2007, 2009
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

#ifndef _MALLOC_H
#define _MALLOC_H 1

#include <features.h>
#include <stddef.h>
#include <stdio.h>
# define __malloc_ptr_t  void *

/* Used by GNU libc internals. */
#define __malloc_size_t size_t
#define __malloc_ptrdiff_t ptrdiff_t

#ifdef __GNUC__

# define __MALLOC_P(args)	args __THROW
/* This macro will be used for functions which might take C++ callback
   functions.  */
# define __MALLOC_PMT(args)	args

#else	/* Not GCC.  */

# define __MALLOC_P(args)	args
# define __MALLOC_PMT(args)	args

#endif	/* GCC.  */


__BEGIN_DECLS

/* Allocate SIZE bytes of memory.  */
extern void *malloc __MALLOC_P ((size_t __size)) __attribute_malloc__ __wur;

/* Allocate NMEMB elements of SIZE bytes each, all initialized to 0.  */
extern void *calloc __MALLOC_P ((size_t __nmemb, size_t __size))
       __attribute_malloc__ __wur;

/* Re-allocate the previously allocated block in __ptr, making the new
   block SIZE bytes long.  */
/* __attribute_malloc__ is not used, because if realloc returns
   the same pointer that was passed to it, aliasing needs to be allowed
   between objects pointed by the old and new pointers.  */
extern void *realloc __MALLOC_P ((void *__ptr, size_t __size))
       __attribute_warn_unused_result__;

/* Free a block allocated by `malloc', `realloc' or `calloc'.  */
extern void free __MALLOC_P ((void *__ptr));

/* Free a block allocated by `calloc'. */
extern void cfree __MALLOC_P ((void *__ptr));

/* Allocate SIZE bytes allocated to ALIGNMENT bytes.  */
extern void *memalign __MALLOC_P ((size_t __alignment, size_t __size))
       __attribute_malloc__ __wur;

/* Allocate SIZE bytes on a page boundary.  */
extern void *valloc __MALLOC_P ((size_t __size))
       __attribute_malloc__ __wur;

/* Equivalent to valloc(minimum-page-that-holds(n)), that is, round up
   __size to nearest pagesize. */
extern void * pvalloc __MALLOC_P ((size_t __size))
       __attribute_malloc__ __wur;

/* Underlying allocation function; successive calls should return
   contiguous pieces of memory.  */
extern void *(*__morecore) __MALLOC_PMT ((ptrdiff_t __size));

/* Default value of `__morecore'.  */
extern void *__default_morecore __MALLOC_P ((ptrdiff_t __size))
       __attribute_malloc__;

/* SVID2/XPG mallinfo structure */

struct mallinfo {
  int arena;    /* non-mmapped space allocated from system */
  int ordblks;  /* number of free chunks */
  int smblks;   /* number of fastbin blocks */
  int hblks;    /* number of mmapped regions */
  int hblkhd;   /* space in mmapped regions */
  int usmblks;  /* maximum total allocated space */
  int fsmblks;  /* space available in freed fastbin blocks */
  int uordblks; /* total allocated space */
  int fordblks; /* total free space */
  int keepcost; /* top-most, releasable (via malloc_trim) space */
};

/* Returns a copy of the updated current mallinfo. */
extern struct mallinfo mallinfo __MALLOC_P ((void));

/* SVID2/XPG mallopt options */
#ifndef M_MXFAST
# define M_MXFAST  1	/* maximum request size for "fastbins" */
#endif
#ifndef M_NLBLKS
# define M_NLBLKS  2	/* UNUSED in this malloc */
#endif
#ifndef M_GRAIN
# define M_GRAIN   3	/* UNUSED in this malloc */
#endif
#ifndef M_KEEP
# define M_KEEP    4	/* UNUSED in this malloc */
#endif

/* mallopt options that actually do something */
#define M_TRIM_THRESHOLD    -1
#define M_TOP_PAD           -2
#define M_MMAP_THRESHOLD    -3
#define M_MMAP_MAX          -4
#define M_CHECK_ACTION      -5
#define M_PERTURB	    -6
#define M_ARENA_TEST	    -7
#define M_ARENA_MAX	    -8

/* General SVID/XPG interface to tunable parameters. */
extern int mallopt __MALLOC_P ((int __param, int __val));

/* Release all but __pad bytes of freed top-most memory back to the
   system. Return 1 if successful, else 0. */
extern int malloc_trim __MALLOC_P ((size_t __pad));

/* Report the number of usable allocated bytes associated with allocated
   chunk __ptr. */
extern size_t malloc_usable_size __MALLOC_P ((void *__ptr));

/* Prints brief summary statistics on stderr. */
extern void malloc_stats __MALLOC_P ((void));

/* Output information about state of allocator to stream FP.  */
extern int malloc_info (int __options, FILE *__fp);

/* Record the state of all malloc variables in an opaque data structure. */
extern void *malloc_get_state __MALLOC_P ((void));

/* Restore the state of all malloc variables from data obtained with
   malloc_get_state(). */
extern int malloc_set_state __MALLOC_P ((void *__ptr));

/* Called once when malloc is initialized; redefining this variable in
   the application provides the preferred way to set up the hook
   pointers. */
extern void (*__malloc_initialize_hook) __MALLOC_PMT ((void));
/* Hooks for debugging and user-defined versions. */
extern void (*__free_hook) __MALLOC_PMT ((void *__ptr,
					__const __malloc_ptr_t));
extern void *(*__malloc_hook) __MALLOC_PMT ((size_t __size,
					     __const __malloc_ptr_t));
extern void *(*__realloc_hook) __MALLOC_PMT ((void *__ptr, size_t __size,
					      __const __malloc_ptr_t));
extern void *(*__memalign_hook) __MALLOC_PMT ((size_t __alignment,
					       size_t __size,
					       __const __malloc_ptr_t));
extern void (*__after_morecore_hook) __MALLOC_PMT ((void));

/* Activate a standard set of debugging hooks. */
extern void __malloc_check_init __MALLOC_P ((void));


__END_DECLS

#endif /* malloc.h */
