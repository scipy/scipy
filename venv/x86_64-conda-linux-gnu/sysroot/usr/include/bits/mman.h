/* Definitions for POSIX memory map interface.  Linux/x86_64 version.
   Copyright (C) 2001-2012 Free Software Foundation, Inc.
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

#ifndef _SYS_MMAN_H
# error "Never use <bits/mman.h> directly; include <sys/mman.h> instead."
#endif

/* The following definitions basically come from the kernel headers.
   But the kernel header is not namespace clean.  */


/* Protections are chosen from these bits, OR'd together.  The
   implementation does not necessarily support PROT_EXEC or PROT_WRITE
   without PROT_READ.  The only guarantees are that no writing will be
   allowed without PROT_WRITE and no access will be allowed for PROT_NONE. */

#define PROT_READ	0x1		/* Page can be read.  */
#define PROT_WRITE	0x2		/* Page can be written.  */
#define PROT_EXEC	0x4		/* Page can be executed.  */
#define PROT_NONE	0x0		/* Page can not be accessed.  */
#define PROT_GROWSDOWN	0x01000000	/* Extend change to start of
					   growsdown vma (mprotect only).  */
#define PROT_GROWSUP	0x02000000	/* Extend change to start of
					   growsup vma (mprotect only).  */

/* Sharing types (must choose one and only one of these).  */
#define MAP_SHARED	0x01		/* Share changes.  */
#define MAP_PRIVATE	0x02		/* Changes are private.  */
#ifdef __USE_MISC
# define MAP_TYPE	0x0f		/* Mask for type of mapping.  */
#endif

/* Other flags.  */
#define MAP_FIXED	0x10		/* Interpret addr exactly.  */
#ifdef __USE_MISC
# define MAP_FILE	0
# define MAP_ANONYMOUS	0x20		/* Don't use a file.  */
# define MAP_ANON	MAP_ANONYMOUS
# define MAP_32BIT	0x40		/* Only give out 32-bit addresses.  */
#endif

/* These are Linux-specific.  */
#ifdef __USE_MISC
# define MAP_GROWSDOWN	0x00100		/* Stack-like segment.  */
# define MAP_DENYWRITE	0x00800		/* ETXTBSY */
# define MAP_EXECUTABLE	0x01000		/* Mark it as an executable.  */
# define MAP_LOCKED	0x02000		/* Lock the mapping.  */
# define MAP_NORESERVE	0x04000		/* Don't check for reservations.  */
# define MAP_POPULATE	0x08000		/* Populate (prefault) pagetables.  */
# define MAP_NONBLOCK	0x10000		/* Do not block on IO.  */
# define MAP_STACK	0x20000		/* Allocation is for a stack.  */
# define MAP_HUGETLB	0x40000		/* Create huge page mapping.  */
#endif

/* Flags to `msync'.  */
#define MS_ASYNC	1		/* Sync memory asynchronously.  */
#define MS_SYNC		4		/* Synchronous memory sync.  */
#define MS_INVALIDATE	2		/* Invalidate the caches.  */

/* Flags for `mlockall'.  */
#define MCL_CURRENT	1		/* Lock all currently mapped pages.  */
#define MCL_FUTURE	2		/* Lock all additions to address
					   space.  */

/* Flags for `mremap'.  */
#ifdef __USE_GNU
# define MREMAP_MAYMOVE	1
# define MREMAP_FIXED	2
#endif

/* Advice to `madvise'.  */
#ifdef __USE_BSD
# define MADV_NORMAL	  0	/* No further special treatment.  */
# define MADV_RANDOM	  1	/* Expect random page references.  */
# define MADV_SEQUENTIAL  2	/* Expect sequential page references.  */
# define MADV_WILLNEED	  3	/* Will need these pages.  */
# define MADV_DONTNEED	  4	/* Don't need these pages.  */
# define MADV_REMOVE	  9	/* Remove these pages and resources.  */
# define MADV_DONTFORK	  10	/* Do not inherit across fork.  */
# define MADV_DOFORK	  11	/* Do inherit across fork.  */
# define MADV_MERGEABLE	  12	/* KSM may merge identical pages.  */
# define MADV_UNMERGEABLE 13	/* KSM may not merge identical pages.  */
# define MADV_DONTDUMP	  16    /* Explicity exclude from the core dump,
                                   overrides the coredump filter bits.  */
# define MADV_DODUMP	  17	/* Clear the MADV_DONTDUMP flag.  */
# define MADV_HWPOISON	  100	/* Poison a page for testing.  */
#endif

/* The POSIX people had to invent similar names for the same things.  */
#ifdef __USE_XOPEN2K
# define POSIX_MADV_NORMAL	0 /* No further special treatment.  */
# define POSIX_MADV_RANDOM	1 /* Expect random page references.  */
# define POSIX_MADV_SEQUENTIAL	2 /* Expect sequential page references.  */
# define POSIX_MADV_WILLNEED	3 /* Will need these pages.  */
# define POSIX_MADV_DONTNEED	4 /* Don't need these pages.  */
#endif
