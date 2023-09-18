/* O_*, F_*, FD_* bit values for Linux/x86-64.
   Copyright (C) 2001,2002,2004,2006,2007,2009,2010
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

#ifndef	_FCNTL_H
# error "Never use <bits/fcntl.h> directly; include <fcntl.h> instead."
#endif

#include <sys/types.h>
#include <bits/wordsize.h>
#ifdef __USE_GNU
# include <bits/uio.h>
#endif


/* open/fcntl - O_SYNC is only implemented on blocks devices and on files
   located on a few file systems.  */
#define O_ACCMODE	   0003
#define O_RDONLY	     00
#define O_WRONLY	     01
#define O_RDWR		     02
#define O_CREAT		   0100	/* not fcntl */
#define O_EXCL		   0200	/* not fcntl */
#define O_NOCTTY	   0400	/* not fcntl */
#define O_TRUNC		  01000	/* not fcntl */
#define O_APPEND	  02000
#define O_NONBLOCK	  04000
#define O_NDELAY	O_NONBLOCK
#define O_SYNC	       04010000
#define O_FSYNC		 O_SYNC
#define O_ASYNC		 020000

#ifdef __USE_XOPEN2K8
# define O_DIRECTORY	0200000	/* Must be a directory.	 */
# define O_NOFOLLOW	0400000	/* Do not follow links.	 */
# define O_CLOEXEC     02000000 /* Set close_on_exec.  */
#endif
#ifdef __USE_GNU
# define O_DIRECT	 040000	/* Direct disk access.	*/
# define O_NOATIME     01000000 /* Do not set atime.  */
#endif

/* For now Linux has synchronisity options for data and read operations.
   We define the symbols here but let them do the same as O_SYNC since
   this is a superset.	*/
#if defined __USE_POSIX199309 || defined __USE_UNIX98
# define O_DSYNC	010000	/* Synchronize data.  */
# define O_RSYNC	O_SYNC	/* Synchronize read operations.	 */
#endif

#ifdef __USE_LARGEFILE64
# if __WORDSIZE == 64
#  define O_LARGEFILE	0
# else
#  define O_LARGEFILE	0100000
# endif
#endif

/* Values for the second argument to `fcntl'.  */
#define F_DUPFD		0	/* Duplicate file descriptor.  */
#define F_GETFD		1	/* Get file descriptor flags.  */
#define F_SETFD		2	/* Set file descriptor flags.  */
#define F_GETFL		3	/* Get file status flags.  */
#define F_SETFL		4	/* Set file status flags.  */
#if __WORDSIZE == 64
# define F_GETLK	5	/* Get record locking info.  */
# define F_SETLK	6	/* Set record locking info (non-blocking).  */
# define F_SETLKW	7	/* Set record locking info (blocking).	*/
/* Not necessary, we always have 64-bit offsets.  */
# define F_GETLK64	5	/* Get record locking info.  */
# define F_SETLK64	6	/* Set record locking info (non-blocking).  */
# define F_SETLKW64	7	/* Set record locking info (blocking).	*/
#else
# ifndef __USE_FILE_OFFSET64
#  define F_GETLK	5	/* Get record locking info.  */
#  define F_SETLK	6	/* Set record locking info (non-blocking).  */
#  define F_SETLKW	7	/* Set record locking info (blocking).	*/
# else
#  define F_GETLK	F_GETLK64  /* Get record locking info.	*/
#  define F_SETLK	F_SETLK64  /* Set record locking info (non-blocking).*/
#  define F_SETLKW	F_SETLKW64 /* Set record locking info (blocking).  */
# endif
# define F_GETLK64	12	/* Get record locking info.  */
# define F_SETLK64	13	/* Set record locking info (non-blocking).  */
# define F_SETLKW64	14	/* Set record locking info (blocking).	*/
#endif

#if defined __USE_BSD || defined __USE_UNIX98 || defined __USE_XOPEN2K8
# define F_SETOWN	8	/* Get owner (process receiving SIGIO).  */
# define F_GETOWN	9	/* Set owner (process receiving SIGIO).  */
#endif

#ifdef __USE_GNU
# define F_SETSIG	10	/* Set number of signal to be sent.  */
# define F_GETSIG	11	/* Get number of signal to be sent.  */
# define F_SETOWN_EX	15	/* Get owner (thread receiving SIGIO).  */
# define F_GETOWN_EX	16	/* Set owner (thread receiving SIGIO).  */
#endif

#ifdef __USE_GNU
# define F_SETLEASE	1024	/* Set a lease.	 */
# define F_GETLEASE	1025	/* Enquire what lease is active.  */
# define F_NOTIFY	1026	/* Request notfications on a directory.	 */
#endif
#ifdef __USE_XOPEN2K8
# define F_DUPFD_CLOEXEC 1030	/* Duplicate file descriptor with
				   close-on-exit set.  */
#endif

/* For F_[GET|SET]FD.  */
#define FD_CLOEXEC	1	/* actually anything with low bit set goes */

/* For posix fcntl() and `l_type' field of a `struct flock' for lockf().  */
#define F_RDLCK		0	/* Read lock.  */
#define F_WRLCK		1	/* Write lock.	*/
#define F_UNLCK		2	/* Remove lock.	 */

/* For old implementation of bsd flock().  */
#define F_EXLCK		4	/* or 3 */
#define F_SHLCK		8	/* or 4 */

#ifdef __USE_BSD
/* Operations for bsd flock(), also used by the kernel implementation.	*/
# define LOCK_SH	1	/* shared lock */
# define LOCK_EX	2	/* exclusive lock */
# define LOCK_NB	4	/* or'd with one of the above to prevent
				   blocking */
# define LOCK_UN	8	/* remove lock */
#endif

#ifdef __USE_GNU
# define LOCK_MAND	32	/* This is a mandatory flock:	*/
# define LOCK_READ	64	/* ... which allows concurrent read operations.	 */
# define LOCK_WRITE	128	/* ... which allows concurrent write operations.  */
# define LOCK_RW	192	/* ... Which allows concurrent read & write operations.	 */
#endif

#ifdef __USE_GNU
/* Types of directory notifications that may be requested with F_NOTIFY.  */
# define DN_ACCESS	0x00000001	/* File accessed.  */
# define DN_MODIFY	0x00000002	/* File modified.  */
# define DN_CREATE	0x00000004	/* File created.  */
# define DN_DELETE	0x00000008	/* File removed.  */
# define DN_RENAME	0x00000010	/* File renamed.  */
# define DN_ATTRIB	0x00000020	/* File changed attibutes.  */
# define DN_MULTISHOT	0x80000000	/* Don't remove notifier.  */
#endif

struct flock
  {
    short int l_type;	/* Type of lock: F_RDLCK, F_WRLCK, or F_UNLCK.	*/
    short int l_whence;	/* Where `l_start' is relative to (like `lseek').  */
#ifndef __USE_FILE_OFFSET64
    __off_t l_start;	/* Offset where the lock begins.  */
    __off_t l_len;	/* Size of the locked area; zero means until EOF.  */
#else
    __off64_t l_start;	/* Offset where the lock begins.  */
    __off64_t l_len;	/* Size of the locked area; zero means until EOF.  */
#endif
    __pid_t l_pid;	/* Process holding the lock.  */
  };

#ifdef __USE_LARGEFILE64
struct flock64
  {
    short int l_type;	/* Type of lock: F_RDLCK, F_WRLCK, or F_UNLCK.	*/
    short int l_whence;	/* Where `l_start' is relative to (like `lseek').  */
    __off64_t l_start;	/* Offset where the lock begins.  */
    __off64_t l_len;	/* Size of the locked area; zero means until EOF.  */
    __pid_t l_pid;	/* Process holding the lock.  */
  };
#endif

#ifdef __USE_GNU
/* Owner types.  */
enum __pid_type
  {
    F_OWNER_TID = 0,		/* Kernel thread.  */
    F_OWNER_PID,		/* Process.  */
    F_OWNER_PGRP,		/* Process group.  */
    F_OWNER_GID = F_OWNER_PGRP	/* Alternative, obsolete name.  */
  };

/* Structure to use with F_GETOWN_EX and F_SETOWN_EX.  */
struct f_owner_ex
  {
    enum __pid_type type;	/* Owner type of ID.  */
    __pid_t pid;		/* ID of owner.  */
  };
#endif

/* Define some more compatibility macros to be backward compatible with
   BSD systems which did not managed to hide these kernel macros.  */
#ifdef	__USE_BSD
# define FAPPEND	O_APPEND
# define FFSYNC		O_FSYNC
# define FASYNC		O_ASYNC
# define FNONBLOCK	O_NONBLOCK
# define FNDELAY	O_NDELAY
#endif /* Use BSD.  */

/* Advise to `posix_fadvise'.  */
#ifdef __USE_XOPEN2K
# define POSIX_FADV_NORMAL	0 /* No further special treatment.  */
# define POSIX_FADV_RANDOM	1 /* Expect random page references.  */
# define POSIX_FADV_SEQUENTIAL	2 /* Expect sequential page references.	 */
# define POSIX_FADV_WILLNEED	3 /* Will need these pages.  */
# define POSIX_FADV_DONTNEED	4 /* Don't need these pages.  */
# define POSIX_FADV_NOREUSE	5 /* Data will be accessed once.  */
#endif


#ifdef __USE_GNU
/* Flags for SYNC_FILE_RANGE.  */
# define SYNC_FILE_RANGE_WAIT_BEFORE	1 /* Wait upon writeout of all pages
					     in the range before performing the
					     write.  */
# define SYNC_FILE_RANGE_WRITE		2 /* Initiate writeout of all those
					     dirty pages in the range which are
					     not presently under writeback.  */
# define SYNC_FILE_RANGE_WAIT_AFTER	4 /* Wait upon writeout of all pages in
					     the range after performing the
					     write.  */

/* Flags for SPLICE and VMSPLICE.  */
# define SPLICE_F_MOVE		1	/* Move pages instead of copying.  */
# define SPLICE_F_NONBLOCK	2	/* Don't block on the pipe splicing
					   (but we may still block on the fd
					   we splice from/to).  */
# define SPLICE_F_MORE		4	/* Expect more data.  */
# define SPLICE_F_GIFT		8	/* Pages passed in are a gift.  */
#endif

__BEGIN_DECLS

#ifdef __USE_GNU

/* Provide kernel hint to read ahead.  */
extern ssize_t readahead (int __fd, __off64_t __offset, size_t __count)
    __THROW;


/* Selective file content synch'ing.  */
extern int sync_file_range (int __fd, __off64_t __offset, __off64_t __count,
			    unsigned int __flags);


/* Splice address range into a pipe.  */
extern ssize_t vmsplice (int __fdout, const struct iovec *__iov,
			 size_t __count, unsigned int __flags);

/* Splice two files together.  */
extern ssize_t splice (int __fdin, __off64_t *__offin, int __fdout,
		       __off64_t *__offout, size_t __len,
		       unsigned int __flags);

/* In-kernel implementation of tee for pipe buffers.  */
extern ssize_t tee (int __fdin, int __fdout, size_t __len,
		    unsigned int __flags);

/* Reserve storage for the data of the file associated with FD.  */
# ifndef __USE_FILE_OFFSET64
extern int fallocate (int __fd, int __mode, __off_t __offset, __off_t __len);
# else
#  ifdef __REDIRECT
extern int __REDIRECT (fallocate, (int __fd, int __mode, __off64_t __offset,
				   __off64_t __len),
		       fallocate64);
#  else
#   define fallocate fallocate64
#  endif
# endif
# ifdef __USE_LARGEFILE64
extern int fallocate64 (int __fd, int __mode, __off64_t __offset,
			__off64_t __len);
# endif
#endif

__END_DECLS
