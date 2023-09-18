/* This just represents the non-kernel parts of <linux/quota.h>.
 *
 * here's the corresponding copyright:
 * Copyright (c) 1982, 1986 Regents of the University of California.
 * All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Robert Elz at The University of Melbourne.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef _SYS_QUOTA_H
#define _SYS_QUOTA_H 1

#include <features.h>
#include <sys/types.h>

/*
 * Select between different incompatible quota versions.
 * Default to the version used by Linux kernel version 2.4.22
 * or later.  */
#ifndef _LINUX_QUOTA_VERSION
# define _LINUX_QUOTA_VERSION 2
#endif

/*
 * Convert diskblocks to blocks and the other way around.
 * currently only to fool the BSD source. :-)
 */
#define dbtob(num) ((num) << 10)
#define btodb(num) ((num) >> 10)

/*
 * Convert count of filesystem blocks to diskquota blocks, meant
 * for filesystems where i_blksize != BLOCK_SIZE
 */
#define fs_to_dq_blocks(num, blksize) (((num) * (blksize)) / BLOCK_SIZE)

/*
 * Definitions for disk quotas imposed on the average user
 * (big brother finally hits Linux).
 *
 * The following constants define the amount of time given a user
 * before the soft limits are treated as hard limits (usually resulting
 * in an allocation failure). The timer is started when the user crosses
 * their soft limit, it is reset when they go below their soft limit.
 */
#define MAX_IQ_TIME  604800	/* (7*24*60*60) 1 week */
#define MAX_DQ_TIME  604800	/* (7*24*60*60) 1 week */

#define MAXQUOTAS 2
#define USRQUOTA  0		/* element used for user quotas */
#define GRPQUOTA  1		/* element used for group quotas */

/*
 * Definitions for the default names of the quotas files.
 */
#define INITQFNAMES { \
   "user",      /* USRQUOTA */ \
   "group",   /* GRPQUOTA */ \
   "undefined", \
};

#define QUOTAFILENAME "quota"
#define QUOTAGROUP "staff"

#define NR_DQHASH 43          /* Just an arbitrary number any suggestions ? */
#define NR_DQUOTS 256         /* Number of quotas active at one time */

/*
 * Command definitions for the 'quotactl' system call.
 * The commands are broken into a main command defined below
 * and a subcommand that is used to convey the type of
 * quota that is being manipulated (see above).
 */
#define SUBCMDMASK  0x00ff
#define SUBCMDSHIFT 8
#define QCMD(cmd, type)  (((cmd) << SUBCMDSHIFT) | ((type) & SUBCMDMASK))

#if _LINUX_QUOTA_VERSION < 2
# define Q_QUOTAON  0x0100	/* enable quotas */
# define Q_QUOTAOFF 0x0200	/* disable quotas */
# define Q_GETQUOTA 0x0300	/* get limits and usage */
# define Q_SETQUOTA 0x0400	/* set limits and usage */
# define Q_SETUSE   0x0500	/* set usage */
# define Q_SYNC     0x0600	/* sync disk copy of a filesystems quotas */
# define Q_SETQLIM  0x0700	/* set limits */
# define Q_GETSTATS 0x0800	/* get collected stats */
# define Q_RSQUASH  0x1000	/* set root_squash option */
#else
# define Q_SYNC     0x800001	/* sync disk copy of a filesystems quotas */
# define Q_QUOTAON  0x800002	/* turn quotas on */
# define Q_QUOTAOFF 0x800003	/* turn quotas off */
# define Q_GETFMT   0x800004	/* get quota format used on given filesystem */
# define Q_GETINFO  0x800005	/* get information about quota files */
# define Q_SETINFO  0x800006	/* set information about quota files */
# define Q_GETQUOTA 0x800007	/* get user quota structure */
# define Q_SETQUOTA 0x800008	/* set user quota structure */
#endif

/*
 * The following structure defines the format of the disk quota file
 * (as it appears on disk) - the file is an array of these structures
 * indexed by user or group number.
 */
#if _LINUX_QUOTA_VERSION < 2
struct dqblk
  {
    u_int32_t dqb_bhardlimit;	/* absolute limit on disk blks alloc */
    u_int32_t dqb_bsoftlimit;	/* preferred limit on disk blks */
    u_int32_t dqb_curblocks;	/* current block count */
    u_int32_t dqb_ihardlimit;	/* maximum # allocated inodes */
    u_int32_t dqb_isoftlimit;	/* preferred inode limit */
    u_int32_t dqb_curinodes;	/* current # allocated inodes */
    time_t dqb_btime;		/* time limit for excessive disk use */
    time_t dqb_itime;		/* time limit for excessive files */
  };
#else

/* Flags that indicate which fields in dqblk structure are valid.  */
#define QIF_BLIMITS	1
#define QIF_SPACE	2
#define QIF_ILIMITS	4
#define QIF_INODES	8
#define QIF_BTIME	16
#define QIF_ITIME	32
#define QIF_LIMITS	(QIF_BLIMITS | QIF_ILIMITS)
#define QIF_USAGE	(QIF_SPACE | QIF_INODES)
#define QIF_TIMES	(QIF_BTIME | QIF_ITIME)
#define QIF_ALL		(QIF_LIMITS | QIF_USAGE | QIF_TIMES)

struct dqblk
  {
    u_int64_t dqb_bhardlimit;	/* absolute limit on disk quota blocks alloc */
    u_int64_t dqb_bsoftlimit;	/* preferred limit on disk quota blocks */
    u_int64_t dqb_curspace;	/* current quota block count */
    u_int64_t dqb_ihardlimit;	/* maximum # allocated inodes */
    u_int64_t dqb_isoftlimit;	/* preferred inode limit */
    u_int64_t dqb_curinodes;	/* current # allocated inodes */
    u_int64_t dqb_btime;	/* time limit for excessive disk use */
    u_int64_t dqb_itime;	/* time limit for excessive files */
    u_int32_t dqb_valid;	/* bitmask of QIF_* constants */
  };
#endif

/*
 * Shorthand notation.
 */
#define	dq_bhardlimit	dq_dqb.dqb_bhardlimit
#define	dq_bsoftlimit	dq_dqb.dqb_bsoftlimit
#if _LINUX_QUOTA_VERSION < 2
# define dq_curblocks	dq_dqb.dqb_curblocks
#else
# define dq_curspace	dq_dqb.dqb_curspace
# define dq_valid	dq_dqb.dqb_valid
#endif
#define	dq_ihardlimit	dq_dqb.dqb_ihardlimit
#define	dq_isoftlimit	dq_dqb.dqb_isoftlimit
#define	dq_curinodes	dq_dqb.dqb_curinodes
#define	dq_btime	dq_dqb.dqb_btime
#define	dq_itime	dq_dqb.dqb_itime

#define dqoff(UID)      ((loff_t)((UID) * sizeof (struct dqblk)))

#if _LINUX_QUOTA_VERSION < 2
struct dqstats
  {
    u_int32_t lookups;
    u_int32_t drops;
    u_int32_t reads;
    u_int32_t writes;
    u_int32_t cache_hits;
    u_int32_t pages_allocated;
    u_int32_t allocated_dquots;
    u_int32_t free_dquots;
    u_int32_t syncs;
  };
#else

/* Flags that indicate which fields in dqinfo structure are valid.  */
# define IIF_BGRACE	1
# define IIF_IGRACE	2
# define IIF_FLAGS	4
# define IIF_ALL	(IIF_BGRACE | IIF_IGRACE | IIF_FLAGS)

struct dqinfo
  {
    u_int64_t dqi_bgrace;
    u_int64_t dqi_igrace;
    u_int32_t dqi_flags;
    u_int32_t dqi_valid;
  };
#endif

__BEGIN_DECLS

extern int quotactl (int __cmd, const char *__special, int __id,
		     caddr_t __addr) __THROW;

__END_DECLS

#endif /* sys/quota.h */
