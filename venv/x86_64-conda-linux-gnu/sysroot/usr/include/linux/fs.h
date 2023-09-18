#ifndef _LINUX_FS_H
#define _LINUX_FS_H

/*
 * This file has definitions for some important file table
 * structures etc.
 */

#include <linux/limits.h>
#include <linux/ioctl.h>
#include <linux/blk_types.h>
#include <linux/types.h>

/*
 * It's silly to have NR_OPEN bigger than NR_FILE, but you can change
 * the file limit at runtime and only root can increase the per-process
 * nr_file rlimit, so it's safe to set up a ridiculously high absolute
 * upper limit on files-per-process.
 *
 * Some programs (notably those using select()) may have to be 
 * recompiled to take full advantage of the new limits..  
 */

/* Fixed constants first: */
#undef NR_OPEN
#define INR_OPEN_CUR 1024	/* Initial setting for nfile rlimits */
#define INR_OPEN_MAX 4096	/* Hard limit for nfile rlimits */

#define BLOCK_SIZE_BITS 10
#define BLOCK_SIZE (1<<BLOCK_SIZE_BITS)

#define SEEK_SET	0	/* seek relative to beginning of file */
#define SEEK_CUR	1	/* seek relative to current file position */
#define SEEK_END	2	/* seek relative to end of file */
#define SEEK_MAX	SEEK_END

struct fstrim_range {
	__u64 start;
	__u64 len;
	__u64 minlen;
};

/* And dynamically-tunable limits and defaults: */
struct files_stat_struct {
	unsigned long nr_files;		/* read only */
	unsigned long nr_free_files;	/* read only */
	unsigned long max_files;		/* tunable */
};

struct inodes_stat_t {
	int nr_inodes;
	int nr_unused;
	int dummy[5];		/* padding for sysctl ABI compatibility */
};


#define NR_FILE  8192	/* this can well be larger on a larger system */

#define MAY_EXEC 1
#define MAY_WRITE 2
#define MAY_READ 4
#define MAY_APPEND 8
#define MAY_ACCESS 16
#define MAY_OPEN 32

/*
 * flags in file.f_mode.  Note that FMODE_READ and FMODE_WRITE must correspond
 * to O_WRONLY and O_RDWR via the strange trick in __dentry_open()
 */

/* file is open for reading */
#define FMODE_READ		((fmode_t)1)
/* file is open for writing */
#define FMODE_WRITE		((fmode_t)2)
/* file is seekable */
#define FMODE_LSEEK		((fmode_t)4)
/* file can be accessed using pread */
#define FMODE_PREAD		((fmode_t)8)
/* file can be accessed using pwrite */
#define FMODE_PWRITE		((fmode_t)16)
/* File is opened for execution with sys_execve / sys_uselib */
#define FMODE_EXEC		((fmode_t)32)
/* File is opened with O_NDELAY (only set for block devices) */
#define FMODE_NDELAY		((fmode_t)64)
/* File is opened with O_EXCL (only set for block devices) */
#define FMODE_EXCL		((fmode_t)128)
/* File is opened using open(.., 3, ..) and is writeable only for ioctls
   (specialy hack for floppy.c) */
#define FMODE_WRITE_IOCTL	((fmode_t)256)
/* 32bit hashes as llseek() offset (for directories) */
#define FMODE_32BITHASH         ((fmode_t)0x200)
/* 64bit hashes as llseek() offset (for directories) */
#define FMODE_64BITHASH         ((fmode_t)0x400)

/*
 * Don't update ctime and mtime.
 *
 * Currently a special hack for the XFS open_by_handle ioctl, but we'll
 * hopefully graduate it to a proper O_CMTIME flag supported by open(2) soon.
 */
#define FMODE_NOCMTIME		((fmode_t)2048)

/* Expect random access pattern */
#define FMODE_RANDOM		((fmode_t)4096)

/* File needs atomic accesses to f_pos */
#define FMODE_ATOMIC_POS	((fmode_t)0x8000)


/*
 * The below are the various read and write types that we support. Some of
 * them include behavioral modifiers that send information down to the
 * block layer and IO scheduler. Terminology:
 *
 *	The block layer uses device plugging to defer IO a little bit, in
 *	the hope that we will see more IO very shortly. This increases
 *	coalescing of adjacent IO and thus reduces the number of IOs we
 *	have to send to the device. It also allows for better queuing,
 *	if the IO isn't mergeable. If the caller is going to be waiting
 *	for the IO, then he must ensure that the device is unplugged so
 *	that the IO is dispatched to the driver.
 *
 *	All IO is handled async in Linux. This is fine for background
 *	writes, but for reads or writes that someone waits for completion
 *	on, we want to notify the block layer and IO scheduler so that they
 *	know about it. That allows them to make better scheduling
 *	decisions. So when the below references 'sync' and 'async', it
 *	is referencing this priority hint.
 *
 * With that in mind, the available types are:
 *
 * READ			A normal read operation. Device will be plugged.
 * READ_SYNC		A synchronous read. Device is not plugged, caller can
 *			immediately wait on this read without caring about
 *			unplugging.
 * READA		Used for read-ahead operations. Lower priority, and the
 *			 block layer could (in theory) choose to ignore this
 *			request if it runs into resource problems.
 * WRITE		A normal async write. Device will be plugged.
 * SWRITE		DEPRECATED; use write_dirty_buffer instead.
 *			Like WRITE, but a special case for ll_rw_block() that
 *			tells it to lock the buffer first. Normally a buffer
 *			must be locked before doing IO.
 * WRITE_SYNC_PLUG	Synchronous write. Identical to WRITE, but passes down
 *			the hint that someone will be waiting on this IO
 *			shortly. The device must still be unplugged explicitly,
 *			WRITE_SYNC_PLUG does not do this as we could be
 *			submitting more writes before we actually wait on any
 *			of them.
 * WRITE_SYNC		Like WRITE_SYNC_PLUG, but also unplugs the device
 *			immediately after submission. The write equivalent
 *			of READ_SYNC.
 * WRITE_ODIRECT_PLUG	Special case write for O_DIRECT only.
 * SWRITE_SYNC
 * SWRITE_SYNC_PLUG	DEPRECATED. Like WRITE_SYNC/WRITE_SYNC_PLUG, but locks
 *			the buffer. See SWRITE.
 * WRITE_BARRIER	DEPRECATED. Always fails. Use FLUSH/FUA instead.
 * WRITE_FLUSH		Like WRITE_SYNC but with preceding cache flush.
 * WRITE_FUA		Like WRITE_SYNC but data is guaranteed to be on
 *			non-volatile media on completion.
 * WRITE_FLUSH_FUA	Combination of WRITE_FLUSH and FUA. The IO is preceded
 *			by a cache flush and data is guaranteed to be on
 *			non-volatile media on completion.
 *
 */
#define RW_MASK			REQ_WRITE
#define RWA_MASK		(1 << BIO_RW_AHEAD)

#define READ			0
#define WRITE			1
#define READA			RWA_MASK
#define SWRITE			(WRITE | READA)

#define READ_SYNC	(READ | (1 << BIO_RW_SYNCIO) | (1 << BIO_RW_UNPLUG))
#define READ_META	(READ | (1 << BIO_RW_META))
#define WRITE_SYNC_PLUG	(WRITE | (1 << BIO_RW_SYNCIO) | (1 << BIO_RW_NOIDLE))
#define WRITE_SYNC	(WRITE_SYNC_PLUG | (1 << BIO_RW_UNPLUG))
#define WRITE_ODIRECT_PLUG	(WRITE | (1 << BIO_RW_SYNCIO))
#define WRITE_META	(WRITE | (1 << BIO_RW_META))
#define SWRITE_SYNC_PLUG						\
			(SWRITE | (1 << BIO_RW_SYNCIO) | (1 << BIO_RW_NOIDLE))
#define SWRITE_SYNC	(SWRITE_SYNC_PLUG | (1 << BIO_RW_UNPLUG))
#define WRITE_BARRIER	(WRITE_SYNC | (1 << BIO_RW_BARRIER))

#define WRITE_FLUSH	(WRITE_SYNC | (1 << BIO_RW_FLUSH))
#define WRITE_FUA	(WRITE_SYNC | (1 << BIO_RW_FUA))
#define WRITE_FLUSH_FUA	(WRITE_FLUSH | WRITE_FUA)

/*
 * These aren't really reads or writes, they pass down information about
 * parts of device that are now unused by the file system.
 * DEPRECATED but preserved for compatibility.
 */
#define DISCARD_NOBARRIER (WRITE | (1 << BIO_RW_DISCARD))
#define DISCARD_BARRIER (DISCARD_NOBARRIER | (1 << BIO_RW_BARRIER))

#define SEL_IN		1
#define SEL_OUT		2
#define SEL_EX		4

/* public flags for file_system_type */
#define FS_REQUIRES_DEV 1 
#define FS_BINARY_MOUNTDATA 2
#define FS_HAS_SUBTYPE 4
#define FS_HAS_NEW_FREEZE 512	/* new freeze mechanism */
#define FS_REVAL_DOT	16384	/* Check the paths ".", ".." for staleness */
#define FS_RENAME_DOES_D_MOVE	32768	/* FS will handle d_move()
					 * during rename() internally.
					 */
#define FS_HANDLE_QUOTA		(1<<16)	/* FS handle quota disable/enable */
#define FS_WEAK_REVALIDATE (1<<17) /* FS has d_op->d_weak_revalidate. Must
					also have FS_REVAL_DOT set! */

 /*
  * the fs is built with the new s_writers member in the superblock
  * and uses all of that associated infrastructure
  */
#define sb_has_new_freeze(sb) ((sb)->s_type->fs_flags & FS_HAS_NEW_FREEZE)

/*
 * These are the fs-independent mount-flags: up to 32 flags are supported
 */
#define MS_RDONLY	 1	/* Mount read-only */
#define MS_NOSUID	 2	/* Ignore suid and sgid bits */
#define MS_NODEV	 4	/* Disallow access to device special files */
#define MS_NOEXEC	 8	/* Disallow program execution */
#define MS_SYNCHRONOUS	16	/* Writes are synced at once */
#define MS_REMOUNT	32	/* Alter flags of a mounted FS */
#define MS_MANDLOCK	64	/* Allow mandatory locks on an FS */
#define MS_DIRSYNC	128	/* Directory modifications are synchronous */
#define MS_NOATIME	1024	/* Do not update access times. */
#define MS_NODIRATIME	2048	/* Do not update directory access times */
#define MS_BIND		4096
#define MS_MOVE		8192
#define MS_REC		16384
#define MS_VERBOSE	32768	/* War is peace. Verbosity is silence.
				   MS_VERBOSE is deprecated. */
#define MS_SILENT	32768
#define MS_POSIXACL	(1<<16)	/* VFS does not apply the umask */
#define MS_UNBINDABLE	(1<<17)	/* change to unbindable */
#define MS_PRIVATE	(1<<18)	/* change to private */
#define MS_SLAVE	(1<<19)	/* change to slave */
#define MS_SHARED	(1<<20)	/* change to shared */
#define MS_RELATIME	(1<<21)	/* Update atime relative to mtime/ctime. */
#define MS_KERNMOUNT	(1<<22) /* this is a kern_mount call */
#define MS_I_VERSION	(1<<23) /* Update inode I_version field */
#define MS_STRICTATIME	(1<<24) /* Always perform atime updates */
#define MS_SNAP_STABLE  (1<<27) /* Snapshot pages during writeback, if needed */
#define MS_BORN		(1<<29)
#define MS_ACTIVE	(1<<30)
#define MS_NOUSER	(1<<31)

/*
 * Superblock flags that can be altered by MS_REMOUNT
 */
#define MS_RMT_MASK	(MS_RDONLY|MS_SYNCHRONOUS|MS_MANDLOCK|MS_I_VERSION)

/*
 * Old magic mount flag and mask
 */
#define MS_MGC_VAL 0xC0ED0000
#define MS_MGC_MSK 0xffff0000

/* Inode flags - they have nothing to superblock flags now */

#define S_SYNC		1	/* Writes are synced at once */
#define S_NOATIME	2	/* Do not update access times */
#define S_APPEND	4	/* Append-only file */
#define S_IMMUTABLE	8	/* Immutable file */
#define S_DEAD		16	/* removed, but still open directory */
#define S_NOQUOTA	32	/* Inode is not counted to quota */
#define S_DIRSYNC	64	/* Directory modifications are synchronous */
#define S_NOCMTIME	128	/* Do not update file c/mtime */
#define S_SWAPFILE	256	/* Do not truncate: swapon got its bmaps */
#define S_PRIVATE	512	/* Inode is fs-internal */
#define S_AUTOMOUNT	2048	/* Automount/referral quasi-directory */
#define S_AOP_EXT	16384 /* fs supports extended aops */

/*
 * Note that nosuid etc flags are inode-specific: setting some file-system
 * flags just means all the inodes inherit those flags by default. It might be
 * possible to override it selectively if you really wanted to with some
 * ioctl() that is not currently implemented.
 *
 * Exception: MS_RDONLY is always applied to the entire file system.
 *
 * Unfortunately, it is possible to change a filesystems flags with it mounted
 * with files in use.  This means that all of the inodes will not have their
 * i_flags updated.  Hence, i_flags no longer inherit the superblock mount
 * flags, so these have to be checked separately. -- rmk@arm.uk.linux.org
 */
#define __IS_FLG(inode,flg) ((inode)->i_sb->s_flags & (flg))

#define IS_RDONLY(inode) ((inode)->i_sb->s_flags & MS_RDONLY)
#define IS_SYNC(inode)		(__IS_FLG(inode, MS_SYNCHRONOUS) || \
					((inode)->i_flags & S_SYNC))
#define IS_DIRSYNC(inode)	(__IS_FLG(inode, MS_SYNCHRONOUS|MS_DIRSYNC) || \
					((inode)->i_flags & (S_SYNC|S_DIRSYNC)))
#define IS_MANDLOCK(inode)	__IS_FLG(inode, MS_MANDLOCK)
#define IS_NOATIME(inode)   __IS_FLG(inode, MS_RDONLY|MS_NOATIME)
#define IS_I_VERSION(inode)   __IS_FLG(inode, MS_I_VERSION)

#define IS_NOQUOTA(inode)	((inode)->i_flags & S_NOQUOTA)
#define IS_APPEND(inode)	((inode)->i_flags & S_APPEND)
#define IS_IMMUTABLE(inode)	((inode)->i_flags & S_IMMUTABLE)
#define IS_POSIXACL(inode)	__IS_FLG(inode, MS_POSIXACL)

#define IS_DEADDIR(inode)	((inode)->i_flags & S_DEAD)
#define IS_NOCMTIME(inode)	((inode)->i_flags & S_NOCMTIME)
#define IS_SWAPFILE(inode)	((inode)->i_flags & S_SWAPFILE)
#define IS_PRIVATE(inode)	((inode)->i_flags & S_PRIVATE)
#define IS_AUTOMOUNT(inode)	((inode)->i_flags & S_AUTOMOUNT)
#define IS_AOP_EXT(inode)	((inode)->i_flags & S_AOP_EXT)

/* the read-only stuff doesn't really belong here, but any other place is
   probably as bad and I don't want to create yet another include file. */

#define BLKROSET   _IO(0x12,93)	/* set device read-only (0 = read-write) */
#define BLKROGET   _IO(0x12,94)	/* get read-only status (0 = read_write) */
#define BLKRRPART  _IO(0x12,95)	/* re-read partition table */
#define BLKGETSIZE _IO(0x12,96)	/* return device size /512 (long *arg) */
#define BLKFLSBUF  _IO(0x12,97)	/* flush buffer cache */
#define BLKRASET   _IO(0x12,98)	/* set read ahead for block device */
#define BLKRAGET   _IO(0x12,99)	/* get current read ahead setting */
#define BLKFRASET  _IO(0x12,100)/* set filesystem (mm/filemap.c) read-ahead */
#define BLKFRAGET  _IO(0x12,101)/* get filesystem (mm/filemap.c) read-ahead */
#define BLKSECTSET _IO(0x12,102)/* set max sectors per request (ll_rw_blk.c) */
#define BLKSECTGET _IO(0x12,103)/* get max sectors per request (ll_rw_blk.c) */
#define BLKSSZGET  _IO(0x12,104)/* get block device sector size */
#if 0
#define BLKPG      _IO(0x12,105)/* See blkpg.h */

/* Some people are morons.  Do not use sizeof! */

#define BLKELVGET  _IOR(0x12,106,size_t)/* elevator get */
#define BLKELVSET  _IOW(0x12,107,size_t)/* elevator set */
/* This was here just to show that the number is taken -
   probably all these _IO(0x12,*) ioctls should be moved to blkpg.h. */
#endif
/* A jump here: 108-111 have been used for various private purposes. */
#define BLKBSZGET  _IOR(0x12,112,size_t)
#define BLKBSZSET  _IOW(0x12,113,size_t)
#define BLKGETSIZE64 _IOR(0x12,114,size_t)	/* return device size in bytes (u64 *arg) */
#define BLKTRACESETUP _IOWR(0x12,115,struct blk_user_trace_setup)
#define BLKTRACESTART _IO(0x12,116)
#define BLKTRACESTOP _IO(0x12,117)
#define BLKTRACETEARDOWN _IO(0x12,118)
#define BLKDISCARD _IO(0x12,119)
#define BLKIOMIN _IO(0x12,120)
#define BLKIOOPT _IO(0x12,121)
#define BLKALIGNOFF _IO(0x12,122)
#define BLKPBSZGET _IO(0x12,123)
#define BLKDISCARDZEROES _IO(0x12,124)

#define BMAP_IOCTL 1		/* obsolete - kept for compatibility */
#define FIBMAP	   _IO(0x00,1)	/* bmap access */
#define FIGETBSZ   _IO(0x00,2)	/* get the block size used for bmap */
#define FIFREEZE	_IOWR('X', 119, int)	/* Freeze */
#define FITHAW		_IOWR('X', 120, int)	/* Thaw */
#define FITRIM		_IOWR('X', 121, struct fstrim_range)	/* Trim */

#define	FS_IOC_GETFLAGS			_IOR('f', 1, long)
#define	FS_IOC_SETFLAGS			_IOW('f', 2, long)
#define	FS_IOC_GETVERSION		_IOR('v', 1, long)
#define	FS_IOC_SETVERSION		_IOW('v', 2, long)
#define FS_IOC_FIEMAP			_IOWR('f', 11, struct fiemap)
#define FS_IOC32_GETFLAGS		_IOR('f', 1, int)
#define FS_IOC32_SETFLAGS		_IOW('f', 2, int)
#define FS_IOC32_GETVERSION		_IOR('v', 1, int)
#define FS_IOC32_SETVERSION		_IOW('v', 2, int)

/*
 * Inode flags (FS_IOC_GETFLAGS / FS_IOC_SETFLAGS)
 */
#define	FS_SECRM_FL			0x00000001 /* Secure deletion */
#define	FS_UNRM_FL			0x00000002 /* Undelete */
#define	FS_COMPR_FL			0x00000004 /* Compress file */
#define FS_SYNC_FL			0x00000008 /* Synchronous updates */
#define FS_IMMUTABLE_FL			0x00000010 /* Immutable file */
#define FS_APPEND_FL			0x00000020 /* writes to file may only append */
#define FS_NODUMP_FL			0x00000040 /* do not dump file */
#define FS_NOATIME_FL			0x00000080 /* do not update atime */
/* Reserved for compression usage... */
#define FS_DIRTY_FL			0x00000100
#define FS_COMPRBLK_FL			0x00000200 /* One or more compressed clusters */
#define FS_NOCOMP_FL			0x00000400 /* Don't compress */
#define FS_ECOMPR_FL			0x00000800 /* Compression error */
/* End compression flags --- maybe not all used */
#define FS_BTREE_FL			0x00001000 /* btree format dir */
#define FS_INDEX_FL			0x00001000 /* hash-indexed directory */
#define FS_IMAGIC_FL			0x00002000 /* AFS directory */
#define FS_JOURNAL_DATA_FL		0x00004000 /* Reserved for ext3 */
#define FS_NOTAIL_FL			0x00008000 /* file tail should not be merged */
#define FS_DIRSYNC_FL			0x00010000 /* dirsync behaviour (directories only) */
#define FS_TOPDIR_FL			0x00020000 /* Top of directory hierarchies*/
#define FS_EXTENT_FL			0x00080000 /* Extents */
#define FS_DIRECTIO_FL			0x00100000 /* Use direct i/o */
#define FS_NOCOW_FL			0x00800000 /* Do not cow file */
#define FS_RESERVED_FL			0x80000000 /* reserved for ext2 lib */

#define FS_FL_USER_VISIBLE		0x0003DFFF /* User visible flags */
#define FS_FL_USER_MODIFIABLE		0x000380FF /* User modifiable flags */


#define SYNC_FILE_RANGE_WAIT_BEFORE	1
#define SYNC_FILE_RANGE_WRITE		2
#define SYNC_FILE_RANGE_WAIT_AFTER	4

#endif /* _LINUX_FS_H */
