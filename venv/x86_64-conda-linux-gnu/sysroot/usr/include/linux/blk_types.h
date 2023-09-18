/*
 * Block data types and constants.  Directly include this file only to
 * break include dependency loop.
 */
#ifndef __LINUX_BLK_TYPES_H
#define __LINUX_BLK_TYPES_H

#ifdef CONFIG_BLOCK

#include <linux/types.h>

/*
 * bio flags
 */
#define BIO_UPTODATE	0	/* ok after I/O completion */
#define BIO_RW_BLOCK	1	/* RW_AHEAD set, and read/write would block */
#define BIO_EOF		2	/* out-out-bounds error */
#define BIO_SEG_VALID	3	/* bi_phys_segments valid */
#define BIO_CLONED	4	/* doesn't own data */
#define BIO_BOUNCED	5	/* bio is a bounce bio */
#define BIO_USER_MAPPED 6	/* contains user pages */
#define BIO_EOPNOTSUPP	7	/* not supported */
#define BIO_CPU_AFFINE	8	/* complete bio on same CPU as submitted */
#define BIO_NULL_MAPPED 9	/* contains invalid user pages */
#define BIO_FS_INTEGRITY 10	/* fs owns integrity data, not block layer */
#define BIO_QUIET	11	/* Make BIO Quiet */
#define BIO_MAPPED_INTEGRITY 12	/* integrity metadata has been remapped */
#define bio_flagged(bio, flag)	((bio)->bi_flags & (1 << (flag)))

/*
 * top 4 bits of bio flags indicate the pool this bio came from
 */
#define BIO_POOL_BITS		(4)
#define BIO_POOL_NONE		((1UL << BIO_POOL_BITS) - 1)
#define BIO_POOL_OFFSET		(BITS_PER_LONG - BIO_POOL_BITS)
#define BIO_POOL_MASK		(1UL << BIO_POOL_OFFSET)
#define BIO_POOL_IDX(bio)	((bio)->bi_flags >> BIO_POOL_OFFSET)	

#endif /* CONFIG_BLOCK */

/*
 * bio bi_rw flags
 *
 * bit 0 -- data direction
 *	If not set, bio is a read from device. If set, it's a write to device.
 * bit 1 -- fail fast device errors
 * bit 2 -- fail fast transport errors
 * bit 3 -- fail fast driver errors
 * bit 4 -- rw-ahead when set
 * bit 5 -- barrier
 *	Insert a serialization point in the IO queue, forcing previously
 *	submitted IO to be completed before this one is issued.
 * bit 6 -- synchronous I/O hint.
 * bit 7 -- Unplug the device immediately after submitting this bio.
 * bit 8 -- metadata request
 *	Used for tracing to differentiate metadata and data IO. May also
 *	get some preferential treatment in the IO scheduler
 * bit 9 -- discard sectors
 *	Informs the lower level device that this range of sectors is no longer
 *	used by the file system and may thus be freed by the device. Used
 *	for flash based storage.
 *	Don't want driver retries for any fast fail whatever the reason.
 * bit 10 -- Tell the IO scheduler not to wait for more requests after this
 *	one has been submitted, even if it is a SYNC request.
 * bit 11 -- FLUSH
 * bit 12 -- FUA
 */
enum bio_rw_flags {
	BIO_RW,
	BIO_RW_FAILFAST_DEV,
	BIO_RW_FAILFAST_TRANSPORT,
	BIO_RW_FAILFAST_DRIVER,
	/* above flags must match REQ_* */
	BIO_RW_AHEAD,
	BIO_RW_BARRIER,
	BIO_RW_SYNCIO,
	BIO_RW_UNPLUG,
	BIO_RW_META,
	BIO_RW_DISCARD,
	BIO_RW_NOIDLE,
	/*
	 * FLUSH and FUA are so bio-based drivers (DM/MD) and FS
	 * have something to work with
	 */
	BIO_RW_FLUSH,
	BIO_RW_FUA,
	/*
	 * This bio has already been subjected to throttling rules. Don't
	 * do it again.
	 */
	BIO_RW_THROTTLED,
	BIO_RW_NOMERGE,
};

#define BIO_WRITE		(1 << BIO_RW)
#define BIO_FLUSH		(1 << BIO_RW_FLUSH)
#define BIO_FUA			(1 << BIO_RW_FUA)
#define BIO_DISCARD		(1 << BIO_RW_DISCARD)
#define BIO_NOMERGE		(1 << BIO_RW_NOMERGE)

/* This mask is used for bio merge checking (RHEL6 doesn't share bio and request flags) */
#define BIO_NOMERGE_FLAGS \
	((1 << BIO_RW_BARRIER) | BIO_FLUSH | BIO_FUA | BIO_NOMERGE)

/*
 * request type modified bits. first four bits match BIO_RW* bits, important
 */
enum rq_flag_bits {
	__REQ_WRITE,		/* was __REQ_RW, not set, read. set, write */
	__REQ_FAILFAST_DEV,	/* no driver retries of device errors */
	__REQ_FAILFAST_TRANSPORT, /* no driver retries of transport errors */
	__REQ_FAILFAST_DRIVER,	/* no driver retries of driver errors */
	/* above flags must match BIO_RW_* */
	__REQ_DISCARD,		/* request to discard sectors */
	__REQ_SORTED,		/* elevator knows about this request */
	__REQ_SOFTBARRIER,	/* may not be passed by ioscheduler */
	__REQ_HARDBARRIER,	/* DEPRECATED: may not be passed by drive either */
	__REQ_FUA,		/* forced unit access */
	__REQ_NOMERGE,		/* don't touch this for merging */
	__REQ_STARTED,		/* drive already may have started this one */
	__REQ_DONTPREP,		/* don't call prep for this one */
	__REQ_QUEUED,		/* uses queueing */
	__REQ_ELVPRIV,		/* elevator private data attached */
	__REQ_FAILED,		/* set if the request failed */
	__REQ_QUIET,		/* don't worry about errors */
	__REQ_PREEMPT,		/* set for "ide_preempt" requests */
	__REQ_ORDERED_COLOR,	/* DEPRECATED: is before or after barrier */
	__REQ_SYNC,		/* was __REQ_RW_SYNC, request is sync (sync write or read) */
	__REQ_ALLOCED,		/* request came from our alloc pool */
	__REQ_META,		/* was __REQ_RW_META, metadata io request */
	__REQ_COPY_USER,	/* contains copies of user pages */
	__REQ_INTEGRITY,	/* DEPRECATED: integrity metadata has been remapped */
	__REQ_NOIDLE,		/* don't anticipate more IO after this one */
	__REQ_IO_STAT,		/* account I/O stat */
	__REQ_MIXED_MERGE,	/* merge of different types, fail separately */
	__REQ_FLUSH,		/* request for cache flush */
	__REQ_FLUSH_SEQ,	/* request for flush sequence */
	__REQ_NR_BITS,		/* stops here */
};

#define REQ_WRITE		(1 << __REQ_WRITE)
#define REQ_RW			REQ_WRITE  /* DEPRECATED */

#define REQ_FAILFAST_DEV	(1 << __REQ_FAILFAST_DEV)
#define REQ_FAILFAST_TRANSPORT	(1 << __REQ_FAILFAST_TRANSPORT)
#define REQ_FAILFAST_DRIVER	(1 << __REQ_FAILFAST_DRIVER)
#define REQ_HARDBARRIER		(1 << __REQ_HARDBARRIER)

#define REQ_SYNC		(1 << __REQ_SYNC)
#define REQ_RW_SYNC		REQ_SYNC  /* DEPRECATED */

#define REQ_META		(1 << __REQ_META)
#define REQ_RW_META		REQ_META  /* DEPRECATED */

#define REQ_DISCARD		(1 << __REQ_DISCARD)
#define REQ_NOIDLE		(1 << __REQ_NOIDLE)

#define REQ_FAILFAST_MASK \
	(REQ_FAILFAST_DEV | REQ_FAILFAST_TRANSPORT | REQ_FAILFAST_DRIVER)
#define REQ_COMMON_MASK \
	(REQ_WRITE | REQ_FAILFAST_MASK | REQ_HARDBARRIER | REQ_SYNC | \
	 REQ_META | REQ_DISCARD | REQ_NOIDLE | REQ_FLUSH | REQ_FUA)
#define REQ_CLONE_MASK		REQ_COMMON_MASK

/* This mask is used for request merge checking */
#define REQ_NOMERGE_FLAGS \
	(REQ_NOMERGE | REQ_STARTED | REQ_HARDBARRIER | REQ_SOFTBARRIER | REQ_FLUSH | REQ_FUA)

#define REQ_SORTED		(1 << __REQ_SORTED)
#define REQ_SOFTBARRIER		(1 << __REQ_SOFTBARRIER)
#define REQ_FUA			(1 << __REQ_FUA)
#define REQ_NOMERGE		(1 << __REQ_NOMERGE)
#define REQ_STARTED		(1 << __REQ_STARTED)
#define REQ_DONTPREP		(1 << __REQ_DONTPREP)
#define REQ_QUEUED		(1 << __REQ_QUEUED)
#define REQ_ELVPRIV		(1 << __REQ_ELVPRIV)
#define REQ_FAILED		(1 << __REQ_FAILED)
#define REQ_QUIET		(1 << __REQ_QUIET)
#define REQ_PREEMPT		(1 << __REQ_PREEMPT)
#define REQ_ORDERED_COLOR	(1 << __REQ_ORDERED_COLOR)  /* DEPRECATED */
#define REQ_ALLOCED		(1 << __REQ_ALLOCED)
#define REQ_COPY_USER		(1 << __REQ_COPY_USER)
#define REQ_INTEGRITY		(1 << __REQ_INTEGRITY)
#define REQ_FLUSH		(1 << __REQ_FLUSH)
#define REQ_FLUSH_SEQ		(1 << __REQ_FLUSH_SEQ)
#define REQ_IO_STAT		(1 << __REQ_IO_STAT)
#define REQ_MIXED_MERGE		(1 << __REQ_MIXED_MERGE)

#define REQ_WRITE_FLUSH		(REQ_WRITE | REQ_SYNC | REQ_NOIDLE | REQ_FLUSH)

#endif /* __LINUX_BLK_TYPES_H */
