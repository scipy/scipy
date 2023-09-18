/*
 * include/linux/nfsd/const.h
 *
 * Various constants related to NFS.
 *
 * Copyright (C) 1995-1997 Olaf Kirch <okir@monad.swb.de>
 */

#ifndef _LINUX_NFSD_CONST_H
#define _LINUX_NFSD_CONST_H

#include <linux/nfs.h>
#include <linux/nfs2.h>
#include <linux/nfs3.h>
#include <linux/nfs4.h>

/*
 * Maximum protocol version supported by knfsd
 */
#define NFSSVC_MAXVERS		3

/*
 * Maximum blocksizes supported by daemon under various circumstances.
 */
#define NFSSVC_MAXBLKSIZE	RPCSVC_MAXPAYLOAD
/* NFSv2 is limited by the protocol specification, see RFC 1094 */
#define NFSSVC_MAXBLKSIZE_V2	(8*1024)


#endif /* _LINUX_NFSD_CONST_H */
