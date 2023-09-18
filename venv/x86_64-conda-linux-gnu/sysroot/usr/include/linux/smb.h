/*
 *  smb.h
 *
 *  Copyright (C) 1995, 1996 by Paal-Kr. Engstad and Volker Lendecke
 *  Copyright (C) 1997 by Volker Lendecke
 *
 */

#ifndef _LINUX_SMB_H
#define _LINUX_SMB_H

#include <linux/types.h>
#include <linux/magic.h>

enum smb_protocol { 
	SMB_PROTOCOL_NONE, 
	SMB_PROTOCOL_CORE, 
	SMB_PROTOCOL_COREPLUS, 
	SMB_PROTOCOL_LANMAN1, 
	SMB_PROTOCOL_LANMAN2, 
	SMB_PROTOCOL_NT1 
};

enum smb_case_hndl {
	SMB_CASE_DEFAULT,
	SMB_CASE_LOWER,
	SMB_CASE_UPPER
};

struct smb_dskattr {
        __u16 total;
        __u16 allocblocks;
        __u16 blocksize;
        __u16 free;
};

struct smb_conn_opt {

        /* The socket */
	unsigned int fd;

	enum smb_protocol protocol;
	enum smb_case_hndl case_handling;

	/* Connection-Options */

	__u32              max_xmit;
	__u16              server_uid;
	__u16              tid;

        /* The following are LANMAN 1.0 options */
        __u16              secmode;
        __u16              maxmux;
        __u16              maxvcs;
        __u16              rawmode;
        __u32              sesskey;

	/* The following are NT LM 0.12 options */
	__u32              maxraw;
	__u32              capabilities;
	__s16              serverzone;
};

#endif
