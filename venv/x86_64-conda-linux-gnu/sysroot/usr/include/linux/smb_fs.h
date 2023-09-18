/*
 *  smb_fs.h
 *
 *  Copyright (C) 1995 by Paal-Kr. Engstad and Volker Lendecke
 *  Copyright (C) 1997 by Volker Lendecke
 *
 */

#ifndef _LINUX_SMB_FS_H
#define _LINUX_SMB_FS_H

#include <linux/smb.h>

/*
 * ioctl commands
 */
#define	SMB_IOC_GETMOUNTUID		_IOR('u', 1, __kernel_old_uid_t)
#define SMB_IOC_NEWCONN                 _IOW('u', 2, struct smb_conn_opt)

/* __kernel_uid_t can never change, so we have to use __kernel_uid32_t */
#define	SMB_IOC_GETMOUNTUID32		_IOR('u', 3, __kernel_uid32_t)



#endif /* _LINUX_SMB_FS_H */
