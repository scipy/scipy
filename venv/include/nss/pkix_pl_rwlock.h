/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_rwlock.h
 *
 * Read/Write Lock Definition
 *
 */

#ifndef _PKIX_PL_RWLOCK_H
#define _PKIX_PL_RWLOCK_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_RWLockStruct {
        PRRWLock* lock;
        PKIX_UInt32 readCount;
        PKIX_Boolean writeLocked;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_RWLock_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_RWLOCK_H */
