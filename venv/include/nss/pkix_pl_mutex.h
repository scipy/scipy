/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_mutex.h
 *
 * Mutual Exclusion (Lock) Object Type Definition
 *
 */

#ifndef _PKIX_PL_MUTEX_H
#define _PKIX_PL_MUTEX_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_MutexStruct {
        PRLock* lock;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_Mutex_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_MUTEX_H */
