/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef CKFWM_H
#define CKFWM_H

/*
 * ckfwm.h
 *
 * This file prototypes the module-private calls of the NSS Cryptoki Framework.
 */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

#ifndef NSSCKT_H
#include "nssckt.h"
#endif /* NSSCKT_H */

#ifndef NSSCKFWT_H
#include "nssckfwt.h"
#endif /* NSSCKFWT_H */

/*
 * nssCKFWHash
 *
 *  nssCKFWHash_Create
 *  nssCKFWHash_Destroy
 *  nssCKFWHash_Add
 *  nssCKFWHash_Remove
 *  nssCKFWHash_Count
 *  nssCKFWHash_Exists
 *  nssCKFWHash_Lookup
 *  nssCKFWHash_Iterate
 */

/*
 * nssCKFWHash_Create
 *
 */
NSS_EXTERN nssCKFWHash *
nssCKFWHash_Create(
    NSSCKFWInstance *fwInstance,
    NSSArena *arena,
    CK_RV *pError);

/*
 * nssCKFWHash_Destroy
 *
 */
NSS_EXTERN void
nssCKFWHash_Destroy(
    nssCKFWHash *hash);

/*
 * nssCKFWHash_Add
 *
 */
NSS_EXTERN CK_RV
nssCKFWHash_Add(
    nssCKFWHash *hash,
    const void *key,
    const void *value);

/*
 * nssCKFWHash_Remove
 *
 */
NSS_EXTERN void
nssCKFWHash_Remove(
    nssCKFWHash *hash,
    const void *it);

/*
 * nssCKFWHash_Count
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWHash_Count(
    nssCKFWHash *hash);

/*
 * nssCKFWHash_Exists
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWHash_Exists(
    nssCKFWHash *hash,
    const void *it);

/*
 * nssCKFWHash_Lookup
 *
 */
NSS_EXTERN void *
nssCKFWHash_Lookup(
    nssCKFWHash *hash,
    const void *it);

/*
 * nssCKFWHash_Iterate
 *
 */
NSS_EXTERN void
nssCKFWHash_Iterate(
    nssCKFWHash *hash,
    nssCKFWHashIterator fcn,
    void *closure);

#endif /* CKFWM_H */
