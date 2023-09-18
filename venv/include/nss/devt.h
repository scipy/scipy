/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef DEVT_H
#define DEVT_H

/*
 * devt.h
 *
 * This file contains definitions for the low-level cryptoki devices.
 */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

#ifndef NSSPKIT_H
#include "nsspkit.h"
#endif /* NSSPKIT_H */

#ifndef NSSDEVT_H
#include "nssdevt.h"
#endif /* NSSDEVT_H */

#ifndef BASET_H
#include "baset.h"
#endif /* BASET_H */

#include "secmodt.h"

PR_BEGIN_EXTERN_C

typedef struct nssSessionStr nssSession;

/* XXX until NSSTokenStr is moved */
struct nssDeviceBaseStr {
    NSSArena *arena;
    PZLock *lock;
    PRInt32 refCount;
    NSSUTF8 *name;
    PRUint32 flags;
};

typedef struct nssTokenObjectCacheStr nssTokenObjectCache;

/* XXX until devobject.c goes away */
struct NSSTokenStr {
    struct nssDeviceBaseStr base;
    NSSSlot *slot;    /* Parent (or peer, if you will) */
    CK_FLAGS ckFlags; /* from CK_TOKEN_INFO.flags */
    PRUint32 flags;
    void *epv;
    nssSession *defaultSession;
    NSSTrustDomain *trustDomain;
    PRIntervalTime lastTime;
    nssTokenObjectCache *cache;
    PK11SlotInfo *pk11slot;
};

typedef enum {
    nssSlotAskPasswordTimes_FirstTime = 0,
    nssSlotAskPasswordTimes_EveryTime = 1,
    nssSlotAskPasswordTimes_Timeout = 2
} nssSlotAskPasswordTimes;

struct nssSlotAuthInfoStr {
    PRTime lastLogin;
    nssSlotAskPasswordTimes askTimes;
    PRIntervalTime askPasswordTimeout;
};

/* values for lastTokenPingState */
typedef enum {
    nssSlotLastPingState_Reset = 0,  /* the state has just been reset, discard
                                      * our cache */
    nssSlotLastPingState_Update = 1, /* we are updating the lastTokenPingTime */
    nssSlotLastPingState_Valid = 2,  /* lastTokenPingTime is valid */
} nssSlotLastPingState;

struct NSSSlotStr {
    struct nssDeviceBaseStr base;
    NSSModule *module; /* Parent */
    CK_SLOT_ID slotID;
    CK_FLAGS ckFlags; /* from CK_SLOT_INFO.flags */
    struct nssSlotAuthInfoStr authInfo;
    PRIntervalTime lastTokenPingTime;
    nssSlotLastPingState lastTokenPingState;
    PZLock *lock;
    void *epv;
    PK11SlotInfo *pk11slot;
    PZLock *isPresentLock;
    PRCondVar *isPresentCondition;
    PRThread *isPresentThread;
};

struct nssSessionStr {
    /* Must not hold slot->lock when taking lock.
     * See ordering in nssSlot_IsTokenPresent.
     */
    PZLock *lock;
    CK_SESSION_HANDLE handle;
    NSSSlot *slot;
    PRBool isRW;
    PRBool ownLock;
};

typedef enum {
    NSSCertificateType_Unknown = 0,
    NSSCertificateType_PKIX = 1
} NSSCertificateType;

typedef enum {
    nssTrustLevel_Unknown = 0,
    nssTrustLevel_NotTrusted = 1,
    nssTrustLevel_Trusted = 2,
    nssTrustLevel_TrustedDelegator = 3,
    nssTrustLevel_MustVerify = 4,
    nssTrustLevel_ValidDelegator = 5
} nssTrustLevel;

typedef struct nssCryptokiInstanceStr nssCryptokiInstance;

struct nssCryptokiInstanceStr {
    CK_OBJECT_HANDLE handle;
    NSSToken *token;
    PRBool isTokenObject;
    NSSUTF8 *label;
};

typedef struct nssCryptokiInstanceStr nssCryptokiObject;

typedef struct nssTokenCertSearchStr nssTokenCertSearch;

typedef enum {
    nssTokenSearchType_AllObjects = 0,
    nssTokenSearchType_SessionOnly = 1,
    nssTokenSearchType_TokenOnly = 2,
    nssTokenSearchType_TokenForced = 3
} nssTokenSearchType;

struct nssTokenCertSearchStr {
    nssTokenSearchType searchType;
    PRStatus (*callback)(NSSCertificate *c, void *arg);
    void *cbarg;
    nssList *cached;
    /* TODO: add a cache query callback if the list would be large
     *       (traversal)
     */
};

struct nssSlotListStr;
typedef struct nssSlotListStr nssSlotList;

struct NSSAlgorithmAndParametersStr {
    CK_MECHANISM mechanism;
};

PR_END_EXTERN_C

#endif /* DEVT_H */
