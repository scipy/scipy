/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PKIT_H
#define PKIT_H

/*
 * pkit.h
 *
 * This file contains definitions for the types of the top-level PKI objects.
 */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

#ifndef BASET_H
#include "baset.h"
#endif /* BASET_H */

#include "certt.h"
#include "pkcs11t.h"

#ifndef NSSPKIT_H
#include "nsspkit.h"
#endif /* NSSPKIT_H */

#ifndef NSSDEVT_H
#include "nssdevt.h"
#endif /* NSSDEVT_H */

#ifndef DEVT_H
#include "devt.h"
#endif /* DEVT_H */

#ifndef nssrwlkt_h__
#include "nssrwlkt.h"
#endif /* nssrwlkt_h__ */

PR_BEGIN_EXTERN_C

/*
 * A note on ephemeral certs
 *
 * The key objects defined here can only be created on tokens, and can only
 * exist on tokens.  Therefore, any instance of a key object must have
 * a corresponding cryptoki instance.  OTOH, certificates created in
 * crypto contexts need not be stored as session objects on the token.
 * There are good performance reasons for not doing so.  The certificate
 * and trust objects have been defined with a cryptoContext field to
 * allow for ephemeral certs, which may have a single instance in a crypto
 * context along with any number (including zero) of cryptoki instances.
 * Since contexts may not share objects, there can be only one context
 * for each object.
 */

typedef enum {
    nssPKILock = 1,
    nssPKIMonitor = 2
} nssPKILockType;

/* nssPKIObject
 *
 * This is the base object class, common to all PKI objects defined in
 * nsspkit.h
 */
struct nssPKIObjectStr {
    /* The arena for all object memory */
    NSSArena *arena;
    /* Atomically incremented/decremented reference counting */
    PRInt32 refCount;
    /* lock protects the array of nssCryptokiInstance's of the object */
    union {
        PZLock *lock;
        PZMonitor *mlock;
    } sync;
    nssPKILockType lockType;
    /* XXX with LRU cache, this cannot be guaranteed up-to-date.  It cannot
     * be compared against the update level of the trust domain, since it is
     * also affected by import/export.  Where is this array needed?
     */
    nssCryptokiObject **instances;
    PRUint32 numInstances;
    /* The object must live in a trust domain */
    NSSTrustDomain *trustDomain;
    /* The object may live in a crypto context */
    NSSCryptoContext *cryptoContext;
    /* XXX added so temp certs can have nickname, think more ... */
    NSSUTF8 *tempName;
};

typedef struct nssDecodedCertStr nssDecodedCert;

typedef struct nssCertificateStoreStr nssCertificateStore;

/* How wide is the scope of this? */
typedef struct nssSMIMEProfileStr nssSMIMEProfile;

typedef struct nssPKIObjectStr nssPKIObject;

struct NSSTrustStr {
    nssPKIObject object;
    NSSCertificate *certificate;
    nssTrustLevel serverAuth;
    nssTrustLevel clientAuth;
    nssTrustLevel emailProtection;
    nssTrustLevel codeSigning;
    PRBool stepUpApproved;
};

struct nssSMIMEProfileStr {
    nssPKIObject object;
    NSSCertificate *certificate;
    NSSASCII7 *email;
    NSSDER *subject;
    NSSItem *profileTime;
    NSSItem *profileData;
};

struct NSSCertificateStr {
    nssPKIObject object;
    NSSCertificateType type;
    NSSItem id;
    NSSBER encoding;
    NSSDER issuer;
    NSSDER subject;
    NSSDER serial;
    NSSASCII7 *email;
    nssDecodedCert *decoding;
};

struct NSSPrivateKeyStr;

struct NSSPublicKeyStr;

struct NSSSymmetricKeyStr;

typedef struct nssTDCertificateCacheStr nssTDCertificateCache;

struct NSSTrustDomainStr {
    PRInt32 refCount;
    NSSArena *arena;
    NSSCallback *defaultCallback;
    nssList *tokenList;
    nssListIterator *tokens;
    nssTDCertificateCache *cache;
    NSSRWLock *tokensLock;
    void *spkDigestInfo;
    CERTStatusConfig *statusConfig;
};

struct NSSCryptoContextStr {
    PRInt32 refCount;
    NSSArena *arena;
    NSSTrustDomain *td;
    NSSToken *token;
    nssSession *session;
    nssCertificateStore *certStore;
};

struct NSSTimeStr {
    PRTime prTime;
};

struct NSSCRLStr {
    nssPKIObject object;
    NSSDER encoding;
    NSSUTF8 *url;
    PRBool isKRL;
};

typedef struct NSSCRLStr NSSCRL;

struct NSSPoliciesStr;

struct NSSAlgorithmAndParametersStr;

struct NSSPKIXCertificateStr;

PR_END_EXTERN_C

#endif /* PKIT_H */
