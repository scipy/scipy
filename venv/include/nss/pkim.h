/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PKIM_H
#define PKIM_H

#ifndef BASE_H
#include "base.h"
#endif /* BASE_H */

#ifndef PKI_H
#include "pki.h"
#endif /* PKI_H */

#ifndef PKITM_H
#include "pkitm.h"
#endif /* PKITM_H */

PR_BEGIN_EXTERN_C

/* nssPKIObject
 *
 * This is the base object class, common to all PKI objects defined in
 * in this module.  Each object can be safely 'casted' to an nssPKIObject,
 * then passed to these methods.
 *
 * nssPKIObject_Create
 * nssPKIObject_Destroy
 * nssPKIObject_AddRef
 * nssPKIObject_AddInstance
 * nssPKIObject_HasInstance
 * nssPKIObject_GetTokens
 * nssPKIObject_GetNicknameForToken
 * nssPKIObject_RemoveInstanceForToken
 * nssPKIObject_DeleteStoredObject
 */

NSS_EXTERN void nssPKIObject_Lock(nssPKIObject *object);
NSS_EXTERN void nssPKIObject_Unlock(nssPKIObject *object);
NSS_EXTERN PRStatus nssPKIObject_NewLock(nssPKIObject *object,
                                         nssPKILockType lockType);
NSS_EXTERN void nssPKIObject_DestroyLock(nssPKIObject *object);

/* nssPKIObject_Create
 *
 * A generic PKI object.  It must live in a trust domain.  It may be
 * initialized with a token instance, or alternatively in a crypto context.
 */
NSS_EXTERN nssPKIObject *
nssPKIObject_Create(
    NSSArena *arenaOpt,
    nssCryptokiObject *instanceOpt,
    NSSTrustDomain *td,
    NSSCryptoContext *ccOpt,
    nssPKILockType lockType);

/* nssPKIObject_AddRef
 */
NSS_EXTERN nssPKIObject *
nssPKIObject_AddRef(nssPKIObject *object);

/* nssPKIObject_Destroy
 *
 * Returns true if object was destroyed.  This notifies the subclass that
 * all references are gone and it should delete any members it owns.
 */
NSS_EXTERN PRBool
nssPKIObject_Destroy(nssPKIObject *object);

/* nssPKIObject_AddInstance
 *
 * Add a token instance to the object, if it does not have it already.
 */
NSS_EXTERN PRStatus
nssPKIObject_AddInstance(
    nssPKIObject *object,
    nssCryptokiObject *instance);

/* nssPKIObject_HasInstance
 *
 * Query the object for a token instance.
 */
NSS_EXTERN PRBool
nssPKIObject_HasInstance(
    nssPKIObject *object,
    nssCryptokiObject *instance);

/* nssPKIObject_GetTokens
 *
 * Get all tokens which have an instance of the object.
 */
NSS_EXTERN NSSToken **
nssPKIObject_GetTokens(
    nssPKIObject *object,
    PRStatus *statusOpt);

/* nssPKIObject_GetNicknameForToken
 *
 * tokenOpt == NULL means take the first available, otherwise return the
 * nickname for the specified token.
 */
NSS_EXTERN NSSUTF8 *
nssPKIObject_GetNicknameForToken(
    nssPKIObject *object,
    NSSToken *tokenOpt);

/* nssPKIObject_RemoveInstanceForToken
 *
 * Remove the instance of the object on the specified token.
 */
NSS_EXTERN PRStatus
nssPKIObject_RemoveInstanceForToken(
    nssPKIObject *object,
    NSSToken *token);

/* nssPKIObject_DeleteStoredObject
 *
 * Delete all token instances of the object, as well as any crypto context
 * instances (TODO).  If any of the instances are read-only, or if the
 * removal fails, the object will keep those instances.  'isFriendly' refers
 * to the object -- can this object be removed from a friendly token without
 * login?  For example, certificates are friendly, private keys are not.
 * Note that if the token is not friendly, authentication will be required
 * regardless of the value of 'isFriendly'.
 */
NSS_EXTERN PRStatus
nssPKIObject_DeleteStoredObject(
    nssPKIObject *object,
    NSSCallback *uhh,
    PRBool isFriendly);

NSS_EXTERN nssCryptokiObject **
nssPKIObject_GetInstances(
    nssPKIObject *object);

NSS_EXTERN NSSCertificate **
nssTrustDomain_FindCertificatesByID(
    NSSTrustDomain *td,
    NSSItem *id,
    NSSCertificate **rvOpt,
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSCRL **
nssTrustDomain_FindCRLsBySubject(
    NSSTrustDomain *td,
    NSSDER *subject);

/* module-private nsspki methods */

NSS_EXTERN NSSCryptoContext *
nssCryptoContext_Create(
    NSSTrustDomain *td,
    NSSCallback *uhhOpt);

/* XXX for the collection */
NSS_EXTERN NSSCertificate *
nssCertificate_Create(nssPKIObject *object);

NSS_EXTERN PRStatus
nssCertificate_SetCertTrust(
    NSSCertificate *c,
    NSSTrust *trust);

NSS_EXTERN nssDecodedCert *
nssCertificate_GetDecoding(NSSCertificate *c);

extern PRIntn
nssCertificate_SubjectListSort(
    void *v1,
    void *v2);

NSS_EXTERN nssDecodedCert *
nssDecodedCert_Create(
    NSSArena *arenaOpt,
    NSSDER *encoding,
    NSSCertificateType type);

NSS_EXTERN PRStatus
nssDecodedCert_Destroy(nssDecodedCert *dc);

NSS_EXTERN NSSTrust *
nssTrust_Create(
    nssPKIObject *object,
    NSSItem *certData);

NSS_EXTERN NSSCRL *
nssCRL_Create(nssPKIObject *object);

NSS_EXTERN NSSCRL *
nssCRL_AddRef(NSSCRL *crl);

NSS_EXTERN PRStatus
nssCRL_Destroy(NSSCRL *crl);

NSS_EXTERN PRStatus
nssCRL_DeleteStoredObject(
    NSSCRL *crl,
    NSSCallback *uhh);

NSS_EXTERN NSSPrivateKey *
nssPrivateKey_Create(nssPKIObject *o);

NSS_EXTERN NSSDER *
nssCRL_GetEncoding(NSSCRL *crl);

NSS_EXTERN NSSPublicKey *
nssPublicKey_Create(nssPKIObject *object);

/* nssCertificateArray
 *
 * These are being thrown around a lot, might as well group together some
 * functionality.
 *
 * nssCertificateArray_Destroy
 * nssCertificateArray_Join
 * nssCertificateArray_FindBestCertificate
 * nssCertificateArray_Traverse
 */

/* nssCertificateArray_Destroy
 *
 * Will destroy the array and the certs within it.  If the array was created
 * in an arena, will *not* (of course) destroy the arena.  However, is safe
 * to call this method on an arena-allocated array.
 */
NSS_EXTERN void
nssCertificateArray_Destroy(NSSCertificate **certs);

/* nssCertificateArray_Join
 *
 * Join two arrays into one.  The two arrays, certs1 and certs2, should
 * be considered invalid after a call to this function (they may be destroyed
 * as part of the join).  certs1 and/or certs2 may be NULL.  Safe to
 * call with arrays allocated in an arena, the result will also be in the
 * arena.
 */
NSS_EXTERN NSSCertificate **
nssCertificateArray_Join(
    NSSCertificate **certs1,
    NSSCertificate **certs2);

/* nssCertificateArray_FindBestCertificate
 *
 * Use the usual { time, usage, policies } to find the best cert in the
 * array.
 */
NSS_EXTERN NSSCertificate *
nssCertificateArray_FindBestCertificate(
    NSSCertificate **certs,
    NSSTime *timeOpt,
    const NSSUsage *usage,
    NSSPolicies *policiesOpt);

/* nssCertificateArray_Traverse
 *
 * Do the callback for each cert, terminate the traversal if the callback
 * fails.
 */
NSS_EXTERN PRStatus
nssCertificateArray_Traverse(
    NSSCertificate **certs,
    PRStatus (*callback)(NSSCertificate *c, void *arg),
    void *arg);

NSS_EXTERN void
nssCRLArray_Destroy(NSSCRL **crls);

/* nssPKIObjectCollection
 *
 * This is a handy way to group objects together and perform operations
 * on them.  It can also handle "proto-objects"-- references to
 * objects instances on tokens, where the actual object hasn't
 * been formed yet.
 *
 * nssCertificateCollection_Create
 * nssPrivateKeyCollection_Create
 * nssPublicKeyCollection_Create
 *
 * If this was a language that provided for inheritance, each type would
 * inherit all of the following methods.  Instead, there is only one
 * type (nssPKIObjectCollection), shared among all.  This may cause
 * confusion; an alternative would be to define all of the methods
 * for each subtype (nssCertificateCollection_Destroy, ...), but that doesn't
 * seem worth the code bloat..  It is left up to the caller to remember
 * what type of collection he/she is dealing with.
 *
 * nssPKIObjectCollection_Destroy
 * nssPKIObjectCollection_Count
 * nssPKIObjectCollection_AddObject
 * nssPKIObjectCollection_AddInstances
 * nssPKIObjectCollection_Traverse
 *
 * Back to type-specific methods.
 *
 * nssPKIObjectCollection_GetCertificates
 * nssPKIObjectCollection_GetCRLs
 * nssPKIObjectCollection_GetPrivateKeys
 * nssPKIObjectCollection_GetPublicKeys
 */

/* nssCertificateCollection_Create
 *
 * Create a collection of certificates in the specified trust domain.
 * Optionally provide a starting set of certs.
 */
NSS_EXTERN nssPKIObjectCollection *
nssCertificateCollection_Create(
    NSSTrustDomain *td,
    NSSCertificate **certsOpt);

/* nssCRLCollection_Create
 *
 * Create a collection of CRLs/KRLs in the specified trust domain.
 * Optionally provide a starting set of CRLs.
 */
NSS_EXTERN nssPKIObjectCollection *
nssCRLCollection_Create(
    NSSTrustDomain *td,
    NSSCRL **crlsOpt);

/* nssPrivateKeyCollection_Create
 *
 * Create a collection of private keys in the specified trust domain.
 * Optionally provide a starting set of keys.
 */
NSS_EXTERN nssPKIObjectCollection *
nssPrivateKeyCollection_Create(
    NSSTrustDomain *td,
    NSSPrivateKey **pvkOpt);

/* nssPublicKeyCollection_Create
 *
 * Create a collection of public keys in the specified trust domain.
 * Optionally provide a starting set of keys.
 */
NSS_EXTERN nssPKIObjectCollection *
nssPublicKeyCollection_Create(
    NSSTrustDomain *td,
    NSSPublicKey **pvkOpt);

/* nssPKIObjectCollection_Destroy
 */
NSS_EXTERN void
nssPKIObjectCollection_Destroy(nssPKIObjectCollection *collection);

/* nssPKIObjectCollection_Count
 */
NSS_EXTERN PRUint32
nssPKIObjectCollection_Count(nssPKIObjectCollection *collection);

NSS_EXTERN PRStatus
nssPKIObjectCollection_AddObject(
    nssPKIObjectCollection *collection,
    nssPKIObject *object);

/* nssPKIObjectCollection_AddInstances
 *
 * Add a set of object instances to the collection.  The instances
 * will be sorted into any existing certs/proto-certs that may be in
 * the collection.  The instances will be absorbed by the collection,
 * the array should not be used after this call (except to free it).
 *
 * Failure means the collection is in an invalid state.
 *
 * numInstances = 0 means the array is NULL-terminated
 */
NSS_EXTERN PRStatus
nssPKIObjectCollection_AddInstances(
    nssPKIObjectCollection *collection,
    nssCryptokiObject **instances,
    PRUint32 numInstances);

/* nssPKIObjectCollection_Traverse
 */
NSS_EXTERN PRStatus
nssPKIObjectCollection_Traverse(
    nssPKIObjectCollection *collection,
    nssPKIObjectCallback *callback);

/* This function is being added for NSS 3.5.  It corresponds to the function
 * nssToken_TraverseCertificates.  The idea is to use the collection during
 * a traversal, creating certs each time a new instance is added for which
 * a cert does not already exist.
 */
NSS_EXTERN PRStatus
nssPKIObjectCollection_AddInstanceAsObject(
    nssPKIObjectCollection *collection,
    nssCryptokiObject *instance);

/* nssPKIObjectCollection_GetCertificates
 *
 * Get all of the certificates in the collection.
 */
NSS_EXTERN NSSCertificate **
nssPKIObjectCollection_GetCertificates(
    nssPKIObjectCollection *collection,
    NSSCertificate **rvOpt,
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSCRL **
nssPKIObjectCollection_GetCRLs(
    nssPKIObjectCollection *collection,
    NSSCRL **rvOpt,
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSPrivateKey **
nssPKIObjectCollection_GetPrivateKeys(
    nssPKIObjectCollection *collection,
    NSSPrivateKey **rvOpt,
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSPublicKey **
nssPKIObjectCollection_GetPublicKeys(
    nssPKIObjectCollection *collection,
    NSSPublicKey **rvOpt,
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSTime *
NSSTime_Now(NSSTime *timeOpt);

NSS_EXTERN NSSTime *
NSSTime_SetPRTime(
    NSSTime *timeOpt,
    PRTime prTime);

NSS_EXTERN PRTime
NSSTime_GetPRTime(
    NSSTime *time);

NSS_EXTERN nssHash *
nssHash_CreateCertificate(
    NSSArena *arenaOpt,
    PRUint32 numBuckets);

/* 3.4 Certificate cache routines */

NSS_EXTERN PRStatus
nssTrustDomain_InitializeCache(
    NSSTrustDomain *td,
    PRUint32 cacheSize);

NSS_EXTERN PRStatus
nssTrustDomain_AddCertsToCache(
    NSSTrustDomain *td,
    NSSCertificate **certs,
    PRUint32 numCerts);

NSS_EXTERN void
nssTrustDomain_RemoveCertFromCacheLOCKED(
    NSSTrustDomain *td,
    NSSCertificate *cert);

NSS_EXTERN void
nssTrustDomain_LockCertCache(NSSTrustDomain *td);

NSS_EXTERN void
nssTrustDomain_UnlockCertCache(NSSTrustDomain *td);

NSS_IMPLEMENT PRStatus
nssTrustDomain_DestroyCache(NSSTrustDomain *td);

/*
 * Remove all certs for the given token from the cache.  This is
 * needed if the token is removed.
 */
NSS_EXTERN PRStatus
nssTrustDomain_RemoveTokenCertsFromCache(
    NSSTrustDomain *td,
    NSSToken *token);

NSS_EXTERN PRStatus
nssTrustDomain_UpdateCachedTokenCerts(
    NSSTrustDomain *td,
    NSSToken *token);

/*
 * Find all cached certs with this nickname (label).
 */
NSS_EXTERN NSSCertificate **
nssTrustDomain_GetCertsForNicknameFromCache(
    NSSTrustDomain *td,
    const NSSUTF8 *nickname,
    nssList *certListOpt);

/*
 * Find all cached certs with this email address.
 */
NSS_EXTERN NSSCertificate **
nssTrustDomain_GetCertsForEmailAddressFromCache(
    NSSTrustDomain *td,
    NSSASCII7 *email,
    nssList *certListOpt);

/*
 * Find all cached certs with this subject.
 */
NSS_EXTERN NSSCertificate **
nssTrustDomain_GetCertsForSubjectFromCache(
    NSSTrustDomain *td,
    NSSDER *subject,
    nssList *certListOpt);

/*
 * Look for a specific cert in the cache.
 */
NSS_EXTERN NSSCertificate *
nssTrustDomain_GetCertForIssuerAndSNFromCache(
    NSSTrustDomain *td,
    NSSDER *issuer,
    NSSDER *serialNum);

/*
 * Look for a specific cert in the cache.
 */
NSS_EXTERN NSSCertificate *
nssTrustDomain_GetCertByDERFromCache(
    NSSTrustDomain *td,
    NSSDER *der);

/* Get all certs from the cache */
/* XXX this is being included to make some old-style calls word, not to
 *     say we should keep it
 */
NSS_EXTERN NSSCertificate **
nssTrustDomain_GetCertsFromCache(
    NSSTrustDomain *td,
    nssList *certListOpt);

NSS_EXTERN void
nssTrustDomain_DumpCacheInfo(
    NSSTrustDomain *td,
    void (*cert_dump_iter)(const void *, void *, void *),
    void *arg);

NSS_EXTERN void
nssCertificateList_AddReferences(
    nssList *certList);

PR_END_EXTERN_C

#endif /* PKIM_H */
