/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PKISTORE_H
#define PKISTORE_H

#ifndef NSSPKIT_H
#include "nsspkit.h"
#endif /* NSSPKIT_H */

#ifndef BASE_H
#include "base.h"
#endif /* BASE_H */

PR_BEGIN_EXTERN_C

/*
 * PKI Stores
 *
 * This is a set of routines for managing local stores of PKI objects.
 * Currently, the only application is in crypto contexts, where the
 * certificate store is used.  In the future, methods should be added
 * here for storing local references to keys.
 */

/*
 * nssCertificateStore
 *
 * Manages local store of certificate, trust, and S/MIME profile objects.
 * Within a crypto context, mappings of cert to trust and cert to S/MIME
 * profile are always 1-1.  Therefore, it is reasonable to store all objects
 * in a single collection, indexed by the certificate.
 */

NSS_EXTERN nssCertificateStore *
nssCertificateStore_Create(
    NSSArena *arenaOpt);

NSS_EXTERN PRStatus
nssCertificateStore_Destroy(
    nssCertificateStore *store);

/* Atomic Find cert in store, or add this cert to the store.
** Ref counts properly maintained.
*/
NSS_EXTERN NSSCertificate *
nssCertificateStore_FindOrAdd(
    nssCertificateStore *store,
    NSSCertificate *c);

NSS_EXTERN void
nssCertificateStore_RemoveCertLOCKED(
    nssCertificateStore *store,
    NSSCertificate *cert);

struct nssCertificateStoreTraceStr {
    nssCertificateStore *store;
    PZLock *lock;
    PRBool locked;
    PRBool unlocked;
};

typedef struct nssCertificateStoreTraceStr nssCertificateStoreTrace;

NSS_EXTERN void
nssCertificateStore_Lock(
    nssCertificateStore *store, nssCertificateStoreTrace *out);

NSS_EXTERN void
nssCertificateStore_Unlock(
    nssCertificateStore *store, const nssCertificateStoreTrace *in,
    nssCertificateStoreTrace *out);

NSS_EXTERN NSSCertificate **
nssCertificateStore_FindCertificatesBySubject(
    nssCertificateStore *store,
    NSSDER *subject,
    NSSCertificate *rvOpt[],
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSCertificate **
nssCertificateStore_FindCertificatesByNickname(
    nssCertificateStore *store,
    const NSSUTF8 *nickname,
    NSSCertificate *rvOpt[],
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSCertificate **
nssCertificateStore_FindCertificatesByEmail(
    nssCertificateStore *store,
    NSSASCII7 *email,
    NSSCertificate *rvOpt[],
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSCertificate *
nssCertificateStore_FindCertificateByIssuerAndSerialNumber(
    nssCertificateStore *store,
    NSSDER *issuer,
    NSSDER *serial);

NSS_EXTERN NSSCertificate *
nssCertificateStore_FindCertificateByEncodedCertificate(
    nssCertificateStore *store,
    NSSDER *encoding);

NSS_EXTERN PRStatus
nssCertificateStore_AddTrust(
    nssCertificateStore *store,
    NSSTrust *trust);

NSS_EXTERN NSSTrust *
nssCertificateStore_FindTrustForCertificate(
    nssCertificateStore *store,
    NSSCertificate *cert);

NSS_EXTERN PRStatus
nssCertificateStore_AddSMIMEProfile(
    nssCertificateStore *store,
    nssSMIMEProfile *profile);

NSS_EXTERN nssSMIMEProfile *
nssCertificateStore_FindSMIMEProfileForCertificate(
    nssCertificateStore *store,
    NSSCertificate *cert);

NSS_EXTERN void
nssCertificateStore_DumpStoreInfo(
    nssCertificateStore *store,
    void (*cert_dump_iter)(const void *, void *, void *),
    void *arg);

PR_END_EXTERN_C

#endif /* PKISTORE_H */
