/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PKI_H
#define PKI_H

#ifndef NSSDEVT_H
#include "nssdevt.h"
#endif /* NSSDEVT_H */

#ifndef NSSPKI_H
#include "nsspki.h"
#endif /* NSSPKI_H */

#ifndef PKIT_H
#include "pkit.h"
#endif /* PKIT_H */

PR_BEGIN_EXTERN_C

NSS_EXTERN NSSCallback *
nssTrustDomain_GetDefaultCallback(
    NSSTrustDomain *td,
    PRStatus *statusOpt);

NSS_EXTERN NSSCertificate **
nssTrustDomain_FindCertificatesBySubject(
    NSSTrustDomain *td,
    NSSDER *subject,
    NSSCertificate *rvOpt[],
    PRUint32 maximumOpt,
    NSSArena *arenaOpt);

NSS_EXTERN NSSTrust *
nssTrustDomain_FindTrustForCertificate(
    NSSTrustDomain *td,
    NSSCertificate *c);

NSS_EXTERN NSSCertificate *
nssCertificate_AddRef(NSSCertificate *c);

NSS_EXTERN PRStatus
nssCertificate_Destroy(NSSCertificate *c);

NSS_EXTERN NSSDER *
nssCertificate_GetEncoding(NSSCertificate *c);

NSS_EXTERN NSSDER *
nssCertificate_GetIssuer(NSSCertificate *c);

NSS_EXTERN NSSDER *
nssCertificate_GetSerialNumber(NSSCertificate *c);

NSS_EXTERN NSSDER *
nssCertificate_GetSubject(NSSCertificate *c);

/* Returns a copy, Caller must free using nss_ZFreeIf */
NSS_EXTERN NSSUTF8 *
nssCertificate_GetNickname(
    NSSCertificate *c,
    NSSToken *tokenOpt);

NSS_EXTERN NSSASCII7 *
nssCertificate_GetEmailAddress(NSSCertificate *c);

NSS_EXTERN PRBool
nssCertificate_IssuerAndSerialEqual(
    NSSCertificate *c1,
    NSSCertificate *c2);

NSS_EXTERN NSSPrivateKey *
nssPrivateKey_AddRef(NSSPrivateKey *vk);

NSS_EXTERN PRStatus
nssPrivateKey_Destroy(NSSPrivateKey *vk);

NSS_EXTERN NSSItem *
nssPrivateKey_GetID(NSSPrivateKey *vk);

NSS_EXTERN NSSUTF8 *
nssPrivateKey_GetNickname(
    NSSPrivateKey *vk,
    NSSToken *tokenOpt);

NSS_EXTERN PRStatus
nssPublicKey_Destroy(NSSPublicKey *bk);

NSS_EXTERN NSSItem *
nssPublicKey_GetID(NSSPublicKey *vk);

NSS_EXTERN NSSCertificate **
nssCryptoContext_FindCertificatesBySubject(
    NSSCryptoContext *cc,
    NSSDER *subject,
    NSSCertificate *rvOpt[],
    PRUint32 maximumOpt, /* 0 for no max */
    NSSArena *arenaOpt);

/* putting here for now, needs more thought */
NSS_EXTERN PRStatus
nssCryptoContext_ImportTrust(
    NSSCryptoContext *cc,
    NSSTrust *trust);

NSS_EXTERN NSSTrust *
nssCryptoContext_FindTrustForCertificate(
    NSSCryptoContext *cc,
    NSSCertificate *cert);

NSS_EXTERN PRStatus
nssCryptoContext_ImportSMIMEProfile(
    NSSCryptoContext *cc,
    nssSMIMEProfile *profile);

NSS_EXTERN nssSMIMEProfile *
nssCryptoContext_FindSMIMEProfileForCertificate(
    NSSCryptoContext *cc,
    NSSCertificate *cert);

NSS_EXTERN NSSTrust *
nssTrust_AddRef(NSSTrust *trust);

NSS_EXTERN PRStatus
nssTrust_Destroy(NSSTrust *trust);

NSS_EXTERN nssSMIMEProfile *
nssSMIMEProfile_AddRef(nssSMIMEProfile *profile);

NSS_EXTERN PRStatus
nssSMIMEProfile_Destroy(nssSMIMEProfile *profile);

NSS_EXTERN nssSMIMEProfile *
nssSMIMEProfile_Create(
    NSSCertificate *cert,
    NSSItem *profileTime,
    NSSItem *profileData);

PR_END_EXTERN_C

#endif /* PKI_H */
