/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PKINSS3HACK_H
#define PKINSS3HACK_H

#ifndef NSSDEVT_H
#include "nssdevt.h"
#endif /* NSSDEVT_H */

#ifndef DEVT_H
#include "devt.h"
#endif /* DEVT_H */

#ifndef NSSPKIT_H
#include "nsspkit.h"
#endif /* NSSPKIT_H */

#include "base.h"

#include "cert.h"

PR_BEGIN_EXTERN_C

#define NSSITEM_FROM_SECITEM(nssit, secit) \
    (nssit)->data = (void *)(secit)->data; \
    (nssit)->size = (PRUint32)(secit)->len;

#define SECITEM_FROM_NSSITEM(secit, nssit)          \
    (secit)->data = (unsigned char *)(nssit)->data; \
    (secit)->len = (unsigned int)(nssit)->size;

NSS_EXTERN NSSTrustDomain *
STAN_GetDefaultTrustDomain();

NSS_EXTERN NSSCryptoContext *
STAN_GetDefaultCryptoContext();

NSS_EXTERN PRStatus
STAN_InitTokenForSlotInfo(NSSTrustDomain *td, PK11SlotInfo *slot);

NSS_EXTERN PRStatus
STAN_ResetTokenInterator(NSSTrustDomain *td);

NSS_EXTERN PRStatus
STAN_LoadDefaultNSS3TrustDomain(void);

NSS_EXTERN PRStatus
STAN_Shutdown();

NSS_EXTERN SECStatus
STAN_AddModuleToDefaultTrustDomain(SECMODModule *module);

NSS_EXTERN SECStatus
STAN_RemoveModuleFromDefaultTrustDomain(SECMODModule *module);

NSS_EXTERN CERTCertificate *
STAN_ForceCERTCertificateUpdate(NSSCertificate *c);

NSS_EXTERN CERTCertificate *
STAN_GetCERTCertificate(NSSCertificate *c);

NSS_EXTERN CERTCertificate *
STAN_GetCERTCertificateOrRelease(NSSCertificate *c);

NSS_EXTERN NSSCertificate *
STAN_GetNSSCertificate(CERTCertificate *c);

NSS_EXTERN CERTCertTrust *
nssTrust_GetCERTCertTrustForCert(NSSCertificate *c, CERTCertificate *cc);

NSS_EXTERN PRStatus
STAN_DeleteCertTrustMatchingSlot(NSSCertificate *c);

NSS_EXTERN PRStatus
STAN_ChangeCertTrust(CERTCertificate *cc, CERTCertTrust *trust);

NSS_EXTERN PRStatus
nssPKIX509_GetIssuerAndSerialFromDER(NSSDER *der,
                                     NSSDER *issuer, NSSDER *serial);

NSS_EXTERN char *
STAN_GetCERTCertificateName(PLArenaPool *arenaOpt, NSSCertificate *c);

NSS_EXTERN char *
STAN_GetCERTCertificateNameForInstance(PLArenaPool *arenaOpt,
                                       NSSCertificate *c,
                                       nssCryptokiInstance *instance);

/* exposing this */
NSS_EXTERN NSSCertificate *
NSSCertificate_Create(NSSArena *arenaOpt);

/* This function is being put here because it is a hack for
 * PK11_FindCertFromNickname.
 */
NSS_EXTERN NSSCertificate *
nssTrustDomain_FindBestCertificateByNicknameForToken(
    NSSTrustDomain *td,
    NSSToken *token,
    NSSUTF8 *name,
    NSSTime *timeOpt, /* NULL for "now" */
    NSSUsage *usage,
    NSSPolicies *policiesOpt /* NULL for none */
);

/* This function is being put here because it is a hack for
 * PK11_FindCertsFromNickname.
 */
NSS_EXTERN NSSCertificate **
nssTrustDomain_FindCertificatesByNicknameForToken(
    NSSTrustDomain *td,
    NSSToken *token,
    NSSUTF8 *name,
    NSSCertificate *rvOpt[],
    PRUint32 maximumOpt, /* 0 for no max */
    NSSArena *arenaOpt);

/* CERT_TraversePermCertsForSubject */
NSS_EXTERN PRStatus
nssTrustDomain_TraverseCertificatesBySubject(
    NSSTrustDomain *td,
    NSSDER *subject,
    PRStatus (*callback)(NSSCertificate *c, void *arg),
    void *arg);

/* CERT_TraversePermCertsForNickname */
NSS_EXTERN PRStatus
nssTrustDomain_TraverseCertificatesByNickname(
    NSSTrustDomain *td,
    NSSUTF8 *nickname,
    PRStatus (*callback)(NSSCertificate *c, void *arg),
    void *arg);

/* SEC_TraversePermCerts */
NSS_EXTERN PRStatus
nssTrustDomain_TraverseCertificates(
    NSSTrustDomain *td,
    PRStatus (*callback)(NSSCertificate *c, void *arg),
    void *arg);

/* CERT_AddTempCertToPerm */
NSS_EXTERN PRStatus
nssTrustDomain_AddTempCertToPerm(NSSCertificate *c);

PR_END_EXTERN_C

#endif /* PKINSS3HACK_H */
