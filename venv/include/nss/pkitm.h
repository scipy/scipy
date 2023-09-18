/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef PKITM_H
#define PKITM_H

/*
 * pkitm.h
 *
 * This file contains PKI-module specific types.
 */

#ifndef BASET_H
#include "baset.h"
#endif /* BASET_H */

#ifndef PKIT_H
#include "pkit.h"
#endif /* PKIT_H */

PR_BEGIN_EXTERN_C

typedef enum nssCertIDMatchEnum {
    nssCertIDMatch_Yes = 0,
    nssCertIDMatch_No = 1,
    nssCertIDMatch_Unknown = 2
} nssCertIDMatch;

/*
 * nssDecodedCert
 *
 * This is an interface to allow the PKI module access to certificate
 * information that can only be found by decoding.  The interface is
 * generic, allowing each certificate type its own way of providing
 * the information
 */
struct nssDecodedCertStr {
    NSSCertificateType type;
    void *data;
    /* returns the unique identifier for the cert */
    NSSItem *(*getIdentifier)(nssDecodedCert *dc);
    /* returns the unique identifier for this cert's issuer */
    void *(*getIssuerIdentifier)(nssDecodedCert *dc);
    /* is id the identifier for this cert? */
    nssCertIDMatch (*matchIdentifier)(nssDecodedCert *dc, void *id);
    /* is this cert a valid CA cert? */
    PRBool (*isValidIssuer)(nssDecodedCert *dc);
    /* returns the cert usage */
    NSSUsage *(*getUsage)(nssDecodedCert *dc);
    /* is time within the validity period of the cert? */
    PRBool (*isValidAtTime)(nssDecodedCert *dc, NSSTime *time);
    /* is the validity period of this cert newer than cmpdc? */
    PRBool (*isNewerThan)(nssDecodedCert *dc, nssDecodedCert *cmpdc);
    /* does the usage for this cert match the requested usage? */
    PRBool (*matchUsage)(nssDecodedCert *dc, const NSSUsage *usage);
    /* is this cert trusted for the requested usage? */
    PRBool (*isTrustedForUsage)(nssDecodedCert *dc,
                                const NSSUsage *usage);
    /* extract the email address */
    NSSASCII7 *(*getEmailAddress)(nssDecodedCert *dc);
    /* extract the DER-encoded serial number */
    PRStatus (*getDERSerialNumber)(nssDecodedCert *dc,
                                   NSSDER *derSerial, NSSArena *arena);
};

struct NSSUsageStr {
    PRBool anyUsage;
    SECCertUsage nss3usage;
    PRBool nss3lookingForCA;
};

typedef struct nssPKIObjectCollectionStr nssPKIObjectCollection;

typedef struct
{
    union {
        PRStatus (*cert)(NSSCertificate *c, void *arg);
        PRStatus (*crl)(NSSCRL *crl, void *arg);
        PRStatus (*pvkey)(NSSPrivateKey *vk, void *arg);
        PRStatus (*pbkey)(NSSPublicKey *bk, void *arg);
    } func;
    void *arg;
} nssPKIObjectCallback;

PR_END_EXTERN_C

#endif /* PKITM_H */
