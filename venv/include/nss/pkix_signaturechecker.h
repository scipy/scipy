/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_signaturechecker.h
 *
 * Header file for validate signature function
 *
 */

#ifndef _PKIX_SIGNATURECHECKER_H
#define _PKIX_SIGNATURECHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pkix_SignatureCheckerState pkix_SignatureCheckerState;

struct pkix_SignatureCheckerState {
        PKIX_Boolean prevCertCertSign;
        PKIX_UInt32 certsRemaining;
        PKIX_PL_PublicKey *prevPublicKey; /* Subject PubKey of last cert */
        PKIX_List *prevPublicKeyList; /* of PKIX_PL_PublicKey */
        PKIX_PL_OID *keyUsageOID;
};

PKIX_Error *
pkix_SignatureChecker_Initialize(
        PKIX_PL_PublicKey *trustedPubKey,
        PKIX_UInt32 certsRemaining,
        PKIX_CertChainChecker **pChecker,
        void *plContext);

PKIX_Error *
pkix_SignatureCheckerState_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_SIGNATURECHECKER_H */
