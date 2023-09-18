/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_validate.h
 *
 * Header file for validateChain function
 *
 */

#ifndef _PKIX_VALIDATE_H
#define _PKIX_VALIDATE_H
#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

PKIX_Error *
pkix_CheckChain(
        PKIX_List *certs,
        PKIX_UInt32 numCerts,
        PKIX_TrustAnchor *anchor,
        PKIX_List *checkers,
        PKIX_RevocationChecker *revChecker,
        PKIX_List *buildCheckedExtOIDs,
        PKIX_ProcessingParams *procParams,
        PKIX_UInt32 *pCertCheckedIndex,
        PKIX_UInt32 *pCheckerIndex,
        PKIX_Boolean *pRevChecking,
        PKIX_UInt32 *pReasonCode,
        void **pNBIOContext,
        PKIX_PL_PublicKey **pFinalSubjPubKey,
        PKIX_PolicyNode **pPolicyTree,
        PKIX_VerifyNode **pVerifyTree,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_VALIDATE_H */
