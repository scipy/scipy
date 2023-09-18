/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_defaultcrlchecker.h
 *
 * Header file for default CRL function
 *
 */

#ifndef _PKIX_CRLCHECKER_H
#define _PKIX_CRLCHECKER_H

#include "pkixt.h"
#include "pkix_revocationmethod.h"
#include "pkix_crlsel.h"

#ifdef __cplusplus
extern "C" {
#endif

/* NOTE: nbio logic removed. Will be replaced later. */

PKIX_Error *
pkix_CrlChecker_CheckLocal(
        PKIX_PL_Cert *cert,
        PKIX_PL_Cert *issuer,
        PKIX_PL_Date *date,
        pkix_RevocationMethod *checkerObject,
        PKIX_ProcessingParams *procParams,
        PKIX_UInt32 methodFlags,
        PKIX_Boolean chainVerificationState,
        PKIX_RevocationStatus *pRevStatus,
        CERTCRLEntryReasonCode *reasonCode,
        void *plContext);

PKIX_Error *
pkix_CrlChecker_CheckExternal(
        PKIX_PL_Cert *cert,
        PKIX_PL_Cert *issuer,
        PKIX_PL_Date *date,
        pkix_RevocationMethod *checkerObject,
        PKIX_ProcessingParams *procParams,
        PKIX_UInt32 methodFlags,
        PKIX_RevocationStatus *pRevStatus,
        CERTCRLEntryReasonCode *reasonCode,
        void **pNBIOContext,
        void *plContext);

PKIX_Error *
pkix_CrlChecker_Create(PKIX_RevocationMethodType methodType,
                       PKIX_UInt32 flags,
                       PKIX_UInt32 priority,
                       pkix_LocalRevocationCheckFn localRevChecker,
                       pkix_ExternalRevocationCheckFn externalRevChecker,
                       PKIX_List *certStores,
                       PKIX_PL_VerifyCallback crlVerifyFn,
                       pkix_RevocationMethod **pChecker,
                       void *plContext);

PKIX_Error *
pkix_CrlChecker_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CRLCHECKER_H */
