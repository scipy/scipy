/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_ocspcertid.h
 *
 * Public Key Object Definitions
 *
 */

#ifndef _PKIX_PL_OCSPCERTID_H
#define _PKIX_PL_OCSPCERTID_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_OcspCertIDStruct {
        CERTOCSPCertID *certID;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_OcspCertID_RegisterSelf(void *plContext);

PKIX_Error *
PKIX_PL_OcspCertID_Create(
        PKIX_PL_Cert *cert,
        PKIX_PL_Date *validity,
        PKIX_PL_OcspCertID **object,
        void *plContext);

PKIX_Error *
PKIX_PL_OcspCertID_GetFreshCacheStatus(
        PKIX_PL_OcspCertID *cid, 
        PKIX_PL_Date *validity,
        PKIX_Boolean *hasFreshStatus,
        PKIX_Boolean *statusIsGood,
        SECErrorCodes *missingResponseError,
        void *plContext);

PKIX_Error *
PKIX_PL_OcspCertID_RememberOCSPProcessingFailure(
        PKIX_PL_OcspCertID *cid, 
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_OCSPCERTID_H */
