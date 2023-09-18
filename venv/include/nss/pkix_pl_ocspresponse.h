/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_ocspresponse.h
 *
 * OcspResponse Object Definitions
 *
 */

#ifndef _PKIX_PL_OCSPRESPONSE_H
#define _PKIX_PL_OCSPRESPONSE_H

#include "pkix_pl_common.h"
#include "pkix_pl_ocspcertid.h"
#include "hasht.h"
#include "cryptohi.h"
#include "ocspti.h"
#include "ocspi.h"
#include "plbase64.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_OCSP_RESPONSE_LEN (64*1024)

struct PKIX_PL_OcspResponseStruct{
        PLArenaPool *arena;
        const PKIX_PL_OcspRequest *request;
        const SEC_HttpClientFcn *httpClient;
        SEC_HTTP_SERVER_SESSION serverSession;
        SEC_HTTP_REQUEST_SESSION sessionRequest;
        PKIX_PL_VerifyCallback verifyFcn;
        SECItem *encodedResponse;
        CERTCertDBHandle *handle;
        PRTime producedAt;
        PKIX_PL_Date *producedAtDate;
        PKIX_PL_Cert *pkixSignerCert;
        CERTOCSPResponse *nssOCSPResponse;
        CERTCertificate *signerCert;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_OcspResponse_RegisterSelf(void *plContext);

PKIX_Error *
pkix_pl_OcspResponse_Create(
        PKIX_PL_OcspRequest *request,
        const char *httpMechanism,
        void *responder,
        PKIX_PL_VerifyCallback verifyFcn,
        void **pNBIOContext,
        PKIX_PL_OcspResponse **pResponse,
        void *plContext);

PKIX_Error *
pkix_pl_OcspResponse_Decode(
        PKIX_PL_OcspResponse *response,
        PKIX_Boolean *passed,
        SECErrorCodes *pReturnCode,
        void *plContext);

PKIX_Error *
pkix_pl_OcspResponse_GetStatus(
        PKIX_PL_OcspResponse *response,
        PKIX_Boolean *passed,
        SECErrorCodes *pReturnCode,
        void *plContext);

PKIX_Error *
pkix_pl_OcspResponse_VerifySignature(
        PKIX_PL_OcspResponse *response,
        PKIX_PL_Cert *cert,
        PKIX_ProcessingParams *procParams,
        PKIX_Boolean *pPassed,
        void **pNBIOContext,
        void *plContext);

PKIX_Error *
pkix_pl_OcspResponse_GetStatusForCert(
        PKIX_PL_OcspCertID *cid,
        PKIX_PL_OcspResponse *response,
        PKIX_Boolean allowCachingOfFailures,
        PKIX_PL_Date *validity,
        PKIX_Boolean *pPassed,
        SECErrorCodes *pReturnCode,
        void *plContext);

PKIX_Error *
PKIX_PL_OcspResponse_UseBuildChain(
        PKIX_PL_Cert *signerCert,
	PKIX_PL_Date *producedAt,
        PKIX_ProcessingParams *procParams,
        void **pNBIOContext,
        void **pState,
        PKIX_BuildResult **pBuildResult,
        PKIX_VerifyNode **pVerifyTree,
	void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_OCSPRESPONSE_H */
