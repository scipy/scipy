/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_httpcertstore.h
 *
 * HTTPCertstore Object Type Definition
 *
 */

#ifndef _PKIX_PL_HTTPCERTSTORE_H
#define _PKIX_PL_HTTPCERTSTORE_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_HttpCertStoreContextStruct {
        const SEC_HttpClientFcn *client;
        SEC_HTTP_SERVER_SESSION serverSession;
        SEC_HTTP_REQUEST_SESSION requestSession;
	char *path;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_HttpCertStoreContext_RegisterSelf(void *plContext);

void pkix_pl_HttpCertStore_Shutdown(void *plContext);

PKIX_Error *
pkix_pl_HttpCertStore_CreateWithAsciiName(
        PKIX_PL_HttpClient *client,
	char *locationAscii,
        PKIX_CertStore **pCertStore,
        void *plContext);

PKIX_Error *
pkix_HttpCertStore_FindSocketConnection(
        PRIntervalTime timeout,
        char *hostname,
        PRUint16 portnum,
        PRErrorCode *pStatus,
        PKIX_PL_Socket **pSocket,
        void *plContext);

PKIX_Error *
pkix_pl_HttpCertStore_ProcessCertResponse(
	PRUint16 responseCode,
	const char *responseContentType,
	const char *responseData,
        PRUint32 responseDataLen,
	PKIX_List **pCertList,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_HTTPCERTSTORE_H */
