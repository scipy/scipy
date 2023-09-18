/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_httpdefaultclient.h
 *
 * HTTPDefaultClient Object Type Definition
 *
 */

#ifndef _PKIX_PL_HTTPDEFAULTCLIENT_H
#define _PKIX_PL_HTTPDEFAULTCLIENT_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define HTTP_DATA_BUFSIZE 4096
#define HTTP_HEADER_BUFSIZE 1024
#define HTTP_MIN_AVAILABLE_BUFFER_SIZE 512

typedef enum {
        HTTP_NOT_CONNECTED,
        HTTP_CONNECT_PENDING,
        HTTP_CONNECTED,
        HTTP_SEND_PENDING,
        HTTP_RECV_HDR,
        HTTP_RECV_HDR_PENDING,
        HTTP_RECV_BODY,
        HTTP_RECV_BODY_PENDING,
        HTTP_COMPLETE,
        HTTP_ERROR
} HttpConnectStatus;

typedef enum {
	HTTP_POST_METHOD,
	HTTP_GET_METHOD
} HttpMethod;

struct PKIX_PL_HttpDefaultClientStruct {
        HttpConnectStatus connectStatus;
        PRUint16 portnum;
        PRIntervalTime timeout;
        PKIX_UInt32 bytesToWrite;
        PKIX_UInt32 send_http_data_len;
        PKIX_UInt32 rcv_http_data_len;
        PKIX_UInt32 capacity;
        PKIX_UInt32 filledupBytes;
        PKIX_UInt32 responseCode;
        PKIX_UInt32 maxResponseLen;
        PKIX_UInt32 GETLen;
        PKIX_UInt32 POSTLen;
        PRUint32 *pRcv_http_data_len;
        PRPollDesc pollDesc;
        void *callbackList; /* cast this to (PKIX_PL_Socket_Callback *) */
        char *GETBuf;
        char *POSTBuf;
        char *rcvBuf;
        char *host;
        char *path;
        char *rcvContentType;
        void *rcvHeaders;
        HttpMethod send_http_method;
        const char *send_http_content_type;
        const char *send_http_data;
        PRUint16 *rcv_http_response_code;
        const char **rcv_http_content_type;
        const char **rcv_http_headers;
        const char **rcv_http_data;
        PKIX_PL_Socket *socket;
        void *plContext;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_HttpDefaultClient_RegisterSelf(void *plContext);

SECStatus
pkix_pl_HttpDefaultClient_CreateSessionFcn(
        const char *host,
        PRUint16 portnum,
        SEC_HTTP_SERVER_SESSION *pSession);

SECStatus
pkix_pl_HttpDefaultClient_KeepAliveSessionFcn(
        SEC_HTTP_SERVER_SESSION session,
        PRPollDesc **pPollDesc);

SECStatus
pkix_pl_HttpDefaultClient_FreeSessionFcn(
        SEC_HTTP_SERVER_SESSION session);

SECStatus
pkix_pl_HttpDefaultClient_RequestCreateFcn(
        SEC_HTTP_SERVER_SESSION session,
        const char *http_protocol_variant, /* usually "http" */
        const char *path_and_query_string,
        const char *http_request_method, 
        const PRIntervalTime timeout, 
        SEC_HTTP_REQUEST_SESSION *pRequest);

SECStatus
pkix_pl_HttpDefaultClient_SetPostDataFcn(
        SEC_HTTP_REQUEST_SESSION request,
        const char *http_data, 
        const PRUint32 http_data_len,
        const char *http_content_type);

SECStatus
pkix_pl_HttpDefaultClient_AddHeaderFcn(
        SEC_HTTP_REQUEST_SESSION request,
        const char *http_header_name, 
        const char *http_header_value);

SECStatus
pkix_pl_HttpDefaultClient_TrySendAndReceiveFcn(
        SEC_HTTP_REQUEST_SESSION request,
        PRPollDesc **pPollDesc,
        PRUint16 *http_response_code, 
        const char **http_response_content_type, 
        const char **http_response_headers, 
        const char **http_response_data, 
        PRUint32 *http_response_data_len); 

SECStatus
pkix_pl_HttpDefaultClient_CancelFcn(
        SEC_HTTP_REQUEST_SESSION request);

SECStatus
pkix_pl_HttpDefaultClient_FreeFcn(
        SEC_HTTP_REQUEST_SESSION request);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_HTTPDEFAULTCLIENT_H */
