/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_ldapresponse.h
 *
 * LdapResponse Object Definitions
 *
 */

#ifndef _PKIX_PL_LDAPRESPONSE_H
#define _PKIX_PL_LDAPRESPONSE_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_LdapResponseStruct{
        LDAPMessage decoded;
        PKIX_UInt32 partialLength;
        PKIX_UInt32 totalLength;
        SECItem derEncoded;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_LdapResponse_Create(
        LDAPMessageType responseType,
        PKIX_UInt32 totalLength,
        PKIX_UInt32 bytesAvailable,
        void *partialData,
        PKIX_UInt32 *pBytesConsumed,
        PKIX_PL_LdapResponse **pResponse,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_Append(
        PKIX_PL_LdapResponse *response,
        PKIX_UInt32 partialLength,
        void *partialData,
        PKIX_UInt32 *bytesConsumed,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_IsComplete(
        PKIX_PL_LdapResponse *response,
        PKIX_Boolean *pIsComplete,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_Decode(
        PLArenaPool *arena,
        PKIX_PL_LdapResponse *response,
        SECStatus *pStatus,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_GetMessage(
        PKIX_PL_LdapResponse *response,
        LDAPMessage **pMessage,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_GetMessageType(
        PKIX_PL_LdapResponse *response,
        LDAPMessageType *pMessageType,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_GetCapacity(
        PKIX_PL_LdapResponse *response,
        PKIX_UInt32 *pCapacity,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_GetResultCode(
        PKIX_PL_LdapResponse *response,
        LDAPResultCode *pResultCode,
        void *plContext);

PKIX_Error *
pkix_pl_LdapResponse_GetAttributes(
        PKIX_PL_LdapResponse *response,
        LDAPSearchResponseAttr ***pAttributes,
        void *plContext);

PKIX_Error *pkix_pl_LdapResponse_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_LDAPRESPONSE_H */
