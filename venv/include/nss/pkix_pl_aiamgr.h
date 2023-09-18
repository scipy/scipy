/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_aiamgr.h
 *
 * AIAMgr Object Definitions
 *
 */

#ifndef _PKIX_PL_AIAMGR_H
#define _PKIX_PL_AIAMGR_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_AIAMgrStruct {
        /* pointer to cert cache */
        /* pointer to crl cache */
        PKIX_UInt32 method;
        PKIX_UInt32 aiaIndex;
        PKIX_UInt32 numAias;
        PKIX_List *aia;
        PKIX_PL_GeneralName *location;
        PKIX_List *results;
	union {
#ifndef NSS_PKIX_NO_LDAP
	        PKIX_PL_LdapClient *ldapClient;
#endif
		struct {
		        const SEC_HttpClientFcn *httpClient;
			SEC_HTTP_SERVER_SESSION serverSession;
			SEC_HTTP_REQUEST_SESSION requestSession;
			char *path;
		} hdata;
	} client;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_AIAMgr_RegisterSelf(void *plContext);

#ifndef NSS_PKIX_NO_LDAP
PKIX_Error *PKIX_PL_LdapClient_InitiateRequest(
        PKIX_PL_LdapClient *client,
        LDAPRequestParams *requestParams,
        void **pPollDesc,
        PKIX_List **pResponse,
        void *plContext);

PKIX_Error *PKIX_PL_LdapClient_ResumeRequest(
        PKIX_PL_LdapClient *client,
        void **pPollDesc,
        PKIX_List **pResponse,
        void *plContext);
#endif /* !NSS_PKIX_NO_LDAP */

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_AIAMGR_H */
