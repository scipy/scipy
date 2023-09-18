/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_ldapdefaultclient.h
 *
 * LDAPDefaultClient Object Type Definition
 *
 */

#ifndef _PKIX_PL_LDAPDEFAULTCLIENT_H
#define _PKIX_PL_LDAPDEFAULTCLIENT_H

#include "pkix_pl_ldapt.h"
#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * At the time of this version, there are unresolved questions about the LDAP
 * protocol. Although RFC1777 describes a BIND and UNBIND message, it is not
 * clear whether they are appropriate to this application. We have tested only
 * using servers that do not expect authentication, and that reject BIND
 * messages. It is not clear what values might be appropriate for the bindname
 * and authentication fields, which are currently implemented as char strings
 * supplied by the caller. (If this changes, the API and possibly the templates
 * will have to change.) Therefore the LDAPClient_Create API contains a
 * BindAPI structure, a union, which will have to be revised and extended when
 * this area of the protocol is better understood.
 *
 */

typedef enum {
        CONNECT_PENDING,
        CONNECTED,
        BIND_PENDING,
        BIND_RESPONSE,
        BIND_RESPONSE_PENDING,
        BOUND,
        SEND_PENDING,
        RECV,
        RECV_PENDING,
        RECV_INITIAL,
        RECV_NONINITIAL,
        ABANDON_PENDING
} LdapClientConnectStatus;

struct PKIX_PL_LdapDefaultClientStruct {
        PKIX_PL_LdapClient vtable;
        LdapClientConnectStatus connectStatus;
        PKIX_UInt32 messageID;
        PKIX_PL_HashTable *cachePtr;
        PKIX_PL_Socket *clientSocket;
        PRPollDesc pollDesc;
        void *callbackList; /* cast this to (PKIX_PL_Socket_Callback *) */
        LDAPBindAPI *bindAPI;
        PLArenaPool *arena;
        PRTime lastIO;
        void *sendBuf;
        PKIX_UInt32 bytesToWrite;
        void *rcvBuf;
        PKIX_UInt32 capacity;
        void *currentInPtr;
        PKIX_UInt32 currentBytesAvailable;
        void *bindMsg;
        PKIX_UInt32 bindMsgLen;
        PKIX_List *entriesFound;
        PKIX_PL_LdapRequest *currentRequest;
        PKIX_PL_LdapResponse *currentResponse;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_LdapDefaultClient_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_LDAPDEFAULTCLIENT_H */
