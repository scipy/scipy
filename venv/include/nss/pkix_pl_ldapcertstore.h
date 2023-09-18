/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_ldapcertstore.h
 *
 * LDAPCertstore Object Type Definition
 *
 */

#ifndef _PKIX_PL_LDAPCERTSTORE_H
#define _PKIX_PL_LDAPCERTSTORE_H

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
 * will have to change.) Therefore the CertStore_Create API contains a BindAPI
 * structure, a union, which will have to be revised and extended when this
 * area of the protocol is better understood.
 *
 * It is further assumed that a given LdapCertStore will connect only to a
 * single server, and that the creation of the socket will initiate the
 * CONNECT. Therefore the LdapCertStore handles only the case of continuing
 * the connection, if nonblocking I/O is being used.
 */

typedef enum {
        LDAP_CONNECT_PENDING,
        LDAP_CONNECTED,
        LDAP_BIND_PENDING,
        LDAP_BIND_RESPONSE,
        LDAP_BIND_RESPONSE_PENDING,
        LDAP_BOUND,
        LDAP_SEND_PENDING,
        LDAP_RECV,
        LDAP_RECV_PENDING,
        LDAP_RECV_INITIAL,
        LDAP_RECV_NONINITIAL,
        LDAP_ABANDON_PENDING
} LDAPConnectStatus;

#define LDAP_CACHEBUCKETS       128
#define RCVBUFSIZE              512

struct PKIX_PL_LdapCertStoreContext {
        PKIX_PL_LdapClient *client;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_LdapCertStoreContext_RegisterSelf(void *plContext);

PKIX_Error *
pkix_pl_LdapCertStore_BuildCertList(
        PKIX_List *responseList,
        PKIX_List **pCerts,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_LDAPCERTSTORE_H */
