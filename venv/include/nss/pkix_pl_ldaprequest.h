/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_ldaprequest.h
 *
 * LdapRequest Object Definitions
 *
 */

#ifndef _PKIX_PL_LDAPREQUEST_H
#define _PKIX_PL_LDAPREQUEST_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
        USER_CERT,
        CA_CERT,
        CROSS_CERT,
        CRL,
        ARL,
        DELTA_CRL
} PKIX_PL_LdapAttr;

struct PKIX_PL_LdapRequestStruct{
        PLArenaPool *arena;
        PKIX_UInt32 msgnum;
        char *issuerDN;
        ScopeType scope;
        DerefType derefAliases;
        PKIX_UInt32 sizeLimit;
        PKIX_UInt32 timeLimit;
        char attrsOnly;
        LDAPFilter *filter;
        LdapAttrMask attrBits;
        SECItem attributes[MAX_LDAPATTRS];
        SECItem **attrArray;
        SECItem *encoded;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_LdapRequest_Create(
        PLArenaPool *arena,
        PKIX_UInt32 msgnum,
        char *issuerDN,
        ScopeType scope,
        DerefType derefAliases,
        PKIX_UInt32 sizeLimit,
        PKIX_UInt32 timeLimit,
        char attrsOnly,
        LDAPFilter *filter,
        LdapAttrMask attrBits,
        PKIX_PL_LdapRequest **pRequestMsg,
        void *plContext);

PKIX_Error *
pkix_pl_LdapRequest_AttrTypeToBit(
        SECItem *attrType,
        LdapAttrMask *pAttrBit,
        void *plContext);

PKIX_Error *
pkix_pl_LdapRequest_AttrStringToBit(
        char *attrString,
        LdapAttrMask *pAttrBit,
        void *plContext);

PKIX_Error *
pkix_pl_LdapRequest_GetEncoded(
        PKIX_PL_LdapRequest *request,
        SECItem **pRequestBuf,
        void *plContext);

PKIX_Error *pkix_pl_LdapRequest_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_LDAPREQUEST_H */
