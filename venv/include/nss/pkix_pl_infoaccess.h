/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_infoaccess.h
 *
 * InfoAccess Object Definitions
 *
 */

#ifndef _PKIX_PL_INFOACCESS_H
#define _PKIX_PL_INFOACCESS_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_InfoAccessStruct{
        PKIX_UInt32 method;
        PKIX_PL_GeneralName *location;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_InfoAccess_RegisterSelf(void *plContext);

PKIX_Error *
pkix_pl_InfoAccess_CreateList(
        CERTAuthInfoAccess **authInfoAccess,
        PKIX_List **pAiaList, /* of PKIX_PL_InfoAccess */
        void *plContext);

#ifndef NSS_PKIX_NO_LDAP
PKIX_Error *
pkix_pl_InfoAccess_ParseLocation(
        PKIX_PL_GeneralName *generalName,
        PLArenaPool *arena,
        LDAPRequestParams *request,
        char **pDomainName,
        void *plContext);
#endif /* !NSS_PKIX_NO_LDAP */

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_INFOACCESS_H */
