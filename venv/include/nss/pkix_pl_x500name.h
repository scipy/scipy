/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_x500name.h
 *
 * X500Name Object Type Definitions
 *
 */

#ifndef _PKIX_PL_X500NAME_H
#define _PKIX_PL_X500NAME_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif


struct PKIX_PL_X500NameStruct{
        PLArenaPool *arena; /* X500Name arena. Shared arena with nssDN
                             * and derName */
        CERTName nssDN;
        SECItem derName;    /* adding DER encoded CERTName to the structure
                             * to avoid unnecessary name encoding when pass
                             * der name to cert finder */
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_X500Name_RegisterSelf(void *plContext);

PKIX_Error *pkix_pl_X500Name_GetDERName(
        PKIX_PL_X500Name *xname,
        PLArenaPool *arena,
        SECItem **pSECName,
        void *plContext);

#ifdef BUILD_LIBPKIX_TESTS
PKIX_Error * pkix_pl_X500Name_CreateFromUtf8(
        char *stringRep,
        PKIX_PL_X500Name **pName,
        void *plContext);
#endif /* BUILD_LIBPKIX_TESTS */

PKIX_Error *pkix_pl_X500Name_GetCommonName(
        PKIX_PL_X500Name *xname,
        unsigned char **pCommonName,
        void *plContext);

PKIX_Error *
pkix_pl_X500Name_GetCountryName(
        PKIX_PL_X500Name *xname,
        unsigned char **pCountryName,
        void *plContext);

PKIX_Error *
pkix_pl_X500Name_GetOrgName(
        PKIX_PL_X500Name *xname,
        unsigned char **pOrgName,
        void *plContext);

PKIX_Error *
pkix_pl_X500Name_GetCERTName(
        PKIX_PL_X500Name *xname,
        CERTName **pCERTName,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_X500NAME_H */
