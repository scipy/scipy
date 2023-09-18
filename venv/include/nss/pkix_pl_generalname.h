/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_generalname.h
 *
 * GeneralName Object Definitions
 *
 */

#ifndef _PKIX_PL_GENERALNAME_H
#define _PKIX_PL_GENERALNAME_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_GeneralNameStruct{
        CERTGeneralNameList *nssGeneralNameList;
        CERTGeneralNameType type;
        PKIX_PL_X500Name *directoryName;
        PKIX_PL_OID *oid;
        OtherName *OthName;
        SECItem *other;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_GeneralName_Create(
        CERTGeneralName *nssAltName,
        PKIX_PL_GeneralName **pGenName,
        void *plContext);

PKIX_Error *
pkix_pl_GeneralName_GetNssGeneralName(
        PKIX_PL_GeneralName *genName,
        CERTGeneralName **pNssGenName,
        void *plContext);

PKIX_Error *pkix_pl_GeneralName_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_GENERALNAME_H */
