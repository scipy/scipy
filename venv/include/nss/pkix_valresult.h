/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_valresult.h
 *
 * ValidateResult Object Type Definition
 *
 */

#ifndef _PKIX_VALIDATERESULT_H
#define _PKIX_VALIDATERESULT_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_ValidateResultStruct {
        PKIX_PL_PublicKey *pubKey;
        PKIX_TrustAnchor *anchor;
        PKIX_PolicyNode *policyTree;
};

/* see source file for function documentation */

PKIX_Error *
pkix_ValidateResult_Create(
        PKIX_PL_PublicKey *pubKey,
        PKIX_TrustAnchor *anchor,
        PKIX_PolicyNode *policyTree,
        PKIX_ValidateResult **pResult,
        void *plContext);

PKIX_Error *pkix_ValidateResult_RegisterSelf(void *plContext);


#ifdef __cplusplus
}
#endif

#endif /* _PKIX_VALIDATERESULT_H */
