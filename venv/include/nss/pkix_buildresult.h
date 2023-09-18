/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_buildresult.h
 *
 * BuildResult Object Type Definition
 *
 */

#ifndef _PKIX_BUILDRESULT_H
#define _PKIX_BUILDRESULT_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_BuildResultStruct {
        PKIX_ValidateResult *valResult;
        PKIX_List *certChain;
};

/* see source file for function documentation */

PKIX_Error *
pkix_BuildResult_Create(
        PKIX_ValidateResult *valResult,
        PKIX_List *certChain,
        PKIX_BuildResult **pResult,
        void *plContext);

PKIX_Error *pkix_BuildResult_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_BUILDRESULT_H */
