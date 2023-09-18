/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_valparams.h
 *
 * ValidateParams Object Type Definition
 *
 */

#ifndef _PKIX_VALIDATEPARAMS_H
#define _PKIX_VALIDATEPARAMS_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_ValidateParamsStruct {
        PKIX_ProcessingParams *procParams;      /* Never NULL */
        PKIX_List *chain;                       /* Never NULL */
};

/* see source file for function documentation */

PKIX_Error *pkix_ValidateParams_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_VALIDATEPARAMS_H */
