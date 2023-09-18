/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_certchainchecker.h
 *
 * CertChainChecker Object Type Definition
 *
 */

#ifndef _PKIX_CERTCHAINCHECKER_H
#define _PKIX_CERTCHAINCHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_CertChainCheckerStruct {
        PKIX_CertChainChecker_CheckCallback checkCallback;
        PKIX_List *extensions;
        PKIX_PL_Object *state;
        PKIX_Boolean forwardChecking;
        PKIX_Boolean isForwardDirectionExpected;
};

/* see source file for function documentation */

PKIX_Error *pkix_CertChainChecker_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CERTCHAINCHECKER_H */
