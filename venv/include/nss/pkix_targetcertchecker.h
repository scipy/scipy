/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_targetcertchecker.h
 *
 * Header file for validate target cert function
 *
 */

#ifndef _PKIX_TARGETCERTCHECKER_H
#define _PKIX_TARGETCERTCHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pkix_TargetCertCheckerState pkix_TargetCertCheckerState;

struct pkix_TargetCertCheckerState {
        PKIX_CertSelector *certSelector;
        PKIX_List *pathToNameList;
        PKIX_List *extKeyUsageList; /* List of PKIX_PL_OID */
        PKIX_List *subjAltNameList;
        PKIX_Boolean subjAltNameMatchAll;
        PKIX_UInt32 certsRemaining;
        PKIX_PL_OID *extKeyUsageOID;
        PKIX_PL_OID *subjAltNameOID;
};

PKIX_Error *
pkix_TargetCertChecker_Initialize(
        PKIX_CertSelector *certSelector,
        PKIX_UInt32 certsRemaining,
        PKIX_CertChainChecker **pChecker,
        void *plContext);

PKIX_Error *
pkix_TargetCertCheckerState_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_TARGETCERTCHECKER_H */
