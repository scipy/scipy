/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_nameconstraintschecker.h
 *
 * Header file for validate Name Constraints Checker function
 *
 */

#ifndef _PKIX_NAMECONSTRAINTSCHECKER_H
#define _PKIX_NAMECONSTRAINTSCHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pkix_NameConstraintsCheckerState \
        pkix_NameConstraintsCheckerState;

struct pkix_NameConstraintsCheckerState {
        PKIX_PL_CertNameConstraints *nameConstraints;
        PKIX_PL_OID *nameConstraintsOID;
        PKIX_UInt32 certsRemaining;
};

PKIX_Error *
pkix_NameConstraintsChecker_Initialize(
        PKIX_PL_CertNameConstraints *trustedNC,
        PKIX_UInt32 numCerts,
        PKIX_CertChainChecker **pChecker,
        void *plContext);

PKIX_Error *
pkix_NameConstraintsCheckerState_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_NAMECONSTRAINTSCHECKER_H */
