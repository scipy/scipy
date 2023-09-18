/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_policychecker.h
 *
 * Header file for policy checker.
 *
 */

#ifndef _PKIX_POLICYCHECKER_H
#define _PKIX_POLICYCHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PKIX_PolicyCheckerStateStruct PKIX_PolicyCheckerState;

struct PKIX_PolicyCheckerStateStruct{
        PKIX_PL_OID *certPoliciesExtension;             /* const */
        PKIX_PL_OID *policyMappingsExtension;           /* const */
        PKIX_PL_OID *policyConstraintsExtension;        /* const */
        PKIX_PL_OID *inhibitAnyPolicyExtension;         /* const */
        PKIX_PL_OID *anyPolicyOID;                      /* const */
        PKIX_Boolean initialIsAnyPolicy;                /* const */
        PKIX_PolicyNode *validPolicyTree;
        PKIX_List *userInitialPolicySet;                /* immutable */
        PKIX_List *mappedUserInitialPolicySet;
        PKIX_Boolean policyQualifiersRejected;
        PKIX_Boolean initialPolicyMappingInhibit;
        PKIX_Boolean initialExplicitPolicy;
        PKIX_Boolean initialAnyPolicyInhibit;
        PKIX_UInt32 explicitPolicy;
        PKIX_UInt32 inhibitAnyPolicy;
        PKIX_UInt32 policyMapping;
        PKIX_UInt32 numCerts;
        PKIX_UInt32 certsProcessed;
        PKIX_PolicyNode *anyPolicyNodeAtBottom;
        PKIX_PolicyNode *newAnyPolicyNode;
        /*
         * The following variables do not survive from one
         * certificate to the next. They are needed at each
         * level of recursive routines, any by placing them
         * in the state object we can pass fewer arguments.
         */
        PKIX_Boolean certPoliciesCritical;
        PKIX_List *mappedPolicyOIDs;
};

PKIX_Error *
pkix_PolicyChecker_Initialize(
        PKIX_List *initialPolicies,
        PKIX_Boolean policyQualifiersRejected,
        PKIX_Boolean initialPolicyMappingInhibit,
        PKIX_Boolean initialExplicitPolicy,
        PKIX_Boolean initialAnyPolicyInhibit,
        PKIX_UInt32 numCerts,
        PKIX_CertChainChecker **pChecker,
        void *plContext);

/* --Private-Functions-------------------------------------------- */

PKIX_Error *
pkix_PolicyCheckerState_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_POLICYCHECKER_H */
