/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_procparams.h
 *
 * ProcessingParams Object Type Definition
 *
 */

#ifndef _PKIX_PROCESSINGPARAMS_H
#define _PKIX_PROCESSINGPARAMS_H

#include "pkix_tools.h"


#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_ProcessingParamsStruct {
        PKIX_List *trustAnchors;        /* Never NULL */
        PKIX_List *hintCerts;	/* user-supplied partial chain, may be NULL */
        PKIX_CertSelector *constraints;
        PKIX_PL_Date *date;
        PKIX_List *initialPolicies;     /* list of PKIX_PL_OID */
        PKIX_Boolean initialPolicyMappingInhibit;
        PKIX_Boolean initialAnyPolicyInhibit;
        PKIX_Boolean initialExplicitPolicy;
        PKIX_Boolean qualifiersRejected;
        PKIX_List *certChainCheckers;
        PKIX_List *certStores;
        PKIX_Boolean isCrlRevocationCheckingEnabled;
        PKIX_Boolean isCrlRevocationCheckingEnabledWithNISTPolicy;
        PKIX_RevocationChecker *revChecker;
        PKIX_ResourceLimits *resourceLimits;
        PKIX_Boolean useAIAForCertFetching;
        PKIX_Boolean qualifyTargetCert;
        PKIX_Boolean useOnlyTrustAnchors;
};

/* see source file for function documentation */

PKIX_Error *pkix_ProcessingParams_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PROCESSINGPARAMS_H */
