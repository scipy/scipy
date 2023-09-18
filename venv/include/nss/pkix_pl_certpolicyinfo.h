/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_certpolicyinfo.h
 *
 * PolicyInfo Type Definitions
 *
 */

#ifndef _PKIX_PL_POLICYINFO_H
#define _PKIX_PL_POLICYINFO_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This structure reflects the contents of the policy info extension as
 * described in Section 4.2.1.5 of RFC3280.
 *
 *  PolicyInformation ::= SEQUENCE {
 *      policyIdentifier        CertPolicyId,
 *      PolicyQualifiers        SEQUENCE SIZE (1..MAX) OF
 *                                  PolicyQualifierInfo OPTIONAL }
 *
 */
struct PKIX_PL_CertPolicyInfoStruct {
        PKIX_PL_OID *cpID;
        PKIX_List *policyQualifiers; /* LIST of PKIX_PL_CertPolicyQualifier */
};

PKIX_Error *
pkix_pl_CertPolicyInfo_Create(
        PKIX_PL_OID *oid,
        PKIX_List *qualifiers,
        PKIX_PL_CertPolicyInfo **pObject,
        void *plContext);

PKIX_Error *
pkix_pl_CertPolicyInfo_RegisterSelf(
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_POLICYINFO_H */
