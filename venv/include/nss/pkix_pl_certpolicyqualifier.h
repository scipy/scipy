/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_certpolicyqualifier.h
 *
 * PolicyQualifier Type Definitions
 *
 */

#ifndef _PKIX_PL_POLICYQUALIFIER_H
#define _PKIX_PL_POLICYQUALIFIER_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This structure reflects the contents of the policy qualifier extension as
 * described in Section 4.2.1.5 of RFC3280.
 *
 *  PolicyQualifierInfo ::= SEQUENCE {
 *      policyQualifierId       PolicyQualifierId,
 *      qualifier               ANY DEFINED BY policyQualifierId }
 *
 *  PolicyQualifierId ::=
 *      OBJECT IDENTIFIER (id-qt-cps | id-qt-unotice)
 *
 */
struct PKIX_PL_CertPolicyQualifierStruct {
        PKIX_PL_OID *policyQualifierId;
        PKIX_PL_ByteArray *qualifier;
};

PKIX_Error *
pkix_pl_CertPolicyQualifier_Create(
        PKIX_PL_OID *oid,
        PKIX_PL_ByteArray *qualifierArray,
        PKIX_PL_CertPolicyQualifier **pObject,
        void *plContext);

PKIX_Error *
pkix_pl_CertPolicyQualifier_RegisterSelf(
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_POLICYQUALIFIER_H */
