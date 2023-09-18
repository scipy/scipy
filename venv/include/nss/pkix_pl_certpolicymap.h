/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_certpolicymap.h
 *
 * CertPolicyMap Object Definitions
 *
 */

#ifndef _PKIX_PL_CERTPOLICYMAP_H
#define _PKIX_PL_CERTPOLICYMAP_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This structure reflects the contents of the policy mapping extension as
 * described in Section 4.2.1.6 of RFC3280.
 *
 *  PolicyMappings ::= SEQUENCE SIZE (1..MAX) OF SEQUENCE {
 *      issuerDomainPolicy      CertPolicyId,
 *      subjectDomainPolicy     CertPolicyId }
 *
 */
struct PKIX_PL_CertPolicyMapStruct {
        PKIX_PL_OID *issuerDomainPolicy;
        PKIX_PL_OID *subjectDomainPolicy;
};

PKIX_Error *
pkix_pl_CertPolicyMap_Create(
        PKIX_PL_OID *issuerDomainPolicy,
        PKIX_PL_OID *subjectDomainPolicy,
        PKIX_PL_CertPolicyMap **pObject,
        void *plContext);

PKIX_Error *
pkix_pl_CertPolicyMap_RegisterSelf(
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_CERTPOLICYMAP_H */
