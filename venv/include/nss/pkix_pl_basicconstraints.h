/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_basicconstraints.h
 *
 * BasicConstraints Object Definitions
 *
 */

#ifndef _PKIX_PL_BASICCONSTRAINTS_H
#define _PKIX_PL_BASICCONSTRAINTS_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This structure reflects the contents of the basic constraints
 * extension as described in Section 4.2.1.10 of RFC 3280.
 * The cA flag indicates whether the public key in this certificate
 * belongs to a certification authority. The pathLen constraint
 * gives the maximum number of non-self-issued intermediate certificates
 * that may follow this certificate in a valid certification path.
 */
struct PKIX_PL_CertBasicConstraintsStruct {
        PKIX_Boolean isCA;
        PKIX_Int32 pathLen;
};

PKIX_Error *
pkix_pl_CertBasicConstraints_Create(
        PKIX_Boolean isCA,
        PKIX_Int32 pathLen,
        PKIX_PL_CertBasicConstraints **object,
        void *plContext);

PKIX_Error *
pkix_pl_CertBasicConstraints_RegisterSelf(
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_BASICCONSTRAINTS_H */
