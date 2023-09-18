/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_trustanchor.h
 *
 * TrustAnchor Object Type Definition
 *
 */

#ifndef _PKIX_TRUSTANCHOR_H
#define _PKIX_TRUSTANCHOR_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_TrustAnchorStruct {
        PKIX_PL_Cert *trustedCert;
        PKIX_PL_X500Name *caName;
        PKIX_PL_PublicKey *caPubKey;
        PKIX_PL_CertNameConstraints *nameConstraints;
};

/* see source file for function documentation */

PKIX_Error *pkix_TrustAnchor_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_TRUSTANCHOR_H */
