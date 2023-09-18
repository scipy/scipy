/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_pk11certstore.h
 *
 * PK11Certstore Object Type Definition
 *
 */

#ifndef _PKIX_PL_PK11CERTSTORE_H
#define _PKIX_PL_PK11CERTSTORE_H

#include "pkix_pl_common.h"
#include "certi.h"

#ifdef __cplusplus
extern "C" {
#endif

/* see source file for function documentation */
PKIX_Error *
PKIX_PL_Pk11CertStore_Create(
        PKIX_CertStore **pCertStore,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_PK11CERTSTORE_H */
