/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_store.h
 *
 * CertStore Object Type Definition
 *
 */

#ifndef _PKIX_STORE_H
#define _PKIX_STORE_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_CertStoreStruct {
        PKIX_CertStore_CertCallback certCallback;
        PKIX_CertStore_CRLCallback crlCallback;
        PKIX_CertStore_CertContinueFunction certContinue;
        PKIX_CertStore_CrlContinueFunction crlContinue;
        PKIX_CertStore_CheckTrustCallback trustCallback;
        PKIX_CertStore_ImportCrlCallback importCrlCallback;
        PKIX_CertStore_CheckRevokationByCrlCallback checkRevByCrlCallback;
        PKIX_PL_Object *certStoreContext;
        PKIX_Boolean cacheFlag;
        PKIX_Boolean localFlag; /* TRUE if CertStore is local */
};

/* see source file for function documentation */

PKIX_Error *pkix_CertStore_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_STORE_H */
