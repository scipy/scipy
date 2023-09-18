/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_nsscontext.h
 *
 * NSSContext Object Type Definition
 *
 */


#ifndef _PKIX_PL_NSSCONTEXT_H
#define _PKIX_PL_NSSCONTEXT_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_NssContextStruct {
        SECCertificateUsage certificateUsage;
        PLArenaPool *arena;
        void *wincx;
        PKIX_UInt32 timeoutSeconds;
        PKIX_UInt32 maxResponseLength;
        PRTime crlReloadDelay;
        PRTime badDerCrlReloadDelay;
        CERTChainVerifyCallback chainVerifyCallback;
        PKIX_Boolean certSignatureCheck;
};

PKIX_Error *
pkix_pl_NssContext_GetCertUsage
        (PKIX_PL_NssContext *nssContext, SECCertificateUsage *pCertUsage);

/* XXX move the setter into the public header. */
PKIX_Error *
pkix_pl_NssContext_SetCertUsage
        (SECCertificateUsage certUsage, PKIX_PL_NssContext *nssContext);

PKIX_Error *
pkix_pl_NssContext_GetCertSignatureCheck
        (PKIX_PL_NssContext *nssContext, PKIX_Boolean *pCheckSig);

PKIX_Error *
pkix_pl_NssContext_SetCertSignatureCheck
        (PKIX_Boolean checkSig, PKIX_PL_NssContext *nssContext);

PKIX_Error *
pkix_pl_NssContext_GetWincx(PKIX_PL_NssContext *nssContext, void **pWincx);

/* XXX move the setter into the public header. */
PKIX_Error *
pkix_pl_NssContext_SetWincx(void *wincx, PKIX_PL_NssContext *nssContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_NSSCONTEXT_H */
