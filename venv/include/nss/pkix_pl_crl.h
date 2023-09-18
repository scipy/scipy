/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_crl.h
 *
 * CRL Object Type Definitions
 *
 */

#ifndef _PKIX_PL_CRL_H
#define _PKIX_PL_CRL_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_CRLStruct {
        CERTSignedCrl *nssSignedCrl;
        PKIX_PL_X500Name *issuer;
        PKIX_PL_OID *signatureAlgId;
        PKIX_PL_BigInt *crlNumber;
        PKIX_Boolean crlNumberAbsent;
        PKIX_List *crlEntryList; /* list of PKIX_PL_CRLEntry */
        PKIX_List *critExtOids;
        SECItem *adoptedDerCrl;
        SECItem *derGenName; /* der of general name which was used
                              * to download the crl. */
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_CRL_RegisterSelf(void *plContext);

PKIX_Error *
pkix_pl_CRL_CreateWithSignedCRL(CERTSignedCrl *nssSignedCrl,
                                SECItem *derCrl,
                                SECItem *derGenName,
                                PKIX_PL_CRL **pCrl,
                                void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_CRL_H */
