/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_crldp.h
 *
 * Crp DP Object Definitions
 *
 */
#include "pkix_pl_common.h"

#ifndef _PKIX_PL_CRLDP_H
#define _PKIX_PL_CRLDP_H

#ifdef __cplusplus
extern "C" {
#endif

/* CRLDP object can not be used without holding a reference
 * to the pkix certificate they belong to. The memory for dp der
 * object is allocated on nssCert certificate - a member of
 * PKIX_PL_Cert struct. */
typedef struct pkix_pl_CrlDpStruct {
    /* reference to decoded crldp that allocated on nssCert arena. */
    const CRLDistributionPoint *nssdp;
    DistributionPointTypes distPointType;
    union {
	CERTGeneralName *fullName;
        /* if dp is a relative name, the issuerName is a merged value
         * of crlIssuer and a relative name. Must be destroyed by CrlDp
         * destructor. */
        CERTName *issuerName;
    } name;
    PKIX_Boolean isPartitionedByReasonCode;
} pkix_pl_CrlDp;


PKIX_Error *
pkix_pl_CrlDp_RegisterSelf(void *plContext);

/* Parses CRLDistributionPoint structure and creaetes
 * pkix_pl_CrlDp object. */
PKIX_Error *
pkix_pl_CrlDp_Create(const CRLDistributionPoint *dp,
                     const CERTName *certIssuerName,
                     pkix_pl_CrlDp **pPkixDP,
                     void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_CRLDP_H */
