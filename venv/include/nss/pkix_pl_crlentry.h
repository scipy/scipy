/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_crlentry.h
 *
 * CRL Entry Type Object Definitions
 *
 */

#ifndef _PKIX_PL_CRLENTRY_H
#define _PKIX_PL_CRLENTRY_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PKIX_PL_CRL_REASONCODE_NOTSET (-1)

struct PKIX_PL_CRLEntryStruct {
        CERTCrlEntry *nssCrlEntry;
        PKIX_PL_BigInt *serialNumber;
        PKIX_List *critExtOids;
        PKIX_Int32 userReasonCode;
        PKIX_Boolean userReasonCodeAbsent;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_CRLEntry_RegisterSelf(void *plContext);

/* following functions are called by CRL only hence not public */

PKIX_Error *
pkix_pl_CRLEntry_Create(
        CERTCrlEntry **nssCrlEntry, /* head of entry list */
        PKIX_List **pCrlEntryList,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_CRLENTRY_H */
