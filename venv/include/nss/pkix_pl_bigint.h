/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_bigint.h
 *
 * Bigint Object Definitions
 *
 */

#ifndef _PKIX_PL_BIGINT_H
#define _PKIX_PL_BIGINT_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_BigIntStruct {
        char *dataRep;
        PKIX_UInt32 length;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_BigInt_CreateWithBytes(
        char *bytes,
        PKIX_UInt32 length,
        PKIX_PL_BigInt **pBigInt,
        void *plContext);

PKIX_Error *pkix_pl_BigInt_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_BIGINT_H */
