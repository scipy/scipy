/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_comcrlselparams.h
 *
 * ComCrlSelParams Object Type Definition
 *
 */

#ifndef _PKIX_COMCRLSELPARAMS_H
#define _PKIX_COMCRLSELPARAMS_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_ComCRLSelParamsStruct {
        PKIX_List *issuerNames; /* list of PKIX_PL_X500Name */
        PKIX_PL_Cert *cert; /* certificate being checked */
        PKIX_List *crldpList;
        PKIX_PL_Date *date;
        PKIX_Boolean nistPolicyEnabled;
        PKIX_PL_BigInt *maxCRLNumber;
        PKIX_PL_BigInt *minCRLNumber;
};

/* see source file for function documentation */

PKIX_Error *pkix_ComCRLSelParams_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_COMCRLSELPARAMS_H */
