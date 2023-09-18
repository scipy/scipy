/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_expirationchecker.h
 *
 * Header file for validate expiration function
 *
 */

#ifndef _PKIX_EXPIRATIONCHECKER_H
#define _PKIX_EXPIRATIONCHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

PKIX_Error *
pkix_ExpirationChecker_Initialize(
        PKIX_PL_Date *testDate,
        PKIX_CertChainChecker **pChecker,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_EXPIRATIONCHECKER_H */
