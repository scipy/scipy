/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_namechainingchecker.h
 *
 * Header file for name chaining checker.
 *
 */

#ifndef _PKIX_NAMECHAININGCHECKER_H
#define _PKIX_NAMECHAININGCHECKER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

PKIX_Error *
pkix_NameChainingChecker_Initialize(
        PKIX_PL_X500Name *trustedCAName,
        PKIX_CertChainChecker **pChecker,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_NAMECHAININGCHECKER_H */
