/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_ekuchecker.h
 *
 * User Defined Object Type Extended Key Usage Definition
 *
 */

#ifndef _PKIX_EKUCHECKER_H
#define _PKIX_EKUCHECKER_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * FUNCTION: PKIX_PL_EkuChecker_Create
 *
 * DESCRIPTION:
 *  Create a CertChainChecker with EkuCheckerState and add it into
 *  PKIX_ProcessingParams object.
 *
 * PARAMETERS
 *  "params"
 *      a PKIX_ProcessingParams links to PKIX_ComCertSelParams where a list of
 *      Extended Key Usage OIDs specified by application can be retrieved for
 *      verification.
 *  "ekuChecker" 
 *      Address of created ekuchecker.
 *  "plContext"
 *      Platform-specific context pointer.
 *
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 *
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a UserDefinedModules Error if the function fails in a non-fatal
 *  way.
 *  Returns a Fatal Error
 */
PKIX_Error *
PKIX_EkuChecker_Create(
        PKIX_ProcessingParams *params,
        PKIX_CertChainChecker **ekuChecker,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_EkuChecker_GetRequiredEku
 *
 * DESCRIPTION:
 *  This function retrieves application specified ExtenedKeyUsage(s) from
 *  ComCertSetparams and converts its OID representations to SECCertUsageEnum.
 *  The result is stored and returned in bit mask at "pRequiredExtKeyUsage".
 *
 * PARAMETERS
 *  "certSelector"
 *      a PKIX_CertSelector links to PKIX_ComCertSelParams where a list of
 *      Extended Key Usage OIDs specified by application can be retrieved for
 *      verification. Must be non-NULL.
 *  "pRequiredExtKeyUsage"
 *      Address where the result is returned. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 *
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 *
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a UserDefinedModules Error if the function fails in a non-fatal
 *  way.
 *  Returns a Fatal Error
 */
PKIX_Error *
pkix_EkuChecker_GetRequiredEku(
        PKIX_CertSelector *certSelector,
        PKIX_UInt32 *pRequiredExtKeyUsage,
        void *plContext);

/* see source file for function documentation */
PKIX_Error *pkix_EkuChecker_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_EKUCHECKER_H */
