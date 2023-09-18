/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_revocationmethod.h
 *
 * RevocationMethod Object
 *
 */

#ifndef _PKIX_REVOCATIONMETHOD_H
#define _PKIX_REVOCATIONMETHOD_H

#include "pkixt.h"
#include "pkix_revocationchecker.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pkix_RevocationMethodStruct pkix_RevocationMethod;

/* Local revocation check function prototype definition.
 * Revocation methods capable of checking revocation though local
 * means(cache) should implement this prototype. */
typedef PKIX_Error *
pkix_LocalRevocationCheckFn(PKIX_PL_Cert *cert, PKIX_PL_Cert *issuer,
                            PKIX_PL_Date *date, 
                            pkix_RevocationMethod *checkerObject,
                            PKIX_ProcessingParams *procParams,
                            PKIX_UInt32 methodFlags,
                            PKIX_Boolean chainVerificationState,
                            PKIX_RevocationStatus *pRevStatus,
                            CERTCRLEntryReasonCode *reasonCode,
                            void *plContext);

/* External revocation check function prototype definition.
 * Revocation methods that required external communications(crldp
 * ocsp) shoult implement this prototype. */
typedef PKIX_Error *
pkix_ExternalRevocationCheckFn(PKIX_PL_Cert *cert, PKIX_PL_Cert *issuer,
                               PKIX_PL_Date *date,
                               pkix_RevocationMethod *checkerObject,
                               PKIX_ProcessingParams *procParams,
                               PKIX_UInt32 methodFlags,
                               PKIX_RevocationStatus *pRevStatus,
                               CERTCRLEntryReasonCode *reasonCode,
                               void **pNBIOContext, void *plContext);

/* Revocation method structure assosiates revocation types with
 * a set of flags on the method, a priority of the method (0
 * corresponds to the highest priority), and method local/external
 * checker functions. */
struct pkix_RevocationMethodStruct {
    PKIX_RevocationMethodType methodType;
    PKIX_UInt32 flags;
    PKIX_UInt32 priority;
    pkix_LocalRevocationCheckFn (*localRevChecker);
    pkix_ExternalRevocationCheckFn (*externalRevChecker);
};

PKIX_Error *
pkix_RevocationMethod_Duplicate(PKIX_PL_Object *object,
                                PKIX_PL_Object *newObject,
                                void *plContext);

PKIX_Error *
pkix_RevocationMethod_Init(pkix_RevocationMethod *method,
                           PKIX_RevocationMethodType methodType,
                           PKIX_UInt32 flags,
                           PKIX_UInt32 priority,
                           pkix_LocalRevocationCheckFn localRevChecker,
                           pkix_ExternalRevocationCheckFn externalRevChecker,
                           void *plContext);


#ifdef __cplusplus
}
#endif

#endif /* _PKIX_REVOCATIONMETHOD_H */
