/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the PKIX_RevocationChecker
 * type.
 *
 */

#ifndef _PKIX_REVCHECKER_H
#define _PKIX_REVCHECKER_H

#include "pkixt.h"
#include "pkix_pl_pki.h"

#ifdef __cplusplus
extern "C" {
#endif

/* General
 *
 * Please refer to the libpkix Programmer's Guide for detailed information
 * about how to use the libpkix library. Certain key warnings and notices from
 * that document are repeated here for emphasis.
 *
 * All identifiers in this file (and all public identifiers defined in
 * libpkix) begin with "PKIX_". Private identifiers only intended for use
 * within the library begin with "pkix_".
 *
 * A function returns NULL upon success, and a PKIX_Error pointer upon failure.
 *
 * Unless otherwise noted, for all accessor (gettor) functions that return a
 * PKIX_PL_Object pointer, callers should assume that this pointer refers to a
 * shared object. Therefore, the caller should treat this shared object as
 * read-only and should not modify this shared object. When done using the
 * shared object, the caller should release the reference to the object by
 * using the PKIX_PL_Object_DecRef function.
 *
 * While a function is executing, if its arguments (or anything referred to by
 * its arguments) are modified, free'd, or destroyed, the function's behavior
 * is undefined.
 *
 */

/* PKIX_RevocationChecker
 *
 * PKIX_RevocationChecker provides a standard way of revocation checking.
 * Caller should configure two set of tests(represented at lists of
 * RevocationMethod objects) to be performed on the leaf and on the rest of
 * the chain certificates.
 *
 * PKIX_RevocationMethods provide a standard way for the caller to insert
 * their own custom revocation checks to verify the revocation status of
 * certificates. This may be useful in many scenarios, including when the
 * caller wishes to use their own revocation checking mechanism instead of (or
 * in addition to) the default revocation checking mechanism provided by
 * libpkix, which uses CRLs and OCSP. 
 *
 * Once the caller has created the RevocationMethod object(s), the caller
 * then specifies the RevocationMethod object(s) in a RevocationCheck object
 * and sets it into a ProcessingParams.
 */

/*
 * FUNCTION: PKIX_RevocationChecker_Create
 * DESCRIPTION:
 *
 * Creates a revocation checker object with the given flags. Revocation will
 * be checked at the current date.
 *
 * PARAMETERS:
 *  "leafMethodListFlags"
 *      Defines a set of method independent flags that will be used to check
 *      revocation of the leaf cert in the chain.
 *  "chainMethodListFlags"
 *      Defines a set of method independent flags that will be used to check
 *      revocation of the remaining certs in the chain.
 *  "pChecker"
 *      The return address of created checker.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a RevocationChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_RevocationChecker_Create(
    PKIX_UInt32 leafMethodListFlags,
    PKIX_UInt32 chainMethodListFlags,
    PKIX_RevocationChecker **pChecker,
    void *plContext);

/*
 * FUNCTION: PKIX_RevocationChecker_CreateAndAddMethod
 * DESCRIPTION:
 *
 * Creates revocation method object with given parameters and adds it
 * to revocation checker method list.
 *
 * PARAMETERS:
 *  "revChecker"
 *      Address of revocation checker structure.
 *  "procParams"
 *      Address of ProcessingParams used to initialize the checker.
 *      Must be non-NULL.
 *  "methodType"
 *      Type of the method. Currently only two types are
 *      supported: crl and ocsp. (See PKIX_RevocationMethodType enum).
 *  "methodFlags"
 *      Set of flags for the method.
 *  "methodPriority"
 *      Method priority. (0 corresponds to the highest priority)
 *  "verificationFn"
 *      User call back function that will perform validation of fetched
 *      revocation information(new crl or ocsp response)
 *  "isLeafMethod"
 *      Boolean flag that if set to true indicates that the method should
 *      should be used for leaf cert revocation test(false for chain set
 *      methods).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a RevocationChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_RevocationChecker_CreateAndAddMethod(
    PKIX_RevocationChecker *revChecker,
    PKIX_ProcessingParams *params,
    PKIX_RevocationMethodType methodType,
    PKIX_UInt32 methodFlags,
    PKIX_UInt32 methodPriority,
    PKIX_PL_VerifyCallback verificationFn,
    PKIX_Boolean isLeafMethod,
    void *plContext);

/*
 * FUNCTION: PKIX_RevocationChecker_Check
 * DESCRIPTION:
 *
 * Verifies revocation status of the certificate. Issuer cert is given to
 * be used in verification of revocation information. Performed verification
 * check depends on configured revocation methods(ocsp, crl. See
 * PKIX_RevocationChecker_CreateAndAddMethod function) and a point of chain
 * building process at which PKIX_RevocationChecker_Check was invoked.
 * For security reasons, the cert status is checked only against cached
 * revocation information during chain building stage(no trust anchor yes has
 * been found). The fresh revocation information fetching is done only at chain
 * verification stage after trust anchor was identified.
 * 
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose revocation status is to be determined.
 *      Must be non-NULL.
 *  "issuer"
 *      Issuer cert that potentially holds public key that will be used
 *      to verify revocation info.
 *  "revChecker"
 *      Address of revocation checker structure.
 *  "procParams"
 *      Address of ProcessingParams used to initialize the checker.
 *      Must be non-NULL.
 *  "chainVerificationState"
 *     Need to be set to true, if the check was called during chain verification
 *     as an opposite to chain building.
 *  "testingLeafCert"
 *     Set to true if verifying revocation status of a leaf cert.
 *  "revStatus"
 *     Address of the returned revocation status of the cert.
 *  "pResultCode"
 *      Address where revocation status will be stored. Must be non-NULL.
 *  "pNBIOContext"
 *      Address at which platform-dependent non-blocking I/O context is stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a RevocationChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_RevocationChecker_Check(PKIX_PL_Cert *cert,
                             PKIX_PL_Cert *issuer,
                             PKIX_RevocationChecker *revChecker,
                             PKIX_ProcessingParams *procParams,
                             PKIX_Boolean chainVerificationState,
                             PKIX_Boolean testingLeafCert,
                             PKIX_RevocationStatus *revStatus,
                             PKIX_UInt32 *pReasonCode,
                             void **pNbioContext,
                             void *plContext);
    
#ifdef __cplusplus
}
#endif

#endif /* _PKIX_REVCHECKER_H */
