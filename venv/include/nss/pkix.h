/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines the public API for libpkix. These are the top-level
 * functions in the library. They perform the primary operations of this
 * library: building and validating chains of X.509 certificates.
 *
 */

#ifndef _PKIX_H
#define _PKIX_H

#include "pkixt.h"
#include "pkix_util.h"
#include "pkix_results.h"
#include "pkix_certstore.h"
#include "pkix_certsel.h"
#include "pkix_crlsel.h"
#include "pkix_checker.h"
#include "pkix_revchecker.h"
#include "pkix_pl_system.h"
#include "pkix_pl_pki.h"
#include "pkix_params.h"

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

/*
 * FUNCTION: PKIX_Initialize
 * DESCRIPTION:
 *
 * No PKIX_* types and functions should be used before this function is called
 * and returns successfully. This function should only be called once. If it
 * is called more than once, the behavior is undefined.
 *
 * NSS applications are expected to call NSS_Init, and need not know that
 * NSS will call this function (with "platformInitNeeded" set to PKIX_FALSE).
 * PKIX applications are expected instead to call this function with
 * "platformInitNeeded" set to PKIX_TRUE.
 *
 * This function initializes data structures critical to the operation of
 * libpkix. It also ensures that the API version (major.minor) desired by the
 * caller (the "desiredMajorVersion", "minDesiredMinorVersion", and
 * "maxDesiredMinorVersion") is compatible with the API version supported by
 * the library. As such, the library must support the "desiredMajorVersion"
 * of the API and must support a minor version that falls between
 * "minDesiredMinorVersion" and "maxDesiredMinorVersion", inclusive. If
 * compatibility exists, the function returns NULL and stores the library's
 * actual minor version at "pActualMinorVersion" (which may be greater than
 * "desiredMinorVersion"). If no compatibility exists, the function returns a
 * PKIX_Error pointer. If the caller wishes to specify that the largest
 * minor version available should be used, then maxDesiredMinorVersion should
 * be set to the macro PKIX_MAX_MINOR_VERSION (defined in pkixt.h).
 *
 * PARAMETERS:
 *  "platformInitNeeded"
 *      Boolean indicating whether the platform layer initialization code
 *      has previously been run, or should be called from this function.
 *  "desiredMajorVersion"
 *      The major version of the libpkix API the application wishes to use.
 *  "minDesiredMinorVersion"
 *      The minimum minor version of the libpkix API the application wishes
 *      to use.
 *  "maxDesiredMinorVersion"
 *      The maximum minor version of the libpkix API the application wishes
 *      to use.
 *  "pActualMinorVersion"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "pPlContext"
 *      Address at which platform-specific context pointer is stored. Must
 *      be non-NULL.
 * THREAD SAFETY:
 *  Not Thread Safe
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Initialize Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Initialize(
        PKIX_Boolean platformInitNeeded,
        PKIX_UInt32 desiredMajorVersion,
        PKIX_UInt32 minDesiredMinorVersion,
        PKIX_UInt32 maxDesiredMinorVersion,
        PKIX_UInt32 *pActualMinorVersion,
        void **pPlContext);

/*
 * FUNCTION: PKIX_Shutdown
 * DESCRIPTION:
 *
 *  This function deallocates any memory used by libpkix and shuts down any
 *  ongoing operations. This function should only be called once. If it is
 *  called more than once, the behavior is undefined.
 *
 *  No PKIX_* types and functions should be used after this function is called
 *  and returns successfully.
 * PARAMETERS:
 *  "plContext" - Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Shutdown(void *plContext);

/*
 * FUNCTION: PKIX_ValidateChain
 * DESCRIPTION:
 *
 *  This function attempts to validate the CertChain that has been set in the
 *  ValidateParams pointed to by "params" using an RFC 3280-compliant
 *  algorithm. If successful, this function returns NULL and stores the
 *  ValidateResult at "pResult", which holds additional information, such as
 *  the policy tree and the target's public key. If unsuccessful, an Error is
 *  returned. Note: This function does not currently support non-blocking I/O.
 *
 *  If "pVerifyTree" is non-NULL, a chain of VerifyNodes is created which
 *  tracks the results of the validation. That is, either each node in the
 *  chain has a NULL Error component, or the last node contains an Error
 *  which indicates why the validation failed.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ValidateParams used to validate CertChain. Must be non-NULL.
 *  "pResult"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "pVerifyTree"
 *      Address where a VerifyTree is stored, if non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (See Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Validate Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateChain(
        PKIX_ValidateParams *params,
        PKIX_ValidateResult **pResult,
	PKIX_VerifyNode **pVerifyTree,
        void *plContext);

/*
 * FUNCTION: PKIX_ValidateChain_NB
 * DESCRIPTION:
 *
 *  This function is the equivalent of PKIX_ValidateChain, except that it
 *  supports non-blocking I/O. When called with "pNBIOContext" pointing to NULL
 *  it initiates a new chain validation as in PKIX_ValidateChain, ignoring the
 *  value in all input variables except "params". If forced to suspend
 *  processing by a WOULDBLOCK return from some operation, such as a CertStore
 *  request, it stores the platform-dependent I/O context at "pNBIOContext" and
 *  stores other intermediate variables at "pCertIndex", "pAnchorIndex", 
 *  "pCheckerIndex", "pRevChecking", and "pCheckers".
 *
 *  When called subsequently with that non-NULL value at "pNBIOContext", it
 *  relies on those intermediate values to be untouched, and it resumes chain
 *  validation where it left off. Its behavior is undefined if any of the
 *  intermediate values was not preserved.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ValidateParams used to validate CertChain. Must be non-NULL.
 *  "pCertIndex"
 *      The UInt32 value of the index to the Cert chain, indicating which Cert
 *      is currently being processed.
 *  "pAnchorIndex"
 *      The UInt32 value of the index to the Anchor chain, indicating which
 *      Trust Anchor is currently being processed.
 *  "pCheckerIndex"
 *      The UInt32 value of the index to the List of CertChainCheckers,
 *      indicating which Checker is currently processing.
 *  "pRevChecking"
 *      The Boolean flag indicating whether normal checking or revocation
 *      checking is occurring for the Cert indicated by "pCertIndex".
 *  "pCheckers"
 *      The address of the List of CertChainCheckers. Must be non-NULL.
 *  "pNBIOContext"
 *      The address of the platform-dependend I/O context. Must be a non-NULL
 *      pointer to a NULL value for the call to initiate chain validation.
 *  "pResult"
 *      Address where ValidateResult object pointer will be stored. Must be
 *      non-NULL.
 *  "pVerifyTree"
 *      Address where a VerifyTree is stored, if non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a VALIDATE Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */PKIX_Error *
PKIX_ValidateChain_NB(
	PKIX_ValidateParams *params,
	PKIX_UInt32 *pCertIndex,
	PKIX_UInt32 *pAnchorIndex,
	PKIX_UInt32 *pCheckerIndex,
	PKIX_Boolean *pRevChecking,
	PKIX_List **pCheckers,
	void **pNBIOContext,
	PKIX_ValidateResult **pResult,
	PKIX_VerifyNode **pVerifyTree,
	void *plContext);

/*
 * FUNCTION: PKIX_BuildChain
 * DESCRIPTION:
 *
 *  If called with a NULL "state", this function attempts to build and validate
 *  a CertChain according to the ProcessingParams pointed to by "params", using
 *  an RFC 3280-compliant validation algorithm. If successful, this function
 *  returns NULL and stores the BuildResult at "pResult", which holds the built
 *  CertChain, as well as additional information, such as the policy tree and
 *  the target's public key. If unsuccessful, an Error is returned.
 *
 *  If the chain building is blocked by a CertStore using non-blocking I/O, this
 *  function stores platform-dependent non-blocking I/O context at
 *  "pNBIOContext", its state at "pState", and NULL at "pResult". The caller
 *  may be able to determine, in a platform-dependent way, when the I/O has
 *  completed. In any case, calling the function again with "pState" containing
 *  the returned value will allow the chain building to resume.
 *
 *  If chain building is completed, either successfully or unsuccessfully, NULL
 *  is stored at "pNBIOContext".
 *
 *  If "pVerifyTree" is non-NULL, a tree of VerifyNodes is created which
 *  tracks the results of the building. That is, each node of the tree either
 *  has a NULL Error component, or it is a leaf node and it contains an Error
 *  which indicates why the chain building could not proceed on this branch.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams used to build and validate CertChain.
 *      Must be non-NULL.
 *  "pNBIOContext"
 *      Address where platform-dependent information is store if the build
 *      is suspended waiting for non-blocking I/O. Must be non-NULL.
 *  "pState"
 *      Address of BuildChain state. Must be NULL on initial call, and the
 *      value previously returned on subsequent calls.
 *  "pResult"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "pVerifyTree"
 *      Address where a VerifyTree is stored, if non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (See Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Build Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_BuildChain(
        PKIX_ProcessingParams *params,
        void **pNBIOContext,
        void **pState,
        PKIX_BuildResult **pResult,
	PKIX_VerifyNode **pVerifyNode,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_H */
