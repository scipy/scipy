/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the PKIX_CertChainChecker type.
 *
 */

#ifndef _PKIX_CHECKER_H
#define _PKIX_CHECKER_H

#include "pkixt.h"

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

/* PKIX_CertChainChecker
 *
 * PKIX_CertChainCheckers provide a standard way for the caller to insert their
 * own custom checks to validate certificates. This may be useful in many
 * scenarios, including when the caller wishes to validate private certificate
 * extensions. The CheckCallback allows custom certificate processing to take
 * place. Additionally, a CertChainChecker can optionally maintain state
 * between successive calls to the CheckCallback. This certChainCheckerState
 * must be an Object (although any object type), allowing it to be
 * reference-counted and allowing it to provide the standard Object functions
 * (Equals, Hashcode, ToString, Compare, Duplicate). If the caller wishes
 * their CertChainChecker to be used during chain building, their
 * certChainCheckerState object must implement an appropriate Duplicate
 * function. The builder uses this Duplicate function when backtracking.
 *
 * Once the caller has created a CertChainChecker object, the caller then
 * specifies a CertChainChecker object in a ProcessingParams object
 * and passes the ProcessingParams object to PKIX_ValidateChain or
 * PKIX_BuildChain, which uses the objects to call the user's callback
 * functions as needed during the validation or building process.
 *
 * A CertChainChecker may be presented certificates in the "reverse" direction
 * (from trust anchor to target) or in the "forward" direction (from target to
 * trust anchor). All CertChainCheckers must support "reverse checking", while
 * support for "forward checking" is optional, but recommended. If "forward
 * checking" is not supported, building chains may be much less efficient. The
 * PKIX_CertChainChecker_IsForwardCheckingSupported function is used to
 * determine whether forward checking is supported, and the
 * PKIX_CertChainChecker_IsForwardDirectionExpected function is used to
 * determine whether the CertChainChecker has been initialized to expect the
 * certificates to be presented in the "forward" direction.
 */

/*
 * FUNCTION: PKIX_CertChainChecker_CheckCallback
 * DESCRIPTION:
 *
 *  This callback function checks whether the specified Cert pointed to by
 *  "cert" is valid using "checker's" internal certChainCheckerState (if any)
 *  and removes the critical extensions that it processes (if any) from the
 *  List of OIDs (possibly empty) pointed to by "unresolvedCriticalExtensions".
 *  If the checker finds that the certificate is not valid, an Error pointer is
 *  returned.
 *
 *  If the checker uses non-blocking I/O, the address of a platform-dependent
 *  non-blocking I/O context ("nbioContext") will be stored at "pNBIOContext",
 *  which the caller may use, in a platform-dependent way, to wait, poll, or
 *  otherwise determine when to try again. If the checker does not use
 *  non-blocking I/O, NULL will always be stored at "pNBIOContext". If a non-NULL
 *  value was stored, on a subsequent call the checker will attempt to complete
 *  the pending I/O and, if successful, NULL will be stored at "pNBIOContext".
 *
 * PARAMETERS:
 *  "checker"
 *      Address of CertChainChecker whose certChainCheckerState and
 *      CheckCallback logic is to be used. Must be non-NULL.
 *  "cert"
 *      Address of Cert that is to be validated using "checker".
 *      Must be non-NULL.
 *  "unresolvedCriticalExtensions"
 *      Address of List of OIDs that represents the critical certificate
 *      extensions that have yet to be resolved. This parameter may be
 *      modified during the function call. Must be non-NULL.
 *  "pNBIOContext"
 *      Address at which is stored a platform-dependent structure indicating
 *      whether checking was suspended for non-blocking I/O. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertChainChecker_CheckCallback)(
        PKIX_CertChainChecker *checker,
        PKIX_PL_Cert *cert,
        PKIX_List *unresolvedCriticalExtensions,  /* list of PKIX_PL_OID */
        void **pNBIOContext,
        void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_Create
 * DESCRIPTION:
 *
 *  Creates a new CertChainChecker and stores it at "pChecker". The new
 *  CertChainChecker uses the CheckCallback pointed to by "callback" as its
 *  callback function. It uses the Object pointed to by "initialState" (if
 *  any) as its initial state. As noted above, the initial state Object must
 *  provide a custom implementation of PKIX_PL_Object_Duplicate if the
 *  CertChainChecker is to be used during certificate chain building.
 *
 *  A CertChainChecker may be presented certificates in the "reverse"
 *  direction (from trust anchor to target) or in the "forward" direction
 *  (from target to trust anchor). All CertChainCheckers must support
 *  "reverse checking", while support for "forward checking" is optional. The
 *  CertChainChecker is initialized with two Boolean flags that deal with this
 *  distinction: "forwardCheckingSupported" and "forwardDirectionExpected".
 *  If the "forwardCheckingSupported" Boolean flag is TRUE, it indicates that
 *  this CertChainChecker is capable of checking certificates in the "forward"
 *  direction (as well as the "reverse" direction, which all CertChainCheckers
 *  MUST support). The "forwardDirectionExpected" Boolean flag indicates in
 *  which direction the CertChainChecker should expect the certificates to be
 *  presented. This is particularly useful for CertChainCheckers that are
 *  capable of checking in either the "forward" direction or the "reverse"
 *  direction, but have different processing steps depending on the direction.
 *
 *  The CertChainChecker also uses the List of OIDs pointed to by "extensions"
 *  as the supported certificate extensions. All certificate extensions that
 *  the CertChainChecker might possibly recognize and be able to process
 *  should be included in the List of supported extensions. If "checker" does
 *  not recognize or process any certificate extensions, "extensions" should
 *  be set to NULL.
 *
 * PARAMETERS:
 *  "callback"
 *      The CheckCallback function to be used. Must be non-NULL.
 *  "forwardCheckingSupported"
 *      A Boolean value indicating whether or not this CertChainChecker is
 *      capable of checking certificates in the "forward" direction.
 *  "forwardDirectionExpected"
 *      A Boolean value indicating whether or not this CertChainChecker should
 *      be used to check in the "forward" direction.
 *  "extensions"
 *      Address of List of OIDs representing the supported extensions.
 *  "initialState"
 *      Address of Object representing the CertChainChecker's initial state
 *      (if any).
 *  "pChecker"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_Create(
    PKIX_CertChainChecker_CheckCallback callback,
    PKIX_Boolean forwardCheckingSupported,
    PKIX_Boolean forwardDirectionExpected,
    PKIX_List *extensions,  /* list of PKIX_PL_OID */
    PKIX_PL_Object *initialState,
    PKIX_CertChainChecker **pChecker,
    void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_GetCheckCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "checker's" Check callback function and puts it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "checker"
 *      The CertChainChecker whose Check callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where Check callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_GetCheckCallback(
        PKIX_CertChainChecker *checker,
        PKIX_CertChainChecker_CheckCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_IsForwardCheckingSupported
 * DESCRIPTION:
 *
 *  Checks whether forward checking is supported by the CertChainChecker
 *  pointed to by "checker" and stores the Boolean result at
 *  "pForwardCheckingSupported".
 *
 *  A CertChainChecker may be presented certificates in the "reverse"
 *  direction (from trust anchor to target) or in the "forward" direction
 *  (from target to trust anchor). All CertChainCheckers must support
 *  "reverse checking", while support for "forward checking" is optional. This
 *  function is used to determine whether forward checking is supported.
 *
 * PARAMETERS:
 *  "checker"
 *      The CertChainChecker whose ability to validate certificates in the
 *      "forward" direction is to be checked. Must be non-NULL.
 *  "pForwardCheckingSupported"
 *      Destination of the Boolean result. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_IsForwardCheckingSupported(
        PKIX_CertChainChecker *checker,
        PKIX_Boolean *pForwardCheckingSupported,
        void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_IsForwardDirectionExpected
 * DESCRIPTION:
 *
 *  Checks whether the CertChainChecker pointed to by "checker" has been
 *  initialized to expect the certificates to be presented in the "forward"
 *  direction and stores the Boolean result at "pForwardDirectionExpected".
 *
 *  A CertChainChecker may be presented certificates in the "reverse"
 *  direction (from trust anchor to target) or in the "forward" direction
 *  (from target to trust anchor). All CertChainCheckers must support
 *  "reverse checking", while support for "forward checking" is optional. This
 *  function is used to determine in which direction the CertChainChecker
 *  expects the certificates to be presented.
 *
 * PARAMETERS:
 *  "checker"
 *      The CertChainChecker that has been initialized to expect certificates
 *      in either the "forward" or "reverse" directions. Must be non-NULL.
 *  "pForwardDirectionExpected"
 *      Destination of the Boolean result. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_IsForwardDirectionExpected(
        PKIX_CertChainChecker *checker,
        PKIX_Boolean *pForwardDirectionExpected,
        void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_GetSupportedExtensions
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a List of OIDs (each OID corresponding to a
 *  certificate extension supported by the CertChainChecker pointed to by
 *  "checker") and stores it at "pExtensions". All certificate extensions that
 *  the CertChainChecker might possibly recognize and be able to process
 *  should be included in the List of supported extensions. If "checker" does
 *  not recognize or process any certificate extensions, this function stores
 *  NULL at "pExtensions".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "checker"
 *      Address of CertChainChecker whose supported extension OIDs are to be
 *      stored. Must be non-NULL.
 *  "pExtensions"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_GetSupportedExtensions(
        PKIX_CertChainChecker *checker,
        PKIX_List **pExtensions, /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_GetCertChainCheckerState
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a PKIX_PL_Object representing the internal state
 *  (if any) of the CertChainChecker pointed to by "checker" and stores it at
 *  "pCertChainCheckerState".
 *
 * PARAMETERS:
 *  "checker"
 *      Address of CertChainChecker whose state is to be stored.
 *      Must be non-NULL.
 *  "pCertChainCheckerState"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_GetCertChainCheckerState(
        PKIX_CertChainChecker *checker,
        PKIX_PL_Object **pCertChainCheckerState,
        void *plContext);

/*
 * FUNCTION: PKIX_CertChainChecker_SetCertChainCheckerState
 * DESCRIPTION:
 *
 *  Sets the internal state of the CertChainChecker pointed to by "checker"
 *  using the Object pointed to by "certChainCheckerState". If "checker" needs
 *  a NULL internal state, "certChainCheckerState" should be set to NULL.
 *
 * PARAMETERS:
 *  "checker"
 *      Address of CertChainChecker whose state is to be set. Must be non-NULL.
 *  "certChainCheckerState"
 *      Address of Object representing internal state.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "checker"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertChainChecker Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertChainChecker_SetCertChainCheckerState(
        PKIX_CertChainChecker *checker,
        PKIX_PL_Object *certChainCheckerState,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CHECKER_H */
