/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the various parameters used
 * by the top-level functions.
 *
 */

#ifndef _PKIX_PARAMS_H
#define _PKIX_PARAMS_H

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

/* PKIX_ProcessingParams
 *
 * PKIX_ProcessingParams are parameters used when validating or building a
 * chain of certificates. Using the parameters, the caller can specify several
 * things, including the various inputs to the PKIX chain validation
 * algorithm (such as trust anchors, initial policies, etc), any customized
 * functionality (such as CertChainCheckers, RevocationCheckers, CertStores),
 * and whether revocation checking should be disabled.
 *
 * Once the caller has created the ProcessingParams object, the caller then
 * passes it to PKIX_ValidateChain or PKIX_BuildChain, which uses it to call
 * the user's callback functions as needed during the validation or building
 * process.
 *
 * If a parameter is not set (or is set to NULL), it will be set to the
 * default value for that parameter. The default value for the Date parameter
 * is NULL, which indicates the current time when the path is validated. The
 * default for the remaining parameters is the least constrained.
 */

/*
 * FUNCTION: PKIX_ProcessingParams_Create
 * DESCRIPTION:
 *
 *  Creates a new ProcessingParams object. Trust anchor list is set to
 *  newly created empty list of trust. In this case trust anchors will
 *  be taken from provided cert store. Pointed to the created
 *  ProcessingParams object is stored in "pParams".
 *
 * PARAMETERS:
 *  "anchors"
 *      Address of List of (non-empty) TrustAnchors to be used.
 *      Must be non-NULL.
 *  "pParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_Create(
        PKIX_ProcessingParams **pParams,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetCertChainCheckers
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of CertChainCheckers (if any) that are set
 *  in the ProcessingParams pointed to by "params" and stores it at
 *  "pCheckers". Each CertChainChecker represents a custom certificate
 *  validation check used by PKIX_ValidateChain or PKIX_BuildChain as needed
 *  during the validation or building process. If "params" does not have any
 *  CertChainCheckers, this function stores an empty List at "pCheckers".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of CertChainCheckers (if any)
 *      are to be stored. Must be non-NULL.
 *  "pCheckers"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetCertChainCheckers(
        PKIX_ProcessingParams *params,
        PKIX_List **pCheckers, /* list of PKIX_CertChainChecker */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetCertChainCheckers
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a List of
 *  CertChainCheckers pointed to by "checkers". Each CertChainChecker
 *  represents a custom certificate validation check used by
 *  PKIX_ValidateChain or PKIX_BuildChain as needed during the validation or
 *  building process. If "checkers" is NULL, no CertChainCheckers will be used.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of CertChainCheckers is to be
 *      set. Must be non-NULL.
 *  "checkers"
 *      Address of List of CertChainCheckers to be set. If NULL, no
 *      CertChainCheckers will be used.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params" and "checkers"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetCertChainCheckers(
        PKIX_ProcessingParams *params,
        PKIX_List *checkers,  /* list of PKIX_CertChainChecker */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_AddCertChainChecker
 * DESCRIPTION:
 *
 *  Adds the CertChainChecker pointed to by "checker" to the ProcessingParams
 *  pointed to by "params". The CertChainChecker represents a custom
 *  certificate validation check used by PKIX_ValidateChain or PKIX_BuildChain
 *  as needed during the validation or building process.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams to be added to. Must be non-NULL.
 *  "checker"
 *      Address of CertChainChecker to be added. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_AddCertChainChecker(
        PKIX_ProcessingParams *params,
        PKIX_CertChainChecker *checker,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetRevocationChecker
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the RevocationChecker that are set
 *  in the ProcessingParams pointed to by "params" and stores it at
 *  "pRevChecker". Each RevocationChecker represents a revocation
 *  check used by PKIX_ValidateChain or PKIX_BuildChain as needed during the
 *  validation or building process. If "params" does not have any
 *  RevocationCheckers, this function stores an empty List at "pRevChecker".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of RevocationCheckers
 *      is to be stored. Must be non-NULL.
 *  "pRevChecker"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetRevocationChecker(
        PKIX_ProcessingParams *params,
        PKIX_RevocationChecker **pChecker,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetRevocationChecker
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a 
 *  RevocationChecker pointed to by "revChecker". Revocation
 *  checker object should be created and assigned to processing
 *  parameters before chain build or validation can begin.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of RevocationCheckers is to be
 *      set. Must be non-NULL.
 *  "revChecker"
 *      Address of RevocationChecker to be set. Must be set before chain
 *      building or validation.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetRevocationChecker(
        PKIX_ProcessingParams *params,
        PKIX_RevocationChecker *revChecker,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetCertStores
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of CertStores (if any) that are set in the
 *  ProcessingParams pointed to by "params" and stores it at "pStores". Each
 *  CertStore represents a particular repository from which certificates and
 *  CRLs can be retrieved by PKIX_ValidateChain or PKIX_BuildChain as needed
 *  during the validation or building process. If "params" does not have any
 *  CertStores, this function stores an empty List at "pStores".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of CertStores (if any) are to
 *      be stored. Must be non-NULL.
 *  "pStores"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetCertStores(
        PKIX_ProcessingParams *params,
        PKIX_List **pStores,  /* list of PKIX_CertStore */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetCertStores
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a List of CertStores
 *  pointed to by "stores". Each CertStore represents a particular repository
 *  from which certificates and CRLs can be retrieved by PKIX_ValidateChain or
 *  PKIX_BuildChain as needed during the validation or building process. If
 *  "stores" is NULL, no CertStores will be used.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of CertStores is to be set.
 *      Must be non-NULL.
 *  "stores"
 *      Address of List of CertStores to be set. If NULL, no CertStores will
 *      be used.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetCertStores(
        PKIX_ProcessingParams *params,
        PKIX_List *stores,  /* list of PKIX_CertStore */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_AddCertStore
 * DESCRIPTION:
 *
 *  Adds the CertStore pointed to by "store" to the ProcessingParams pointed
 *  to by "params". The CertStore represents a particular repository from
 *  which certificates and CRLs can be retrieved by PKIX_ValidateChain or
 *  PKIX_BuildChain as needed during the validation or building process.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams to be added to. Must be non-NULL.
 *  "store"
 *      Address of CertStore to be added.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_AddCertStore(
        PKIX_ProcessingParams *params,
        PKIX_CertStore *store,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetDate
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Date (if any) that is set in the
 *  ProcessingParams pointed to by "params" and stores it at "pDate". The
 *  Date represents the time for which the validation of the certificate chain
 *  should be determined. If "params" does not have any Date set, this function
 *  stores NULL at "pDate".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose Date (if any) is to be stored.
 *      Must be non-NULL.
 *  "pDate"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetDate(
        PKIX_ProcessingParams *params,
        PKIX_PL_Date **pDate,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetDate
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a Date pointed to by
 *  "date". The Date represents the time for which the validation of the
 *  certificate chain should be determined. If "date" is NULL, the current
 *  time is used during validation.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose Date is to be set. Must be non-NULL.
 *  "date"
 *      Address of Date to be set. If NULL, current time is used.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetDate(
        PKIX_ProcessingParams *params,
        PKIX_PL_Date *date,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetInitialPolicies
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of OIDs (if any) that are set in the
 *  ProcessingParams pointed to by "params" and stores it at "pInitPolicies".
 *  Each OID represents an initial policy identifier, indicating that any
 *  one of these policies would be acceptable to the certificate user for
 *  the purposes of certification path processing. If "params" does not have
 *  any initial policies, this function stores an empty List at
 *  "pInitPolicies".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of OIDs (if any) are to be
 *      stored. Must be non-NULL.
 *  "pInitPolicies"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetInitialPolicies(
        PKIX_ProcessingParams *params,
        PKIX_List **pInitPolicies,    /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetInitialPolicies
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a List of OIDs
 *  pointed to by "initPolicies".
 *
 *  Each OID represents an initial policy identifier, indicating that any
 *  one of these policies would be acceptable to the certificate user for
 *  the purposes of certification path processing. By default, any policy
 *  is acceptable (i.e. all policies), so a user that wants to allow any
 *  policy as acceptable does not need to call this method. Similarly, if
 *  initPolicies is NULL or points to an empty List, all policies are
 *  acceptable.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of OIDs is to be set.
 *      Must be non-NULL.
 *  "initPolicies"
 *      Address of List of OIDs to be set. If NULL or if pointing to an empty
 *      List, all policies are acceptable.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetInitialPolicies(
        PKIX_ProcessingParams *params,
        PKIX_List *initPolicies,    /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetPolicyQualifiersRejected
 * DESCRIPTION:
 *
 *  Checks whether the ProcessingParams pointed to by "params" indicate that
 *  policy qualifiers should be rejected and stores the Boolean result at
 *  "pRejected".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams used to determine whether or not policy
 *      qualifiers should be rejected. Must be non-NULL.
 *  "pRejected"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetPolicyQualifiersRejected(
        PKIX_ProcessingParams *params,
        PKIX_Boolean *pRejected,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetPolicyQualifiersRejected
 * DESCRIPTION:
 *
 *  Specifies in the ProcessingParams pointed to by "params" whether policy
 *  qualifiers are rejected using the Boolean value of "rejected".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams to be set. Must be non-NULL.
 *  "rejected"
 *      Boolean value indicating whether policy qualifiers are to be rejected.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetPolicyQualifiersRejected(
        PKIX_ProcessingParams *params,
        PKIX_Boolean rejected,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetTargetCertConstraints
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the CertSelector (if any) that is set in the
 *  ProcessingParams pointed to by "params" and stores it at "pConstraints".
 *  The CertSelector represents the constraints to be placed on the target
 *  certificate. If "params" does not have any CertSelector set, this function
 *  stores NULL at "pConstraints".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose CertSelector (if any) is to be
 *      stored. Must be non-NULL.
 *  "pConstraints"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetTargetCertConstraints(
        PKIX_ProcessingParams *params,
        PKIX_CertSelector **pConstraints,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetTargetCertConstraints
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a CertSelector
 *  pointed to by "constraints". The CertSelector represents the constraints
 *  to be placed on the target certificate. If "constraints" is NULL, no
 *  constraints are defined.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose CertSelector is to be set.
 *      Must be non-NULL.
 *  "constraints"
 *      Address of CertSelector to be set. If NULL, no constraints are defined.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetTargetCertConstraints(
        PKIX_ProcessingParams *params,
        PKIX_CertSelector *constraints,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetTrustAnchors
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of TrustAnchors that are set in
 *  the ProcessingParams pointed to by "params" and stores it at "pAnchors".
 *  If the function succeeds, the pointer to the List is guaranteed to be
 *  non-NULL and the List is guaranteed to be non-empty.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of TrustAnchors are to
 *      be stored. Must be non-NULL.
 *  "pAnchors"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetTrustAnchors(
        PKIX_ProcessingParams *params,
        PKIX_List **pAnchors,  /* list of TrustAnchor */
        void *plContext);
/*
 * FUNCTION: PKIX_ProcessingParams_SetTrustAnchors
 * DESCRIPTION:
 *
 * Sets user defined set of trust anchors. The handling of the trust anchors
 * may be furthered alter via PKIX_ProcessingParams_SetUseOnlyTrustAnchors.
 * By default, a certificate will be considered invalid if it does not chain
 * to a trusted anchor from this list.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of TrustAnchors are to
 *      be stored. Must be non-NULL.
 *  "anchors"
 *      Address of the trust anchors list object. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetTrustAnchors(
        PKIX_ProcessingParams *params,
        PKIX_List *pAnchors,  /* list of TrustAnchor */
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetUseOnlyTrustAnchors
 * DESCRIPTION:
 *
 * Retrieves a pointer to the Boolean. The boolean value represents
 * the switch value that is used to identify whether trust anchors, if
 * specified, should be the exclusive source of trust information.
 * If the function succeeds, the pointer to the Boolean is guaranteed to be
 * non-NULL.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams. Must be non-NULL.
 *  "pUseOnlyTrustAnchors"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetUseOnlyTrustAnchors(
        PKIX_ProcessingParams *params,
        PKIX_Boolean *pUseOnlyTrustAnchors,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetUseOnlyTrustAnchors
 * DESCRIPTION:
 *
 * Configures whether trust anchors are used as the exclusive source of trust.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams. Must be non-NULL.
 *  "useOnlyTrustAnchors"
 *      If true, indicates that trust anchors should be used exclusively when
 *      they have been specified via PKIX_ProcessingParams_SetTrustAnchors. A
 *      certificate will be considered invalid if it does not chain to a
 *      trusted anchor from that list.
 *      If false, indicates that the trust anchors are additive to whatever
 *      existing trust stores are configured. A certificate is considered
 *      valid if it chains to EITHER a trusted anchor from that list OR a
 *      certificate marked trusted in a trust store.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetUseOnlyTrustAnchors(
        PKIX_ProcessingParams *params,
        PKIX_Boolean useOnlyTrustAnchors,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetUseAIAForCertFetching
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Boolean. The boolean value represents
 *  the switch value that is used to identify if url in cert AIA extension
 *  may be used for cert fetching.
 *  If the function succeeds, the pointer to the Boolean is guaranteed to be
 *  non-NULL.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams. Must be non-NULL.
 *  "pUseAIA"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetUseAIAForCertFetching(
        PKIX_ProcessingParams *params,
        PKIX_Boolean *pUseAIA,  /* list of TrustAnchor */
        void *plContext);
/*
 * FUNCTION: PKIX_ProcessingParams_SetTrustAnchors
 * DESCRIPTION:
 *
 * Sets switch value that defines if url in cert AIA extension
 * may be used for cert fetching.
 * 
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams.
 *  "useAIA"
 *      Address of the trust anchors list object. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetUseAIAForCertFetching(
        PKIX_ProcessingParams *params,
        PKIX_Boolean useAIA,  
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetQualifyTargetCert
 * DESCRIPTION:
 *
 * Sets a boolean value that tells if libpkix needs to check that
 * the target certificate satisfies the conditions set in processing
 * parameters. Includes but not limited to date, ku and eku checks.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of TrustAnchors are to
 *      be stored. Must be non-NULL.
 *  "qualifyTargetCert"
 *      boolean value if set to true will trigger qualification of the
 *      target certificate.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetQualifyTargetCert(
        PKIX_ProcessingParams *params,
        PKIX_Boolean qualifyTargetCert,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetHintCerts
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a List of Certs supplied by the user as a suggested
 *  partial CertChain (subject to verification), that are set in the
 *  ProcessingParams pointed to by "params", and stores it at "pHintCerts".
 *  The List returned may be empty or NULL.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of TrustAnchors are to
 *      be stored. Must be non-NULL.
 *  "pHintCerts"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetHintCerts(
        PKIX_ProcessingParams *params,
        PKIX_List **pHintCerts,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetHintCerts
 * DESCRIPTION:
 *
 *  Stores a pointer to a List of Certs supplied by the user as a suggested
 *  partial CertChain (subject to verification), as an element in the
 *  ProcessingParams pointed to by "params". The List may be empty or NULL.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose List of HintCerts is to be stored.
 *      Must be non-NULL.
 *  "hintCerts"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetHintCerts(
        PKIX_ProcessingParams *params,
        PKIX_List *hintCerts,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_GetResourceLimits
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ResourceLimits (if any) that is set in the
 *  ProcessingParams pointed to by "params" and stores it at "pResourceLimits".
 *  The ResourceLimits represent the maximum resource usage that the caller 
 *  desires (such as MaxTime). The ValidateChain or BuildChain call will not
 *  exceed these maximum limits. If "params" does not have any ResourceLimits
 *  set, this function stores NULL at "pResourceLimits".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose ResourceLimits (if any) are to be
 *      stored. Must be non-NULL.
 *  "pResourceLimits"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_GetResourceLimits(
        PKIX_ProcessingParams *params,
        PKIX_ResourceLimits **pResourceLimits,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetResourceLimits
 * DESCRIPTION:
 *
 *  Sets the ProcessingParams pointed to by "params" with a ResourceLimits
 *  object pointed to by "resourceLimits". The ResourceLimits represent the
 *  maximum resource usage that the caller desires (such as MaxTime). The
 *  ValidateChain or BuildChain call will not exceed these maximum limits.
 *  If "resourceLimits" is NULL, no ResourceLimits are defined.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams whose ResourceLimits are to be set.
 *      Must be non-NULL.
 *  "resourceLimits"
 *      Address of ResourceLimits to be set. If NULL, no limits are defined.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetResourceLimits(
        PKIX_ProcessingParams *params,
        PKIX_ResourceLimits *resourceLimits,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_IsAnyPolicyInhibited
 * DESCRIPTION:
 *
 *  Checks whether the ProcessingParams pointed to by "params" indicate that
 *  anyPolicy is inhibited and stores the Boolean result at "pInhibited".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams used to determine whether or not anyPolicy
 *      inhibited. Must be non-NULL.
 *  "pInhibited"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_IsAnyPolicyInhibited(
        PKIX_ProcessingParams *params,
        PKIX_Boolean *pInhibited,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetAnyPolicyInhibited
 * DESCRIPTION:
 *
 *  Specifies in the ProcessingParams pointed to by "params" whether anyPolicy
 *  is inhibited using the Boolean value of "inhibited".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams to be set. Must be non-NULL.
 *  "inhibited"
 *      Boolean value indicating whether anyPolicy is to be inhibited.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetAnyPolicyInhibited(
        PKIX_ProcessingParams *params,
        PKIX_Boolean inhibited,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_IsExplicitPolicyRequired
 * DESCRIPTION:
 *
 *  Checks whether the ProcessingParams pointed to by "params" indicate that
 *  explicit policies are required and stores the Boolean result at
 *  "pRequired".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams used to determine whether or not explicit
 *      policies are required. Must be non-NULL.
 *  "pRequired"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_IsExplicitPolicyRequired(
        PKIX_ProcessingParams *params,
        PKIX_Boolean *pRequired,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetExplicitPolicyRequired
 * DESCRIPTION:
 *
 *  Specifies in the ProcessingParams pointed to by "params" whether explicit
 *  policies are required using the Boolean value of "required".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams to be set. Must be non-NULL.
 *  "required"
 *      Boolean value indicating whether explicit policies are to be required.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetExplicitPolicyRequired(
        PKIX_ProcessingParams *params,
        PKIX_Boolean required,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_IsPolicyMappingInhibited
 * DESCRIPTION:
 *
 *  Checks whether the ProcessingParams pointed to by "params" indicate that
 *  policyMapping is inhibited and stores the Boolean result at "pInhibited".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams used to determine whether or not policy
 *      mappings are inhibited. Must be non-NULL.
 *  "pInhibited"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_IsPolicyMappingInhibited(
        PKIX_ProcessingParams *params,
        PKIX_Boolean *pInhibited,
        void *plContext);

/*
 * FUNCTION: PKIX_ProcessingParams_SetPolicyMappingInhibited
 * DESCRIPTION:
 *
 *  Specifies in the ProcessingParams pointed to by "params" whether policy
 *  mapping is inhibited using the Boolean value of "inhibited".
 *
 * PARAMETERS:
 *  "params"
 *      Address of ProcessingParams to be set. Must be non-NULL.
 *  "inhibited"
 *      Boolean value indicating whether policy mapping is to be inhibited.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ProcessingParams_SetPolicyMappingInhibited(
        PKIX_ProcessingParams *params,
        PKIX_Boolean inhibited,
        void *plContext);


/* PKIX_ValidateParams
 *
 * PKIX_ValidateParams consists of a ProcessingParams object as well as the
 * List of Certs (certChain) that the caller is trying to validate.
 */

/*
 * FUNCTION: PKIX_ValidateParams_Create
 * DESCRIPTION:
 *
 *  Creates a new ValidateParams object and stores it at "pParams".
 *
 * PARAMETERS:
 *  "procParams"
 *      Address of ProcessingParams to be used. Must be non-NULL.
 *  "chain"
 *      Address of List of Certs (certChain) to be validated. Must be non-NULL.
 *  "pParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateParams_Create(
        PKIX_ProcessingParams *procParams,
        PKIX_List *chain,
        PKIX_ValidateParams **pParams,
        void *plContext);

/*
 * FUNCTION: PKIX_ValidateParams_GetProcessingParams
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ProcessingParams that represent the basic
 *  certificate processing parameters used during chain validation and chain
 *  building from the ValidateParams pointed to by "valParams" and stores it
 *  at "pProcParams". If the function succeeds, the pointer to the
 *  ProcessingParams is guaranteed to be non-NULL.
 *
 * PARAMETERS:
 *  "valParams"
 *      Address of ValidateParams whose ProcessingParams are to be stored.
 *      Must be non-NULL.
 *  "pProcParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateParams_GetProcessingParams(
        PKIX_ValidateParams *valParams,
        PKIX_ProcessingParams **pProcParams,
        void *plContext);

/*
 * FUNCTION: PKIX_ValidateParams_GetCertChain
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of Certs (certChain) that is set in the
 *  ValidateParams pointed to by "valParams" and stores it at "pChain". If the
 *  function succeeds, the pointer to the CertChain is guaranteed to be
 *  non-NULL.
 *
 * PARAMETERS:
 *  "valParams"
 *      Address of ValidateParams whose CertChain is to be stored.
 *      Must be non-NULL.
 *  "pChain"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateParams_GetCertChain(
        PKIX_ValidateParams *valParams,
        PKIX_List **pChain,
        void *plContext);

/* PKIX_TrustAnchor
 *
 * A PKIX_TrustAnchor represents a trusted entity and can be specified using a
 * self-signed certificate or using the trusted CA's name and public key. In
 * order to limit the trust in the trusted entity, name constraints can also
 * be imposed on the trust anchor.
 */

/*
 * FUNCTION: PKIX_TrustAnchor_CreateWithCert
 * DESCRIPTION:
 *
 *  Creates a new TrustAnchor object using the Cert pointed to by "cert" as
 *  the trusted certificate and stores it at "pAnchor". Once created, a
 *  TrustAnchor is immutable.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert to use as trusted certificate. Must be non-NULL.
 *  "pAnchor"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_TrustAnchor_CreateWithCert(
        PKIX_PL_Cert *cert,
        PKIX_TrustAnchor **pAnchor,
        void *plContext);

/*
 * FUNCTION: PKIX_TrustAnchor_CreateWithNameKeyPair
 * DESCRIPTION:
 *
 *  Creates a new TrustAnchor object using the X500Name pointed to by "name",
 *  and the PublicKey pointed to by "pubKey" and stores it at "pAnchor". The
 *  CertNameConstraints pointed to by "nameConstraints" (if any) are used to
 *  limit the trust placed in this trust anchor. To indicate that name
 *  constraints don't apply, set "nameConstraints" to NULL. Once created, a
 *  TrustAnchor is immutable.
 *
 * PARAMETERS:
 *  "name"
 *      Address of X500Name to use as name of trusted CA. Must be non-NULL.
 *  "pubKey"
 *      Address of PublicKey to use as trusted public key. Must be non-NULL.
 *  "nameConstraints"
 *      Address of CertNameConstraints to use as initial name constraints.
 *      If NULL, no name constraints are applied.
 *  "pAnchor"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_TrustAnchor_CreateWithNameKeyPair(
        PKIX_PL_X500Name *name,
        PKIX_PL_PublicKey *pubKey,
        PKIX_PL_CertNameConstraints *nameConstraints,
        PKIX_TrustAnchor **pAnchor,
        void *plContext);

/*
 * FUNCTION: PKIX_TrustAnchor_GetTrustedCert
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Cert that is set in the TrustAnchor pointed to
 *  by "anchor" and stores it at "pCert". If "anchor" does not have a Cert
 *  set, this function stores NULL at "pCert".
 *
 * PARAMETERS:
 *  "anchor"
 *      Address of TrustAnchor whose Cert is to be stored. Must be non-NULL.
 *  "pChain"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_TrustAnchor_GetTrustedCert(
        PKIX_TrustAnchor *anchor,
        PKIX_PL_Cert **pCert,
        void *plContext);

/*
 * FUNCTION: PKIX_TrustAnchor_GetCAName
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the CA's X500Name (if any) that is set in the
 *  TrustAnchor pointed to by "anchor" and stores it at "pCAName". If "anchor"
 *  does not have an X500Name set, this function stores NULL at "pCAName".
 *
 * PARAMETERS:
 *  "anchor"
 *      Address of TrustAnchor whose CA Name is to be stored. Must be non-NULL.
 *  "pCAName"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_TrustAnchor_GetCAName(
        PKIX_TrustAnchor *anchor,
        PKIX_PL_X500Name **pCAName,
        void *plContext);

/*
 * FUNCTION: PKIX_TrustAnchor_GetCAPublicKey
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the CA's PublicKey (if any) that is set in the
 *  TrustAnchor pointed to by "anchor" and stores it at "pPubKey". If "anchor"
 *  does not have a PublicKey set, this function stores NULL at "pPubKey".
 *
 * PARAMETERS:
 *  "anchor"
 *      Address of TrustAnchor whose CA PublicKey is to be stored.
 *      Must be non-NULL.
 *  "pPubKey"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_TrustAnchor_GetCAPublicKey(
        PKIX_TrustAnchor *anchor,
        PKIX_PL_PublicKey **pPubKey,
        void *plContext);

/*
 * FUNCTION: PKIX_TrustAnchor_GetNameConstraints
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the CertNameConstraints (if any) set in the
 *  TrustAnchor pointed to by "anchor" and stores it at "pConstraints". If
 *  "anchor" does not have any CertNameConstraints set, this function stores
 *  NULL at "pConstraints".
 *
 * PARAMETERS:
 *  "anchor"
 *      Address of TrustAnchor whose CertNameConstraints are to be stored.
 *      Must be non-NULL.
 *  "pConstraints"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Params Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_TrustAnchor_GetNameConstraints(
        PKIX_TrustAnchor *anchor,
        PKIX_PL_CertNameConstraints **pNameConstraints,
        void *plContext);

/* PKIX_ResourceLimits
 *
 *  A PKIX_ResourceLimits object represents the maximum resource usage that
 *  the caller desires. The ValidateChain or BuildChain call
 *  will not exceed these maximum limits. For example, the caller may want
 *  a timeout value of 1 minute, meaning that if the ValidateChain or
 *  BuildChain function is unable to finish in 1 minute, it should abort
 *  with an Error.
 */

/*
 * FUNCTION: PKIX_ResourceLimits_Create
 * DESCRIPTION:
 *
 *  Creates a new ResourceLimits object and stores it at "pResourceLimits".
 *
 * PARAMETERS:
 *  "pResourceLimits"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_Create(
        PKIX_ResourceLimits **pResourceLimits,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_GetMaxTime
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the maximum time that is
 *  set in the ResourceLimits object pointed to by "resourceLimits" and stores
 *  it at "pMaxTime". This maximum time (in seconds) should not be exceeded
 *  by the function whose ProcessingParams contain this ResourceLimits object
 *  (typically ValidateChain or BuildChain). It essentially functions as a
 *  time-out value and is only appropriate if blocking I/O is being used.
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum time (in seconds) is
 *      to be stored. Must be non-NULL.
 *  "pMaxTime"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_GetMaxTime(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 *pMaxTime,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_SetMaxTime
 * DESCRIPTION:
 *
 *  Sets the maximum time of the ResourceLimits object pointed to by
 *  "resourceLimits" using the PKIX_UInt32 value of "maxTime". This
 *  maximum time (in seconds) should not be exceeded by the function
 *  whose ProcessingParams contain this ResourceLimits object
 *  (typically ValidateChain or BuildChain). It essentially functions as a
 *  time-out value and is only appropriate if blocking I/O is being used.
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum time (in seconds) is
 *      to be set. Must be non-NULL.
 *  "maxTime"
 *      Value of PKIX_UInt32 representing the maximum time (in seconds)
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_SetMaxTime(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 maxTime,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_GetMaxFanout
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the maximum fanout that is
 *  set in the ResourceLimits object pointed to by "resourceLimits" and stores
 *  it at "pMaxFanout". This maximum fanout (number of certs) should not be
 *  exceeded by the function whose ProcessingParams contain this ResourceLimits
 *  object (typically ValidateChain or BuildChain). If the builder encounters
 *  more than this maximum number of certificates when searching for the next
 *  candidate certificate, it should abort and return an error. This
 *  parameter is only relevant for ValidateChain if it needs to internally call
 *  BuildChain (e.g. in order to build the chain to a CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum fanout (number of certs)
 *      is to be stored. Must be non-NULL.
 *  "pMaxFanout"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_GetMaxFanout(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 *pMaxFanout,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_SetMaxFanout
 * DESCRIPTION:
 *
 *  Sets the maximum fanout of the ResourceLimits object pointed to by
 *  "resourceLimits" using the PKIX_UInt32 value of "maxFanout". This maximum
 *  fanout (number of certs) should not be exceeded by the function whose
 *  ProcessingParams contain this ResourceLimits object (typically ValidateChain
 *  or BuildChain). If the builder encounters more than this maximum number of
 *  certificates when searching for the next candidate certificate, it should
 *  abort and return an Error. This parameter is only relevant for ValidateChain
 *  if it needs to internally call BuildChain (e.g. in order to build the
 *  chain to a CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum fanout (number of certs)
 *      is to be set. Must be non-NULL.
 *  "maxFanout"
 *      Value of PKIX_UInt32 representing the maximum fanout (number of certs)
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_SetMaxFanout(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 maxFanout,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_GetMaxDepth
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the maximum depth that is
 *  set in the ResourceLimits object pointed to by "resourceLimits" and stores
 *  it at "pMaxDepth". This maximum depth (number of certs) should not be
 *  exceeded by the function whose ProcessingParams contain this ResourceLimits
 *  object (typically ValidateChain or BuildChain). If the builder encounters
 *  more than this maximum number of certificates when searching for the next
 *  candidate certificate, it should abort and return an error. This
 *  parameter is only relevant for ValidateChain if it needs to internally call
 *  BuildChain (e.g. in order to build the chain to a CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum depth (number of certs)
 *      is to be stored. Must be non-NULL.
 *  "pMaxDepth"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_GetMaxDepth(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 *pMaxDepth,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_SetMaxDepth
 * DESCRIPTION:
 *
 *  Sets the maximum depth of the ResourceLimits object pointed to by
 *  "resourceLimits" using the PKIX_UInt32 value of "maxDepth". This maximum
 *  depth (number of certs) should not be exceeded by the function whose
 *  ProcessingParams contain this ResourceLimits object (typically ValidateChain
 *  or BuildChain). If the builder encounters more than this maximum number of
 *  certificates when searching for the next candidate certificate, it should
 *  abort and return an Error. This parameter is only relevant for ValidateChain
 *  if it needs to internally call BuildChain (e.g. in order to build the
 *  chain to a CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum depth (number of certs)
 *      is to be set. Must be non-NULL.
 *  "maxDepth"
 *      Value of PKIX_UInt32 representing the maximum depth (number of certs)
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_SetMaxDepth(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 maxDepth,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_GetMaxNumberOfCerts
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the maximum number of traversed
 *  certs that is set in the ResourceLimits object pointed to by "resourceLimits"
 *  and stores it at "pMaxNumber". This maximum number of traversed certs should
 *  not be exceeded by the function whose ProcessingParams contain this ResourceLimits
 *  object (typically ValidateChain or BuildChain). If the builder traverses more
 *  than this number of certs during the build process, it should abort and
 *  return an Error. This parameter is only relevant for ValidateChain if it
 *  needs to internally call BuildChain (e.g. in order to build the chain to a
 *  CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum number of traversed certs
 *      is to be stored. Must be non-NULL.
 *  "pMaxNumber"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_GetMaxNumberOfCerts(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 *pMaxNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_SetMaxNumberOfCerts
 * DESCRIPTION:
 *
 *  Sets the maximum number of traversed certs of the ResourceLimits object
 *  pointed to by "resourceLimits" using the PKIX_UInt32 value of "maxNumber".
 *  This maximum number of traversed certs should not be exceeded by the function 
 *  whose ProcessingParams contain this ResourceLimits object (typically ValidateChain
 *  or BuildChain). If the builder traverses more than this number of certs
 *  during the build process, it should abort and return an Error. This parameter
 *  is only relevant for ValidateChain if it needs to internally call BuildChain
 *  (e.g. in order to build the chain to a CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum number of traversed certs
 *      is to be set. Must be non-NULL.
 *  "maxNumber"
 *      Value of PKIX_UInt32 representing the maximum number of traversed certs
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_SetMaxNumberOfCerts(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 maxNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_GetMaxNumberOfCRLs
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the maximum number of traversed
 *  CRLs that is set in the ResourceLimits object pointed to by "resourceLimits"
 *  and stores it at "pMaxNumber". This maximum number of traversed CRLs should
 *  not be exceeded by the function whose ProcessingParams contain this ResourceLimits
 *  object (typically ValidateChain or BuildChain). If the builder traverses more
 *  than this number of CRLs during the build process, it should abort and
 *  return an Error. This parameter is only relevant for ValidateChain if it
 *  needs to internally call BuildChain (e.g. in order to build the chain to a
 *  CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum number of traversed CRLs
 *      is to be stored. Must be non-NULL.
 *  "pMaxNumber"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_GetMaxNumberOfCRLs(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 *pMaxNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ResourceLimits_SetMaxNumberOfCRLs
 * DESCRIPTION:
 *
 *  Sets the maximum number of traversed CRLs of the ResourceLimits object
 *  pointed to by "resourceLimits" using the PKIX_UInt32 value of "maxNumber".
 *  This maximum number of traversed CRLs should not be exceeded by the function 
 *  whose ProcessingParams contain this ResourceLimits object (typically ValidateChain
 *  or BuildChain). If the builder traverses more than this number of CRLs
 *  during the build process, it should abort and return an Error. This parameter
 *  is only relevant for ValidateChain if it needs to internally call BuildChain
 *  (e.g. in order to build the chain to a CRL's issuer).
 *
 * PARAMETERS:
 *  "resourceLimits"
 *      Address of ResourceLimits object whose maximum number of traversed CRLs
 *      is to be set. Must be non-NULL.
 *  "maxNumber"
 *      Value of PKIX_UInt32 representing the maximum number of traversed CRLs
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a ResourceLimits Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ResourceLimits_SetMaxNumberOfCRLs(
        PKIX_ResourceLimits *resourceLimits,
        PKIX_UInt32 maxNumber,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PARAMS_H */
