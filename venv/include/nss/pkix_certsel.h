/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the PKIX_CertSelector and the
 * PKIX_ComCertSelParams types.
 *
 */

#ifndef _PKIX_CERTSEL_H
#define _PKIX_CERTSEL_H

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

/* PKIX_CertSelector
 *
 * PKIX_CertSelectors provide a standard way for the caller to select
 * certificates based on particular criteria. A CertSelector is typically used
 * by the caller to specify the constraints they wish to impose on the target
 * certificate in a chain. (see pkix_params.h) A CertSelector is also often
 * used to retrieve certificates from a CertStore that match the selector's
 * criteria. (See pkix_certstore.h) For example, the caller may wish to only
 * select those certificates that have a particular Subject Distinguished Name
 * and a particular value for a private certificate extension. The
 * MatchCallback allows the caller to specify the custom matching logic to be
 * used by a CertSelector.
 *
 * By default, the MatchCallback is set to point to the default implementation
 * provided by libpkix, which understands how to process the most common
 * parameters. If the default implementation is used, the caller should set
 * these common parameters using PKIX_CertSelector_SetCommonCertSelectorParams.
 * Any common parameter that is not set is assumed to be disabled, which means
 * the default MatchCallback implementation will select all certificates
 * without regard to that particular disabled parameter. For example, if the
 * SerialNumber parameter is not set, MatchCallback will not filter out any
 * certificate based on its serial number. As such, if no parameters are set,
 * all are disabled and any certificate will match. If a parameter is
 * disabled, its associated PKIX_ComCertSelParams_Get* function returns a
 * default value of NULL, or -1 for PKIX_ComCertSelParams_GetBasicConstraints
 * and PKIX_ComCertSelParams_GetVersion, or 0 for
 * PKIX_ComCertSelParams_GetKeyUsage.
 *
 * If a custom implementation is desired, the default implementation can be
 * overridden by calling PKIX_CertSelector_SetMatchCallback. In this case, the
 * CertSelector can be initialized with a certSelectorContext, which is where
 * the caller can specify the desired parameters the caller wishes to match
 * against. Note that this certSelectorContext must be an Object (although any
 * object type), allowing it to be reference-counted and allowing it to
 * provide the standard Object functions (Equals, Hashcode, ToString, Compare,
 * Duplicate).
 *
 */

/*
 * FUNCTION: PKIX_CertSelector_MatchCallback
 * DESCRIPTION:
 *
 *  This callback function determines whether the specified Cert pointed to by
 *  "cert" matches the criteria of the CertSelector pointed to by "selector".
 *  If the Cert does not matches the CertSelector's criteria, an exception will
 *  be thrown.
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CertSelector whose MatchCallback logic and parameters are
 *      to be used. Must be non-NULL.
 *  "cert"
 *      Address of Cert that is to be matched using "selector".
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertSelector_MatchCallback)(
        PKIX_CertSelector *selector,
        PKIX_PL_Cert *cert,
        void *plContext);

/*
 * FUNCTION: PKIX_CertSelector_Create
 * DESCRIPTION:
 *
 *  Creates a new CertSelector using the Object pointed to by
 *  "certSelectorContext" (if any) and stores it at "pSelector". As noted
 *  above, by default, the MatchCallback is set to point to the default
 *  implementation provided by libpkix, which understands how to process
 *  ComCertSelParams objects. This is overridden if the MatchCallback pointed
 *  to by "callback" is not NULL, in which case the parameters are specified
 *  using the certSelectorContext.
 *
 * PARAMETERS:
 *  "callback"
 *      The MatchCallback function to be used.
 *  "certSelectorContext"
 *      Address of Object representing the CertSelector's context (if any).
 *  "pSelector"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertSelector_Create(
        PKIX_CertSelector_MatchCallback callback,
        PKIX_PL_Object *certSelectorContext,
        PKIX_CertSelector **pSelector,
        void *plContext);

/*
 * FUNCTION: PKIX_CertSelector_GetMatchCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "selector's" Match callback function and puts it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "selector"
 *      The CertSelector whose Match callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where Match callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertSelector_GetMatchCallback(
        PKIX_CertSelector *selector,
        PKIX_CertSelector_MatchCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertSelector_GetCertSelectorContext
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a PKIX_PL_Object representing the context (if any)
 *  of the CertSelector pointed to by "selector" and stores it at
 *  "pCertSelectorContext".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CertSelector whose context is to be stored.
 *      Must be non-NULL.
 *  "pCertSelectorContext"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertSelector_GetCertSelectorContext(
        PKIX_CertSelector *selector,
        PKIX_PL_Object **pCertSelectorContext,
        void *plContext);

/*
 * FUNCTION: PKIX_CertSelector_GetCommonCertSelectorParams
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ComCertSelParams object that represent the
 *  common parameters of the CertSelector pointed to by "selector" and stores
 *  it at "pCommonCertSelectorParams". If there are no common parameters
 *  stored with the CertSelector, this function stores NULL at
 *  "pCommonCertSelectorParams".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CertSelector whose ComCertSelParams object is to be stored.
 *      Must be non-NULL.
 *  "pCommonCertSelectorParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertSelector_GetCommonCertSelectorParams(
        PKIX_CertSelector *selector,
        PKIX_ComCertSelParams **pCommonCertSelectorParams,
        void *plContext);

/*
 * FUNCTION: PKIX_CertSelector_SetCommonCertSelectorParams
 * DESCRIPTION:
 *
 *  Sets the common parameters for the CertSelector pointed to by "selector"
 *  using the ComCertSelParams object pointed to by "commonCertSelectorParams".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CertSelector whose common parameters are to be set.
 *      Must be non-NULL.
 *  "commonCertSelectorParams"
 *      Address of ComCertSelParams object representing the common parameters.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "selector"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertSelector_SetCommonCertSelectorParams(
        PKIX_CertSelector *selector,
        PKIX_ComCertSelParams *commonCertSelectorParams,
        void *plContext);

/* PKIX_ComCertSelParams
 *
 * PKIX_ComCertSelParams objects are X.509 parameters commonly used with
 * CertSelectors, especially when enforcing constraints on a target
 * certificate or determining which certificates to retrieve from a CertStore.
 * ComCertSelParams objects are typically used with those CertSelectors that
 * use the default implementation of MatchCallback, which understands how to
 * process ComCertSelParams objects.
 */

/*
 * FUNCTION: PKIX_ComCertSelParams_Create
 * DESCRIPTION:
 *
 *  Creates a new ComCertSelParams object and stores it at "pParams".
 *
 * PARAMETERS:
 *  "pParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_Create(
        PKIX_ComCertSelParams **pParams,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSubjAltNames
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of GeneralNames (if any) representing the
 *  subject alternative names criterion that is set in the ComCertSelParams
 *  object pointed to by "params" and stores it at "pNames". In order to match
 *  against this criterion, a certificate must contain all or at least one of
 *  the criterion's subject alternative names (depending on the result of
 *  PKIX_ComCertSelParams_GetMatchAllSubjAltNames). The default behavior
 *  requires a certificate to contain all of the criterion's subject
 *  alternative names in order to match.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pNames", in which case all certificates are considered to match this
 *  criterion.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject alternative names
 *      criterion (if any) is to be stored. Must be non-NULL.
 *  "pNames"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSubjAltNames(
        PKIX_ComCertSelParams *params,
        PKIX_List **pNames, /* list of PKIX_PL_GeneralName */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSubjAltNames
 * DESCRIPTION:
 *
 *  Sets the subject alternative names criterion of the ComCertSelParams object
 *  pointed to by "params" using a List of GeneralNames pointed to by "names".
 *  In order to match against this criterion, a certificate must contain all or
 *  at least one of the criterion's subject alternative names (depending on the
 *  result of PKIX_ComCertSelParams_GetMatchAllSubjAltNames). The default
 *  behavior requires a certificate to contain all of the criterion's subject
 *  alternative names in order to match.
 *
 *  If "names" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject alternative
 *      names criterion is to be set. Must be non-NULL.
 *  "names"
 *      Address of List of GeneralNames used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSubjAltNames(
        PKIX_ComCertSelParams *params,
        PKIX_List *names,  /* list of PKIX_PL_GeneralName */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_AddSubjAltName
 * DESCRIPTION:
 *
 *  Adds to the subject alternative names criterion of the ComCertSelParams
 *  object pointed to by "params" using the GeneralName pointed to by "name".
 *  In order to match against this criterion, a certificate must contain all
 *  or at least one of the criterion's subject alternative names (depending on
 *  the result of PKIX_ComCertSelParams_GetMatchAllSubjAltNames). The default
 *  behavior requires a certificate to contain all of the criterion's subject
 *  alternative names in order to match.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject alternative names
 *      criterion is to be added to. Must be non-NULL.
 *  "name"
 *      Address of GeneralName to be added.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_AddSubjAltName(
        PKIX_ComCertSelParams *params,
        PKIX_PL_GeneralName *name,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetPathToNames
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of GeneralNames (if any) representing the
 *  path to names criterion that is set in the ComCertSelParams object pointed
 *  to by "params" and stores it at "pNames". In order to match against this
 *  criterion, a certificate must not include name constraints that would
 *  prohibit building a path to the criterion's specified names.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pNames", in which case all certificates are considered to match this
 *  criterion.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose path to names criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pNames"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetPathToNames(
        PKIX_ComCertSelParams *params,
        PKIX_List **pNames,  /* list of PKIX_PL_GeneralName */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetPathToNames
 * DESCRIPTION:
 *
 *  Sets the path to names criterion of the ComCertSelParams object pointed to
 *  by "params" using a List of GeneralNames pointed to by "names". In order to
 *  match against this criterion, a certificate must not include name
 *  constraints that would prohibit building a path to the criterion's
 *  specified names.
 *
 *  If "names" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose path to names criterion
 *      is to be set. Must be non-NULL.
 *  "names"
 *      Address of List of GeneralNames used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetPathToNames(
        PKIX_ComCertSelParams *params,
        PKIX_List *names,    /* list of PKIX_PL_GeneralName */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_AddPathToName
 * DESCRIPTION:
 *
 *  Adds to the path to names criterion of the ComCertSelParams object pointed
 *  to by "params" using the GeneralName pointed to by "pathToName". In order
 *  to match against this criterion, a certificate must not include name
 *  constraints that would prohibit building a path to the criterion's
 *  specified names.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose path to names criterion is to
 *      be added to. Must be non-NULL.
 *  "pathToName"
 *      Address of GeneralName to be added.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_AddPathToName(
        PKIX_ComCertSelParams *params,
        PKIX_PL_GeneralName *pathToName,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetAuthorityKeyIdentifier
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ByteArray (if any) representing the authority
 *  key identifier criterion that is set in the ComCertSelParams object
 *  pointed to by "params" and stores it at "pAuthKeyId". In order to match
 *  against this criterion, a certificate must contain an
 *  AuthorityKeyIdentifier extension whose value matches the criterion's
 *  authority key identifier value.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pAuthKeyId", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose authority key identifier
 *      criterion (if any) is to be stored. Must be non-NULL.
 *  "pAuthKeyId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetAuthorityKeyIdentifier(
        PKIX_ComCertSelParams *params,
        PKIX_PL_ByteArray **pAuthKeyId,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetAuthorityKeyIdentifier
 * DESCRIPTION:
 *
 *  Sets the authority key identifier criterion of the ComCertSelParams object
 *  pointed to by "params" to the ByteArray pointed to by "authKeyId". In
 *  order to match against this criterion, a certificate must contain an
 *  AuthorityKeyIdentifier extension whose value matches the criterion's
 *  authority key identifier value.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose authority key identifier
 *      criterion is to be set. Must be non-NULL.
 *  "authKeyId"
 *      Address of ByteArray used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetAuthorityKeyIdentifier(
        PKIX_ComCertSelParams *params,
        PKIX_PL_ByteArray *authKeyId,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSubjKeyIdentifier
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ByteArray (if any) representing the subject key
 *  identifier criterion that is set in the ComCertSelParams object pointed to
 *  by "params" and stores it at "pSubjKeyId". In order to match against this
 *  criterion, a certificate must contain a SubjectKeyIdentifier extension
 *  whose value matches the criterion's subject key identifier value.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pSubjKeyId", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject key identifier
 *      criterion (if any) is to be stored. Must be non-NULL.
 *  "pSubjKeyId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSubjKeyIdentifier(
        PKIX_ComCertSelParams *params,
        PKIX_PL_ByteArray **pSubjKeyId,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSubjKeyIdentifier
 * DESCRIPTION:
 *
 *  Sets the subject key identifier criterion of the ComCertSelParams object
 *  pointed to by "params" using a ByteArray pointed to by "subjKeyId". In
 *  order to match against this criterion, a certificate must contain an
 *  SubjectKeyIdentifier extension whose value matches the criterion's subject
 *  key identifier value.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject key identifier
 *      criterion is to be set. Must be non-NULL.
 *  "subjKeyId"
 *      Address of ByteArray used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSubjKeyIdentifier(
        PKIX_ComCertSelParams *params,
        PKIX_PL_ByteArray *subKeyId,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSubjPubKey
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the PublicKey (if any) representing the subject
 *  public key criterion that is set in the ComCertSelParams object pointed to
 *  by "params" and stores it at "pPubKey". In order to match against this
 *  criterion, a certificate must contain a SubjectPublicKey that matches the
 *  criterion's public key.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pPubKey", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject public key criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pPubKey"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSubjPubKey(
        PKIX_ComCertSelParams *params,
        PKIX_PL_PublicKey **pPubKey,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSubjPubKey
 * DESCRIPTION:
 *
 *  Sets the subject public key criterion of the ComCertSelParams object
 *  pointed to by "params" using a PublicKey pointed to by "pubKey". In order
 *  to match against this criterion, a certificate must contain a
 *  SubjectPublicKey that matches the criterion's public key.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject public key
 *      criterion is to be set. Must be non-NULL.
 *  "pubKey"
 *      Address of PublicKey used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSubjPubKey(
        PKIX_ComCertSelParams *params,
        PKIX_PL_PublicKey *pubKey,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSubjPKAlgId
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the OID (if any) representing the subject public key
 *  algorithm identifier criterion that is set in the ComCertSelParams object
 *  pointed to by "params" and stores it at "pPubKey". In order to match
 *  against this criterion, a certificate must contain a SubjectPublicKey with
 *  an algorithm that matches the criterion's algorithm.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pAlgId", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject public key algorithm
 *      identifier (if any) is to be stored. Must be non-NULL.
 *  "pAlgId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSubjPKAlgId(
        PKIX_ComCertSelParams *params,
        PKIX_PL_OID **pAlgId,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSubjPKAlgId
 * DESCRIPTION:
 *
 *  Sets the subject public key algorithm identifier criterion of the
 *  ComCertSelParams object pointed to by "params" using an OID pointed to by
 *  "algId". In order to match against this criterion, a certificate must
 *  contain a SubjectPublicKey with an algorithm that matches the criterion's
 *  algorithm.
 *
 *  If "algId" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject public key
 *      algorithm identifier criterion is to be set. Must be non-NULL.
 *  "algId"
 *      Address of OID used to set criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSubjPKAlgId(
        PKIX_ComCertSelParams *params,
        PKIX_PL_OID *algId,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetBasicConstraints
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the minimum path length (if any) representing the
 *  basic constraints criterion that is set in the ComCertSelParams object
 *  pointed to by "params" and stores it at "pMinPathLength". In order to
 *  match against this criterion, there are several possibilities.
 *
 *  1) If the criterion's minimum path length is greater than or equal to zero,
 *  a certificate must include a BasicConstraints extension with a pathLen of
 *  at least this value.
 *
 *  2) If the criterion's minimum path length is -2, a certificate must be an
 *  end-entity certificate.
 *
 *  3) If the criterion's minimum path length is -1, no basic constraints check
 *  is done and all certificates are considered to match this criterion.
 *
 *  The semantics of other values of the criterion's minimum path length are
 *  undefined but may be defined in future versions of the API.
 *
 *  If "params" does not have this criterion set, this function stores -1 at
 *  "pMinPathLength", in which case all certificates are considered to match
 *  this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose basic constraints criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pMinPathLength"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetBasicConstraints(
        PKIX_ComCertSelParams *params,
        PKIX_Int32 *pMinPathLength,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetBasicConstraints
 * DESCRIPTION:
 *
 *  Sets the basic constraints criterion of the ComCertSelParams object
 *  pointed to by "params" using the integer value of "minPathLength". In
 *  order to match against this criterion, there are several possibilities.
 *
 *  1) If the criterion's minimum path length is greater than or equal to zero,
 *  a certificate must include a BasicConstraints extension with a pathLen of
 *  at least this value.
 *
 *  2) If the criterion's minimum path length is -2, a certificate must be an
 *  end-entity certificate.
 *
 *  3) If the criterion's minimum path length is -1, no basic constraints check
 *  is done and all certificates are considered to match this criterion.
 *
 *  The semantics of other values of the criterion's minimum path length are
 *  undefined but may be defined in future versions of the API.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose basic constraints
 *      criterion is to be set. Must be non-NULL.
 *  "minPathLength"
 *      Value of PKIX_Int32 used to set the criterion
 *      (or -1 to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetBasicConstraints(
        PKIX_ComCertSelParams *params,
        PKIX_Int32 minPathLength,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetCertificate
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Cert (if any) representing the certificate
 *  criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pCert". In order to match against this
 *  criterion, a certificate must be equal to the criterion's certificate. If
 *  this criterion is specified, it is usually not necessary to specify any
 *  other criteria, since this criterion requires an exact certificate match.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pCert", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose certificate criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pCert"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetCertificate(
        PKIX_ComCertSelParams *params,
        PKIX_PL_Cert **pCert,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetCertificate
 * DESCRIPTION:
 *
 *  Sets the certificate criterion of the ComCertSelParams object pointed to by
 * "params" using a Cert pointed to by "cert". In order to match against this
 *  criterion, a certificate must be equal to the criterion's certificate.
 *  If this criterion is specified, it is usually not necessary to specify
 *  any other criteria, since this criterion requires an exact certificate
 *  match.
 *
 *  If "cert" is NULL, all certificates are considered to match this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose certificate criterion is to be
 *      set. Must be non-NULL.
 *  "cert"
 *      Address of Cert used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetCertificate(
        PKIX_ComCertSelParams *params,
        PKIX_PL_Cert *cert,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetCertificateValid
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Date (if any) representing the certificate
 *  validity criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pDate". In order to match against this
 *  criterion, a certificate's validity period must include the criterion's
 *  Date.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pDate", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose certificate validity criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pDate"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetCertificateValid(
        PKIX_ComCertSelParams *params,
        PKIX_PL_Date **pDate,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetCertificateValid
 * DESCRIPTION:
 *
 *  Sets the certificate validity criterion of the ComCertSelParams object
 *  pointed to by "params" using a Date pointed to by "date". In order to
 *  match against this criterion, a certificate's validity period must include
 *  the criterion's Date.
 *
 *  If "date" is NULL, all certificates are considered to match this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose certificate validity criterion
 *      is to be set. Must be non-NULL.
 *  "date"
 *      Address of Date used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetCertificateValid(
        PKIX_ComCertSelParams *params,
        PKIX_PL_Date *date,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSerialNumber
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the BigInt (if any) representing the serial number
 *  criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pSerialNumber". In order to match against this
 *  criterion, a certificate must have a serial number equal to the
 *  criterion's serial number.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pSerialNumber", in which case all certificates are considered to match
 *  this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose serial number criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pSerialNumber"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSerialNumber(
        PKIX_ComCertSelParams *params,
        PKIX_PL_BigInt **pSerialNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSerialNumber
 * DESCRIPTION:
 *
 *  Sets the serial number criterion of the ComCertSelParams object pointed to
 *  by "params" using a BigInt pointed to by "serialNumber". In order to match
 *  against this criterion, a certificate must have a serial number equal to
 *  the criterion's serial number.
 *
 *  If "serialNumber" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose serial number criterion is to
 *      be set. Must be non-NULL.
 *  "serialNumber"
 *      Address of BigInt used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSerialNumber(
        PKIX_ComCertSelParams *params,
        PKIX_PL_BigInt *serialNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetVersion
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the version criterion that is
 *  set in the ComCertSelParams object pointed to by "params" and stores it at
 *  "pVersion". In order to match against this criterion, a certificate's
 *  version must be equal to the criterion's version.
 *
 *  The version number will either be 0, 1, or 2 (corresponding to
 *  v1, v2, or v3, respectively).
 *
 *  If "params" does not have this criterion set, this function stores
 *  0xFFFFFFFF at "pVersion", in which case all certificates are considered
 *  to match this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose version criterion (if any) is
 *      to be stored. Must be non-NULL.
 *  "pVersion"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetVersion(
        PKIX_ComCertSelParams *params,
        PKIX_UInt32 *pVersion,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetVersion
 * DESCRIPTION:
 *
 *  Sets the version criterion of the ComCertSelParams object pointed to by
 *  "params" using the integer value of "version". In order to match against
 *  this criterion, a certificate's version must be equal to the criterion's
 *  version. If the criterion's version is -1, no version check is done and
 *  all certificates are considered to match this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose version criterion is to be
 *      set. Must be non-NULL.
 *  "version"
 *      Value of PKIX_Int32 used to set the criterion
 *      (or -1 to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetVersion(
        PKIX_ComCertSelParams *params,
        PKIX_Int32 version,
        void *plContext);


/*
 * FUNCTION: PKIX_ComCertSelParams_GetKeyUsage
 * DESCRIPTION:
 *
 *  Retrieves a PKIX_UInt32 (if any) representing the key usage criterion that
 *  is set in the ComCertSelParams object pointed to by "params" and stores it
 *  at "pKeyUsage". In order to match against this criterion, a certificate
 *  must allow the criterion's key usage values. Note that a certificate that
 *  has no KeyUsage extension implicity allows all key usages. Note also that
 *  this functions supports a maximum of 32 key usage bits.
 *
 *  If "params" does not have this criterion set, this function stores zero at
 *  "pKeyUsage", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose key usage criterion (if any)
 *      is to be stored. Must be non-NULL.
 *  "pKeyUsage"
 *      Address where PKIX_UInt32 will be stored. Must not be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetKeyUsage(
        PKIX_ComCertSelParams *params,
        PKIX_UInt32 *pKeyUsage,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetKeyUsage
 * DESCRIPTION:
 *
 *  Sets the key usage criterion of the ComCertSelParams object pointed to by
 *  "params" using the integer value of "keyUsage". In order to match against
 *  this criterion, a certificate must allow the criterion's key usage values.
 *  Note that a certificate that has no KeyUsage extension implicity allows
 *  all key usages. Note also that this functions supports a maximum of 32 key
 *  usage bits.
 *
 *  If the criterion's key usage value is zero, no key usage check is done and
 *  all certificates are considered to match this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose key usage criterion is to be
 *      set. Must be non-NULL.
 *  "keyUsage"
 *      Value of PKIX_Int32 used to set the criterion
 *      (or zero to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetKeyUsage(
        PKIX_ComCertSelParams *params,
        PKIX_UInt32 keyUsage,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetExtendedKeyUsage
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of OIDs (if any) representing the extended
 *  key usage criterion that is set in the ComCertSelParams object pointed to
 *  by "params" and stores it at "pExtKeyUsage". In order to match against this
 *  criterion, a certificate's ExtendedKeyUsage extension must allow the
 *  criterion's extended key usages. Note that a certificate that has no
 *  ExtendedKeyUsage extension implicity allows all key purposes.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pExtKeyUsage", in which case all certificates are considered to match
 *  this criterion.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose extended key usage criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pExtKeyUsage"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetExtendedKeyUsage(
        PKIX_ComCertSelParams *params,
        PKIX_List **pExtKeyUsage, /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetExtendedKeyUsage
 * DESCRIPTION:
 *
 *  Sets the extended key usage criterion of the ComCertSelParams object
 *  pointed to by "params" using a List of OIDs pointed to by "extKeyUsage".
 *  In order to match against this criterion, a certificate's ExtendedKeyUsage
 *  extension must allow the criterion's extended key usages. Note that a
 *  certificate that has no ExtendedKeyUsage extension implicitly allows all
 *  key purposes.
 *
 *  If "extKeyUsage" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose extended key usage criterion
 *      is to be set. Must be non-NULL.
 *  "extKeyUsage"
 *      Address of List of OIDs used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetExtendedKeyUsage(
        PKIX_ComCertSelParams *params,
        PKIX_List *extKeyUsage,  /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetPolicy
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of OIDs (if any) representing the policy
 *  criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pPolicy". In order to match against this
 *  criterion, a certificate's CertificatePolicies extension must include at
 *  least one of the criterion's policies. If "params" has this criterion set,
 *  but the List of OIDs is empty, then a certificate's CertificatePolicies
 *  extension must include at least some policy.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pPolicy", in which case all certificates are considered to match this
 *  criterion.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose policy criterion (if any) is
 *      to be stored. Must be non-NULL.
 *  "pPolicy"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetPolicy(
        PKIX_ComCertSelParams *params,
        PKIX_List **pPolicy,  /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetPolicy
 * DESCRIPTION:
 *
 *  Sets the policy criterion of the ComCertSelParams object pointed to by
 *  "params" using a List of OIDs pointed to by "policy". In order to match
 *  against this criterion, a certificate's CertificatePolicies extension must
 *  include at least one of the criterion's policies. If "params" has this
 *  criterion set, but the List of OIDs is empty, then a certificate's
 *  CertificatePolicies extension must include at least some policy.
 *
 *  If "policy" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose policy criterion is to be set.
 *      Must be non-NULL.
 *  "policy"
 *      Address of List of OIDs used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetPolicy(
        PKIX_ComCertSelParams *params,
        PKIX_List *policy,    /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetIssuer
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the X500Name (if any) representing the issuer
 *  criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pIssuer". In order to match against this
 *  criterion, a certificate's IssuerName must match the criterion's issuer
 *  name.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pIssuer", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose issuer criterion (if any) is
 *      to be stored. Must be non-NULL.
 *  "pIssuer"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetIssuer(
        PKIX_ComCertSelParams *params,
        PKIX_PL_X500Name **pIssuer,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetIssuer
 * DESCRIPTION:
 *
 *  Sets the issuer criterion of the ComCertSelParams object pointed to by
 *  "params" using an X500Name pointed to by "issuer". In order to match
 *  against this criterion, a certificate's IssuerName must match the
 *  criterion's issuer name.
 *
 *  If "issuer" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose issuer criterion is to be set.
 *      Must be non-NULL.
 *  "issuer"
 *      Address of X500Name used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetIssuer(
        PKIX_ComCertSelParams *params,
        PKIX_PL_X500Name *issuer,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSubject
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the X500Name (if any) representing the subject
 *  criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pSubject". In order to match against this
 *  criterion, a certificate's SubjectName must match the criterion's subject
 *  name.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pSubject", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject criterion (if any) is
 *      to be stored. Must be non-NULL.
 *  "pSubject"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSubject(
        PKIX_ComCertSelParams *params,
        PKIX_PL_X500Name **pSubject,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSubject
 * DESCRIPTION:
 *
 *  Sets the subject criterion of the ComCertSelParams object pointed to by
 *  "params" using an X500Name pointed to by "subject". In order to match
 *  against this criterion, a certificate's SubjectName must match the
 *  criterion's subject name.
 *
 *  If "subject" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject criterion is to be
 *      set. Must be non-NULL.
 *  "subject"
 *      Address of X500Name used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSubject(
        PKIX_ComCertSelParams *params,
        PKIX_PL_X500Name *subject,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetSubjectAsByteArray
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ByteArray (if any) representing the subject
 *  criterion that is set in the ComCertSelParams object pointed to by
 *  "params" and stores it at "pSubject". In order to match against this
 *  criterion, a certificate's SubjectName must match the criterion's subject
 *  name.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pSubject", in which case all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject criterion (if any) is
 *      to be stored. Must be non-NULL.
 *  "pSubject"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetSubjectAsByteArray(
        PKIX_ComCertSelParams *params,
        PKIX_PL_ByteArray **pSubject,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetSubjectAsByteArray
 * DESCRIPTION:
 *
 *  Sets the subject criterion of the ComCertSelParams object pointed to by
 *  "params" using a ByteArray pointed to by "subject". In order to match
 *  against this criterion, a certificate's SubjectName must match the
 *  criterion's subject name.
 *
 *  If "subject" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose subject criterion is to be
 *      set. Must be non-NULL.
 *  "subject"
 *      Address of ByteArray used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetSubjectAsByteArray(
        PKIX_ComCertSelParams *params,
        PKIX_PL_ByteArray *subject,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetNameConstraints
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the X500Name (if any) representing the name
 *  constraints criterion that is set in the ComCertSelParams object pointed
 *  to by "params" and stores it at "pConstraints". In order to match against
 *  this criterion, a certificate's subject and subject alternative names must
 *  be allowed by the criterion's name constraints.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pConstraints", in which case all certificates are considered to match
 *  this criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose name constraints criterion
 *      (if any) is to be stored. Must be non-NULL.
 *  "pConstraints"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetNameConstraints(
        PKIX_ComCertSelParams *params,
        PKIX_PL_CertNameConstraints **pConstraints,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetNameConstraints
 * DESCRIPTION:
 *
 *  Sets the name constraints criterion of the ComCertSelParams object pointed
 *  to by "params" using the CertNameConstraints pointed to by "constraints".
 *  In order to match against this criterion, a certificate's subject and
 *  subject alternative names must be allowed by the criterion's name
 *  constraints.
 *
 *  If "constraints" is NULL, all certificates are considered to match this
 *  criterion.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose name constraints criterion is
 *      to be set. Must be non-NULL.
 *  "constraints"
 *      Address of CertNameConstraints used to set the criterion
 *      (or NULL to disable the criterion).
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetNameConstraints(
        PKIX_ComCertSelParams *params,
        PKIX_PL_CertNameConstraints *constraints,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetMatchAllSubjAltNames
 * DESCRIPTION:
 *
 *  Checks whether the ComCertSelParams object pointed to by "params" indicate
 *  that all subject alternative names are to be matched and stores the Boolean
 *  result at "pMatch". This Boolean value determines the behavior of the
 *  subject alternative names criterion.
 *
 *  In order to match against the subject alternative names criterion, if the
 *  Boolean value at "pMatch" is PKIX_TRUE, a certificate must contain all of
 *  the criterion's subject alternative names. If the Boolean value at
 *  "pMatch" is PKIX_FALSE, a certificate must contain at least one of the
 *  criterion's subject alternative names. The default behavior is as if the
 *  Boolean value at "pMatch" is PKIX_TRUE.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object used to determine whether all
 *      subject alternative names must be matched. Must be non-NULL.
 *  "pMatch"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_GetMatchAllSubjAltNames(
        PKIX_ComCertSelParams *params,
        PKIX_Boolean *pMatch,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetMatchAllSubjAltNames
 * DESCRIPTION:
 *
 *  Sets the match flag of the ComCertSelParams object pointed to by "params"
 *  using the Boolean value of "match". This Boolean value determines the
 *  behavior of the subject alternative names criterion.
 *
 *  In order to match against the subject alternative names criterion, if the
 *  "match" is PKIX_TRUE, a certificate must contain all of the criterion's
 *  subject alternative names. If the "match" is PKIX_FALSE, a certificate
 *  must contain at least one of the criterion's subject alternative names.
 *  The default behavior is as if "match" is PKIX_TRUE.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose match flag is to be set.
 *      Must be non-NULL.
 *  "match"
 *      Boolean value used to set the match flag.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetMatchAllSubjAltNames(
        PKIX_ComCertSelParams *params,
        PKIX_Boolean match,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_GetLeafCertFlag
 * DESCRIPTION:
 *
 * Return "leafCert" flag of the ComCertSelParams structure. If set to true,
 * the flag indicates that a selector should filter out all cert that are not
 * qualified to be a leaf cert according to the specified key/ekey usages.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object used to determine whether all
 *      subject alternative names must be matched. Must be non-NULL.
 *  "pLeafFlag"
 *      Address of returned value.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error*
PKIX_ComCertSelParams_GetLeafCertFlag(
        PKIX_ComCertSelParams *params,
        PKIX_Boolean *pLeafFlag,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCertSelParams_SetLeafCertFlag
 * DESCRIPTION:
 *
 * Sets a flag that if its value is true, indicates that the selector
 * should only pick certs that qualifies to be leaf for this cert path
 * validation.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCertSelParams object whose match flag is to be set.
 *      Must be non-NULL.
 *  "leafFlag"
 *      Boolean value used to set the leaf flag.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCertSelParams_SetLeafCertFlag(
        PKIX_ComCertSelParams *params,
        PKIX_Boolean leafFlag,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CERTSEL_H */
