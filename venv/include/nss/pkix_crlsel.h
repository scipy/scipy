/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the PKIX_CRLSelector and the
 * PKIX_ComCRLSelParams types.
 *
 */


#ifndef _PKIX_CRLSEL_H
#define _PKIX_CRLSEL_H

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

/* PKIX_CRLSelector
 *
 * PKIX_CRLSelectors provide a standard way for the caller to select CRLs
 * based on particular criteria. A CRLSelector is typically used by libpkix
 * to retrieve CRLs from a CertStore during certificate chain validation or
 * building. (see pkix_certstore.h) For example, the caller may wish to only
 * select those CRLs that have a particular issuer or a particular value for a
 * private CRL extension. The MatchCallback allows the caller to specify the
 * custom matching logic to be used by a CRLSelector.

 * By default, the MatchCallback is set to point to the default implementation
 * provided by libpkix, which understands how to process the most common
 * parameters. If the default implementation is used, the caller should set
 * these common parameters using PKIX_CRLSelector_SetCommonCRLSelectorParams.
 * Any common parameter that is not set is assumed to be disabled, which means
 * the default MatchCallback implementation will select all CRLs without
 * regard to that particular disabled parameter. For example, if the
 * MaxCRLNumber parameter is not set, MatchCallback will not filter out any
 * CRL based on its CRL number. As such, if no parameters are set, all are
 * disabled and any CRL will match. If a parameter is disabled, its associated
 * PKIX_ComCRLSelParams_Get* function returns a default value of NULL.
 *
 * If a custom implementation is desired, the default implementation can be
 * overridden by calling PKIX_CRLSelector_SetMatchCallback. In this case, the
 * CRLSelector can be initialized with a crlSelectorContext, which is where
 * the caller can specify the desired parameters the caller wishes to match
 * against. Note that this crlSelectorContext must be a PKIX_PL_Object,
 * allowing it to be reference-counted and allowing it to provide the standard
 * PKIX_PL_Object functions (Equals, Hashcode, ToString, Compare, Duplicate).
 *
 */

/*
 * FUNCTION: PKIX_CRLSelector_MatchCallback
 * DESCRIPTION:
 *
 *  This callback function determines whether the specified CRL pointed to by
 *  "crl" matches the criteria of the CRLSelector pointed to by "selector".
 *  If the CRL matches the CRLSelector's criteria, PKIX_TRUE is stored at
 *  "pMatch". Otherwise PKIX_FALSE is stored at "pMatch".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CRLSelector whose MatchCallback logic and parameters are
 *      to be used. Must be non-NULL.
 *  "crl"
 *      Address of CRL that is to be matched using "selector". Must be non-NULL.
 *  "pMatch"
 *      Address at which Boolean result is stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CRLSelector_MatchCallback)(
        PKIX_CRLSelector *selector,
        PKIX_PL_CRL *crl,
        PKIX_Boolean *pMatch,
        void *plContext);

/*
 * FUNCTION: PKIX_CRLSelector_Create
 * DESCRIPTION:
 *
 *  Creates a new CRLSelector using the Object pointed to by
 *  "crlSelectorContext" (if any) and stores it at "pSelector". As noted
 *  above, by default, the MatchCallback is set to point to the default
 *  implementation provided by libpkix, which understands how to process
 *  ComCRLSelParams. This is overridden if the MatchCallback pointed to by
 *  "callback" is not NULL, in which case the parameters are specified using
 *  the Object pointed to by "crlSelectorContext".
 *
 * PARAMETERS:
 *  "issue"
 *      crl issuer.
 *  "crlDpList"
 *      distribution points list
 *  "callback"
 *      The MatchCallback function to be used.
 *  "pSelector"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CRLSelector_Create(
        PKIX_PL_Cert *issuer,
        PKIX_List *crlDpList,
        PKIX_PL_Date *date,
        PKIX_CRLSelector **pSelector,
        void *plContext);

/*
 * FUNCTION: PKIX_CRLSelector_GetMatchCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "selector's" Match callback function and puts it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "selector"
 *      The CRLSelector whose Match callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where Match callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CRLSelector_GetMatchCallback(
        PKIX_CRLSelector *selector,
        PKIX_CRLSelector_MatchCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CRLSelector_GetCRLSelectorContext
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a PKIX_PL_Object representing the context (if any)
 *  of the CRLSelector pointed to by "selector" and stores it at
 *  "pCRLSelectorContext".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CRLSelector whose context is to be stored. Must be non-NULL.
 *  "pCRLSelectorContext"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CRLSelector_GetCRLSelectorContext(
        PKIX_CRLSelector *selector,
        void **pCRLSelectorContext,
        void *plContext);

/*
 * FUNCTION: PKIX_CRLSelector_GetCommonCRLSelectorParams
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the ComCRLSelParams object that represent the common
 *  parameters of the CRLSelector pointed to by "selector" and stores it at
 *  "pCommonCRLSelectorParams". If there are no common parameters stored with
 *  the CRLSelector, this function stores NULL at "pCommonCRLSelectorParams".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CRLSelector whose ComCRLSelParams are to be stored.
 *      Must be non-NULL.
 *  "pCommonCRLSelectorParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CRLSelector_GetCommonCRLSelectorParams(
        PKIX_CRLSelector *selector,
        PKIX_ComCRLSelParams **pCommonCRLSelectorParams,
        void *plContext);

/*
 * FUNCTION: PKIX_CRLSelector_SetCommonCRLSelectorParams
 * DESCRIPTION:
 *
 *  Sets the common parameters for the CRLSelector pointed to by "selector"
 *  using the ComCRLSelParams pointed to by "commonCRLSelectorParams".
 *
 * PARAMETERS:
 *  "selector"
 *      Address of CRLSelector whose common parameters are to be set.
 *      Must be non-NULL.
 *  "commonCRLSelectorParams"
 *      Address of ComCRLSelParams representing the common parameters.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "selector"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CRLSelector_SetCommonCRLSelectorParams(
        PKIX_CRLSelector *selector,
        PKIX_ComCRLSelParams *commonCRLSelectorParams,
        void *plContext);

/* PKIX_ComCRLSelParams
 *
 * PKIX_ComCRLSelParams are X.509 parameters commonly used with CRLSelectors,
 * especially determining which CRLs to retrieve from a CertStore.
 * PKIX_ComCRLSelParams are typically used with those CRLSelectors that use
 * the default implementation of MatchCallback, which understands how to
 * process ComCRLSelParams.
 */

/*
 * FUNCTION: PKIX_ComCRLSelParams_Create
 * DESCRIPTION:
 *
 *  Creates a new ComCRLSelParams object and stores it at "pParams".
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
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_Create(
        PKIX_ComCRLSelParams **pParams,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_GetIssuerNames
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of X500Names (if any) representing the
 *  issuer names criterion that is set in the ComCRLSelParams pointed to by
 *  "params" and stores it at "pNames". In order to match against this
 *  criterion, a CRL's IssuerName must match at least one of the criterion's
 *  issuer names.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pNames", in which case all CRLs are considered to match.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose issuer names criterion (if any) is to
 *      be stored. Must be non-NULL.
 *  "pNames"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_GetIssuerNames(
        PKIX_ComCRLSelParams *params,
        PKIX_List **pNames,  /* list of PKIX_PL_X500Name */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetIssuerNames
 * DESCRIPTION:
 *
 *  Sets the issuer names criterion of the ComCRLSelParams pointed to by
 *  "params" using a List of X500Names pointed to by "names". In order to match
 *  against this criterion, a CRL's IssuerName must match at least one of the
 *  criterion's issuer names.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParamsParams whose issuer names criterion is to be
 *      set. Must be non-NULL.
 *  "names"
 *      Address of List of X500Names used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_SetIssuerNames(
        PKIX_ComCRLSelParams *params,
        PKIX_List *names,   /* list of PKIX_PL_X500Name */
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_AddIssuerName
 * DESCRIPTION:
 *
 *  Adds to the issuer names criterion of the ComCRLSelParams pointed to by
 *  "params" using the X500Name pointed to by "name". In order to match
 *  against this criterion, a CRL's IssuerName must match at least one of the
 *  criterion's issuer names.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose issuer names criterion is to be added
 *      to. Must be non-NULL.
 *  "name"
 *      Address of X500Name to be added.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_AddIssuerName(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_X500Name *name,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_GetCertificateChecking
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Cert (if any) representing the certificate whose
 *  revocation status is being checked. This is not a criterion. It is simply
 *  optional information that may help a CertStore find relevant CRLs.
 *
 *  If "params" does not have a certificate set, this function stores NULL at
 *  "pCert", in which case there is no optional information to provide.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose certificate being checked (if any) is
 *      to be stored. Must be non-NULL.
 *  "pCert"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_GetCertificateChecking(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_Cert **pCert,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetCertificateChecking
 * DESCRIPTION:
 *
 *  Sets the ComCRLSelParams pointed to by "params" with the certificate
 *  (pointed to by "cert") whose revocation status is being checked. This is
 *  not a criterion. It is simply optional information that may help a
 *  CertStore find relevant CRLs.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose certificate being checked is to be
 *      set. Must be non-NULL.
 *  "cert"
 *      Address of Cert whose revocation status is being checked
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_SetCertificateChecking(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_Cert *cert,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_GetDateAndTime
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Date (if any) representing the dateAndTime
 *  criterion that is set in the ComCRLSelParams pointed to by "params" and
 *  stores it at "pDate". In order to match against this criterion, a CRL's
 *  thisUpdate component must be less than or equal to the criterion's
 *  dateAndTime and the CRL's nextUpdate component must be later than the
 *  criterion's dateAndTime. There is no match if the CRL does not contain a
 *  nextUpdate component.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pDate", in which case all CRLs are considered to match.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose dateAndTime criterion (if any) is to
 *      be stored. Must be non-NULL.
 *  "pDate"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_GetDateAndTime(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_Date **pDate,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetDateAndTime
 * DESCRIPTION:
 *
 *  Sets the dateAndTime criterion of the ComCRLSelParams pointed to by
 *  "params" using a Date pointed to by "date". In order to match against this
 *  criterion, a CRL's thisUpdate component must be less than or equal to the
 *  criterion's dateAndTime and the CRL's nextUpdate component must be later
 *  than the criterion's dateAndTime. There is no match if the CRL does not
 *  contain a nextUpdate component.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParamsParams whose dateAndTime criterion is to be
 *      set. Must be non-NULL.
 *  "date"
 *      Address of Date used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_SetDateAndTime(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_Date *date,
        void *plContext);

/* 
 * FUNCTION: PKIX_ComCRLSelParams_GetNISTPolicyEnabled
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Boolean representing the NIST CRL policy
 *  activation flag that is set in the ComCRLSelParams pointed to by "params"
 *  and stores it at "enabled". If enabled, a CRL must have nextUpdate field.
 *
 *  Default value for this flag is TRUE.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose NIST CRL policy criterion  is to
 *      be stored. Must be non-NULL.
 *  "pEnabled"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_GetNISTPolicyEnabled(
        PKIX_ComCRLSelParams *params,
        PKIX_Boolean *pEnabled,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetNISTPolicyEnabled
 * DESCRIPTION:
 *
 *  Sets the NIST crl policy criterion of the ComCRLSelParams pointed to by
 *  "params" using a "enabled" flag. In order to match against this
 *  criterion, a CRL's nextUpdate must be available and criterion's
 *  dataAndTime must be within thisUpdate and nextUpdate time period.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParamsParams whose NIST CRL policy criterion
 *      is to be set. Must be non-NULL.
 *  "enabled"
 *      Address of Bollean used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_SetNISTPolicyEnabled(
        PKIX_ComCRLSelParams *params,
        PKIX_Boolean enabled,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_GetMaxCRLNumber
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the BigInt (if any) representing the maxCRLNumber
 *  criterion that is set in the ComCRLSelParams pointed to by "params" and
 *  stores it at "pNumber". In order to match against this criterion, a CRL
 *  must have a CRL number extension whose value is less than or equal to the
 *  criterion's value.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pNumber", in which case all CRLs are considered to match.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose maxCRLNumber criterion (if any) is to
 *      be stored. Must be non-NULL.
 *  "pNumber"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_GetMaxCRLNumber(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_BigInt **pNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetMaxCRLNumber
 * DESCRIPTION:
 *
 *  Sets the maxCRLNumber criterion of the ComCRLSelParams pointed to by
 *  "params" using a BigInt pointed to by "number". In order to match against
 *  this criterion, a CRL must have a CRL number extension whose value is less
 *  than or equal to the criterion's value.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParamsParams whose maxCRLNumber criterion is to be
 *      set. Must be non-NULL.
 *  "number"
 *      Address of BigInt used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_SetMaxCRLNumber(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_BigInt *number,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_GetMinCRLNumber
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the BigInt (if any) representing the minCRLNumber
 *  criterion that is set in the ComCRLSelParams pointed to by "params" and
 *  stores it at "pNumber". In order to match against this criterion, a CRL
 *  must have a CRL number extension whose value is greater than or equal to
 *  the criterion's value.
 *
 *  If "params" does not have this criterion set, this function stores NULL at
 *  "pNumber", in which case all CRLs are considered to match.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParams whose minCRLNumber criterion (if any) is to
 *      be stored. Must be non-NULL.
 *  "pNumber"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_GetMinCRLNumber(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_BigInt **pNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetMinCRLNumber
 * DESCRIPTION:
 *
 *  Sets the minCRLNumber criterion of the ComCRLSelParams pointed to by
 *  "params" using a BigInt pointed to by "number". In order to match against
 *  this criterion, a CRL must have a CRL number extension whose value is
 *  greater than or equal to the criterion's value.
 *
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParamsParams whose minCRLNumber criterion is to be
 *      set. Must be non-NULL.
 *  "number"
 *      Address of BigInt used to set the criterion
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ComCRLSelParams_SetMinCRLNumber(
        PKIX_ComCRLSelParams *params,
        PKIX_PL_BigInt *number,
        void *plContext);

/*
 * FUNCTION: PKIX_ComCRLSelParams_SetCrlDp
 * DESCRIPTION:
 *
 * Sets crldp list that can be used to download a crls.
 * 
 * PARAMETERS:
 *  "params"
 *      Address of ComCRLSelParamsParams whose minCRLNumber criterion is to be
 *      set. Must be non-NULL.
 *  "crldpList"
 *      A list of CRLDPs. Can be an emptry list.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "params"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRLSelector Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error*
PKIX_ComCRLSelParams_SetCrlDp(
         PKIX_ComCRLSelParams *params,
         PKIX_List *crldpList,
         void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CRLSEL_H */
