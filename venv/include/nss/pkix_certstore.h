/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the PKIX_CertStore type.
 *
 */

#ifndef _PKIX_CERTSTORE_H
#define _PKIX_CERTSTORE_H

#include "pkixt.h"
#include "certt.h"

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

/* PKIX_CertStore
 *
 * A PKIX_CertStore provides a standard way for the caller to retrieve
 * certificates and CRLs from a particular repository (or "store") of
 * certificates and CRLs, including LDAP directories, flat files, local
 * databases, etc. The CertCallback allows custom certificate retrieval logic
 * to be used while the CRLCallback allows custom CRL retrieval logic to be
 * used. Additionally, a CertStore can be initialized with a certStoreContext,
 * which is where the caller can specify configuration data such as the host
 * name of an LDAP server. Note that this certStoreContext must be an
 * Object (although any object type), allowing it to be reference-counted and
 * allowing it to provide the standard Object functions (Equals, Hashcode,
 * ToString, Compare, Duplicate). Please note that each certStoreContext must
 * provide Equals and Hashcode functions in order for the caching (on Cert and
 * CertChain) to work correctly. When providing those two functions, it is not
 * required that all the components of the object be hashed or checked for 
 * equality, but merely that the functions distinguish between unique
 * instances of the certStoreContext.
 *
 * Once the caller has created the CertStore object, the caller then specifies
 * these CertStore objects in a ProcessingParams object and passes that object
 * to PKIX_ValidateChain or PKIX_BuildChain, which uses the objects to call the
 * user's callback functions as needed during the validation or building
 * process.
 *
 * The order of CertStores stored (as a list) at ProcessingParams determines
 * the order in which certificates are retrieved. Trusted CertStores should
 * precede non-trusted ones on the list of CertStores so their certificates
 * are evaluated ahead of other certificates selected on the basis of the same
 * selector criteria.
 *
 * The CheckTrustCallback function is used when the CertStore object
 * supports trust status, which means a Cert's trust status can be altered
 * dynamically. When a CertStore object is created, if the
 * CheckTrustCallback is initialized to be non-NULL, this CertStore is
 * defaulted as supporting trust. Then whenever a Cert needs to (re)check its
 * trust status, this callback can be invoked. When a Cert is retrieved by
 * a CertStore supports trust, at its GetCertCallback, the CertStore
 * information should be updated in Cert's data structure so the link between
 * the Cert and CertStore exists.
 *
 */

/*
 * FUNCTION: PKIX_CertStore_CertCallback
 * DESCRIPTION:
 *
 *  This callback function retrieves from the CertStore pointed to by "store"
 *  all the certificates that match the CertSelector pointed to by "selector".
 *  It places these certificates in a List and stores a pointer to the List at
 *  "pCerts". If no certificates are found which match the CertSelector's
 *  criteria, this function stores an empty List at "pCerts". In either case, if
 *  the operation is completed, NULL is stored at "pNBIOContext".
 *
 *  A CertStore which uses non-blocking I/O may store platform-dependent
 *  information at "pNBIOContext" and NULL at "pCerts" to indicate that I/O is
 *  pending. A subsequent call to PKIX_CertStore_CertContinue is required to
 *  finish the operation and to obtain the List of Certs.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which Certs are to be retrieved.
 *      Must be non-NULL.
 *  "selector"
 *      Address of CertSelector whose criteria must be satisfied.
 *      Must be non-NULL.
 *  "verifyNode"
 *      Parent log node for tracking of filtered out certs.
 *  "pNBIOContext"
 *      Address at which platform-dependent information is stored if the
 *      operation is suspended for non-blocking I/O. Must be non-NULL.
 *  "pCerts"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertStore_CertCallback)(
        PKIX_CertStore *store,
        PKIX_CertSelector *selector,
        PKIX_VerifyNode *verifyNode,
        void **pNBIOContext,
        PKIX_List **pCerts,  /* list of PKIX_PL_Cert */
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_CertContinue
 * DESCRIPTION:
 *
 *  This function continues the non-blocking operation initiated by an earlier
 *  call to the CertCallback function, for the CertStore pointed to by "store". 
 *  If an earlier call did not terminate with the WOULDBLOCK indication (non-NULL
 *  value returned in "pNBIOContext") calling this function will return a fatal
 *  error. If the operation is completed the certificates found are placed in a
 *  List, a pointer to which is stored at "pCerts". If no certificates are found
 *  which match the CertSelector's criteria, this function stores an empty List
 *  at "pCerts". In either case, if the operation is completed, NULL is stored
 *  at "pNBIOContext".
 *
 *  If non-blocking I/O is still pending this function stores platform-dependent
 *  information at "pNBIOContext" and NULL at "pCerts". A subsequent call to
 *  PKIX_CertStore_CertContinue is required to finish the operation and to
 *  obtain the List of Certs.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which Certs are to be retrieved.
 *      Must be non-NULL.
 *  "selector"
 *      Address of CertSelector whose criteria must be satisfied.
 *      Must be non-NULL.
 *  "verifyNode"
 *      Parent log node for tracking of filtered out certs.
 *  "pNBIOContext"
 *      Address at which platform-dependent information is stored if the
 *      operation is suspended for non-blocking I/O. Must be non-NULL.
 *  "pCerts"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_CertContinue(
        PKIX_CertStore *store,
        PKIX_CertSelector *selector,
        PKIX_VerifyNode *verifyNode,
        void **pNBIOContext,
        PKIX_List **pCerts,  /* list of PKIX_PL_Cert */
        void *plContext);

typedef PKIX_Error *
(*PKIX_CertStore_CertContinueFunction)(
        PKIX_CertStore *store,
        PKIX_CertSelector *selector,
        PKIX_VerifyNode *verifyNode,
        void **pNBIOContext,
        PKIX_List **pCerts,  /* list of PKIX_PL_Cert */
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_CRLCallback
 * DESCRIPTION:
 *
 *  This callback function retrieves from the CertStore pointed to by "store"
 *  all the CRLs that match the CRLSelector pointed to by "selector". It
 *  places these CRLs in a List and stores a pointer to the List at "pCRLs".
 *  If no CRLs are found which match the CRLSelector's criteria, this function
 *  stores an empty List at "pCRLs". In either case, if the operation is
 *  completed, NULL is stored at "pNBIOContext".
 *
 *  A CertStore which uses non-blocking I/O may store platform-dependent
 *  information at "pNBIOContext" and NULL at "pCrls" to indicate that I/O is
 *  pending. A subsequent call to PKIX_CertStore_CRLContinue is required to
 *  finish the operation and to obtain the List of Crls.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which CRLs are to be retrieved.
 *      Must be non-NULL.
 *  "selector"
 *      Address of CRLSelector whose criteria must be satisfied.
 *      Must be non-NULL.
 *  "pCrls"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertStore_CRLCallback)(
        PKIX_CertStore *store,
        PKIX_CRLSelector *selector,
        void **pNBIOContext,
        PKIX_List **pCrls,  /* list of PKIX_PL_CRL */
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_ImportCrlCallback
 * DESCRIPTION:
 *
 * The function imports crl list into a cert store. Stores that
 * have local cache may only have that function defined.
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which CRLs are to be retrieved.
 *      Must be non-NULL.
 *  "issuerName"
 *      Name of the issuer that will be used to track bad der crls.
 *  "crlList"
 *      Address on the importing crl list.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertStore_ImportCrlCallback)(
        PKIX_CertStore *store,
        PKIX_PL_X500Name *issuerName,
        PKIX_List *crlList,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_CheckRevokationByCrlCallback
 * DESCRIPTION:
 *
 * The function checks revocation status of a cert with specified
 * issuer, date. It returns revocation status of a cert and
 * a reason code(if any) if a cert was revoked.
 * 
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which CRLs are to be retrieved.
 *      Must be non-NULL.
 *  "cert"
 *      Certificate which revocation status will be checked.
 *  "issuer"
 *      Issuer certificate of the "crl".
 *  "date"
 *      Date of the revocation check.
 *  "crlDownloadDone"
 *      Indicates, that all needed crl downloads are done by the time of
 *      the revocation check.
 *  "reasonCode"
 *      If cert is revoked, returned reason code for  which a cert was revoked.
 *  "revStatus"
 *      Returned revocation status of the cert. See PKIX_RevocationStatus
 *      for more details
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertStore_CheckRevokationByCrlCallback)(
        PKIX_CertStore *store,
        PKIX_PL_Cert *cert,
        PKIX_PL_Cert *issuer,
        PKIX_PL_Date *date,
        PKIX_Boolean  crlDownloadDone,
        CERTCRLEntryReasonCode *reasonCode,
        PKIX_RevocationStatus *revStatus,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_CrlContinue
 * DESCRIPTION:
 *
 *  This function continues the non-blocking operation initiated by an earlier
 *  call to the CRLCallback function, for the CertStore pointed to by "store". 
 *  If an earlier call did not terminate with the WOULDBLOCK indication (non-NULL
 *  value returned in "pNBIOContext") calling this function will return a fatal
 *  error. If the operation is completed the crls found are placed in a List, a
 *  pointer to which is stored at "pCrls". If no crls are found which match the
 *  CRLSelector's criteria, this function stores an empty List at "pCrls". In
 *  either case, if the operation is completed, NULL is stored at "pNBIOContext".
 *
 *  If non-blocking I/O is still pending this function stores platform-dependent
 *  information at "pNBIOContext" and NULL at "pCrls". A subsequent call to
 *  PKIX_CertStore_CrlContinue is required to finish the operation and to
 *  obtain the List of Crls.
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which Crls are to be retrieved.
 *      Must be non-NULL.
 *  "selector"
 *      Address of CRLSelector whose criteria must be satisfied.
 *      Must be non-NULL.
 *  "pNBIOContext"
 *      Address at which platform-dependent information is stored if the
 *      operation is suspended for non-blocking I/O. Must be non-NULL.
 *  "pCrls"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_CrlContinue(
        PKIX_CertStore *store,
        PKIX_CRLSelector *selector,
        void **pNBIOContext,
        PKIX_List **pCrls,  /* list of PKIX_PL_CRL */
        void *plContext);

typedef PKIX_Error *
(*PKIX_CertStore_CrlContinueFunction)(
        PKIX_CertStore *store,
        PKIX_CRLSelector *selector,
        void **pNBIOContext,
        PKIX_List **pCrls,  /* list of PKIX_PL_CRL */
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_CheckTrustCallback
 * DESCRIPTION:
 *
 *  This callback function rechecks "cert's" trust status from the CertStore
 *  pointed to by "store".
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore from which Certs are to be checked.
 *      Must be non-NULL.
 *  "cert"
 *      Address of Cert whose trust status needs to be rechecked.
 *      Must be non-NULL.
 *  "pTrusted"
 *      Address of PKIX_Boolean where the trust status is returned.
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
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_CertStore_CheckTrustCallback)(
        PKIX_CertStore *store,
        PKIX_PL_Cert *cert,
        PKIX_Boolean *pTrusted,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_Create
 * DESCRIPTION:
 *
 *  Creates a new CertStore and stores it at "pStore". The new CertStore uses
 *  the CertCallback pointed to by "certCallback" and the CRLCallback pointed
 *  to by "crlCallback" as its callback functions and uses the Object pointed
 *  to by "certStoreContext" as its context . Note that this certStoreContext
 *  must be an Object (although any object type), allowing it to be
 *  reference-counted and allowing it to provide the standard Object functions
 *  (Equals, Hashcode, ToString, Compare, Duplicate). Once created, a
 *  CertStore object is immutable, although the underlying repository can
 *  change. For example, a CertStore will often be a front-end for a database
 *  or directory. The contents of that directory can change after the
 *  CertStore object is created, but the CertStore object remains immutable.
 *
 * PARAMETERS:
 *  "certCallback"
 *      The CertCallback function to be used. Must be non-NULL.
 *  "crlCallback"
 *      The CRLCallback function to be used. Must be non-NULL.
 *  "certContinue"
 *      The function to be used to resume a certCallback that returned with a
 *      WOULDBLOCK condition. Must be non-NULL if certStore supports non-blocking
 *      I/O.
 *  "crlContinue"
 *      The function to be used to resume a crlCallback that returned with a
 *      WOULDBLOCK condition. Must be non-NULL if certStore supports non-blocking
 *      I/O.
 *  "trustCallback"
 *      Address of PKIX_CertStore_CheckTrustCallback which is called to
 *      verify the trust status of Certs in this CertStore.
 *  "certStoreContext"
 *      Address of Object representing the CertStore's context (if any).
 *  "cachedFlag"
 *      If TRUE indicates data retrieved from CertStore should be cached.
 *  "localFlag"
 *      Boolean value indicating whether this CertStore is local.
 *  "pStore"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CertStore Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_Create(
        PKIX_CertStore_CertCallback certCallback,
        PKIX_CertStore_CRLCallback crlCallback,
        PKIX_CertStore_CertContinueFunction certContinue,
        PKIX_CertStore_CrlContinueFunction crlContinue,
        PKIX_CertStore_CheckTrustCallback trustCallback,
        PKIX_CertStore_ImportCrlCallback importCrlCallback,
        PKIX_CertStore_CheckRevokationByCrlCallback checkRevByCrlCallback,
        PKIX_PL_Object *certStoreContext,
        PKIX_Boolean cachedFlag,
        PKIX_Boolean localFlag,
        PKIX_CertStore **pStore,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetCertCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "store's" Cert callback function and put it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "store"
 *      The CertStore whose Cert callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where Cert callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetCertCallback(
        PKIX_CertStore *store,
        PKIX_CertStore_CertCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetCRLCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "store's" CRL callback function and put it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "store"
 *      The CertStore whose CRL callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where CRL callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetCRLCallback(
        PKIX_CertStore *store,
        PKIX_CertStore_CRLCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetImportCrlCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "store's" Import CRL callback function and put it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "store"
 *      The CertStore whose CRL callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where CRL callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetImportCrlCallback(
        PKIX_CertStore *store,
        PKIX_CertStore_ImportCrlCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetCheckRevByCrl
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "store's" CRL revocation checker callback function
 *  and put it in "pCallback".
 *
 * PARAMETERS:
 *  "store"
 *      The CertStore whose CRL callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where CRL callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetCrlCheckerFn(
        PKIX_CertStore *store,
        PKIX_CertStore_CheckRevokationByCrlCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetTrustCallback
 * DESCRIPTION:
 *
 *  Retrieves the function pointer to the CheckTrust callback function of the
 *  CertStore pointed to by "store" and stores it at "pCallback".
 *
 * PARAMETERS:
 *  "store"
 *      The CertStore whose CheckTrust callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where CheckTrust callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetTrustCallback(
        PKIX_CertStore *store,
        PKIX_CertStore_CheckTrustCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetCertStoreContext
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Object representing the context (if any)
 *  of the CertStore pointed to by "store" and stores it at
 *  "pCertStoreContext".
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore whose context is to be stored. Must be non-NULL.
 *  "pCertStoreContext"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetCertStoreContext(
        PKIX_CertStore *store,
        PKIX_PL_Object **pCertStoreContext,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetCertStoreCacheFlag
 * DESCRIPTION:
 *
 *  Retrieves the Boolean cache flag of the CertStore pointed to by "store" and
 *  stores it at "pCachedFlag".
 *
 * PARAMETERS:
 *  "store"
 *      Address of CertStore whose cache flag is to be stored. Must be non-NULL.
 *  "pCacheFlag"
 *      Address where the result will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetCertStoreCacheFlag(
        PKIX_CertStore *store,
        PKIX_Boolean *pCacheFlag,
        void *plContext);

/*
 * FUNCTION: PKIX_CertStore_GetLocalFlag
 * DESCRIPTION:
 *
 *  Retrieves the Boolean localFlag for the CertStore pointed to by "store" and
 *  stores it at "pLocalFlag". The localFlag is TRUE if the CertStore can
 *  fulfill a request without performing network I/O.
 *
 * PARAMETERS:
 *  "store"
 *      The CertStore whose Local flag is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where the Boolean LocalFlag will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_CertStore_GetLocalFlag(
        PKIX_CertStore *store,
        PKIX_Boolean *pLocalFlag,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CERTSTORE_H */
