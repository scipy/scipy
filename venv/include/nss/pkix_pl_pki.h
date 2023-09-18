/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines several platform independent functions to
 * manipulate certificates and CRLs in a portable manner.
 *
 */

#ifndef _PKIX_PL_PKI_H
#define _PKIX_PL_PKI_H

#include "pkixt.h"
#include "seccomon.h"
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

/*
 * Cert
 *
 * A Cert represents an X.509 certificate. It can be created using the bytes
 * of a valid ASN.1 DER encoding. Once created, a Cert is immutable. The
 * following functions include accessors (gettors) for the various components
 * of an X.509 certificate. Also included are functions to perform various
 * checks on a certificate, including name constraints, key usage, validity
 * (expiration), and signature verification.
 */

/*
 * FUNCTION: PKIX_PL_Cert_Create
 * DESCRIPTION:
 *
 *  Creates a new certificate using the bytes in the ByteArray pointed to by
 *  "byteArray" and stores it at "pCert". If the bytes are not a valid ASN.1
 *  DER encoding of a certificate, a PKIX_Error pointer is returned. Once
 *  created, a Cert is immutable.
 *
 *  Certificate  ::=  SEQUENCE  {
 *      tbsCertificate          TBSCertificate,
 *      signatureAlgorithm      AlgorithmIdentifier,
 *      signatureValue          BIT STRING  }
 *
 *  AlgorithmIdentifier  ::=  SEQUENCE  {
 *      algorithm               OBJECT IDENTIFIER,
 *      parameters              ANY DEFINED BY algorithm OPTIONAL  }
 *
 *  TBSCertificate  ::=  SEQUENCE  {
 *      version         [0]  EXPLICIT Version DEFAULT v1,
 *      serialNumber    CertificateSerialNumber,
 *      signature       AlgorithmIdentifier,
 *      issuer          Name,
 *      validity        Validity,
 *      subject         Name,
 *      subjectPublicKeyInfo SubjectPublicKeyInfo,
 *      issuerUniqueID  [1]  IMPLICIT UniqueIdentifier OPTIONAL,
 *                          -- If present, version MUST be v2 or v3
 *      subjectUniqueID [2]  IMPLICIT UniqueIdentifier OPTIONAL,
 *                              -- If present, version MUST be v2 or v3
 *      extensions      [3]  EXPLICIT Extensions OPTIONAL
 *                              -- If present, version MUST be v3
 *      }
 *
 *  Version  ::=  INTEGER  {  v1(0), v2(1), v3(2)  }
 *
 *  CertificateSerialNumber  ::=  INTEGER
 *
 *  Validity ::= SEQUENCE {
 *      notBefore       Time,
 *      notAfter        Time }
 *
 *  Time ::= CHOICE {
 *      utcTime         UTCTime,
 *      generalTime     GeneralizedTime }
 *
 *  UniqueIdentifier  ::=  BIT STRING
 *
 *  SubjectPublicKeyInfo  ::=  SEQUENCE  {
 *      algorithm               AlgorithmIdentifier,
 *      subjectPublicKey        BIT STRING  }
 *
 *  Extensions  ::=  SEQUENCE SIZE (1..MAX) OF Extension
 *
 *  Extension  ::=  SEQUENCE  {
 *      extnID          OBJECT IDENTIFIER,
 *      critical        BOOLEAN DEFAULT FALSE,
 *      extnValue       OCTET STRING  }
 *
 * PARAMETERS:
 *  "byteArray"
 *      Address of ByteArray representing the CERT's DER encoding.
 *      Must be non-NULL.
 *  "pCert"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_Create(
        PKIX_PL_ByteArray *byteArray,
        PKIX_PL_Cert **pCert,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_CreateFromCERTCertificate
 * DESCRIPTION:
 *
 * Creates a new certificate using passed in CERTCertificate object.
 *
 * PARAMETERS:
 *  "nssCert"
 *      The object that will be used to create new PKIX_PL_Cert.
 *  "pCert"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_CreateFromCERTCertificate(
        const CERTCertificate *nssCert,
        PKIX_PL_Cert **pCert,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetCERTCertificate
 * DESCRIPTION:
 *
 * Returns underlying CERTCertificate structure. Return CERTCertificate
 * object is duplicated and should be destroyed by caller.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of PKIX_PL_Cert. Must be non-NULL.
 *  "pCert"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetCERTCertificate(
        PKIX_PL_Cert *cert,
        CERTCertificate **pnssCert, 
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetVersion
 * DESCRIPTION:
 *
 *  Retrieves the version of the Cert pointed to by "cert" and stores it at
 *  "pVersion". The version number will either be 0, 1, or 2 (corresponding to
 *  v1, v2, or v3, respectively).
 *
 *  Version  ::=  INTEGER  {  v1(0), v2(1), v3(2)  }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose version is to be stored. Must be non-NULL.
 *  "pVersion"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetVersion(
        PKIX_PL_Cert *cert,
        PKIX_UInt32 *pVersion,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetSerialNumber
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the BigInt that represents the serial number of the
 *  Cert pointed to by "cert" and stores it at "pSerialNumber".
 *
 *  CertificateSerialNumber  ::=  INTEGER
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose serial number is to be stored. Must be non-NULL.
 *  "pSerial"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSerialNumber(
        PKIX_PL_Cert *cert,
        PKIX_PL_BigInt **pSerial,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetIssuer
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the X500Name that represents the issuer DN of the
 *  Cert pointed to by "cert" and stores it at "pIssuer".
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose issuer is to be stored. Must be non-NULL.
 *  "pIssuer"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetIssuer(
        PKIX_PL_Cert *cert,
        PKIX_PL_X500Name **pIssuer,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetSubject
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the X500Name that represents the subject DN of the
 *  Cert pointed to by "cert" and stores it at "pSubject". If the Cert does not
 *  have a subject DN, this function stores NULL at "pSubject".
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subject is to be stored. Must be non-NULL.
 *  "pSubject"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSubject(
        PKIX_PL_Cert *cert,
        PKIX_PL_X500Name **pSubject,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetSubjectPublicKeyAlgId
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the OID that represents the subject public key
 *  algorithm of the Cert pointed to by "cert" and stores it at
 *  "pSubjKeyAlgId".
 *
 *  SubjectPublicKeyInfo  ::=  SEQUENCE  {
 *      algorithm               AlgorithmIdentifier,
 *      subjectPublicKey        BIT STRING  }
 *
 *  AlgorithmIdentifier  ::=  SEQUENCE  {
 *      algorithm               OBJECT IDENTIFIER,
 *      parameters              ANY DEFINED BY algorithm OPTIONAL  }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subject public key algorithm OID is to be stored.
 *      Must be non-NULL.
 *  "pSubjKeyAlgId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSubjectPublicKeyAlgId(
        PKIX_PL_Cert *cert,
        PKIX_PL_OID **pSubjKeyAlgId,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetSubjectPublicKey
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the PublicKey that represents the subject public key
 *  of the Cert pointed to by "cert" and stores it at "pPublicKey".
 *
 *  SubjectPublicKeyInfo  ::=  SEQUENCE  {
 *      algorithm               AlgorithmIdentifier,
 *      subjectPublicKey        BIT STRING  }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subject public key is to be stored.
 *      Must be non-NULL.
 *  "pPublicKey"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSubjectPublicKey(
        PKIX_PL_Cert *cert,
        PKIX_PL_PublicKey **pPublicKey,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_PublicKey_NeedsDSAParameters
 * DESCRIPTION:
 *
 * Determines if the PublicKey pointed to by "pubKey" is a DSA Key with null
 * parameters and stores the result at "pNeedsParams". 
 *
 * PARAMETERS:
 *  "pubKey"
 *      Address of the Public Key of interest. Must be non-NULL.
 *  "pNeedsParams"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a PublicKey Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_PublicKey_NeedsDSAParameters(
        PKIX_PL_PublicKey *pubKey,
        PKIX_Boolean *pNeedsParams,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_PublicKey_MakeInheritedDSAPublicKey
 * DESCRIPTION:
 *
 * This function is used for DSA key parameter inheritance, which allows a
 * first DSA key with omitted parameters (pointed to by "firstKey") to inherit
 * the PQG parameters of a second DSA key that does have parameters. (pointed
 * to by "secondKey"). Once created, a PublicKey is immutable.
 *
 * Specifically, the algorithm used by the function is:
 *
 * If the first PublicKey is not a DSA public key with omitted parameters,
 *      the function stores NULL at "pResultKey". (No Error is returned)
 * Else if the second PublicKey is not a DSA public key with non-NULL,
 *      parameters, the function returns an Error.
 * Else
 *      the function creates a third PublicKey with a "Y" value from the
 *      first PublicKey and the DSA parameters from the second PublicKey,
 *      and stores it at "pResultKey".
 *
 * PARAMETERS:
 *  "firstKey"
 *      Address of a Public Key that needs to inherit DSA parameters.
 *      Must be non-NULL.
 *  "secondKey"
 *      Address of a Public Key that has DSA parameters that will be inherited
 *      by "firstKey". Must be non-NULL.
 *  "pResultKey"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a PublicKey Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_PublicKey_MakeInheritedDSAPublicKey(
        PKIX_PL_PublicKey *firstKey,
        PKIX_PL_PublicKey *secondKey,
        PKIX_PL_PublicKey **pResultKey,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetCriticalExtensionOIDs
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of OIDs (each OID corresponding to a
 *  critical extension of the Cert pointed to by "cert") and stores it at
 *  "pExtensions". If "cert" does not have any critical extensions, this
 *  function stores an empty List at "pExtensions".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose critical extension OIDs are to be stored.
 *      Must be non-NULL.
 *  "pExtensions"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetCriticalExtensionOIDs(
        PKIX_PL_Cert *cert,
        PKIX_List **pExtensions,  /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetAuthorityKeyIdentifier
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a ByteArray representing the authority key
 *  identifier extension of the Cert pointed to by "cert" and stores it at
 *  "pAuthKeyId".
 *
 *  Note that this function only retrieves the keyIdentifier component
 *  (OCTET STRING) of the AuthorityKeyIdentifier extension, when present.
 *
 *  If "cert" does not have an AuthorityKeyIdentifier extension or if the
 *  keyIdentifier component of the AuthorityKeyIdentifier extension is not
 *  present, this function stores NULL at "pAuthKeyId".
 *
 *  AuthorityKeyIdentifier ::= SEQUENCE {
 *      keyIdentifier                   [0] KeyIdentifier           OPTIONAL,
 *      authorityCertIssuer             [1] GeneralNames            OPTIONAL,
 *      authorityCertSerialNumber       [2] CertificateSerialNumber OPTIONAL  }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose authority key identifier is to be stored.
 *      Must be non-NULL.
 *  "pAuthKeyId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetAuthorityKeyIdentifier(
        PKIX_PL_Cert *cert,
        PKIX_PL_ByteArray **pAuthKeyId,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetSubjectKeyIdentifier
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a ByteArray representing the subject key identifier
 *  extension of the Cert pointed to by "cert" and stores it at "pSubjKeyId".
 *  If "cert" does not have a SubjectKeyIdentifier extension, this function
 *  stores NULL at "pSubjKeyId".
 *
 *  SubjectKeyIdentifier ::= KeyIdentifier
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subject key identifier is to be stored.
 *      Must be non-NULL.
 *  "pSubjKeyId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSubjectKeyIdentifier(
        PKIX_PL_Cert *cert,
        PKIX_PL_ByteArray **pSubjKeyId,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetSubjectAltNames
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of GeneralNames (each GeneralName
 *  representing a subject alternative name found in the subject alternative
 *  names extension of the Cert pointed to by "cert") and stores it at
 *  "pSubjectAltNames". If "cert" does not have a SubjectAlternativeNames
 *  extension, this function stores NULL at "pSubjectAltNames".
 *
 *  Note that the List returned by this function is immutable.
 *
 *  SubjectAltName ::= GeneralNames
 *
 *  GeneralNames ::= SEQUENCE SIZE (1..MAX) OF GeneralName
 *
 *  GeneralName ::= CHOICE {
 *      otherName                       [0]     OtherName,
 *      rfc822Name                      [1]     IA5String,
 *      dNSName                         [2]     IA5String,
 *      x400Address                     [3]     ORAddress,
 *      directoryName                   [4]     Name,
 *      ediPartyName                    [5]     EDIPartyName,
 *      uniformResourceIdentifier       [6]     IA5String,
 *      iPAddress                       [7]     OCTET STRING,
 *      registeredID                    [8]     OBJECT IDENTIFIER }
 *
 *  OtherName ::= SEQUENCE {
 *      type-id                         OBJECT IDENTIFIER,
 *      value                           [0] EXPLICIT ANY DEFINED BY type-id }
 *
 *  EDIPartyName ::= SEQUENCE {
 *      nameAssigner                    [0]     DirectoryString OPTIONAL,
 *      partyName                       [1]     DirectoryString }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subjectAltNames are to be stored.
 *      Must be non-NULL.
 *  "pSubjectAltNames"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSubjectAltNames(
        PKIX_PL_Cert *cert,
        PKIX_List **pSubjectAltNames,  /* list of PKIX_PL_GeneralName */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetAllSubjectNames
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of GeneralNames (each GeneralName
 *  representing a subject DN or a subject alternative name found in the
 *  subject alternative names extension of the Cert pointed to by "cert") and
 *  stores it at "pAllSubjectNames".If the Subject DN of "cert" is empty and
 *  it does not have a SubjectAlternativeNames extension, this function stores
 *  NULL at "pAllSubjectNames".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subject DN and subjectAltNames are to be stored.
 *      Must be non-NULL.
 *  "pAllSubjectNames"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetAllSubjectNames(
        PKIX_PL_Cert *cert,
        PKIX_List **pAllSubjectNames,  /* list of PKIX_PL_GeneralName */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetExtendedKeyUsage
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a List of OIDs (each OID corresponding to an
 *  extended key usage of the Cert pointed to by "cert") and stores it at
 *  "pKeyUsage". If "cert" does not have an extended key usage extension, this
 *  function stores a NULL at "pKeyUsage".
 *
 *  Note that the List returned by this function is immutable.
 *
 *  ExtKeyUsageSyntax ::= SEQUENCE SIZE (1..MAX) OF KeyPurposeId
 *
 *  KeyPurposeId ::= OBJECT IDENTIFIER
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose extended key usage OIDs are to be stored.
 *      Must be non-NULL.
 *  "pKeyUsage"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetExtendedKeyUsage(
        PKIX_PL_Cert *cert,
        PKIX_List **pKeyUsage,  /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetNameConstraints
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a CertNameConstraints object representing the name
 *  constraints extension of the Cert pointed to by "cert" and stores it at
 *  "pNameConstraints".
 *
 *  If "cert" does not have a name constraints extension, this function stores
 *  NULL at "pNameConstraints".
 *
 *  NameConstraints ::= SEQUENCE {
 *      permittedSubtrees       [0]     GeneralSubtrees OPTIONAL,
 *      excludedSubtrees        [1]     GeneralSubtrees OPTIONAL }
 *
 *  GeneralSubtrees ::= SEQUENCE SIZE (1..MAX) OF GeneralSubtree
 *
 *  GeneralSubtree ::= SEQUENCE {
 *      base                    GeneralName,
 *      minimum         [0]     BaseDistance DEFAULT 0,
 *      maximum         [1]     BaseDistance OPTIONAL }
 *
 *  BaseDistance ::= INTEGER (0..MAX)
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose name constraints extension is to be stored.
 *      Must be non-NULL.
 *  "pNameConstraints"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetNameConstraints(
        PKIX_PL_Cert *cert,
        PKIX_PL_CertNameConstraints **pNameConstraints,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetBasicConstraints
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a CertBasicConstraints object representing the basic
 *  constraints extension of the Cert pointed to by "cert" and stores it at
 *  "pBasicConstraints".
 *
 *  If "cert" does not have a basic constraints extension, this function stores
 *  NULL at "pBasicConstraints". Once created, a CertBasicConstraints object
 *  is immutable.
 *
 *  BasicConstraints ::= SEQUENCE {
 *      cA                      BOOLEAN DEFAULT FALSE,
 *      pathLenConstraint       INTEGER (0..MAX) OPTIONAL }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose basic constraints extension is to be stored.
 *      Must be non-NULL.
 *  "pBasicConstraints"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetBasicConstraints(
        PKIX_PL_Cert *cert,
        PKIX_PL_CertBasicConstraints **pBasicConstraints,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_BasicConstraints_GetCAFlag
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a Boolean value representing the cA Flag component
 *  of the CertBasicConstraints object pointed to by "basicConstraints" and
 *  stores it at "pResult".
 *
 *  BasicConstraints ::= SEQUENCE {
 *      cA                      BOOLEAN DEFAULT FALSE,
 *      pathLenConstraint       INTEGER (0..MAX) OPTIONAL }
 *
 * PARAMETERS:
 *  "basicConstraints"
 *      Address of CertBasicConstraints whose cA Flag is to be stored.
 *      Must be non-NULL.
 *  "pResult"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_BasicConstraints_GetCAFlag(
        PKIX_PL_CertBasicConstraints *basicConstraints,
        PKIX_Boolean *pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_BasicConstraints_GetPathLenConstraint
 * DESCRIPTION:
 *
 *  Retrieves a pointer to an integer value representing the pathLenConstraint
 *  component of the CertBasicConstraints object pointed to by
 *  "basicConstraints" and stores it at "pPathLenConstraint". If the
 *  pathLenConstraint component is not present, this function stores -1 at
 *  "pPathLenConstraint".
 *
 * PARAMETERS:
 *  "basicConstraints"
 *      Address of CertBasicConstraints whose pathLen is to be stored.
 *      Must be non-NULL.
 *  "pPathLenConstraint"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_BasicConstraints_GetPathLenConstraint(
        PKIX_PL_CertBasicConstraints *basicConstraints,
        PKIX_Int32 *pPathLenConstraint,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetPolicyInformation
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a List of CertPolicyInfos found in the certificate
 *  policies extension of the Cert pointed to by "cert" and stores it at
 *  "pPolicyInfo". If "cert" does not have a certificate policies extension,
 *  this function stores NULL at "pPolicyInfo". Once created, a CertPolicyInfo
 *  object is immutable.
 *
 *  Note that the List returned by this function is immutable.
 *
 *  certificatePolicies ::= SEQUENCE SIZE (1..MAX) OF PolicyInformation
 *
 *  PolicyInformation ::= SEQUENCE {
 *      policyIdentifier   CertPolicyId,
 *      policyQualifiers   SEQUENCE SIZE (1..MAX) OF
 *                              PolicyQualifierInfo OPTIONAL }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose CertPolicyInfos are to be stored.
 *      Must be non-NULL.
 *  "pPolicyInfo"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetPolicyInformation(
        PKIX_PL_Cert *cert,
        PKIX_List **pPolicyInfo, /* list of PKIX_PL_CertPolicyInfo */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CertPolicyInfo_GetPolicyId
 * DESCRIPTION:
 *
 *  Retrieves a pointer to an OID representing the policyIdentifier of the
 *  CertPolicyInfo pointed to by "policyInfo" and stores it at "pCertPolicyId".
 *
 *  certificatePolicies ::= SEQUENCE SIZE (1..MAX) OF PolicyInformation
 *
 *  PolicyInformation ::= SEQUENCE {
 *      policyIdentifier   CertPolicyId,
 *      policyQualifiers   SEQUENCE SIZE (1..MAX) OF
 *                              PolicyQualifierInfo OPTIONAL }
 *
 *  CertPolicyId ::= OBJECT IDENTIFIER
 *
 * PARAMETERS:
 *  "policyInfo"
 *      Address of CertPolicyInfo whose policy identifier is to be stored.
 *      Must be non-NULL.
 *  "pCertPolicyId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CertPolicyInfo_GetPolicyId(
        PKIX_PL_CertPolicyInfo *policyInfo,
        PKIX_PL_OID **pCertPolicyId,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CertPolicyInfo_GetPolQualifiers
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a List of the CertPolicyQualifiers representing
 *  the policyQualifiers of the CertPolicyInfo pointed to by "policyInfo" and
 *  stores it at "pPolicyQualifiers". If "policyInfo" does not have any
 *  policyQualifiers, this function stores NULL at "pPolicyQualifiers". Once
 *  created, a CertPolicyQualifier is immutable.
 *
 *  Note that the List returned by this function is immutable.
 *
 *  certificatePolicies ::= SEQUENCE SIZE (1..MAX) OF PolicyInformation
 *
 *  PolicyInformation ::= SEQUENCE {
 *      policyIdentifier   CertPolicyId,
 *      policyQualifiers   SEQUENCE SIZE (1..MAX) OF
 *                              PolicyQualifierInfo OPTIONAL }
 *
 *  PolicyQualifierInfo ::= SEQUENCE {
 *      policyQualifierId  PolicyQualifierId,
 *      qualifier       ANY DEFINED BY policyQualifierId }
 *
 * PARAMETERS:
 *  "policyInfo"
 *      Address of CertPolicyInfo whose policy qualifiers List is to be stored.
 *      Must be non-NULL.
 *  "pPolicyQualifiers"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CertPolicyInfo_GetPolQualifiers(
        PKIX_PL_CertPolicyInfo *policyInfo,
        PKIX_List **pPolicyQualifiers, /* list of PKIX_PL_CertPolicyQualifier */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_PolicyQualifier_GetPolicyQualifierId
 * DESCRIPTION:
 *
 *  Retrieves a pointer to an OID representing the policyQualifierId of the
 *  CertPolicyQualifier pointed to by "policyQualifier" and stores it at
 *  "pPolicyQualifierId".
 *
 *  PolicyQualifierInfo ::= SEQUENCE {
 *      policyQualifierId       PolicyQualifierId,
 *      qualifier               ANY DEFINED BY policyQualifierId }
 *
 *  PolicyQualifierId ::=
 *      OBJECT IDENTIFIER ( id-qt-cps | id-qt-unotice )
 *
 * PARAMETERS:
 *  "policyQualifier"
 *      Address of CertPolQualifier whose policyQualifierId is to be stored.
 *      Must be non-NULL.
 *  "pPolicyQualifierId"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_PolicyQualifier_GetPolicyQualifierId(
        PKIX_PL_CertPolicyQualifier *policyQualifier,
        PKIX_PL_OID **pPolicyQualifierId,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_PolicyQualifier_GetQualifier
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a ByteArray representing the qualifier of the
 *  CertPolicyQualifier pointed to by "policyQualifier" and stores it at
 *  "pQualifier".
 *
 *  PolicyQualifierInfo ::= SEQUENCE {
 *      policyQualifierId       PolicyQualifierId,
 *      qualifier               ANY DEFINED BY policyQualifierId }
 *
 * PARAMETERS:
 *  "policyQualifier"
 *      Address of CertPolicyQualifier whose qualifier is to be stored.
 *      Must be non-NULL.
 *  "pQualifier"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_PolicyQualifier_GetQualifier(
        PKIX_PL_CertPolicyQualifier *policyQualifier,
        PKIX_PL_ByteArray **pQualifier,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetPolicyMappings
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a List of CertPolicyMaps found in the policy
 *  mappings extension of the Cert pointed to by "cert" and stores it at
 *  "pPolicyMappings". If "cert" does not have a policy mappings extension,
 *  this function stores NULL at "pPolicyMappings". Once created, a
 *  CertPolicyMap is immutable.
 *
 *  Note that the List returned by this function is immutable.
 *
 *  PolicyMappings ::= SEQUENCE SIZE (1..MAX) OF SEQUENCE {
 *      issuerDomainPolicy      CertPolicyId,
 *      subjectDomainPolicy     CertPolicyId }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose CertPolicyMaps are to be stored.
 *      Must be non-NULL.
 *  "pPolicyMappings"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetPolicyMappings(
        PKIX_PL_Cert *cert,
        PKIX_List **pPolicyMappings, /* list of PKIX_PL_CertPolicyMap */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CertPolicyMap_GetIssuerDomainPolicy
 * DESCRIPTION:
 *
 *  Retrieves a pointer to an OID representing the issuerDomainPolicy of the
 *  CertPolicyMap pointed to by "policyMapping" and stores it at
 *  "pIssuerDomainPolicy".
 *
 *  PolicyMappings ::= SEQUENCE SIZE (1..MAX) OF SEQUENCE {
 *      issuerDomainPolicy      CertPolicyId,
 *      subjectDomainPolicy     CertPolicyId }
 *
 * PARAMETERS:
 *  "policyMapping"
 *      Address of CertPolicyMap whose issuerDomainPolicy is to be stored.
 *      Must be non-NULL.
 *  "pIssuerDomainPolicy"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CertPolicyMap_GetIssuerDomainPolicy(
        PKIX_PL_CertPolicyMap *policyMapping,
        PKIX_PL_OID **pIssuerDomainPolicy,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CertPolicyMap_GetSubjectDomainPolicy
 * DESCRIPTION:
 *
 *  Retrieves a pointer to an OID representing the subjectDomainPolicy of the
 *  CertPolicyMap pointed to by "policyMapping" and stores it at
 *  "pSubjectDomainPolicy".
 *
 *  PolicyMappings ::= SEQUENCE SIZE (1..MAX) OF SEQUENCE {
 *      issuerDomainPolicy      CertPolicyId,
 *      subjectDomainPolicy     CertPolicyId }
 *
 * PARAMETERS:
 *  "policyMapping"
 *      Address of CertPolicyMap whose subjectDomainPolicy is to be stored.
 *      Must be non-NULL.
 *  "pSubjectDomainPolicy"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CertPolicyMap_GetSubjectDomainPolicy(
        PKIX_PL_CertPolicyMap *policyMapping,
        PKIX_PL_OID **pSubjectDomainPolicy,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetRequireExplicitPolicy
 * DESCRIPTION:
 *
 *  Retrieves the requireExplicitPolicy value of the policy constraints
 *  extension of the Cert pointed to by "cert" and stores it at "pSkipCerts".
 *  If "cert" does not have a policy constraints extension or the
 *  requireExplicitPolicy component is not populated, this function stores -1
 *  at "pSkipCerts".
 *
 *  PolicyConstraints ::= SEQUENCE {
 *      requireExplicitPolicy   [0] SkipCerts OPTIONAL,
 *      inhibitPolicyMapping    [1] SkipCerts OPTIONAL }
 *
 *  SkipCerts ::= INTEGER (0..MAX)
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose requireExplicitPolicy value is to be stored.
 *      Must be non-NULL.
 *  "pSkipCerts"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetRequireExplicitPolicy(
        PKIX_PL_Cert *cert,
        PKIX_Int32 *pSkipCerts,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetPolicyMappingInhibited
 * DESCRIPTION:
 *
 *  Retrieves the inhibitPolicyMapping value of the policy constraints
 *  extension of the Cert pointed to by "cert" and stores it at "pSkipCerts".
 *  If "cert" does not have a policy constraints extension or the
 *  inhibitPolicyMapping component is not populated, this function stores -1
 *  at "pSkipCerts".
 *
 *  PolicyConstraints ::= SEQUENCE {
 *      requireExplicitPolicy   [0] SkipCerts OPTIONAL,
 *      inhibitPolicyMapping    [1] SkipCerts OPTIONAL }
 *
 *  SkipCerts ::= INTEGER (0..MAX)
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose requireExplicitPolicy value is to be stored.
 *      Must be non-NULL.
 *  "pSkipCerts"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetPolicyMappingInhibited(
        PKIX_PL_Cert *cert,
        PKIX_Int32 *pSkipCerts,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetInhibitAnyPolicy
 * DESCRIPTION:
 *
 *  Retrieves the value of the inhibit any-policy extension of the Cert
 *  pointed to by "cert" and stores it at "pSkipCerts". If "cert" does not have
 *  an inhibit any-policy extension, this function stores -1 at "pSkipCerts".
 *
 *  InhibitAnyPolicy ::= SkipCerts
 *
 *  SkipCerts ::= INTEGER (0..MAX)
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose inhibit any-policy extensions value is to be
 *      stored. Must be non-NULL.
 *  "pSkipCerts"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetInhibitAnyPolicy(
        PKIX_PL_Cert *cert,
        PKIX_Int32 *pSkipCerts,
        void *plContext);

/* policy processing functions */

/*
 * FUNCTION: PKIX_PL_Cert_AreCertPoliciesCritical
 * DESCRIPTION:
 *
 *  Checks whether the certificate policies extension of the Cert pointed to
 *  by "cert" is critical and stores the Boolean result at "pCritical". If
 *  "cert" does not have a certificate policies extension, this function
 *  stores NULL at "pCritical".
 *
 *  XXX what distinguishes NULL from PKIX_FALSE?
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose certificate policies extension's criticality is
 *      to be determined. Must be non-NULL.
 *  "pCritical"
 *      Address where PKIX_Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_AreCertPoliciesCritical(
        PKIX_PL_Cert *cert,
        PKIX_Boolean *pCritical,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_CheckNameConstraints
 * DESCRIPTION:
 *
 *  Checks whether the subject distinguished name and subject alternative names
 *  of the Cert pointed to by "cert" satisfy the CertNameConstraints pointed
 *  to by "nameConstraints". If the CertNameConstraints are not satisfied, a
 *  PKIX_Error pointer is returned. If "nameConstraints" is NULL, the function
 *  does nothing.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose subject names are to be checked.
 *      Must be non-NULL.
 *  "nameConstraints"
 *      Address of CertNameConstraints that need to be satisfied.
 *  "treatCommonNameAsDNSName"
 *      PKIX_TRUE if the subject common name should be considered a dNSName
 *      when evaluating name constraints.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_CheckNameConstraints(
        PKIX_PL_Cert *cert,
        PKIX_PL_CertNameConstraints *nameConstraints,
        PKIX_Boolean treatCommonNameAsDNSName,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_MergeNameConstraints
 * DESCRIPTION:
 *
 *  Merges the CertNameConstraints pointed to by "firstNC" and the
 *  CertNameConstraints pointed to by "secondNC" and stores the merged
 *  CertNameConstraints at "pResultNC". If "secondNC" is NULL, the
 *  CertNameConstraints pointed to by "firstNC" is stored at "pResultNC".
 *
 *  Once created, a CertNameConstraints object is immutable.
 *
 * PARAMETERS:
 *  "firstNC"
 *      Address of first CertNameConstraints to be merged. Must be non-NULL.
 *  "secondNC"
 *      Address of second CertNameConstraints to be merged
 *  "pResultNC"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_MergeNameConstraints(
        PKIX_PL_CertNameConstraints *firstNC,
        PKIX_PL_CertNameConstraints *secondNC,
        PKIX_PL_CertNameConstraints **pResultNC,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_VerifyKeyUsage
 * DESCRIPTION:
 *
 *  Verifies that the keyUsage bit(s) specified by "keyUsage" appear in the
 *  keyUsage extension of the Cert pointed to by "cert". The keyUsage bit
 *  values specified in pkixt.h are supported, and can be bitwise or'ed if
 *  multiple bit values are to be verified. If the keyUsages do not all appear
 *  in the keyUsage extension of "cert", a PKIX_Error pointer is returned.
 *
 *  KeyUsage ::= BIT STRING {
 *      digitalSignature        (0),
 *      nonRepudiation          (1),
 *      keyEncipherment         (2),
 *      dataEncipherment        (3),
 *      keyAgreement            (4),
 *      keyCertSign             (5),
 *      cRLSign                 (6),
 *      encipherOnly            (7),
 *      decipherOnly            (8) }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose keyUsage bits are to be verified.
 *      Must be non-NULL.
 *  "keyUsage"
 *      Constant representing keyUsage bit(s) that all must appear in keyUsage
 *      extension of "cert".
 *  "plContext" - Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_VerifyKeyUsage(
        PKIX_PL_Cert *cert,
        PKIX_UInt32 keyUsage,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_VerifyCertAndKeyType
 * DESCRIPTION:
 *
 * Verifies cert and key types against certificate usage that is
 * a part of plContext(pkix_pl_nsscontext) structure. Throws an error
 * if cert or key types does not match.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose keyUsage bits are to be verified.
 *      Must be non-NULL.
 *  "isLeafCert"
 *      What type of a cert has been verified.
 *  "plContext" - Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_VerifyCertAndKeyType(
        PKIX_PL_Cert *cert,
        PKIX_Boolean isChainCert,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_CheckValidity
 * DESCRIPTION:
 *
 *  Checks whether the Cert pointed to by "cert" would be valid at the time
 *  represented by the Date pointed to by "date". If "date" is NULL, then this
 *  function checks whether the Cert would be valid at the current time. If the
 *  Cert would not be valid at the specified Date, a PKIX_Error pointer is
 *  returned.
 *
 *  Validity ::= SEQUENCE {
 *      notBefore       Time,
 *      notAfter        Time }
 *
 *  Time ::= CHOICE {
 *      utcTime         UTCTime,
 *      generalTime     GeneralizedTime }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose validity is to be checked. Must be non-NULL.
 *  "date"
 *      Address of Date at which the Cert is being checked for validity.
 *      If NULL, the current time is used for the Date.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_CheckValidity(
        PKIX_PL_Cert *cert,
        PKIX_PL_Date *date,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetValidityNotAfter
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the Date that represents the notAfter time of the
 *  Certificate pointed to by "cert" and stores it at "pDate".
 *
 *  Validity ::= SEQUENCE {
 *      notBefore       Time,
 *      notAfter        Time }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose validity time is to be retrieved. Must be
 *      non-NULL.
 *  "date"
 *      Address of Date at which the Cert's notAfter time is being retrieved.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetValidityNotAfter(
        PKIX_PL_Cert *cert,
        PKIX_PL_Date **pDate,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_VerifySignature
 * DESCRIPTION:
 *
 *  Verifies the signature on the Cert pointed to by "cert" using the
 *  PublicKey pointed to by "pubKey". If the signature doesn't verify, an
 *  Error pointer is returned.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose signature is to be verified. Must be non-NULL.
 *  "pubKey"
 *      Address of a Public Key used to verify the signature. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_VerifySignature(
        PKIX_PL_Cert *cert,
        PKIX_PL_PublicKey *pubKey,
        void *plContext);

/* A set of flags to indicate how explicitly configured trust anchors should be
 * handled by PKIX_PL_Cert_IsCertTrusted
 */
typedef enum PKIX_PL_TrustAnchorModeEnum {
        /* Indicates trust anchors should be ignored; only the underlying
         * platform's trust settings should be used.
         */
        PKIX_PL_TrustAnchorMode_Ignore,

        /* Indicates that explicitly configured trust anchors may be considered
         * trustworthy, if present.
         * Note: If the underlying platform supports marking a certificate as
         *       explicitly untrustworthy, explicitly configured trust anchors
         *       MAY be ignored/rejected.
         */
        PKIX_PL_TrustAnchorMode_Additive,

        /* Indicates that ONLY trust anchors should be considered as
         * trustworthy.
         * Note: If the underlying platform supports marking a certificate as
         *       explicitly untrustworthy, explicitly configured trust anchors
         *       MAY be ignored/rejected.
         */
        PKIX_PL_TrustAnchorMode_Exclusive
} PKIX_PL_TrustAnchorMode;

/*
 * FUNCTION: PKIX_PL_Cert_IsCertTrusted
 * DESCRIPTION:
 *
 *  Checks the Cert specified by "cert" to determine, in a manner that depends
 *  on the underlying platform, whether it is trusted, and stores the result in
 *  "pTrusted". If a certificate is trusted it means that a chain built to that
 *  certificate, and satisfying all the usage, policy, validity, and other
 *  tests, is a valid chain and the End Entity certificate from which it was
 *  built can be trusted.
 *
 *  If the Certificate is not intrinsically trustworthy, it still might end up a
 *  component in a successful chain.
 *
 *  If the Certificate is intrinsically untrustworthy, this function will return
 *  an error. 
 *
 * PARAMETERS
 *  "cert"
 *      Address of Cert whose trustworthiness is to be determined. Must be
 *      non-NULL.
 *  "trustAnchorMode"
 *      A PKIX_PL_TrustAnchorMode that indicates how explicitly defined user
 *      trust anchors should be handled.
 *  "pTrusted"
 *      Address where the Boolean value will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CERT Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_IsCertTrusted(
        PKIX_PL_Cert *cert,
        PKIX_PL_TrustAnchorMode trustAnchorMode,
        PKIX_Boolean *pTrusted,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_IsLeafCertTrusted
 * DESCRIPTION:
 *
 *  Checks the Leaf Cert specified by "cert" to determine, in a manner that 
 *  depends on the underlying platform, whether it is trusted, and stores the 
 *  result in "pTrusted". If a certificate is trusted it means that this
 *  End Entify certificate has been marked as trusted for the requested usage,
 *  policy, validity, and other tests.
 *
 *  If the Certificate is not intrinsically trustworthy, we can still try to 
 *  build a successful chain.
 *
 *  If the Certificate is intrinsically untrustworthy, this function will return
 *  an error. 
 *
 * PARAMETERS
 *  "cert"
 *      Address of Cert whose trustworthiness is to be determined. Must be
 *      non-NULL.
 *  "pTrusted"
 *      Address where the Boolean value will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CERT Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_IsLeafCertTrusted(
        PKIX_PL_Cert *cert,
        PKIX_Boolean *pTrusted,
        void *plContext);

/* FUNCTION: PKIX_PL_Cert_SetAsTrustAnchor */
PKIX_Error*
PKIX_PL_Cert_SetAsTrustAnchor(PKIX_PL_Cert *cert, 
                              void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetCacheFlag
 * DESCRIPTION:
 *
 *  Retrieves the value of the cache flag in "cert" and return it at address
 *  pointed by "pCacheFlag". The initila cache flag is determined by the
 *  CertStore this "cert" is fetched from. When CertStore is created, user
 *  need to specify if the data should be cached.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose cache flag is fetched. Must be non-NULL.
 *  "pCacheFlag"
 *      Address where PKIX_Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetCacheFlag(
        PKIX_PL_Cert *cert,
        PKIX_Boolean *pCacheFlag,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_SetCacheFlag
 * DESCRIPTION:
 *
 *  Set the value of the cache flag in "cert" base on the boolean value stored
 *  at "cacheFlag". This function is meant to be used by CertStore after a
 *  Cert is created.
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert where "cacheFlag" is stored. Must be non-NULL.
 *  "cacheFlag"
 *      PKIX_Boolean flag for cache flag.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_SetCacheFlag(
        PKIX_PL_Cert *cert,
        PKIX_Boolean cacheFlag,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_GetTrustCertStore
 * DESCRIPTION:
 *
 *  Retrieves the value of the CertStore in "cert" and return it at address
 *  pointed by "pCertStore".
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose CertStore is fetched. Must be non-NULL.
 *  "pTrustCertStore"
 *      Address where CertStore will be stored and returned. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetTrustCertStore(
        PKIX_PL_Cert *cert,
        PKIX_CertStore **pTrustCertStore,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Cert_SetTrustCertStore
 * DESCRIPTION:
 *
 *  Set the value of the CertStore "certStore" in "cert".
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert where "certStore" will be stored. Must be non-NULL.
 *  "trustCertStore"
 *      Address where the CertStore is. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_SetTrustCertStore(
        PKIX_PL_Cert *cert,
        PKIX_CertStore *trustCertStore,
        void *plContext);


/*
 * FUNCTION: PKIX_PL_Cert_GetAuthorityInfoAccess
 * DESCRIPTION:
 *
 *  Retrieves the value(s) of the Authority Information Access in "cert" and
 *  returns it in a list at address pointed by "pAuthorityInfoAccess".
 *
 *  SubjectInfoAccess ::=
 *    SEQUENCE SIZE (1..MAX) of AccessDescription
 *    AccessDescription ::= SEQUENCE {
 *        accessMethod     OBJECT IDENTIFIER,
 *        accessLocation   GeneralName
 *    }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose Authority Information Access is fetched.
 *      Must be non-NULL.
 *  "pAuthorityInfoAccess"
 *      Address where Authority InfoAccess will be stored and returned.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetAuthorityInfoAccess(
        PKIX_PL_Cert *cert,
        PKIX_List **pAiaList, /* of PKIX_PL_InfoAccess */
        void *plContext);


/*
 * FUNCTION: PKIX_PL_Cert_GetSubjectInfoAccess
 * DESCRIPTION:
 *
 *  Retrieves the value(s) of the Subject Information Access in "cert" and
 *  returns it in a list at address pointed by "pSubjectInfoAccess".
 *
 *  SubjectInfoAccess ::=
 *    SEQUENCE SIZE (1..MAX) of AccessDescription
 *    AccessDescription ::= SEQUENCE {
 *        accessMethod     OBJECT IDENTIFIER,
 *        accessLocation   GeneralName
 *    }
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose Subject Information Access is fetched.
 *      Must be non-NULL.
 *  "pSubjectInfoAccess"
 *      Address where Subject InfoAccess will be stored and returned.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetSubjectInfoAccess(
        PKIX_PL_Cert *cert,
        PKIX_List **pSiaList, /* of PKIX_PL_InfoAccess */
        void *plContext);



/*
 * FUNCTION: PKIX_PL_Cert_GetCrlDp
 * DESCRIPTION:
 *
 *  Retrieves the value(s) of the CRL Distribution Point Extension and
 *  returns it in a list at address pointed by "pDpList".
 *
 * PARAMETERS:
 *  "cert"
 *      Address of Cert whose Subject Information Access is fetched.
 *      Must be non-NULL.
 *  "pDpList"
 *      Address where CRL DP will be stored and returned.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Cert_GetCrlDp(PKIX_PL_Cert *cert,
                      PKIX_List **pDpList,
                      void *plContext);


/*
 * InfoAccess 
 *
 * To hold Authority Information Access or Subject Information Access
 * retrieved from a Certificate.
 */

#define PKIX_INFOACCESS_OCSP          1
#define PKIX_INFOACCESS_CA_ISSUERS    2
#define PKIX_INFOACCESS_TIMESTAMPING  3
#define PKIX_INFOACCESS_CA_REPOSITORY 5

#define PKIX_INFOACCESS_LOCATION_UNKNOWN 0
#define PKIX_INFOACCESS_LOCATION_HTTP    1
#ifndef NSS_PKIX_NO_LDAP
#define PKIX_INFOACCESS_LOCATION_LDAP    2
#endif

/*
 * FUNCTION: PKIX_PL_InfoAccess_GetMethod
 * DESCRIPTION:
 *
 *  Stores the method of the Information Access from "infoAccess" and
 *  returns in "pMethod".
 *
 *  SubjectInfoAccess ::=
 *    AccessDescription ::= SEQUENCE {
 *        accessMethod     OBJECT IDENTIFIER,
 *        accessLocation   GeneralName
 *    }
 *
 * PARAMETERS:
 *  "infoAccess"
 *      Address of PKIX_PL_InfoAccess that has the access data.
 *      Must be non-NULL.
 *  "pMethod"
 *      Address where access method will be stored and returned.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_InfoAccess_GetMethod(
        PKIX_PL_InfoAccess *infoAccess,
        PKIX_UInt32 *pMethod,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_InfoAccess_GetLocation
 * DESCRIPTION:
 *
 *  Stores the location of the Information Access from "infoAccess" and
 *  returns in "pLocation".
 *
 *  SubjectInfoAccess ::=
 *    AccessDescription ::= SEQUENCE {
 *        accessMethod     OBJECT IDENTIFIER,
 *        accessLocation   GeneralName
 *    }
 *
 * PARAMETERS:
 *  "infoAccess"
 *      Address of PKIX_PL_InfoAccess that has the access data.
 *      Must be non-NULL.
 *  "pLocation"
 *      Address where access location will be stored and returned.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_InfoAccess_GetLocation(
        PKIX_PL_InfoAccess *infoAccess,
        PKIX_PL_GeneralName **pLocation,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_InfoAccess_GetLocationType
 * DESCRIPTION:
 *
 *  Stores the type of location of the Information Access from "infoAccess" and
 *  returns in "pType".
 *
 *  SubjectInfoAccess ::=
 *    AccessDescription ::= SEQUENCE {
 *        accessMethod     OBJECT IDENTIFIER,
 *        accessLocation   GeneralName
 *    }
 *
 * PARAMETERS:
 *  "infoAccess"
 *      Address of PKIX_PL_InfoAccess that has the access data.
 *      Must be non-NULL.
 *  "pType"
 *      Address where access location type will be stored and returned.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Cert Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_InfoAccess_GetLocationType(
        PKIX_PL_InfoAccess *infoAccess,
        PKIX_UInt32 *pType,
        void *plContext);

PKIX_Error *
pkix_pl_InfoAccess_GetAIACerts(
        PKIX_PL_InfoAccess *ia,
        void **pNBIOContext,
        void **pHandle,
        PKIX_List **pCerts,
        void *plContext);

/*
 * CRL
 *
 * A CRL represents an X.509 certificate revocation list. It can be created
 * using the bytes of a valid ASN.1 DER encoding. Once created, a CRL is
 * immutable. The following functions include accessors (gettors) for the
 * various components of an X.509 CRL, as well as a function for signature
 * verification.
 */

/*
 * FUNCTION: PKIX_PL_CRL_Create
 * DESCRIPTION:
 *
 *  Creates a new CRL using the bytes in the ByteArray pointed to by
 *  "byteArray" and stores it at "pCRL". If the bytes are not a valid ASN.1
 *  DER encoding of a CRL, a PKIX_Error pointer is returned. Once created, a
 *  CRL is immutable.
 *
 *  CertificateList  ::=  SEQUENCE  {
 *      tbsCertList             TBSCertList,
 *      signatureAlgorithm      AlgorithmIdentifier,
 *      signatureValue          BIT STRING  }
 *
 *  TBSCertList  ::=  SEQUENCE  {
 *      version                 Version OPTIONAL,
 *                              -- if present, MUST be v2
 *      signature               AlgorithmIdentifier,
 *      issuer                  Name,
 *      thisUpdate              Time,
 *      nextUpdate              Time OPTIONAL,
 *      revokedCertificates     SEQUENCE OF SEQUENCE  {
 *              userCertificate         CertificateSerialNumber,
 *              revocationDate          Time,
 *              crlEntryExtensions      Extensions OPTIONAL
 *                                      -- if present, MUST be v2
 *                              }  OPTIONAL,
 *      crlExtensions           [0]  EXPLICIT Extensions OPTIONAL
 *                                      -- if present, MUST be v2
 *      }
 *
 * PARAMETERS:
 *  "byteArray"
 *      Address of ByteArray representing the CRL's DER encoding.
 *      Must be non-NULL.
 *  "pCRL"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_Create(
        PKIX_PL_ByteArray *byteArray,
        PKIX_PL_CRL **pCRL,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_GetIssuer
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the X500Name that represents the issuer of the CRL
 *  pointed to by "crl" and stores it at "pCRLIssuer".
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose issuer is to be stored. Must be non-NULL.
 *  "pCRLIssuer"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_GetIssuer(
        PKIX_PL_CRL *crl,
        PKIX_PL_X500Name **pCRLIssuer,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_GetCriticalExtensionOIDs
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of OIDs (each OID corresponding to a
 *  critical extension of the CRL pointed to by "crl") and stores it at
 *  "pExtensions". If "crl" does not have any critical extensions, this
 *  function stores an empty List at "pExtensions".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose critical extension OIDs are to be stored.
 *      Must be non-NULL.
 *  "pExtensions"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_GetCriticalExtensionOIDs(
        PKIX_PL_CRL *crl,
        PKIX_List **pExtensions,   /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_GetCRLEntryForSerialNumber
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the CRLEntry (found in the CRL pointed to by "crl")
 *  corresponding to the BigInt pointed to by "serialNumber" and stores it at
 *  "pCRLEntry". If there is no such CRLEntry, this functions stores NULL at
 *  "pCRLEntry". Once created, a CRLEntry is immutable.
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose CRL Entries are to be searched. Must be non-NULL.
 *  "serialNumber"
 *      Address of BigInt representing serial number of certificate whose
 *      CRLEntry is to be found. Must be non-NULL.
 *  "pCRLEntry"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_GetCRLEntryForSerialNumber(
        PKIX_PL_CRL *crl,
        PKIX_PL_BigInt *serialNumber,
        PKIX_PL_CRLEntry **pCRLEntry,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_GetCRLNumber
 * DESCRIPTION:
 *  Retrieves the CRL Number from extension. This is non-critical extension.
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose version is to be stored. Must be non-NULL.
 *  "pCrlNumber"
 *      Address where a CRL Number will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_GetCRLNumber(
        PKIX_PL_CRL *crl,
        PKIX_PL_BigInt **pCrlNumber,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_VerifyUpdateTime
 * DESCRIPTION:
 *
 *  Checks whether the CRL pointed to by "crl" would be valid at the time
 *  represented by the Date pointed to by "date" and stores the Boolean result
 *  at "pResult". This check is done only when NIST policy is enforced.
 *
 *  Time ::= CHOICE {
 *      utcTime         UTCTime,
 *      generalTime     GeneralizedTime }
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose validity is to be checked. Must be non-NULL.
 *  "date"
 *      Address of Date at which the CRL is being checked for validity.
 *      Must be non-NULL.
 *  "pResult"
 *      Address of Boolean result. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_VerifyUpdateTime(
        PKIX_PL_CRL *crl,
        PKIX_PL_Date *date,
        PKIX_Boolean *pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_VerifySignature
 * DESCRIPTION:
 *
 *  Verifies the signature on the CRL pointed to by "crl" using the PublicKey
 *  pointed to by "pubKey". If the signature doesn't verify, a PKIX_Error
 *  pointer is returned.
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose signature is to be verified. Must be non-NULL.
 *  "pubKey"
 *      Address of a Public Key used to verify the signature. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_VerifySignature(
        PKIX_PL_CRL *crl,
        PKIX_PL_PublicKey *pubKey,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRL_ReleaseDerCrl
 * DESCRIPTION:
 *
 * Relinguish the ownership for the crl der. The operation will succeed if
 * a crl owns the der. If the crl was created from existing crl and does not
 * own the der, then the function will return null.
 *
 * PARAMETERS:
 *  "crl"
 *      Address of CRL whose signature is to be verified. Must be non-NULL.
 *  "derCrl"
 *      Pointer to a SECItem that has der crl.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_ReleaseDerCrl(PKIX_PL_CRL *crl,
                         SECItem **derCrl,
                         void *plContext);
/*
 * FUNCTION: PKIX_PL_CRL_AdoptDerCrl
 * DESCRIPTION:
 *
 * Adopt memory of the der. The secItem that contains der will be
 * freed with destruction of parent pkix crl structure.
 *
 * * PARAMETERS:
 *  "crl"
 *      Address of CRL whose signature is to be verified. Must be non-NULL.
 *  "derCrl"
 *      Pointer to a SECItem that has der crl.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRL_AdoptDerCrl(PKIX_PL_CRL *crl,
                        SECItem *derCrl,
                        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRLEntry_GetCRLEntryReasonCode
 * DESCRIPTION:
 *
 *  Retrieves the value of the reason code extension of the CRLEntry pointed
 *  to by "crlEntry" and stores it at "pReason". If the "crlEntry" has no
 *  reason code extension, this function stores -1 at "pReason".
 *
 *  CRLReason ::= ENUMERATED {
 *      unspecified             (0),
 *      keyCompromise           (1),
 *      cACompromise            (2),
 *      affiliationChanged      (3),
 *      superseded              (4),
 *      cessationOfOperation    (5),
 *      certificateHold         (6),
 *      removeFromCRL           (8),
 *      privilegeWithdrawn      (9),
 *      aACompromise            (10) }
 *
 * PARAMETERS:
 *  "crlEntry"
 *      Address of CRLEntry whose reason code bit values are to be returned
 *      at "pReason". Must be non-NULL.
 *  "pReason"
 *      Address of PKIX_Int32 where reason code is stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRLEntry_GetCRLEntryReasonCode(
        PKIX_PL_CRLEntry *crlEntry,
        PKIX_Int32 *pReason,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_CRLEntry_GetCriticalExtensionOIDs
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of OIDs (each OID corresponding to a
 *  critical extension of the CRLEntry pointed to by "crlEntry") and stores it
 *  at "pExtensions". If "crlEntry" does not have any critical extensions, this
 *  function stores an empty List at "pExtensions".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "crlEntry"
 *      Address of CRLEntry whose critical extension OIDs are to be stored.
 *      Must be non-NULL.
 *  "pExtensions"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a CRL Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CRLEntry_GetCriticalExtensionOIDs(
        PKIX_PL_CRLEntry *crlEntry,
        PKIX_List **pExtensions,  /* list of PKIX_PL_OID */
        void *plContext);

#ifdef BUILD_LIBPKIX_TESTS
/*
 * FUNCTION: PKIX_PL_X500Name_Create
 * DESCRIPTION:
 *
 *  Creates a new X500Name using the UTF8 string representation pointed to by
 *  "stringRep" and stores it at "pName". Once created, an X500Name is
 *  immutable.
 *
 *  Name ::= CHOICE {
 *    RDNSequence }
 *
 *  RDNSequence ::= SEQUENCE OF RelativeDistinguishedName
 *
 *  RelativeDistinguishedName ::=
 *    SET OF AttributeTypeAndValue
 *
 *  AttributeTypeAndValue ::= SEQUENCE {
 *      type    AttributeType,
 *      value   AttributeValue }
 *
 *  AttributeType ::= OBJECT IDENTIFIER
 *
 *  AttributeValue ::= ANY DEFINED BY AttributeType
 *
 *  DirectoryString ::= CHOICE {
 *      teletexString           TeletexString (SIZE (1..MAX)),
 *      printableString         PrintableString (SIZE (1..MAX)),
 *      universalString         UniversalString (SIZE (1..MAX)),
 *      utf8String              UTF8String (SIZE (1..MAX)),
 *      bmpString               BMPString (SIZE (1..MAX)) }
 *
 * PARAMETERS:
 *  "stringRep"
 *      Address of UTF8 String representation of X500Name. Must be non-NULL.
 *  "pName"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an X500Name Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_X500Name_Create (
        PKIX_PL_String *stringRep,
        PKIX_PL_X500Name **pName,
        void *plContext);

#endif /* BUILD_LIBPKIX_TESTS */

/*
 * FUNCTION: PKIX_PL_X500Name_CreateFromCERTName
 * DESCRIPTION:
 * 
 * The function creates x500Name using der encoded DN and/or pointer to
 * CERTName. If arument "name" is NULL, but derName is supplied when
 * the function generates nssDN(CERTName type) from der data. If derName
 * is not supplied, CERTName *name will not be used to generate DN DER
 * encoding.
 *
 * PARAMETERS:
 *  "derName"
 *      Address of DER representation of X500Name. Can be NULL
 *  "name"
 *      Address of CERTName representation of X500Name. Can be NULL
 *  "pName"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an X500Name Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_X500Name_CreateFromCERTName(
        SECItem *derName,
        CERTName *name,
        PKIX_PL_X500Name **pName,
        void *plContext);


/*
 * TYPE: PKIX_PL_X500Name_Match
 * DESCRIPTION:
 *  Checks whether the X500Name pointed to by "firstX500Name" MATCHES the
 *  X500Name pointed to by "secondX500Name" and stores the boolean result at
 *  "pResult". Two X500Names MATCH if they meet the conditions specified by
 *  RFC 3280 (section 4.1.2.4). Namely:
 *
 *      "This specification requires only a subset of the name comparison
 *      functionality specified in the X.500 series of specifications.
 *      Conforming implementations are REQUIRED to implement the following
 *      name comparison rules:
 *
 *      (a)  attribute values encoded in different types (e.g., PrintableString
 *      and BMPString) MAY be assumed to represent different strings;
 *
 *      (b) attribute values in types other than PrintableString are case
 *      sensitive (this permits matching of attribute values as binary objects)
 *
 *      (c)  attribute values in PrintableString are not case sensitive
 *      (e.g., "Marianne Swanson" is the same as "MARIANNE SWANSON"); and
 *
 *      (d)  attribute values in PrintableString are compared after removing
 *      leading and trailing white space and converting internal substrings of
 *      one or more consecutive white space characters to a single space."
 *
 * PARAMETERS:
 *  "firstX500Name"
 *      Address of first X500Name to compare. Must be non-NULL.
 *  "secondX500Name"
 *      Address of second X500Name to compare. Must be non-NULL.
 *  "pResult"
 *      Address of Boolean result. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an X500Name Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_X500Name_Match(
        PKIX_PL_X500Name *firstX500Name,
        PKIX_PL_X500Name *secondX500Name,
        PKIX_Boolean *pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Date_Create_UTCTime
 * DESCRIPTION:
 *  Creates a new Date of type UTCTime using the string representation pointed
 *  to by "stringRep" and stores it at "pDate". The UTCTime restriction means
 *  that the year can only be specified by the least significant two digits
 *  (YY). As such, Only the years 1950-2049 can be represented. If "stringRep"
 *  is NULL, this function creates a new Date representing the current time
 *  and stores it at "pDate". Once created, a Date is immutable.
 *
 *  If YY is greater than or equal to 50, the year is interpreted as 19YY.
 *  If YY is less than 50, the year is interpreted as 20YY.
 *
 *  The string representation of the date must be in the following form:
 *      "YYMMDDhhmmssZ" where:
 *
 *  YY is the least significant two digits of the year
 *  MM is the month (01 to 12)
 *  DD is the day (01 to 31)
 *  hh is the hour (00 to 23)
 *  mm are the minutes (00 to 59)
 *  ss are the seconds (00 to 59)
 *  Z indicates that local time is GMT
 *
 * PARAMETERS:
 *  "stringRep"
 *      Address of String representation of Date.
 *      If NULL, current time is used.
 *  "pDate"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Date Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Date_Create_UTCTime (
        PKIX_PL_String *stringRep,
        PKIX_PL_Date **pDate,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Date_Create_UTCTime
 * DESCRIPTION:
 *  Creates a new Date from PRTime data.
 *
 * PARAMETERS:
 *  "time"
 *      Represented time in PRTime type.
 *  "pDate"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Date Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Date_CreateFromPRTime(
        PRTime time,
        PKIX_PL_Date **pDate,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Date_Create_CurrentOffBySeconds
 * DESCRIPTION:
 *  Creates a new Date of type UTCTime for current time with seconds off by
 *  "secondsOffset" and returns it at "pDate".
 *
 * PARAMETERS:
 *  "secondsOffset"
 *      A PKIX_Int32 indicates the time offset from current. If "secondsOffset"
 *      is negative, the time is in past.
 *  "pDate"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Date Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Date_Create_CurrentOffBySeconds(
        PKIX_Int32 secondsOffset,
        PKIX_PL_Date **pDate,
        void *plContext);

#ifdef BUILD_LIBPKIX_TESTS
/*
 * FUNCTION: PKIX_PL_GeneralName_Create
 * DESCRIPTION:
 *
 *  Creates a new GeneralName of type "nameType" using the string
 *  representation pointed to by "stringRep" and stores it at "pGName".
 *  All of the GeneralName type format values specified in pkixt.h are
 *  supported, with the exception of PKIX_OTHER_NAME, PKIX_EDIPARTY_NAME,
 *  PKIX_IP_NAME, and PKIX_X400_ADDRESS. A PKIX_ESCASCII string representation
 *  should be used for all supported nameTypes, with the exception of
 *  registeredID and directoryName. For registeredID, the string representation
 *  should be the same as that used by PKIX_PL_OID_Create. For directoryName,
 *  the string representation should be the same as that used by
 *  PKIX_PL_X500Name_Create. If an unsupported name type is used, an Error is
 *  returned. Once created, a GeneralName is immutable.
 *
 *  GeneralName ::= CHOICE {
 *      otherName                       [0]     OtherName,
 *      rfc822Name                      [1]     IA5String,
 *      dNSName                         [2]     IA5String,
 *      x400Address                     [3]     ORAddress,
 *      directoryName                   [4]     Name,
 *      ediPartyName                    [5]     EDIPartyName,
 *      uniformResourceIdentifier       [6]     IA5String,
 *      iPAddress                       [7]     OCTET STRING,
 *      registeredID                    [8]     OBJECT IDENTIFIER }
 *
 *
 * NOTE: This function is allowed to be called only by pkix tests programs.
 * 
 * PARAMETERS:
 *  "nameType"
 *      Type of GeneralName to be created. This must be one of the GeneralName
 *      type format values specified in pkixt.h
 *  "stringRep"
 *      Address of String representation of GeneralName. Must be non-NULL.
 *  "pGName"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a GeneralName Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_GeneralName_Create (
        PKIX_UInt32 nameType,
        PKIX_PL_String *stringRep,
        PKIX_PL_GeneralName **pGName,
        void *plContext);
#endif /* BUILD_LIBPKIX_TESTS */

/*
 * FUNCTION: PKIX_PL_CertNameConstraints_CheckNamesInNameSpace
 * DESCRIPTION:
 *
 *  This function checks whether names in "nameList" comply with
 *  "nameConstraints". It stores PKIX_TRUE at "pCheckPass" if the names meet the
 *  requirement of the NameConstraints, PKIX_FALSE otherwise.
 *
 * PARAMETERS
 *  "nameList"
 *      List of GeneralNames that are checked for compliance. May be empty
 *      or NULL.
 *  "nameConstraints"
 *      Address of CertNameConstraints that provides lists of permitted
 *      and excluded names. Must be non-NULL.
 *  "pCheckPass"
 *      Address where PKIX_TRUE is returned if the all names in "nameList" are
 *      valid. Must be non-NULL.
 *  "plContext" - Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a NameConstraints Error if the function fails in a
 *  non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_CertNameConstraints_CheckNamesInNameSpace(
        PKIX_List *nameList, /* List of PKIX_PL_GeneralName */
        PKIX_PL_CertNameConstraints *nameConstraints,
        PKIX_Boolean *pCheckPass,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_AIAMgr_Create
 * DESCRIPTION:
 *
 *  This function creates an AIAMgr to handle retrieval of Certs and CRLs
 *  from servers given by AIA Certificate extensions. It manages connections
 *  and caches. The manager created is stored at "pAIAMgr".
 *
 * PARAMETERS:
 *  "pAIAMgr"
 *      The address at which the result is stored. Must be non-NULL.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an AIAMgr Error if the function fails in a non-fatal way
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_AIAMgr_Create(
        PKIX_PL_AIAMgr **pAIAMgr,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_AIAMgr_GetAIACerts
 * DESCRIPTION:
 *
 *  This function uses the AIAMgr pointed to by "aiaMgr" to retrieve the Certs
 *  specified by an AIA certificate extension, if any, in the Cert pointed to by
 *  "prevCert", storing the results at "pCerts". If the certificate has no such
 *  extension, this function stores NULL at "pCerts".
 *
 *  If the request is suspended for non-blocking I/O, a platform-dependent
 *  context is stored at "pNBIOContext" and NULL is stored at "pCerts". This
 *  return is referred to as the WOULDBLOCK state. Note that the caller must
 *  check for a non-NULL value at "pNBIOContext", to distinguish this state from
 *  the "no such extension" return described in the first paragraph. (The
 *  alternative would be to return an empty List, but it seemed wrong to incur
 *  the overhead of creating and destroying an empty List for the most common
 *  situation.)
 *
 *  After a WOULDBLOCK return, the user may continue the operation by calling
 *  pkix_AIAMgr_GetAIACerts (possibly more than once, if the function again
 *  returns in the WOULDBLOCK state) with the previously-returned non-NULL
 *  value of "pNBIOContext". When results are complete, NULL is stored at
 *  "pNBIOContext", and the results (which may be NULL) are stored at "pCerts".
 *
 * PARAMETERS:
 *  "aiaMgr"
 *      The AIAMgr which controls the retrieval of certificates. Must be
 *      non-NULL.
 *  "prevCert"
 *      Address of PKIX_PL_Cert which may provide an AIA or SIA extension. Must
 *      be non-NULL.
 *  "pNBIOContext"
 *      Address at which platform-dependent information is returned if request
 *      is suspended for non-blocking I/O. Must be non-NULL.
 *  "pCerts"
 *      Address at which the returned List is stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an AIAMgr Error if the function fails in a non-fatal way
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_AIAMgr_GetAIACerts(
        PKIX_PL_AIAMgr *aiaMgr,
        PKIX_PL_Cert *prevCert,
        void **pNBIOContext,
        PKIX_List **pCerts,
        void *plContext);

typedef PKIX_Error *
(*PKIX_PL_VerifyCallback)(
        PKIX_PL_Object *signedObject,
        PKIX_PL_Cert *signerCert, /* can be unknown */
        PKIX_PL_Date *producedAt,
        PKIX_ProcessingParams *procParams,
        void **pNBIOContext,
        void **pState,
        PKIX_BuildResult **pBuildResult,
        PKIX_VerifyNode **pVerifyTree,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_PKI_H */
