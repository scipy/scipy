/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * ocspi.h - NSS internal interfaces to OCSP code
 */

#ifndef _OCSPI_H_
#define _OCSPI_H_

SECStatus OCSP_InitGlobal(void);
SECStatus OCSP_ShutdownGlobal(void);

ocspResponseData *
ocsp_GetResponseData(CERTOCSPResponse *response, SECItem **tbsResponseDataDER);

ocspSignature *
ocsp_GetResponseSignature(CERTOCSPResponse *response);

SECItem *
ocsp_DigestValue(PLArenaPool *arena, SECOidTag digestAlg,
                 SECItem *fill, const SECItem *src);

PRBool
ocsp_CertIsOCSPDefaultResponder(CERTCertDBHandle *handle, CERTCertificate *cert);

CERTCertificate *
ocsp_GetSignerCertificate(CERTCertDBHandle *handle, ocspResponseData *tbsData,
                          ocspSignature *signature, CERTCertificate *issuer);

SECStatus
ocsp_VerifyResponseSignature(CERTCertificate *signerCert,
                             ocspSignature *signature,
                             SECItem *tbsResponseDataDER,
                             void *pwArg);

CERTOCSPRequest *
cert_CreateSingleCertOCSPRequest(CERTOCSPCertID *certID,
                                 CERTCertificate *singleCert,
                                 PRTime time,
                                 PRBool addServiceLocator,
                                 CERTCertificate *signerCert);

typedef enum { ocspMissing,
               ocspFresh,
               ocspStale } OCSPFreshness;

SECStatus
ocsp_GetCachedOCSPResponseStatus(CERTOCSPCertID *certID,
                                 PRTime time,
                                 PRBool ignoreOcspFailureMode,
                                 SECStatus *rvOcsp,
                                 SECErrorCodes *missingResponseError,
                                 OCSPFreshness *freshness);

/*
 * FUNCTION: cert_ProcessOCSPResponse
 *  Same behavior and basic parameters as CERT_GetOCSPStatusForCertID.
 *  In addition it can update the OCSP cache (using information
 *  available internally to this function).
 * INPUTS:
 *  CERTCertDBHandle *handle
 *    certificate DB of the cert that is being checked
 *  CERTOCSPResponse *response
 *    the OCSP response we want to retrieve status from.
 *  CERTOCSPCertID *certID
 *    the ID we want to look for from the response.
 *  CERTCertificate *signerCert
 *    the certificate that was used to sign the OCSP response.
 *    must be obtained via a call to CERT_VerifyOCSPResponseSignature.
 *  PRTime time
 *    The time at which we're checking the status for.
 *  PRBool *certIDWasConsumed
 *    In and Out parameter.
 *    If certIDWasConsumed is NULL on input,
 *    this function might produce a deep copy of cert ID
 *    for storing it in the cache.
 *    If out value is true, ownership of parameter certID was
 *    transferred to the OCSP cache.
 *  SECStatus *cacheUpdateStatus
 *    This optional out parameter will contain the result
 *    of the cache update operation (if requested).
 *  RETURN:
 *    The return value is not influenced by the cache operation,
 *    it matches the documentation for CERT_CheckOCSPStatus
 */

SECStatus
cert_ProcessOCSPResponse(CERTCertDBHandle *handle,
                         CERTOCSPResponse *response,
                         CERTOCSPCertID *certID,
                         CERTCertificate *signerCert,
                         PRTime time,
                         PRBool *certIDWasConsumed,
                         SECStatus *cacheUpdateStatus);

/*
 * FUNCTION: cert_RememberOCSPProcessingFailure
 *  If an application notices a failure during OCSP processing,
 *  it should finally call this function. The failure will be recorded
 *  in the OCSP cache in order to avoid repetitive failures.
 * INPUTS:
 *  CERTOCSPCertID *certID
 *    the ID that was used for the failed OCSP processing
 *  PRBool *certIDWasConsumed
 *    Out parameter, if set to true, ownership of parameter certID was
 *    transferred to the OCSP cache.
 *  RETURN:
 *    Status of the cache update operation.
 */

SECStatus
cert_RememberOCSPProcessingFailure(CERTOCSPCertID *certID,
                                   PRBool *certIDWasConsumed);

/*
 * FUNCTION: ocsp_GetResponderLocation
 *  Check ocspx context for user-designated responder URI first. If not
 *  found, checks cert AIA extension.
 * INPUTS:
 *  CERTCertDBHandle *handle
 *    certificate DB of the cert that is being checked
 *  CERTCertificate *cert
 *     The certificate being examined.
 *  PRBool *certIDWasConsumed
 *    Out parameter, if set to true, URI of default responder is
 *    returned.
 *  RETURN:
 *    Responder URI.
 */
char *
ocsp_GetResponderLocation(CERTCertDBHandle *handle,
                          CERTCertificate *cert,
                          PRBool canUseDefaultLocation,
                          PRBool *isDefault);

/* FUNCTION: ocsp_FetchingFailureIsVerificationFailure
 * The function checks the global ocsp settings and
 * tells how to treat an ocsp response fetching failure.
 * RETURNS:
 *   if PR_TRUE is returned, then treat fetching as a
 *   revoked cert status.
 */
PRBool
ocsp_FetchingFailureIsVerificationFailure(void);

size_t
ocsp_UrlEncodeBase64Buf(const char *base64Buf, char *outputBuf);

SECStatus
ocsp_GetVerifiedSingleResponseForCertID(CERTCertDBHandle *handle,
                                        CERTOCSPResponse *response,
                                        CERTOCSPCertID *certID,
                                        CERTCertificate *signerCert,
                                        PRTime time,
                                        CERTOCSPSingleResponse **pSingleResponse);

SECStatus
ocsp_CertHasGoodStatus(ocspCertStatus *status, PRTime time);

void
ocsp_CacheSingleResponse(CERTOCSPCertID *certID,
                         CERTOCSPSingleResponse *single,
                         PRBool *certIDWasConsumed);

#endif /* _OCSPI_H_ */
