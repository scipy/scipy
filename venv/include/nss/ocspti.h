/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Private header defining OCSP types.
 */

#ifndef _OCSPTI_H_
#define _OCSPTI_H_

#include "ocspt.h"

#include "certt.h"
#include "plarena.h"
#include "seccomon.h"
#include "secoidt.h"

/*
 * Some notes about naming conventions...
 *
 * The public data types all start with "CERTOCSP" (e.g. CERTOCSPRequest).
 * (Even the public types are opaque, however.  Only their names are
 * "exported".)
 *
 * Internal-only data types drop the "CERT" prefix and use only the
 * lower-case "ocsp" (e.g. ocspTBSRequest), for brevity sake.
 *
 * In either case, the base/suffix of the type name usually matches the
 * name as defined in the OCSP specification.  The exceptions to this are:
 *  - When there is overlap between the "OCSP" or "ocsp" prefix and
 *    the name used in the standard.  That is, you cannot strip off the
 *    "CERTOCSP" or "ocsp" prefix and necessarily get the name of the
 *    type as it is defined in the standard; the "real" name will be
 *    *either* "OCSPSuffix" or just "Suffix".
 *  - When the name in the standard was a little too generic.  (e.g. The
 *    standard defines "Request" but we call it a "SingleRequest".)
 *    In this case a comment above the type definition calls attention
 *    to the difference.
 *
 * The definitions laid out in this header file are intended to follow
 * the same order as the definitions in the OCSP specification itself.
 * With the OCSP standard in hand, you should be able to move through
 * this file and follow along.  To future modifiers of this file: please
 * try to keep it that way.  The only exceptions are the few cases where
 * we need to define a type before it is referenced (e.g. enumerations),
 * whereas in the OCSP specification these are usually defined the other
 * way around (reference before definition).
 */

/*
 * Forward-declarations of internal-only data structures.
 *
 * These are in alphabetical order (case-insensitive); please keep it that way!
 */
typedef struct ocspBasicOCSPResponseStr ocspBasicOCSPResponse;
typedef struct ocspCertStatusStr ocspCertStatus;
typedef struct ocspResponderIDStr ocspResponderID;
typedef struct ocspResponseBytesStr ocspResponseBytes;
typedef struct ocspResponseDataStr ocspResponseData;
typedef struct ocspRevokedInfoStr ocspRevokedInfo;
typedef struct ocspServiceLocatorStr ocspServiceLocator;
typedef struct ocspSignatureStr ocspSignature;
typedef struct ocspSingleRequestStr ocspSingleRequest;
typedef struct ocspSingleResponseStr ocspSingleResponse;
typedef struct ocspTBSRequestStr ocspTBSRequest;

/*
 * An OCSPRequest; this is what is sent (encoded) to an OCSP responder.
 */
struct CERTOCSPRequestStr {
    PLArenaPool *arena; /* local; not part of encoding */
    ocspTBSRequest *tbsRequest;
    ocspSignature *optionalSignature;
};

/*
 * A TBSRequest; when an OCSPRequest is signed, the encoding of this
 * is what the signature is actually applied to.  ("TBS" == To Be Signed)
 * Whether signed or not, however, this structure will be present, and
 * is the "meat" of the OCSPRequest.
 *
 * Note that the "requestorName" field cannot be encoded/decoded in the
 * same pass as the entire request -- it needs to be handled with a special
 * call to convert to/from our internal form of a GeneralName.  Thus the
 * "derRequestorName" field, which is the actual DER-encoded bytes.
 *
 * The "extensionHandle" field is used on creation only; it holds
 * in-progress extensions as they are optionally added to the request.
 */
struct ocspTBSRequestStr {
    SECItem version;                    /* an INTEGER */
    SECItem *derRequestorName;          /* encoded GeneralName; see above */
    CERTGeneralNameList *requestorName; /* local; not part of encoding */
    ocspSingleRequest **requestList;
    CERTCertExtension **requestExtensions;
    void *extensionHandle; /* local; not part of encoding */
};

/*
 * This is the actual signature information for an OCSPRequest (applied to
 * the TBSRequest structure) or for a BasicOCSPResponse (applied to a
 * ResponseData structure).
 *
 * Note that the "signature" field itself is a BIT STRING; operations on
 * it need to keep that in mind, converting the length to bytes as needed
 * and back again afterward (so that the length is usually expressing bits).
 *
 * The "cert" field is the signer's certificate.  In the case of a received
 * signature, it will be filled in when the signature is verified.  In the
 * case of a created signature, it is filled in on creation and will be the
 * cert used to create the signature when the signing-and-encoding occurs,
 * as well as the cert (and its chain) to fill in derCerts if requested.
 *
 * The extra fields cache information about the signature after we have
 * attempted a verification.  "wasChecked", if true, means the signature
 * has been checked against the appropriate data and thus that "status"
 * contains the result of that verification.  If "status" is not SECSuccess,
 * "failureReason" is a copy of the error code that was set at the time;
 * presumably it tells why the signature verification failed.
 */
struct ocspSignatureStr {
    SECAlgorithmID signatureAlgorithm;
    SECItem signature;     /* a BIT STRING */
    SECItem **derCerts;    /* a SEQUENCE OF Certificate */
    CERTCertificate *cert; /* local; not part of encoding */
    PRBool wasChecked;     /* local; not part of encoding */
    SECStatus status;      /* local; not part of encoding */
    int failureReason;     /* local; not part of encoding */
};

/*
 * An OCSPRequest contains a SEQUENCE OF these, one for each certificate
 * whose status is being checked.
 *
 * Note that in the OCSP specification this is just called "Request",
 * but since that seemed confusing (vs. an OCSPRequest) and to be more
 * consistent with the parallel type "SingleResponse", I called it a
 * "SingleRequest".
 *
 * XXX figure out how to get rid of that arena -- there must be a way
 */
struct ocspSingleRequestStr {
    PLArenaPool *arena; /* just a copy of the response arena,
					 * needed here for extension handling
					 * routines, on creation only */
    CERTOCSPCertID *reqCert;
    CERTCertExtension **singleRequestExtensions;
};

/*
 * A CertID is the means of identifying a certificate, used both in requests
 * and in responses.
 *
 * When in a SingleRequest it specifies the certificate to be checked.
 * When in a SingleResponse it is the cert whose status is being given.
 */
struct CERTOCSPCertIDStr {
    SECAlgorithmID hashAlgorithm;
    SECItem issuerNameHash;     /* an OCTET STRING */
    SECItem issuerKeyHash;      /* an OCTET STRING */
    SECItem serialNumber;       /* an INTEGER */
    SECItem issuerSHA1NameHash; /* keep other hashes around when */
    SECItem issuerMD5NameHash;  /* we have them */
    SECItem issuerMD2NameHash;
    SECItem issuerSHA1KeyHash; /* keep other hashes around when */
    SECItem issuerMD5KeyHash;  /* we have them */
    SECItem issuerMD2KeyHash;
    PLArenaPool *poolp;
};

/*
 * This describes the value of the responseStatus field in an OCSPResponse.
 * The corresponding ASN.1 definition is:
 *
 * OCSPResponseStatus	::=	ENUMERATED {
 *	successful		(0),	--Response has valid confirmations
 *	malformedRequest	(1),	--Illegal confirmation request
 *	internalError		(2),	--Internal error in issuer
 *	tryLater		(3),	--Try again later
 *					--(4) is not used
 *	sigRequired		(5),	--Must sign the request
 *	unauthorized		(6),	--Request unauthorized
 * }
 */
typedef enum {
    ocspResponse_min = 0,
    ocspResponse_successful = 0,
    ocspResponse_malformedRequest = 1,
    ocspResponse_internalError = 2,
    ocspResponse_tryLater = 3,
    ocspResponse_unused = 4,
    ocspResponse_sigRequired = 5,
    ocspResponse_unauthorized = 6,
    ocspResponse_max = 6 /* Please update max when adding values.
                          * Remember to also update arrays, e.g.
                          * "responseStatusNames" in ocspclnt.c
                          * and potentially other places. */
} ocspResponseStatus;

/*
 * An OCSPResponse is what is sent (encoded) by an OCSP responder.
 *
 * The field "responseStatus" is the ASN.1 encoded value; the field
 * "statusValue" is simply that same value translated into our local
 * type ocspResponseStatus.
 */
struct CERTOCSPResponseStr {
    PLArenaPool *arena;               /* local; not part of encoding */
    SECItem responseStatus;           /* an ENUMERATED, see above */
    ocspResponseStatus statusValue;   /* local; not part of encoding */
    ocspResponseBytes *responseBytes; /* only when status is successful */
};

/*
 * A ResponseBytes (despite appearances) is what contains the meat
 * of a successful response -- but still in encoded form.  The type
 * given as "responseType" tells you how to decode the string.
 *
 * We look at the OID and translate it into our local OID representation
 * "responseTypeTag", and use that value to tell us how to decode the
 * actual response itself.  For now the only kind of OCSP response we
 * know about is a BasicOCSPResponse.  However, the intention in the
 * OCSP specification is to allow for other response types, so we are
 * building in that flexibility from the start and thus put a pointer
 * to that data structure inside of a union.  Whenever OCSP adds more
 * response types, just add them to the union.
 */
struct ocspResponseBytesStr {
    SECItem responseType;      /* an OBJECT IDENTIFIER */
    SECOidTag responseTypeTag; /* local; not part of encoding */
    SECItem response;          /* an OCTET STRING */
    union {
        ocspBasicOCSPResponse *basic; /* when type is id-pkix-ocsp-basic */
    } decodedResponse;                /* local; not part of encoding */
};

/*
 * A BasicOCSPResponse -- when the responseType in a ResponseBytes is
 * id-pkix-ocsp-basic, the "response" OCTET STRING above is the DER
 * encoding of one of these.
 *
 * Note that in the OCSP specification, the signature fields are not
 * part of a separate sub-structure.  But since they are the same fields
 * as we define for the signature in a request, it made sense to share
 * the C data structure here and in some shared code to operate on them.
 */
struct ocspBasicOCSPResponseStr {
    SECItem tbsResponseDataDER;
    ocspResponseData *tbsResponseData; /* "tbs" == To Be Signed */
    ocspSignature responseSignature;
};

/*
 * A ResponseData is the part of a BasicOCSPResponse that is signed
 * (after it is DER encoded).  It contains the real details of the response
 * (a per-certificate status).
 */
struct ocspResponseDataStr {
    SECItem version; /* an INTEGER */
    SECItem derResponderID;
    ocspResponderID *responderID; /* local; not part of encoding */
    SECItem producedAt;           /* a GeneralizedTime */
    CERTOCSPSingleResponse **responses;
    CERTCertExtension **responseExtensions;
};

struct ocspResponderIDStr {
    CERTOCSPResponderIDType responderIDType; /* local; not part of encoding */
    union {
        CERTName name;   /* when ocspResponderID_byName */
        SECItem keyHash; /* when ocspResponderID_byKey */
        SECItem other;   /* when ocspResponderID_other */
    } responderIDValue;
};

/*
 * The ResponseData in a BasicOCSPResponse contains a SEQUENCE OF
 * SingleResponse -- one for each certificate whose status is being supplied.
 *
 * XXX figure out how to get rid of that arena -- there must be a way
 */
struct CERTOCSPSingleResponseStr {
    PLArenaPool *arena; /* just a copy of the response arena,
					 * needed here for extension handling
					 * routines, on creation only */
    CERTOCSPCertID *certID;
    SECItem derCertStatus;
    ocspCertStatus *certStatus; /* local; not part of encoding */
    SECItem thisUpdate;         /* a GeneralizedTime */
    SECItem *nextUpdate;        /* a GeneralizedTime */
    CERTCertExtension **singleExtensions;
};

/*
 * A CertStatus is the actual per-certificate status.  Its ASN.1 definition:
 *
 * CertStatus	::=	CHOICE {
 *	good			[0] IMPLICIT NULL,
 *	revoked			[1] IMPLICIT RevokedInfo,
 *	unknown			[2] IMPLICIT UnknownInfo }
 *
 * (where for now UnknownInfo is defined to be NULL but in the
 * future may be replaced with an enumeration).
 *
 * Because it is CHOICE, the status value and its associated information
 * (if any) are actually encoded together.  To represent this same
 * information internally, we explicitly define a type and save it,
 * along with the value, into a data structure.
 */

typedef enum {
    ocspCertStatus_good,    /* cert is not revoked */
    ocspCertStatus_revoked, /* cert is revoked */
    ocspCertStatus_unknown, /* cert was unknown to the responder */
    ocspCertStatus_other    /* status was not an expected value */
} ocspCertStatusType;

/*
 * This is the actual per-certificate status.
 *
 * The "goodInfo" and "unknownInfo" items are only place-holders for a NULL.
 * (Though someday OCSP may replace UnknownInfo with an enumeration that
 * gives more detailed information.)
 */
struct ocspCertStatusStr {
    ocspCertStatusType certStatusType; /* local; not part of encoding */
    union {
        SECItem *goodInfo;            /* when ocspCertStatus_good */
        ocspRevokedInfo *revokedInfo; /* when ocspCertStatus_revoked */
        SECItem *unknownInfo;         /* when ocspCertStatus_unknown */
        SECItem *otherInfo;           /* when ocspCertStatus_other */
    } certStatusInfo;
};

/*
 * A RevokedInfo gives information about a revoked certificate -- when it
 * was revoked and why.
 */
struct ocspRevokedInfoStr {
    SECItem revocationTime;    /* a GeneralizedTime */
    SECItem *revocationReason; /* a CRLReason; ignored for now */
};

/*
 * ServiceLocator can be included as one of the singleRequestExtensions.
 * When added, it specifies the (name of the) issuer of the cert being
 * checked, and optionally the value of the AuthorityInfoAccess extension
 * if the cert has one.
 */
struct ocspServiceLocatorStr {
    CERTName *issuer;
    SECItem locator; /* DER encoded authInfoAccess extension from cert */
};

#endif /* _OCSPTI_H_ */
