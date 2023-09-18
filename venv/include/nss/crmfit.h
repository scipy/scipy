/* -*- Mode: C; tab-width: 8 -*-*/
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _CRMFIT_H_
#define _CRMFIT_H_

struct CRMFCertReqMessagesStr {
    CRMFCertReqMsg **messages;
    PLArenaPool *poolp;
};

struct CRMFCertExtensionStr {
    SECItem id;
    SECItem critical;
    SECItem value;
};

struct CRMFOptionalValidityStr {
    SECItem notBefore;
    SECItem notAfter;
};

struct CRMFCertTemplateStr {
    SECItem version;
    SECItem serialNumber;
    SECAlgorithmID *signingAlg;
    CERTName *issuer;
    CRMFOptionalValidity *validity;
    CERTName *subject;
    CERTSubjectPublicKeyInfo *publicKey;
    SECItem issuerUID;
    SECItem subjectUID;
    CRMFCertExtension **extensions;
    int numExtensions;
};

struct CRMFCertIDStr {
    SECItem issuer;       /* General Name */
    SECItem serialNumber; /*INTEGER*/
};

struct CRMFEncryptedValueStr {
    SECAlgorithmID *intendedAlg;
    SECAlgorithmID *symmAlg;
    SECItem encSymmKey; /*BIT STRING   */
    SECAlgorithmID *keyAlg;
    SECItem valueHint; /*OCTET STRING */
    SECItem encValue;  /*BIT STRING   */
};

/*
 * The field derValue will contain the actual der
 * to include in the encoding or that was read in
 * from a der blob.
 */
struct CRMFEncryptedKeyStr {
    union {
        SEC_PKCS7ContentInfo *envelopedData;
        CRMFEncryptedValue encryptedValue;
    } value;
    CRMFEncryptedKeyChoice encKeyChoice;
    SECItem derValue;
};

/* ASN1 must only have one of the following 3 options. */
struct CRMFPKIArchiveOptionsStr {
    union {
        CRMFEncryptedKey encryptedKey;
        SECItem keyGenParameters;
        SECItem archiveRemGenPrivKey; /* BOOLEAN */
    } option;
    CRMFPKIArchiveOptionsType archOption;
};

struct CRMFPKIPublicationInfoStr {
    SECItem action; /* Possible values                    */
                    /* dontPublish (0), pleasePublish (1) */
    CRMFSinglePubInfo **pubInfos;
};

struct CRMFControlStr {
    SECOidTag tag;
    SECItem derTag;
    SECItem derValue;
    /* These will be C structures used to represent the various
     * options.  Values that can't be stored as der right away.
     * After creating these structures, we'll place their der
     * encoding in derValue so the encoder knows how to get to
     * it.
     */
    union {
        CRMFCertID oldCertId;
        CRMFPKIArchiveOptions archiveOptions;
        CRMFPKIPublicationInfo pubInfo;
        CRMFProtocolEncrKey protEncrKey;
    } value;
};

struct CRMFCertRequestStr {
    SECItem certReqId;
    CRMFCertTemplate certTemplate;
    CRMFControl **controls;
    /* The following members are used by the internal implementation, but
     * are not part of the encoding.
     */
    PLArenaPool *poolp;
    PRUint32 requestID; /* This is the value that will be encoded into
                         * the certReqId field.
                         */
};

struct CRMFAttributeStr {
    SECItem derTag;
    SECItem derValue;
};

struct CRMFCertReqMsgStr {
    CRMFCertRequest *certReq;
    CRMFProofOfPossession *pop;
    CRMFAttribute **regInfo;
    SECItem derPOP;
    /* This arena will be used for allocating memory when decoding.
     */
    PLArenaPool *poolp;
    PRBool isDecoded;
};

struct CRMFPOPOSigningKeyInputStr {
    /* ASN1 must have only one of the next 2 options */
    union {
        SECItem sender; /*General Name*/
        CRMFPKMACValue *publicKeyMAC;
    } authInfo;
    CERTSubjectPublicKeyInfo publicKey;
};

struct CRMFPOPOSigningKeyStr {
    SECItem derInput; /*If in the future we support
                       *POPOSigningKeyInput, this will
                       *a C structure representation
                       *instead.
                       */
    SECAlgorithmID *algorithmIdentifier;
    SECItem signature; /* This is a BIT STRING. Remember */
};                     /* that when interpreting.        */

/* ASN1 must only choose one of these members */
struct CRMFPOPOPrivKeyStr {
    union {
        SECItem thisMessage;       /* BIT STRING */
        SECItem subsequentMessage; /*INTEGER*/
        SECItem dhMAC;             /*BIT STRING*/
    } message;
    CRMFPOPOPrivKeyChoice messageChoice;
};

/* ASN1 must only have one of these options. */
struct CRMFProofOfPossessionStr {
    union {
        SECItem raVerified;
        CRMFPOPOSigningKey signature;
        CRMFPOPOPrivKey keyEncipherment;
        CRMFPOPOPrivKey keyAgreement;
    } popChoice;
    CRMFPOPChoice popUsed; /*Not part of encoding*/
};

struct CRMFPKMACValueStr {
    SECAlgorithmID algID;
    SECItem value; /*BIT STRING*/
};

struct CRMFSinglePubInfoStr {
    SECItem pubMethod;            /* Possible Values:
                                   *   dontCare (0)
                                   *   x500     (1)
                                   *   web      (2)
                                   *   ldap     (3)
                                   */
    CERTGeneralName *pubLocation; /* General Name */
};

#endif /* _CRMFIT_H_ */
