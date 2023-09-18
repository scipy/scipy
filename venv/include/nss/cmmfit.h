/* -*- Mode: C; tab-width: 8 -*-*/
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _CMMFIT_H_
#define _CMMFIT_H_

/*
 * All fields marked by a PKIStausInfo in comments is an integer
 * with the following possible values.
 *
 *  Integer Value          Meaning
 *  -------------          -------
 *         0               granted- got exactly what you asked for.
 *
 *         1               grantedWithMods-got something like what you asked
 *                          for;requester is responsible for ascertainging the
 *                          differences.
 *
 *         2               rejection-you don't get what you asked for; more
 *                          information elsewhere in the message
 *
 *         3               waiting-the request body part has not yet been
 *                          processed, expect to hear more later.
 *
 *         4               revocationWarning-this message contains a warning
 *                          that a revocation is imminent.
 *
 *         5               revocationNotification-notification that a
 *                          revocation has occurred.
 *
 *         6               keyUpdateWarning-update already done for the
 *                          oldCertId specified in FullCertTemplate.
 */

struct CMMFPKIStatusInfoStr {
    SECItem status;
    SECItem statusString;
    SECItem failInfo;
};

struct CMMFCertOrEncCertStr {
    union {
        CERTCertificate *certificate;
        CRMFEncryptedValue *encryptedCert;
    } cert;
    CMMFCertOrEncCertChoice choice;
    SECItem derValue;
};

struct CMMFCertifiedKeyPairStr {
    CMMFCertOrEncCert certOrEncCert;
    CRMFEncryptedValue *privateKey;
    SECItem derPublicationInfo; /* We aren't creating
                                 * PKIPublicationInfo's, so
                                 * we'll store away the der
                                 * here if we decode one that
                                 * does have pubInfo.
                                 */
    SECItem unwrappedPrivKey;
};

struct CMMFCertResponseStr {
    SECItem certReqId;
    CMMFPKIStatusInfo status; /*PKIStatusInfo*/
    CMMFCertifiedKeyPair *certifiedKeyPair;
};

struct CMMFCertRepContentStr {
    CERTCertificate **caPubs;
    CMMFCertResponse **response;
    PLArenaPool *poolp;
    PRBool isDecoded;
};

struct CMMFChallengeStr {
    SECAlgorithmID *owf;
    SECItem witness;
    SECItem senderDER;
    SECItem key;
    SECItem challenge;
    SECItem randomNumber;
};

struct CMMFRandStr {
    SECItem integer;
    SECItem senderHash;
    CERTGeneralName *sender;
};

struct CMMFPOPODecKeyChallContentStr {
    CMMFChallenge **challenges;
    PLArenaPool *poolp;
    int numChallenges;
    int numAllocated;
};

struct CMMFPOPODecKeyRespContentStr {
    SECItem **responses;
    PLArenaPool *poolp;
};

struct CMMFKeyRecRepContentStr {
    CMMFPKIStatusInfo status; /* PKIStatusInfo */
    CERTCertificate *newSigCert;
    CERTCertificate **caCerts;
    CMMFCertifiedKeyPair **keyPairHist;
    PLArenaPool *poolp;
    int numKeyPairs;
    int allocKeyPairs;
    PRBool isDecoded;
};

#endif /* _CMMFIT_H_ */
