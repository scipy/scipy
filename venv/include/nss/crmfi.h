/* -*- Mode: C; tab-width: 8 -*-*/
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _CRMFI_H_
#define _CRMFI_H_
/* This file will contain all declarations common to both
 * encoding and decoding of CRMF Cert Requests.  This header
 * file should only be included internally by CRMF implementation
 * files.
 */
#include "secasn1.h"
#include "crmfit.h"
#include "secerr.h"
#include "blapit.h"

#define CRMF_DEFAULT_ARENA_SIZE 1024

/*
 * Explanation for the definition of MAX_WRAPPED_KEY_LEN:
 *
 * It's used for internal buffers to transport a wrapped private key.
 * The value is in BYTES.
 * We want to define a reasonable upper bound for this value.
 * Ideally this could be calculated, but in order to simplify the code
 * we want to estimate the maximum requires size.
 * See also bug 655850 for the full explanation.
 *
 * We know the largest wrapped keys are RSA keys.
 * We'll estimate the maximum size needed for wrapped RSA keys,
 * and assume it's sufficient for wrapped keys of any type we support.
 *
 * The maximum size of RSA keys in bits is defined elsewhere as
 *   RSA_MAX_MODULUS_BITS
 *
 * The idea is to define MAX_WRAPPED_KEY_LEN based on the above.
 *
 * A wrapped RSA key requires about
 *   ( ( RSA_MAX_MODULUS_BITS / 8 ) * 5.5) + 65
 * bytes.
 *
 * Therefore, a safe upper bound is:
 *   ( ( RSA_MAX_MODULUS_BITS / 8 ) *8 ) = RSA_MAX_MODULUS_BITS
 *
 */
#define MAX_WRAPPED_KEY_LEN RSA_MAX_MODULUS_BITS

#define CRMF_BITS_TO_BYTES(bits) (((bits) + 7) / 8)
#define CRMF_BYTES_TO_BITS(bytes) ((bytes)*8)

struct crmfEncoderArg {
    SECItem *buffer;
    unsigned long allocatedLen;
};

struct crmfEncoderOutput {
    CRMFEncoderOutputCallback fn;
    void *outputArg;
};

/*
 * This function is used by the API for encoding functions that are
 * exposed through the API, ie all of the CMMF_Encode* and CRMF_Encode*
 * functions.
 */
extern void
crmf_encoder_out(void *arg, const char *buf, unsigned long len,
                 int depth, SEC_ASN1EncodingPart data_kind);

/*
 * This function is used when we want to encode something locally within
 * the library, ie the CertRequest so that we can produce its signature.
 */
extern SECStatus
crmf_init_encoder_callback_arg(struct crmfEncoderArg *encoderArg,
                               SECItem *derDest);

/*
 * This is the callback function we feed to the ASN1 encoder when doing
 * internal DER-encodings.  ie, encoding the cert request so we can
 * produce a signature.
 */
extern void
crmf_generic_encoder_callback(void *arg, const char *buf, unsigned long len,
                              int depth, SEC_ASN1EncodingPart data_kind);

/* The ASN1 templates that need to be seen by internal files
 * in order to implement CRMF.
 */
extern const SEC_ASN1Template CRMFCertReqMsgTemplate[];
extern const SEC_ASN1Template CRMFRAVerifiedTemplate[];
extern const SEC_ASN1Template CRMFPOPOSigningKeyTemplate[];
extern const SEC_ASN1Template CRMFPOPOKeyEnciphermentTemplate[];
extern const SEC_ASN1Template CRMFPOPOKeyAgreementTemplate[];
extern const SEC_ASN1Template CRMFThisMessageTemplate[];
extern const SEC_ASN1Template CRMFSubsequentMessageTemplate[];
extern const SEC_ASN1Template CRMFDHMACTemplate[];
extern const SEC_ASN1Template CRMFEncryptedKeyWithEncryptedValueTemplate[];
extern const SEC_ASN1Template CRMFEncryptedValueTemplate[];

/*
 * Use these two values for encoding Boolean values.
 */
extern const unsigned char hexTrue;
extern const unsigned char hexFalse;
/*
 * Prototypes for helper routines used internally by multiple files.
 */
extern SECStatus crmf_encode_integer(PLArenaPool *poolp, SECItem *dest,
                                     long value);
extern SECStatus crmf_make_bitstring_copy(PLArenaPool *arena, SECItem *dest,
                                          SECItem *src);

extern SECStatus crmf_copy_pkiarchiveoptions(PLArenaPool *poolp,
                                             CRMFPKIArchiveOptions *destOpt,
                                             CRMFPKIArchiveOptions *srcOpt);
extern SECStatus
crmf_destroy_pkiarchiveoptions(CRMFPKIArchiveOptions *inArchOptions,
                               PRBool freeit);
extern const SEC_ASN1Template *
crmf_get_pkiarchiveoptions_subtemplate(CRMFControl *inControl);

extern SECStatus crmf_copy_encryptedkey(PLArenaPool *poolp,
                                        CRMFEncryptedKey *srcEncrKey,
                                        CRMFEncryptedKey *destEncrKey);
extern SECStatus
crmf_copy_encryptedvalue(PLArenaPool *poolp,
                         CRMFEncryptedValue *srcValue,
                         CRMFEncryptedValue *destValue);

extern SECStatus
crmf_copy_encryptedvalue_secalg(PLArenaPool *poolp,
                                SECAlgorithmID *srcAlgId,
                                SECAlgorithmID **destAlgId);

extern SECStatus crmf_template_copy_secalg(PLArenaPool *poolp,
                                           SECAlgorithmID **dest,
                                           SECAlgorithmID *src);

extern SECStatus crmf_copy_cert_name(PLArenaPool *poolp, CERTName **dest,
                                     CERTName *src);

extern SECStatus crmf_template_add_public_key(PLArenaPool *poolp,
                                              CERTSubjectPublicKeyInfo **dest,
                                              CERTSubjectPublicKeyInfo *pubKey);

extern CRMFCertExtension *crmf_create_cert_extension(PLArenaPool *poolp,
                                                     SECOidTag tag,
                                                     PRBool isCritical,
                                                     SECItem *data);
extern CRMFCertRequest *
crmf_copy_cert_request(PLArenaPool *poolp, CRMFCertRequest *srcReq);

extern SECStatus crmf_destroy_encrypted_value(CRMFEncryptedValue *inEncrValue,
                                              PRBool freeit);

extern CRMFEncryptedValue *
crmf_create_encrypted_value_wrapped_privkey(SECKEYPrivateKey *inPrivKey,
                                            SECKEYPublicKey *inPubKey,
                                            CRMFEncryptedValue *destValue);

extern CK_MECHANISM_TYPE
crmf_get_mechanism_from_public_key(SECKEYPublicKey *inPubKey);

extern SECStatus
crmf_encrypted_value_unwrap_priv_key(PLArenaPool *poolp,
                                     CRMFEncryptedValue *encValue,
                                     SECKEYPrivateKey *privKey,
                                     SECKEYPublicKey *newPubKey,
                                     SECItem *nickname,
                                     PK11SlotInfo *slot,
                                     unsigned char keyUsage,
                                     SECKEYPrivateKey **unWrappedKey,
                                     void *wincx);

extern SECItem *
crmf_get_public_value(SECKEYPublicKey *pubKey, SECItem *dest);

extern CRMFCertExtension *
crmf_copy_cert_extension(PLArenaPool *poolp, CRMFCertExtension *inExtension);

extern SECStatus
crmf_create_prtime(SECItem *src, PRTime **dest);
#endif /*_CRMFI_H_*/
