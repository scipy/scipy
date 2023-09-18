/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Support routines for CMS implementation, none of which are exported.
 *
 * Do not export this file!  If something in here is really needed outside
 * of smime code, first try to add a CMS interface which will do it for
 * you.  If that has a problem, then just move out what you need, changing
 * its name as appropriate!
 */

#ifndef _CMSLOCAL_H_
#define _CMSLOCAL_H_

#include "cms.h"
#include "cmsreclist.h"
#include "secasn1t.h"

extern const SEC_ASN1Template NSSCMSContentInfoTemplate[];

struct NSSCMSContentInfoPrivateStr {
    NSSCMSCipherContext *ciphcx;
    NSSCMSDigestContext *digcx;
    PRBool dontStream;
};

/************************************************************************/
SEC_BEGIN_PROTOS

/*
 * private content Info stuff
 */

/* initialize the private content info field. If this returns
 * SECSuccess, the cinfo->private field is safe to dereference.
 */
SECStatus NSS_CMSContentInfo_Private_Init(NSSCMSContentInfo *cinfo);

/***********************************************************************
 * cmscipher.c - en/decryption routines
 ***********************************************************************/

/*
 * NSS_CMSCipherContext_StartDecrypt - create a cipher context to do decryption
 * based on the given bulk * encryption key and algorithm identifier (which may include an iv).
 */
extern NSSCMSCipherContext *
NSS_CMSCipherContext_StartDecrypt(PK11SymKey *key, SECAlgorithmID *algid);

/*
 * NSS_CMSCipherContext_StartEncrypt - create a cipher object to do encryption,
 * based on the given bulk encryption key and algorithm tag.  Fill in the algorithm
 * identifier (which may include an iv) appropriately.
 */
extern NSSCMSCipherContext *
NSS_CMSCipherContext_StartEncrypt(PLArenaPool *poolp, PK11SymKey *key, SECAlgorithmID *algid);

extern void
NSS_CMSCipherContext_Destroy(NSSCMSCipherContext *cc);

/*
 * NSS_CMSCipherContext_DecryptLength - find the output length of the next call to decrypt.
 *
 * cc - the cipher context
 * input_len - number of bytes used as input
 * final - true if this is the final chunk of data
 *
 * Result can be used to perform memory allocations.  Note that the amount
 * is exactly accurate only when not doing a block cipher or when final
 * is false, otherwise it is an upper bound on the amount because until
 * we see the data we do not know how many padding bytes there are
 * (always between 1 and bsize).
 */
extern unsigned int
NSS_CMSCipherContext_DecryptLength(NSSCMSCipherContext *cc, unsigned int input_len, PRBool final);

/*
 * NSS_CMSCipherContext_EncryptLength - find the output length of the next call to encrypt.
 *
 * cc - the cipher context
 * input_len - number of bytes used as input
 * final - true if this is the final chunk of data
 *
 * Result can be used to perform memory allocations.
 */
extern unsigned int
NSS_CMSCipherContext_EncryptLength(NSSCMSCipherContext *cc, unsigned int input_len, PRBool final);

/*
 * NSS_CMSCipherContext_Decrypt - do the decryption
 *
 * cc - the cipher context
 * output - buffer for decrypted result bytes
 * output_len_p - number of bytes in output
 * max_output_len - upper bound on bytes to put into output
 * input - pointer to input bytes
 * input_len - number of input bytes
 * final - true if this is the final chunk of data
 *
 * Decrypts a given length of input buffer (starting at "input" and
 * containing "input_len" bytes), placing the decrypted bytes in
 * "output" and storing the output length in "*output_len_p".
 * "cc" is the return value from NSS_CMSCipher_StartDecrypt.
 * When "final" is true, this is the last of the data to be decrypted.
 */
extern SECStatus
NSS_CMSCipherContext_Decrypt(NSSCMSCipherContext *cc, unsigned char *output,
                             unsigned int *output_len_p, unsigned int max_output_len,
                             const unsigned char *input, unsigned int input_len,
                             PRBool final);

/*
 * NSS_CMSCipherContext_Encrypt - do the encryption
 *
 * cc - the cipher context
 * output - buffer for decrypted result bytes
 * output_len_p - number of bytes in output
 * max_output_len - upper bound on bytes to put into output
 * input - pointer to input bytes
 * input_len - number of input bytes
 * final - true if this is the final chunk of data
 *
 * Encrypts a given length of input buffer (starting at "input" and
 * containing "input_len" bytes), placing the encrypted bytes in
 * "output" and storing the output length in "*output_len_p".
 * "cc" is the return value from NSS_CMSCipher_StartEncrypt.
 * When "final" is true, this is the last of the data to be encrypted.
 */
extern SECStatus
NSS_CMSCipherContext_Encrypt(NSSCMSCipherContext *cc, unsigned char *output,
                             unsigned int *output_len_p, unsigned int max_output_len,
                             const unsigned char *input, unsigned int input_len,
                             PRBool final);

/************************************************************************
 * cmspubkey.c - public key operations
 ************************************************************************/

/*
 * NSS_CMSUtil_EncryptSymKey_RSA - wrap a symmetric key with RSA
 *
 * this function takes a symmetric key and encrypts it using an RSA public key
 * according to PKCS#1 and RFC2633 (S/MIME)
 */
extern SECStatus
NSS_CMSUtil_EncryptSymKey_RSA(PLArenaPool *poolp, CERTCertificate *cert,
                              PK11SymKey *key,
                              SECItem *encKey);

extern SECStatus
NSS_CMSUtil_EncryptSymKey_RSAPubKey(PLArenaPool *poolp,
                                    SECKEYPublicKey *publickey,
                                    PK11SymKey *bulkkey, SECItem *encKey);

/*
 * NSS_CMSUtil_DecryptSymKey_RSA - unwrap a RSA-wrapped symmetric key
 *
 * this function takes an RSA-wrapped symmetric key and unwraps it, returning a symmetric
 * key handle. Please note that the actual unwrapped key data may not be allowed to leave
 * a hardware token...
 */
extern PK11SymKey *
NSS_CMSUtil_DecryptSymKey_RSA(SECKEYPrivateKey *privkey, SECItem *encKey, SECOidTag bulkalgtag);

extern SECStatus
NSS_CMSUtil_EncryptSymKey_ESDH(PLArenaPool *poolp, CERTCertificate *cert, PK11SymKey *key,
                               SECItem *encKey, SECItem **ukm, SECAlgorithmID *keyEncAlg,
                               SECItem *originatorPubKey);

extern PK11SymKey *
NSS_CMSUtil_DecryptSymKey_ESDH(SECKEYPrivateKey *privkey, SECItem *encKey,
                               SECAlgorithmID *keyEncAlg, SECOidTag bulkalgtag, void *pwfn_arg);

/************************************************************************
 * cmsreclist.c - recipient list stuff
 ************************************************************************/
extern NSSCMSRecipient **nss_cms_recipient_list_create(NSSCMSRecipientInfo **recipientinfos);
extern void nss_cms_recipient_list_destroy(NSSCMSRecipient **recipient_list);
extern NSSCMSRecipientEncryptedKey *NSS_CMSRecipientEncryptedKey_Create(PLArenaPool *poolp);

/************************************************************************
 * cmsarray.c - misc array functions
 ************************************************************************/
/*
 * NSS_CMSArray_Alloc - allocate an array in an arena
 */
extern void **
NSS_CMSArray_Alloc(PLArenaPool *poolp, int n);

/*
 * NSS_CMSArray_Add - add an element to the end of an array
 */
extern SECStatus
NSS_CMSArray_Add(PLArenaPool *poolp, void ***array, void *obj);

/*
 * NSS_CMSArray_IsEmpty - check if array is empty
 */
extern PRBool
NSS_CMSArray_IsEmpty(void **array);

/*
 * NSS_CMSArray_Count - count number of elements in array
 */
extern int
NSS_CMSArray_Count(void **array);

/*
 * NSS_CMSArray_Sort - sort an array ascending, in place
 *
 * If "secondary" is not NULL, the same reordering gets applied to it.
 * If "tertiary" is not NULL, the same reordering gets applied to it.
 * "compare" is a function that returns
 *  < 0 when the first element is less than the second
 *  = 0 when the first element is equal to the second
 *  > 0 when the first element is greater than the second
 */
extern void
NSS_CMSArray_Sort(void **primary, int (*compare)(void *, void *), void **secondary, void **tertiary);

/************************************************************************
 * cmsattr.c - misc attribute functions
 ************************************************************************/
/*
 * NSS_CMSAttribute_Create - create an attribute
 *
 * if value is NULL, the attribute won't have a value. It can be added later
 * with NSS_CMSAttribute_AddValue.
 */
extern NSSCMSAttribute *
NSS_CMSAttribute_Create(PLArenaPool *poolp, SECOidTag oidtag, SECItem *value, PRBool encoded);

/*
 * NSS_CMSAttribute_AddValue - add another value to an attribute
 */
extern SECStatus
NSS_CMSAttribute_AddValue(PLArenaPool *poolp, NSSCMSAttribute *attr, SECItem *value);

/*
 * NSS_CMSAttribute_GetType - return the OID tag
 */
extern SECOidTag
NSS_CMSAttribute_GetType(NSSCMSAttribute *attr);

/*
 * NSS_CMSAttribute_GetValue - return the first attribute value
 *
 * We do some sanity checking first:
 * - Multiple values are *not* expected.
 * - Empty values are *not* expected.
 */
extern SECItem *
NSS_CMSAttribute_GetValue(NSSCMSAttribute *attr);

/*
 * NSS_CMSAttribute_CompareValue - compare the attribute's first value against data
 */
extern PRBool
NSS_CMSAttribute_CompareValue(NSSCMSAttribute *attr, SECItem *av);

/*
 * NSS_CMSAttributeArray_Encode - encode an Attribute array as SET OF Attributes
 *
 * If you are wondering why this routine does not reorder the attributes
 * first, and might be tempted to make it do so, see the comment by the
 * call to ReorderAttributes in cmsencode.c.  (Or, see who else calls this
 * and think long and hard about the implications of making it always
 * do the reordering.)
 */
extern SECItem *
NSS_CMSAttributeArray_Encode(PLArenaPool *poolp, NSSCMSAttribute ***attrs, SECItem *dest);

/*
 * NSS_CMSAttributeArray_Reorder - sort attribute array by attribute's DER encoding
 *
 * make sure that the order of the attributes guarantees valid DER (which must be
 * in lexigraphically ascending order for a SET OF); if reordering is necessary it
 * will be done in place (in attrs).
 */
extern SECStatus
NSS_CMSAttributeArray_Reorder(NSSCMSAttribute **attrs);

/*
 * NSS_CMSAttributeArray_FindAttrByOidTag - look through a set of attributes and
 * find one that matches the specified object ID.
 *
 * If "only" is true, then make sure that there is not more than one attribute
 * of the same type.  Otherwise, just return the first one found. (XXX Does
 * anybody really want that first-found behavior?  It was like that when I found it...)
 */
extern NSSCMSAttribute *
NSS_CMSAttributeArray_FindAttrByOidTag(NSSCMSAttribute **attrs, SECOidTag oidtag, PRBool only);

/*
 * NSS_CMSAttributeArray_AddAttr - add an attribute to an
 * array of attributes.
 */
extern SECStatus
NSS_CMSAttributeArray_AddAttr(PLArenaPool *poolp, NSSCMSAttribute ***attrs, NSSCMSAttribute *attr);

/*
 * NSS_CMSAttributeArray_SetAttr - set an attribute's value in a set of attributes
 */
extern SECStatus
NSS_CMSAttributeArray_SetAttr(PLArenaPool *poolp, NSSCMSAttribute ***attrs,
                              SECOidTag type, SECItem *value, PRBool encoded);

/*
 * NSS_CMSSignedData_AddTempCertificate - add temporary certificate references.
 * They may be needed for signature verification on the data, for example.
 */
extern SECStatus
NSS_CMSSignedData_AddTempCertificate(NSSCMSSignedData *sigd, CERTCertificate *cert);

/*
 * local function to handle compatibility issues
 * by mapping a signature algorithm back to a digest.
 */
SECOidTag NSS_CMSUtil_MapSignAlgs(SECOidTag signAlg);

/************************************************************************/

/*
 * local functions to handle user defined S/MIME content types
 */

PRBool NSS_CMSType_IsWrapper(SECOidTag type);
PRBool NSS_CMSType_IsData(SECOidTag type);
size_t NSS_CMSType_GetContentSize(SECOidTag type);
const SEC_ASN1Template *NSS_CMSType_GetTemplate(SECOidTag type);

void NSS_CMSGenericWrapperData_Destroy(SECOidTag type,
                                       NSSCMSGenericWrapperData *gd);
SECStatus NSS_CMSGenericWrapperData_Decode_BeforeData(SECOidTag type,
                                                      NSSCMSGenericWrapperData *gd);
SECStatus NSS_CMSGenericWrapperData_Decode_AfterData(SECOidTag type,
                                                     NSSCMSGenericWrapperData *gd);
SECStatus NSS_CMSGenericWrapperData_Decode_AfterEnd(SECOidTag type,
                                                    NSSCMSGenericWrapperData *gd);
SECStatus NSS_CMSGenericWrapperData_Encode_BeforeStart(SECOidTag type,
                                                       NSSCMSGenericWrapperData *gd);
SECStatus NSS_CMSGenericWrapperData_Encode_BeforeData(SECOidTag type,
                                                      NSSCMSGenericWrapperData *gd);
SECStatus NSS_CMSGenericWrapperData_Encode_AfterData(SECOidTag type,
                                                     NSSCMSGenericWrapperData *gd);

SEC_END_PROTOS

#endif /* _CMSLOCAL_H_ */
