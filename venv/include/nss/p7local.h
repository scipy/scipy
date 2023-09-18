/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Support routines for PKCS7 implementation, none of which are exported.
 * This file should only contain things that are needed by both the
 * encoding/creation side *and* the decoding/decryption side.  Anything
 * else should just be static routines in the appropriate file.
 *
 * Do not export this file!  If something in here is really needed outside
 * of pkcs7 code, first try to add a PKCS7 interface which will do it for
 * you.  If that has a problem, then just move out what you need, changing
 * its name as appropriate!
 */

#ifndef _P7LOCAL_H_
#define _P7LOCAL_H_

#include "secpkcs7.h"
#include "secasn1t.h"

extern const SEC_ASN1Template sec_PKCS7ContentInfoTemplate[];

/* opaque objects */
typedef struct sec_pkcs7_cipher_object sec_PKCS7CipherObject;

/************************************************************************/
SEC_BEGIN_PROTOS

/*
 * Look through a set of attributes and find one that matches the
 * specified object ID.  If "only" is true, then make sure that
 * there is not more than one attribute of the same type.  Otherwise,
 * just return the first one found. (XXX Does anybody really want
 * that first-found behavior?  It was like that when I found it...)
 */
extern SEC_PKCS7Attribute *sec_PKCS7FindAttribute(SEC_PKCS7Attribute **attrs,
                                                  SECOidTag oidtag,
                                                  PRBool only);
/*
 * Return the single attribute value, doing some sanity checking first:
 * - Multiple values are *not* expected.
 * - Empty values are *not* expected.
 */
extern SECItem *sec_PKCS7AttributeValue(SEC_PKCS7Attribute *attr);

/*
 * Encode a set of attributes (found in "src").
 */
extern SECItem *sec_PKCS7EncodeAttributes(PLArenaPool *poolp,
                                          SECItem *dest, void *src);

/*
 * Make sure that the order of the attributes guarantees valid DER
 * (which must be in lexigraphically ascending order for a SET OF);
 * if reordering is necessary it will be done in place (in attrs).
 */
extern SECStatus sec_PKCS7ReorderAttributes(SEC_PKCS7Attribute **attrs);

/*
 * Create a context for decrypting, based on the given key and algorithm.
 */
extern sec_PKCS7CipherObject *
sec_PKCS7CreateDecryptObject(PK11SymKey *key, SECAlgorithmID *algid);

/*
 * Create a context for encrypting, based on the given key and algorithm,
 * and fill in the algorithm id.
 */
extern sec_PKCS7CipherObject *
sec_PKCS7CreateEncryptObject(PLArenaPool *poolp, PK11SymKey *key,
                             SECOidTag algtag, SECAlgorithmID *algid);

/*
 * Destroy the given decryption or encryption object.
 */
extern void sec_PKCS7DestroyDecryptObject(sec_PKCS7CipherObject *obj);
extern void sec_PKCS7DestroyEncryptObject(sec_PKCS7CipherObject *obj);

/*
 * What will be the output length of the next call to encrypt/decrypt?
 * Result can be used to perform memory allocations.  Note that the amount
 * is exactly accurate only when not doing a block cipher or when final
 * is false, otherwise it is an upper bound on the amount because until
 * we see the data we do not know how many padding bytes there are
 * (always between 1 and the cipher block size).
 *
 * Note that this can return zero, which does not mean that the cipher
 * operation can be skipped!  (It simply means that there are not enough
 * bytes to make up an entire block; the bytes will be reserved until
 * there are enough to encrypt/decrypt at least one block.)  However,
 * if zero is returned it *does* mean that no output buffer need be
 * passed in to the subsequent cipher operation, as no output bytes
 * will be stored.
 */
extern unsigned int sec_PKCS7DecryptLength(sec_PKCS7CipherObject *obj,
                                           unsigned int input_len,
                                           PRBool final);
extern unsigned int sec_PKCS7EncryptLength(sec_PKCS7CipherObject *obj,
                                           unsigned int input_len,
                                           PRBool final);

/*
 * Decrypt a given length of input buffer (starting at "input" and
 * containing "input_len" bytes), placing the decrypted bytes in
 * "output" and storing the output length in "*output_len_p".
 * "obj" is the return value from sec_PKCS7CreateDecryptObject.
 * When "final" is true, this is the last of the data to be decrypted.
 */
extern SECStatus sec_PKCS7Decrypt(sec_PKCS7CipherObject *obj,
                                  unsigned char *output,
                                  unsigned int *output_len_p,
                                  unsigned int max_output_len,
                                  const unsigned char *input,
                                  unsigned int input_len,
                                  PRBool final);

/*
 * Encrypt a given length of input buffer (starting at "input" and
 * containing "input_len" bytes), placing the encrypted bytes in
 * "output" and storing the output length in "*output_len_p".
 * "obj" is the return value from sec_PKCS7CreateEncryptObject.
 * When "final" is true, this is the last of the data to be encrypted.
 */
extern SECStatus sec_PKCS7Encrypt(sec_PKCS7CipherObject *obj,
                                  unsigned char *output,
                                  unsigned int *output_len_p,
                                  unsigned int max_output_len,
                                  const unsigned char *input,
                                  unsigned int input_len,
                                  PRBool final);

/************************************************************************/
SEC_END_PROTOS

#endif /* _P7LOCAL_H_ */
