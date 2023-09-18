/*
 * blapi.h - public prototypes for the freebl library
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _BLAPI_H_
#define _BLAPI_H_

#include "blapit.h"
#include "hasht.h"
#include "cmac.h"
#include "alghmac.h"

SEC_BEGIN_PROTOS

/*
** RSA encryption/decryption. When encrypting/decrypting the output
** buffer must be at least the size of the public key modulus.
*/

extern SECStatus BL_Init(void);

/*
** Generate and return a new RSA public and private key.
**  Both keys are encoded in a single RSAPrivateKey structure.
**  "cx" is the random number generator context
**  "keySizeInBits" is the size of the key to be generated, in bits.
**     512, 1024, etc.
**  "publicExponent" when not NULL is a pointer to some data that
**     represents the public exponent to use. The data is a byte
**     encoded integer, in "big endian" order.
*/
extern RSAPrivateKey *RSA_NewKey(int keySizeInBits,
                                 SECItem *publicExponent);

/*
** Perform a raw public-key operation
**  Length of input and output buffers are equal to key's modulus len.
*/
extern SECStatus RSA_PublicKeyOp(RSAPublicKey *key,
                                 unsigned char *output,
                                 const unsigned char *input);

/*
** Perform a raw private-key operation
**  Length of input and output buffers are equal to key's modulus len.
*/
extern SECStatus RSA_PrivateKeyOp(RSAPrivateKey *key,
                                  unsigned char *output,
                                  const unsigned char *input);

/*
** Perform a raw private-key operation, and check the parameters used in
** the operation for validity by performing a test operation first.
**  Length of input and output buffers are equal to key's modulus len.
*/
extern SECStatus RSA_PrivateKeyOpDoubleChecked(RSAPrivateKey *key,
                                               unsigned char *output,
                                               const unsigned char *input);

/*
** Perform a check of private key parameters for consistency.
*/
extern SECStatus RSA_PrivateKeyCheck(const RSAPrivateKey *key);

/*
** Given only minimal private key parameters, fill in the rest of the
** parameters.
**
**
** All the entries, including those supplied by the caller, will be
** overwritten with data alocated out of the arena.
**
** If no arena is supplied, one will be created.
**
** The following fields must be supplied in order for this function
** to succeed:
**   one of either publicExponent or privateExponent
**   two more of the following 5 parameters (not counting the above).
**      modulus (n)
**      prime1  (p)
**      prime2  (q)
**      publicExponent (e)
**      privateExponent (d)
**
** NOTE: if only the publicExponent, privateExponent, and one prime is given,
** then there may be more than one RSA key that matches that combination. If
** we find 2 possible valid keys that meet this criteria, we return an error.
** If we return the wrong key, and the original modulus is compared to the
** new modulus, both can be factored by calculateing gcd(n_old,n_new) to get
** the common prime.
**
** NOTE: in some cases the publicExponent must be less than 2^23 for this
** function to work correctly. (The case where we have only one of: modulus
** prime1 and prime2).
**
** All parameters will be replaced in the key structure with new parameters
** allocated out of the arena. There is no attempt to free the old structures.
** prime1 will always be greater than prime2 (even if the caller supplies the
** smaller prime as prime1 or the larger prime as prime2). The parameters are
** not overwritten on failure.
**
** While the remaining Chinese remainder theorem parameters (dp,dp, and qinv)
** can also be used in reconstructing the private key, they are currently
** ignored in this implementation.
*/
extern SECStatus RSA_PopulatePrivateKey(RSAPrivateKey *key);

/********************************************************************
** RSA algorithm
*/

/********************************************************************
** Raw signing/encryption/decryption operations.
**
** No padding or formatting will be applied.
** inputLen MUST be equivalent to the modulus size (in bytes).
*/
extern SECStatus
RSA_SignRaw(RSAPrivateKey *key,
            unsigned char *output,
            unsigned int *outputLen,
            unsigned int maxOutputLen,
            const unsigned char *input,
            unsigned int inputLen);

extern SECStatus
RSA_CheckSignRaw(RSAPublicKey *key,
                 const unsigned char *sig,
                 unsigned int sigLen,
                 const unsigned char *hash,
                 unsigned int hashLen);

extern SECStatus
RSA_CheckSignRecoverRaw(RSAPublicKey *key,
                        unsigned char *data,
                        unsigned int *dataLen,
                        unsigned int maxDataLen,
                        const unsigned char *sig,
                        unsigned int sigLen);

extern SECStatus
RSA_EncryptRaw(RSAPublicKey *key,
               unsigned char *output,
               unsigned int *outputLen,
               unsigned int maxOutputLen,
               const unsigned char *input,
               unsigned int inputLen);

extern SECStatus
RSA_DecryptRaw(RSAPrivateKey *key,
               unsigned char *output,
               unsigned int *outputLen,
               unsigned int maxOutputLen,
               const unsigned char *input,
               unsigned int inputLen);

/********************************************************************
** RSAES-OAEP encryption/decryption, as defined in RFC 3447, Section 7.1.
**
** Note: Only MGF1 is supported as the mask generation function. It will be
** used with maskHashAlg as the inner hash function.
**
** Unless performing Known Answer Tests, "seed" should be NULL, indicating that
** freebl should generate a random value. Otherwise, it should be an octet
** string of seedLen bytes, which should be the same size as the output of
** hashAlg.
*/
extern SECStatus
RSA_EncryptOAEP(RSAPublicKey *key,
                HASH_HashType hashAlg,
                HASH_HashType maskHashAlg,
                const unsigned char *label,
                unsigned int labelLen,
                const unsigned char *seed,
                unsigned int seedLen,
                unsigned char *output,
                unsigned int *outputLen,
                unsigned int maxOutputLen,
                const unsigned char *input,
                unsigned int inputLen);

extern SECStatus
RSA_DecryptOAEP(RSAPrivateKey *key,
                HASH_HashType hashAlg,
                HASH_HashType maskHashAlg,
                const unsigned char *label,
                unsigned int labelLen,
                unsigned char *output,
                unsigned int *outputLen,
                unsigned int maxOutputLen,
                const unsigned char *input,
                unsigned int inputLen);

/********************************************************************
** RSAES-PKCS1-v1_5 encryption/decryption, as defined in RFC 3447, Section 7.2.
*/
extern SECStatus
RSA_EncryptBlock(RSAPublicKey *key,
                 unsigned char *output,
                 unsigned int *outputLen,
                 unsigned int maxOutputLen,
                 const unsigned char *input,
                 unsigned int inputLen);

extern SECStatus
RSA_DecryptBlock(RSAPrivateKey *key,
                 unsigned char *output,
                 unsigned int *outputLen,
                 unsigned int maxOutputLen,
                 const unsigned char *input,
                 unsigned int inputLen);

/********************************************************************
** RSASSA-PSS signing/verifying, as defined in RFC 3447, Section 8.1.
**
** Note: Only MGF1 is supported as the mask generation function. It will be
** used with maskHashAlg as the inner hash function.
**
** Unless performing Known Answer Tests, "salt" should be NULL, indicating that
** freebl should generate a random value.
*/
extern SECStatus
RSA_SignPSS(RSAPrivateKey *key,
            HASH_HashType hashAlg,
            HASH_HashType maskHashAlg,
            const unsigned char *salt,
            unsigned int saltLen,
            unsigned char *output,
            unsigned int *outputLen,
            unsigned int maxOutputLen,
            const unsigned char *input,
            unsigned int inputLen);

extern SECStatus
RSA_CheckSignPSS(RSAPublicKey *key,
                 HASH_HashType hashAlg,
                 HASH_HashType maskHashAlg,
                 unsigned int saltLen,
                 const unsigned char *sig,
                 unsigned int sigLen,
                 const unsigned char *hash,
                 unsigned int hashLen);

/********************************************************************
** RSASSA-PKCS1-v1_5 signing/verifying, as defined in RFC 3447, Section 8.2.
**
** These functions expect as input to be the raw value to be signed. For most
** cases using PKCS1-v1_5, this should be the value of T, the DER-encoded
** DigestInfo structure defined in Section 9.2, Step 2.
** Note: This can also be used for signatures that use PKCS1-v1_5 padding, such
** as the signatures used in SSL/TLS, which sign a raw hash.
*/
extern SECStatus
RSA_Sign(RSAPrivateKey *key,
         unsigned char *output,
         unsigned int *outputLen,
         unsigned int maxOutputLen,
         const unsigned char *data,
         unsigned int dataLen);

extern SECStatus
RSA_CheckSign(RSAPublicKey *key,
              const unsigned char *sig,
              unsigned int sigLen,
              const unsigned char *data,
              unsigned int dataLen);

extern SECStatus
RSA_CheckSignRecover(RSAPublicKey *key,
                     unsigned char *output,
                     unsigned int *outputLen,
                     unsigned int maxOutputLen,
                     const unsigned char *sig,
                     unsigned int sigLen);

/********************************************************************
** DSA signing algorithm
*/

/* Generate a new random value within the interval [2, q-1].
*/
extern SECStatus DSA_NewRandom(PLArenaPool *arena, const SECItem *q,
                               SECItem *random);

/*
** Generate and return a new DSA public and private key pair,
**  both of which are encoded into a single DSAPrivateKey struct.
**  "params" is a pointer to the PQG parameters for the domain
**  Uses a random seed.
*/
extern SECStatus DSA_NewKey(const PQGParams *params,
                            DSAPrivateKey **privKey);

/* signature is caller-supplied buffer of at least 20 bytes.
** On input,  signature->len == size of buffer to hold signature.
**            digest->len    == size of digest.
** On output, signature->len == size of signature in buffer.
** Uses a random seed.
*/
extern SECStatus DSA_SignDigest(DSAPrivateKey *key,
                                SECItem *signature,
                                const SECItem *digest);

/* signature is caller-supplied buffer of at least 20 bytes.
** On input,  signature->len == size of buffer to hold signature.
**            digest->len    == size of digest.
*/
extern SECStatus DSA_VerifyDigest(DSAPublicKey *key,
                                  const SECItem *signature,
                                  const SECItem *digest);

/* For FIPS compliance testing. Seed must be exactly 20 bytes long */
extern SECStatus DSA_NewKeyFromSeed(const PQGParams *params,
                                    const unsigned char *seed,
                                    DSAPrivateKey **privKey);

/* For FIPS compliance testing. Seed must be exactly 20 bytes. */
extern SECStatus DSA_SignDigestWithSeed(DSAPrivateKey *key,
                                        SECItem *signature,
                                        const SECItem *digest,
                                        const unsigned char *seed);

/******************************************************
** Diffie Helman key exchange algorithm
*/

/* Generates parameters for Diffie-Helman key generation.
**  primeLen is the length in bytes of prime P to be generated.
*/
extern SECStatus DH_GenParam(int primeLen, DHParams **params);

/* Generates a public and private key, both of which are encoded in a single
**  DHPrivateKey struct. Params is input, privKey are output.
**  This is Phase 1 of Diffie Hellman.
*/
extern SECStatus DH_NewKey(DHParams *params,
                           DHPrivateKey **privKey);

/*
** DH_Derive does the Diffie-Hellman phase 2 calculation, using the
** other party's publicValue, and the prime and our privateValue.
** maxOutBytes is the requested length of the generated secret in bytes.
** A zero value means produce a value of any length up to the size of
** the prime.   If successful, derivedSecret->data is set
** to the address of the newly allocated buffer containing the derived
** secret, and derivedSecret->len is the size of the secret produced.
** The size of the secret produced will depend on the value of outBytes.
** If outBytes is 0, the key length will be all the significant bytes of
** the derived secret (leading zeros are dropped). This length could be less
** than the length of the prime. If outBytes is nonzero, the length of the
** produced key will be outBytes long. If the key is truncated, the most
** significant bytes are truncated. If it is expanded, zero bytes are added
** at the beginning.
** It is the caller's responsibility to free the allocated buffer
** containing the derived secret.
*/
extern SECStatus DH_Derive(SECItem *publicValue,
                           SECItem *prime,
                           SECItem *privateValue,
                           SECItem *derivedSecret,
                           unsigned int outBytes);

/*
** KEA_CalcKey returns octet string with the private key for a dual
** Diffie-Helman  key generation as specified for government key exchange.
*/
extern SECStatus KEA_Derive(SECItem *prime,
                            SECItem *public1,
                            SECItem *public2,
                            SECItem *private1,
                            SECItem *private2,
                            SECItem *derivedSecret);

/*
 * verify that a KEA or DSA public key is a valid key for this prime and
 * subprime domain.
 */
extern PRBool KEA_Verify(SECItem *Y, SECItem *prime, SECItem *subPrime);

/* verify a value is prime */
PRBool KEA_PrimeCheck(SECItem *prime);

/****************************************
 * J-PAKE key transport
 */

/* Given gx == g^x, create a Schnorr zero-knowledge proof for the value x
 * using the specified hash algorithm and signer ID. The signature is
 * returned in the values gv and r. testRandom must be NULL for a PRNG
 * generated random committment to be used in the sigature. When testRandom
 * is non-NULL, that value must contain a value in the subgroup q; that
 * value will be used instead of a PRNG-generated committment in order to
 * facilitate known-answer tests.
 *
 * If gxIn is non-NULL then it must contain a pre-computed value of g^x that
 * will be used by the function; in this case, the gxOut parameter must be NULL.
 * If the gxIn parameter is NULL then gxOut must be non-NULL; in this case
 * gxOut will contain the value g^x on output.
 *
 * gx (if not supplied by the caller), gv, and r will be allocated in the arena.
 * The arena is *not* optional so do not pass NULL for the arena parameter.
 * The arena should be zeroed when it is freed.
 */
SECStatus
JPAKE_Sign(PLArenaPool *arena, const PQGParams *pqg, HASH_HashType hashType,
           const SECItem *signerID, const SECItem *x,
           const SECItem *testRandom, const SECItem *gxIn, SECItem *gxOut,
           SECItem *gv, SECItem *r);

/* Given gx == g^x, verify the Schnorr zero-knowledge proof (gv, r) for the
 * value x using the specified hash algorithm and signer ID.
 *
 * The arena is *not* optional so do not pass NULL for the arena parameter.
 */
SECStatus
JPAKE_Verify(PLArenaPool *arena, const PQGParams *pqg,
             HASH_HashType hashType, const SECItem *signerID,
             const SECItem *peerID, const SECItem *gx,
             const SECItem *gv, const SECItem *r);

/* Call before round 2 with x2, s, and x2s all non-NULL. This will calculate
 * base = g^(x1+x3+x4) (mod p) and x2s = x2*s (mod q). The values to send in
 * round 2 (A and the proof of knowledge of x2s) can then be calculated with
 * JPAKE_Sign using pqg->base = base and x = x2s.
 *
 * Call after round 2 with x2, s, and x2s all NULL, and passing (gx1, gx2, gx3)
 * instead of (gx1, gx3, gx4). This will calculate base = g^(x1+x2+x3). Then call
 * JPAKE_Verify with pqg->base = base and then JPAKE_Final.
 *
 * base and x2s will be allocated in the arena. The arena is *not* optional so
 * do not pass NULL for the arena parameter. The arena should be zeroed when it
 * is freed.
*/
SECStatus
JPAKE_Round2(PLArenaPool *arena, const SECItem *p, const SECItem *q,
             const SECItem *gx1, const SECItem *gx3, const SECItem *gx4,
             SECItem *base, const SECItem *x2, const SECItem *s, SECItem *x2s);

/* K = (B/g^(x2*x4*s))^x2 (mod p)
 *
 * K will be allocated in the arena. The arena is *not* optional so do not pass
 * NULL for the arena parameter. The arena should be zeroed when it is freed.
 */
SECStatus
JPAKE_Final(PLArenaPool *arena, const SECItem *p, const SECItem *q,
            const SECItem *x2, const SECItem *gx4, const SECItem *x2s,
            const SECItem *B, SECItem *K);

/******************************************************
** Elliptic Curve algorithms
*/

/* Generates a public and private key, both of which are encoded
** in a single ECPrivateKey struct. Params is input, privKey are
** output.
*/
extern SECStatus EC_NewKey(ECParams *params,
                           ECPrivateKey **privKey);

extern SECStatus EC_NewKeyFromSeed(ECParams *params,
                                   ECPrivateKey **privKey,
                                   const unsigned char *seed,
                                   int seedlen);

/* Validates an EC public key as described in Section 5.2.2 of
 * X9.62. Such validation prevents against small subgroup attacks
 * when the ECDH primitive is used with the cofactor.
 */
extern SECStatus EC_ValidatePublicKey(ECParams *params,
                                      SECItem *publicValue);

/*
** ECDH_Derive performs a scalar point multiplication of a point
** representing a (peer's) public key and a large integer representing
** a private key (its own). Both keys must use the same elliptic curve
** parameters. If the withCofactor parameter is true, the
** multiplication also uses the cofactor associated with the curve
** parameters.  The output of this scheme is the x-coordinate of the
** resulting point. If successful, derivedSecret->data is set to the
** address of the newly allocated buffer containing the derived
** secret, and derivedSecret->len is the size of the secret
** produced. It is the caller's responsibility to free the allocated
** buffer containing the derived secret.
*/
extern SECStatus ECDH_Derive(SECItem *publicValue,
                             ECParams *params,
                             SECItem *privateValue,
                             PRBool withCofactor,
                             SECItem *derivedSecret);

/* On input,  signature->len == size of buffer to hold signature.
**            digest->len    == size of digest.
** On output, signature->len == size of signature in buffer.
** Uses a random seed.
*/
extern SECStatus ECDSA_SignDigest(ECPrivateKey *key,
                                  SECItem *signature,
                                  const SECItem *digest);

/* On input,  signature->len == size of buffer to hold signature.
**            digest->len    == size of digest.
*/
extern SECStatus ECDSA_VerifyDigest(ECPublicKey *key,
                                    const SECItem *signature,
                                    const SECItem *digest);

/* Uses the provided seed. */
extern SECStatus ECDSA_SignDigestWithSeed(ECPrivateKey *key,
                                          SECItem *signature,
                                          const SECItem *digest,
                                          const unsigned char *seed,
                                          const int seedlen);

/******************************************/
/*
** RC4 symmetric stream cypher
*/

/*
** Create a new RC4 context suitable for RC4 encryption/decryption.
**  "key" raw key data
**  "len" the number of bytes of key data
*/
extern RC4Context *RC4_CreateContext(const unsigned char *key, int len);

extern RC4Context *RC4_AllocateContext(void);
extern SECStatus RC4_InitContext(RC4Context *cx,
                                 const unsigned char *key,
                                 unsigned int keylen,
                                 const unsigned char *,
                                 int,
                                 unsigned int,
                                 unsigned int);

/*
** Destroy an RC4 encryption/decryption context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void RC4_DestroyContext(RC4Context *cx, PRBool freeit);

/*
** Perform RC4 encryption.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus RC4_Encrypt(RC4Context *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/*
** Perform RC4 decryption.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus RC4_Decrypt(RC4Context *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** RC2 symmetric block cypher
*/

/*
** Create a new RC2 context suitable for RC2 encryption/decryption.
**  "key" raw key data
**  "len" the number of bytes of key data
**  "iv" is the CBC initialization vector (if mode is NSS_RC2_CBC)
**  "mode" one of NSS_RC2 or NSS_RC2_CBC
**  "effectiveKeyLen" is the effective key length (as specified in
**      RFC 2268) in bytes (not bits).
**
** When mode is set to NSS_RC2_CBC the RC2 cipher is run in "cipher block
** chaining" mode.
*/
extern RC2Context *RC2_CreateContext(const unsigned char *key, unsigned int len,
                                     const unsigned char *iv, int mode,
                                     unsigned effectiveKeyLen);
extern RC2Context *RC2_AllocateContext(void);
extern SECStatus RC2_InitContext(RC2Context *cx,
                                 const unsigned char *key,
                                 unsigned int keylen,
                                 const unsigned char *iv,
                                 int mode,
                                 unsigned int effectiveKeyLen,
                                 unsigned int);

/*
** Destroy an RC2 encryption/decryption context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void RC2_DestroyContext(RC2Context *cx, PRBool freeit);

/*
** Perform RC2 encryption.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus RC2_Encrypt(RC2Context *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/*
** Perform RC2 decryption.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus RC2_Decrypt(RC2Context *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** RC5 symmetric block cypher -- 64-bit block size
*/

/*
** Create a new RC5 context suitable for RC5 encryption/decryption.
**      "key" raw key data
**      "len" the number of bytes of key data
**      "iv" is the CBC initialization vector (if mode is NSS_RC5_CBC)
**      "mode" one of NSS_RC5 or NSS_RC5_CBC
**
** When mode is set to NSS_RC5_CBC the RC5 cipher is run in "cipher block
** chaining" mode.
*/
extern RC5Context *RC5_CreateContext(const SECItem *key, unsigned int rounds,
                                     unsigned int wordSize, const unsigned char *iv, int mode);
extern RC5Context *RC5_AllocateContext(void);
extern SECStatus RC5_InitContext(RC5Context *cx,
                                 const unsigned char *key,
                                 unsigned int keylen,
                                 const unsigned char *iv,
                                 int mode,
                                 unsigned int rounds,
                                 unsigned int wordSize);

/*
** Destroy an RC5 encryption/decryption context.
**      "cx" the context
**      "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void RC5_DestroyContext(RC5Context *cx, PRBool freeit);

/*
** Perform RC5 encryption.
**      "cx" the context
**      "output" the output buffer to store the encrypted data.
**      "outputLen" how much data is stored in "output". Set by the routine
**         after some data is stored in output.
**      "maxOutputLen" the maximum amount of data that can ever be
**         stored in "output"
**      "input" the input data
**      "inputLen" the amount of input data
*/
extern SECStatus RC5_Encrypt(RC5Context *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/*
** Perform RC5 decryption.
**      "cx" the context
**      "output" the output buffer to store the decrypted data.
**      "outputLen" how much data is stored in "output". Set by the routine
**         after some data is stored in output.
**      "maxOutputLen" the maximum amount of data that can ever be
**         stored in "output"
**      "input" the input data
**      "inputLen" the amount of input data
*/

extern SECStatus RC5_Decrypt(RC5Context *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** DES symmetric block cypher
*/

/*
** Create a new DES context suitable for DES encryption/decryption.
**  "key" raw key data
**  "len" the number of bytes of key data
**  "iv" is the CBC initialization vector (if mode is NSS_DES_CBC or
**     mode is DES_EDE3_CBC)
**  "mode" one of NSS_DES, NSS_DES_CBC, NSS_DES_EDE3 or NSS_DES_EDE3_CBC
**  "encrypt" is PR_TRUE if the context will be used for encryption
**
** When mode is set to NSS_DES_CBC or NSS_DES_EDE3_CBC then the DES
** cipher is run in "cipher block chaining" mode.
*/
extern DESContext *DES_CreateContext(const unsigned char *key,
                                     const unsigned char *iv,
                                     int mode, PRBool encrypt);
extern DESContext *DES_AllocateContext(void);
extern SECStatus DES_InitContext(DESContext *cx,
                                 const unsigned char *key,
                                 unsigned int keylen,
                                 const unsigned char *iv,
                                 int mode,
                                 unsigned int encrypt,
                                 unsigned int);

/*
** Destroy an DES encryption/decryption context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void DES_DestroyContext(DESContext *cx, PRBool freeit);

/*
** Perform DES encryption.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
**
** NOTE: the inputLen must be a multiple of DES_KEY_LENGTH
*/
extern SECStatus DES_Encrypt(DESContext *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/*
** Perform DES decryption.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
**
** NOTE: the inputLen must be a multiple of DES_KEY_LENGTH
*/
extern SECStatus DES_Decrypt(DESContext *cx, unsigned char *output,
                             unsigned int *outputLen, unsigned int maxOutputLen,
                             const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** SEED symmetric block cypher
*/
extern SEEDContext *
SEED_CreateContext(const unsigned char *key, const unsigned char *iv,
                   int mode, PRBool encrypt);
extern SEEDContext *SEED_AllocateContext(void);
extern SECStatus SEED_InitContext(SEEDContext *cx,
                                  const unsigned char *key,
                                  unsigned int keylen,
                                  const unsigned char *iv,
                                  int mode, unsigned int encrypt,
                                  unsigned int);
extern void SEED_DestroyContext(SEEDContext *cx, PRBool freeit);
extern SECStatus
SEED_Encrypt(SEEDContext *cx, unsigned char *output,
             unsigned int *outputLen, unsigned int maxOutputLen,
             const unsigned char *input, unsigned int inputLen);
extern SECStatus
SEED_Decrypt(SEEDContext *cx, unsigned char *output,
             unsigned int *outputLen, unsigned int maxOutputLen,
             const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** AES symmetric block cypher (Rijndael)
*/

/*
** Create a new AES context suitable for AES encryption/decryption.
**  "key" raw key data
**  "keylen" the number of bytes of key data (16, 24, or 32)
**  "blocklen" is the blocksize to use. NOTE: only 16 is supported!
*/
extern AESContext *
AES_CreateContext(const unsigned char *key, const unsigned char *iv,
                  int mode, int encrypt,
                  unsigned int keylen, unsigned int blocklen);
extern AESContext *AES_AllocateContext(void);
extern SECStatus AES_InitContext(AESContext *cx,
                                 const unsigned char *key,
                                 unsigned int keylen,
                                 const unsigned char *iv,
                                 int mode,
                                 unsigned int encrypt,
                                 unsigned int blocklen);

/*
** Destroy a AES encryption/decryption context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void
AES_DestroyContext(AESContext *cx, PRBool freeit);

/*
** Perform AES encryption.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
AES_Encrypt(AESContext *cx, unsigned char *output,
            unsigned int *outputLen, unsigned int maxOutputLen,
            const unsigned char *input, unsigned int inputLen);

/*
** Perform AES decryption.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
AES_Decrypt(AESContext *cx, unsigned char *output,
            unsigned int *outputLen, unsigned int maxOutputLen,
            const unsigned char *input, unsigned int inputLen);
/*
** Perform AES AEAD operation (either encrypt or decrypt), controlled by
** the context.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
**  "params" pointer to an AEAD specific param PKCS #11 param structure
**  "paramsLen" length of the param structure pointed to by params
**  "aad" addition authenticated data
**  "aadLen" the amount of additional authenticated data.
*/
extern SECStatus
AES_AEAD(AESContext *cx, unsigned char *output,
         unsigned int *outputLen, unsigned int maxOutputLen,
         const unsigned char *input, unsigned int inputLen,
         void *params, unsigned int paramsLen,
         const unsigned char *aad, unsigned int aadLen);

/******************************************/
/*
** AES key wrap algorithm, RFC 3394
*/

/*
** Create a new AES context suitable for AES encryption/decryption.
**  "key" raw key data
**      "iv"  The 8 byte "initial value"
**      "encrypt", a boolean, true for key wrapping, false for unwrapping.
**  "keylen" the number of bytes of key data (16, 24, or 32)
*/
extern AESKeyWrapContext *
AESKeyWrap_CreateContext(const unsigned char *key, const unsigned char *iv,
                         int encrypt, unsigned int keylen);
extern AESKeyWrapContext *AESKeyWrap_AllocateContext(void);
extern SECStatus
AESKeyWrap_InitContext(AESKeyWrapContext *cx,
                       const unsigned char *key,
                       unsigned int keylen,
                       const unsigned char *iv,
                       int,
                       unsigned int encrypt,
                       unsigned int);

/*
** Destroy a AES KeyWrap context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void
AESKeyWrap_DestroyContext(AESKeyWrapContext *cx, PRBool freeit);

/*
** Perform AES key wrap.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
AESKeyWrap_Encrypt(AESKeyWrapContext *cx, unsigned char *output,
                   unsigned int *outputLen, unsigned int maxOutputLen,
                   const unsigned char *input, unsigned int inputLen);

/*
** Perform AES key unwrap.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
AESKeyWrap_Decrypt(AESKeyWrapContext *cx, unsigned char *output,
                   unsigned int *outputLen, unsigned int maxOutputLen,
                   const unsigned char *input, unsigned int inputLen);

/*
** Perform AES padded key wrap.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
AESKeyWrap_EncryptKWP(AESKeyWrapContext *cx, unsigned char *output,
                      unsigned int *outputLen, unsigned int maxOutputLen,
                      const unsigned char *input, unsigned int inputLen);

/*
** Perform AES padded key unwrap.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
AESKeyWrap_DecryptKWP(AESKeyWrapContext *cx, unsigned char *output,
                      unsigned int *outputLen, unsigned int maxOutputLen,
                      const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** Camellia symmetric block cypher
*/

/*
** Create a new Camellia context suitable for Camellia encryption/decryption.
**  "key" raw key data
**  "keylen" the number of bytes of key data (16, 24, or 32)
*/
extern CamelliaContext *
Camellia_CreateContext(const unsigned char *key, const unsigned char *iv,
                       int mode, int encrypt, unsigned int keylen);

extern CamelliaContext *Camellia_AllocateContext(void);
extern SECStatus Camellia_InitContext(CamelliaContext *cx,
                                      const unsigned char *key,
                                      unsigned int keylen,
                                      const unsigned char *iv,
                                      int mode,
                                      unsigned int encrypt,
                                      unsigned int unused);
/*
** Destroy a Camellia encryption/decryption context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void
Camellia_DestroyContext(CamelliaContext *cx, PRBool freeit);

/*
** Perform Camellia encryption.
**  "cx" the context
**  "output" the output buffer to store the encrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
Camellia_Encrypt(CamelliaContext *cx, unsigned char *output,
                 unsigned int *outputLen, unsigned int maxOutputLen,
                 const unsigned char *input, unsigned int inputLen);

/*
** Perform Camellia decryption.
**  "cx" the context
**  "output" the output buffer to store the decrypted data.
**  "outputLen" how much data is stored in "output". Set by the routine
**     after some data is stored in output.
**  "maxOutputLen" the maximum amount of data that can ever be
**     stored in "output"
**  "input" the input data
**  "inputLen" the amount of input data
*/
extern SECStatus
Camellia_Decrypt(CamelliaContext *cx, unsigned char *output,
                 unsigned int *outputLen, unsigned int maxOutputLen,
                 const unsigned char *input, unsigned int inputLen);

/******************************************/
/*
** ChaCha20 block cipher
*/

extern SECStatus ChaCha20_InitContext(ChaCha20Context *ctx,
                                      const unsigned char *key,
                                      unsigned int keyLen,
                                      const unsigned char *nonce,
                                      unsigned int nonceLen,
                                      PRUint32 ctr);

extern ChaCha20Context *ChaCha20_CreateContext(const unsigned char *key,
                                               unsigned int keyLen,
                                               const unsigned char *nonce,
                                               unsigned int nonceLen,
                                               PRUint32 ctr);

extern void ChaCha20_DestroyContext(ChaCha20Context *ctx, PRBool freeit);

/******************************************/
/*
** ChaCha20+Poly1305 AEAD
*/

extern SECStatus ChaCha20Poly1305_InitContext(ChaCha20Poly1305Context *ctx,
                                              const unsigned char *key,
                                              unsigned int keyLen,
                                              unsigned int tagLen);

extern ChaCha20Poly1305Context *ChaCha20Poly1305_CreateContext(
    const unsigned char *key, unsigned int keyLen, unsigned int tagLen);

extern void ChaCha20Poly1305_DestroyContext(ChaCha20Poly1305Context *ctx,
                                            PRBool freeit);

extern SECStatus ChaCha20Poly1305_Seal(
    const ChaCha20Poly1305Context *ctx, unsigned char *output,
    unsigned int *outputLen, unsigned int maxOutputLen,
    const unsigned char *input, unsigned int inputLen,
    const unsigned char *nonce, unsigned int nonceLen,
    const unsigned char *ad, unsigned int adLen);

extern SECStatus ChaCha20Poly1305_Open(
    const ChaCha20Poly1305Context *ctx, unsigned char *output,
    unsigned int *outputLen, unsigned int maxOutputLen,
    const unsigned char *input, unsigned int inputLen,
    const unsigned char *nonce, unsigned int nonceLen,
    const unsigned char *ad, unsigned int adLen);

extern SECStatus ChaCha20Poly1305_Encrypt(
    const ChaCha20Poly1305Context *ctx, unsigned char *output,
    unsigned int *outputLen, unsigned int maxOutputLen,
    const unsigned char *input, unsigned int inputLen,
    const unsigned char *nonce, unsigned int nonceLen,
    const unsigned char *ad, unsigned int adLen, unsigned char *tagOut);

extern SECStatus ChaCha20Poly1305_Decrypt(
    const ChaCha20Poly1305Context *ctx, unsigned char *output,
    unsigned int *outputLen, unsigned int maxOutputLen,
    const unsigned char *input, unsigned int inputLen,
    const unsigned char *nonce, unsigned int nonceLen,
    const unsigned char *ad, unsigned int adLen, unsigned char *tagIn);

extern SECStatus ChaCha20_Xor(
    unsigned char *output, const unsigned char *block, unsigned int len,
    const unsigned char *k, const unsigned char *nonce, PRUint32 ctr);

/******************************************/
/*
** MD5 secure hash function
*/

/*
** Hash a null terminated string "src" into "dest" using MD5
*/
extern SECStatus MD5_Hash(unsigned char *dest, const char *src);

/*
** Hash a non-null terminated string "src" into "dest" using MD5
*/
extern SECStatus MD5_HashBuf(unsigned char *dest, const unsigned char *src,
                             PRUint32 src_length);

/*
** Create a new MD5 context
*/
extern MD5Context *MD5_NewContext(void);

/*
** Destroy an MD5 secure hash context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void MD5_DestroyContext(MD5Context *cx, PRBool freeit);

/*
** Reset an MD5 context, preparing it for a fresh round of hashing
*/
extern void MD5_Begin(MD5Context *cx);

/*
** Update the MD5 hash function with more data.
**  "cx" the context
**  "input" the data to hash
**  "inputLen" the amount of data to hash
*/
extern void MD5_Update(MD5Context *cx,
                       const unsigned char *input, unsigned int inputLen);

/*
** Finish the MD5 hash function. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 16 bytes of digest data are stored
**  "digestLen" where the digest length (16) is stored
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void MD5_End(MD5Context *cx, unsigned char *digest,
                    unsigned int *digestLen, unsigned int maxDigestLen);

/*
** Export the current state of the MD5 hash without appending the standard
** padding and length bytes. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 16 bytes of digest data are stored
**  "digestLen" where the digest length (16) is stored (optional)
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void MD5_EndRaw(MD5Context *cx, unsigned char *digest,
                       unsigned int *digestLen, unsigned int maxDigestLen);

/*
 * Return the the size of a buffer needed to flatten the MD5 Context into
 *    "cx" the context
 *  returns size;
 */
extern unsigned int MD5_FlattenSize(MD5Context *cx);

/*
 * Flatten the MD5 Context into a buffer:
 *    "cx" the context
 *    "space" the buffer to flatten to
 *  returns status;
 */
extern SECStatus MD5_Flatten(MD5Context *cx, unsigned char *space);

/*
 * Resurrect a flattened context into a MD5 Context
 *    "space" the buffer of the flattend buffer
 *    "arg" ptr to void used by cryptographic resurrect
 *  returns resurected context;
 */
extern MD5Context *MD5_Resurrect(unsigned char *space, void *arg);
extern void MD5_Clone(MD5Context *dest, MD5Context *src);

/*
** trace the intermediate state info of the MD5 hash.
*/
extern void MD5_TraceState(MD5Context *cx);

/******************************************/
/*
** MD2 secure hash function
*/

/*
** Hash a null terminated string "src" into "dest" using MD2
*/
extern SECStatus MD2_Hash(unsigned char *dest, const char *src);

/*
** Create a new MD2 context
*/
extern MD2Context *MD2_NewContext(void);

/*
** Destroy an MD2 secure hash context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void MD2_DestroyContext(MD2Context *cx, PRBool freeit);

/*
** Reset an MD2 context, preparing it for a fresh round of hashing
*/
extern void MD2_Begin(MD2Context *cx);

/*
** Update the MD2 hash function with more data.
**  "cx" the context
**  "input" the data to hash
**  "inputLen" the amount of data to hash
*/
extern void MD2_Update(MD2Context *cx,
                       const unsigned char *input, unsigned int inputLen);

/*
** Finish the MD2 hash function. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 16 bytes of digest data are stored
**  "digestLen" where the digest length (16) is stored
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void MD2_End(MD2Context *cx, unsigned char *digest,
                    unsigned int *digestLen, unsigned int maxDigestLen);

/*
 * Return the the size of a buffer needed to flatten the MD2 Context into
 *    "cx" the context
 *  returns size;
 */
extern unsigned int MD2_FlattenSize(MD2Context *cx);

/*
 * Flatten the MD2 Context into a buffer:
 *    "cx" the context
 *    "space" the buffer to flatten to
 *  returns status;
 */
extern SECStatus MD2_Flatten(MD2Context *cx, unsigned char *space);

/*
 * Resurrect a flattened context into a MD2 Context
 *    "space" the buffer of the flattend buffer
 *    "arg" ptr to void used by cryptographic resurrect
 *  returns resurected context;
 */
extern MD2Context *MD2_Resurrect(unsigned char *space, void *arg);
extern void MD2_Clone(MD2Context *dest, MD2Context *src);

/******************************************/
/*
** SHA-1 secure hash function
*/

/*
** Hash a null terminated string "src" into "dest" using SHA-1
*/
extern SECStatus SHA1_Hash(unsigned char *dest, const char *src);

/*
** Hash a non-null terminated string "src" into "dest" using SHA-1
*/
extern SECStatus SHA1_HashBuf(unsigned char *dest, const unsigned char *src,
                              PRUint32 src_length);

/*
** Create a new SHA-1 context
*/
extern SHA1Context *SHA1_NewContext(void);

/*
** Destroy a SHA-1 secure hash context.
**  "cx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void SHA1_DestroyContext(SHA1Context *cx, PRBool freeit);

/*
** Reset a SHA-1 context, preparing it for a fresh round of hashing
*/
extern void SHA1_Begin(SHA1Context *cx);

/*
** Update the SHA-1 hash function with more data.
**  "cx" the context
**  "input" the data to hash
**  "inputLen" the amount of data to hash
*/
extern void SHA1_Update(SHA1Context *cx, const unsigned char *input,
                        unsigned int inputLen);

/*
** Finish the SHA-1 hash function. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 16 bytes of digest data are stored
**  "digestLen" where the digest length (20) is stored
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void SHA1_End(SHA1Context *cx, unsigned char *digest,
                     unsigned int *digestLen, unsigned int maxDigestLen);

/*
** Export the current state of the SHA-1 hash without appending the standard
** padding and length bytes. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 20 bytes of digest data are stored
**  "digestLen" where the digest length (20) is stored (optional)
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void SHA1_EndRaw(SHA1Context *cx, unsigned char *digest,
                        unsigned int *digestLen, unsigned int maxDigestLen);

/*
** trace the intermediate state info of the SHA1 hash.
*/
extern void SHA1_TraceState(SHA1Context *cx);

/*
 * Return the the size of a buffer needed to flatten the SHA-1 Context into
 *    "cx" the context
 *  returns size;
 */
extern unsigned int SHA1_FlattenSize(SHA1Context *cx);

/*
 * Flatten the SHA-1 Context into a buffer:
 *    "cx" the context
 *    "space" the buffer to flatten to
 *  returns status;
 */
extern SECStatus SHA1_Flatten(SHA1Context *cx, unsigned char *space);

/*
 * Resurrect a flattened context into a SHA-1 Context
 *    "space" the buffer of the flattend buffer
 *    "arg" ptr to void used by cryptographic resurrect
 *  returns resurected context;
 */
extern SHA1Context *SHA1_Resurrect(unsigned char *space, void *arg);
extern void SHA1_Clone(SHA1Context *dest, SHA1Context *src);

/******************************************/

/******************************************/
/*
** SHA-2 secure hash function
** The SHA-2 family includes SHA224, SHA256, SHA384, and SHA512
*/

extern SHA224Context *SHA224_NewContext(void);
extern void SHA224_DestroyContext(SHA224Context *cx, PRBool freeit);
extern void SHA224_Begin(SHA224Context *cx);
extern void SHA224_Update(SHA224Context *cx, const unsigned char *input,
                          unsigned int inputLen);
extern void SHA224_End(SHA224Context *cx, unsigned char *digest,
                       unsigned int *digestLen, unsigned int maxDigestLen);
/*
** Export the current state of the SHA-224 hash without appending the standard
** padding and length bytes. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 28 bytes of digest data are stored
**  "digestLen" where the digest length (28) is stored (optional)
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void SHA224_EndRaw(SHA224Context *cx, unsigned char *digest,
                          unsigned int *digestLen, unsigned int maxDigestLen);
extern SECStatus SHA224_HashBuf(unsigned char *dest, const unsigned char *src,
                                PRUint32 src_length);
extern SECStatus SHA224_Hash(unsigned char *dest, const char *src);
extern void SHA224_TraceState(SHA224Context *cx);
extern unsigned int SHA224_FlattenSize(SHA224Context *cx);
extern SECStatus SHA224_Flatten(SHA224Context *cx, unsigned char *space);
extern SHA224Context *SHA224_Resurrect(unsigned char *space, void *arg);
extern void SHA224_Clone(SHA224Context *dest, SHA224Context *src);

/******************************************/

extern SHA256Context *SHA256_NewContext(void);
extern void SHA256_DestroyContext(SHA256Context *cx, PRBool freeit);
extern void SHA256_Begin(SHA256Context *cx);
extern void SHA256_Update(SHA256Context *cx, const unsigned char *input,
                          unsigned int inputLen);
extern void SHA256_End(SHA256Context *cx, unsigned char *digest,
                       unsigned int *digestLen, unsigned int maxDigestLen);
/*
** Export the current state of the SHA-256 hash without appending the standard
** padding and length bytes. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 32 bytes of digest data are stored
**  "digestLen" where the digest length (32) is stored (optional)
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void SHA256_EndRaw(SHA256Context *cx, unsigned char *digest,
                          unsigned int *digestLen, unsigned int maxDigestLen);
extern SECStatus SHA256_HashBuf(unsigned char *dest, const unsigned char *src,
                                PRUint32 src_length);
extern SECStatus SHA256_Hash(unsigned char *dest, const char *src);
extern void SHA256_TraceState(SHA256Context *cx);
extern unsigned int SHA256_FlattenSize(SHA256Context *cx);
extern SECStatus SHA256_Flatten(SHA256Context *cx, unsigned char *space);
extern SHA256Context *SHA256_Resurrect(unsigned char *space, void *arg);
extern void SHA256_Clone(SHA256Context *dest, SHA256Context *src);

/******************************************/

extern SHA512Context *SHA512_NewContext(void);
extern void SHA512_DestroyContext(SHA512Context *cx, PRBool freeit);
extern void SHA512_Begin(SHA512Context *cx);
extern void SHA512_Update(SHA512Context *cx, const unsigned char *input,
                          unsigned int inputLen);
/*
** Export the current state of the SHA-512 hash without appending the standard
** padding and length bytes. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 64 bytes of digest data are stored
**  "digestLen" where the digest length (64) is stored (optional)
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void SHA512_EndRaw(SHA512Context *cx, unsigned char *digest,
                          unsigned int *digestLen, unsigned int maxDigestLen);
extern void SHA512_End(SHA512Context *cx, unsigned char *digest,
                       unsigned int *digestLen, unsigned int maxDigestLen);
extern SECStatus SHA512_HashBuf(unsigned char *dest, const unsigned char *src,
                                PRUint32 src_length);
extern SECStatus SHA512_Hash(unsigned char *dest, const char *src);
extern void SHA512_TraceState(SHA512Context *cx);
extern unsigned int SHA512_FlattenSize(SHA512Context *cx);
extern SECStatus SHA512_Flatten(SHA512Context *cx, unsigned char *space);
extern SHA512Context *SHA512_Resurrect(unsigned char *space, void *arg);
extern void SHA512_Clone(SHA512Context *dest, SHA512Context *src);

/******************************************/

extern SHA384Context *SHA384_NewContext(void);
extern void SHA384_DestroyContext(SHA384Context *cx, PRBool freeit);
extern void SHA384_Begin(SHA384Context *cx);
extern void SHA384_Update(SHA384Context *cx, const unsigned char *input,
                          unsigned int inputLen);
extern void SHA384_End(SHA384Context *cx, unsigned char *digest,
                       unsigned int *digestLen, unsigned int maxDigestLen);
/*
** Export the current state of the SHA-384 hash without appending the standard
** padding and length bytes. Produce the digested results in "digest"
**  "cx" the context
**  "digest" where the 48 bytes of digest data are stored
**  "digestLen" where the digest length (48) is stored (optional)
**  "maxDigestLen" the maximum amount of data that can ever be
**     stored in "digest"
*/
extern void SHA384_EndRaw(SHA384Context *cx, unsigned char *digest,
                          unsigned int *digestLen, unsigned int maxDigestLen);
extern SECStatus SHA384_HashBuf(unsigned char *dest, const unsigned char *src,
                                PRUint32 src_length);
extern SECStatus SHA384_Hash(unsigned char *dest, const char *src);
extern void SHA384_TraceState(SHA384Context *cx);
extern unsigned int SHA384_FlattenSize(SHA384Context *cx);
extern SECStatus SHA384_Flatten(SHA384Context *cx, unsigned char *space);
extern SHA384Context *SHA384_Resurrect(unsigned char *space, void *arg);
extern void SHA384_Clone(SHA384Context *dest, SHA384Context *src);

/******************************************/
/*
** SHA-3 secure hash function
** The SHA-3 family includes SHA3_224, SHA3_256, SHA3_384, and SHA3_512
*/

extern SHA3_224Context *SHA3_224_NewContext(void);
extern void SHA3_224_DestroyContext(SHA3_224Context *cx, PRBool freeit);
extern unsigned int SHA3_224_FlattenSize(SHA3_224Context *cx);
extern void SHA3_224_Begin(SHA3_224Context *cx);
extern void SHA3_224_Update(SHA3_224Context *cx, const unsigned char *input,
                            unsigned int inputLen);
extern void SHA3_224_End(SHA3_224Context *cx, unsigned char *digest,
                         unsigned int *digestLen, unsigned int maxDigestLen);

extern SECStatus SHA3_224_HashBuf(unsigned char *dest, const unsigned char *src,
                                  PRUint32 src_length);
extern SECStatus SHA3_224_Hash(unsigned char *dest, const char *src);

/******************************************/

extern SHA3_256Context *SHA3_256_NewContext(void);
extern void SHA3_256_DestroyContext(SHA3_256Context *cx, PRBool freeit);
extern unsigned int SHA3_256_FlattenSize(SHA3_256Context *cx);
extern void SHA3_256_Begin(SHA3_256Context *cx);
extern void SHA3_256_Update(SHA3_256Context *cx, const unsigned char *input,
                            unsigned int inputLen);
extern void SHA3_256_End(SHA3_256Context *cx, unsigned char *digest,
                         unsigned int *digestLen, unsigned int maxDigestLen);

extern SECStatus SHA3_256_HashBuf(unsigned char *dest, const unsigned char *src,
                                  PRUint32 src_length);
extern SECStatus SHA3_256_Hash(unsigned char *dest, const char *src);

/******************************************/

extern SHA3_384Context *SHA3_384_NewContext(void);
extern void SHA3_384_DestroyContext(SHA3_384Context *cx, PRBool freeit);
extern unsigned int SHA3_384_FlattenSize(SHA3_384Context *cx);
extern void SHA3_384_Begin(SHA3_384Context *cx);
extern void SHA3_384_Update(SHA3_384Context *cx, const unsigned char *input,
                            unsigned int inputLen);
extern void SHA3_384_End(SHA3_384Context *cx, unsigned char *digest,
                         unsigned int *digestLen, unsigned int maxDigestLen);

extern SECStatus SHA3_384_HashBuf(unsigned char *dest, const unsigned char *src,
                                  PRUint32 src_length);
extern SECStatus SHA3_384_Hash(unsigned char *dest, const char *src);

/******************************************/

extern SHA3_512Context *SHA3_512_NewContext(void);
extern void SHA3_512_DestroyContext(SHA3_512Context *cx, PRBool freeit);
extern unsigned int SHA3_512_FlattenSize(SHA3_512Context *cx);
extern void SHA3_512_Begin(SHA3_512Context *cx);
extern void SHA3_512_Update(SHA3_512Context *cx, const unsigned char *input,
                            unsigned int inputLen);
extern void SHA3_512_End(SHA3_512Context *cx, unsigned char *digest,
                         unsigned int *digestLen, unsigned int maxDigestLen);

extern SECStatus SHA3_512_HashBuf(unsigned char *dest, const unsigned char *src,
                                  PRUint32 src_length);
extern SECStatus SHA3_512_Hash(unsigned char *dest, const char *src);

/******************************************/
/*
** SHAKE XOF functions from SHA-3
** The SHAKE family includes SHAKE_128 and SHAKE_256
*/

extern SHAKE_128Context *SHAKE_128_NewContext(void);
extern void SHAKE_128_DestroyContext(SHAKE_128Context *cx, PRBool freeit);
extern void SHAKE_128_Begin(SHAKE_128Context *cx);
extern void SHAKE_128_Absorb(SHAKE_128Context *cx, const unsigned char *input,
                             unsigned int inputLen);
extern void SHAKE_128_SqueezeEnd(SHAKE_128Context *cx, unsigned char *digest,
                                 unsigned int digestLen);
extern SECStatus SHAKE_128_HashBuf(unsigned char *dest, unsigned int dest_len,
                                   const unsigned char *src, PRUint32 src_length);
extern SECStatus SHAKE_128_Hash(unsigned char *dest, unsigned int dest_len, const char *src);

/******************************************/

extern SHAKE_256Context *SHAKE_256_NewContext(void);
extern void SHAKE_256_DestroyContext(SHAKE_256Context *cx, PRBool freeit);
extern void SHAKE_256_Begin(SHAKE_256Context *cx);
extern void SHAKE_256_Absorb(SHAKE_256Context *cx, const unsigned char *input,
                             unsigned int inputLen);
extern void SHAKE_256_SqueezeEnd(SHAKE_256Context *cx, unsigned char *digest,
                                 unsigned int digestLen);
extern SECStatus SHAKE_256_HashBuf(unsigned char *dest, unsigned int dest_len,
                                   const unsigned char *src, PRUint32 src_length);
extern SECStatus SHAKE_256_Hash(unsigned char *dest, unsigned int dest_len, const char *src);

/****************************************
 * implement TLS 1.0 Pseudo Random Function (PRF) and TLS P_hash function
 */

extern SECStatus
TLS_PRF(const SECItem *secret, const char *label, SECItem *seed,
        SECItem *result, PRBool isFIPS);

extern SECStatus
TLS_P_hash(HASH_HashType hashAlg, const SECItem *secret, const char *label,
           SECItem *seed, SECItem *result, PRBool isFIPS);

/******************************************/
/*
** Implements the Blake2b hash function.
*/

/*
** Hash a null terminated string "src" into "dest" using Blake2b
*/
extern SECStatus BLAKE2B_Hash(unsigned char *dest, const char *src);

/*
** Hash a non-null terminated string "src" into "dest" using Blake2b
*/
extern SECStatus BLAKE2B_HashBuf(unsigned char *output,
                                 const unsigned char *input, PRUint32 inlen);

extern SECStatus BLAKE2B_MAC_HashBuf(unsigned char *output,
                                     const unsigned char *input,
                                     unsigned int inlen,
                                     const unsigned char *key,
                                     unsigned int keylen);

/*
** Create a new Blake2b context
*/
extern BLAKE2BContext *BLAKE2B_NewContext(void);

/*
** Destroy a Blake2b secure hash context.
**  "ctx" the context
**  "freeit" if PR_TRUE then free the object as well as its sub-objects
*/
extern void BLAKE2B_DestroyContext(BLAKE2BContext *ctx, PRBool freeit);

/*
** Reset a Blake2b context, preparing it for a fresh round of hashing
*/
extern SECStatus BLAKE2B_Begin(BLAKE2BContext *ctx);

extern SECStatus BLAKE2B_MAC_Begin(BLAKE2BContext *ctx, const PRUint8 *key,
                                   const size_t keylen);

/*
** Update the Blake hash function with more data.
*/
extern SECStatus BLAKE2B_Update(BLAKE2BContext *ctx, const unsigned char *in,
                                unsigned int inlen);

/*
** Finish the Blake hash function. Produce the digested results in "digest"
*/
extern SECStatus BLAKE2B_End(BLAKE2BContext *ctx, unsigned char *out,
                             unsigned int *digestLen, size_t maxDigestLen);

/*
 * Return the size of a buffer needed to flatten the Blake2b Context into
 *    "ctx" the context
 *  returns size;
 */
extern unsigned int BLAKE2B_FlattenSize(BLAKE2BContext *ctx);

/*
 * Flatten the Blake2b Context into a buffer:
 *    "ctx" the context
 *    "space" the buffer to flatten to
 *  returns status;
 */
extern SECStatus BLAKE2B_Flatten(BLAKE2BContext *ctx, unsigned char *space);

/*
 * Resurrect a flattened context into a Blake2b Context
 *    "space" the buffer of the flattend buffer
 *    "arg" ptr to void used by cryptographic resurrect
 *  returns resurected context
 */
extern BLAKE2BContext *BLAKE2B_Resurrect(unsigned char *space, void *arg);
extern void BLAKE2B_Clone(BLAKE2BContext *dest, BLAKE2BContext *src);

/******************************************/
/*
** Pseudo Random Number Generation.  FIPS compliance desirable.
*/

/*
** Initialize the global RNG context and give it some seed input taken
** from the system.  This function is thread-safe and will only allow
** the global context to be initialized once.  The seed input is likely
** small, so it is imperative that RNG_RandomUpdate() be called with
** additional seed data before the generator is used.  A good way to
** provide the generator with additional entropy is to call
** RNG_SystemInfoForRNG().  Note that NSS_Init() does exactly that.
*/
extern SECStatus RNG_RNGInit(void);

/*
** Update the global random number generator with more seeding
** material
*/
extern SECStatus RNG_RandomUpdate(const void *data, size_t bytes);

/*
** Generate some random bytes, using the global random number generator
** object.
*/
extern SECStatus RNG_GenerateGlobalRandomBytes(void *dest, size_t len);

/* Destroy the global RNG context.  After a call to RNG_RNGShutdown()
** a call to RNG_RNGInit() is required in order to use the generator again,
** along with seed data (see the comment above RNG_RNGInit()).
*/
extern void RNG_RNGShutdown(void);

extern void RNG_SystemInfoForRNG(void);

/*
 * FIPS 186-2 Change Notice 1 RNG Algorithm 1, used both to
 * generate the DSA X parameter and as a generic purpose RNG.
 *
 * The following two FIPS186Change functions are needed for
 * NIST RNG Validation System.
 */

/*
 * FIPS186Change_GenerateX is now deprecated. It will return SECFailure with
 * the error set to PR_NOT_IMPLEMENTED_ERROR.
 */
extern SECStatus
FIPS186Change_GenerateX(unsigned char *XKEY,
                        const unsigned char *XSEEDj,
                        unsigned char *x_j);

/*
 * When generating the DSA X parameter, we generate 2*GSIZE bytes
 * of random output and reduce it mod q.
 *
 * Input: w, 2*GSIZE bytes
 *        q, DSA_SUBPRIME_LEN bytes
 * Output: xj, DSA_SUBPRIME_LEN bytes
 */
extern SECStatus
FIPS186Change_ReduceModQForDSA(const unsigned char *w,
                               const unsigned char *q,
                               unsigned char *xj);

/* To allow NIST KAT tests */
extern SECStatus
PRNGTEST_Instantiate_Kat(const PRUint8 *entropy, unsigned int entropy_len,
                         const PRUint8 *nonce, unsigned int nonce_len,
                         const PRUint8 *personal_string, unsigned int ps_len);

/*
 * The following functions are for FIPS poweron self test and FIPS algorithm
 * testing.
 */
extern SECStatus
PRNGTEST_Instantiate(const PRUint8 *entropy, unsigned int entropy_len,
                     const PRUint8 *nonce, unsigned int nonce_len,
                     const PRUint8 *personal_string, unsigned int ps_len);

extern SECStatus
PRNGTEST_Reseed(const PRUint8 *entropy, unsigned int entropy_len,
                const PRUint8 *additional, unsigned int additional_len);

extern SECStatus
PRNGTEST_Generate(PRUint8 *bytes, unsigned int bytes_len,
                  const PRUint8 *additional, unsigned int additional_len);

extern SECStatus
PRNGTEST_Uninstantiate(void);

extern SECStatus
PRNGTEST_RunHealthTests(void);

/* Generate PQGParams and PQGVerify structs.
 * Length of seed and length of h both equal length of P.
 * All lengths are specified by "j", according to the table above.
 *
 * The verify parameters will conform to FIPS186-1.
 */
extern SECStatus
PQG_ParamGen(unsigned int j,      /* input : determines length of P. */
             PQGParams **pParams, /* output: P Q and G returned here */
             PQGVerify **pVfy);   /* output: counter and seed. */

/* Generate PQGParams and PQGVerify structs.
 * Length of P specified by j.  Length of h will match length of P.
 * Length of SEED in bytes specified in seedBytes.
 * seedBbytes must be in the range [20..255] or an error will result.
 *
 * The verify parameters will conform to FIPS186-1.
 */
extern SECStatus
PQG_ParamGenSeedLen(
    unsigned int j,         /* input : determines length of P. */
    unsigned int seedBytes, /* input : length of seed in bytes.*/
    PQGParams **pParams,    /* output: P Q and G returned here */
    PQGVerify **pVfy);      /* output: counter and seed. */

/* Generate PQGParams and PQGVerify structs.
 * Length of P specified by L in bits.
 * Length of Q specified by N in bits.
 * Length of SEED in bytes specified in seedBytes.
 * seedBbytes must be in the range [N..L*2] or an error will result.
 *
 * Not that J uses the above table, L is the length exact. L and N must
 * match the table below or an error will result:
 *
 *  L            N
 * 1024         160
 * 2048         224
 * 2048         256
 * 3072         256
 *
 * If N or seedBytes are set to zero, then PQG_ParamGenSeedLen will
 * pick a default value (typically the smallest secure value for these
 * variables).
 *
 * The verify parameters will conform to FIPS186-3 using the smallest
 * permissible hash for the key strength.
 */
extern SECStatus
PQG_ParamGenV2(
    unsigned int L,         /* input : determines length of P. */
    unsigned int N,         /* input : determines length of Q. */
    unsigned int seedBytes, /* input : length of seed in bytes.*/
    PQGParams **pParams,    /* output: P Q and G returned here */
    PQGVerify **pVfy);      /* output: counter and seed. */

/*  Test PQGParams for validity as DSS PQG values.
 *  If vfy is non-NULL, test PQGParams to make sure they were generated
 *       using the specified seed, counter, and h values.
 *
 *  Return value indicates whether Verification operation ran successfully
 *  to completion, but does not indicate if PQGParams are valid or not.
 *  If return value is SECSuccess, then *pResult has these meanings:
 *       SECSuccess: PQGParams are valid.
 *       SECFailure: PQGParams are invalid.
 *
 * Verify the PQG againts the counter, SEED and h.
 * These tests are specified in FIPS 186-3 Appendix A.1.1.1, A.1.1.3, and A.2.2
 * PQG_VerifyParams will automatically choose the appropriate test.
 */

extern SECStatus PQG_VerifyParams(const PQGParams *params,
                                  const PQGVerify *vfy, SECStatus *result);

extern void PQG_DestroyParams(PQGParams *params);

extern void PQG_DestroyVerify(PQGVerify *vfy);

/*
 * clean-up any global tables freebl may have allocated after it starts up.
 * This function is not thread safe and should be called only after the
 * library has been quiessed.
 */
extern void BL_Cleanup(void);

/* unload freebl shared library from memory */
extern void BL_Unload(void);

/**************************************************************************
 *  Verify a given Shared library signature                               *
 **************************************************************************/
PRBool BLAPI_SHVerify(const char *name, PRFuncPtr addr);

/**************************************************************************
 *  Verify a given filename's signature                               *
 **************************************************************************/
PRBool BLAPI_SHVerifyFile(const char *shName);

/**************************************************************************
 *  Verify Are Own Shared library signature                               *
 **************************************************************************/
PRBool BLAPI_VerifySelf(const char *name);

/*********************************************************************/
extern const SECHashObject *HASH_GetRawHashObject(HASH_HashType hashType);

extern void BL_SetForkState(PRBool forked);

/*
** pepare an ECParam structure from DEREncoded params
 */
extern SECStatus EC_FillParams(PLArenaPool *arena,
                               const SECItem *encodedParams, ECParams *params);
extern SECStatus EC_DecodeParams(const SECItem *encodedParams,
                                 ECParams **ecparams);
extern SECStatus EC_CopyParams(PLArenaPool *arena, ECParams *dstParams,
                               const ECParams *srcParams);

/*
 * use the internal table to get the size in bytes of a single EC point
 */
extern int EC_GetPointSize(const ECParams *params);

SEC_END_PROTOS

#endif /* _BLAPI_H_ */
