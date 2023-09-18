/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _ALGHMAC_H_
#define _ALGHMAC_H_

typedef struct HMACContextStr HMACContext;

SEC_BEGIN_PROTOS

/* destroy HMAC context */
extern void
HMAC_Destroy(HMACContext *cx, PRBool freeit);

/* create HMAC context
 *  hash_obj    hash object from SECRawHashObjects[]
 *  secret      the secret with which the HMAC is performed.
 *  secret_len  the length of the secret.
 *  isFIPS      true if conforming to FIPS 198.
 *
 * NULL is returned if an error occurs.
 */
extern HMACContext *
HMAC_Create(const SECHashObject *hash_obj, const unsigned char *secret,
            unsigned int secret_len, PRBool isFIPS);

/* like HMAC_Create, except caller allocates HMACContext. */
SECStatus
HMAC_Init(HMACContext *cx, const SECHashObject *hash_obj,
          const unsigned char *secret, unsigned int secret_len, PRBool isFIPS);

/* like HMAC_Init, except caller passes in an existing context
 * previously used by either HMAC_Create or HMAC_Init. */
SECStatus
HMAC_ReInit(HMACContext *cx, const SECHashObject *hash_obj,
            const unsigned char *secret, unsigned int secret_len, PRBool isFIPS);

/* reset HMAC for a fresh round */
extern void
HMAC_Begin(HMACContext *cx);

/* update HMAC
 *  cx          HMAC Context
 *  data        the data to perform HMAC on
 *  data_len    the length of the data to process
 */
extern void
HMAC_Update(HMACContext *cx, const unsigned char *data, unsigned int data_len);

/* Finish HMAC -- place the results within result
 *  cx          HMAC context
 *  result      buffer for resulting hmac'd data
 *  result_len  where the resultant hmac length is stored
 *  max_result_len  maximum possible length that can be stored in result
 */
extern SECStatus
HMAC_Finish(HMACContext *cx, unsigned char *result, unsigned int *result_len,
            unsigned int max_result_len);

/* clone a copy of the HMAC state.  this is usefult when you would
 * need to keep a running hmac but also need to extract portions
 * partway through the process.
 */
extern HMACContext *
HMAC_Clone(HMACContext *cx);

SEC_END_PROTOS

#endif
