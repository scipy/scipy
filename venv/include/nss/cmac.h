/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _CMAC_H_
#define _CMAC_H_

typedef struct CMACContextStr CMACContext;

SEC_BEGIN_PROTOS

/* Enum for identifying the underlying block cipher we're using internally. */
typedef enum {
    CMAC_AES = 0
} CMACCipher;

/* Initialize an existing CMACContext struct. */
SECStatus CMAC_Init(CMACContext *ctx, CMACCipher type,
                    const unsigned char *key, unsigned int key_len);

/* Allocate and initialize a new CMAC context with the specified cipher and
 * key. */
CMACContext *CMAC_Create(CMACCipher type, const unsigned char *key,
                         unsigned int key_len);

/* Called automatically by CMAC_*{Create,Init}(...). Only useful for restarting
 * an already-started CMAC instance. */
SECStatus CMAC_Begin(CMACContext *ctx);

/* Add the specified bytes into the CMAC state. */
SECStatus CMAC_Update(CMACContext *ctx, const unsigned char *data,
                      unsigned int data_len);

/* Finalize the CMAC state and return the result. */
SECStatus CMAC_Finish(CMACContext *ctx, unsigned char *result,
                      unsigned int *result_len,
                      unsigned int max_result_len);

/* Note: CMAC_Clone isn't implemented here because AES doesn't expose a
 * context-cloning operation. */

/* Destroy a CMAC context, optionally freeing it. */
void CMAC_Destroy(CMACContext *ctx, PRBool free_it);

SEC_END_PROTOS

#endif
