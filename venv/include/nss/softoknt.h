/*
 * softoknt.h - public data structures for the software token library
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _SOFTOKNT_H_
#define _SOFTOKNT_H_

#define NSS_SOFTOKEN_DEFAULT_CHUNKSIZE 2048
#define DES_BLOCK_SIZE 8     /* bytes */
#define MAX_DES3_KEY_SIZE 24 /* DES_BLOCK_SIZE * 3 */
#define SFTK_MAX_DERIVE_KEY_SIZE 64

/*
 * FIPS 140-2 auditing
 */
typedef enum {
    NSS_AUDIT_ERROR = 3,   /* errors */
    NSS_AUDIT_WARNING = 2, /* warning messages */
    NSS_AUDIT_INFO = 1     /* informational messages */
} NSSAuditSeverity;

typedef enum {
    NSS_AUDIT_ACCESS_KEY = 0,
    NSS_AUDIT_CHANGE_KEY,
    NSS_AUDIT_COPY_KEY,
    NSS_AUDIT_CRYPT,
    NSS_AUDIT_DERIVE_KEY,
    NSS_AUDIT_DESTROY_KEY,
    NSS_AUDIT_DIGEST_KEY,
    NSS_AUDIT_FIPS_STATE,
    NSS_AUDIT_GENERATE_KEY,
    NSS_AUDIT_INIT_PIN,
    NSS_AUDIT_INIT_TOKEN,
    NSS_AUDIT_LOAD_KEY,
    NSS_AUDIT_LOGIN,
    NSS_AUDIT_LOGOUT,
    NSS_AUDIT_SELF_TEST,
    NSS_AUDIT_SET_PIN,
    NSS_AUDIT_UNWRAP_KEY,
    NSS_AUDIT_WRAP_KEY
} NSSAuditType;

#endif /* _SOFTOKNT_H_ */
