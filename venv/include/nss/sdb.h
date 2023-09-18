/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file implements PKCS 11 on top of our existing security modules
 *
 * For more information about PKCS 11 See PKCS 11 Token Inteface Standard.
 *   This implementation has two slots:
 *      slot 1 is our generic crypto support. It does not require login.
 *   It supports Public Key ops, and all they bulk ciphers and hashes.
 *   It can also support Private Key ops for imported Private keys. It does
 *   not have any token storage.
 *      slot 2 is our private key support. It requires a login before use. It
 *   can store Private Keys and Certs as token objects. Currently only private
 *   keys and their associated Certificates are saved on the token.
 *
 *   In this implementation, session objects are only visible to the session
 *   that created or generated them.
 */

/*
 * the following data structures should be moved to a 'rdb.h'.
 */

#ifndef _SDB_H
#define _SDB_H 1
#include "pkcs11t.h"
#include "secitem.h"
#include "sftkdbt.h"

#define STATIC_CMD_SIZE 2048

typedef struct SDBFindStr SDBFind;
typedef struct SDBStr SDB;

struct SDBStr {
    void *private;
    int version;
    int reserved;
    int sdb_flags;
    void *app_private;
    CK_RV(*sdb_FindObjectsInit)
    (SDB *sdb, const CK_ATTRIBUTE *template,
     CK_ULONG count, SDBFind **find);
    CK_RV(*sdb_FindObjects)
    (SDB *sdb, SDBFind *find, CK_OBJECT_HANDLE *ids,
     CK_ULONG arraySize, CK_ULONG *count);
    CK_RV(*sdb_FindObjectsFinal)
    (SDB *sdb, SDBFind *find);
    CK_RV(*sdb_GetAttributeValue)
    (SDB *sdb, CK_OBJECT_HANDLE object,
     CK_ATTRIBUTE *template, CK_ULONG count);
    CK_RV(*sdb_SetAttributeValue)
    (SDB *sdb, CK_OBJECT_HANDLE object,
     const CK_ATTRIBUTE *template, CK_ULONG count);
    CK_RV(*sdb_CreateObject)
    (SDB *sdb, CK_OBJECT_HANDLE *object,
     const CK_ATTRIBUTE *template, CK_ULONG count);
    CK_RV(*sdb_DestroyObject)
    (SDB *sdb, CK_OBJECT_HANDLE object);
    CK_RV(*sdb_GetMetaData)
    (SDB *sdb, const char *id,
     SECItem *item1, SECItem *item2);
    CK_RV(*sdb_PutMetaData)
    (SDB *sdb, const char *id,
     const SECItem *item1, const SECItem *item2);
    CK_RV(*sdb_Begin)
    (SDB *sdb);
    CK_RV(*sdb_Commit)
    (SDB *sdb);
    CK_RV(*sdb_Abort)
    (SDB *sdb);
    CK_RV(*sdb_Reset)
    (SDB *sdb);
    CK_RV(*sdb_Close)
    (SDB *sdb);
    void (*sdb_SetForkState)(PRBool forked);
    CK_RV(*sdb_GetNewObjectID)
    (SDB *db, CK_OBJECT_HANDLE *object);
    CK_RV(*sdb_DestroyMetaData)
    (SDB *db, const char *id);
};

CK_RV s_open(const char *directory, const char *certPrefix,
             const char *keyPrefix,
             int cert_version, int key_version,
             int flags, SDB **certdb, SDB **keydb, int *newInit);
CK_RV s_shutdown();

#if defined(_WIN32)
wchar_t *sdb_UTF8ToWide(const char *buf);
#endif

/* flags */
#define SDB_RDONLY 1
#define SDB_RDWR 2
#define SDB_CREATE 4
#define SDB_HAS_META 8
#define SDB_FIPS 0x10

#endif
