/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This code defines the glue layer between softoken and the legacy DB library
 */
#include "sdb.h"

/*
 * function prototypes for the callbacks into softoken from the legacyDB
 */

typedef SECStatus (*LGEncryptFunc)(PLArenaPool *arena, SDB *sdb,
                                   SECItem *plainText, SECItem **cipherText);
typedef SECStatus (*LGDecryptFunc)(SDB *sdb, SECItem *cipherText,
                                   SECItem **plainText);

/*
 * function prototypes for the exported functions.
 */
typedef CK_RV (*LGOpenFunc)(const char *dir, const char *certPrefix,
                            const char *keyPrefix,
                            int certVersion, int keyVersion, int flags,
                            SDB **certDB, SDB **keyDB);
typedef char **(*LGReadSecmodFunc)(const char *appName,
                                   const char *filename,
                                   const char *dbname, char *params, PRBool rw);
typedef SECStatus (*LGReleaseSecmodFunc)(const char *appName,
                                         const char *filename,
                                         const char *dbname, char **params, PRBool rw);
typedef SECStatus (*LGDeleteSecmodFunc)(const char *appName,
                                        const char *filename,
                                        const char *dbname, char *params, PRBool rw);
typedef SECStatus (*LGAddSecmodFunc)(const char *appName,
                                     const char *filename,
                                     const char *dbname, char *params, PRBool rw);
typedef SECStatus (*LGShutdownFunc)(PRBool forked);
typedef void (*LGSetForkStateFunc)(PRBool);
typedef void (*LGSetCryptFunc)(LGEncryptFunc, LGDecryptFunc);

/*
 * Softoken Glue Functions
 */
CK_RV sftkdbCall_open(const char *dir, const char *certPrefix,
                      const char *keyPrefix,
                      int certVersion, int keyVersion, int flags,
                      SDB **certDB, SDB **keyDB);
char **sftkdbCall_ReadSecmodDB(const char *appName, const char *filename,
                               const char *dbname, char *params, PRBool rw);
SECStatus sftkdbCall_ReleaseSecmodDBData(const char *appName,
                                         const char *filename, const char *dbname,
                                         char **moduleSpecList, PRBool rw);
SECStatus sftkdbCall_DeleteSecmodDB(const char *appName,
                                    const char *filename, const char *dbname,
                                    char *args, PRBool rw);
SECStatus sftkdbCall_AddSecmodDB(const char *appName,
                                 const char *filename, const char *dbname,
                                 char *module, PRBool rw);
CK_RV sftkdbCall_Shutdown(void);
