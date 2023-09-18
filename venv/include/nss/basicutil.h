/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
#ifndef _BASIC_UTILS_H_
#define _BASIC_UTILS_H_

#include "seccomon.h"
#include "secitem.h"
#include "secoid.h"
#include "secoidt.h"
#include "secport.h"
#include "prerror.h"
#include "base64.h"
#include "secasn1.h"
#include "secder.h"
#include "ecl-exp.h"
#include <stdio.h>

#ifdef SECUTIL_NEW
typedef int (*SECU_PPFunc)(PRFileDesc *out, SECItem *item,
                           char *msg, int level);
#else
typedef int (*SECU_PPFunc)(FILE *out, SECItem *item, char *msg, int level);
#endif

/* print out an error message */
extern void SECU_PrintError(const char *progName, const char *msg, ...);

/* print out a system error message */
extern void SECU_PrintSystemError(const char *progName, const char *msg, ...);

/* print a formatted error message */
extern void SECU_PrintErrMsg(FILE *out, int level, const char *progName,
                             const char *msg, ...);

/* Read the contents of a file into a SECItem */
extern SECStatus SECU_FileToItem(SECItem *dst, PRFileDesc *src);
extern SECStatus SECU_TextFileToItem(SECItem *dst, PRFileDesc *src);

/* Indent based on "level" */
extern void SECU_Indent(FILE *out, int level);

/* Print a newline to out */
extern void SECU_Newline(FILE *out);

/* Print integer value and hex */
extern void SECU_PrintInteger(FILE *out, const SECItem *i, const char *m,
                              int level);

/* Print SECItem as hex */
extern void SECU_PrintAsHex(FILE *out, const SECItem *i, const char *m,
                            int level);

/* dump a buffer in hex and ASCII */
extern void SECU_PrintBuf(FILE *out, const char *msg, const void *vp, int len);

/* Dump contents of private key */
extern int SECU_PrintPrivateKey(FILE *out, SECItem *der, char *m, int level);

/* Init PKCS11 stuff */
extern SECStatus SECU_PKCS11Init(PRBool readOnly);

/* Dump contents of signed data */
extern int SECU_PrintSignedData(FILE *out, SECItem *der, const char *m,
                                int level, SECU_PPFunc inner);

extern void SECU_PrintString(FILE *out, const SECItem *si, const char *m,
                             int level);
extern void SECU_PrintAny(FILE *out, const SECItem *i, const char *m, int level);

extern void SECU_PrintPRandOSError(const char *progName);

/* Caller ensures that dst is at least item->len*2+1 bytes long */
void
SECU_SECItemToHex(const SECItem *item, char *dst);

/* Requires 0x prefix. Case-insensitive. Will do in-place replacement if
 * successful */
SECStatus
SECU_SECItemHexStringToBinary(SECItem *srcdest);

/*
** Read a hex string into a SecItem.
*/
extern SECItem *SECU_HexString2SECItem(PLArenaPool *arena, SECItem *item,
                                       const char *str);

extern SECStatus SECU_ecName2params(ECCurveName curve, SECItem *params);

/*
 *
 *  Utilities for parsing security tools command lines
 *
 */

/*  A single command flag  */
typedef struct {
    char flag;
    PRBool needsArg;
    char *arg;
    PRBool activated;
    char *longform;
} secuCommandFlag;

/*  A full array of command/option flags  */
typedef struct
{
    int numCommands;
    int numOptions;

    secuCommandFlag *commands;
    secuCommandFlag *options;
} secuCommand;

/*  fill the "arg" and "activated" fields for each flag  */
SECStatus
SECU_ParseCommandLine(int argc, char **argv, char *progName,
                      const secuCommand *cmd);
char *
SECU_GetOptionArg(const secuCommand *cmd, int optionNum);

/*
 *
 *  Error messaging
 *
 */

void printflags(char *trusts, unsigned int flags);

#if !defined(XP_UNIX) && !defined(XP_OS2)
extern int ffs(unsigned int i);
#endif

#include "secerr.h"

extern const char *hex;
extern const char printable[];

#endif /* _BASIC_UTILS_H_ */
