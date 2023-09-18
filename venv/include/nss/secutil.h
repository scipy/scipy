/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
#ifndef _SEC_UTIL_H_
#define _SEC_UTIL_H_

#include "seccomon.h"
#include "secitem.h"
#include "secport.h"
#include "prerror.h"
#include "base64.h"
#include "keyhi.h"
#include "secpkcs7.h"
#include "secasn1.h"
#include "secder.h"
#include <stdio.h>

#include "basicutil.h"
#include "sslerr.h"
#include "sslt.h"
#include "blapi.h"

#define SEC_CT_PRIVATE_KEY "private-key"
#define SEC_CT_PUBLIC_KEY "public-key"
#define SEC_CT_CERTIFICATE "certificate"
#define SEC_CT_CERTIFICATE_REQUEST "certificate-request"
#define SEC_CT_CERTIFICATE_ID "certificate-identity"
#define SEC_CT_PKCS7 "pkcs7"
#define SEC_CT_PKCS12 "pkcs12"
#define SEC_CT_CRL "crl"
#define SEC_CT_NAME "name"

#define NS_CERTREQ_HEADER "-----BEGIN NEW CERTIFICATE REQUEST-----"
#define NS_CERTREQ_TRAILER "-----END NEW CERTIFICATE REQUEST-----"

#define NS_CERT_HEADER "-----BEGIN CERTIFICATE-----"
#define NS_CERT_TRAILER "-----END CERTIFICATE-----"

#define NS_CRL_HEADER "-----BEGIN CRL-----"
#define NS_CRL_TRAILER "-----END CRL-----"

#define SECU_Strerror PORT_ErrorToString

typedef struct {
    enum {
        PW_NONE = 0,
        PW_FROMFILE = 1,
        PW_PLAINTEXT = 2,
        PW_EXTERNAL = 3
    } source;
    char *data;
} secuPWData;

/*
** Change a password on a token, or initialize a token with a password
** if it does not already have one.
** Use passwd to send the password in plaintext, pwFile to specify a
** file containing the password, or NULL for both to prompt the user.
*/
SECStatus SECU_ChangePW(PK11SlotInfo *slot, char *passwd, char *pwFile);

/*
** Change a password on a token, or initialize a token with a password
** if it does not already have one.
** In this function, you can specify both the old and new passwords
** as either a string or file. NOTE: any you don't specify will
** be prompted for
*/
SECStatus SECU_ChangePW2(PK11SlotInfo *slot, char *oldPass, char *newPass,
                         char *oldPwFile, char *newPwFile);

/*  These were stolen from the old sec.h... */
/*
** Check a password for legitimacy. Passwords must be at least 8
** characters long and contain one non-alphabetic. Return DSTrue if the
** password is ok, DSFalse otherwise.
*/
extern PRBool SEC_CheckPassword(char *password);

/*
** Blind check of a password. Complement to SEC_CheckPassword which
** ignores length and content type, just retuning DSTrue is the password
** exists, DSFalse if NULL
*/
extern PRBool SEC_BlindCheckPassword(char *password);

/*
** Get a password.
** First prompt with "msg" on "out", then read the password from "in".
** The password is then checked using "chkpw".
*/
extern char *SEC_GetPassword(FILE *in, FILE *out, char *msg,
                             PRBool (*chkpw)(char *));

char *SECU_FilePasswd(PK11SlotInfo *slot, PRBool retry, void *arg);

char *SECU_GetPasswordString(void *arg, char *prompt);

/*
** Write a dongle password.
** Uses MD5 to hash constant system data (hostname, etc.), and then
** creates RC4 key to encrypt a password "pw" into a file "fd".
*/
extern SECStatus SEC_WriteDongleFile(int fd, char *pw);

/*
** Get a dongle password.
** Uses MD5 to hash constant system data (hostname, etc.), and then
** creates RC4 key to decrypt and return a password from file "fd".
*/
extern char *SEC_ReadDongleFile(int fd);

/* End stolen headers */

/* Just sticks the two strings together with a / if needed */
char *SECU_AppendFilenameToDir(char *dir, char *filename);

/* Returns result of PR_GetEnvSecure("SSL_DIR") or NULL */
extern char *SECU_DefaultSSLDir(void);

/*
** Should be called once during initialization to set the default
**    directory for looking for cert.db, key.db, and cert-nameidx.db files
** Removes trailing '/' in 'base'
** If 'base' is NULL, defaults to set to .netscape in home directory.
*/
extern char *SECU_ConfigDirectory(const char *base);

/*
** Basic callback function for SSL_GetClientAuthDataHook
*/
extern int
SECU_GetClientAuthData(void *arg, PRFileDesc *fd,
                       struct CERTDistNamesStr *caNames,
                       struct CERTCertificateStr **pRetCert,
                       struct SECKEYPrivateKeyStr **pRetKey);

extern PRBool SECU_GetWrapEnabled(void);
extern void SECU_EnableWrap(PRBool enable);

extern PRBool SECU_GetUtf8DisplayEnabled(void);
extern void SECU_EnableUtf8Display(PRBool enable);

/* revalidate the cert and print information about cert verification
 * failure at time == now */
extern void
SECU_printCertProblems(FILE *outfile, CERTCertDBHandle *handle,
                       CERTCertificate *cert, PRBool checksig,
                       SECCertificateUsage certUsage, void *pinArg, PRBool verbose);

/* revalidate the cert and print information about cert verification
 * failure at specified time */
extern void
SECU_printCertProblemsOnDate(FILE *outfile, CERTCertDBHandle *handle,
                             CERTCertificate *cert, PRBool checksig, SECCertificateUsage certUsage,
                             void *pinArg, PRBool verbose, PRTime datetime);

/* print out CERTVerifyLog info. */
extern void
SECU_displayVerifyLog(FILE *outfile, CERTVerifyLog *log,
                      PRBool verbose);

/* Read in a DER from a file, may be ascii  */
extern SECStatus
SECU_ReadDERFromFile(SECItem *der, PRFileDesc *inFile, PRBool ascii,
                     PRBool warnOnPrivateKeyInAsciiFile);

/* Print integer value and hex */
extern void SECU_PrintInteger(FILE *out, const SECItem *i, const char *m,
                              int level);

/* Print ObjectIdentifier symbolically */
extern SECOidTag SECU_PrintObjectID(FILE *out, const SECItem *oid,
                                    const char *m, int level);

/* Print AlgorithmIdentifier symbolically */
extern void SECU_PrintAlgorithmID(FILE *out, SECAlgorithmID *a, char *m,
                                  int level);

/*
 * Format and print the UTC Time "t".  If the tag message "m" is not NULL,
 * do indent formatting based on "level" and add a newline afterward;
 * otherwise just print the formatted time string only.
 */
extern void SECU_PrintUTCTime(FILE *out, const SECItem *t, const char *m,
                              int level);

/*
 * Format and print the Generalized Time "t".  If the tag message "m"
 * is not NULL, * do indent formatting based on "level" and add a newline
 * afterward; otherwise just print the formatted time string only.
 */
extern void SECU_PrintGeneralizedTime(FILE *out, const SECItem *t,
                                      const char *m, int level);

/*
 * Format and print the UTC or Generalized Time "t".  If the tag message
 * "m" is not NULL, do indent formatting based on "level" and add a newline
 * afterward; otherwise just print the formatted time string only.
 */
extern void SECU_PrintTimeChoice(FILE *out, const SECItem *t, const char *m,
                                 int level);

/* callback for listing certs through pkcs11 */
extern SECStatus SECU_PrintCertNickname(CERTCertListNode *cert, void *data);

/* Dump all certificate nicknames in a database */
extern SECStatus
SECU_PrintCertificateNames(CERTCertDBHandle *handle, PRFileDesc *out,
                           PRBool sortByName, PRBool sortByTrust);

/* See if nickname already in database. Return 1 true, 0 false, -1 error */
int SECU_CheckCertNameExists(CERTCertDBHandle *handle, char *nickname);

/* Dump contents of cert req */
extern int SECU_PrintCertificateRequest(FILE *out, SECItem *der, char *m,
                                        int level);

/* Dump contents of certificate */
extern int SECU_PrintCertificate(FILE *out, const SECItem *der, const char *m,
                                 int level);

extern int SECU_PrintCertificateBasicInfo(FILE *out, const SECItem *der, const char *m,
                                          int level);

extern int SECU_PrintDumpDerIssuerAndSerial(FILE *out, SECItem *der, char *m,
                                            int level);

/* Dump contents of a DER certificate name (issuer or subject) */
extern int SECU_PrintDERName(FILE *out, SECItem *der, const char *m, int level);

/* print trust flags on a cert */
extern void SECU_PrintTrustFlags(FILE *out, CERTCertTrust *trust, char *m,
                                 int level);

extern int SECU_PrintSubjectPublicKeyInfo(FILE *out, SECItem *der, char *m,
                                          int level);

/* Dump contents of private key */
extern int SECU_PrintPrivateKey(FILE *out, SECItem *der, char *m, int level);

/* Dump contents of an RSA public key */
extern void SECU_PrintRSAPublicKey(FILE *out, SECKEYPublicKey *pk, char *m, int level);

/* Dump contents of a DSA public key */
extern void SECU_PrintDSAPublicKey(FILE *out, SECKEYPublicKey *pk, char *m, int level);

/* Print the MD5 and SHA1 fingerprints of a cert */
extern int SECU_PrintFingerprints(FILE *out, SECItem *derCert, char *m,
                                  int level);

/* Pretty-print any PKCS7 thing */
extern int SECU_PrintPKCS7ContentInfo(FILE *out, SECItem *der, char *m,
                                      int level);
/* Pretty-print a pkcs12 file */
extern SECStatus SECU_PrintPKCS12(FILE *out, const SECItem *der, char *m, int level);
/* Init PKCS11 stuff */
extern SECStatus SECU_PKCS11Init(PRBool readOnly);

/* Dump contents of signed data */
extern int SECU_PrintSignedData(FILE *out, SECItem *der, const char *m,
                                int level, SECU_PPFunc inner);

/* Dump contents of signed data, excluding the signature */
extern int SECU_PrintSignedContent(FILE *out, SECItem *der, char *m, int level,
                                   SECU_PPFunc inner);

/* Print cert data and its trust flags */
extern SECStatus SEC_PrintCertificateAndTrust(CERTCertificate *cert,
                                              const char *label,
                                              CERTCertTrust *trust);

extern int SECU_PrintCrl(FILE *out, SECItem *der, char *m, int level);

extern void
SECU_PrintCRLInfo(FILE *out, CERTCrl *crl, char *m, int level);

extern void SECU_PrintString(FILE *out, const SECItem *si, const char *m,
                             int level);
extern void SECU_PrintAny(FILE *out, const SECItem *i, const char *m, int level);

extern void SECU_PrintPolicy(FILE *out, SECItem *value, char *msg, int level);
extern void SECU_PrintPrivKeyUsagePeriodExtension(FILE *out, SECItem *value,
                                                  char *msg, int level);

extern void SECU_PrintExtensions(FILE *out, CERTCertExtension **extensions,
                                 char *msg, int level);

extern void SECU_PrintNameQuotesOptional(FILE *out, CERTName *name,
                                         const char *msg, int level,
                                         PRBool quotes);
extern void SECU_PrintName(FILE *out, CERTName *name, const char *msg,
                           int level);
extern void SECU_PrintRDN(FILE *out, CERTRDN *rdn, const char *msg, int level);

#ifdef SECU_GetPassword
/* Convert a High public Key to a Low public Key */
extern SECKEYLowPublicKey *SECU_ConvHighToLow(SECKEYPublicKey *pubHighKey);
#endif

extern char *SECU_GetModulePassword(PK11SlotInfo *slot, PRBool retry, void *arg);

extern SECStatus DER_PrettyPrint(FILE *out, const SECItem *it, PRBool raw);

extern char *SECU_SECModDBName(void);

/* Fetch and register an oid if it hasn't been done already */
extern void SECU_cert_fetchOID(SECOidTag *data, const SECOidData *src);

extern SECStatus SECU_RegisterDynamicOids(void);

/* Identifies hash algorithm tag by its string representation. */
extern SECOidTag SECU_StringToSignatureAlgTag(const char *alg);

/* Store CRL in output file or pk11 db. Also
 * encodes with base64 and exports to file if ascii flag is set
 * and file is not NULL. */
extern SECStatus SECU_StoreCRL(PK11SlotInfo *slot, SECItem *derCrl,
                               PRFileDesc *outFile, PRBool ascii, char *url);

/*
** DER sign a single block of data using private key encryption and the
** MD5 hashing algorithm. This routine first computes a digital signature
** using SEC_SignData, then wraps it with an CERTSignedData and then der
** encodes the result.
**     "arena" is the memory arena to use to allocate data from
**     "sd" returned CERTSignedData
**     "result" the final der encoded data (memory is allocated)
**     "buf" the input data to sign
**     "len" the amount of data to sign
**     "pk" the private key to encrypt with
*/
extern SECStatus SECU_DerSignDataCRL(PLArenaPool *arena, CERTSignedData *sd,
                                     unsigned char *buf, int len,
                                     SECKEYPrivateKey *pk, SECOidTag algID);

typedef enum {
    noKeyFound = 1,
    noSignatureMatch = 2,
    failToEncode = 3,
    failToSign = 4,
    noMem = 5
} SignAndEncodeFuncExitStat;

extern SECStatus
SECU_SignAndEncodeCRL(CERTCertificate *issuer, CERTSignedCrl *signCrl,
                      SECOidTag hashAlgTag, SignAndEncodeFuncExitStat *resCode);

extern SECStatus
SECU_CopyCRL(PLArenaPool *destArena, CERTCrl *destCrl, CERTCrl *srcCrl);

/*
** Finds the crl Authority Key Id extension. Returns NULL if no such extension
** was found.
*/
CERTAuthKeyID *
SECU_FindCRLAuthKeyIDExten(PLArenaPool *arena, CERTSignedCrl *crl);

/*
 * Find the issuer of a crl. Cert usage should be checked before signing a crl.
 */
CERTCertificate *
SECU_FindCrlIssuer(CERTCertDBHandle *dbHandle, SECItem *subject,
                   CERTAuthKeyID *id, PRTime validTime);

/* call back function used in encoding of an extension. Called from
 * SECU_EncodeAndAddExtensionValue */
typedef SECStatus (*EXTEN_EXT_VALUE_ENCODER)(PLArenaPool *extHandleArena,
                                             void *value, SECItem *encodedValue);

/* Encodes and adds extensions to the CRL or CRL entries. */
SECStatus
SECU_EncodeAndAddExtensionValue(PLArenaPool *arena, void *extHandle,
                                void *value, PRBool criticality, int extenType,
                                EXTEN_EXT_VALUE_ENCODER EncodeValueFn);

/* Caller ensures that dst is at least item->len*2+1 bytes long */
void
SECU_SECItemToHex(const SECItem *item, char *dst);

/* Requires 0x prefix. Case-insensitive. Will do in-place replacement if
 * successful */
SECStatus
SECU_SECItemHexStringToBinary(SECItem *srcdest);

/* Parse a version range string, with "min" and "max" version numbers,
 * separated by colon (":"), and return the result in vr and v2.
 *
 * Both min and max values are optional.
 * The following syntax is used to specify the enabled protocol versions:
 * A string with only a max value is expected as ":{max}",
 * and all implemented versions less than or equal to max will be enabled.
 * A string with only a min value is expected as "{min}:",
 * and all implemented versions greater than or equal to min will be enabled.
 * A string consisting of a colon only means "all versions enabled".
 *
 * In order to avoid a link dependency from libsectool to libssl,
 * the caller must provide the desired default values for the min/max values,
 * by providing defaultVersionRange (which can be obtained from libssl by
 * calling SSL_VersionRangeGetSupported).
 */
SECStatus
SECU_ParseSSLVersionRangeString(const char *input,
                                const SSLVersionRange defaultVersionRange,
                                SSLVersionRange *vrange);

SECStatus parseGroupList(const char *arg, SSLNamedGroup **enabledGroups,
                         unsigned int *enabledGroupsCount);
SECStatus parseSigSchemeList(const char *arg,
                             const SSLSignatureScheme **enabledSigSchemes,
                             unsigned int *enabledSigSchemeCount);
typedef struct {
    SECItem label;
    PRBool hasContext;
    SECItem context;
    unsigned int outputLength;
} secuExporter;

SECStatus parseExporters(const char *arg,
                         const secuExporter **enabledExporters,
                         unsigned int *enabledExporterCount);

SECStatus exportKeyingMaterials(PRFileDesc *fd,
                                const secuExporter *exporters,
                                unsigned int exporterCount);

SECStatus readPSK(const char *arg, SECItem *psk, SECItem *label);

/*
 *
 *  Error messaging
 *
 */

void printflags(char *trusts, unsigned int flags);

#if !defined(XP_UNIX) && !defined(XP_OS2)
extern int ffs(unsigned int i);
#endif

/* Finds certificate by searching it in the DB or by examinig file
 * in the local directory. */
CERTCertificate *
SECU_FindCertByNicknameOrFilename(CERTCertDBHandle *handle,
                                  char *name, PRBool ascii,
                                  void *pwarg);
#include "secerr.h"
#include "sslerr.h"

#endif /* _SEC_UTIL_H_ */
