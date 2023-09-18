/*
 * softoken.h - private data structures and prototypes for the softoken lib
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _SOFTOKEN_H_
#define _SOFTOKEN_H_

#include "blapi.h"
#include "lowkeyti.h"
#include "softoknt.h"
#include "secoidt.h"

#include "pkcs11t.h"

SEC_BEGIN_PROTOS

/*
 * Convenience wrapper for doing a single PKCS#1 v1.5 RSA operations where the
 * encoded digest info is computed internally, rather than by the caller.
 *
 * The HashSign variants expect as input the value of H, the computed hash
 * from RFC 3447, Section 9.2, Step 1, and will compute the DER-encoded
 * DigestInfo structure internally prior to signing/verifying.
 */
extern SECStatus
RSA_HashSign(SECOidTag hashOid, NSSLOWKEYPrivateKey *key,
             unsigned char *sig, unsigned int *sigLen, unsigned int maxLen,
             const unsigned char *hash, unsigned int hashLen);

extern SECStatus
RSA_HashCheckSign(SECOidTag hashOid, NSSLOWKEYPublicKey *key,
                  const unsigned char *sig, unsigned int sigLen,
                  const unsigned char *hash, unsigned int hashLen);

/*
** Prepare a buffer for padded CBC encryption, growing to the appropriate
** boundary, filling with the appropriate padding.
**
** blockSize must be a power of 2.
**
** We add from 1 to blockSize bytes -- we *always* grow.
** The extra bytes contain the value of the length of the padding:
** if we have 2 bytes of padding, then the padding is "0x02, 0x02".
**
** NOTE: If arena is non-NULL, we re-allocate from there, otherwise
** we assume (and use) PR memory (re)allocation.
*/
extern unsigned char *CBC_PadBuffer(PLArenaPool *arena, unsigned char *inbuf,
                                    unsigned int inlen, unsigned int *outlen,
                                    int blockSize);

/****************************************/
/*
** Power-Up selftests are required for FIPS.
*/
/* make sure Power-up selftests have been run. */
extern CK_RV sftk_FIPSEntryOK(PRBool rerun);

/*
** make known fixed PKCS #11 key types to their sizes in bytes
*/
unsigned long sftk_MapKeySize(CK_KEY_TYPE keyType);

/*
** FIPS 140-2 auditing
*/
extern PRBool sftk_audit_enabled;

extern void sftk_LogAuditMessage(NSSAuditSeverity severity,
                                 NSSAuditType, const char *msg);

extern void sftk_AuditCreateObject(CK_SESSION_HANDLE hSession,
                                   CK_ATTRIBUTE_PTR pTemplate, CK_ULONG ulCount,
                                   CK_OBJECT_HANDLE_PTR phObject, CK_RV rv);

extern void sftk_AuditCopyObject(CK_SESSION_HANDLE hSession,
                                 CK_OBJECT_HANDLE hObject,
                                 CK_ATTRIBUTE_PTR pTemplate, CK_ULONG ulCount,
                                 CK_OBJECT_HANDLE_PTR phNewObject, CK_RV rv);

extern void sftk_AuditDestroyObject(CK_SESSION_HANDLE hSession,
                                    CK_OBJECT_HANDLE hObject, CK_RV rv);

extern void sftk_AuditGetObjectSize(CK_SESSION_HANDLE hSession,
                                    CK_OBJECT_HANDLE hObject, CK_ULONG_PTR pulSize,
                                    CK_RV rv);

extern void sftk_AuditGetAttributeValue(CK_SESSION_HANDLE hSession,
                                        CK_OBJECT_HANDLE hObject, CK_ATTRIBUTE_PTR pTemplate,
                                        CK_ULONG ulCount, CK_RV rv);

extern void sftk_AuditSetAttributeValue(CK_SESSION_HANDLE hSession,
                                        CK_OBJECT_HANDLE hObject, CK_ATTRIBUTE_PTR pTemplate,
                                        CK_ULONG ulCount, CK_RV rv);

extern void sftk_AuditCryptInit(const char *opName,
                                CK_SESSION_HANDLE hSession,
                                CK_MECHANISM_PTR pMechanism,
                                CK_OBJECT_HANDLE hKey, CK_RV rv);

extern void sftk_AuditGenerateKey(CK_SESSION_HANDLE hSession,
                                  CK_MECHANISM_PTR pMechanism,
                                  CK_ATTRIBUTE_PTR pTemplate, CK_ULONG ulCount,
                                  CK_OBJECT_HANDLE_PTR phKey, CK_RV rv);

extern void sftk_AuditGenerateKeyPair(CK_SESSION_HANDLE hSession,
                                      CK_MECHANISM_PTR pMechanism,
                                      CK_ATTRIBUTE_PTR pPublicKeyTemplate,
                                      CK_ULONG ulPublicKeyAttributeCount,
                                      CK_ATTRIBUTE_PTR pPrivateKeyTemplate,
                                      CK_ULONG ulPrivateKeyAttributeCount,
                                      CK_OBJECT_HANDLE_PTR phPublicKey,
                                      CK_OBJECT_HANDLE_PTR phPrivateKey, CK_RV rv);

extern void sftk_AuditWrapKey(CK_SESSION_HANDLE hSession,
                              CK_MECHANISM_PTR pMechanism,
                              CK_OBJECT_HANDLE hWrappingKey, CK_OBJECT_HANDLE hKey,
                              CK_BYTE_PTR pWrappedKey,
                              CK_ULONG_PTR pulWrappedKeyLen, CK_RV rv);

extern void sftk_AuditUnwrapKey(CK_SESSION_HANDLE hSession,
                                CK_MECHANISM_PTR pMechanism,
                                CK_OBJECT_HANDLE hUnwrappingKey,
                                CK_BYTE_PTR pWrappedKey, CK_ULONG ulWrappedKeyLen,
                                CK_ATTRIBUTE_PTR pTemplate, CK_ULONG ulAttributeCount,
                                CK_OBJECT_HANDLE_PTR phKey, CK_RV rv);

extern void sftk_AuditDeriveKey(CK_SESSION_HANDLE hSession,
                                CK_MECHANISM_PTR pMechanism,
                                CK_OBJECT_HANDLE hBaseKey,
                                CK_ATTRIBUTE_PTR pTemplate, CK_ULONG ulAttributeCount,
                                CK_OBJECT_HANDLE_PTR phKey, CK_RV rv);

extern void sftk_AuditDigestKey(CK_SESSION_HANDLE hSession,
                                CK_OBJECT_HANDLE hKey, CK_RV rv);

/*
** FIPS 140-2 Error state
*/
extern PRBool sftk_fatalError;

/*
** macros to check for forked child process after C_Initialize
*/
/* for PKCS #11 3.0, default is NO_FORK_CHECK, if you want it, now you
 * need to define DO_FORK_CHECK */
#if defined(XP_UNIX) && defined(DO_FORK_CHECK)

#ifdef DEBUG

#define FORK_ASSERT()                                            \
    {                                                            \
        char *forkAssert = PR_GetEnvSecure("NSS_STRICT_NOFORK"); \
        if ((!forkAssert) || (0 == strcmp(forkAssert, "1"))) {   \
            PORT_Assert(0);                                      \
        }                                                        \
    }

#else

#define FORK_ASSERT()

#endif

/* we have 3 methods of implementing the fork checks :
 * - Solaris "mixed" method
 * - pthread_atfork method
 * - getpid method
 */

#if !defined(CHECK_FORK_MIXED) && !defined(CHECK_FORK_PTHREAD) && \
    !defined(CHECK_FORK_GETPID)

/* Choose fork check method automatically unless specified
 * This section should be updated as more platforms get pthread fixes
 * to unregister fork handlers in dlclose.
 */

#ifdef SOLARIS

/* Solaris 8, s9 use PID checks, s10 uses pthread_atfork */

#define CHECK_FORK_MIXED

#elif defined(LINUX) || defined(__GLIBC__) || defined(FREEBSD) || defined(OPENBSD)

#define CHECK_FORK_PTHREAD

#else

/* Other Unix platforms use only PID checks. Even if pthread_atfork is
 * available, the behavior of dlclose isn't guaranteed by POSIX to
 * unregister the fork handler. */

#define CHECK_FORK_GETPID

#endif

#endif

#if defined(CHECK_FORK_MIXED)

extern PRBool usePthread_atfork;
#include <unistd.h>
extern pid_t myPid;
extern PRBool forked;

#define PARENT_FORKED() (usePthread_atfork ? forked : (myPid && myPid != getpid()))

#elif defined(CHECK_FORK_PTHREAD)

extern PRBool forked;

#define PARENT_FORKED() forked

#elif defined(CHECK_FORK_GETPID)

#include <unistd.h>
extern pid_t myPid;

#define PARENT_FORKED() (myPid && myPid != getpid())

#endif

extern PRBool parentForkedAfterC_Initialize;
extern PRBool sftkForkCheckDisabled;

#define CHECK_FORK()                                     \
    do {                                                 \
        if (!sftkForkCheckDisabled && PARENT_FORKED()) { \
            FORK_ASSERT();                               \
            return CKR_DEVICE_ERROR;                     \
        }                                                \
    } while (0)

#define SKIP_AFTER_FORK(x)              \
    if (!parentForkedAfterC_Initialize) \
    x

#define ENABLE_FORK_CHECK()                                       \
    {                                                             \
        char *doForkCheck = PR_GetEnvSecure("NSS_STRICT_NOFORK"); \
        if (doForkCheck && !strcmp(doForkCheck, "DISABLED")) {    \
            sftkForkCheckDisabled = PR_TRUE;                      \
        }                                                         \
    }

#else

/* non-Unix platforms, or fork check disabled */

#define CHECK_FORK()
#define SKIP_AFTER_FORK(x) x
#define ENABLE_FORK_CHECK()

#ifndef NO_FORK_CHECK
#define NO_FORK_CHECK
#endif

#endif

/*
 * If we were trying to be complete, we would have both FORK_SAFE
 * and non-Fork safe interfaces here. That would require doubling
 * the functions in our function list for both this and the FIPS
 * interface. Since NSS now always asks for a FORK_SAFE interface,
 * and can fall back to a non-FORK_SAFE interface, we set only
 * export one set of interfaces here */
#ifdef NO_FORK_CHECK
#define NSS_INTERFACE_FLAGS CKF_INTERFACE_FORK_SAFE
#else
#define NSS_INTERFACE_FLAGS 0
#endif

SEC_END_PROTOS

#endif /* _SOFTOKEN_H_ */
