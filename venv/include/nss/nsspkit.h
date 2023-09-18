/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef NSSPKIT_H
#define NSSPKIT_H

/*
 * nsspkit.h
 *
 * This file defines the types of the top-level PKI objects.
 */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

PR_BEGIN_EXTERN_C

/*
 * NSSCertificate
 *
 * This is the public representation of a Certificate.  The certificate
 * may be one found on a smartcard or other token, one decoded from data
 * received as part of a protocol, one constructed from constituent
 * parts, etc.  Usually it is associated with ("in") a trust domain; as
 * it can be verified only within a trust domain.  The underlying type
 * of certificate may be of any supported standard, e.g. PKIX, PGP, etc.
 *
 * People speak of "verifying (with) the server's, or correspondant's,
 * certificate"; for simple operations we support that simplification
 * by implementing public-key crypto operations as methods on this type.
 */

struct NSSCertificateStr;
typedef struct NSSCertificateStr NSSCertificate;

/*
 * NSSUserCertificate
 *
 * A ``User'' certificate is one for which the private key is available.
 * People speak of "using my certificate to sign my email" and "using
 * my certificate to authenticate to (or login to) the server"; for
 * simple operations, we support that simplification by implementing
 * private-key crypto operations as methods on this type.
 *
 * The current design only weakly distinguishes between certificates
 * and user certificates: as far as the compiler goes they're
 * interchangeable; debug libraries only have one common pointer-tracker;
 * etc.  However, attempts to do private-key operations on a certificate
 * for which the private key is not available will fail.
 *
 * Open design question: should these types be more firmly separated?
 */

typedef NSSCertificate NSSUserCertificate;

/*
 * NSSPrivateKey
 *
 * This is the public representation of a Private Key.  In general,
 * the actual value of the key is not available, but operations may
 * be performed with it.
 */

struct NSSPrivateKeyStr;
typedef struct NSSPrivateKeyStr NSSPrivateKey;

/*
 * NSSPublicKey
 *
 */

struct NSSPublicKeyStr;
typedef struct NSSPublicKeyStr NSSPublicKey;

/*
 * NSSSymmetricKey
 *
 */

struct NSSSymmetricKeyStr;
typedef struct NSSSymmetricKeyStr NSSSymmetricKey;

/*
 * NSSTrustDomain
 *
 * A Trust Domain is the field in which certificates may be validated.
 * A trust domain will generally have one or more cryptographic modules
 * open; these modules perform the cryptographic operations, and
 * provide the basic "root" trust information from which the trust in
 * a specific certificate or key depends.
 *
 * A client program, or a simple server, would typically have one
 * trust domain.  A server supporting multiple "virtual servers" might
 * have a separate trust domain for each virtual server.  The separate
 * trust domains might share some modules (e.g., a hardware crypto
 * accelerator) but not others (e.g., the tokens storing the different
 * servers' private keys, or the databases with each server's trusted
 * root certificates).
 *
 * This object descends from the "permananet database" in the old code.
 */

struct NSSTrustDomainStr;
typedef struct NSSTrustDomainStr NSSTrustDomain;

/*
 * NSSCryptoContext
 *
 * A Crypto Context is a short-term, "helper" object which is used
 * for the lifetime of one ongoing "crypto operation."  Such an
 * operation may be the creation of a signed message, the use of an
 * TLS socket connection, etc.  Each crypto context is "in" a
 * specific trust domain, and it may have associated with it a
 * distinguished certificate, public key, private key, and/or
 * symmetric key.  It can also temporarily hold and use temporary
 * data (e.g. intermediate certificates) which is not stored
 * permanently in the trust domain.
 *
 * In OO terms, this interface inherits interfaces from the trust
 * domain, the certificates, and the keys.  It also provides
 * streaming crypto operations.
 *
 * This object descends from the "temporary database" concept in the
 * old code, but it has changed a lot as a result of what we've
 * learned.
 */

typedef struct NSSCryptoContextStr NSSCryptoContext;

/*
 * fgmr others
 */

/*
 * OBJECT IDENTIFIER
 *
 * This is the basic OID that crops up everywhere.
 */

struct NSSOIDStr; /* unused opaque structure */
typedef struct NSSOIDStr NSSOID;

/*
 * NSSTime
 *
 * Unfortunately, we need an "exceptional" value to indicate
 * an error upon return, or "no value" on input.  Note that zero
 * is a perfectly valid value for both time_t and PRTime.
 *
 * If we were to create a "range" object, with two times for
 * Not Before and Not After, we would have an obvious place for
 * the somewhat arbitrary logic involved in comparing them.
 *
 * Failing that, let's have an NSSTime_CompareRanges function.
 */

struct NSSTimeStr;
typedef struct NSSTimeStr NSSTime;

struct NSSTrustStr;
typedef struct NSSTrustStr NSSTrust;

/*
 * NSSUsage
 *
 * This is trickier than originally planned; I'll write up a
 * doc on it.
 *
 * We'd still like nsspki.h to have a list of common usages,
 * e.g.:
 *
 *  extern const NSSUsage *NSSUsage_ClientAuth;
 *  extern const NSSUsage *NSSUsage_ServerAuth;
 *  extern const NSSUsage *NSSUsage_SignEmail;
 *  extern const NSSUsage *NSSUsage_EncryptEmail;
 *  etc.
 */

struct NSSUsageStr;
typedef struct NSSUsageStr NSSUsage;

/*
 * NSSPolicies
 *
 * Placeholder, for now.
 */

struct NSSPoliciesStr;
typedef struct NSSPoliciesStr NSSPolicies;

/*
 * NSSAlgorithmAndParameters
 *
 * Algorithm is an OID
 * Parameters depend on the algorithm
 */

struct NSSAlgorithmAndParametersStr;
typedef struct NSSAlgorithmAndParametersStr NSSAlgorithmAndParameters;

/*
 * NSSCallback
 *
 * At minimum, a "challenge" method and a closure argument.
 * Usually the challenge will just be prompting for a password.
 * How OO do we want to make it?
 */

typedef struct NSSCallbackStr NSSCallback;

struct NSSCallbackStr {
    /* Prompt for a password to initialize a slot.  */
    PRStatus (*getInitPW)(NSSUTF8 *slotName, void *arg,
                          NSSUTF8 **ssoPW, NSSUTF8 **userPW);
    /* Prompt for oldPW and newPW in order to change the
     * password on a slot.
     */
    PRStatus (*getNewPW)(NSSUTF8 *slotName, PRUint32 *retries, void *arg,
                         NSSUTF8 **oldPW, NSSUTF8 **newPW);
    /* Prompt for slot password.  */
    PRStatus (*getPW)(NSSUTF8 *slotName, PRUint32 *retries, void *arg,
                      NSSUTF8 **password);
    void *arg;
};

/* set errors - user cancelled, ... */

typedef PRUint32 NSSOperations;
/* 1) Do we want these to be preprocessor definitions or constants? */
/* 2) What is the correct and complete list? */

#define NSSOperations_ENCRYPT 0x0001
#define NSSOperations_DECRYPT 0x0002
#define NSSOperations_WRAP 0x0004
#define NSSOperations_UNWRAP 0x0008
#define NSSOperations_SIGN 0x0010
#define NSSOperations_SIGN_RECOVER 0x0020
#define NSSOperations_VERIFY 0x0040
#define NSSOperations_VERIFY_RECOVER 0x0080

struct NSSPKIXCertificateStr;

PR_END_EXTERN_C

#endif /* NSSPKIT_H */
