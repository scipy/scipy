/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * certi.h - private data structures for the certificate library
 */
#ifndef _CERTI_H_
#define _CERTI_H_

#include "certt.h"
#include "nssrwlkt.h"

/*
#define GLOBAL_RWLOCK 1
*/

#define DPC_RWLOCK 1

/* all definitions in this file are subject to change */

typedef struct OpaqueCRLFieldsStr OpaqueCRLFields;
typedef struct CRLEntryCacheStr CRLEntryCache;
typedef struct CRLDPCacheStr CRLDPCache;
typedef struct CRLIssuerCacheStr CRLIssuerCache;
typedef struct CRLCacheStr CRLCache;
typedef struct CachedCrlStr CachedCrl;
typedef struct NamedCRLCacheStr NamedCRLCache;
typedef struct NamedCRLCacheEntryStr NamedCRLCacheEntry;

struct OpaqueCRLFieldsStr {
    PRBool partial;
    PRBool decodingError;
    PRBool badEntries;
    PRBool badDER;
    PRBool badExtensions;
    PRBool heapDER;
};

typedef struct PreAllocatorStr PreAllocator;

struct PreAllocatorStr {
    PRSize len;
    void* data;
    PRSize used;
    PLArenaPool* arena;
    PRSize extra;
};

/*  CRL entry cache.
    This is the same as an entry plus the next/prev pointers for the hash table
*/

struct CRLEntryCacheStr {
    CERTCrlEntry entry;
    CRLEntryCache *prev, *next;
};

#define CRL_CACHE_INVALID_CRLS 0x0001      /* this state will be set             \
                 if we have CRL objects with an invalid DER or signature. Can be \
                 cleared if the invalid objects are deleted from the token */
#define CRL_CACHE_LAST_FETCH_FAILED 0x0002 /* this state will be set        \
            if the last CRL fetch encountered an error. Can be cleared if a \
            new fetch succeeds */

#define CRL_CACHE_OUT_OF_MEMORY 0x0004 /* this state will be set \
            if we don't have enough memory to build the hash table of entries */

typedef enum {
    CRL_OriginToken = 0,   /* CRL came from PKCS#11 token */
    CRL_OriginExplicit = 1 /* CRL was explicitly added to the cache, from RAM */
} CRLOrigin;

typedef enum {
    dpcacheNoEntry = 0,           /* no entry found for this SN */
    dpcacheFoundEntry = 1,        /* entry found for this SN */
    dpcacheCallerError = 2,       /* invalid args */
    dpcacheInvalidCacheError = 3, /* CRL in cache may be bad DER */
                                  /* or unverified */
    dpcacheEmpty = 4,             /* no CRL in cache */
    dpcacheLookupError = 5        /* internal error */
} dpcacheStatus;

struct CachedCrlStr {
    CERTSignedCrl* crl;
    CRLOrigin origin;
    /* hash table of entries. We use a PLHashTable and pre-allocate the
       required amount of memory in one shot, so that our allocator can
       simply pass offsets into it when hashing.

       This won't work anymore when we support delta CRLs and iCRLs, because
       the size of the hash table will vary over time. At that point, the best
       solution will be to allocate large CRLEntry structures by modifying
       the DER decoding template. The extra space would be for next/prev
       pointers. This would allow entries from different CRLs to be mixed in
       the same hash table.
    */
    PLHashTable* entries;
    PreAllocator* prebuffer; /* big pre-allocated buffer mentioned above */
    PRBool sigChecked;       /* this CRL signature has already been checked */
    PRBool sigValid;         /* signature verification status .
                                Only meaningful if checked is PR_TRUE . */
    PRBool unbuildable;      /* Avoid using assosiated CRL is it fails
                              * a decoding step */
};

/*  CRL distribution point cache object
    This is a cache of CRL entries for a given distribution point of an issuer
    It is built from a collection of one full and 0 or more delta CRLs.
*/

struct CRLDPCacheStr {
#ifdef DPC_RWLOCK
    NSSRWLock* lock;
#else
    PRLock* lock;
#endif
    SECItem* issuerDERCert; /* issuer DER cert. Don't hold a reference
                               to the actual cert so the trust can be
                               updated on the cert automatically.
                               XXX there may be multiple issuer certs,
                               with different validity dates. Also
                               need to deal with SKID/AKID . See
                               bugzilla 217387, 233118 */

    CERTCertDBHandle* dbHandle;

    SECItem* subject;           /* DER of issuer subject */
    SECItem* distributionPoint; /* DER of distribution point. This may be
                                   NULL when distribution points aren't
                                   in use (ie. the CA has a single CRL).
                                   Currently not used. */

    /* array of full CRLs matching this distribution point */
    PRUint32 ncrls;   /* total number of CRLs in crls */
    CachedCrl** crls; /* array of all matching CRLs */
    /* XCRL With iCRLs and multiple DPs, the CRL can be shared accross several
       issuers. In the future, we'll need to globally recycle the CRL in a
       separate list in order to avoid extra lookups, decodes, and copies */

    /* pointers to good decoded CRLs used to build the cache */
    CachedCrl* selected; /* full CRL selected for use in the cache */
#if 0
    /* for future use */
    PRInt32 numdeltas;      /* number of delta CRLs used for the cache */
    CachedCrl** deltas;     /* delta CRLs used for the cache */
#endif
    /* cache invalidity bitflag */
    PRUint16 invalid;  /* this state will be set if either
        CRL_CACHE_INVALID_CRLS or CRL_CACHE_LAST_FETCH_FAILED is set.
        In those cases, all certs are considered to have unknown status.
        The invalid state can only be cleared during an update if all
        error states are cleared */
    PRBool refresh;    /* manual refresh from tokens has been forced */
    PRBool mustchoose; /* trigger reselection algorithm, for case when
                          RAM CRL objects are dropped from the cache */
    PRTime lastfetch;  /* time a CRL token fetch was last performed */
    PRTime lastcheck;  /* time CRL token objects were last checked for
                          existence */
};

/*  CRL issuer cache object
    This object tracks all the distribution point caches for a given issuer.
    XCRL once we support multiple issuing distribution points, this object
    will be a hash table. For now, it just holds the single CRL distribution
    point cache structure.
*/

struct CRLIssuerCacheStr {
    SECItem* subject; /* DER of issuer subject */
    CRLDPCache* dpp;
};

/*  CRL revocation cache object
    This object tracks all the issuer caches
*/

struct CRLCacheStr {
#ifdef GLOBAL_RWLOCK
    NSSRWLock* lock;
#else
    PRLock* lock;
#endif
    /* hash table of issuer to CRLIssuerCacheStr,
       indexed by issuer DER subject */
    PLHashTable* issuers;
};

SECStatus InitCRLCache(void);
SECStatus ShutdownCRLCache(void);

/* Returns a pointer to an environment-like string, a series of
** null-terminated strings, terminated by a zero-length string.
** This function is intended to be internal to NSS.
*/
extern char* cert_GetCertificateEmailAddresses(CERTCertificate* cert);

/*
 * These functions are used to map subjectKeyID extension values to certs
 * and to keep track of the checks for user certificates in each slot
 */
SECStatus cert_CreateSubjectKeyIDHashTable(void);

SECStatus cert_AddSubjectKeyIDMapping(SECItem* subjKeyID,
                                      CERTCertificate* cert);

SECStatus cert_UpdateSubjectKeyIDSlotCheck(SECItem* slotid, int series);

int cert_SubjectKeyIDSlotCheckSeries(SECItem* slotid);

/*
 * Call this function to remove an entry from the mapping table.
 */
SECStatus cert_RemoveSubjectKeyIDMapping(SECItem* subjKeyID);

SECStatus cert_DestroySubjectKeyIDHashTable(void);

SECItem* cert_FindDERCertBySubjectKeyID(SECItem* subjKeyID);

/* return maximum length of AVA value based on its type OID tag. */
extern int cert_AVAOidTagToMaxLen(SECOidTag tag);

/* Make an AVA, allocated from pool, from OID and DER encoded value */
extern CERTAVA* CERT_CreateAVAFromRaw(PLArenaPool* pool, const SECItem* OID,
                                      const SECItem* value);

/* Make an AVA from binary input specified by SECItem */
extern CERTAVA* CERT_CreateAVAFromSECItem(PLArenaPool* arena, SECOidTag kind,
                                          int valueType, SECItem* value);

/*
 * get a DPCache object for the given issuer subject and dp
 * Automatically creates the cache object if it doesn't exist yet.
 */
SECStatus AcquireDPCache(CERTCertificate* issuer, const SECItem* subject,
                         const SECItem* dp, PRTime t, void* wincx,
                         CRLDPCache** dpcache, PRBool* writeLocked);

/* check if a particular SN is in the CRL cache and return its entry */
dpcacheStatus DPCache_Lookup(CRLDPCache* cache, const SECItem* sn,
                             CERTCrlEntry** returned);

/* release a DPCache object that was previously acquired */
void ReleaseDPCache(CRLDPCache* dpcache, PRBool writeLocked);

/*
 * map Stan errors into NSS errors
 * This function examines the stan error stack and automatically sets
 * PORT_SetError(); to the appropriate SEC_ERROR value.
 */
void CERT_MapStanError();

/* Like CERT_VerifyCert, except with an additional argument, flags. The
 * flags are defined immediately below.
 */
SECStatus cert_VerifyCertWithFlags(CERTCertDBHandle* handle,
                                   CERTCertificate* cert, PRBool checkSig,
                                   SECCertUsage certUsage, PRTime t,
                                   PRUint32 flags, void* wincx,
                                   CERTVerifyLog* log);

/* Use the default settings.
 * cert_VerifyCertWithFlags(..., CERT_VERIFYCERT_USE_DEFAULTS, ...) is
 * equivalent to CERT_VerifyCert(...);
 */
#define CERT_VERIFYCERT_USE_DEFAULTS 0

/* Skip all the OCSP checks during certificate verification, regardless of
 * the global OCSP settings. By default, certificate |cert| will have its
 * revocation status checked via OCSP according to the global OCSP settings.
 *
 * OCSP checking is always skipped when certUsage is certUsageStatusResponder.
 */
#define CERT_VERIFYCERT_SKIP_OCSP 1

/* Interface function for libpkix cert validation engine:
 * cert_verify wrapper. */
SECStatus cert_VerifyCertChainPkix(CERTCertificate* cert, PRBool checkSig,
                                   SECCertUsage requiredUsage, PRTime time,
                                   void* wincx, CERTVerifyLog* log,
                                   PRBool* sigError, PRBool* revoked);

SECStatus cert_InitLocks(void);

SECStatus cert_DestroyLocks(void);

/*
 * fill in nsCertType field of the cert based on the cert extension
 */
extern SECStatus cert_GetCertType(CERTCertificate* cert);

/*
 * compute and return the value of nsCertType for cert, but do not
 * update the CERTCertificate.
 */
extern PRUint32 cert_ComputeCertType(CERTCertificate* cert);

extern PRBool cert_EKUAllowsIPsecIKE(CERTCertificate* cert,
                                     PRBool* isCritical);

void cert_AddToVerifyLog(CERTVerifyLog* log, CERTCertificate* cert,
                         long errorCode, unsigned int depth, void* arg);

/* Insert a DER CRL into the CRL cache, and take ownership of it.
 *
 * cert_CacheCRLByGeneralName takes ownership of the memory in crl argument
 * completely.  crl must be freeable by SECITEM_FreeItem. It will be freed
 * immediately if it is rejected from the CRL cache, or later during cache
 * updates when a new crl is available, or at shutdown time.
 *
 * canonicalizedName represents the source of the CRL, a GeneralName.
 * The format of the encoding is not restricted, but all callers of
 * cert_CacheCRLByGeneralName and cert_FindCRLByGeneralName must use
 * the same encoding. To facilitate X.500 name matching, a canonicalized
 * encoding of the GeneralName should be used, if available.
 */

SECStatus cert_CacheCRLByGeneralName(CERTCertDBHandle* dbhandle, SECItem* crl,
                                     const SECItem* canonicalizedName);

struct NamedCRLCacheStr {
    PRLock* lock;
    PLHashTable* entries;
};

/* NamedCRLCacheEntryStr is filled in by cert_CacheCRLByGeneralName,
 * and read by cert_FindCRLByGeneralName */
struct NamedCRLCacheEntryStr {
    SECItem* canonicalizedName;
    SECItem* crl; /* DER, kept only if CRL
                   * is successfully cached */
    PRBool inCRLCache;
    PRTime successfulInsertionTime; /* insertion time */
    PRTime lastAttemptTime;         /* time of last call to
                              cert_CacheCRLByGeneralName with this name */
    PRBool badDER;                  /* ASN.1 error */
    PRBool dupe;                    /* matching DER CRL already in CRL cache */
    PRBool unsupported;             /* IDP, delta, any other reason */
};

typedef enum {
    certRevocationStatusRevoked = 0,
    certRevocationStatusValid = 1,
    certRevocationStatusUnknown = 2
} CERTRevocationStatus;

/* Returns detailed status of the cert(revStatus variable). Tells if
 * issuer cache has OriginFetchedWithTimeout crl in it. */
SECStatus cert_CheckCertRevocationStatus(CERTCertificate* cert,
                                         CERTCertificate* issuer,
                                         const SECItem* dp, PRTime t,
                                         void* wincx,
                                         CERTRevocationStatus* revStatus,
                                         CERTCRLEntryReasonCode* revReason);

SECStatus cert_AcquireNamedCRLCache(NamedCRLCache** returned);

/* cert_FindCRLByGeneralName must be called only while the named cache is
 * acquired, and the entry is only valid until cache is released.
 */
SECStatus cert_FindCRLByGeneralName(NamedCRLCache* ncc,
                                    const SECItem* canonicalizedName,
                                    NamedCRLCacheEntry** retEntry);

SECStatus cert_ReleaseNamedCRLCache(NamedCRLCache* ncc);

/* This is private for now.  Maybe shoule be public. */
CERTGeneralName* cert_GetSubjectAltNameList(const CERTCertificate* cert,
                                            PLArenaPool* arena);

/* Count DNS names and IP addresses in a list of GeneralNames */
PRUint32 cert_CountDNSPatterns(CERTGeneralName* firstName);

/*
 * returns the trust status of the leaf certificate based on usage.
 * If the leaf is explicitly untrusted, this function will fail and
 * failedFlags will be set to the trust bit value that lead to the failure.
 * If the leaf is trusted, isTrusted is set to true and the function returns
 * SECSuccess. This function does not check if the cert is fit for a
 * particular usage.
 */
SECStatus cert_CheckLeafTrust(CERTCertificate* cert, SECCertUsage usage,
                              unsigned int* failedFlags, PRBool* isTrusted);

/*
 * Acquire the cert temp/perm lock
 */
void CERT_LockCertTempPerm(const CERTCertificate* cert);

/*
 * Release the temp/perm lock
 */
void CERT_UnlockCertTempPerm(const CERTCertificate* cert);

/*
 * Acquire the cert trust lock
 * There is currently one global lock for all certs, but I'm putting a cert
 * arg here so that it will be easy to make it per-cert in the future if
 * that turns out to be necessary.
 */
void CERT_LockCertTrust(const CERTCertificate* cert);

/*
 * Release the cert trust lock
 */
void CERT_UnlockCertTrust(const CERTCertificate* cert);

#endif /* _CERTI_H_ */
