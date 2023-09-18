/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * Internal header file included only by files in pkcs11 dir, or in
 * pkcs11 specific client and server files.
 */
#ifndef _SECMODI_H_
#define _SECMODI_H_ 1

#include <stddef.h>

#include "pkcs11.h"
#include "nssilock.h"
#include "secoidt.h"
#include "secdert.h"
#include "certt.h"
#include "secmodt.h"
#include "keythi.h"

SEC_BEGIN_PROTOS

/* proto-types */
extern SECStatus SECMOD_DeletePermDB(SECMODModule *module);
extern SECStatus SECMOD_AddPermDB(SECMODModule *module);
extern SECStatus SECMOD_Shutdown(void);
void nss_DumpModuleLog(void);

extern int secmod_PrivateModuleCount;

extern void SECMOD_Init(void);
SECStatus secmod_ModuleInit(SECMODModule *mod, SECMODModule **oldModule,
                            PRBool *alreadyLoaded);

/* list managment */
extern SECStatus SECMOD_AddModuleToList(SECMODModule *newModule);
extern SECStatus SECMOD_AddModuleToDBOnlyList(SECMODModule *newModule);
extern SECStatus SECMOD_AddModuleToUnloadList(SECMODModule *newModule);
extern void SECMOD_RemoveList(SECMODModuleList **, SECMODModuleList *);
extern void SECMOD_AddList(SECMODModuleList *, SECMODModuleList *, SECMODListLock *);
extern SECMODListLock *SECMOD_NewListLock(void);
extern void SECMOD_DestroyListLock(SECMODListLock *);
extern void SECMOD_GetWriteLock(SECMODListLock *);
extern void SECMOD_ReleaseWriteLock(SECMODListLock *);

/* Operate on modules by name */
extern SECMODModule *SECMOD_FindModuleByID(SECMODModuleID);
extern SECMODModule *secmod_FindModuleByFuncPtr(void *funcPtr);

/* database/memory management */
extern SECMODModuleList *SECMOD_NewModuleListElement(void);
extern SECMODModuleList *SECMOD_DestroyModuleListElement(SECMODModuleList *);
extern void SECMOD_DestroyModuleList(SECMODModuleList *);
extern SECStatus SECMOD_AddModule(SECMODModule *newModule);

extern unsigned long SECMOD_InternaltoPubCipherFlags(unsigned long internalFlags);

/* Library functions */
SECStatus secmod_LoadPKCS11Module(SECMODModule *, SECMODModule **oldModule);
SECStatus SECMOD_UnloadModule(SECMODModule *);
void SECMOD_SetInternalModule(SECMODModule *);
PRBool secmod_IsInternalKeySlot(SECMODModule *);
void secmod_SetInternalKeySlotFlag(SECMODModule *mod, PRBool val);

/* tools for checking if we are loading the same database twice */
typedef struct SECMODConfigListStr SECMODConfigList;
/* collect all the databases in a given spec */
SECMODConfigList *secmod_GetConfigList(PRBool isFIPS, char *spec, int *count);
/* see is a spec matches a database on the list */
PRBool secmod_MatchConfigList(const char *spec,
                              SECMODConfigList *conflist, int count);
/* returns the slot id from a module and modulespec */
CK_SLOT_ID secmod_GetSlotIDFromModuleSpec(const char *moduleSpec, SECMODModule *module);
/* free our list of databases */
void secmod_FreeConfigList(SECMODConfigList *conflist, int count);

/* parsing parameters */
/* returned char * must be freed by caller with PORT_Free */
/* children and ids are null terminated arrays which must be freed with
 * secmod_FreeChildren */
char *secmod_ParseModuleSpecForTokens(PRBool convert,
                                      PRBool isFIPS,
                                      const char *moduleSpec,
                                      char ***children,
                                      CK_SLOT_ID **ids);
void secmod_FreeChildren(char **children, CK_SLOT_ID *ids);
char *secmod_MkAppendTokensList(PLArenaPool *arena, char *origModuleSpec,
                                char *newModuleSpec, CK_SLOT_ID newID,
                                char **children, CK_SLOT_ID *ids);

void SECMOD_SlotDestroyModule(SECMODModule *module, PRBool fromSlot);
CK_RV pk11_notify(CK_SESSION_HANDLE session, CK_NOTIFICATION event,
                  CK_VOID_PTR pdata);
void pk11_SignedToUnsigned(CK_ATTRIBUTE *attrib);
CK_OBJECT_HANDLE pk11_FindObjectByTemplate(PK11SlotInfo *slot,
                                           CK_ATTRIBUTE *inTemplate, size_t tsize);
CK_OBJECT_HANDLE *pk11_FindObjectsByTemplate(PK11SlotInfo *slot,
                                             CK_ATTRIBUTE *inTemplate, size_t tsize, int *objCount);

#define PK11_GETTAB(x) ((CK_FUNCTION_LIST_3_0_PTR)((x)->functionList))
#define PK11_SETATTRS(x, id, v, l) \
    (x)->type = (id);              \
    (x)->pValue = (v);             \
    (x)->ulValueLen = (l);
SECStatus PK11_CreateNewObject(PK11SlotInfo *slot, CK_SESSION_HANDLE session,
                               const CK_ATTRIBUTE *theTemplate, int count,
                               PRBool token, CK_OBJECT_HANDLE *objectID);

SECStatus pbe_PK11AlgidToParam(SECAlgorithmID *algid, SECItem *mech);
SECStatus PBE_PK11ParamToAlgid(SECOidTag algTag, SECItem *param,
                               PLArenaPool *arena, SECAlgorithmID *algId);

PK11SymKey *pk11_TokenKeyGenWithFlagsAndKeyType(PK11SlotInfo *slot,
                                                CK_MECHANISM_TYPE type, SECItem *param, CK_KEY_TYPE keyType,
                                                int keySize, SECItem *keyId, CK_FLAGS opFlags,
                                                PK11AttrFlags attrFlags, void *wincx);

CK_MECHANISM_TYPE pk11_GetPBECryptoMechanism(SECAlgorithmID *algid,
                                             SECItem **param, SECItem *pwd, PRBool faulty3DES);

extern void pk11sdr_Init(void);
extern void pk11sdr_Shutdown(void);

/*
 * Private to pk11wrap.
 */

PRBool pk11_LoginStillRequired(PK11SlotInfo *slot, void *wincx);
CK_SESSION_HANDLE pk11_GetNewSession(PK11SlotInfo *slot, PRBool *owner);
void pk11_CloseSession(PK11SlotInfo *slot, CK_SESSION_HANDLE sess, PRBool own);
PK11SymKey *pk11_ForceSlot(PK11SymKey *symKey, CK_MECHANISM_TYPE type,
                           CK_ATTRIBUTE_TYPE operation);
/* Convert key operation flags to PKCS #11 attributes. */
unsigned int pk11_OpFlagsToAttributes(CK_FLAGS flags,
                                      CK_ATTRIBUTE *attrs, CK_BBOOL *ckTrue);
/* Check for bad (conflicting) attribute flags */
PRBool pk11_BadAttrFlags(PK11AttrFlags attrFlags);
/* Convert key attribute flags to PKCS #11 attributes. */
unsigned int pk11_AttrFlagsToAttributes(PK11AttrFlags attrFlags,
                                        CK_ATTRIBUTE *attrs, CK_BBOOL *ckTrue, CK_BBOOL *ckFalse);
PRBool pk11_FindAttrInTemplate(CK_ATTRIBUTE *attr, unsigned int numAttrs,
                               CK_ATTRIBUTE_TYPE target);

CK_MECHANISM_TYPE pk11_mapWrapKeyType(KeyType keyType);
PK11SymKey *pk11_KeyExchange(PK11SlotInfo *slot, CK_MECHANISM_TYPE type,
                             CK_ATTRIBUTE_TYPE operation, CK_FLAGS flags, PRBool isPerm,
                             PK11SymKey *symKey);

PRBool pk11_HandleTrustObject(PK11SlotInfo *slot, CERTCertificate *cert,
                              CERTCertTrust *trust);
CK_OBJECT_HANDLE pk11_FindPubKeyByAnyCert(CERTCertificate *cert,
                                          PK11SlotInfo **slot, void *wincx);
SECStatus pk11_AuthenticateUnfriendly(PK11SlotInfo *slot, PRBool loadCerts,
                                      void *wincx);
int PK11_NumberObjectsFor(PK11SlotInfo *slot, CK_ATTRIBUTE *findTemplate,
                          int templateCount);
SECItem *pk11_GetLowLevelKeyFromHandle(PK11SlotInfo *slot,
                                       CK_OBJECT_HANDLE handle);
SECStatus PK11_TraverseSlot(PK11SlotInfo *slot, void *arg);
CK_OBJECT_HANDLE pk11_FindPrivateKeyFromCertID(PK11SlotInfo *slot,
                                               SECItem *keyID);
SECKEYPrivateKey *PK11_MakePrivKey(PK11SlotInfo *slot, KeyType keyType,
                                   PRBool isTemp, CK_OBJECT_HANDLE privID, void *wincx);
CERTCertificate *PK11_MakeCertFromHandle(PK11SlotInfo *slot,
                                         CK_OBJECT_HANDLE certID, CK_ATTRIBUTE *privateLabel);

SECItem *pk11_GenerateNewParamWithKeyLen(CK_MECHANISM_TYPE type, int keyLen);
SECItem *pk11_ParamFromIVWithLen(CK_MECHANISM_TYPE type,
                                 SECItem *iv, int keyLen);
SECItem *pk11_mkcertKeyID(CERTCertificate *cert);

SEC_END_PROTOS

#endif
