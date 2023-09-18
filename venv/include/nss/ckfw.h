/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef CKFW_H
#define CKFW_H

/*
 * ckfw.h
 *
 * This file prototypes the private calls of the NSS Cryptoki Framework.
 */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

#ifndef NSSCKT_H
#include "nssckt.h"
#endif /* NSSCKT_H */

#ifndef NSSCKFWT_H
#include "nssckfwt.h"
#endif /* NSSCKFWT_H */

#ifndef NSSCKMDT_H
#include "nssckmdt.h"
#endif /* NSSCKMDT_H */

/*
 * NSSCKFWInstance
 *
 *  -- create/destroy --
 *  nssCKFWInstance_Create
 *  nssCKFWInstance_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWInstance_GetMDInstance
 *  nssCKFWInstance_GetArena
 *  nssCKFWInstance_MayCreatePthreads
 *  nssCKFWInstance_CreateMutex
 *  nssCKFWInstance_GetConfigurationData
 *  nssCKFWInstance_GetInitArgs
 *
 *  -- private accessors --
 *  nssCKFWInstance_CreateSessionHandle
 *  nssCKFWInstance_ResolveSessionHandle
 *  nssCKFWInstance_DestroySessionHandle
 *  nssCKFWInstance_FindSessionHandle
 *  nssCKFWInstance_CreateObjectHandle
 *  nssCKFWInstance_ResolveObjectHandle
 *  nssCKFWInstance_DestroyObjectHandle
 *  nssCKFWInstance_FindObjectHandle
 *
 *  -- module fronts --
 *  nssCKFWInstance_GetNSlots
 *  nssCKFWInstance_GetCryptokiVersion
 *  nssCKFWInstance_GetManufacturerID
 *  nssCKFWInstance_GetFlags
 *  nssCKFWInstance_GetLibraryDescription
 *  nssCKFWInstance_GetLibraryVersion
 *  nssCKFWInstance_GetModuleHandlesSessionObjects
 *  nssCKFWInstance_GetSlots
 *  nssCKFWInstance_WaitForSlotEvent
 *
 *  -- debugging versions only --
 *  nssCKFWInstance_verifyPointer
 */

/*
 * nssCKFWInstance_Create
 *
 */
NSS_EXTERN NSSCKFWInstance *
nssCKFWInstance_Create(
    CK_C_INITIALIZE_ARGS_PTR pInitArgs,
    CryptokiLockingState LockingState,
    NSSCKMDInstance *mdInstance,
    CK_RV *pError);

/*
 * nssCKFWInstance_Destroy
 *
 */
NSS_EXTERN CK_RV
nssCKFWInstance_Destroy(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetMDInstance
 *
 */
NSS_EXTERN NSSCKMDInstance *
nssCKFWInstance_GetMDInstance(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetArena
 *
 */
NSS_EXTERN NSSArena *
nssCKFWInstance_GetArena(
    NSSCKFWInstance *fwInstance,
    CK_RV *pError);

/*
 * nssCKFWInstance_MayCreatePthreads
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWInstance_MayCreatePthreads(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_CreateMutex
 *
 */
NSS_EXTERN NSSCKFWMutex *
nssCKFWInstance_CreateMutex(
    NSSCKFWInstance *fwInstance,
    NSSArena *arena,
    CK_RV *pError);

/*
 * nssCKFWInstance_GetConfigurationData
 *
 */
NSS_EXTERN NSSUTF8 *
nssCKFWInstance_GetConfigurationData(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetInitArgs
 *
 */
NSS_EXTERN CK_C_INITIALIZE_ARGS_PTR
nssCKFWInstance_GetInitArgs(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_CreateSessionHandle
 *
 */
NSS_EXTERN CK_SESSION_HANDLE
nssCKFWInstance_CreateSessionHandle(
    NSSCKFWInstance *fwInstance,
    NSSCKFWSession *fwSession,
    CK_RV *pError);

/*
 * nssCKFWInstance_ResolveSessionHandle
 *
 */
NSS_EXTERN NSSCKFWSession *
nssCKFWInstance_ResolveSessionHandle(
    NSSCKFWInstance *fwInstance,
    CK_SESSION_HANDLE hSession);

/*
 * nssCKFWInstance_DestroySessionHandle
 *
 */
NSS_EXTERN void
nssCKFWInstance_DestroySessionHandle(
    NSSCKFWInstance *fwInstance,
    CK_SESSION_HANDLE hSession);

/*
 * nssCKFWInstance_FindSessionHandle
 *
 */
NSS_EXTERN CK_SESSION_HANDLE
nssCKFWInstance_FindSessionHandle(
    NSSCKFWInstance *fwInstance,
    NSSCKFWSession *fwSession);

/*
 * nssCKFWInstance_CreateObjectHandle
 *
 */
NSS_EXTERN CK_OBJECT_HANDLE
nssCKFWInstance_CreateObjectHandle(
    NSSCKFWInstance *fwInstance,
    NSSCKFWObject *fwObject,
    CK_RV *pError);

/*
 * nssCKFWInstance_ResolveObjectHandle
 *
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWInstance_ResolveObjectHandle(
    NSSCKFWInstance *fwInstance,
    CK_OBJECT_HANDLE hObject);

/*
 * nssCKFWInstance_ReassignObjectHandle
 *
 */
NSS_EXTERN CK_RV
nssCKFWInstance_ReassignObjectHandle(
    NSSCKFWInstance *fwInstance,
    CK_OBJECT_HANDLE hObject,
    NSSCKFWObject *fwObject);

/*
 * nssCKFWInstance_DestroyObjectHandle
 *
 */
NSS_EXTERN void
nssCKFWInstance_DestroyObjectHandle(
    NSSCKFWInstance *fwInstance,
    CK_OBJECT_HANDLE hObject);

/*
 * nssCKFWInstance_FindObjectHandle
 *
 */
NSS_EXTERN CK_OBJECT_HANDLE
nssCKFWInstance_FindObjectHandle(
    NSSCKFWInstance *fwInstance,
    NSSCKFWObject *fwObject);

/*
 * nssCKFWInstance_GetNSlots
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWInstance_GetNSlots(
    NSSCKFWInstance *fwInstance,
    CK_RV *pError);

/*
 * nssCKFWInstance_GetCryptokiVersion
 *
 */
NSS_EXTERN CK_VERSION
nssCKFWInstance_GetCryptokiVersion(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetManufacturerID
 *
 */
NSS_EXTERN CK_RV
nssCKFWInstance_GetManufacturerID(
    NSSCKFWInstance *fwInstance,
    CK_CHAR manufacturerID[32]);

/*
 * nssCKFWInstance_GetFlags
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWInstance_GetFlags(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetLibraryDescription
 *
 */
NSS_EXTERN CK_RV
nssCKFWInstance_GetLibraryDescription(
    NSSCKFWInstance *fwInstance,
    CK_CHAR libraryDescription[32]);

/*
 * nssCKFWInstance_GetLibraryVersion
 *
 */
NSS_EXTERN CK_VERSION
nssCKFWInstance_GetLibraryVersion(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetModuleHandlesSessionObjects
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWInstance_GetModuleHandlesSessionObjects(
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWInstance_GetSlots
 *
 */
NSS_EXTERN NSSCKFWSlot **
nssCKFWInstance_GetSlots(
    NSSCKFWInstance *fwInstance,
    CK_RV *pError);

/*
 * nssCKFWInstance_WaitForSlotEvent
 *
 */
NSS_EXTERN NSSCKFWSlot *
nssCKFWInstance_WaitForSlotEvent(
    NSSCKFWInstance *fwInstance,
    CK_BBOOL block,
    CK_RV *pError);

/*
 * nssCKFWInstance_verifyPointer
 *
 */
NSS_EXTERN CK_RV
nssCKFWInstance_verifyPointer(
    const NSSCKFWInstance *fwInstance);

/*
 * NSSCKFWSlot
 *
 *  -- create/destroy --
 *  nssCKFWSlot_Create
 *  nssCKFWSlot_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWSlot_GetMDSlot
 *  nssCKFWSlot_GetFWInstance
 *  nssCKFWSlot_GetMDInstance
 *
 *  -- private accessors --
 *  nssCKFWSlot_GetSlotID
 *
 *  -- module fronts --
 *  nssCKFWSlot_GetSlotDescription
 *  nssCKFWSlot_GetManufacturerID
 *  nssCKFWSlot_GetTokenPresent
 *  nssCKFWSlot_GetRemovableDevice
 *  nssCKFWSlot_GetHardwareSlot
 *  nssCKFWSlot_GetHardwareVersion
 *  nssCKFWSlot_GetFirmwareVersion
 *  nssCKFWSlot_GetToken
 */

/*
 * nssCKFWSlot_Create
 *
 */
NSS_EXTERN NSSCKFWSlot *
nssCKFWSlot_Create(
    NSSCKFWInstance *fwInstance,
    NSSCKMDSlot *mdSlot,
    CK_SLOT_ID slotID,
    CK_RV *pError);

/*
 * nssCKFWSlot_Destroy
 *
 */
NSS_EXTERN CK_RV
nssCKFWSlot_Destroy(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetMDSlot
 *
 */
NSS_EXTERN NSSCKMDSlot *
nssCKFWSlot_GetMDSlot(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetFWInstance
 *
 */

NSS_EXTERN NSSCKFWInstance *
nssCKFWSlot_GetFWInstance(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetMDInstance
 *
 */

NSS_EXTERN NSSCKMDInstance *
nssCKFWSlot_GetMDInstance(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetSlotID
 *
 */
NSS_EXTERN CK_SLOT_ID
nssCKFWSlot_GetSlotID(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetSlotDescription
 *
 */
NSS_EXTERN CK_RV
nssCKFWSlot_GetSlotDescription(
    NSSCKFWSlot *fwSlot,
    CK_CHAR slotDescription[64]);

/*
 * nssCKFWSlot_GetManufacturerID
 *
 */
NSS_EXTERN CK_RV
nssCKFWSlot_GetManufacturerID(
    NSSCKFWSlot *fwSlot,
    CK_CHAR manufacturerID[32]);

/*
 * nssCKFWSlot_GetTokenPresent
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWSlot_GetTokenPresent(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetRemovableDevice
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWSlot_GetRemovableDevice(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetHardwareSlot
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWSlot_GetHardwareSlot(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetHardwareVersion
 *
 */
NSS_EXTERN CK_VERSION
nssCKFWSlot_GetHardwareVersion(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetFirmwareVersion
 *
 */
NSS_EXTERN CK_VERSION
nssCKFWSlot_GetFirmwareVersion(
    NSSCKFWSlot *fwSlot);

/*
 * nssCKFWSlot_GetToken
 *
 */
NSS_EXTERN NSSCKFWToken *
nssCKFWSlot_GetToken(
    NSSCKFWSlot *fwSlot,
    CK_RV *pError);

/*
 * nssCKFWSlot_ClearToken
 *
 */
NSS_EXTERN void
nssCKFWSlot_ClearToken(
    NSSCKFWSlot *fwSlot);

/*
 * NSSCKFWToken
 *
 *  -- create/destroy --
 *  nssCKFWToken_Create
 *  nssCKFWToken_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWToken_GetMDToken
 *  nssCKFWToken_GetFWSlot
 *  nssCKFWToken_GetMDSlot
 *  nssCKFWToken_GetSessionState
 *
 *  -- private accessors --
 *  nssCKFWToken_SetSessionState
 *  nssCKFWToken_RemoveSession
 *  nssCKFWToken_CloseAllSessions
 *  nssCKFWToken_GetSessionCount
 *  nssCKFWToken_GetRwSessionCount
 *  nssCKFWToken_GetRoSessionCount
 *  nssCKFWToken_GetSessionObjectHash
 *  nssCKFWToken_GetMDObjectHash
 *  nssCKFWToken_GetObjectHandleHash
 *
 *  -- module fronts --
 *  nssCKFWToken_InitToken
 *  nssCKFWToken_GetLabel
 *  nssCKFWToken_GetManufacturerID
 *  nssCKFWToken_GetModel
 *  nssCKFWToken_GetSerialNumber
 *  nssCKFWToken_GetHasRNG
 *  nssCKFWToken_GetIsWriteProtected
 *  nssCKFWToken_GetLoginRequired
 *  nssCKFWToken_GetUserPinInitialized
 *  nssCKFWToken_GetRestoreKeyNotNeeded
 *  nssCKFWToken_GetHasClockOnToken
 *  nssCKFWToken_GetHasProtectedAuthenticationPath
 *  nssCKFWToken_GetSupportsDualCryptoOperations
 *  nssCKFWToken_GetMaxSessionCount
 *  nssCKFWToken_GetMaxRwSessionCount
 *  nssCKFWToken_GetMaxPinLen
 *  nssCKFWToken_GetMinPinLen
 *  nssCKFWToken_GetTotalPublicMemory
 *  nssCKFWToken_GetFreePublicMemory
 *  nssCKFWToken_GetTotalPrivateMemory
 *  nssCKFWToken_GetFreePrivateMemory
 *  nssCKFWToken_GetHardwareVersion
 *  nssCKFWToken_GetFirmwareVersion
 *  nssCKFWToken_GetUTCTime
 *  nssCKFWToken_OpenSession
 *  nssCKFWToken_GetMechanismCount
 *  nssCKFWToken_GetMechanismTypes
 *  nssCKFWToken_GetMechanism
 */

/*
 * nssCKFWToken_Create
 *
 */
NSS_EXTERN NSSCKFWToken *
nssCKFWToken_Create(
    NSSCKFWSlot *fwSlot,
    NSSCKMDToken *mdToken,
    CK_RV *pError);

/*
 * nssCKFWToken_Destroy
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_Destroy(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMDToken
 *
 */
NSS_EXTERN NSSCKMDToken *
nssCKFWToken_GetMDToken(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetArena
 *
 */
NSS_EXTERN NSSArena *
nssCKFWToken_GetArena(
    NSSCKFWToken *fwToken,
    CK_RV *pError);

/*
 * nssCKFWToken_GetFWSlot
 *
 */
NSS_EXTERN NSSCKFWSlot *
nssCKFWToken_GetFWSlot(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMDSlot
 *
 */
NSS_EXTERN NSSCKMDSlot *
nssCKFWToken_GetMDSlot(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetSessionState
 *
 */
NSS_EXTERN CK_STATE
nssCKFWToken_GetSessionState(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_InitToken
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_InitToken(
    NSSCKFWToken *fwToken,
    NSSItem *pin,
    NSSUTF8 *label);

/*
 * nssCKFWToken_GetLabel
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_GetLabel(
    NSSCKFWToken *fwToken,
    CK_CHAR label[32]);

/*
 * nssCKFWToken_GetManufacturerID
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_GetManufacturerID(
    NSSCKFWToken *fwToken,
    CK_CHAR manufacturerID[32]);

/*
 * nssCKFWToken_GetModel
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_GetModel(
    NSSCKFWToken *fwToken,
    CK_CHAR model[16]);

/*
 * nssCKFWToken_GetSerialNumber
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_GetSerialNumber(
    NSSCKFWToken *fwToken,
    CK_CHAR serialNumber[16]);

/*
 * nssCKFWToken_GetHasRNG
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetHasRNG(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetIsWriteProtected
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetIsWriteProtected(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetLoginRequired
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetLoginRequired(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetUserPinInitialized
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetUserPinInitialized(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetRestoreKeyNotNeeded
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetRestoreKeyNotNeeded(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetHasClockOnToken
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetHasClockOnToken(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetHasProtectedAuthenticationPath
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetHasProtectedAuthenticationPath(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetSupportsDualCryptoOperations
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWToken_GetSupportsDualCryptoOperations(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMaxSessionCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetMaxSessionCount(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMaxRwSessionCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetMaxRwSessionCount(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMaxPinLen
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetMaxPinLen(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMinPinLen
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetMinPinLen(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetTotalPublicMemory
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetTotalPublicMemory(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetFreePublicMemory
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetFreePublicMemory(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetTotalPrivateMemory
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetTotalPrivateMemory(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetFreePrivateMemory
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetFreePrivateMemory(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetHardwareVersion
 *
 */
NSS_EXTERN CK_VERSION
nssCKFWToken_GetHardwareVersion(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetFirmwareVersion
 *
 */
NSS_EXTERN CK_VERSION
nssCKFWToken_GetFirmwareVersion(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetUTCTime
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_GetUTCTime(
    NSSCKFWToken *fwToken,
    CK_CHAR utcTime[16]);

/*
 * nssCKFWToken_OpenSession
 *
 */
NSS_EXTERN NSSCKFWSession *
nssCKFWToken_OpenSession(
    NSSCKFWToken *fwToken,
    CK_BBOOL rw,
    CK_VOID_PTR pApplication,
    CK_NOTIFY Notify,
    CK_RV *pError);

/*
 * nssCKFWToken_GetMechanismCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetMechanismCount(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMechanismTypes
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_GetMechanismTypes(
    NSSCKFWToken *fwToken,
    CK_MECHANISM_TYPE types[]);

/*
 * nssCKFWToken_GetMechanism
 *
 */
NSS_EXTERN NSSCKFWMechanism *
nssCKFWToken_GetMechanism(
    NSSCKFWToken *fwToken,
    CK_MECHANISM_TYPE which,
    CK_RV *pError);

/*
 * nssCKFWToken_SetSessionState
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_SetSessionState(
    NSSCKFWToken *fwToken,
    CK_STATE newState);

/*
 * nssCKFWToken_RemoveSession
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_RemoveSession(
    NSSCKFWToken *fwToken,
    NSSCKFWSession *fwSession);

/*
 * nssCKFWToken_CloseAllSessions
 *
 */
NSS_EXTERN CK_RV
nssCKFWToken_CloseAllSessions(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetSessionCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetSessionCount(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetRwSessionCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetRwSessionCount(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetRoSessionCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWToken_GetRoSessionCount(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetSessionObjectHash
 *
 */
NSS_EXTERN nssCKFWHash *
nssCKFWToken_GetSessionObjectHash(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetMDObjectHash
 *
 */
NSS_EXTERN nssCKFWHash *
nssCKFWToken_GetMDObjectHash(
    NSSCKFWToken *fwToken);

/*
 * nssCKFWToken_GetObjectHandleHash
 *
 */
NSS_EXTERN nssCKFWHash *
nssCKFWToken_GetObjectHandleHash(
    NSSCKFWToken *fwToken);

/*
 * NSSCKFWMechanism
 *
 *  -- create/destroy --
 *  nssCKFWMechanism_Create
 *  nssCKFWMechanism_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWMechanism_GetMDMechanism
 *
 *  -- private accessors --
 *
 *  -- module fronts --
 *  nssCKFWMechanism_GetMinKeySize
 *  nssCKFWMechanism_GetMaxKeySize
 *  nssCKFWMechanism_GetInHardware
 *  nssCKFWMechanism_GetCanEncrypt
 *  nssCKFWMechanism_GetCanDecrypt
 *  nssCKFWMechanism_GetCanDigest
 *  nssCKFWMechanism_GetCanSignRecover
 *  nssCKFWMechanism_GetCanVerify
 *  nssCKFWMechanism_GetCanVerifyRecover
 *  nssCKFWMechanism_GetCanGenerate
 *  nssCKFWMechanism_GetCanGenerateKeyPair
 *  nssCKFWMechanism_GetCanWrap
 *  nssCKFWMechanism_GetCanUnwrap
 *  nssCKFWMechanism_GetCanDerive
 *  nssCKFWMechanism_EncryptInit
 *  nssCKFWMechanism_DecryptInit
 *  nssCKFWMechanism_DigestInit
 *  nssCKFWMechanism_SignInit
 *  nssCKFWMechanism_SignRecoverInit
 *  nssCKFWMechanism_VerifyInit
 *  nssCKFWMechanism_VerifyRecoverInit
 *  nssCKFWMechanism_GenerateKey
 *  nssCKFWMechanism_GenerateKeyPair
 *  nssCKFWMechanism_GetWrapKeyLength
 *  nssCKFWMechanism_WrapKey
 *  nssCKFWMechanism_UnwrapKey
 *  nssCKFWMechanism_DeriveKey
 */

/*
 * nssCKFWMechanism_Create
 *
 */
NSS_EXTERN NSSCKFWMechanism *
nssCKFWMechanism_Create(
    NSSCKMDMechanism *mdMechanism,
    NSSCKMDToken *mdToken,
    NSSCKFWToken *fwToken,
    NSSCKMDInstance *mdInstance,
    NSSCKFWInstance *fwInstance);

/*
 * nssCKFWMechanism_Destroy
 *
 */
NSS_EXTERN void
nssCKFWMechanism_Destroy(
    NSSCKFWMechanism *fwMechanism);

/*
 * nssCKFWMechanism_GetMDMechanism
 *
 */

NSS_EXTERN NSSCKMDMechanism *
nssCKFWMechanism_GetMDMechanism(
    NSSCKFWMechanism *fwMechanism);

/*
 * nssCKFWMechanism_GetMinKeySize
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWMechanism_GetMinKeySize(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetMaxKeySize
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWMechanism_GetMaxKeySize(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetInHardware
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetInHardware(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * the following are determined automatically by which of the cryptographic
 * functions are defined for this mechanism.
 */
/*
 * nssCKFWMechanism_GetCanEncrypt
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanEncrypt(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanDecrypt
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanDecrypt(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanDigest
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanDigest(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanSign
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanSign(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanSignRecover
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanSignRecover(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanVerify
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanVerify(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanVerifyRecover
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanVerifyRecover(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanGenerate
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanGenerate(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanGenerateKeyPair
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanGenerateKeyPair(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanWrap
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanWrap(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanUnwrap
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanUnwrap(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GetCanDerive
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWMechanism_GetCanDerive(
    NSSCKFWMechanism *fwMechanism,
    CK_RV *pError);

/*
 *  nssCKFWMechanism_EncryptInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_EncryptInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 *  nssCKFWMechanism_DecryptInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_DecryptInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 *  nssCKFWMechanism_DigestInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_DigestInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession);

/*
 *  nssCKFWMechanism_SignInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_SignInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 *  nssCKFWMechanism_SignRecoverInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_SignRecoverInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 *  nssCKFWMechanism_VerifyInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_VerifyInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 *  nssCKFWMechanism_VerifyRecoverInit
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_VerifyRecoverInit(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM *pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 * nssCKFWMechanism_GenerateKey
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWMechanism_GenerateKey(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM_PTR pMechanism,
    NSSCKFWSession *fwSession,
    CK_ATTRIBUTE_PTR pTemplate,
    CK_ULONG ulAttributeCount,
    CK_RV *pError);

/*
 * nssCKFWMechanism_GenerateKeyPair
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_GenerateKeyPair(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM_PTR pMechanism,
    NSSCKFWSession *fwSession,
    CK_ATTRIBUTE_PTR pPublicKeyTemplate,
    CK_ULONG ulPublicKeyAttributeCount,
    CK_ATTRIBUTE_PTR pPrivateKeyTemplate,
    CK_ULONG ulPrivateKeyAttributeCount,
    NSSCKFWObject **fwPublicKeyObject,
    NSSCKFWObject **fwPrivateKeyObject);

/*
 * nssCKFWMechanism_GetWrapKeyLength
 */
NSS_EXTERN CK_ULONG
nssCKFWMechanism_GetWrapKeyLength(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM_PTR pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwWrappingKeyObject,
    NSSCKFWObject *fwObject,
    CK_RV *pError);

/*
 * nssCKFWMechanism_WrapKey
 */
NSS_EXTERN CK_RV
nssCKFWMechanism_WrapKey(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM_PTR pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwWrappingKeyObject,
    NSSCKFWObject *fwObject,
    NSSItem *wrappedKey);

/*
 * nssCKFWMechanism_UnwrapKey
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWMechanism_UnwrapKey(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM_PTR pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwWrappingKeyObject,
    NSSItem *wrappedKey,
    CK_ATTRIBUTE_PTR pTemplate,
    CK_ULONG ulAttributeCount,
    CK_RV *pError);

/*
 * nssCKFWMechanism_DeriveKey
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWMechanism_DeriveKey(
    NSSCKFWMechanism *fwMechanism,
    CK_MECHANISM_PTR pMechanism,
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwBaseKeyObject,
    CK_ATTRIBUTE_PTR pTemplate,
    CK_ULONG ulAttributeCount,
    CK_RV *pError);

/*
 * NSSCKFWCryptoOperation
 *
 *  -- create/destroy --
 *  nssCKFWCryptoOperation_Create
 *  nssCKFWCryptoOperation_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWCryptoOperation_GetMDCryptoOperation
 *  nssCKFWCryptoOperation_GetType
 *
 *  -- private accessors --
 *
 *  -- module fronts --
 * nssCKFWCryptoOperation_GetFinalLength
 * nssCKFWCryptoOperation_GetOperationLength
 * nssCKFWCryptoOperation_Final
 * nssCKFWCryptoOperation_Update
 * nssCKFWCryptoOperation_DigestUpdate
 * nssCKFWCryptoOperation_DigestKey
 * nssCKFWCryptoOperation_UpdateFinal
 */

/*
 *  nssCKFWCrytoOperation_Create
 */
NSS_EXTERN NSSCKFWCryptoOperation *
nssCKFWCryptoOperation_Create(
    NSSCKMDCryptoOperation *mdOperation,
    NSSCKMDSession *mdSession,
    NSSCKFWSession *fwSession,
    NSSCKMDToken *mdToken,
    NSSCKFWToken *fwToken,
    NSSCKMDInstance *mdInstance,
    NSSCKFWInstance *fwInstance,
    NSSCKFWCryptoOperationType type,
    CK_RV *pError);

/*
 *  nssCKFWCryptoOperation_Destroy
 */
NSS_EXTERN void
nssCKFWCryptoOperation_Destroy(
    NSSCKFWCryptoOperation *fwOperation);

/*
 *  nssCKFWCryptoOperation_GetMDCryptoOperation
 */
NSS_EXTERN NSSCKMDCryptoOperation *
nssCKFWCryptoOperation_GetMDCryptoOperation(
    NSSCKFWCryptoOperation *fwOperation);

/*
 *  nssCKFWCryptoOperation_GetType
 */
NSS_EXTERN NSSCKFWCryptoOperationType
nssCKFWCryptoOperation_GetType(
    NSSCKFWCryptoOperation *fwOperation);

/*
 * nssCKFWCryptoOperation_GetFinalLength
 */
NSS_EXTERN CK_ULONG
nssCKFWCryptoOperation_GetFinalLength(
    NSSCKFWCryptoOperation *fwOperation,
    CK_RV *pError);

/*
 * nssCKFWCryptoOperation_GetOperationLength
 */
NSS_EXTERN CK_ULONG
nssCKFWCryptoOperation_GetOperationLength(
    NSSCKFWCryptoOperation *fwOperation,
    NSSItem *inputBuffer,
    CK_RV *pError);

/*
 * nssCKFWCryptoOperation_Final
 */
NSS_EXTERN CK_RV
nssCKFWCryptoOperation_Final(
    NSSCKFWCryptoOperation *fwOperation,
    NSSItem *outputBuffer);

/*
 * nssCKFWCryptoOperation_Update
 */
NSS_EXTERN CK_RV
nssCKFWCryptoOperation_Update(
    NSSCKFWCryptoOperation *fwOperation,
    NSSItem *inputBuffer,
    NSSItem *outputBuffer);

/*
 * nssCKFWCryptoOperation_DigestUpdate
 */
NSS_EXTERN CK_RV
nssCKFWCryptoOperation_DigestUpdate(
    NSSCKFWCryptoOperation *fwOperation,
    NSSItem *inputBuffer);

/*
 * nssCKFWCryptoOperation_DigestKey
 */
NSS_EXTERN CK_RV
nssCKFWCryptoOperation_DigestKey(
    NSSCKFWCryptoOperation *fwOperation,
    NSSCKFWObject *fwKey);

/*
 * nssCKFWCryptoOperation_UpdateFinal
 */
NSS_EXTERN CK_RV
nssCKFWCryptoOperation_UpdateFinal(
    NSSCKFWCryptoOperation *fwOperation,
    NSSItem *inputBuffer,
    NSSItem *outputBuffer);

/*
 * nssCKFWCryptoOperation_UpdateCombo
 */
NSS_EXTERN CK_RV
nssCKFWCryptoOperation_UpdateCombo(
    NSSCKFWCryptoOperation *fwOperation,
    NSSCKFWCryptoOperation *fwPeerOperation,
    NSSItem *inputBuffer,
    NSSItem *outputBuffer);

/*
 * NSSCKFWSession
 *
 *  -- create/destroy --
 *  nssCKFWSession_Create
 *  nssCKFWSession_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWSession_GetMDSession
 *  nssCKFWSession_GetArena
 *  nssCKFWSession_CallNotification
 *  nssCKFWSession_IsRWSession
 *  nssCKFWSession_IsSO
 *  nssCKFWSession_GetCurrentCryptoOperation
 *
 *  -- private accessors --
 *  nssCKFWSession_GetFWSlot
 *  nssCKFWSession_GetSessionState
 *  nssCKFWSession_SetFWFindObjects
 *  nssCKFWSession_GetFWFindObjects
 *  nssCKFWSession_SetMDSession
 *  nssCKFWSession_SetHandle
 *  nssCKFWSession_GetHandle
 *  nssCKFWSession_RegisterSessionObject
 *  nssCKFWSession_DeregisterSessionObject
 *  nssCKFWSession_SetCurrentCryptoOperation
 *
 *  -- module fronts --
 *  nssCKFWSession_GetDeviceError
 *  nssCKFWSession_Login
 *  nssCKFWSession_Logout
 *  nssCKFWSession_InitPIN
 *  nssCKFWSession_SetPIN
 *  nssCKFWSession_GetOperationStateLen
 *  nssCKFWSession_GetOperationState
 *  nssCKFWSession_SetOperationState
 *  nssCKFWSession_CreateObject
 *  nssCKFWSession_CopyObject
 *  nssCKFWSession_FindObjectsInit
 *  nssCKFWSession_SeedRandom
 *  nssCKFWSession_GetRandom
 *  nssCKFWSession_Final
 *  nssCKFWSession_Update
 *  nssCKFWSession_DigestUpdate
 *  nssCKFWSession_DigestKey
 *  nssCKFWSession_UpdateFinal
 *  nssCKFWSession_UpdateCombo
 */

/*
 * nssCKFWSession_Create
 *
 */
NSS_EXTERN NSSCKFWSession *
nssCKFWSession_Create(
    NSSCKFWToken *fwToken,
    CK_BBOOL rw,
    CK_VOID_PTR pApplication,
    CK_NOTIFY Notify,
    CK_RV *pError);

/*
 * nssCKFWSession_Destroy
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_Destroy(
    NSSCKFWSession *fwSession,
    CK_BBOOL removeFromTokenHash);

/*
 * nssCKFWSession_GetMDSession
 *
 */
NSS_EXTERN NSSCKMDSession *
nssCKFWSession_GetMDSession(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_GetArena
 *
 */
NSS_EXTERN NSSArena *
nssCKFWSession_GetArena(
    NSSCKFWSession *fwSession,
    CK_RV *pError);

/*
 * nssCKFWSession_CallNotification
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_CallNotification(
    NSSCKFWSession *fwSession,
    CK_NOTIFICATION event);

/*
 * nssCKFWSession_IsRWSession
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWSession_IsRWSession(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_IsSO
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWSession_IsSO(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_GetFWSlot
 *
 */
NSS_EXTERN NSSCKFWSlot *
nssCKFWSession_GetFWSlot(
    NSSCKFWSession *fwSession);

/*
 * nssCFKWSession_GetSessionState
 *
 */
NSS_EXTERN CK_STATE
nssCKFWSession_GetSessionState(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_SetFWFindObjects
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_SetFWFindObjects(
    NSSCKFWSession *fwSession,
    NSSCKFWFindObjects *fwFindObjects);

/*
 * nssCKFWSession_GetFWFindObjects
 *
 */
NSS_EXTERN NSSCKFWFindObjects *
nssCKFWSession_GetFWFindObjects(
    NSSCKFWSession *fwSesssion,
    CK_RV *pError);

/*
 * nssCKFWSession_SetMDSession
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_SetMDSession(
    NSSCKFWSession *fwSession,
    NSSCKMDSession *mdSession);

/*
 * nssCKFWSession_SetHandle
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_SetHandle(
    NSSCKFWSession *fwSession,
    CK_SESSION_HANDLE hSession);

/*
 * nssCKFWSession_GetHandle
 *
 */
NSS_EXTERN CK_SESSION_HANDLE
nssCKFWSession_GetHandle(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_RegisterSessionObject
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_RegisterSessionObject(
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 * nssCKFWSession_DeregisterSessionObject
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_DeregisterSessionObject(
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwObject);

/*
 * nssCKFWSession_GetDeviceError
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWSession_GetDeviceError(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_Login
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_Login(
    NSSCKFWSession *fwSession,
    CK_USER_TYPE userType,
    NSSItem *pin);

/*
 * nssCKFWSession_Logout
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_Logout(
    NSSCKFWSession *fwSession);

/*
 * nssCKFWSession_InitPIN
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_InitPIN(
    NSSCKFWSession *fwSession,
    NSSItem *pin);

/*
 * nssCKFWSession_SetPIN
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_SetPIN(
    NSSCKFWSession *fwSession,
    const NSSItem *oldPin,
    NSSItem *newPin);

/*
 * nssCKFWSession_GetOperationStateLen
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWSession_GetOperationStateLen(
    NSSCKFWSession *fwSession,
    CK_RV *pError);

/*
 * nssCKFWSession_GetOperationState
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_GetOperationState(
    NSSCKFWSession *fwSession,
    NSSItem *buffer);

/*
 * nssCKFWSession_SetOperationState
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_SetOperationState(
    NSSCKFWSession *fwSession,
    NSSItem *state,
    NSSCKFWObject *encryptionKey,
    NSSCKFWObject *authenticationKey);

/*
 * nssCKFWSession_CreateObject
 *
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWSession_CreateObject(
    NSSCKFWSession *fwSession,
    CK_ATTRIBUTE_PTR pTemplate,
    CK_ULONG ulAttributeCount,
    CK_RV *pError);

/*
 * nssCKFWSession_CopyObject
 *
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWSession_CopyObject(
    NSSCKFWSession *fwSession,
    NSSCKFWObject *object,
    CK_ATTRIBUTE_PTR pTemplate,
    CK_ULONG ulAttributeCount,
    CK_RV *pError);

/*
 * nssCKFWSession_FindObjectsInit
 *
 */
NSS_EXTERN NSSCKFWFindObjects *
nssCKFWSession_FindObjectsInit(
    NSSCKFWSession *fwSession,
    CK_ATTRIBUTE_PTR pTemplate,
    CK_ULONG ulAttributeCount,
    CK_RV *pError);

/*
 * nssCKFWSession_SetCurrentCryptoOperation
 */
NSS_IMPLEMENT void
nssCKFWSession_SetCurrentCryptoOperation(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperation *fwOperation,
    NSSCKFWCryptoOperationState state);

/*
 * nssCKFWSession_GetCurrentCryptoOperation
 */
NSS_IMPLEMENT NSSCKFWCryptoOperation *
nssCKFWSession_GetCurrentCryptoOperation(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperationState state);

/*
 * nssCKFWSession_Final
 * (terminate a cryptographic operation and get the result)
 */
NSS_IMPLEMENT CK_RV
nssCKFWSession_Final(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperationType type,
    NSSCKFWCryptoOperationState state,
    CK_BYTE_PTR outBuf,
    CK_ULONG_PTR outBufLen);

/*
 * nssCKFWSession_Update
 * (get the next step of an encrypt/decrypt operation)
 */
NSS_IMPLEMENT CK_RV
nssCKFWSession_Update(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperationType type,
    NSSCKFWCryptoOperationState state,
    CK_BYTE_PTR inBuf,
    CK_ULONG inBufLen,
    CK_BYTE_PTR outBuf,
    CK_ULONG_PTR outBufLen);

/*
 * nssCKFWSession_DigestUpdate
 * (do the next step of an digest/sign/verify operation)
 */
NSS_IMPLEMENT CK_RV
nssCKFWSession_DigestUpdate(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperationType type,
    NSSCKFWCryptoOperationState state,
    CK_BYTE_PTR inBuf,
    CK_ULONG inBufLen);

/*
 * nssCKFWSession_DigestKey
 * (do the next step of an digest/sign/verify operation)
 */
NSS_IMPLEMENT CK_RV
nssCKFWSession_DigestKey(
    NSSCKFWSession *fwSession,
    NSSCKFWObject *fwKey);

/*
 * nssCKFWSession_UpdateFinal
 * (do a single-step of a cryptographic operation and get the result)
 */
NSS_IMPLEMENT CK_RV
nssCKFWSession_UpdateFinal(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperationType type,
    NSSCKFWCryptoOperationState state,
    CK_BYTE_PTR inBuf,
    CK_ULONG inBufLen,
    CK_BYTE_PTR outBuf,
    CK_ULONG_PTR outBufLen);

/*
 * nssCKFWSession_UpdateCombo
 * (do a combination encrypt/decrypt and sign/digest/verify operation)
 */
NSS_IMPLEMENT CK_RV
nssCKFWSession_UpdateCombo(
    NSSCKFWSession *fwSession,
    NSSCKFWCryptoOperationType encryptType,
    NSSCKFWCryptoOperationType digestType,
    NSSCKFWCryptoOperationState digestState,
    CK_BYTE_PTR inBuf,
    CK_ULONG inBufLen,
    CK_BYTE_PTR outBuf,
    CK_ULONG_PTR outBufLen);

/*
 * nssCKFWSession_SeedRandom
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_SeedRandom(
    NSSCKFWSession *fwSession,
    NSSItem *seed);

/*
 * nssCKFWSession_GetRandom
 *
 */
NSS_EXTERN CK_RV
nssCKFWSession_GetRandom(
    NSSCKFWSession *fwSession,
    NSSItem *buffer);

/*
 * NSSCKFWObject
 *
 * -- create/destroy --
 *  nssCKFWObject_Create
 *  nssCKFWObject_Finalize
 *  nssCKFWObject_Destroy
 *
 * -- implement public accessors --
 *  nssCKFWObject_GetMDObject
 *  nssCKFWObject_GetArena
 *
 * -- private accessors --
 *  nssCKFWObject_SetHandle
 *  nssCKFWObject_GetHandle
 *
 * -- module fronts --
 *  nssCKFWObject_IsTokenObject
 *  nssCKFWObject_GetAttributeCount
 *  nssCKFWObject_GetAttributeTypes
 *  nssCKFWObject_GetAttributeSize
 *  nssCKFWObject_GetAttribute
 *  nssCKFWObject_SetAttribute
 *  nssCKFWObject_GetObjectSize
 */

/*
 * nssCKFWObject_Create
 *
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWObject_Create(
    NSSArena *arena,
    NSSCKMDObject *mdObject,
    NSSCKFWSession *fwSession,
    NSSCKFWToken *fwToken,
    NSSCKFWInstance *fwInstance,
    CK_RV *pError);

/*
 * nssCKFWObject_Finalize
 *
 */
NSS_EXTERN void
nssCKFWObject_Finalize(
    NSSCKFWObject *fwObject,
    PRBool removeFromHash);

/*
 * nssCKFWObject_Destroy
 *
 */
NSS_EXTERN void
nssCKFWObject_Destroy(
    NSSCKFWObject *fwObject);

/*
 * nssCKFWObject_GetMDObject
 *
 */
NSS_EXTERN NSSCKMDObject *
nssCKFWObject_GetMDObject(
    NSSCKFWObject *fwObject);

/*
 * nssCKFWObject_GetArena
 *
 */
NSS_EXTERN NSSArena *
nssCKFWObject_GetArena(
    NSSCKFWObject *fwObject,
    CK_RV *pError);

/*
 * nssCKFWObject_SetHandle
 *
 */
NSS_EXTERN CK_RV
nssCKFWObject_SetHandle(
    NSSCKFWObject *fwObject,
    CK_OBJECT_HANDLE hObject);

/*
 * nssCKFWObject_GetHandle
 *
 */
NSS_EXTERN CK_OBJECT_HANDLE
nssCKFWObject_GetHandle(
    NSSCKFWObject *fwObject);

/*
 * nssCKFWObject_IsTokenObject
 *
 */
NSS_EXTERN CK_BBOOL
nssCKFWObject_IsTokenObject(
    NSSCKFWObject *fwObject);

/*
 * nssCKFWObject_GetAttributeCount
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWObject_GetAttributeCount(
    NSSCKFWObject *fwObject,
    CK_RV *pError);

/*
 * nssCKFWObject_GetAttributeTypes
 *
 */
NSS_EXTERN CK_RV
nssCKFWObject_GetAttributeTypes(
    NSSCKFWObject *fwObject,
    CK_ATTRIBUTE_TYPE_PTR typeArray,
    CK_ULONG ulCount);

/*
 * nssCKFWObject_GetAttributeSize
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWObject_GetAttributeSize(
    NSSCKFWObject *fwObject,
    CK_ATTRIBUTE_TYPE attribute,
    CK_RV *pError);

/*
 * nssCKFWObject_GetAttribute
 *
 * Usual NSS allocation rules:
 * If itemOpt is not NULL, it will be returned; otherwise an NSSItem
 * will be allocated.  If itemOpt is not NULL but itemOpt->data is,
 * the buffer will be allocated; otherwise, the buffer will be used.
 * Any allocations will come from the optional arena, if one is
 * specified.
 */
NSS_EXTERN NSSItem *
nssCKFWObject_GetAttribute(
    NSSCKFWObject *fwObject,
    CK_ATTRIBUTE_TYPE attribute,
    NSSItem *itemOpt,
    NSSArena *arenaOpt,
    CK_RV *pError);

/*
 * nssCKFWObject_SetAttribute
 *
 */
NSS_EXTERN CK_RV
nssCKFWObject_SetAttribute(
    NSSCKFWObject *fwObject,
    NSSCKFWSession *fwSession,
    CK_ATTRIBUTE_TYPE attribute,
    NSSItem *value);

/*
 * nssCKFWObject_GetObjectSize
 *
 */
NSS_EXTERN CK_ULONG
nssCKFWObject_GetObjectSize(
    NSSCKFWObject *fwObject,
    CK_RV *pError);

/*
 * NSSCKFWFindObjects
 *
 *  -- create/destroy --
 *  nssCKFWFindObjects_Create
 *  nssCKFWFindObjects_Destroy
 *
 *  -- implement public accessors --
 *  nssCKFWFindObjects_GetMDFindObjects
 *
 *  -- private accessors --
 *
 *  -- module fronts --
 *  nssCKFWFindObjects_Next
 */

/*
 * nssCKFWFindObjects_Create
 *
 */
NSS_EXTERN NSSCKFWFindObjects *
nssCKFWFindObjects_Create(
    NSSCKFWSession *fwSession,
    NSSCKFWToken *fwToken,
    NSSCKFWInstance *fwInstance,
    NSSCKMDFindObjects *mdFindObjects1,
    NSSCKMDFindObjects *mdFindObjects2,
    CK_RV *pError);

/*
 * nssCKFWFindObjects_Destroy
 *
 */
NSS_EXTERN void
nssCKFWFindObjects_Destroy(
    NSSCKFWFindObjects *fwFindObjects);

/*
 * nssCKFWFindObjects_GetMDFindObjects
 *
 */
NSS_EXTERN NSSCKMDFindObjects *
nssCKFWFindObjects_GetMDFindObjects(
    NSSCKFWFindObjects *fwFindObjects);

/*
 * nssCKFWFindObjects_Next
 *
 */
NSS_EXTERN NSSCKFWObject *
nssCKFWFindObjects_Next(
    NSSCKFWFindObjects *fwFindObjects,
    NSSArena *arenaOpt,
    CK_RV *pError);

/*
 * NSSCKFWMutex
 *
 *  nssCKFWMutex_Create
 *  nssCKFWMutex_Destroy
 *  nssCKFWMutex_Lock
 *  nssCKFWMutex_Unlock
 *
 */

/*
 * nssCKFWMutex_Create
 *
 */
NSS_EXTERN NSSCKFWMutex *
nssCKFWMutex_Create(
    CK_C_INITIALIZE_ARGS_PTR pInitArgs,
    CryptokiLockingState LockingState,
    NSSArena *arena,
    CK_RV *pError);

/*
 * nssCKFWMutex_Destroy
 *
 */
NSS_EXTERN CK_RV
nssCKFWMutex_Destroy(
    NSSCKFWMutex *mutex);

/*
 * nssCKFWMutex_Lock
 *
 */
NSS_EXTERN CK_RV
nssCKFWMutex_Lock(
    NSSCKFWMutex *mutex);

/*
 * nssCKFWMutex_Unlock
 *
 */
NSS_EXTERN CK_RV
nssCKFWMutex_Unlock(
    NSSCKFWMutex *mutex);

#endif /* CKFW_H */
