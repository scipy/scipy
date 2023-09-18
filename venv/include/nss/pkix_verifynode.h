/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_verifynode.h
 *
 * VerifyNode Type Definitions
 *
 */

#ifndef _PKIX_VERIFYNODE_H
#define _PKIX_VERIFYNODE_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This structure reflects the contents of a verify node...
 */
struct PKIX_VerifyNodeStruct {
    PKIX_PL_Cert *verifyCert;
    PKIX_List *children;                /* VerifyNodes */
    PKIX_UInt32 depth;
    PKIX_Error *error;
};

PKIX_Error *
pkix_SingleVerifyNode_ToString(
        PKIX_VerifyNode *node,
        PKIX_PL_String **pString,
        void *plContext);

PKIX_Error *
pkix_VerifyNode_Create(
        PKIX_PL_Cert *verifyCert,
        PKIX_UInt32 depth,
        PKIX_Error *error,
        PKIX_VerifyNode **pObject,
        void *plContext);

PKIX_Error *
pkix_VerifyNode_AddToChain(
        PKIX_VerifyNode *parentNode,
        PKIX_VerifyNode *child,
        void *plContext);

PKIX_Error *
pkix_VerifyNode_AddToTree(
        PKIX_VerifyNode *parentNode,
        PKIX_VerifyNode *child,
        void *plContext);

PKIX_Error *
pkix_VerifyNode_SetError(
        PKIX_VerifyNode *node,
        PKIX_Error *error,
        void *plContext);

PKIX_Error *
pkix_VerifyNode_RegisterSelf(
        void *plContext);

PKIX_Error *
pkix_VerifyNode_FindError(
        PKIX_VerifyNode *node,
        PKIX_Error **error,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_VERIFYNODE_H */
