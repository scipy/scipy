/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_policynode.h
 *
 * PolicyNode Type Definitions
 *
 */

#ifndef _PKIX_POLICYNODE_H
#define _PKIX_POLICYNODE_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This structure reflects the contents of a policy node...
 */
struct PKIX_PolicyNodeStruct {
    PKIX_PL_OID *validPolicy;
    PKIX_List *qualifierSet;    /* CertPolicyQualifiers */
    PKIX_Boolean criticality;
    PKIX_List *expectedPolicySet;       /* OIDs */
    PKIX_PolicyNode *parent;
    PKIX_List *children;                /* PolicyNodes */
    PKIX_UInt32 depth;
};

PKIX_Error *
pkix_SinglePolicyNode_ToString(
        PKIX_PolicyNode *node,
        PKIX_PL_String **pString,
        void *plContext);

PKIX_Error *
pkix_PolicyNode_GetChildrenMutable(
        PKIX_PolicyNode *node,
        PKIX_List **pChildren,  /* PolicyNodes */
        void *plContext);

PKIX_Error *
pkix_PolicyNode_Create(
        PKIX_PL_OID *validPolicy,
        PKIX_List *qualifierSet,        /* CertPolicyQualifiers */
        PKIX_Boolean criticality,
        PKIX_List *expectedPolicySet,   /* OIDs */
        PKIX_PolicyNode **pObject,
        void *plContext);

PKIX_Error *
pkix_PolicyNode_AddToParent(
        PKIX_PolicyNode *parentNode,
        PKIX_PolicyNode *child,
        void *plContext);

PKIX_Error *
pkix_PolicyNode_Prune(
        PKIX_PolicyNode *node,
        PKIX_UInt32 depth,
        PKIX_Boolean *pDelete,
        void *plContext);

PKIX_Error *
pkix_PolicyNode_RegisterSelf(
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_POLICYNODE_H */
