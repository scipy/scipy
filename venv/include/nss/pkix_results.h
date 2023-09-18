/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines functions associated with the results used
 * by the top-level functions.
 *
 */

#ifndef _PKIX_RESULTS_H
#define _PKIX_RESULTS_H

#include "pkixt.h"

#ifdef __cplusplus
extern "C" {
#endif

/* General
 *
 * Please refer to the libpkix Programmer's Guide for detailed information
 * about how to use the libpkix library. Certain key warnings and notices from
 * that document are repeated here for emphasis.
 *
 * All identifiers in this file (and all public identifiers defined in
 * libpkix) begin with "PKIX_". Private identifiers only intended for use
 * within the library begin with "pkix_".
 *
 * A function returns NULL upon success, and a PKIX_Error pointer upon failure.
 *
 * Unless otherwise noted, for all accessor (gettor) functions that return a
 * PKIX_PL_Object pointer, callers should assume that this pointer refers to a
 * shared object. Therefore, the caller should treat this shared object as
 * read-only and should not modify this shared object. When done using the
 * shared object, the caller should release the reference to the object by
 * using the PKIX_PL_Object_DecRef function.
 *
 * While a function is executing, if its arguments (or anything referred to by
 * its arguments) are modified, free'd, or destroyed, the function's behavior
 * is undefined.
 *
 */
/* PKIX_ValidateResult
 *
 * PKIX_ValidateResult represents the result of a PKIX_ValidateChain call. It
 * consists of the valid policy tree and public key resulting from validation,
 * as well as the trust anchor used for this chain. Once created, a
 * ValidateResult object is immutable.
 */

/*
 * FUNCTION: PKIX_ValidateResult_GetPolicyTree
 * DESCRIPTION:
 *
 *  Retrieves the PolicyNode component (representing the valid_policy_tree)
 *  from the ValidateResult object pointed to by "result" and stores it at
 *  "pPolicyTree".
 *
 * PARAMETERS:
 *  "result"
 *      Address of ValidateResult whose policy tree is to be stored. Must be
 *      non-NULL.
 *  "pPolicyTree"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateResult_GetPolicyTree(
        PKIX_ValidateResult *result,
        PKIX_PolicyNode **pPolicyTree,
        void *plContext);

/*
 * FUNCTION: PKIX_ValidateResult_GetPublicKey
 * DESCRIPTION:
 *
 *  Retrieves the PublicKey component (representing the valid public_key) of
 *  the ValidateResult object pointed to by "result" and stores it at
 *  "pPublicKey".
 *
 * PARAMETERS:
 *  "result"
 *      Address of ValidateResult whose public key is to be stored.
 *      Must be non-NULL.
 *  "pPublicKey"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateResult_GetPublicKey(
        PKIX_ValidateResult *result,
        PKIX_PL_PublicKey **pPublicKey,
        void *plContext);

/*
 * FUNCTION: PKIX_ValidateResult_GetTrustAnchor
 * DESCRIPTION:
 *
 *  Retrieves the TrustAnchor component (representing the trust anchor used
 *  during chain validation) of the ValidateResult object pointed to by
 *  "result" and stores it at "pTrustAnchor".
 *
 * PARAMETERS:
 *  "result"
 *      Address of ValidateResult whose trust anchor is to be stored.
 *      Must be non-NULL.
 *  "pTrustAnchor"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_ValidateResult_GetTrustAnchor(
        PKIX_ValidateResult *result,
        PKIX_TrustAnchor **pTrustAnchor,
        void *plContext);

/* PKIX_BuildResult
 *
 * PKIX_BuildResult represents the result of a PKIX_BuildChain call. It
 * consists of a ValidateResult object, as well as the built and validated
 * CertChain. Once created, a BuildResult object is immutable.
 */

/*
 * FUNCTION: PKIX_BuildResult_GetValidateResult
 * DESCRIPTION:
 *
 *  Retrieves the ValidateResult component (representing the build's validate
 *  result) of the BuildResult object pointed to by "result" and stores it at
 *  "pResult".
 *
 * PARAMETERS:
 *  "result"
 *      Address of BuildResult whose ValidateResult component is to be stored.
 *      Must be non-NULL.
 *  "pResult"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_BuildResult_GetValidateResult(
        PKIX_BuildResult *result,
        PKIX_ValidateResult **pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_BuildResult_GetCertChain
 * DESCRIPTION:
 *
 *  Retrieves the List of Certs (certChain) component (representing the built
 *  and validated CertChain) of the BuildResult object pointed to by "result"
 *  and stores it at "pChain".
 *
 * PARAMETERS:
 *  "result"
 *      Address of BuildResult whose CertChain component is to be stored.
 *      Must be non-NULL.
 *  "pChain"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_BuildResult_GetCertChain(
        PKIX_BuildResult *result,
        PKIX_List **pChain,
        void *plContext);

/* PKIX_PolicyNode
 *
 * PKIX_PolicyNode represents a node in the policy tree returned in
 * ValidateResult. The policy tree is the same length as the validated
 * certificate chain and the nodes are associated with a particular depth
 * (corresponding to a particular certificate in the chain).
 * PKIX_ValidateResult_GetPolicyTree returns the root node of the valid policy
 * tree. Other nodes can be accessed using the getChildren and getParents
 * functions, and individual elements of a node can be accessed with the
 * appropriate gettors. Once created, a PolicyNode is immutable.
 */

/*
 * FUNCTION: PKIX_PolicyNode_GetChildren
 * DESCRIPTION:
 *
 *  Retrieves the List of PolicyNodes representing the child nodes of the
 *  Policy Node pointed to by "node" and stores it at "pChildren". If "node"
 *  has no child nodes, this function stores an empty List at "pChildren".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose child nodes are to be stored.
 *      Must be non-NULL.
 *  "pChildren"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_GetChildren(
        PKIX_PolicyNode *node,
        PKIX_List **pChildren,  /* list of PKIX_PolicyNode */
        void *plContext);

/*
 * FUNCTION: PKIX_PolicyNode_GetParent
 * DESCRIPTION:
 *
 *  Retrieves the PolicyNode representing the parent node of the PolicyNode
 *  pointed to by "node" and stores it at "pParent". If "node" has no parent
 *  node, this function stores NULL at "pParent".
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose parent node is to be stored.
 *      Must be non-NULL.
 *  "pParent"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_GetParent(
        PKIX_PolicyNode *node,
        PKIX_PolicyNode **pParent,
        void *plContext);

/*
 * FUNCTION: PKIX_PolicyNode_GetValidPolicy
 * DESCRIPTION:
 *
 *  Retrieves the OID representing the valid policy of the PolicyNode pointed
 *  to by "node" and stores it at "pValidPolicy".
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose valid policy is to be stored.
 *      Must be non-NULL.
 *  "pValidPolicy"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_GetValidPolicy(
        PKIX_PolicyNode *node,
        PKIX_PL_OID **pValidPolicy,
        void *plContext);

/*
 * FUNCTION: PKIX_PolicyNode_GetPolicyQualifiers
 * DESCRIPTION:
 *
 *  Retrieves the List of CertPolicyQualifiers representing the policy
 *  qualifiers associated with the PolicyNode pointed to by "node" and stores
 *  it at "pQualifiers". If "node" has no policy qualifiers, this function
 *  stores an empty List at "pQualifiers".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose policy qualifiers are to be stored.
 *      Must be non-NULL.
 *  "pQualifiers"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_GetPolicyQualifiers(
        PKIX_PolicyNode *node,
        PKIX_List **pQualifiers,  /* list of PKIX_PL_CertPolicyQualifier */
        void *plContext);

/*
 * FUNCTION: PKIX_PolicyNode_GetExpectedPolicies
 * DESCRIPTION:
 *
 *  Retrieves the List of OIDs representing the expected policies associated
 *  with the PolicyNode pointed to by "node" and stores it at "pExpPolicies".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose expected policies are to be stored.
 *      Must be non-NULL.
 *  "pExpPolicies"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_GetExpectedPolicies(
        PKIX_PolicyNode *node,
        PKIX_List **pExpPolicies,  /* list of PKIX_PL_OID */
        void *plContext);

/*
 * FUNCTION: PKIX_PolicyNode_IsCritical
 * DESCRIPTION:
 *
 *  Checks the criticality field of the PolicyNode pointed to by "node" and
 *  stores the Boolean result at "pCritical".
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose criticality field is examined.
 *      Must be non-NULL.
 *  "pCritical"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_IsCritical(
        PKIX_PolicyNode *node,
        PKIX_Boolean *pCritical,
        void *plContext);

/*
 * FUNCTION: PKIX_PolicyNode_GetDepth
 * DESCRIPTION:
 *
 *  Retrieves the depth component of the PolicyNode pointed to by "node" and
 *  stores it at "pDepth".
 *
 * PARAMETERS:
 *  "node"
 *      Address of PolicyNode whose depth component is to be stored.
 *      Must be non-NULL.
 *  "pDepth"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Result Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PolicyNode_GetDepth(
        PKIX_PolicyNode *node,
        PKIX_UInt32 *pDepth,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_RESULTS_H */
