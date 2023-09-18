/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _GENAME_H_
#define _GENAME_H_

#include "plarena.h"
#include "seccomon.h"
#include "secoidt.h"
#include "secasn1.h"
#include "secder.h"
#include "certt.h"

/************************************************************************/
SEC_BEGIN_PROTOS

extern const SEC_ASN1Template CERT_GeneralNamesTemplate[];

extern SECItem **cert_EncodeGeneralNames(PLArenaPool *arena,
                                         CERTGeneralName *names);

extern CERTGeneralName *cert_DecodeGeneralNames(PLArenaPool *arena,
                                                SECItem **encodedGenName);

extern SECStatus cert_DestroyGeneralNames(CERTGeneralName *name);

extern SECStatus cert_EncodeNameConstraints(CERTNameConstraints *constraints,
                                            PLArenaPool *arena, SECItem *dest);

extern CERTNameConstraints *cert_DecodeNameConstraints(
    PLArenaPool *arena, const SECItem *encodedConstraints);

extern CERTGeneralName *cert_CombineNamesLists(CERTGeneralName *list1,
                                               CERTGeneralName *list2);

extern CERTNameConstraint *cert_CombineConstraintsLists(
    CERTNameConstraint *list1, CERTNameConstraint *list2);

/*********************************************************************/
/* A thread safe implementation of General Names                     */
/*********************************************************************/

/* Destroy a Single CERTGeneralName */
void CERT_DestroyGeneralName(CERTGeneralName *name);

SECStatus CERT_CompareGeneralName(CERTGeneralName *a, CERTGeneralName *b);

SECStatus CERT_CopyGeneralName(PLArenaPool *arena, CERTGeneralName *dest,
                               CERTGeneralName *src);

/* General Name Lists are a thread safe, reference counting layer to
 * general names */

/* Destroys a CERTGeneralNameList */
void CERT_DestroyGeneralNameList(CERTGeneralNameList *list);

/* Creates a CERTGeneralNameList */
CERTGeneralNameList *CERT_CreateGeneralNameList(CERTGeneralName *name);

/* Compares two CERTGeneralNameList */
SECStatus CERT_CompareGeneralNameLists(CERTGeneralNameList *a,
                                       CERTGeneralNameList *b);

/* returns a copy of the first name of the type requested */
void *CERT_GetGeneralNameFromListByType(CERTGeneralNameList *list,
                                        CERTGeneralNameType type,
                                        PLArenaPool *arena);

/* Adds a name to the tail of the list */
void CERT_AddGeneralNameToList(CERTGeneralNameList *list,
                               CERTGeneralNameType type, void *data,
                               SECItem *oid);

/* returns a duplicate of the CERTGeneralNameList */
CERTGeneralNameList *CERT_DupGeneralNameList(CERTGeneralNameList *list);

/* returns the number of CERTGeneralName objects in the  doubly linked
** list of which *names is a member.
*/
extern int CERT_GetNamesLength(CERTGeneralName *names);

/************************************************************************/

SECStatus CERT_CompareNameSpace(CERTCertificate *cert,
                                CERTGeneralName *namesList,
                                CERTCertificate **certsList,
                                PLArenaPool *reqArena,
                                CERTCertificate **pBadCert);

SEC_END_PROTOS

#endif
