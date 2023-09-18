/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_nameconstraints.h
 *
 * Name Constraints Object Definitions
 *
 */

#ifndef _PKIX_PL_NAMECONSTRAINTS_H
#define _PKIX_PL_NAMECONSTRAINTS_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_CertNameConstraintsStruct {
        PLArenaPool *arena;
        CERTNameConstraints **nssNameConstraintsList;
        PKIX_UInt32 numNssNameConstraints;
        PKIX_List *permittedList; /* list of PKIX_PL_GeneralName */
        PKIX_List *excludedList; /* list of PKIX_PL_GeneralName */
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_CertNameConstraints_RegisterSelf(void *plContext);

PKIX_Error *pkix_pl_CertNameConstraints_Create(
        CERTCertificate *nssCert,
        PKIX_PL_CertNameConstraints **pNameConstraints,
        void *plContext);

PKIX_Error *
pkix_pl_CertNameConstraints_CreateWithNames(
        PKIX_List *names, /* List of PKIX_PL_GeneralName */
        PKIX_PL_CertNameConstraints **pNameConstraints,
        void *plContext);

PKIX_Error *
pkix_pl_CertNameConstraints_CheckNameSpaceNssNames(
        CERTGeneralName *nssSubjectNames,
        PKIX_PL_CertNameConstraints *nameConstraints,
        PKIX_Boolean *pCheckPass,
        void *plContext);

PKIX_Error *
pkix_pl_CertNameConstraints_CheckNameSpacePkixNames(
        PKIX_List *nameList,
        PKIX_PL_CertNameConstraints *nameConstraints,
        PKIX_Boolean *pCheckPass,
        void *plContext);


PKIX_Error *pkix_pl_CertNameConstraints_Merge(
        PKIX_PL_CertNameConstraints *firstNC,
        PKIX_PL_CertNameConstraints *secondNC,
        PKIX_PL_CertNameConstraints **pMergedNC,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_NAMECONSTRAINTS_H */
