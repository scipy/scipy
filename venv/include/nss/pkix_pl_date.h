/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_date.h
 *
 * Date Object Definitions
 *
 */

#ifndef _PKIX_PL_DATE_H
#define _PKIX_PL_DATE_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_DateStruct{
        PRTime nssTime;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_Date_ToString_Helper(
        SECItem *nssTime,
        PKIX_PL_String **pString,
        void *plContext);

PKIX_Error *pkix_pl_Date_RegisterSelf(void *plContext);

PKIX_Error *
pkix_pl_Date_GetPRTime(
        PKIX_PL_Date *date,
        PRTime *pPRTime,
        void *plContext);

PKIX_Error *
pkix_pl_Date_CreateFromPRTime(
        PRTime prtime,
        PKIX_PL_Date **pDate,
        void *plContext);

PKIX_Error *
PKIX_PL_Date_CreateFromPRTime(
        PRTime prtime,
        PKIX_PL_Date **pDate,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_DATE_H */
