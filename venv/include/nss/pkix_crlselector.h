/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_crlselector.h
 *
 * CrlSelector Object Type Definition
 *
 */

#ifndef _PKIX_CRLSELECTOR_H
#define _PKIX_CRLSELECTOR_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_CRLSelectorStruct {
        PKIX_CRLSelector_MatchCallback matchCallback;
        PKIX_ComCRLSelParams *params;
        PKIX_PL_Object *context;
};

/* see source file for function documentation */

PKIX_Error *pkix_CRLSelector_RegisterSelf(void *plContext);

PKIX_Error *
pkix_CRLSelector_Select(
	PKIX_CRLSelector *selector,
	PKIX_List *before,
	PKIX_List **pAfter,
	void *plContext);
#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CRLSELECTOR_H */
