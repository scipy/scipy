/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_certselector.h
 *
 * CertSelector Object Type Definition
 *
 */

#ifndef _PKIX_CERTSELECTOR_H
#define _PKIX_CERTSELECTOR_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_CertSelectorStruct {
        PKIX_CertSelector_MatchCallback matchCallback;
        PKIX_ComCertSelParams *params;
        PKIX_PL_Object *context;
};

/* see source file for function documentation */

PKIX_Error *
pkix_CertSelector_Select(
	PKIX_CertSelector *selector,
	PKIX_List *before,
	PKIX_List **pAfter,
	void *plContext);

PKIX_Error *pkix_CertSelector_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_CERTSELECTOR_H */
