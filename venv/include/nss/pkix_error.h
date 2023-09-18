/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_error.h
 *
 * Error Object Type Definition
 *
 */

#ifndef _PKIX_ERROR_H
#define _PKIX_ERROR_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_ErrorStruct {
        PKIX_ERRORCODE errCode;
        PKIX_ERRORCLASS errClass;  /* was formerly "code" */
        PKIX_Int32 plErr;
        PKIX_Error *cause;
        PKIX_PL_Object *info;
};

/* see source file for function documentation */

extern PKIX_Error * pkix_Error_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_ERROR_H */
