/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_string.h
 *
 * String Object Definitions
 *
 */

#ifndef _PKIX_PL_STRING_H
#define _PKIX_PL_STRING_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_StringStruct {
        void* utf16String;
        PKIX_UInt32 utf16Length;
        /* XXX For Debugging  */
        char* escAsciiString;
        PKIX_UInt32 escAsciiLength;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_String_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_STRING_H */
