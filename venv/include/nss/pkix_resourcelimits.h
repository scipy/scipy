/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_resourcelimits.h
 *
 * ResourceLimits Object Type Definition
 *
 */

#ifndef _PKIX_RESOURCELIMITS_H
#define _PKIX_RESOURCELIMITS_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_ResourceLimitsStruct {
        PKIX_UInt32 maxTime;
        PKIX_UInt32 maxFanout;
        PKIX_UInt32 maxDepth;
        PKIX_UInt32 maxCertsNumber;
        PKIX_UInt32 maxCrlsNumber;
};

/* see source file for function documentation */

PKIX_Error *pkix_ResourceLimits_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_RESOURCELIMITS_H */
