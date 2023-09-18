/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_colcertstore.h
 *
 * CollectionCertstore Object Type Definition
 *
 */

#ifndef _PKIX_PL_COLCERTSTORE_H
#define _PKIX_PL_COLCERTSTORE_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PKIX_PL_CollectionCertStoreContext {
        PKIX_PL_String *storeDir;
        PKIX_List *crlList;
        PKIX_List *certList;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_CollectionCertStoreContext_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_COLCERTSTORE_H */
