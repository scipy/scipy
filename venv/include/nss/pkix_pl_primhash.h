/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_primhash.h
 *
 * Primitive Hashtable Definition
 *
 */

#ifndef _PKIX_PL_PRIMHASH_H
#define _PKIX_PL_PRIMHASH_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pkix_pl_HT_Elem pkix_pl_HT_Elem;

typedef struct pkix_pl_PrimHashTable pkix_pl_PrimHashTable;

typedef struct pkix_pl_Integer pkix_pl_Integer;

struct pkix_pl_Integer{
        PKIX_UInt32 ht_int;
};

struct pkix_pl_HT_Elem {
        void *key;
        void *value;
        PKIX_UInt32 hashCode;
        pkix_pl_HT_Elem *next;
};

struct pkix_pl_PrimHashTable {
        pkix_pl_HT_Elem **buckets;
        PKIX_UInt32 size;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_PrimHashTable_Create(
        PKIX_UInt32 numBuckets,
        pkix_pl_PrimHashTable **pResult,
        void *plContext);

PKIX_Error *
pkix_pl_PrimHashTable_Add(
        pkix_pl_PrimHashTable *ht,
        void *key,
        void *value,
        PKIX_UInt32 hashCode,
        PKIX_PL_EqualsCallback keyComp,
        void *plContext);

PKIX_Error *
pkix_pl_PrimHashTable_Remove(
        pkix_pl_PrimHashTable *ht,
        void *key,
        PKIX_UInt32 hashCode,
        PKIX_PL_EqualsCallback keyComp,
        void **pKey,
        void **pValue,
        void *plContext);

PKIX_Error *
pkix_pl_PrimHashTable_Lookup(
        pkix_pl_PrimHashTable *ht,
        void *key,
        PKIX_UInt32 hashCode,
        PKIX_PL_EqualsCallback keyComp,
        void **pResult,
        void *plContext);

PKIX_Error*
pkix_pl_PrimHashTable_Destroy(
        pkix_pl_PrimHashTable *ht,
        void *plContext);

PKIX_Error *
pkix_pl_PrimHashTable_GetBucketSize(
        pkix_pl_PrimHashTable *ht,
        PKIX_UInt32 hashCode,
        PKIX_UInt32 *pBucketSize,
        void *plContext);

PKIX_Error *
pkix_pl_PrimHashTable_RemoveFIFO(
        pkix_pl_PrimHashTable *ht,
        PKIX_UInt32 hashCode,
        void **pKey,
        void **pValue,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_PRIMHASH_H */
