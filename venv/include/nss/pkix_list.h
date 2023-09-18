/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_list.h
 *
 * List Object Type Definition
 *
 */

#ifndef _PKIX_LIST_H
#define _PKIX_LIST_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef PKIX_Error *
(*PKIX_List_SortComparatorCallback)(
        PKIX_PL_Object *obj1,
        PKIX_PL_Object *obj2,
        PKIX_Int32 *pResult,
        void *plContext);

struct PKIX_ListStruct {
        PKIX_PL_Object *item;
        PKIX_List *next;
        PKIX_Boolean immutable;
        PKIX_UInt32 length;
        PKIX_Boolean isHeader;
};

/* see source file for function documentation */

PKIX_Error *pkix_List_RegisterSelf(void *plContext);

PKIX_Error *
pkix_List_Contains(
        PKIX_List *list,
        PKIX_PL_Object *object,
        PKIX_Boolean *pFound,
        void *plContext);

PKIX_Error *
pkix_List_Remove(
        PKIX_List *list,
        PKIX_PL_Object *target,
        void *plContext);

PKIX_Error *
pkix_List_MergeLists(
        PKIX_List *firstList,
        PKIX_List *secondList,
        PKIX_List **pMergedList,
        void *plContext);

PKIX_Error *
pkix_List_AppendList(
        PKIX_List *toList,
        PKIX_List *fromList,
        void *plContext);

PKIX_Error *
pkix_List_AppendUnique(
        PKIX_List *toList,
        PKIX_List *fromList,
        void *plContext);

PKIX_Error *
pkix_List_RemoveItems(
        PKIX_List *list,
        PKIX_List *deleteList,
        void *plContext);

PKIX_Error *
pkix_List_QuickSort(
        PKIX_List *fromList,
        PKIX_List_SortComparatorCallback comparator,
        PKIX_List **pSortedList,
        void *plContext);

PKIX_Error *
pkix_List_BubbleSort(
        PKIX_List *fromList,
        PKIX_List_SortComparatorCallback comparator,
        PKIX_List **pSortedList,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_LIST_H */
