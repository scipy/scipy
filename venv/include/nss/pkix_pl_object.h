/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_object.h
 *
 * Object Construction, Destruction and Callback Definitions
 *
 */

#ifndef _PKIX_PL_OBJECT_H
#define _PKIX_PL_OBJECT_H

#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Object Implementation Notes:
 *
 * Allocating a new object creates an object header and a block of
 * uninitialized user data. A pointer to this uninitialized data is
 * returned to the user. The structure looks as follows:
 *
 * +--------------------+
 * | MAGIC HEADER       |
 * | (object header)    |
 * +--------------------+
 * | user data          | -- pointer returned from PKIX_PL_Object_Alloc
 * +--------------------+
 *
 * Object operations receive a pointer to raw user data as an argument.
 * The macro HEADER(object) returns a pointer to the object header.
 * An assertion then verifies that the first field is the MAGIC_HEADER.
 */

/* PKIX_PL_Object Structure Definition */
struct PKIX_PL_ObjectStruct {
        PRUint64    magicHeader;
        PKIX_UInt32 type;
        PKIX_Int32 references;
        PRLock *lock;
        PKIX_PL_String *stringRep;
        PKIX_UInt32 hashcode;
        PKIX_Boolean hashcodeCached;
};

/* see source file for function documentation */

PKIX_Error *
pkix_pl_Object_RetrieveEqualsCallback(
        PKIX_PL_Object *object,
        PKIX_PL_EqualsCallback *equalsCallback,
        void *plContext);

extern PKIX_Boolean initializing;
extern PKIX_Boolean initialized;

#ifdef PKIX_USER_OBJECT_TYPE

extern PRLock *classTableLock;

#endif

extern pkix_ClassTable_Entry systemClasses[PKIX_NUMTYPES];

PKIX_Error *
pkix_pl_Object_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_OBJECT_H */
