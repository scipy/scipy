/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * This file defines several platform independent functions to make system
 * calls in a portable manner.
 *
 */

#ifndef _PKIX_PL_SYSTEM_H
#define _PKIX_PL_SYSTEM_H

#include "pkixt.h"

#ifdef __cplusplus
extern "C" {
#endif

/* General
 *
 * Please refer to the libpkix Programmer's Guide for detailed information
 * about how to use the libpkix library. Certain key warnings and notices from
 * that document are repeated here for emphasis.
 *
 * All identifiers in this file (and all public identifiers defined in
 * libpkix) begin with "PKIX_". Private identifiers only intended for use
 * within the library begin with "pkix_".
 *
 * A function returns NULL upon success, and a PKIX_Error pointer upon failure.
 *
 * Unless otherwise noted, for all accessor (gettor) functions that return a
 * PKIX_PL_Object pointer, callers should assume that this pointer refers to a
 * shared object. Therefore, the caller should treat this shared object as
 * read-only and should not modify this shared object. When done using the
 * shared object, the caller should release the reference to the object by
 * using the PKIX_PL_Object_DecRef function.
 *
 * While a function is executing, if its arguments (or anything referred to by
 * its arguments) are modified, free'd, or destroyed, the function's behavior
 * is undefined.
 *
 */

/*
 * FUNCTION: PKIX_PL_Initialize
 * DESCRIPTION:
 *
 *  XXX If this function is really only meant to be used by PKIX_Initialize,
 *  why don't we just put it in a private header file rather than the public
 *  API. I think it may confuse users.
 *
 *  This function should NOT be called by applications. It is only meant to
 *  be used internally. The application needs only to call PKIX_Initialize,
 *  which in turn will call this function.
 *
 *  This function initializes data structures critical to the operation of
 *  libpkix. If initialization is not successful, an Error pointer is
 *  returned. This function should only be called once. If it is called more
 *  than once, the behavior is undefined.
 *
 *  No PKIX_* types and functions should be used before this function is
 *  called and returns successfully.
 *
 * PARAMETERS:
 *  "platformInitNeeded"
 *      Boolean indicating whether platform initialization is to be called
 *  "useArenas"
 *      Boolean indicating whether allocation is to be done using arenas or
 *      individual allocation (malloc).
 *  "pPlContext"
 *      Address at which platform-specific context pointer is stored. Must be
 *      non-NULL.
 * THREAD SAFETY:
 *  Not Thread Safe
 *
 *  This function assumes that no other thread is calling this function while
 *  it is executing.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Initialize(
        PKIX_Boolean platformInitNeeded,
        PKIX_Boolean useArenas,
        void **pPlContext);

/*
 * FUNCTION: PKIX_PL_Shutdown
 * DESCRIPTION:
 *
 *  XXX If this function is really only meant to be used by PKIX_Shutdown,
 *  why don't we just put it in a private header file rather than the public
 *  API. I think it may confuse users.
 *
 *  This function should NOT be called by applications. It is only meant to
 *  be used internally. The application needs only to call PKIX_Shutdown,
 *  which in turn will call this function.
 *
 *  This function deallocates any memory used by the Portability Layer (PL)
 *  component of the libpkix library and shuts down any ongoing operations.
 *  This function should only be called once. If it is called more than once,
 *  the behavior is undefined.
 *
 *  No PKIX_* types and functions should be used after this function is called
 *  and returns successfully.
 *
 * PARAMETERS:
 *  "platformInitNeeded"
 *      Boolean value of whether PKIX initialized NSS: PKIX_TRUE if we
 *      called nssInit, PKIX_FALSE otherwise
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe
 *
 *  This function makes use of global variables and should only be called once.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Shutdown(void *plContext);

/* standard memory management operations (not reference-counted) */

/*
 * FUNCTION: PKIX_PL_Malloc
 * DESCRIPTION:
 *
 *  Allocates a block of "size" bytes. The bytes are not initialized. A
 *  pointer to the newly allocated memory will be stored at "pMemory". The
 *  memory allocated by PKIX_PL_Malloc() may only be freed by PKIX_PL_Free().
 *  If "size" equals zero, this function stores NULL at "pMemory".
 *
 * PARAMETERS:
 *  "size"
 *      Number of bytes to allocate.
 *  "pMemory"
 *      Address where newly allocated pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on underlying thread safety of platform used by PL.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Malloc(
        PKIX_UInt32 size,
        void **pMemory,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Calloc
 * DESCRIPTION:
 *
 *  Allocates memory for an array of "nElem" elements, with each element
 *  requiring "elSize" bytes, and with all the bits initialized to zero. A
 *  pointer to the newly allocated memory will be stored at "pMemory". The
 *  memory allocated by PKIX_PL_Calloc() may only be freed by PKIX_PL_Free().
 *  If "nElem" equals zero or "elSize" equals zero, this function stores NULL
 *  at "pMemory".
 *
 * PARAMETERS:
 *  "nElem"
 *      Number of elements needed.
 *  "elSize"
 *      Number of bytes needed per element.
 *  "pMemory"
 *      Address where newly allocated pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on underlying thread safety of platform used by PL.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Calloc(
        PKIX_UInt32 nElem,
        PKIX_UInt32 elSize,
        void **pMemory,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Realloc
 * DESCRIPTION:
 *
 *  Resizes an existing block of memory (pointed to by "ptr") to "size" bytes.
 *  Stores a pointer to the resized memory at "pNewPtr". The "ptr" must
 *  originate from either PKIX_PL_Malloc(), PKIX_PL_Realloc(), or
 *  PKIX_PL_Calloc(). If "ptr" is NULL, this function behaves as if
 *  PKIX_PL_Malloc were called. If "ptr" is not NULL and "size" equals zero,
 *  the memory pointed to by "ptr" is deallocated and this function stores
 *  NULL at "pPtr".
 *
 * PARAMETERS:
 *  "ptr"
 *      A pointer to an existing block of memory.
 *  "size"
 *      New size in bytes.
 *  "pPtr"
 *      Address where newly allocated pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on underlying thread safety of platform used by PL.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Realloc(
        void *ptr,
        PKIX_UInt32 size,
        void **pNewPtr,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Free
 * DESCRIPTION:
 *
 *  Frees a block of memory pointed to by "ptr". This value must originate with
 *  either PKIX_PL_Malloc(), PKIX_PL_Calloc, or PKIX_PL_Realloc(). If "ptr" is
 *  NULL, the function has no effect.
 *
 * PARAMETERS:
 *  "ptr"
 *      A pointer to an existing block of memory.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on underlying thread safety of platform used by PL.
 * RETURNS:
 *  Returns NULL always.
 */
PKIX_Error *
PKIX_PL_Free(
        void *ptr,
        void *plContext);

/* Callback Types
 *
 * The next few typedefs define function pointer types for the standard
 * functions associated with every object type. See the Implementation
 * Guidelines or the comments below for more information.
 */

/*
 * TYPE: PKIX_PL_DestructorCallback
 * DESCRIPTION:
 *
 *  This callback function destroys (or DecRef's) any pointers contained in
 *  the user data for the Object pointed to by "object" before the Object is
 *  destroyed.
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object to destroy. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts (as long as they're not operating on the same
 *  object and nobody else is performing an operation on the object at the
 *  same time). Both of these conditions should be guaranteed by the fact that
 *  the object's ref count was reduced to 0 in a lock that's still held when
 *  this callback is called.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_PL_DestructorCallback)(
        PKIX_PL_Object *object,
        void *plContext);

/*
 * TYPE: PKIX_PL_EqualsCallback
 * DESCRIPTION:
 *
 *  This callback function compares the Object pointed to by "firstObject" with
 *  the Object pointed to by "secondObject" for equality and stores the result
 *  at "pResult" (PKIX_TRUE if equal; PKIX_FALSE if not).
 *
 * PARAMETERS:
 *  "firstObject"
 *      Address of first object to compare. Must be non-NULL.
 *  "secondObject"
 *      Address of second object to compare. Must be non-NULL.
 *  "pResult"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_PL_EqualsCallback)(
        PKIX_PL_Object *firstObject,
        PKIX_PL_Object *secondObject,
        PKIX_Boolean *pResult,
        void *plContext);

/*
 * TYPE: PKIX_PL_HashcodeCallback
 * DESCRIPTION:
 *
 *  This callback function computes the hashcode of the Object pointed to by
 *  "object" and stores the result at "pValue".
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose hashcode is desired. Must be non-NULL.
 *  "pValue"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_PL_HashcodeCallback)(
        PKIX_PL_Object *object,
        PKIX_UInt32 *pValue,
        void *plContext);

/*
 * TYPE: PKIX_PL_ToStringCallback
 * DESCRIPTION:
 *
 *  This callback function converts the Object pointed to by "object" to a
 *  string representation and stores the result at "pString".
 *
 * PARAMETERS:
 *  "object"
 *      Object to get a string representation from. Must be non-NULL.
 *  "pString"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_PL_ToStringCallback)(
        PKIX_PL_Object *object,
        PKIX_PL_String **pString,
        void *plContext);

/*
 * TYPE: PKIX_PL_ComparatorCallback
 * DESCRIPTION:
 *
 *  This callback function determines how the Object pointed to by
 *  "firstObject" compares to the Object pointed to by "secondObject" and
 *  stores the result at "pResult".
 *
 *  Result is less than 0 if firstObject < secondObject
 *  Result equals 0 if firstObject = secondObject
 *  Result is greater than 0 if firstObject > secondObject
 *
 * PARAMETERS:
 *  "firstObject"
 *      Address of the first Object to compare. Must be non-NULL.
 *  "secondObject"
 *      Address of the second Object to compare. Must be non-NULL.
 *  "pResult"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_PL_ComparatorCallback)(
        PKIX_PL_Object *firstObject,
        PKIX_PL_Object *secondObject,
        PKIX_Int32 *pResult,
        void *plContext);

/*
 * TYPE: PKIX_PL_DuplicateCallback
 * DESCRIPTION:
 *
 *  This callback function creates a copy of the Object pointed to by "object"
 *  and stores it at "pNewObject". Changes to the copy will not affect the
 *  original and vice versa.
 *
 *  Note that if "object" is immutable, the Duplicate callback function simply
 *  needs to increment the reference count on "object" and return a reference
 *  to "object".
 *
 * PARAMETERS:
 *  "object"
 *      Address of the object to be copied. Must be non-NULL.
 *  "pNewObject"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_PL_DuplicateCallback)(
        PKIX_PL_Object *object,
        PKIX_PL_Object **pNewObject,
        void *plContext);

/* reference-counted objects */

/*
 * FUNCTION: PKIX_PL_Object_Alloc
 * DESCRIPTION:
 *
 *  Allocates a new Object of type "type" with "size" bytes and stores the
 *  resulting pointer at "pObject". The reference count of the newly
 *  allocated object will be initialized to 1. To improve performance, each
 *  object maintains a small cache for the results of Hashcode and ToString.
 *  Mutable objects should call InvalidateCache whenever changes are made to
 *  the object's state (after creation). If an error occurs during allocation,
 *  "pObject" will be set to NULL. If "size" equals zero, this function creates
 *  an Object with a reference count of 1, and places a pointer to unallocated
 *  memory at "pMemory".
 *
 * PARAMETERS:
 *  "type"
 *      The type code of this object. See pkixt.h for codes. The type code
 *      must be previously registered with PKIX_PL_Object_RegisterType().
 *  "size"
 *      The number of bytes needed for this object.
 *  "pMemory"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Alloc(
        PKIX_TYPENUM type,
        PKIX_UInt32 size,
        PKIX_PL_Object **pObject,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_IsTypeRegistered
 * DESCRIPTION:
 *
 *  Checks whether "type" has been registered by a previous call to
 *  PKIX_PL_Object_RegisterType() and stores the Boolean result at "pBool".
 *  This function will typically only be called by constructors for specific
 *  types.
 *
 * PARAMETERS:
 *  "type"
 *      The type code to check if valid.
 *  "pBool"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_IsTypeRegistered(
        PKIX_UInt32 type,
        PKIX_Boolean *pBool,
        void *plContext);

#ifdef PKIX_USER_OBJECT_TYPE
/*
 * FUNCTION: PKIX_PL_Object_RegisterType
 * DESCRIPTION:
 *
 *  Registers a new Object with type value "type" and associates it with a set
 *  of functions ("destructor", "equalsFunction", "hashcodeFunction",
 *  "toStringFunction", "comparator", "duplicateFunction"). The new type value
 *  is also associated with a string pointed to by "description", which is used
 *  by the default ToStringCallback. This function may only be called with a
 *  particular "type" value once. If "destructor", "equalsFunction",
 *  "hashcodeFunction", or "toStringFunction" are NULL, default functions will
 *  be registered. However, if "comparator" and "duplicateFunction" are NULL,
 *  no functions will be registered and calls to PKIX_PL_Object_Compare and
 *  PKIX_PL_Object_Duplicate will result in an error.
 *
 * PARAMETERS:
 *  "type"
 *      The type code.
 *  "description"
 *      The string used by the default ToStringCallback. Default used if NULL.
 *  "destructor"
 *      The DestructorCallback function to be set. Default used if NULL.
 *  "equalsFunction"
 *      The EqualsCallback function to be set. Default used if NULL.
 *  "hashcodeFunction"
 *      The HashcodeCallback function to be set. Default used if NULL.
 *  "toStringFunction"
 *      The ToStringCallback function to be set. Default used if NULL.
 *  "comparator"
 *      The ComparatorCallback function to be set. None set if NULL. If no
 *      callback function is set in this field, calls to
 *      PKIX_PL_Object_Compare() will result in an error.
 *  "duplicateFunction"
 *      The DuplicateCallback function to be set. None set if NULL. If no
 *      callback function is set in this field, calls to
 *      PKIX_PL_Object_Duplicate() will result in an error.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Object Error if "type" is already registered.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_RegisterType(
        PKIX_UInt32 type,
        char *description,
        PKIX_PL_DestructorCallback destructor,
        PKIX_PL_EqualsCallback equalsFunction,
        PKIX_PL_HashcodeCallback hashcodeFunction,
        PKIX_PL_ToStringCallback toStringFunction,
        PKIX_PL_ComparatorCallback comparator,
        PKIX_PL_DuplicateCallback duplicateFunction,
        void *plContext);

#endif
/*
 * FUNCTION: PKIX_PL_Object_InvalidateCache
 * DESCRIPTION:
 *
 *  Invalidates the cache of the Object pointed to by "object". The cache
 *  contains results of Hashcode and ToString. This function should be used by
 *  mutable objects whenever changes are made to the Object's state (after
 *  creation).
 *
 *  For example, if ToString is called on a mutable Object, the result will be
 *  computed, cached, and returned. If the Object's state does not change, a
 *  subsequent call to ToString will recognize that the relevant result is
 *  cached and will simply return the result (without calling the Object's
 *  ToStringCallback to recompute it). However, when the Object's state
 *  changes, the cache needs to be invalidated in order to force a subsequent
 *  call to ToString to recompute the result.
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose cache is to be invalidated. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 *
 * THREAD SAFETY
 *  Thread Safe - Object Type Table is locked during modification.
 *
 *  Multiple threads can safely call this function without worrying about
 *  conflicts, even if they're operating on the same object.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_InvalidateCache(
        PKIX_PL_Object *object,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_IncRef
 * DESCRIPTION:
 *
 *  Increments the reference count of the Object pointed to by "object".
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose reference count is to be incremented.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_IncRef(
        PKIX_PL_Object *object,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_DecRef
 * DESCRIPTION:
 *
 *  Decrements the reference count of the Object pointed to by "object". If the
 *  resulting reference count is zero, the destructor (if any) registered for
 *  the Object's type (by PKIX_PL_RegisterType) will be called and then the
 *  Object will be destroyed.
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose reference count is to be decremented.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  If destructor is not called, multiple threads can safely call this function
 *  without worrying about conflicts, even if they're operating on the same
 *  object. If destructor is called, thread safety depends on the callback
 *  defined by PKIX_PL_RegisterType().
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_DecRef(
        PKIX_PL_Object *object,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_Equals
 * DESCRIPTION:
 *
 *  Compares the Object pointed to by "firstObject" with the Object pointed to
 *  by "secondObject" for equality using the callback function registered for
 *  "firstObject"'s type, and stores the Boolean result at "pResult". While
 *  typical callbacks will return PKIX_FALSE if the objects are of different
 *  types, other callbacks may be capable of comparing objects of different
 *  types [which may correctly result in cases where Equals(first, second)
 *  differs from Equals(second, first)].
 *
 * PARAMETERS:
 *  "firstObject"
 *      Address of the first Object to compare. Must be non-NULL.
 *      The EqualsCallback for this Object will be called.
 *  "secondObject"
 *      Address of the second Object to compare. Must be non-NULL.
 *  "pResult"
 *      Address where Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on the callback defined by PKIX_PL_RegisterType().
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Object Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Equals(
        PKIX_PL_Object *firstObject,
        PKIX_PL_Object *secondObject,
        PKIX_Boolean *pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_Hashcode
 * DESCRIPTION:
 *
 *  Computes a hashcode of the Object pointed to by "object" using the
 *  callback registered for "object"'s type and stores it at "pValue". Two
 *  objects which are equal should have the same hashcode. Once a call to
 *  Hashcode has been made, the results are cached and subsequent calls to
 *  Hashcode will return the cached value. For mutable objects, an
 *  InvalidateCache function is provided, which should be called whenever
 *  changes are made to the object's state (after creation).
 *
 * PARAMETERS:
 *  "object"
 *      Address of the Object whose hashcode is desired. Must be non-NULL.
 *      The HashcodeCallback for this object will be called.
 *  "pValue"
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 *
 * THREAD SAFETY:
 *  Thread safety depends on the callback defined by PKIX_PL_RegisterType().
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Object Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Hashcode(
        PKIX_PL_Object *object,
        PKIX_UInt32 *pValue,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_ToString
 * DESCRIPTION:
 *
 *  Creates a string representation of the Object pointed to by "object" using
 *  the callback registered for "object"'s type and stores it at "pString".
 *  Once a call to ToString has been made, the results are cached and
 *  subsequent calls to ToString will return the cached value. For mutable
 *  objects, an InvalidateCache function is provided, which should be called
 *  whenever changes are made to the object's state (after creation).
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose string representation is desired.
 *      Must be non-NULL. The ToStringCallback for this object will be called.
 *  "pString"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on the callback defined by PKIX_PL_RegisterType().
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Object Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_ToString(
        PKIX_PL_Object *object,
        PKIX_PL_String **pString,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_Compare
 * DESCRIPTION:
 *
 *  Compares the Object pointed to by "firstObject" and the Object pointed to
 *  by "secondObject" using the comparator registered for "firstObject"'s type
 *  and stores the result at "pResult". Different types may be compared. This
 *  may correctly result in cases where Compare(first, second) is not the
 *  opposite of Compare(second, first). The PKIX_Int32 value stored at
 *  "pResult" will be:
 *   Less than 0 if "firstObject" < "secondObject"
 *   Equals to 0 if "firstObject" = "secondObject"
 *   Greater than 0 if "firstObject" > "secondObject"
 *
 * PARAMETERS:
 *  "firstObject"
 *      Address of first Object to compare. Must be non-NULL.
 *      The ComparatorCallback for this object will be called.
 *  "secondObject"
 *      Address of second object to compare. Must be non-NULL.
 *  "pResult
 *      Address where PKIX_Int32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on the comparator defined by PKIX_PL_RegisterType().
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Object Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Compare(
        PKIX_PL_Object *firstObject,
        PKIX_PL_Object *secondObject,
        PKIX_Int32 *pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_Duplicate
 * DESCRIPTION:
 *
 *  Creates a duplicate copy of the Object pointed to by "object" using the
 *  callback registered for "object"'s type and stores the copy at
 *  "pNewObject". Changes to the new object will not affect the original and
 *  vice versa.
 *
 *  Note that if "object" is immutable, the Duplicate callback function simply
 *  needs to increment the reference count on "object" and return a reference
 *  to "object".
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object to be duplicated. Must be non-NULL.
 *      The DuplicateCallback for this Object will be called.
 *  "pNewObject"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread safety depends on the callback defined by PKIX_PL_RegisterType().
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Object Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Duplicate(
        PKIX_PL_Object *object,
        PKIX_PL_Object **pNewObject,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_GetType
 * DESCRIPTION:
 *
 *  Retrieves the type code of the Object pointed to by "object" and stores it
 *  at "pType". See pkixt.h for type codes.
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose type is desired. Must be non-NULL.
 *  "pType"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_GetType(
        PKIX_PL_Object *object,
        PKIX_UInt32 *pType,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_Lock
 * DESCRIPTION:
 *
 *  Locks the Mutex associated with the Object pointed to by "object". When an
 *  object is created, it is associated with an object-specific Mutex to allow
 *  for synchronization when the fields of the object are modified.
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose Mutex is to be locked. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Lock(
        PKIX_PL_Object *object,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Object_Unlock
 * DESCRIPTION:
 *
 *  Unlocks the Mutex associated with the Object pointed to by "object". When
 *  an object is created, it is associated with an object-specific Mutex to
 *  allow for synchronization when the fields of the object are modified.
 *
 * PARAMETERS:
 *  "object"
 *      Address of Object whose Mutex is to be unlocked. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Object_Unlock(
        PKIX_PL_Object *object,
        void *plContext);

/* mutexes (locks) */

/*
 * FUNCTION: PKIX_PL_Mutex_Create
 * DESCRIPTION:
 *
 *  Creates a new Mutex and stores it at "pNewLock".
 *
 * PARAMETERS:
 *  "pNewLock"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Mutex_Create(
        PKIX_PL_Mutex **pNewLock,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Mutex_Lock
 * DESCRIPTION:
 *
 *  Locks the Mutex pointed to by "lock". If the Mutex is already locked, this
 *  function will block the current thread until the mutex can be locked by
 *  this thread.
 *
 * PARAMETERS:
 *  "lock"
 *      Address of Mutex to lock. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Mutex_Lock(
        PKIX_PL_Mutex *lock,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Mutex_Unlock
 * DESCRIPTION:
 *
 *  Unlocks the Mutex pointed to by "lock" if the current thread holds the
 *  Mutex.
 *
 * PARAMETERS:
 *  "lock"
 *      Address of Mutex to unlock. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Mutex_Unlock(
        PKIX_PL_Mutex *lock,
        void *plContext);

/* monitor (locks) */

/*
 * FUNCTION: PKIX_PL_MonitorLock_Create
 * DESCRIPTION:
 *
 *  Creates a new PKIX_PL_MonitorLock and stores it at "pNewLock".
 *
 * PARAMETERS:
 *  "pNewLock"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_MonitorLock_Create(
        PKIX_PL_MonitorLock **pNewLock,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_MonitorLock_Enter
 * DESCRIPTION:
 *
 *  Locks the MonitorLock pointed to by "lock". If the MonitorLock is already
 *  locked by other thread, this function will block the current thread. If
 *  the "lock" had been locked by current thread, this function will NOT block.
 *
 * PARAMETERS:
 *  "lock"
 *      Address of MonitorLock to lock. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_MonitorLock_Enter(
        PKIX_PL_MonitorLock *lock,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_MonitorLock_Exit
 * DESCRIPTION:
 *
 *  Unlocks the MonitorLock pointed to by "lock" if the lock counter of 
 *  current thread holds the MonitorLock reach 0, the lock is released.
 *
 * PARAMETERS:
 *  "lock"
 *      Address of MonitorLock to unlock. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_MonitorLock_Exit(
        PKIX_PL_MonitorLock *lock,
        void *plContext);

/* strings and formatted printing */

/*
 * FUNCTION: PKIX_PL_String_Create
 * DESCRIPTION:
 *
 *  Creates a new String using the data pointed to by "pString", the
 *  PKIX_UInt32 pointed to by "stringLen", and the PKIX_UInt32 pointed to by
 *  "fmtIndicator" and stores it at "pString". If the format is PKIX_ESCASCII
 *  the "stringLen" parameter is ignored and the string extends until a zero
 *  byte is  found. Once created, a String object is immutable.
 *
 *  Valid formats are:
 *      PKIX_ESCASCII
 *      PKIX_ESCASCII_DEBUG
 *      PKIX_UTF8
 *      PKIX_UTF8_NULL_TERM
 *      PKIX_UTF16
 *
 * PARAMETERS:
 *  "fmtIndicator"
 *      Format that "stringRep" is encoded with. Must be non-NULL.
 *  "stringRep"
 *      Address of encoded string representation. Must be non-NULL.
 *  "stringLen"
 *      Length of data stored at stringRep.
 *  "pString"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a String Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_String_Create(
        PKIX_UInt32 fmtIndicator,
        const void *stringRep,
        PKIX_UInt32 stringLen,
        PKIX_PL_String **pString,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_Sprintf
 * DESCRIPTION:
 *
 *  Creates a formatted string at "pOut" using the given format "fmt" and a
 *  variable length list of arguments. The format flags are identical to
 *  standard C with the exception that %s expects a PKIX_PL_String*, rather
 *  than a char *, and that {%d, %i, %o, %u, %x, %X} expect PKIX_UInt32 or
 *  PKIX_Int32 instead of int or unsigned int.
 *
 * PARAMETERS:
 *  "pOut"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 *  "fmt"
 *      Address of format string. Must be non-NULL.
 * THREAD SAFETY:
 *  Not Thread Safe - Caller must have exclusive access to all arguments.
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a String Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_Sprintf(
        PKIX_PL_String **pOut,
        void *plContext,
        const PKIX_PL_String *fmt, ...);

/*
 * FUNCTION: PKIX_PL_GetString
 * DESCRIPTION:
 *
 *  Retrieves the String associated with the value of "stringID" (if any) and
 *  stores it at "pString". If no such string is associated with "stringID",
 *  this function uses "defaultString" to create a String and stores it at
 *  "pString".
 *
 * PARAMETERS:
 *  "stringID"
 *      PKIX_UInt32 valud of string identifier.
 *  "defaultString"
 *      Address of a PKIX_ESCASCII encoded string representation.
 *      Must be non-NULL.
 *  "pString"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a String Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_GetString(
        PKIX_UInt32 stringID,
        char *defaultString,
        PKIX_PL_String **pString,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_String_GetEncoded
 * DESCRIPTION:
 *
 *  Retrieves the value of the String pointed to by "string" in the encoding
 *  specified by "fmtIndicator" and stores the result in "pStringRep" and
 *  "pLength", respectively. Note that "pStringRep" is not reference counted
 *  and will need to be freed with PKIX_PL_Free().
 *
 * PARAMETERS:
 *  "string"
 *      Address of String whose encoded value is desired. Must be non-NULL.
 *  "fmtIndicator"
 *      Format of encoding. Supported formats are:
 *      PKIX_ESCASCII, PKIX_ESCASII_DEBUG, PKIX_UTF8, PKIX_UTF8_NULL_TERM, and
 *      PKIX_UTF16. XXX Where are these documented?
 *  "pStringRep"
 *      Address where pointer to encoded value will be stored.
 *      Must be non-NULL.
 *  "pLength"
 *      Address where byte length of encoded value will be stored.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a String Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_String_GetEncoded(
        PKIX_PL_String *string,
        PKIX_UInt32 fmtIndicator,
        void **pStringRep,
        PKIX_UInt32 *pLength,
        void *plContext);

/*
 * Hashtable
 *
 * A hashtable is a very efficient data structure used for mapping keys to
 * values. Any non-null PKIX_PL_Object can be used as a key or as a value,
 * provided that it correctly implements the PKIX_PL_EqualsCallback and the
 * PKIX_PL_HashcodeCallback. A hashtable consists of several buckets, with
 * each bucket capable of holding a linked list of key/value mappings. When
 * adding, retrieving, or deleting a value, the hashcode of the key is used to
 * determine which bucket's linked list is relevant. The corresponding
 * key/value pair is then appended, retrieved, or deleted.
 */

/*
 * FUNCTION: PKIX_PL_HashTable_Create
 * DESCRIPTION:
 *
 *  Creates a new Hashtable with an initial capacity of "numBuckets" buckets
 *  and "maxEntriesPerBucket" of entries limit for each bucket and stores it
 *  at "pResult".
 *
 * PARAMETERS:
 *  "numBuckets"
 *      The initial number of hash table buckets. Must be non-zero.
 *  "maxEntriesPerBucket"
 *      The limit of entries per bucket. Zero means no limit.
 *  "pResult"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_HashTable_Create(
        PKIX_UInt32 numBuckets,
        PKIX_UInt32 maxEntriesPerBucket,
        PKIX_PL_HashTable **pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_HashTable_Add
 * DESCRIPTION:
 *
 *  Adds a key/value mapping using the Objects pointed to by "key" and "value"
 *  to the Hashtable pointed to by "ht".
 *
 *  Function increments key/value reference counts. Caller is responsible to
 *  to decrement(destroy) key/value ref counts(objects). 
 *
 * PARAMETERS:
 *  "ht"
 *      Address of Hashtable to be added to. Must be non-NULL.
 *  "key"
 *      Address of Object to be associated with "value". Must be non-NULL.
 *  "value"
 *      Address of Object to be added to Hashtable. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "ht"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Hashtable Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_HashTable_Add(
        PKIX_PL_HashTable *ht,
        PKIX_PL_Object *key,
        PKIX_PL_Object *value,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_HashTable_Remove
 * DESCRIPTION:
 *
 *  Removes the Object value whose key is equal to the Object pointed to by
 *  "key" from the Hashtable pointed to by "ht". If no such object exists,
 *  this function throws an Error.
 *
 *  Function frees "value" object. Caller is responsible to free "key"
 *  object.
 *
 * PARAMETERS:
 *  "ht"
 *      Address of Hashtable to remove object from. Must be non-NULL.
 *  "key"
 *      Address of Object used for lookup. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "ht"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Hashtable Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_HashTable_Remove(
        PKIX_PL_HashTable *ht,
        PKIX_PL_Object *key,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_HashTable_Lookup
 * DESCRIPTION:
 *
 *  Retrieves the Object whose key equals the Object pointed to by "key" from
 *  the Hashtable associated with "ht" and stores it at "pResult". If no
 *  Object is found, this function stores NULL at "pResult".
 *
 * PARAMETERS:
 *  "ht"
 *      Address of Hashtable to lookup Object from. Must be non-NULL.
 *  "key"
 *      Address of key Object used for lookup. Must be non-NULL.
 *  "pResult"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Hashtable Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_HashTable_Lookup(
        PKIX_PL_HashTable *ht,
        PKIX_PL_Object *key,
        PKIX_PL_Object **pResult,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_ByteArray_Create
 * DESCRIPTION:
 *
 *  Creates a new ByteArray using "length" bytes of data pointed to by "array"
 *  and stores it at "pByteArray". Once created, a ByteArray is immutable.
 *
 * PARAMETERS:
 *  "array"
 *      Address of source data.
 *  "length"
 *      Number of bytes to copy.
 *  "pByteArray"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_ByteArray_Create(
        void *array,
        PKIX_UInt32 length,
        PKIX_PL_ByteArray **pByteArray,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_ByteArray_GetPointer
 * DESCRIPTION:
 *
 *  Allocates enough memory to hold the contents of the ByteArray pointed to
 *  by "byteArray", copies the data from the ByteArray pointed to by
 *  "byteArray" into the newly allocated memory, and stores a pointer to the
 *  memory at "pArray". Note that "pArray" is not reference counted. It will
 *  need to be freed with PKIX_PL_Free().
 *
 * PARAMETERS:
 *  "byteArray"
 *      Address of ByteArray whose data is desired. Must be non-NULL.
 *  "pArray"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_ByteArray_GetPointer(
        PKIX_PL_ByteArray *byteArray,
        void **pArray,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_ByteArray_GetLength
 * DESCRIPTION:
 *
 *  Retrieves the length of the ByteArray pointed to by "byteArray" and stores
 *  the length at "pLength".
 *
 * PARAMETERS:
 *  "byteArray"
 *      Address of ByteArray whose length is desired. Must be non-NULL.
 *  "pLength"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_ByteArray_GetLength(
        PKIX_PL_ByteArray *byteArray,
        PKIX_UInt32 *pLength,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_OID_Create
 * DESCRIPTION:
 *
 *  Creates a new OID using NSS oid tag.
 *
 * PARAMETERS:
 *  "idtag"
 *      nss oid id tag.
 *  "pOID"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an OID Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_OID_Create(
        SECOidTag idtag,
        PKIX_PL_OID **pOID,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_OID_CreateBySECItem
 * DESCRIPTION:
 *
 *  Creates a new OID using a DER encoded OID stored as SECItem.
 *
 * PARAMETERS:
 *  "derOid"
 *      Address of SECItem that holds DER encoded OID.
 *  "pOID"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an OID Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_OID_CreateBySECItem(
        SECItem *derOid,
        PKIX_PL_OID **pOID,
        void *plContext);

/*
 * FUNCTION: PKIX_PL_BigInt_Create
 * DESCRIPTION:
 *
 *  Creates a new BigInt using the source String pointed to by "stringRep" and
 *  stores it at "pBigInt". Valid source Strings consist of an even number of
 *  hexadecimal digits, which are always interpreted as a positive number.
 *  Once created, a BigInt is immutable.
 *
 *  The regexp format is:
 *      HexDigit  ::= [0-9] | [A-F] | [a-f]
 *      DoubleHex ::= HexDigit HexDigit
 *      BigIntSrc ::= (DoubleHex)+
 *
 *  Note that since we are using DoubleHex, the number of characters in the
 *  source MUST be even. Additionally, the first DoubleHex MUST NOT be "00"
 *  unless it is the only DoubleHex.
 *
 *  Valid  :    "09"
 *  Valid  :    "00"    (special case where first and only DoubleHex is "00")
 *  Invalid:    "9"     (not DoubleHex: odd number of characters)
 *  Invalid:    "0009"  (first DoubleHex is "00")
 *
 *  XXX Why does this take a String object while OID_Create takes a char* ?
 *  Perhaps because OID_Create is often used with constant strings and
 *  this function isn't. That's a good reason, but we should explain it
 *  (if it's right)
 * PARAMETERS:
 *  "stringRep"
 *      Address of String representing a BigInt. Must be non-NULL.
 *  "pBigInt"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a BigInt Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_PL_BigInt_Create(
        PKIX_PL_String *stringRep,
        PKIX_PL_BigInt **pBigInt,
        void *plContext);

#ifdef __cplusplus
}
#endif

/*
 * FUNCTION: PKIX_PL_GetPLErrorCode
 * DESCRIPTION:
 *
 *  Returns error code from PL layer.
 *
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  PL layer error code. 
 */
int
PKIX_PL_GetPLErrorCode();

#endif /* _LIBPKIX_SYSTEM_H */
