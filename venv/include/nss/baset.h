/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BASET_H
#define BASET_H

/*
 * baset.h
 *
 * This file contains definitions for the basic types used throughout
 * nss but not available publicly.
 */

#ifndef NSSBASET_H
#include "nssbaset.h"
#endif /* NSSBASET_H */

#include "plhash.h"

PR_BEGIN_EXTERN_C

/*
 * nssArenaMark
 *
 * This type is used to mark the current state of an NSSArena.
 */

struct nssArenaMarkStr;
typedef struct nssArenaMarkStr nssArenaMark;

#ifdef DEBUG
/*
 * ARENA_THREADMARK
 *
 * Optionally, this arena implementation can be compiled with some
 * runtime checking enabled, which will catch the situation where
 * one thread "marks" the arena, another thread allocates memory,
 * and then the mark is released.  Usually this is a surprise to
 * the second thread, and this leads to weird runtime errors.
 * Define ARENA_THREADMARK to catch these cases; we define it for all
 * (internal and external) debug builds.
 */
#define ARENA_THREADMARK

/*
 * ARENA_DESTRUCTOR_LIST
 *
 * Unfortunately, our pointer-tracker facility, used in debug
 * builds to agressively fight invalid pointers, requries that
 * pointers be deregistered when objects are destroyed.  This
 * conflicts with the standard arena usage where "memory-only"
 * objects (that don't hold onto resources outside the arena)
 * can be allocated in an arena, and never destroyed other than
 * when the arena is destroyed.  Therefore we have added a
 * destructor-registratio facility to our arenas.  This was not
 * a simple decision, since we're getting ever-further away from
 * the original arena philosophy.  However, it was felt that
 * adding this in debug builds wouldn't be so bad; as it would
 * discourage them from being used for "serious" purposes.
 * This facility requires ARENA_THREADMARK to be defined.
 */
#ifdef ARENA_THREADMARK
#define ARENA_DESTRUCTOR_LIST
#endif /* ARENA_THREADMARK */

#endif /* DEBUG */

typedef struct nssListStr nssList;
typedef struct nssListIteratorStr nssListIterator;
typedef PRBool (*nssListCompareFunc)(void *a, void *b);
typedef PRIntn (*nssListSortFunc)(void *a, void *b);
typedef void (*nssListElementDestructorFunc)(void *el);

typedef struct nssHashStr nssHash;
typedef void(PR_CALLBACK *nssHashIterator)(const void *key, void *value,
                                           void *arg);

/*
 * nssPointerTracker
 *
 * This type is used in debug builds (both external and internal) to
 * track our object pointers.  Objects of this type must be statically
 * allocated, which means the structure size must be available to the
 * compiler.  Therefore we must expose the contents of this structure.
 * But please don't access elements directly; use the accessors.
 */

#ifdef DEBUG
struct nssPointerTrackerStr {
    PRCallOnceType once;
    PZLock *lock;
    PLHashTable *table;
};
typedef struct nssPointerTrackerStr nssPointerTracker;
#endif /* DEBUG */

/*
 * nssStringType
 *
 * There are several types of strings in the real world.  We try to
 * use only UTF8 and avoid the rest, but that's not always possible.
 * So we have a couple converter routines to go to and from the other
 * string types.  We have to be able to specify those string types,
 * so we have this enumeration.
 */

enum nssStringTypeEnum {
    nssStringType_DirectoryString,
    nssStringType_TeletexString, /* Not "teletext" with trailing 't' */
    nssStringType_PrintableString,
    nssStringType_UniversalString,
    nssStringType_BMPString,
    nssStringType_UTF8String,
    nssStringType_PHGString,
    nssStringType_GeneralString,

    nssStringType_Unknown = -1
};
typedef enum nssStringTypeEnum nssStringType;

PR_END_EXTERN_C

#endif /* BASET_H */
