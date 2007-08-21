/*
 * phash.h -- $Id$
 * portability layer unique key hashing functions
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "extern_c.h"

typedef struct p_hashtab p_hashtab;
typedef unsigned long p_hashkey;

/* randomize the low order 32 bits of an address or integer
 *   such that P_IHASH(x)==P_IHASH(y) if and only if x==y */
#define P_IHASH(x) ((x)^p_hmasks[(((p_hashkey)(x))>>4)&0x3f])
#define P_PHASH(x) P_IHASH((char*)(x)-(char*)0)

extern p_hashkey p_hmasks[64];  /* for P_IHASH, P_PHASH macros */

/* unique key hash tables are basis for all hashing */
extern p_hashtab *p_halloc(p_hashkey size);
extern void p_hfree(p_hashtab *tab, void (*func)(void *));
extern int p_hinsert(p_hashtab *tab, p_hashkey hkey, void *value);
extern void *p_hfind(p_hashtab *tab, p_hashkey hkey);
extern void p_hiter(p_hashtab *tab,
                    void (*func)(void *val, p_hashkey key, void *ctx),
                    void *ctx);

/* global name string to id number correspondence
 *   p_id returns id number of an existing string, or 0
 *   p_idmake returns id number valid until matching p_idfree, never 0
 *   p_idstatic returns id number for statically allocated input name
 *     - name not copied, subsequent calls to p_idfree will be ignored
 *   p_idfree decrements the use counter for the given id number,
 *     freeing the number if there are no more uses
 *     - p_idmake increments use counter if name already exists
 *   p_idname returns 0 if no such id number has been made */
extern p_hashkey p_id(const char *name, int len);
extern p_hashkey p_idmake(const char *name, int len);
extern p_hashkey p_idstatic(char *name);
extern void p_idfree(p_hashkey id);
extern char *p_idname(p_hashkey id);

/* global pointer-to-pointer correspondence
 *   p_getctx returns context of a pointer or 0
 *   p_setctx sets context of a pointer, or deletes it if 0 */
extern void p_setctx(void *ptr, void *context);
extern void *p_getctx(void *ptr);

END_EXTERN_C
