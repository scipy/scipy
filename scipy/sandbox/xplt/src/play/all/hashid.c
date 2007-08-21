/*
 * hashid.c -- $Id$
 * global name string to id number correspondence
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "phash.h"
#include "pstdlib.h"

#include <string.h>

typedef struct id_name id_name;
struct id_name {
  union {
    char *ame;
    id_name *ext;
  } n;
  /* uses >= 0 indicates number of additional uses, ordinary case
   * uses ==-1 indicates idstatic, n.ame will never be freed
   * uses <=-2 indicates locked id_name struct, but
   *           n.ame can still be freed (true uses = -2-uses) */
  long uses;
};

static p_hashtab *id_table = 0;
static id_name *id_freelist = 0;
static id_name id_null;
extern int p_id_collisions; 
int p_id_collisions = 0;  /* expect to be tiny */

static p_hashkey id_hash(const char *name, int len);
static id_name *id_block(void);
static int idnm_compare(const char *n1, const char *n2, int len);

p_hashkey
p_id(const char *name, int len)
{
  p_hashkey h = id_hash(name, len);
  if (id_table) {
    id_name *idnm;
    p_hashkey rehash = h&0xfff;

    for (;;) {
      if (!h) h = 1;
      idnm = p_hfind(id_table, h);
      if (!idnm || !idnm->n.ame) return 0;
      if (!idnm_compare(name, idnm->n.ame, len)) return h;
      if (!rehash) rehash = 3691;
      h += rehash;
    }
  }
  return 0;
}

p_hashkey
p_idmake(const char *name, int len)
{
  id_name *idnm = 0;
  p_hashkey h = id_hash(name, len);
  p_hashkey rehash = h&0xfff;

  if (!id_table) {
    id_table = p_halloc(64);
    id_null.n.ame = 0;
    id_null.uses = -1;
    p_hinsert(id_table, 0, &id_null);
  }

  for (;;) {
    if (!h) h = 1;
    idnm = p_hfind(id_table, h);
    if (!idnm || !idnm->n.ame) break;
    if (!idnm_compare(name, idnm->n.ame, len)) {
      if (idnm->uses>=0) idnm->uses++;
      else if (idnm->uses<-1) idnm->uses--;
      return h;
    }
    if (!rehash) rehash = 3691;
    h += rehash;
    /* a collision locks the hashkey and id_name struct forever */
    if (idnm->uses>=0) {
      idnm->uses = -2-idnm->uses;
      p_id_collisions++;
    }
  }

  if (!idnm) {
    idnm = id_freelist;
    if (!idnm) idnm = id_block();
    idnm->uses = 0;
    id_freelist = idnm->n.ext;
    idnm->n.ame = len>0? p_strncat(0, name, len) : p_strcpy(name);
    p_hinsert(id_table, h, idnm);
  } else {
    idnm->n.ame = len>0? p_strncat(0, name, len) : p_strcpy(name);
  }

  return h;
}

static int
idnm_compare(const char *n1, const char *n2, int len)
{
  if (len) return strncmp(n1, n2, len) || n2[len];
  else     return strcmp(n1, n2);
}

p_hashkey
p_idstatic(char *name)
{
  p_hashkey id = p_idmake(name, 0);
  id_name *idnm = p_hfind(id_table, id);  /* never 0 or idnm->n.ame==0 */
  if (idnm->uses>=0) {
    char *nm = idnm->n.ame;
    idnm->uses = -1;
    idnm->n.ame = name;
    p_free(nm);
  }
  return id;
}

void
p_idfree(p_hashkey id)
{
  if (id_table) {
    id_name *idnm = p_hfind(id_table, id);
    if (idnm && idnm->n.ame) {
      if (!idnm->uses) {
        char *name = idnm->n.ame;
        p_hinsert(id_table, id, (void *)0);
        idnm->n.ext = id_freelist;
        id_freelist = idnm;
        p_free(name);
      } else if (idnm->uses>0) {
        idnm->uses--;
      } else if (idnm->uses==-2) {
        char *name = idnm->n.ame;
        idnm->n.ame = 0;
        p_free(name);
      } else if (idnm->uses<-2) {
        idnm->uses++;
      }
    }
  }
}

char *
p_idname(p_hashkey id)
{
  id_name *idnm = id_table? p_hfind(id_table, id) : 0;
  return idnm? idnm->n.ame : 0;
}

static id_name *
id_block(void)
{
  int i, n = 128;
  id_name *idnm = p_malloc(sizeof(id_name)*n);
  for (i=0 ; i<n-1 ; i++) idnm[i].n.ext = &idnm[i+1];
  idnm[i].n.ext = 0;
  return id_freelist = idnm;
}

static p_hashkey
id_hash(const char *name, int len)
{
  p_hashkey h = 0x34da3723;  /* random bits */
  if (len) for (; len && name[0] ; len--,name++)
    h = ((h<<8) + (h>>1)) ^ (unsigned char)name[0];
  else while (*name)
    h = ((h<<8) + (h>>1)) ^ (unsigned char)(*name++);
  return h;
}
