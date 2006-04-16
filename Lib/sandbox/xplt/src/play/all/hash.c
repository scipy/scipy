/*
 * hash.c -- $Id$
 * unique key hashing functions
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "phash.h"
#include "pstdlib.h"

typedef struct p_hashent p_hashent;

struct p_hashtab {
  p_hashkey mask;       /* nslots-1 */
  p_hashent **slots;
  p_hashent *freelist;
  p_hashent *entries;
  long nitems;          /* informational only */
};

struct p_hashent {
  p_hashent *next;
  p_hashkey hkey;   /* could keep key, but pexpand faster with this */
  void *value;
};

static p_hashent *p_hexpand(p_hashtab *tab);
static int p_hremove(p_hashtab *tab, p_hashkey hkey);

p_hashtab *
p_halloc(p_hashkey size)
{
  p_hashtab *tab;
  p_hashent *e;
  p_hashkey i;
  p_hashkey n = 4;
  if (size>100000) size = 100000;
  while (n<size) n<<= 1;
  n<<= 1;
  tab = p_malloc(sizeof(p_hashtab));
  tab->nitems = 0;
  tab->mask = n-1;
  tab->slots = p_malloc(sizeof(p_hashent *)*n);
  for (i=0 ; i<n ; i++) tab->slots[i] = 0;
  n>>= 1;
  e = p_malloc(sizeof(p_hashent)*n);
  for (i=0 ; i<n-1 ; i++) e[i].next = &e[i+1];
  e[i].next = 0;
  tab->entries = tab->freelist = e;
  return tab;
}

void
p_hfree(p_hashtab *tab, void (*func)(void *))
{
  p_hashent **slots = tab->slots;
  p_hashent *entries = tab->entries;
  if (func) {
    p_hashkey n = tab->mask+1;
    p_hashent *e;
    p_hashkey i;
    for (i=0 ; i<n ; i++)
      for (e=tab->slots[i] ; e ; e=e->next) func(e->value);
  }
  tab->slots = 0;
  tab->freelist = tab->entries = 0;
  p_free(slots);
  p_free(entries);
  p_free(tab);
}

void
p_hiter(p_hashtab *tab,
        void (*func)(void *val, p_hashkey key, void *ctx), void *ctx)
{
  p_hashkey n = tab->mask+1;
  p_hashent *e;
  p_hashkey i;
  for (i=0 ; i<n ; i++)
    for (e=tab->slots[i] ; e ; e=e->next)
      func(e->value, e->hkey, ctx);
}

int
p_hinsert(p_hashtab *tab, p_hashkey hkey, void *value)
{
  p_hashent *e;
  if (p_signalling) return 1;
  if (!value) return p_hremove(tab, hkey);
  for (e=tab->slots[hkey&tab->mask] ; e && e->hkey!=hkey ; e=e->next);
  if (!e) {
    e = tab->freelist;
    if (!e) {
      e = p_hexpand(tab);
      if (!e) return 1;
    }
    /* CRITICAL SECTION BEGIN */
    e->hkey = hkey;
    hkey &= tab->mask;
    tab->freelist = e->next;
    e->next = tab->slots[hkey];
    tab->slots[hkey] = e;
    /* CRITICAL SECTION END */
    tab->nitems++;
  }
  e->value = value;
  return 0;
}

void *
p_hfind(p_hashtab *tab, p_hashkey hkey)
{
  p_hashent *e;
  for (e=tab->slots[hkey&tab->mask] ; e ; e=e->next)
    if (e->hkey==hkey) return e->value;
  return 0;
}

static int
p_hremove(p_hashtab *tab, p_hashkey hkey)
{
  p_hashent *e, **pe = &tab->slots[hkey & tab->mask];
  for (e=*pe ; e ; e=*pe) {
    if (e->hkey==hkey) {
      /* CRITICAL SECTION BEGIN */
      *pe = e->next;
      e->next = tab->freelist;
      tab->freelist = e;
      /* CRITICAL SECTION END */
      tab->nitems--;
      break;
    }
    pe = &e->next;
  }
  return 0;
}

static p_hashent *
p_hexpand(p_hashtab *tab)
{
  p_hashkey n = tab->mask + 1;
  p_hashent **slots = p_malloc(sizeof(p_hashent *)*(n<<1));
  if (slots) {
    p_hashent *e = p_malloc(sizeof(p_hashent)*n);
    if (e) {
      p_hashent *e0, **pe, **qe;
      p_hashkey i;
      for (i=0 ; i<n ; i++) {
        pe = &slots[i];
        qe = &pe[n];
        for (e0=tab->slots[i] ; e0 ; e0=e0->next,e++) {
          /* exactly n/2 passes through this loop body */
          e->value = e0->value;
          if ((e->hkey = e0->hkey) & n) {
            *qe = e;
            qe = &e->next;
          } else {
            *pe = e;
            pe = &e->next;
          }
        }
        *pe = *qe = 0;
      }
      n>>= 1;
      for (i=0 ; i<n-1 ; i++) e[i].next = &e[i+1];
      e[i].next = 0;
      pe = tab->slots;
      e0 = tab->entries;
      /* CRITICAL SECTION BEGIN */
      tab->mask = (tab->mask<<1) | 1;
      tab->slots = slots;
      tab->entries = e-n;
      tab->freelist = e;
      /* CRITICAL SECTION END */
      p_free(pe);
      p_free(e0);
      return e;
    } else {
      p_free(slots);
    }
  }
  return 0;
}
