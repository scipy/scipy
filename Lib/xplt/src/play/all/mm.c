/*
 * mm.c -- $Id$
 *
 * load default function pointers for pstdlib.h interface
 * p_mminit and p_mmdebug override these defaults
 * purpose of this file is so if p_mminit or p_mmdebug never called,
 *   code should run with system memory manager without loading bmm
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "pstdlib.h"

static void *p__malloc(size_t);
static void  p__free(void *);
static void *p__realloc(void *, size_t);

void *(*p_malloc)(size_t)= &p__malloc;
void  (*p_free)(void *)=   &p__free;
void *(*p_realloc)(void *, size_t)= &p__realloc;

static void *p__mmfail(unsigned long n);
void *(*p_mmfail)(unsigned long n)= &p__mmfail;

long p_nallocs = 0;
long p_nfrees = 0;
long p_nsmall = 0;
long p_asmall = 0;

p_twkspc p_wkspc;

static void *
p__malloc(size_t n)
{
  void *p = malloc(n>0? n : 1);
  if (!p) return p_mmfail(n>0? n : 1);
  p_nallocs++;
  return p;
}

static void
p__free(void *p)
{
  if (p) {
    p_nfrees++;
    free(p);
  }
}

static void *
p__realloc(void *p, size_t n)
{
  if (n<=0) n = 1;
  p = p? realloc(p,n) : malloc(n);
  return p? p : p_mmfail(n);
}

/* ARGSUSED */
static void *
p__mmfail(unsigned long n)
{
  return 0;
}
