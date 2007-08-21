/*
 * mminit.c -- $Id$
 * basic memory management interface
 *
 * compile this file twice:
 * P_DEBUG not defined makes production version (p_mminit)
 * P_DEBUG defined makes debug version (p_mmdebug)
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "pstdlib.h"

typedef union mm_block mm_block;
typedef struct mm_arena mm_arena;

union mm_block {
  mm_block *next_free;  /* when not in use */
  mm_arena *arena;      /* when in use */
  /* remainder assure most restrictive alignment */
  long l;
  double d;   /* long double tickles bug for linux gcc 2.7.2 */
  void (*f)(void);
};

#ifndef P_DEBUG
# define MM_OFFSET 1
# define MM_EXTRA 1
# define MM_GUARD(block, n)
# define MM_CHECK(p, n)
# define MM_VCHECK(p)
#else
# define MM_OFFSET 3
# define MM_EXTRA 4
# define MM_GUARD(block, n) p_mmguard(block, n);
# define MM_CHECK(p, n) if (p_mmcheck(p)) return p_mmfail(n);
# define MM_VCHECK(p) if (p_mmcheck(p)) { p_mmfail(0); return; }
long p_mmoffset = MM_OFFSET*sizeof(mm_block);
long p_mmextra = MM_EXTRA*sizeof(mm_block);
#endif

static void *bmalloc(size_t);
static void  bfree(void *);
static void *brealloc(void *, size_t);

struct mm_arena {
  mm_block *next_free;
  mm_block *blocks;     /* current chunk */
  long n_blocks;        /* total number of mm_blocks per chunk */
  long size;            /* number of mm_blocks per unit */
};

/* the numbers are the number of allocation units per chunk
 * and the number of mm_blocks per unit (augmented by MM_EXTRA) */
static mm_arena arenas[6] = {
  { 0, 0, 512*( 1+MM_EXTRA),  1+MM_EXTRA },
  { 0, 0, 384*( 2+MM_EXTRA),  2+MM_EXTRA },
  { 0, 0, 256*( 4+MM_EXTRA),  4+MM_EXTRA },
  { 0, 0, 128*( 8+MM_EXTRA),  8+MM_EXTRA },
  { 0, 0, 64 *(16+MM_EXTRA), 16+MM_EXTRA },
  { 0, 0, 32 *(32+MM_EXTRA), 32+MM_EXTRA }};

static mm_arena *arena_list[32] = {
  arenas,  arenas+1,arenas+2,arenas+2,arenas+3,arenas+3,arenas+3,arenas+3,
  arenas+4,arenas+4,arenas+4,arenas+4,arenas+4,arenas+4,arenas+4,arenas+4,
  arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,
  arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,arenas+5,arenas+5 };

static mm_block *expand_arena(mm_arena *arena);

void
p_mminit(void)  /* actually p_mmdebug if P_DEBUG */
{
  p_malloc =  &bmalloc;
  p_free =    &bfree;
  p_realloc = &brealloc;
}

static void *
bmalloc(size_t n)
{
  mm_block *block;
  if (n<=0) n = 1;
  if (n <= 32*sizeof(mm_block)) {
    mm_arena *arena = arena_list[(n-1)/sizeof(mm_block)];
    block = arena->next_free;
    if (!block) {
      block = expand_arena(arena);
      if (!block) return p_mmfail(n);
    }
    arena->next_free = block->next_free;
    block->arena = arena;
    p_nsmall++;
  } else {
    block = (void *)malloc(MM_EXTRA*sizeof(mm_block)+n);
    if (!block) return p_mmfail(n);
    block->arena = 0;
  }
  MM_GUARD(block, n)
  p_nallocs++;
  return block+MM_OFFSET;
}

static void
bfree(void *p)
{
  if (p) {
    mm_block *block = ((mm_block *)p) - MM_OFFSET;
    mm_arena *arena = block->arena;
    MM_VCHECK(p)
    if (arena) {
      block->next_free = arena->next_free;
      arena->next_free = block;
      p_nsmall--;
    } else {
      free(block);
    }
    p_nfrees++;
  }
}

static mm_block *
expand_arena(mm_arena *arena)
{
  long n = arena->n_blocks;
  long size = sizeof(mm_block)*(n+1);
  mm_block *block = (void *)malloc(size);
  if (!block) return p_mmfail(size);
  block[0].next_free = arena->blocks;
  arena->blocks = block++;
  p_asmall += size;
  for (size=arena->size ; n-=size ; block+=size)
    block->next_free = block + size;
  block->next_free = 0;
  return (arena->next_free = arena->blocks+1);
}

static void *
brealloc(void *p, size_t n)
{
  if (p) {
    mm_block *block = ((mm_block *)p) - MM_OFFSET;
    mm_arena *arena = block->arena;
    MM_CHECK(p, n<=0? 1 : n)
    if (arena) {
      mm_block *new_block;
      long i, old_n = arena->size-MM_EXTRA;
      if (n <= sizeof(mm_block)*old_n) {
        MM_GUARD(block, n)
        return p;
      }
      new_block = bmalloc(n);
      if (!new_block) return p_mmfail(n<=0? 1 : n);
      for (i=0 ; i<old_n ; i++) new_block[i] = block[i+MM_OFFSET];
      block->next_free = arena->next_free;
      arena->next_free = block;
      p_nfrees++;
      return new_block;
    } else {
      /* don't bother trying to put shrinking blocks back into
       * the small block arenas */
      if (n<=0) n = 1;
      block = (void *)realloc(block, MM_EXTRA*sizeof(mm_block)+n);
      if (!block) return p_mmfail(n);
      MM_GUARD(block, n)
      return block+MM_OFFSET;
    }
  } else {
    return bmalloc(n);
  }
}

#ifdef P_DEBUG
/* provide primitive checking for memory overrun and underrun
 * checked on every free or realloc operation, or p_mmcheck can
 * be called anytime by debugger or debugging code */

static mm_block ur_pattern;
static char *pattern = 0;

void
p_mmguard(void *b, unsigned long n)
{
  mm_block *block = b;
  char *p = (char *)&block[2];
  int i;
  if (!pattern) {
    pattern = (char *)&ur_pattern;
    for (i=0 ; i<sizeof(mm_block) ; i++) pattern[i] = 255-i;
  }
  block[1].l = n;
  for (i=0 ; i<sizeof(mm_block) ; i++) p[i] = pattern[i];
  p += sizeof(mm_block)+n;
  for (i=0 ; i<sizeof(mm_block) ; i++) p[i] = pattern[i];
}

int
p_mmcheck(void *p)
{
  mm_block *block = p;
  char *chk = (char *)&block[-1];
  unsigned long n = block[-2].l;
  int i;
  for (i=0 ; i<sizeof(mm_block) ; i++) if (chk[i]!=pattern[i]) break;
  if (i<sizeof(mm_block)) return 1;  /* signal memory underrun */
  chk += sizeof(mm_block)+n;
  for (i=0 ; i<sizeof(mm_block) ; i++) if (chk[i]!=pattern[i]) break;
  if (i<sizeof(mm_block)) return 2;  /* signal memory overrun */
  return 0;
}
#endif
