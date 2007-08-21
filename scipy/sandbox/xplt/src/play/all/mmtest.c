/*
 * mmtest.c -- $Id$
 * test pstdlib.h memory manager
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "pstdlib.h"
#include <string.h>

union mm_block_fake {
  void *next_free, *arena;
  long l;
  double d;
  void (*f)(void);
};

extern long le_random(long scale);

static void *mmfail(void);
static long fail_count = 0;

int
main(int argc, char *argv[])
{
  long n = argc>1? strtol(argv[1], (char **)0, 10) : 0;
  long mx = argc>2? strtol(argv[2], (char **)0, 10) : 0;
  long count = argc>3? strtol(argv[3], (char **)0, 10) : 0;
  long i, len = 0;
  void **ptrs;

  if (!n) n = 10000;
  else if (n<64) n = 64;
  if (mx<34*sizeof(union mm_block_fake)) mx = 34*sizeof(union mm_block_fake);
  if (!count) count = 20;

  if (count>0) p_mminit();
  else count = -count;
  p_mmfail = &mmfail;

  ptrs = p_malloc(sizeof(void *)*n);

  /* first malloc all n and set their values */
  for (i=0 ; i<n ; i++) {
    len = le_random(mx-5) + 6;
    ptrs[i] = p_malloc(len);
    memset(ptrs[i], 0x41, len);
  }

  /* next free half at random */
  for (i=n ; i>n/2 ; ) {
    len = le_random(n);
    if (len<0 || len>=n || !ptrs[len]) continue;
    p_free(ptrs[len]);
    ptrs[len] = 0;
    i--;
  }

  /* next do realloc, malloc and free in random order */
  while ((count--) > 0) {
    for (i=0 ; i<n ; i++) {
      if (i>5 && ptrs[i-5] && !ptrs[i-1] && i<n-1 && ptrs[i+1]) {
        len = le_random(mx-5) + 6;
        ptrs[i] = p_realloc(ptrs[i], len);
        memset(ptrs[i], 0x51, len);
      } else if (!ptrs[i]) {
        len = le_random(mx-5) + 6;
        ptrs[i] = p_malloc(len);
        memset(ptrs[i], 0x61, len);
      } else {
        p_free(ptrs[i]);
        ptrs[i] = 0;
      }
    }
  }

  /* finally free all remaining ones */
  for (i=0 ; i<n ; i++) p_free(ptrs[i]);

  p_free(ptrs);

  return fail_count!=0;
}

static void *
mmfail(void)
{
  fail_count++;
  return 0;
}

static unsigned long s1 = 0x9ad5e0d7;
static unsigned long s2 = 0xf08ff6d8;
static unsigned long s3 = 0xe7ed3a47;

long
le_random(long scale)
{
  unsigned long b = ((s1 << 13) ^ s1) >> 19;
  s1 = ((s1 & 0xfffffffe) << 12) ^ b;
  b = ((s2 << 2) ^ s2) >> 25;
  s2 = ((s2 & 0xfffffff8) << 4) ^ b;
  b = ((s3 << 3) ^ s3) >> 11;
  s3 = ((s3 & 0xfffffff0) << 17) ^ b;
  b = s1 ^ s2 ^ s3;
  return (long)(2.3283064370807974e-10*(b-0.5001)*scale);
}
