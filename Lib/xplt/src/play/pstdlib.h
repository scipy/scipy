/*
 * pstdlib.h -- $Id$
 * portability layer basic memory management interface
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include <stdlib.h>

#include "extern_c.h"

extern void *(*p_malloc)(size_t);
extern void  (*p_free)(void *);
extern void *(*p_realloc)(void *, size_t);

/* above data loaded to system malloc, free, and realloc
 * -- call p_mminit to get mm version
 */
#ifdef P_DEBUG
#define p_mminit p_mmdebug
extern int p_mmcheck(void *p);
extern void p_mmguard(void *b, unsigned long n);
extern long p_mmextra, p_mmoffset;
#endif
extern void p_mminit(void);

/* make trivial memory statistics globally available
 * -- counts total number of allocations, frees, and
 *    current number of large blocks */
extern long p_nallocs;
extern long p_nfrees;
extern long p_nsmall;
extern long p_asmall;

/* define this to get control when mm functions fail
 * -- if it returns, must return 0 */
extern void *(*p_mmfail)(unsigned long n);

/* temporary space */
#define P_WKSIZ 2048
typedef union {
  char c[P_WKSIZ+8];
  int i[P_WKSIZ/8];
  long l[P_WKSIZ/8];
  double d[P_WKSIZ/8];
} p_twkspc;
extern p_twkspc p_wkspc;

/* similar to the string.h functions, but p_malloc destination
 * - 0 src is acceptable */
extern void *p_memcpy(const void *, size_t);
extern char *p_strcpy(const char *);
extern char *p_strncat(const char *, const char *, size_t);

/* dont do anything critical if this is set -- signal an error */
extern volatile int p_signalling;

END_EXTERN_C
