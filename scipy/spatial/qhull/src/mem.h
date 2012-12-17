/*<html><pre>  -<a                             href="qh-mem.htm"
  >-------------------------------</a><a name="TOP">-</a>

   mem.h
     prototypes for memory management functions

   see qh-mem.htm, mem.c and qset.h

   for error handling, writes message and calls
     qh_errexit(qhmem_ERRmem, NULL, NULL) if insufficient memory
       and
     qh_errexit(qhmem_ERRqhull, NULL, NULL) otherwise

   Copyright (c) 1993-2012 The Geometry Center.
   $Id: //main/2011/qhull/src/libqhull/mem.h#4 $$Change: 1464 $
   $DateTime: 2012/01/25 22:58:41 $$Author: bbarber $
*/

#ifndef qhDEFmem
#define qhDEFmem 1

#include <stdio.h>

/*-<a                             href="qh-mem.htm#TOC"
  >-------------------------------</a><a name="NOmem">-</a>

  qh_NOmem
    turn off quick-fit memory allocation

  notes:
    mem.c implements Quickfit memory allocation for about 20% time
    savings.  If it fails on your machine, try to locate the
    problem, and send the answer to qhull@qhull.org.  If this can
    not be done, define qh_NOmem to use malloc/free instead.

   #define qh_NOmem
*/

/*-<a                             href="qh-mem.htm#TOC"
>-------------------------------</a><a name="TRACEshort">-</a>

qh_TRACEshort
Trace short and quick memory allocations at T5

*/
#define qh_TRACEshort

/*-------------------------------------------
    to avoid bus errors, memory allocation must consider alignment requirements.
    malloc() automatically takes care of alignment.   Since mem.c manages
    its own memory, we need to explicitly specify alignment in
    qh_meminitbuffers().

    A safe choice is sizeof(double).  sizeof(float) may be used if doubles
    do not occur in data structures and pointers are the same size.  Be careful
    of machines (e.g., DEC Alpha) with large pointers.  If gcc is available,
    use __alignof__(double) or fmax_(__alignof__(float), __alignof__(void *)).

   see <a href="user.h#MEMalign">qh_MEMalign</a> in user.h for qhull's alignment
*/

#define qhmem_ERRmem 4    /* matches qh_ERRmem in libqhull.h */
#define qhmem_ERRqhull 5  /* matches qh_ERRqhull in libqhull.h */

/*-<a                             href="qh-mem.htm#TOC"
  >--------------------------------</a><a name="ptr_intT">-</a>

  ptr_intT
    for casting a void * to an integer-type that holds a pointer
    Used for integer expressions (e.g., computing qh_gethash() in poly.c)

  notes:
    WARN64 -- these notes indicate 64-bit issues
    On 64-bit machines, a pointer may be larger than an 'int'.
    qh_meminit()/mem.c checks that 'ptr_intT' holds a 'void*'
    ptr_intT is typically a signed value, but not necessarily so
    size_t is typically unsigned, but should match the parameter type
    Qhull uses int instead of size_t except for system calls such as malloc, qsort, qh_malloc, etc.
    This matches Qt convention and is easier to work with.  
*/
#if _MSC_VER && defined(_WIN64)
typedef long long ptr_intT;
#else
typedef long ptr_intT;
#endif

/*-<a                             href="qh-mem.htm#TOC"
  >--------------------------------</a><a name="qhmemT">-</a>

  qhmemT
    global memory structure for mem.c

 notes:
   users should ignore qhmem except for writing extensions
   qhmem is allocated in mem.c

   qhmem could be swapable like qh and qhstat, but then
   multiple qh's and qhmem's would need to keep in synch.
   A swapable qhmem would also waste memory buffers.  As long
   as memory operations are atomic, there is no problem with
   multiple qh structures being active at the same time.
   If you need separate address spaces, you can swap the
   contents of qhmem.
*/
typedef struct qhmemT qhmemT;
extern qhmemT qhmem;

#ifndef DEFsetT
#define DEFsetT 1
typedef struct setT setT;          /* defined in qset.h */
#endif

/* Update qhmem in mem.c if add or remove fields */
struct qhmemT {               /* global memory management variables */
  int      BUFsize;           /* size of memory allocation buffer */
  int      BUFinit;           /* initial size of memory allocation buffer */
  int      TABLEsize;         /* actual number of sizes in free list table */
  int      NUMsizes;          /* maximum number of sizes in free list table */
  int      LASTsize;          /* last size in free list table */
  int      ALIGNmask;         /* worst-case alignment, must be 2^n-1 */
  void   **freelists;          /* free list table, linked by offset 0 */
  int     *sizetable;         /* size of each freelist */
  int     *indextable;        /* size->index table */
  void    *curbuffer;         /* current buffer, linked by offset 0 */
  void    *freemem;           /*   free memory in curbuffer */
  int      freesize;          /*   size of freemem in bytes */
  setT    *tempstack;         /* stack of temporary memory, managed by users */
  FILE    *ferr;              /* file for reporting errors, only user is qh_fprintf() */
  int      IStracing;         /* =5 if tracing memory allocations */
  int      cntquick;          /* count of quick allocations */
                              /* Note: removing statistics doesn't effect speed */
  int      cntshort;          /* count of short allocations */
  int      cntlong;           /* count of long allocations */
  int      freeshort;         /* count of short memfrees */
  int      freelong;          /* count of long memfrees */
  int      totbuffer;         /* total short memory buffers minus buffer links */
  int      totdropped;        /* total dropped memory at end of short memory buffers (e.g., freesize) */
  int      totfree;           /* total size of free, short memory on freelists */
  int      totlong;           /* total size of long memory in use */
  int      maxlong;           /*   maximum totlong */
  int      totshort;          /* total size of short memory in use */
  int      totunused;         /* total unused short memory (estimated, short size - request size of first allocations) */
  int      cntlarger;         /* count of setlarger's */
  int      totlarger;         /* total copied by setlarger */
};


/*==================== -macros ====================*/

/*-<a                             href="qh-mem.htm#TOC"
  >--------------------------------</a><a name="memalloc_">-</a>

  qh_memalloc_(insize, object, type)
    returns object of size bytes
        assumes size<=qhmem.LASTsize and void **freelistp is a temp
*/

#if defined qh_NOmem
#define qh_memalloc_(insize, freelistp, object, type) {\
  object= (type*)qh_memalloc(insize); }
#elif defined qh_TRACEshort
#define qh_memalloc_(insize, freelistp, object, type) {\
    freelistp= NULL; /* Avoid warnings */ \
    object= (type*)qh_memalloc(insize); }
#else /* !qh_NOmem */

#define qh_memalloc_(insize, freelistp, object, type) {\
  freelistp= qhmem.freelists + qhmem.indextable[insize];\
  if ((object= (type*)*freelistp)) {\
    qhmem.totshort += qhmem.sizetable[qhmem.indextable[insize]]; \
    qhmem.totfree -= qhmem.sizetable[qhmem.indextable[insize]]; \
    qhmem.cntquick++;  \
    *freelistp= *((void **)*freelistp);\
  }else object= (type*)qh_memalloc(insize);}
#endif

/*-<a                             href="qh-mem.htm#TOC"
  >--------------------------------</a><a name="memfree_">-</a>

  qh_memfree_(object, insize)
    free up an object

  notes:
    object may be NULL
    assumes size<=qhmem.LASTsize and void **freelistp is a temp
*/
#if defined qh_NOmem
#define qh_memfree_(object, insize, freelistp) {\
  qh_memfree(object, insize); }
#elif defined qh_TRACEshort
#define qh_memfree_(object, insize, freelistp) {\
    freelistp= NULL; /* Avoid warnings */ \
    qh_memfree(object, insize); }
#else /* !qh_NOmem */

#define qh_memfree_(object, insize, freelistp) {\
  if (object) { \
    qhmem .freeshort++;\
    freelistp= qhmem.freelists + qhmem.indextable[insize];\
    qhmem.totshort -= qhmem.sizetable[qhmem.indextable[insize]]; \
    qhmem.totfree += qhmem.sizetable[qhmem.indextable[insize]]; \
    *((void **)object)= *freelistp;\
    *freelistp= object;}}
#endif

/*=============== prototypes in alphabetical order ============*/

void *qh_memalloc(int insize);
void qh_memfree(void *object, int insize);
void qh_memfreeshort(int *curlong, int *totlong);
void qh_meminit(FILE *ferr);
void qh_meminitbuffers(int tracelevel, int alignment, int numsizes,
                        int bufsize, int bufinit);
void qh_memsetup(void);
void qh_memsize(int size);
void qh_memstatistics(FILE *fp);
void qh_memtotal(int *totlong, int *curlong, int *totshort, int *curshort, int *maxlong, int *totbuffer);

#endif /* qhDEFmem */
