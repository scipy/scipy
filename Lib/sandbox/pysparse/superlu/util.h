#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Macros */
#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#define USER_MALLOC(size) superlu_malloc(size)
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)


#define SUPERLU_MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define SUPERLU_MIN(x, y) 	( (x) < (y) ? (x) : (y) )

/* 
 * Constants 
 */
#define EMPTY	(-1)
#define NO	(-1)
#define FALSE	0
#define TRUE	1

/*
 * Type definitions
 */
typedef float    flops_t;
typedef unsigned char Logical;

/* 
 * The following enumerate type is used by the statistics variable 
 * SuperLUStat, to keep track of flop count and time spent at various stages.
 *
 * Note that not all of the fields are disjoint.
 */
typedef enum {
    COLPERM, /* find a column ordering that minimizes fills */
    RELAX,   /* find artificial supernodes */
    ETREE,   /* compute column etree */
    EQUIL,   /* equilibrate the original matrix */
    FACT,    /* perform LU factorization */
    RCOND,   /* estimate reciprocal condition number */
    SOLVE,   /* forward and back solves */
    REFINE,  /* perform iterative refinement */
    FLOAT,   /* time spent in floating-point operations */
    TRSV,    /* fraction of FACT spent in xTRSV */
    GEMV,    /* fraction of FACT spent in xGEMV */
    FERR,    /* estimate error bounds after iterative refinement */
    NPHASES  /* total number of phases */
} PhaseType;

typedef struct {
    int     *panel_histo; /* histogram of panel size distribution */
    double  *utime;       /* running time at various phases */
    flops_t *ops;         /* operation count at various phases */
} SuperLUStat_t;

/* Macros */
#define FIRSTCOL_OF_SNODE(i)	(xsup[i])


#ifdef __cplusplus
extern "C" {
#endif

extern void    superlu_abort_and_exit(char*);
extern void    *superlu_malloc (size_t);
extern int     *intMalloc (int);
extern int     *intCalloc (int);
extern void    superlu_free (void*);
extern void    SetIWork (int, int, int, int *, int **, int **, int **,
                         int **, int **, int **, int **);
extern void    StatInit(int, int);
extern void    StatFree();
extern int     sp_coletree (int *, int *, int *, int, int, int *);
extern void    relax_snode (const int, int *, const int, int *, int *);
extern void    resetrep_col (const int, const int *, int *);
extern int     spcoletree (int *, int *, int *, int, int, int *);
extern int     *TreePostorder (int, int *);
extern double  SuperLU_timer_ ();
extern int     sp_ienv (int);
extern int     lsame_ (char *, char *);
extern int     xerbla_ (char *, int *);
extern void    ifill (int *, int, int);
extern void    snode_profile (int, int *);
extern void    super_stats (int, int *);
extern void    PrintSumm (char *, int, int, int);
extern void    PrintStat (SuperLUStat_t *);
extern void    print_panel_seg(int, int, int, int, int *, int *);
extern void    check_repfnz(int, int, int, int *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_UTIL */
