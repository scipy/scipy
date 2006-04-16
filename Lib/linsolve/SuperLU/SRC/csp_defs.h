
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#ifndef __SUPERLU_cSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_cSP_DEFS

/*
 * File name:		csp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#include <string.h>
#endif

/* Define my integer type int_t */
typedef int int_t; /* default */

#include "Cnames.h"
#include "supermatrix.h"
#include "util.h"
#include "scomplex.h"


/*
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 */
typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    complex  *lusup;   /* L supernodes */
    int     *xlusup;
    complex  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} GlobalLU_t;

typedef struct {
    float for_lu;
    float total_needed;
    int   expansions;
} mem_usage_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Driver routines */
extern void
cgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
cgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       mem_usage_t *, SuperLUStat_t *, int *);

/* Supernodal LU factor related */
extern void
cCreate_CompCol_Matrix(SuperMatrix *, int, int, int, complex *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_CompRow_Matrix(SuperMatrix *, int, int, int, complex *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
cCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
cCreate_Dense_Matrix(SuperMatrix *, int, int, complex *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, complex *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
cCopy_Dense_Matrix(int, int, complex *, int, complex *, int);

extern void    countnz (const int, int *, int *, int *, GlobalLU_t *);
extern void    fixupL (const int, const int *, GlobalLU_t *);

extern void    callocateA (int, int, complex **, int **, int **);
extern void    cgstrf (superlu_options_t*, SuperMatrix*, float, 
                       int, int, int*, void *, int, int *, int *, 
                       SuperMatrix *, SuperMatrix *, SuperLUStat_t*, int *);
extern int     csnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, GlobalLU_t *);
extern int     csnode_bmod (const int, const int, const int, complex *,
                              complex *, GlobalLU_t *, SuperLUStat_t*);
extern void    cpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, complex *, int *, int *, int *,
			   int *, int *, int *, int *, GlobalLU_t *);
extern void    cpanel_bmod (const int, const int, const int, const int,
                           complex *, complex *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int     ccolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, GlobalLU_t *);
extern int     ccolumn_bmod (const int, const int, complex *,
			   complex *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
extern int     ccopy_to_ucol (int, int, int *, int *, int *,
                              complex *, GlobalLU_t *);         
extern int     cpivotL (const int, const float, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
extern void    cpruneL (const int, const int *, const int, const int,
			  const int *, const int *, int *, GlobalLU_t *);
extern void    creadmt (int *, int *, int *, complex **, int **, int **);
extern void    cGenXtrue (int, int, complex *, int);
extern void    cFillRHS (trans_t, int, complex *, int, SuperMatrix *,
			  SuperMatrix *);
extern void    cgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                        SuperMatrix *, SuperLUStat_t*, int *);


/* Driver related */

extern void    cgsequ (SuperMatrix *, float *, float *, float *,
			float *, float *, int *);
extern void    claqgs (SuperMatrix *, float *, float *, float,
                        float, float, char *);
extern void    cgscon (char *, SuperMatrix *, SuperMatrix *, 
		         float, float *, SuperLUStat_t*, int *);
extern float   cPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
extern void    cgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, float *, 
                       float *, SuperMatrix *, SuperMatrix *,
                       float *, float *, SuperLUStat_t*, int *);

extern int     sp_ctrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, complex *, SuperLUStat_t*, int *);
extern int     sp_cgemv (char *, complex, SuperMatrix *, complex *,
			int, complex, complex *, int);

extern int     sp_cgemm (char *, char *, int, int, int, complex,
			SuperMatrix *, complex *, int, complex, 
			complex *, int);

/* Memory-related */
extern int     cLUMemInit (fact_t, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     GlobalLU_t *, int **, complex **);
extern void    cSetRWork (int, int, complex *, complex **, complex **);
extern void    cLUWorkFree (int *, complex *, GlobalLU_t *);
extern int     cLUMemXpand (int, int, MemType, int *, GlobalLU_t *);

extern complex  *complexMalloc(int);
extern complex  *complexCalloc(int);
extern float  *floatMalloc(int);
extern float  *floatCalloc(int);
extern int     cmemory_usage(const int, const int, const int, const int);
extern int     cQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/* Auxiliary routines */
extern void    creadhb(int *, int *, int *, complex **, int **, int **);
extern void    cCompRow_to_CompCol(int, int, int, complex*, int*, int*,
		                   complex **, int **, int **);
extern void    cfill (complex *, int, complex);
extern void    cinf_norm_error (int, SuperMatrix *, complex *);
extern void    PrintPerf (SuperMatrix *, SuperMatrix *, mem_usage_t *,
			 complex, complex, complex *, complex *, char *);

/* Routines for debugging */
extern void    cPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    cPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    cPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    print_lu_col(char *, int, int, int *, GlobalLU_t *);
extern void    check_tempv(int, complex *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_cSP_DEFS */

