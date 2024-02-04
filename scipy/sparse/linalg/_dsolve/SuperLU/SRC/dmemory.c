/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file dmemory.c
 * \brief Memory details
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory.
 * June 30, 2009
 * </pre>
 */
#include "slu_ddefs.h"


/* Internal prototypes */
void  *dexpand (int_t *, MemType, int_t, int, GlobalLU_t *);
int   dLUWorkInit (int, int, int, int **, double **, GlobalLU_t *);
void  copy_mem_double (int_t, void *, void *);
void  dStackCompress (GlobalLU_t *);
void  dSetupSpace (void *, int_t, GlobalLU_t *);
void  *duser_malloc (int, int, GlobalLU_t *);
void  duser_free (int, int, GlobalLU_t *);

/* External prototypes (in memory.c - prec-independent) */
extern void    copy_mem_int    (int, void *, void *);
extern void    user_bcopy      (char *, char *, int);


/* Macros to manipulate stack */
#define StackFull(x)         ( x + Glu->stack.used >= Glu->stack.size )
#define NotDoubleAlign(addr) ( (intptr_t)addr & 7 )
#define DoubleAlign(addr)    ( ((intptr_t)addr + 7) & ~7L )	
#define TempSpace(m, w)      ( (2*w + 4 + NO_MARKER) * m * sizeof(int) + \
			      (w + 1) * m * sizeof(double) )
#define Reduce(alpha)        ((alpha + 1) / 2)  /* i.e. (alpha-1)/2 + 1 */




/*! \brief Setup the memory model to be used for factorization.
 *  
 *    lwork = 0: use system malloc;
 *    lwork > 0: use user-supplied work[] space.
 */
void dSetupSpace(void *work, int_t lwork, GlobalLU_t *Glu)
{
    if ( lwork == 0 ) {
	Glu->MemModel = SYSTEM; /* malloc/free */
    } else if ( lwork > 0 ) {
	Glu->MemModel = USER;   /* user provided space */
	Glu->stack.used = 0;
	Glu->stack.top1 = 0;
	Glu->stack.top2 = (lwork/4)*4; /* must be word addressable */
	Glu->stack.size = Glu->stack.top2;
	Glu->stack.array = (void *) work;
    }
}



void *duser_malloc(int bytes, int which_end, GlobalLU_t *Glu)
{
    void *buf;
    
    if ( StackFull(bytes) ) return (NULL);

    if ( which_end == HEAD ) {
	buf = (char*) Glu->stack.array + Glu->stack.top1;
	Glu->stack.top1 += bytes;
    } else {
	Glu->stack.top2 -= bytes;
	buf = (char*) Glu->stack.array + Glu->stack.top2;
    }
    
    Glu->stack.used += bytes;
    return buf;
}


void duser_free(int bytes, int which_end, GlobalLU_t *Glu)
{
    if ( which_end == HEAD ) {
	Glu->stack.top1 -= bytes;
    } else {
	Glu->stack.top2 += bytes;
    }
    Glu->stack.used -= bytes;
}



/*!
 * Calculate memory usage
 *
 * \param mem_usage consists of the following fields:
 *    - <tt>for_lu (float)</tt>
 *      The amount of space used in bytes for the L\\U data structures.
 *    - <tt>total_needed (float)</tt>
 *      The amount of space needed in bytes to perform factorization.
 */
int dQuerySpace(SuperMatrix *L, SuperMatrix *U, mem_usage_t *mem_usage)
{
    SCformat *Lstore;
    NCformat *Ustore;
    register int n, iword, dword, panel_size = sp_ienv(1);

    Lstore = L->Store;
    Ustore = U->Store;
    n = L->ncol;
    iword = sizeof(int);
    dword = sizeof(double);

    /* For LU factors */
    mem_usage->for_lu = (float)( (4.0*n + 3.0) * iword +
                                 Lstore->nzval_colptr[n] * dword +
                                 Lstore->rowind_colptr[n] * iword );
    mem_usage->for_lu += (float)( (n + 1.0) * iword +
				 Ustore->colptr[n] * (dword + iword) );

    /* Working storage to support factorization */
    mem_usage->total_needed = mem_usage->for_lu +
	(float)( (2.0 * panel_size + 4.0 + NO_MARKER) * n * iword +
		(panel_size + 1.0) * n * dword );

    return 0;
} /* dQuerySpace */


/*!
 * Calculate memory usage
 *
 * \param mem_usage consists of the following fields:
 *    - <tt>for_lu (float)</tt>
 *      The amount of space used in bytes for the L\\U data structures.
 *    - <tt>total_needed (float)</tt>
 *      The amount of space needed in bytes to perform factorization.
 *
 */
int ilu_dQuerySpace(SuperMatrix *L, SuperMatrix *U, mem_usage_t *mem_usage)
{
    SCformat *Lstore;
    NCformat *Ustore;
    register int n, panel_size = sp_ienv(1);
    register float iword, dword;

    Lstore = L->Store;
    Ustore = U->Store;
    n = L->ncol;
    iword = sizeof(int);
    dword = sizeof(double);

    /* For LU factors */
    mem_usage->for_lu = (float)( (4.0f * n + 3.0f) * iword +
				 Lstore->nzval_colptr[n] * dword +
				 Lstore->rowind_colptr[n] * iword );
    mem_usage->for_lu += (float)( (n + 1.0f) * iword +
				 Ustore->colptr[n] * (dword + iword) );

    /* Working storage to support factorization.
       ILU needs 5*n more integers than LU */
    mem_usage->total_needed = mem_usage->for_lu +
	(float)( (2.0f * panel_size + 9.0f + NO_MARKER) * n * iword +
		(panel_size + 1.0f) * n * dword );

    return 0;
} /* ilu_dQuerySpace */


/*! \brief Allocate storage for the data structures common to all factor routines.
 *
 * <pre>
 * For those unpredictable size, estimate as fill_ratio * nnz(A).
 * Return value:
 *     If lwork = -1, return the estimated amount of space required, plus n;
 *     otherwise, return the amount of space actually allocated when
 *     memory allocation failure occurred.
 * </pre> 
 */
int_t
dLUMemInit(fact_t fact, void *work, int_t lwork, int m, int n, int_t annz,
	  int panel_size, double fill_ratio, SuperMatrix *L, SuperMatrix *U,
          GlobalLU_t *Glu, int **iwork, double **dwork)
{
    int      info, iword, dword;
    SCformat *Lstore;
    NCformat *Ustore;
    int      *xsup, *supno;
    int_t    *lsub, *xlsub;
    double   *lusup;
    int_t    *xlusup;
    double   *ucol;
    int_t    *usub, *xusub;
    int_t    nzlmax, nzumax, nzlumax;
    
    iword     = sizeof(int);
    dword     = sizeof(double);
    Glu->n    = n;
    Glu->num_expansions = 0;

    Glu->expanders = (ExpHeader *) SUPERLU_MALLOC( NO_MEMTYPE *
                                                     sizeof(ExpHeader) );
    if ( !Glu->expanders ) ABORT("SUPERLU_MALLOC fails for expanders");
    
    if ( fact != SamePattern_SameRowPerm ) {
	/* Guess for L\U factors */
	nzumax = nzlumax = nzlmax = fill_ratio * annz;
	//nzlmax = SUPERLU_MAX(1, fill_ratio/4.) * annz;

	if ( lwork == -1 ) {
	    return ( GluIntArray(n) * iword + TempSpace(m, panel_size)
		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword + n );
        } else {
	    dSetupSpace(work, lwork, Glu);
	}
	
#if ( PRNTlevel >= 1 )
	printf("dLUMemInit() called: fill_ratio %.0f, nzlmax %lld, nzumax %lld\n", 
	       fill_ratio, (long long) nzlmax, (long long) nzumax);
	fflush(stdout);
#endif	
	
	/* Integer pointers for L\U factors */
	if ( Glu->MemModel == SYSTEM ) {
	    xsup   = int32Malloc(n+1);
	    supno  = int32Malloc(n+1);
	    xlsub  = intMalloc(n+1);
	    xlusup = intMalloc(n+1);
	    xusub  = intMalloc(n+1);
	} else {
	    xsup   = (int *)duser_malloc((n+1) * iword, HEAD, Glu);
	    supno  = (int *)duser_malloc((n+1) * iword, HEAD, Glu);
	    xlsub  = duser_malloc((n+1) * iword, HEAD, Glu);
	    xlusup = duser_malloc((n+1) * iword, HEAD, Glu);
	    xusub  = duser_malloc((n+1) * iword, HEAD, Glu);
	}

	lusup = (double *) dexpand( &nzlumax, LUSUP, 0, 0, Glu );
	ucol  = (double *) dexpand( &nzumax, UCOL, 0, 0, Glu );
	lsub  = (int_t *) dexpand( &nzlmax, LSUB, 0, 0, Glu );
	usub  = (int_t *) dexpand( &nzumax, USUB, 0, 1, Glu );

	while ( !lusup || !ucol || !lsub || !usub ) {
	    if ( Glu->MemModel == SYSTEM ) {
		SUPERLU_FREE(lusup); 
		SUPERLU_FREE(ucol); 
		SUPERLU_FREE(lsub); 
		SUPERLU_FREE(usub);
	    } else {
		duser_free((nzlumax+nzumax)*dword+(nzlmax+nzumax)*iword,
                            HEAD, Glu);
	    }
	    nzlumax /= 2;
	    nzumax /= 2;
	    nzlmax /= 2;
	    if ( nzlumax < annz ) {
		printf("Not enough memory to perform factorization.\n");
		return (dmemory_usage(nzlmax, nzumax, nzlumax, n) + n);
	    }
#if ( PRNTlevel >= 1)
	    printf("dLUMemInit() reduce size: nzlmax %ld, nzumax %ld\n", 
		   (long) nzlmax, (long) nzumax);
	    fflush(stdout);
#endif
	    lusup = (double *) dexpand( &nzlumax, LUSUP, 0, 0, Glu );
	    ucol  = (double *) dexpand( &nzumax, UCOL, 0, 0, Glu );
	    lsub  = (int_t *) dexpand( &nzlmax, LSUB, 0, 0, Glu );
	    usub  = (int_t *) dexpand( &nzumax, USUB, 0, 1, Glu );
	}
	
    } else {
	/* fact == SamePattern_SameRowPerm */
	Lstore   = L->Store;
	Ustore   = U->Store;
	xsup     = Lstore->sup_to_col;
	supno    = Lstore->col_to_sup;
	xlsub    = Lstore->rowind_colptr;
	xlusup   = Lstore->nzval_colptr;
	xusub    = Ustore->colptr;
	nzlmax   = Glu->nzlmax;    /* max from previous factorization */
	nzumax   = Glu->nzumax;
	nzlumax  = Glu->nzlumax;
	
	if ( lwork == -1 ) {
	    return ( GluIntArray(n) * iword + TempSpace(m, panel_size)
		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword + n );
        } else if ( lwork == 0 ) {
	    Glu->MemModel = SYSTEM;
	} else {
	    Glu->MemModel = USER;
	    Glu->stack.top2 = (lwork/4)*4; /* must be word-addressable */
	    Glu->stack.size = Glu->stack.top2;
	}
	
	lsub  = Glu->expanders[LSUB].mem  = Lstore->rowind;
	lusup = Glu->expanders[LUSUP].mem = Lstore->nzval;
	usub  = Glu->expanders[USUB].mem  = Ustore->rowind;
	ucol  = Glu->expanders[UCOL].mem  = Ustore->nzval;;
	Glu->expanders[LSUB].size         = nzlmax;
	Glu->expanders[LUSUP].size        = nzlumax;
	Glu->expanders[USUB].size         = nzumax;
	Glu->expanders[UCOL].size         = nzumax;	
    }

    Glu->xsup    = xsup;
    Glu->supno   = supno;
    Glu->lsub    = lsub;
    Glu->xlsub   = xlsub;
    Glu->lusup   = (void *) lusup;
    Glu->xlusup  = xlusup;
    Glu->ucol    = (void *) ucol;
    Glu->usub    = usub;
    Glu->xusub   = xusub;
    Glu->nzlmax  = nzlmax;
    Glu->nzumax  = nzumax;
    Glu->nzlumax = nzlumax;
    
    info = dLUWorkInit(m, n, panel_size, iwork, dwork, Glu);
    if ( info )
	return ( info + dmemory_usage(nzlmax, nzumax, nzlumax, n) + n);
    
    ++Glu->num_expansions;
    return 0;
    
} /* dLUMemInit */

/*! \brief Allocate known working storage.
 * Returns 0 if success, otherwise
 * returns the number of bytes allocated so far when failure occurred.
 */
int
dLUWorkInit(int m, int n, int panel_size, int **iworkptr, 
            double **dworkptr, GlobalLU_t *Glu)
{
    int    isize, dsize, extra;
    double *old_ptr;
    int    maxsuper = SUPERLU_MAX( sp_ienv(3), sp_ienv(7) ),
           rowblk   = sp_ienv(4);

    /* xplore[m] and xprune[n] can be 64-bit; they are allocated separately */
    //isize = ( (2 * panel_size + 3 + NO_MARKER ) * m + n ) * sizeof(int);
    isize = ( (2 * panel_size + 2 + NO_MARKER ) * m ) * sizeof(int);
    dsize = (m * panel_size +
	     NUM_TEMPV(m,panel_size,maxsuper,rowblk)) * sizeof(double);
    
    if ( Glu->MemModel == SYSTEM ) 
	*iworkptr = (int *) int32Calloc(isize/sizeof(int));
    else
	*iworkptr = (int *) duser_malloc(isize, TAIL, Glu);
    if ( ! *iworkptr ) {
	fprintf(stderr, "dLUWorkInit: malloc fails for local iworkptr[]\n");
	return (isize + n);
    }

    if ( Glu->MemModel == SYSTEM )
	*dworkptr = (double *) SUPERLU_MALLOC(dsize);
    else {
	*dworkptr = (double *) duser_malloc(dsize, TAIL, Glu);
	if ( NotDoubleAlign(*dworkptr) ) {
	    old_ptr = *dworkptr;
	    *dworkptr = (double*) DoubleAlign(*dworkptr);
	    *dworkptr = (double*) ((double*)*dworkptr - 1);
	    extra = (char*)old_ptr - (char*)*dworkptr;
#if ( DEBUGlevel>=1 )
	    printf("dLUWorkInit: not aligned, extra %d\n", extra); fflush(stdout);
#endif	    
	    Glu->stack.top2 -= extra;
	    Glu->stack.used += extra;
	}
    }
    if ( ! *dworkptr ) {
	fprintf(stderr, "malloc fails for local dworkptr[].");
	return (isize + dsize + n);
    }
	
    return 0;
} /* end dLUWorkInit */


/*! \brief Set up pointers for real working arrays.
 */
void
dSetRWork(int m, int panel_size, double *dworkptr,
	 double **dense, double **tempv)
{
    double zero = 0.0;

    int maxsuper = SUPERLU_MAX( sp_ienv(3), sp_ienv(7) ),
        rowblk   = sp_ienv(4);
    *dense = dworkptr;
    *tempv = *dense + panel_size*m;
    dfill (*dense, m * panel_size, zero);
    dfill (*tempv, NUM_TEMPV(m,panel_size,maxsuper,rowblk), zero);     
}
	
/*! \brief Free the working storage used by factor routines.
 */
void dLUWorkFree(int *iwork, double *dwork, GlobalLU_t *Glu)
{
    if ( Glu->MemModel == SYSTEM ) {
	SUPERLU_FREE (iwork);
	SUPERLU_FREE (dwork);
    } else {
	Glu->stack.used -= (Glu->stack.size - Glu->stack.top2);
	Glu->stack.top2 = Glu->stack.size;
/*	dStackCompress(Glu);  */
    }
    
    SUPERLU_FREE (Glu->expanders);	
    Glu->expanders = NULL;
}

/*! \brief Expand the data structures for L and U during the factorization.
 * 
 * <pre>
 * Return value:   0 - successful return
 *               > 0 - number of bytes allocated when run out of space
 * </pre>
 */
int_t
dLUMemXpand(int jcol,
	   int_t next,          /* number of elements currently in the factors */
	   MemType mem_type,  /* which type of memory to expand  */
	   int_t *maxlen,       /* modified - maximum length of a data structure */
	   GlobalLU_t *Glu    /* modified - global LU data structures */
	   )
{
    void   *new_mem;
    
#if ( DEBUGlevel>=1 ) 
    printf("dLUMemXpand[1]: jcol %d, next %lld, maxlen %lld, MemType %d\n",
	   jcol, (long long) next, (long long) *maxlen, mem_type);
#endif    

    if (mem_type == USUB) 
    	new_mem = dexpand(maxlen, mem_type, next, 1, Glu);
    else
	new_mem = dexpand(maxlen, mem_type, next, 0, Glu);
    
    if ( !new_mem ) {
	int_t    nzlmax  = Glu->nzlmax;
	int_t    nzumax  = Glu->nzumax;
	int_t    nzlumax = Glu->nzlumax;
    	fprintf(stderr, "Can't expand MemType %d: jcol %d\n", mem_type, jcol);
    	return (dmemory_usage(nzlmax, nzumax, nzlumax, Glu->n) + Glu->n);
    }

    switch ( mem_type ) {
      case LUSUP:
	Glu->lusup   = (void *) new_mem;
	Glu->nzlumax = *maxlen;
	break;
      case UCOL:
	Glu->ucol   = (void *) new_mem;
	Glu->nzumax = *maxlen;
	break;
      case LSUB:
	Glu->lsub   = (int_t *) new_mem;
	Glu->nzlmax = *maxlen;
	break;
      case USUB:
	Glu->usub   = (int_t *) new_mem;
	Glu->nzumax = *maxlen;
	break;
      default: break;
    }
    
    return 0;
    
}


void
copy_mem_double(int_t howmany, void *old, void *new)
{
    register int_t i;
    double *dold = old;
    double *dnew = new;
    for (i = 0; i < howmany; i++) dnew[i] = dold[i];
}

/*! \brief Expand the existing storage to accommodate more fill-ins.
 */
void
*dexpand (
	 int_t *prev_len,   /* length used from previous call */
	 MemType type,    /* which part of the memory to expand */
	 int_t len_to_copy, /* size of the memory to be copied to new store */
	 int keep_prev,   /* = 1: use prev_len;
			     = 0: compute new_len to expand */
	 GlobalLU_t *Glu  /* modified - global LU data structures */
	)
{
    float    EXPAND = 1.5;
    float    alpha;
    void     *new_mem, *old_mem;
    int_t    new_len, bytes_to_copy;
    int      tries, lword, extra;
    ExpHeader *expanders = Glu->expanders; /* Array of 4 types of memory */

    alpha = EXPAND;

    if ( Glu->num_expansions == 0 || keep_prev ) {
        /* First time allocate requested */
        new_len = *prev_len;
    } else {
	new_len = alpha * *prev_len;
    }
    
    if ( type == LSUB || type == USUB ) lword = sizeof(int_t);
    else lword = sizeof(double);

    if ( Glu->MemModel == SYSTEM ) {
	new_mem = (void *) SUPERLU_MALLOC((size_t)new_len * lword);
	if ( Glu->num_expansions != 0 ) {
	    tries = 0;
	    if ( keep_prev ) {
		if ( !new_mem ) return (NULL);
	    } else {
		while ( !new_mem ) {
		    if ( ++tries > 10 ) return (NULL);
		    alpha = Reduce(alpha);
		    new_len = alpha * *prev_len;
		    new_mem = (void *) SUPERLU_MALLOC((size_t)new_len * lword);
		}
	    }
	    if ( type == LSUB || type == USUB ) {
		copy_mem_int(len_to_copy, expanders[type].mem, new_mem);
	    } else {
		copy_mem_double(len_to_copy, expanders[type].mem, new_mem);
	    }
	    SUPERLU_FREE (expanders[type].mem);
	}
	expanders[type].mem = (void *) new_mem;
	
    } else { /* MemModel == USER */
    
	if ( Glu->num_expansions == 0 ) { /* First time initialization */
	
	    new_mem = duser_malloc(new_len * lword, HEAD, Glu);
	    if ( NotDoubleAlign(new_mem) &&
		(type == LUSUP || type == UCOL) ) {
		old_mem = new_mem;
		new_mem = (void *)DoubleAlign(new_mem);
		extra = (char*)new_mem - (char*)old_mem;
#if ( DEBUGlevel>=1 )
		printf("expand(): not aligned, extra %d\n", extra);
#endif		
		Glu->stack.top1 += extra;
		Glu->stack.used += extra;
	    }
	    
	    expanders[type].mem = (void *) new_mem;
	    
	} else { /* CASE: num_expansions != 0 */
	
	    tries = 0;
	    extra = (new_len - *prev_len) * lword;
	    if ( keep_prev ) {
		if ( StackFull(extra) ) return (NULL);
	    } else {
		while ( StackFull(extra) ) {
		    if ( ++tries > 10 ) return (NULL);
		    alpha = Reduce(alpha);
		    new_len = alpha * *prev_len;
		    extra = (new_len - *prev_len) * lword;	    
		}
	    }

	      /* Need to expand the memory: moving the content after the current MemType
	      	 to make extra room for the current MemType.
              	 Memory layout: [ LUSUP || UCOL || LSUB || USUB ]
	      */
  	    if ( type != USUB ) {
		new_mem = (void*)((char*)expanders[type + 1].mem + extra);
		bytes_to_copy = (char*)Glu->stack.array + Glu->stack.top1
		    - (char*)expanders[type + 1].mem;
		user_bcopy(expanders[type+1].mem, new_mem, bytes_to_copy);

		if ( type < USUB ) {
		    Glu->usub = expanders[USUB].mem =
			(void*)((char*)expanders[USUB].mem + extra);
		}
		if ( type < LSUB ) {
		    Glu->lsub = expanders[LSUB].mem =
			(void*)((char*)expanders[LSUB].mem + extra);
		}
		if ( type < UCOL ) {
		    Glu->ucol = expanders[UCOL].mem =
			(void*)((char*)expanders[UCOL].mem + extra);
		}
		Glu->stack.top1 += extra;
		Glu->stack.used += extra;
		if ( type == UCOL ) {
		    Glu->stack.top1 += extra;   /* Add same amount for USUB */
		    Glu->stack.used += extra;
		}
		
	    } /* end expansion */

	} /* else ... */
    }

    expanders[type].size = new_len;
    *prev_len = new_len;
    if ( Glu->num_expansions ) ++Glu->num_expansions;
    
    return (void *) expanders[type].mem;
    
} /* dexpand */


/*! \brief Compress the work[] array to remove fragmentation.
 */
void
dStackCompress(GlobalLU_t *Glu)
{
    register int iword, dword, ndim;
    char     *last, *fragment;
    int_t    *ifrom, *ito;
    double   *dfrom, *dto;
    int_t    *xlsub, *lsub, *xusub, *usub, *xlusup;
    double   *ucol, *lusup;
    
    iword = sizeof(int);
    dword = sizeof(double);
    ndim = Glu->n;

    xlsub  = Glu->xlsub;
    lsub   = Glu->lsub;
    xusub  = Glu->xusub;
    usub   = Glu->usub;
    xlusup = Glu->xlusup;
    ucol   = Glu->ucol;
    lusup  = Glu->lusup;
    
    dfrom = ucol;
    dto = (double *)((char*)lusup + xlusup[ndim] * dword);
    copy_mem_double(xusub[ndim], dfrom, dto);
    ucol = dto;

    ifrom = lsub;
    ito = (int_t *) ((char*)ucol + xusub[ndim] * iword);
    copy_mem_int(xlsub[ndim], ifrom, ito);
    lsub = ito;
    
    ifrom = usub;
    ito = (int_t *) ((char*)lsub + xlsub[ndim] * iword);
    copy_mem_int(xusub[ndim], ifrom, ito);
    usub = ito;
    
    last = (char*)usub + xusub[ndim] * iword;
    fragment = (char*) (((char*)Glu->stack.array + Glu->stack.top1) - last);
    Glu->stack.used -= (long int) fragment;
    Glu->stack.top1 -= (long int) fragment;

    Glu->ucol = ucol;
    Glu->lsub = lsub;
    Glu->usub = usub;
    
#if ( DEBUGlevel>=1 )
    printf("dStackCompress: fragment %lld\n", (long long) fragment);
    /* for (last = 0; last < ndim; ++last)
	print_lu_col("After compress:", last, 0);*/
#endif    
    
}

/*! \brief Allocate storage for original matrix A
 */
void
dallocateA(int n, int_t nnz, double **a, int_t **asub, int_t **xa)
{
    *a    = (double *) doubleMalloc(nnz);
    *asub = (int_t *) intMalloc(nnz);
    *xa   = (int_t *) intMalloc(n+1);
}


double *doubleMalloc(size_t n)
{
    double *buf;
    buf = (double *) SUPERLU_MALLOC(n * (size_t) sizeof(double)); 
    if ( !buf ) {
	ABORT("SUPERLU_MALLOC failed for buf in doubleMalloc()\n");
    }
    return (buf);
}

double *doubleCalloc(size_t n)
{
    double *buf;
    register size_t i;
    double zero = 0.0;
    buf = (double *) SUPERLU_MALLOC(n * (size_t) sizeof(double));
    if ( !buf ) {
	ABORT("SUPERLU_MALLOC failed for buf in doubleCalloc()\n");
    }
    for (i = 0; i < n; ++i) buf[i] = zero;
    return (buf);
}


int_t dmemory_usage(const int_t nzlmax, const int_t nzumax,
		  const int_t nzlumax, const int n)
{
    register int iword, liword, dword;

    iword   = sizeof(int);
    liword  = sizeof(int_t);
    dword   = sizeof(double);
    
    return (10 * n * iword +
	    nzlmax * liword + nzumax * (liword + dword) + nzlumax * dword);

}
