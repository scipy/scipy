/* 
 * lu_internal.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Functions to load/save/reset struct lu objects
 *
 */

#include "lu_internal.h"

/* private entries in xstore */
#define BASICLU_TASK 256
#define BASICLU_FTCOLUMN_IN 257
#define BASICLU_FTCOLUMN_OUT 258
#define BASICLU_PIVOT_ROW 259
#define BASICLU_PIVOT_COL 260
#define BASICLU_RANKDEF 261
#define BASICLU_MIN_COLNZ 262
#define BASICLU_MIN_ROWNZ 263
#define BASICLU_MARKER 266
#define BASICLU_UPDATE_COST_NUMER 267
#define BASICLU_UPDATE_COST_DENOM 268
#define BASICLU_PIVOTLEN 269


/* ==========================================================================
 * lu_load
 *
 * Initialize @this from @istore, @xstore if these are a valid BASICLU
 * instance. The remaining arguments are copied only and can be NULL.
 *
 * Return BASICLU_OK or BASICLU_ERROR_invalid_store
 * ========================================================================== */

lu_int lu_load(
    struct lu *this, lu_int *istore, double *xstore, lu_int *Li, double *Lx,
    lu_int *Ui, double *Ux, lu_int *Wi, double *Wx)
{
    lu_int m, *iptr;
    double *xptr;

    if (!istore || istore[0] != BASICLU_HASH ||
        !xstore || xstore[0] != BASICLU_HASH )
    {
        return BASICLU_ERROR_invalid_store;
    }

    /* user parameters */
    this->Lmem                  = xstore[BASICLU_MEMORYL];
    this->Umem                  = xstore[BASICLU_MEMORYU];
    this->Wmem                  = xstore[BASICLU_MEMORYW];
    this->droptol               = xstore[BASICLU_DROP_TOLERANCE];
    this->abstol                = xstore[BASICLU_ABS_PIVOT_TOLERANCE];
    this->reltol                = xstore[BASICLU_REL_PIVOT_TOLERANCE];
    this->reltol                = fmin(this->reltol, 1.0);
    this->nzbias                = xstore[BASICLU_BIAS_NONZEROS];
    this->maxsearch             = xstore[BASICLU_MAXN_SEARCH_PIVOT];
    this->pad                   = xstore[BASICLU_PAD];
    this->stretch               = xstore[BASICLU_STRETCH];
    this->compress_thres        = xstore[BASICLU_COMPRESSION_THRESHOLD];
    this->sparse_thres          = xstore[BASICLU_SPARSE_THRESHOLD];
    this->search_rows           = xstore[BASICLU_SEARCH_ROWS] != 0;

    /* user readable */
    this->m = m                 = xstore[BASICLU_DIM];
    this->addmemL               = 0;
    this->addmemU               = 0;
    this->addmemW               = 0;

    this->nupdate               = xstore[BASICLU_NUPDATE];
    this->nforrest              = xstore[BASICLU_NFORREST];
    this->nfactorize            = xstore[BASICLU_NFACTORIZE];
    this->nupdate_total         = xstore[BASICLU_NUPDATE_TOTAL];
    this->nforrest_total        = xstore[BASICLU_NFORREST_TOTAL];
    this->nsymperm_total        = xstore[BASICLU_NSYMPERM_TOTAL];
    this->Lnz                   = xstore[BASICLU_LNZ];
    this->Unz                   = xstore[BASICLU_UNZ];
    this->Rnz                   = xstore[BASICLU_RNZ];
    this->min_pivot             = xstore[BASICLU_MIN_PIVOT];
    this->max_pivot             = xstore[BASICLU_MAX_PIVOT];
    this->max_eta               = xstore[BASICLU_MAX_ETA];
    this->update_cost_numer     = xstore[BASICLU_UPDATE_COST_NUMER];
    this->update_cost_denom     = xstore[BASICLU_UPDATE_COST_DENOM];
    this->time_factorize        = xstore[BASICLU_TIME_FACTORIZE];
    this->time_solve            = xstore[BASICLU_TIME_SOLVE];
    this->time_update           = xstore[BASICLU_TIME_UPDATE];
    this->time_factorize_total  = xstore[BASICLU_TIME_FACTORIZE_TOTAL];
    this->time_solve_total      = xstore[BASICLU_TIME_SOLVE_TOTAL];
    this->time_update_total     = xstore[BASICLU_TIME_UPDATE_TOTAL];
    this->Lflops                = xstore[BASICLU_LFLOPS];
    this->Uflops                = xstore[BASICLU_UFLOPS];
    this->Rflops                = xstore[BASICLU_RFLOPS];
    this->condestL              = xstore[BASICLU_CONDEST_L];
    this->condestU              = xstore[BASICLU_CONDEST_U];
    this->normL                 = xstore[BASICLU_NORM_L];
    this->normU                 = xstore[BASICLU_NORM_U];
    this->normestLinv           = xstore[BASICLU_NORMEST_LINV];
    this->normestUinv           = xstore[BASICLU_NORMEST_UINV];
    this->onenorm               = xstore[BASICLU_MATRIX_ONENORM];
    this->infnorm               = xstore[BASICLU_MATRIX_INFNORM];
    this->residual_test         = xstore[BASICLU_RESIDUAL_TEST];

    this->matrix_nz             = xstore[BASICLU_MATRIX_NZ];
    this->rank                  = xstore[BASICLU_RANK];
    this->bump_size             = xstore[BASICLU_BUMP_SIZE];
    this->bump_nz               = xstore[BASICLU_BUMP_NZ];
    this->nsearch_pivot         = xstore[BASICLU_NSEARCH_PIVOT];
    this->nexpand               = xstore[BASICLU_NEXPAND];
    this->ngarbage              = xstore[BASICLU_NGARBAGE];
    this->factor_flops          = xstore[BASICLU_FACTOR_FLOPS];
    this->time_singletons       = xstore[BASICLU_TIME_SINGLETONS];
    this->time_search_pivot     = xstore[BASICLU_TIME_SEARCH_PIVOT];
    this->time_elim_pivot       = xstore[BASICLU_TIME_ELIM_PIVOT];

    this->pivot_error           = xstore[BASICLU_PIVOT_ERROR];

    /* private */
    this->task                  = xstore[BASICLU_TASK];
    this->pivot_row             = xstore[BASICLU_PIVOT_ROW];
    this->pivot_col             = xstore[BASICLU_PIVOT_COL];
    this->ftran_for_update      = xstore[BASICLU_FTCOLUMN_IN];
    this->btran_for_update      = xstore[BASICLU_FTCOLUMN_OUT];
    this->marker                = xstore[BASICLU_MARKER];
    this->pivotlen              = xstore[BASICLU_PIVOTLEN];
    this->rankdef               = xstore[BASICLU_RANKDEF];
    this->min_colnz             = xstore[BASICLU_MIN_COLNZ];
    this->min_rownz             = xstore[BASICLU_MIN_ROWNZ];

    /* aliases to user arrays */
    this->Lindex = Li; this->Lvalue = Lx;
    this->Uindex = Ui; this->Uvalue = Ux;
    this->Windex = Wi; this->Wvalue = Wx;

    /* partition istore for factorize */
    iptr = istore + 1;
    this->colcount_flink        = iptr; iptr += 2*m+2;
    this->colcount_blink        = iptr; iptr += 2*m+2;
    this->rowcount_flink        = iptr; iptr += 2*m+2;
    this->rowcount_blink        = iptr; iptr += 2*m+2;
    this->Wbegin                = iptr; iptr += 2*m+1;
    this->Wend                  = iptr; iptr += 2*m+1;
    this->Wflink                = iptr; iptr += 2*m+1;
    this->Wblink                = iptr; iptr += 2*m+1;
    this->pinv                  = iptr; iptr += m;
    this->qinv                  = iptr; iptr += m;
    this->Lbegin_p              = iptr; iptr += m+1;
    this->Ubegin                = iptr; iptr += m+1;
    this->iwork0                = iptr; iptr += m;

    /* share istore memory for solve/update */
    this->pivotcol              = this->colcount_flink;
    this->pivotrow              = this->colcount_blink;
    this->Rbegin                = this->rowcount_flink;
    this->eta_row               = this->rowcount_flink + m+1;
    this->iwork1                = this->rowcount_blink;
    this->Lbegin                = this->Wbegin + m+1;
    this->Ltbegin               = this->Wend + m+1;
    this->Ltbegin_p             = this->Wflink + m+1;
    this->p                     = this->Wblink + m+1;
    this->pmap                  = this->pinv;
    this->qmap                  = this->qinv;
    this->marked                = this->iwork0;

    /* partition xstore for factorize and update */
    xptr = xstore + 512;
    this->work0                 = xptr; xptr += m;
    this->work1                 = xptr; xptr += m;
    this->col_pivot             = xptr; xptr += m;
    this->row_pivot             = xptr; xptr += m;

    /*
     * Reset @marked if increasing @marker by four causes overflow.
     */
    if (this->marker > LU_INT_MAX-4)
    {
        memset(this->marked, 0, m*sizeof(lu_int));
        this->marker = 0;
    }

    /*
     * One past the final position in @Wend must hold the file size.
     * The file has 2*m lines while factorizing and m lines otherwise.
     */
    if (this->nupdate >= 0)
        this->Wend[m] = this->Wmem;
    else
        this->Wend[2*m] = this->Wmem;

    return BASICLU_OK;
}


/* ==========================================================================
 * lu_save
 *
 * Copy scalar entries (except for user parameters) from @this to @istore,
 * @xstore. Store status code.
 *
 * Return @status
 * ========================================================================== */

lu_int lu_save(
    const struct lu *this, lu_int *istore, double *xstore, lu_int status)
{
    /* user readable */
    xstore[BASICLU_STATUS]                  = status;
    xstore[BASICLU_ADD_MEMORYL]             = this->addmemL;
    xstore[BASICLU_ADD_MEMORYU]             = this->addmemU;
    xstore[BASICLU_ADD_MEMORYW]             = this->addmemW;

    xstore[BASICLU_NUPDATE]                 = this->nupdate;
    xstore[BASICLU_NFORREST]                = this->nforrest;
    xstore[BASICLU_NFACTORIZE]              = this->nfactorize;
    xstore[BASICLU_NUPDATE_TOTAL]           = this->nupdate_total;
    xstore[BASICLU_NFORREST_TOTAL]          = this->nforrest_total;
    xstore[BASICLU_NSYMPERM_TOTAL]          = this->nsymperm_total;
    xstore[BASICLU_LNZ]                     = this->Lnz;
    xstore[BASICLU_UNZ]                     = this->Unz;
    xstore[BASICLU_RNZ]                     = this->Rnz;
    xstore[BASICLU_MIN_PIVOT]               = this->min_pivot;
    xstore[BASICLU_MAX_PIVOT]               = this->max_pivot;
    xstore[BASICLU_MAX_ETA]                 = this->max_eta;
    xstore[BASICLU_UPDATE_COST_NUMER]       = this->update_cost_numer;
    xstore[BASICLU_UPDATE_COST_DENOM]       = this->update_cost_denom;
    xstore[BASICLU_UPDATE_COST]             =
        this->update_cost_numer / this->update_cost_denom;
    xstore[BASICLU_TIME_FACTORIZE]          = this->time_factorize;
    xstore[BASICLU_TIME_SOLVE]              = this->time_solve;
    xstore[BASICLU_TIME_UPDATE]             = this->time_update;
    xstore[BASICLU_TIME_FACTORIZE_TOTAL]    = this->time_factorize_total;
    xstore[BASICLU_TIME_SOLVE_TOTAL]        = this->time_solve_total;
    xstore[BASICLU_TIME_UPDATE_TOTAL]       = this->time_update_total;
    xstore[BASICLU_LFLOPS]                  = this->Lflops;
    xstore[BASICLU_UFLOPS]                  = this->Uflops;
    xstore[BASICLU_RFLOPS]                  = this->Rflops;
    xstore[BASICLU_CONDEST_L]               = this->condestL;
    xstore[BASICLU_CONDEST_U]               = this->condestU;
    xstore[BASICLU_NORM_L]                  = this->normL;
    xstore[BASICLU_NORM_U]                  = this->normU;
    xstore[BASICLU_NORMEST_LINV]            = this->normestLinv;
    xstore[BASICLU_NORMEST_UINV]            = this->normestUinv;
    xstore[BASICLU_MATRIX_ONENORM]          = this->onenorm;
    xstore[BASICLU_MATRIX_INFNORM]          = this->infnorm;
    xstore[BASICLU_RESIDUAL_TEST]           = this->residual_test;

    xstore[BASICLU_MATRIX_NZ]               = this->matrix_nz;
    xstore[BASICLU_RANK]                    = this->rank;
    xstore[BASICLU_BUMP_SIZE]               = this->bump_size;
    xstore[BASICLU_BUMP_NZ]                 = this->bump_nz;
    xstore[BASICLU_NSEARCH_PIVOT]           = this->nsearch_pivot;
    xstore[BASICLU_NEXPAND]                 = this->nexpand;
    xstore[BASICLU_NGARBAGE]                = this->ngarbage;
    xstore[BASICLU_FACTOR_FLOPS]            = this->factor_flops;
    xstore[BASICLU_TIME_SINGLETONS]         = this->time_singletons;
    xstore[BASICLU_TIME_SEARCH_PIVOT]       = this->time_search_pivot;
    xstore[BASICLU_TIME_ELIM_PIVOT]         = this->time_elim_pivot;

    xstore[BASICLU_PIVOT_ERROR]             = this->pivot_error;

    /* private */
    xstore[BASICLU_TASK]                    = this->task;
    xstore[BASICLU_PIVOT_ROW]               = this->pivot_row;
    xstore[BASICLU_PIVOT_COL]               = this->pivot_col;
    xstore[BASICLU_FTCOLUMN_IN]             = this->ftran_for_update;
    xstore[BASICLU_FTCOLUMN_OUT]            = this->btran_for_update;
    xstore[BASICLU_MARKER]                  = this->marker;
    xstore[BASICLU_PIVOTLEN]                = this->pivotlen;
    xstore[BASICLU_RANKDEF]                 = this->rankdef;
    xstore[BASICLU_MIN_COLNZ]               = this->min_colnz;
    xstore[BASICLU_MIN_ROWNZ]               = this->min_rownz;

    return status;
}


/* ==========================================================================
 * lu_reset
 *
 * Reset @this for a new factorization. Invalidate current factorization.
 * ========================================================================== */

void lu_reset(struct lu *this)
{
    /* user readable */
    this->nupdate = -1;         /* invalidate factorization */
    this->nforrest = 0;
    this->Lnz = 0;
    this->Unz = 0;
    this->Rnz = 0;
    this->min_pivot = 0;
    this->max_pivot = 0;
    this->max_eta = 0;
    this->update_cost_numer = 0;
    this->update_cost_denom = 1;
    this->time_factorize = 0;
    this->time_solve = 0;
    this->time_update = 0;
    this->Lflops = 0;
    this->Uflops = 0;
    this->Rflops = 0;
    this->condestL = 0;
    this->condestU = 0;
    this->normL = 0;
    this->normU = 0;
    this->normestLinv = 0;
    this->normestUinv = 0;
    this->onenorm = 0;
    this->infnorm = 0;
    this->residual_test = 0;

    this->matrix_nz = 0;
    this->rank = 0;
    this->bump_size = 0;
    this->bump_nz = 0;
    this->nsearch_pivot = 0;
    this->nexpand = 0;
    this->ngarbage = 0;
    this->factor_flops = 0;
    this->time_singletons = 0;
    this->time_search_pivot = 0;
    this->time_elim_pivot = 0;

    this->pivot_error = 0;

    /* private */
    this->task = NO_TASK;
    this->pivot_row = -1;
    this->pivot_col = -1;
    this->ftran_for_update = -1;
    this->btran_for_update = -1;
    this->marker = 0;
    this->pivotlen = 0;
    this->rankdef = 0;
    this->min_colnz = 1;
    this->min_rownz = 1;

    /*
     * One past the final position in @Wend must hold the file size.
     * The file has 2*m lines during factorization.
     */
    this->Wend[2*this->m] = this->Wmem;

    /*
     * The integer workspace iwork0 must be zeroed for a new factorization.
     * The double workspace work0 actually needs only be zeroed once in the
     * initialization of xstore. However, it is easier and more consistent
     * to do that here as well.
     */
    memset(this->iwork0, 0, this->m * sizeof(lu_int));
    memset(this->work0, 0, this->m * sizeof(lu_int));
}
