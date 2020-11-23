#ifndef _LU_INTERNAL_H
#define _LU_INTERNAL_H

#include "lu_def.h"

/* -------------------------------------------------------------------------- */
/* struct lu */
/* -------------------------------------------------------------------------- */

/*
    This data structure provides access to istore, xstore.

    lu_* routines do not access istore, xstore directly. Instead, they operate
    on a struct lu object. Scalar quantities stored in istore, xstore are copied
    to a struct lu object by lu_load() and copied back by lu_save(). Subarrays
    of istore, xstore and the user arrays Li, Lx, Ui, Ux, Wi, Wx are aliased by
    pointers in struct lu.
*/

struct lu
{
    /* user parameters, not modified */
    lu_int Lmem;
    lu_int Umem;
    lu_int Wmem;
    double droptol;
    double abstol;
    double reltol;
    lu_int nzbias;
    lu_int maxsearch;
    lu_int pad;
    double stretch;
    double compress_thres;
    double sparse_thres;
    lu_int search_rows;

    /* user readable */
    lu_int m;
    lu_int addmemL;
    lu_int addmemU;
    lu_int addmemW;

    lu_int nupdate;
    lu_int nforrest;
    lu_int nfactorize;
    lu_int nupdate_total;
    lu_int nforrest_total;
    lu_int nsymperm_total;
    lu_int Lnz;                  /* nz in L excluding diagonal */
    lu_int Unz;                  /* nz in U excluding diagonal */
    lu_int Rnz;                  /* nz in update etas excluding diagonal */
    double min_pivot;
    double max_pivot;
    double max_eta;
    double update_cost_numer;
    double update_cost_denom;
    double time_factorize;
    double time_solve;
    double time_update;
    double time_factorize_total;
    double time_solve_total;
    double time_update_total;
    lu_int Lflops;
    lu_int Uflops;
    lu_int Rflops;
    double condestL;
    double condestU;
    double normL;
    double normU;
    double normestLinv;
    double normestUinv;
    double onenorm;             /* 1-norm and inf-norm of matrix after fresh  */
    double infnorm;             /* factorization with dependent cols replaced */
    double residual_test;       /* computed by lu_residual_test() */

    lu_int matrix_nz;           /* nz in basis matrix when factorized */
    lu_int rank;                /* rank of basis matrix when factorized */
    lu_int bump_size;
    lu_int bump_nz;
    lu_int nsearch_pivot;       /* # rows/cols searched for pivot */
    lu_int nexpand;             /* # rows/cols expanded in factorize */
    lu_int ngarbage;            /* # garbage collections in factorize */
    lu_int factor_flops;        /* # flops in factorize */
    double time_singletons;
    double time_search_pivot;
    double time_elim_pivot;

    double pivot_error;         /* error estimate for pivot in last update */

    /* private */
    lu_int task;                /* the part of factorization in progress */
    lu_int pivot_row;           /* chosen pivot row */
    lu_int pivot_col;           /* chosen pivot column */
    lu_int ftran_for_update;    /* >= 0 if FTRAN prepared for update */
    lu_int btran_for_update;    /* >= 0 if BTRAN prepared for update */
    lu_int marker;              /* see @marked, below */
    lu_int pivotlen;            /* length of @pivotcol, @pivotrow; <= 2*m */
    lu_int rankdef;             /* # columns removed from active submatrix
                                   because maximum was 0 or < abstol */
    lu_int min_colnz;           /* colcount lists 1..min_colnz-1 are empty */
    lu_int min_rownz;           /* rowcount lists 1..min_rownz-1 are empty */

    /* aliases to user arrays */
    lu_int *Lindex, *Uindex, *Windex;
    double *Lvalue, *Uvalue, *Wvalue;

    /*
     * pointers into istore
     *
     * When two declaration lists are on one line, then the arrays from the
     * second list share memory with the array from the first list. The arrays
     * from the first lists are used during factorization, the arrays from the
     * second lists are used during solves/updates.
     */

    /* documented in lu_singletons.c, lu_setup_bump.c, lu_build_factors.c */
    lu_int *colcount_flink;     lu_int *pivotcol;
    lu_int *colcount_blink;     lu_int *pivotrow;
    lu_int *rowcount_flink;     lu_int *Rbegin, *eta_row;
    lu_int *rowcount_blink;     lu_int *iwork1;
    lu_int *Wbegin;             lu_int *Lbegin;    /* + Wbegin reused */
    lu_int *Wend;               lu_int *Ltbegin;   /* + Wend   reused */
    lu_int *Wflink;             lu_int *Ltbegin_p; /* + Wflink reused */
    lu_int *Wblink;             lu_int *p;         /* + Wblink reused */
    lu_int *pinv;               lu_int *pmap;
    lu_int *qinv;               lu_int *qmap;
    lu_int *Lbegin_p;           /* Lbegin_p reused */
    lu_int *Ubegin;             /* Ubegin   reused */

    lu_int *iwork0;             lu_int *marked;
    /* iwork0: size m workspace, zeroed */
    /* marked: size m workspace, 0 <= marked[i] <= @marker */

    /* pointers into xstore */
    double *work0;              /* size m workspace, zeroed */
    double *work1;              /* size m workspace, uninitialized */
    double *col_pivot;          /* pivot elements by column index */
    double *row_pivot;          /* pivot elements by row index */
};


/* -------------------------------------------------------------------------- */
/* Internal function prototypes */
/* -------------------------------------------------------------------------- */

lu_int lu_load(
    struct lu *this, lu_int *istore, double *xstore, lu_int *Li, double *Lx,
    lu_int *Ui, double *Ux, lu_int *Wi, double *Wx);

lu_int lu_save(
    const struct lu *this, lu_int *istore, double *xstore, lu_int status);

void lu_reset(struct lu *this);

void lu_initialize(lu_int m, lu_int *istore, double *xstore);

lu_int lu_factorize_bump (struct lu *this);

lu_int lu_build_factors(struct lu *this);

void lu_garbage_perm(struct lu *this);

lu_int lu_markowitz(struct lu *this);

lu_int lu_pivot(struct lu *this);

lu_int lu_setup_bump(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx);

lu_int lu_singletons(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx);

void lu_solve_dense(
    struct lu *this, const double *rhs, double *lhs, char trans);

lu_int lu_solve_for_update(
    struct lu *this, const lu_int nrhs, const lu_int *irhs, const double *xrhs,
    lu_int *nlhs, lu_int *ilhs, double *xlhs, char trans);

void lu_solve_sparse(
    struct lu *this, const lu_int nrhs, const lu_int *irhs, const double *xrhs,
    lu_int *nlhs, lu_int *ilhs, double *xlhs, char trans);

lu_int lu_dfs(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M);

lu_int lu_solve_symbolic(
    const lu_int m, const lu_int *begin, const lu_int *end, const lu_int *index,
    const lu_int nrhs, const lu_int *irhs, lu_int *ilhs, lu_int *pstack,
    lu_int *marked, const lu_int M);

lu_int lu_solve_triangular(
    const lu_int nz_symb, const lu_int *pattern_symb, const lu_int *begin,
    const lu_int *end, const lu_int *index, const double *value,
    const double *pivot, const double droptol, double *lhs, lu_int *pattern,
    lu_int *flops);

lu_int lu_update(struct lu *this, double xtbl);

double lu_condest(
    lu_int m, const lu_int *Ubegin, const lu_int *Ui, const double *Ux,
    const double *pivot, const lu_int *perm, int upper, double *work,
    double *norm, double *norminv);

double lu_normest(
    lu_int m, const lu_int *Ubegin, const lu_int *Ui, const double *Ux,
    const double *pivot, const lu_int *perm, int upper, double *work);

void lu_matrix_norm(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx);

void lu_residual_test(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx);

#endif
