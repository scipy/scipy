/*
 * basiclu_factorize.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

lu_int basiclu_factorize
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    const lu_int Bbegin[],
    const lu_int Bend[],
    const lu_int Bi[],
    const double Bx[],
    lu_int c0ntinue
)
{
    struct lu this;
    lu_int status;
    double factor_cost;

    status = lu_load(&this, istore, xstore, Li, Lx, Ui, Ux, Wi, Wx);
    if (status != BASICLU_OK)
        return status;

    if (! (Li && Lx && Ui && Ux && Wi && Wx && Bbegin && Bend && Bi && Bx))
    {
        status = BASICLU_ERROR_argument_missing;
        return lu_save(&this, istore, xstore, status);
    }
    if (!c0ntinue)
    {
        lu_reset(&this);
        this.task = SINGLETONS;
    }

    /* continue factorization */
    switch (this.task)
    {
    case SINGLETONS: goto singletons;
    case SETUP_BUMP: goto setup_bump;
    case FACTORIZE_BUMP: goto factorize_bump;
    case BUILD_FACTORS: goto build_factors;
    }
    status = BASICLU_ERROR_invalid_call;
    return lu_save(&this, istore, xstore, status);

    /*
     * Each of the following four parts of the factorization calls a routine
     * lu_do_something() which may request reallocation. In this case return
     * to the caller immediately, keeping the entry point in this.task.
     */

singletons:
    this.task = SINGLETONS;
    status = lu_singletons(&this, Bbegin, Bend, Bi, Bx);
    if (status != BASICLU_OK)
        goto return_to_caller;

setup_bump:
    this.task = SETUP_BUMP;
    status = lu_setup_bump(&this, Bbegin, Bend, Bi, Bx);
    if (status != BASICLU_OK)
        goto return_to_caller;

factorize_bump:
    this.task = FACTORIZE_BUMP;
    status = lu_factorize_bump(&this);
    if (status != BASICLU_OK)
        goto return_to_caller;

build_factors:
    this.task = BUILD_FACTORS;
    status = lu_build_factors(&this);
    if (status != BASICLU_OK)
        goto return_to_caller;

    /* factorization successfully finished */
    this.task               = NO_TASK;
    this.nupdate            =  0; /* make factorization valid */
    this.ftran_for_update   = -1;
    this.btran_for_update   = -1;
    this.nfactorize++;

    this.condestL = lu_condest(this.m, this.Lbegin, this.Lindex, this.Lvalue,
                               NULL, this.p, 0, this.work1,
                               &this.normL, &this.normestLinv);
    this.condestU = lu_condest(this.m, this.Ubegin, this.Uindex, this.Uvalue,
                               this.row_pivot, this.p, 1, this.work1,
                               &this.normU, &this.normestUinv);

    /* measure numerical stability of the factorization */
    lu_residual_test(&this, Bbegin, Bend, Bi, Bx);

    /*
     * factor_cost is a deterministic measure of the factorization cost.
     * The parameters have been adjusted such that (on my computer)
     * 1e-6 * factor_cost =~ time_factorize.
     *
     * update_cost measures the accumulated cost of updates/solves compared
     * to the last factorization. It is computed from
     *
     *   update_cost = update_cost_numer / update_cost_denom.
     *
     * update_cost_denom is fixed here.
     * update_cost_numer is zero here and increased by solves/updates.
     *
     */
    factor_cost =
        0.04 * this.m +
        0.07 * this.matrix_nz +
        0.20 * this.bump_nz +
        0.20 * this.nsearch_pivot +
        0.008 * this.factor_flops;

    this.update_cost_denom = factor_cost * 250;

    if (this.rank < this.m)
        status = BASICLU_WARNING_singular_matrix;

return_to_caller:
    return lu_save(&this, istore, xstore, status);
}
