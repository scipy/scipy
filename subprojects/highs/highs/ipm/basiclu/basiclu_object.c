/*
 * basiclu_object.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

/*
 * lu_free()
 *
 * p = lu_free(p) deallocates p if not NULL and sets p to NULL.
 */
static void *lu_free(void *p)
{
    if (p) free(p);
    return NULL;
}

/*
 * lu_reallocix() - reallocate two arrays
 *
 * @nz new number of elements per array
 * @p_Ai location of pointer to integer memory
 * @p_Ax location of pointer to floating point memory
 *
 * If a reallocation fails, then do not overwrite the old pointer.
 *
 * Return: BASICLU_OK or BASICLU_ERROR_out_of_memory
 */
static lu_int lu_reallocix(lu_int nz, lu_int **p_Ai, double **p_Ax)
{
    lu_int *Ainew;
    double *Axnew;

    Ainew = realloc(*p_Ai, nz * sizeof(lu_int));
    if (Ainew)
        *p_Ai = Ainew;
    Axnew = realloc(*p_Ax, nz * sizeof(double));
    if (Axnew)
        *p_Ax = Axnew;

    return Ainew && Axnew ? BASICLU_OK : BASICLU_ERROR_out_of_memory;
}

/*
 * lu_realloc_obj()
 *
 * reallocate Li,Lx and/or Ui,Ux and/or Wi,Wx as requested in xstore
 *
 * Return: BASICLU_OK or BASICLU_ERROR_out_of_memory
 */
static lu_int lu_realloc_obj(struct basiclu_object *obj)
{
    double *xstore = obj->xstore;
    lu_int addmemL = xstore[BASICLU_ADD_MEMORYL];
    lu_int addmemU = xstore[BASICLU_ADD_MEMORYU];
    lu_int addmemW = xstore[BASICLU_ADD_MEMORYW];
    double realloc_factor = fmax(1.0, obj->realloc_factor);
    lu_int nelem;
    lu_int status = BASICLU_OK;

    if (status == BASICLU_OK && addmemL > 0)
    {
        nelem = xstore[BASICLU_MEMORYL] + addmemL;
        nelem *= realloc_factor;
        status = lu_reallocix(nelem, &obj->Li, &obj->Lx);
        if (status == BASICLU_OK)
            xstore[BASICLU_MEMORYL] = nelem;
    }
    if (status == BASICLU_OK && addmemU > 0)
    {
        nelem = xstore[BASICLU_MEMORYU] + addmemU;
        nelem *= realloc_factor;
        status = lu_reallocix(nelem, &obj->Ui, &obj->Ux);
        if (status == BASICLU_OK)
            xstore[BASICLU_MEMORYU] = nelem;
    }
    if (status == BASICLU_OK && addmemW > 0)
    {
        nelem = xstore[BASICLU_MEMORYW] + addmemW;
        nelem *= realloc_factor;
        status = lu_reallocix(nelem, &obj->Wi, &obj->Wx);
        if (status == BASICLU_OK)
            xstore[BASICLU_MEMORYW] = nelem;
    }
    return status;
}

/* 
 * isvalid() - test if @obj is an allocated BASICLU object
 */
static lu_int isvalid(struct basiclu_object *obj)
{
    return obj && obj->istore && obj->xstore;
}

/*
 * lu_clear_lhs() - reset contents of lhs to zero
 */
static void lu_clear_lhs(struct basiclu_object *obj)
{
    lu_int m = obj->xstore[BASICLU_DIM];
    lu_int nzsparse = obj->xstore[BASICLU_SPARSE_THRESHOLD] * m;
    lu_int nz = obj->nzlhs;
    lu_int p;

    if (nz)
    {
        if (nz <= nzsparse)
            for (p = 0; p < nz; p++)
                obj->lhs[obj->ilhs[p]] = 0;
        else
            memset(obj->lhs, 0, m*sizeof(double));
        obj->nzlhs = 0;
    }
}

/*
 * basiclu_obj_initialize()
 */
lu_int basiclu_obj_initialize(struct basiclu_object *obj, lu_int m)
{
    lu_int imemsize, xmemsize, fmemsize;

    if (!obj)
        return BASICLU_ERROR_argument_missing;
    if (m < 0)
        return BASICLU_ERROR_invalid_argument;

    if (m == 0)
    {
        obj->istore = NULL;
        obj->xstore = NULL;
        obj->Li = NULL;
        obj->Lx = NULL;
        obj->Ui = NULL;
        obj->Ux = NULL;
        obj->Wi = NULL;
        obj->Wx = NULL;
        obj->lhs = NULL;
        obj->ilhs = NULL;
        obj->nzlhs = 0;
        return BASICLU_OK;
    }

    imemsize = BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * m;
    xmemsize = BASICLU_SIZE_XSTORE_1 + BASICLU_SIZE_XSTORE_M * m;
    fmemsize = m;               /* initial length of Li, Lx, Ui, Ux, Wi, Wx */

    obj->istore = malloc(imemsize * sizeof(lu_int));
    obj->xstore = malloc(xmemsize * sizeof(double));
    obj->Li = malloc(fmemsize * sizeof(lu_int));
    obj->Lx = malloc(fmemsize * sizeof(double));
    obj->Ui = malloc(fmemsize * sizeof(lu_int));
    obj->Ux = malloc(fmemsize * sizeof(double));
    obj->Wi = malloc(fmemsize * sizeof(lu_int));
    obj->Wx = malloc(fmemsize * sizeof(double));
    obj->lhs = calloc(m, sizeof(double));
    obj->ilhs = malloc(m * sizeof(lu_int));
    obj->nzlhs = 0;
    obj->realloc_factor = 1.5;

    if (! (obj->istore && obj->xstore && obj->Li && obj->Lx && obj->Ui &&
           obj->Ux && obj->Wi && obj->Wx && obj->lhs && obj->ilhs))
    {
        basiclu_obj_free(obj);
        return BASICLU_ERROR_out_of_memory;
    }

    lu_initialize(m, obj->istore, obj->xstore);
    obj->xstore[BASICLU_MEMORYL] = fmemsize;
    obj->xstore[BASICLU_MEMORYU] = fmemsize;
    obj->xstore[BASICLU_MEMORYW] = fmemsize;
    return BASICLU_OK;
}

/*
 * basiclu_obj_free()
 */
void basiclu_obj_free(struct basiclu_object *obj)
{
    if (obj)
    {
        obj->istore = lu_free(obj->istore);
        obj->xstore = lu_free(obj->xstore);
        obj->Li = lu_free(obj->Li);
        obj->Lx = lu_free(obj->Lx);
        obj->Ui = lu_free(obj->Ui);
        obj->Ux = lu_free(obj->Ux);
        obj->Wi = lu_free(obj->Wi);
        obj->Wx = lu_free(obj->Wx);
        obj->lhs = lu_free(obj->lhs);
        obj->ilhs = lu_free(obj->ilhs);
        obj->nzlhs = -1;
    }
}

/*
 * basiclu_obj_factorize()
 */
lu_int basiclu_obj_factorize(struct basiclu_object *obj, const lu_int *Bbegin,
                             const lu_int *Bend, const lu_int *Bi,
                             const double *Bx)
{
    lu_int status;

    if (!isvalid(obj))
        return BASICLU_ERROR_invalid_object;

    status = basiclu_factorize(obj->istore, obj->xstore, obj->Li, obj->Lx,
                               obj->Ui, obj->Ux, obj->Wi, obj->Wx, Bbegin, Bend,
                               Bi, Bx, 0);

    while (status == BASICLU_REALLOCATE)
    {
        status = lu_realloc_obj(obj);
        if (status != BASICLU_OK)
            break;
        status = basiclu_factorize(obj->istore, obj->xstore, obj->Li, obj->Lx,
                                   obj->Ui, obj->Ux, obj->Wi, obj->Wx, Bbegin,
                                   Bend, Bi, Bx, 1);
    }

    return status;
}

/*
 * basiclu_obj_get_factors()
 */
lu_int basiclu_obj_get_factors(struct basiclu_object *obj,
                               lu_int rowperm[], lu_int colperm[],
                               lu_int Lcolptr[], lu_int Lrowidx[],
                               double Lvalue[], lu_int Ucolptr[],
                               lu_int Urowidx[], double Uvalue[])
{
    if (!isvalid(obj))
        return BASICLU_ERROR_invalid_object;

    return basiclu_get_factors(obj->istore, obj->xstore, obj->Li, obj->Lx,
                               obj->Ui, obj->Ux, obj->Wi, obj->Wx, rowperm,
                               colperm, Lcolptr, Lrowidx, Lvalue, Ucolptr,
                               Urowidx, Uvalue);
}

/*
 * basiclu_obj_solve_dense()
 */
lu_int basiclu_obj_solve_dense(struct basiclu_object *obj, const double rhs[],
                               double lhs[], char trans)
{
    if (!isvalid(obj))
        return BASICLU_ERROR_invalid_object;

    return basiclu_solve_dense(obj->istore, obj->xstore, obj->Li, obj->Lx,
                               obj->Ui, obj->Ux, obj->Wi, obj->Wx, rhs, lhs,
                               trans);
}

/*
 * basiclu_obj_solve_sparse()
 */
lu_int basiclu_obj_solve_sparse(struct basiclu_object *obj, lu_int nzrhs,
                                const lu_int irhs[], const double xrhs[],
                                char trans)
{
    if (!isvalid(obj))
        return BASICLU_ERROR_invalid_object;

    lu_clear_lhs(obj);
    return basiclu_solve_sparse(obj->istore, obj->xstore, obj->Li, obj->Lx,
                                obj->Ui, obj->Ux, obj->Wi, obj->Wx, nzrhs, irhs,
                                xrhs, &obj->nzlhs, obj->ilhs, obj->lhs, trans);
}

/*
 * basiclu_obj_solve_for_update()
 */
lu_int basiclu_obj_solve_for_update(struct basiclu_object *obj, lu_int nzrhs,
                                    const lu_int irhs[], const double xrhs[],
                                    char trans, lu_int want_solution)
{
    lu_int status = BASICLU_OK;

    if (!isvalid(obj))
        return BASICLU_ERROR_invalid_object;

    lu_clear_lhs(obj);
    while (status == BASICLU_OK)
    {
        status = basiclu_solve_for_update(obj->istore, obj->xstore, obj->Li,
                                          obj->Lx, obj->Ui, obj->Ux, obj->Wi,
                                          obj->Wx, nzrhs, irhs, xrhs,
                                          want_solution ? &obj->nzlhs : NULL,
                                          obj->ilhs, obj->lhs, trans);
        if (status != BASICLU_REALLOCATE)
            break;
        status = lu_realloc_obj(obj);
    }

    return status;
}

/*
 * basiclu_obj_update()
 */
lu_int basiclu_obj_update(struct basiclu_object *obj, double xtbl)
{
    lu_int status = BASICLU_OK;

    if (!isvalid(obj))
        return BASICLU_ERROR_invalid_object;

    while (status == BASICLU_OK)
    {
        status = basiclu_update(obj->istore, obj->xstore, obj->Li, obj->Lx,
                                obj->Ui, obj->Ux, obj->Wi, obj->Wx, xtbl);
        if (status != BASICLU_REALLOCATE)
            break;
        status = lu_realloc_obj(obj);
    }

    return status;
}
