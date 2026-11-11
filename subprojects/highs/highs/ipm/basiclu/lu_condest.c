/*
 * lu_condest.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * LINPACK condition number estimate
 *
 */

#include "ipm/basiclu/lu_internal.h"

/*
 * lu_condest()
 *
 * Given m-by-m matrix U such that U[perm,perm] is upper triangular,
 * return estimate for 1-norm condition number of U.
 * If @norm is not NULL, it holds the 1-norm of the matrix on return.
 * If @norminv is not NULL, it holds the estimated 1-norm of the inverse on
 * return.
 *
 * The other function arguments are the same as in lu_normest().
 *
 */
double lu_condest(
    lu_int m, const lu_int *Ubegin, const lu_int *Ui, const double *Ux,
    const double *pivot, const lu_int *perm, int upper, double *work,
    double *norm, double *norminv)
{
    lu_int j, p;
    double Unorm, Uinvnorm;

    /* compute 1-norm of U */
    Unorm = 0;
    for (j = 0; j < m; j++)
    {
        double colsum = pivot ? fabs(pivot[j]) : 1;
        for (p = Ubegin[j]; Ui[p] >= 0; p++)
            colsum += fabs(Ux[p]);
        Unorm = fmax(Unorm, colsum);
    }

    /* estimate 1-norm of U^{-1} */
    Uinvnorm = lu_normest(m, Ubegin, Ui, Ux, pivot, perm, upper, work);

    if (norm) *norm = Unorm;
    if (norminv) *norminv = Uinvnorm;

    return Unorm * Uinvnorm;
}

/*
 * lu_normest()
 *
 * Given m-by-m matrix U such that U[perm,perm] is triangular,
 * estimate 1-norm of U^{-1} by computing
 *
 *   U'x = b, Uy = x, normest = max{norm(y)_1/norm(x)_1, norm(x)_inf},
 *
 * where the entries of b are +/-1 chosen dynamically to make x large.
 * The method is described in [1].
 *
 * @Ubegin, @Ui, @Ux matrix U in compressed column format without pivots,
 *                   columns are terminated by a negative index
 * @pivot pivot elements by column index of U; NULL if unit pivots
 * @perm permutation to triangular form; NULL if identity
 * @upper nonzero if permuted matrix is upper triangular; zero if lower
 * @work size m workspace, uninitialized on entry/return
 *
 * Return: estimate for 1-norm of U^{-1}
 *
 * [1] I. Duff, A. Erisman, J. Reid, "Direct Methods for Sparse Matrices"
 *
 */
double lu_normest(
    lu_int m, const lu_int *Ubegin, const lu_int *Ui, const double *Ux,
    const double *pivot, const lu_int *perm, int upper, double *work)
{
    lu_int i, j, k, kbeg, kend, kinc, p;
    double x1norm, xinfnorm, y1norm, temp;

    x1norm = 0;
    xinfnorm = 0;
    if (upper)
    {
        kbeg = 0; kend = m; kinc = 1;
    }
    else
    {
        kbeg = m-1; kend = -1; kinc = -1;
    }
    for (k = kbeg; k != kend; k += kinc)
    {
        j = perm ? perm[k] : k;
        temp = 0;
        for (p = Ubegin[j]; (i = Ui[p]) >= 0; p++)
            temp -= work[i] * Ux[p];
        temp += temp >= 0 ? 1 : -1; /* choose b[i] = 1 or b[i] = -1 */
        if (pivot) temp /= pivot[j];
        work[j] = temp;
        x1norm += fabs(temp);
        xinfnorm = fmax(xinfnorm, fabs(temp));
    }

    y1norm = 0;
    if (upper)
    {
        kbeg = m-1; kend = -1; kinc = -1;
    }
    else
    {
        kbeg = 0; kend = m; kinc = 1;
    }
    for (k = kbeg; k != kend; k += kinc)
    {
        j = perm ? perm[k] : k;
        if (pivot) work[j] /= pivot[j];
        temp = work[j];
        for (p = Ubegin[j]; (i = Ui[p]) >= 0; p++)
            work[i] -= temp * Ux[p];
        y1norm += fabs(temp);
    }

    return fmax(y1norm/x1norm, xinfnorm);
}
