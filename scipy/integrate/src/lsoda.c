#include "lsoda.h"
#include <stdio.h>

// Forward declarations (these functions are implemented later in this file)
static void ewset(const int n, const int itol, double* rtol, double* atol, double* ycur, double* ewt);
static double vmnorm(const int n, double* restrict v, double* restrict w);
static void intdy(const double t, const int k, double* yh, const int nyh, double* dky, int* iflag, lsoda_common_struct_t* S);


/**
 * @brief Computes the norm of a banded n by n matrix consistent with the weighted max-norm on vectors.
 *
 * This function computes the norm of a banded matrix stored in the array a,
 * using the weights in the array w. The matrix is stored in LAPACK banded format,
 * with nra as the leading dimension, and bandwidths ml (lower) and mu (upper).
 *
 * The norm is defined as:
 *   bnorm = max_{i=0,...,n-1} ( w[i] * sum_{j=jlo}^{jhi} |a[(i1-j-1) + j*nra]| / w[j] )
 * where:
 *   i1 = i + mu + 1
 *   jlo = max(i - ml, 0)
 *   jhi = min(i + mu, n - 1)
 *
 * @param n    Number of rows/columns of the matrix.
 * @param a    Pointer to the banded matrix data (LAPACK format).
 * @param nra  Leading dimension of the matrix a (nra >= ml + mu + 1).
 * @param ml   Lower half-bandwidth.
 * @param mu   Upper half-bandwidth.
 * @param w    Pointer to the weights array of length n.
 * @return     The computed weighted max-norm of the banded matrix.
 */
static double
bnorm(const int n, double* restrict a, const int nra, const int ml, const int mu,
      double* restrict w)
{
    int i1, jlo, jhi;
    double an = 0.0, sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum = 0.0;
        i1 = i + mu + 1;
        jlo = (i - ml < 0 ? 0 : i - ml);
        jhi = (i + mu > n - 1 ? n - 1 : i + mu);
        for (int j = jlo; j <= jhi; j++)
        {
            sum = sum + fabs(a[(i1 - j - 1) + j*nra]) / w[j];  // Original formula was correct
        }
        an = fmax(an, sum * w[i]);
    }
    return an;
}


/**
 * @brief Computes the norm of a matrix using the weighted max-norm.
 *
 * This function computes the norm of a full n by n matrix, stored in
 * the array a, that is consistent with the weighted max-norm on vectors,
 * with weights stored in the array w..
 *
 * @param n Size of the square matrix.
 * @param a Input matrix of size n by n.
 * @param w Weights vector of length n.
 * @return The computed norm of the matrix.
 *
 */
static double
fnorm(const int n, double* restrict a, double* restrict w)
{
    double an = 0.0;
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum = sum + fabs(a[i * n + j]) / w[j];
        }
        // 10
        an = fmax(an, sum*w[i]);
    }
    // 20
    return an;
}


/**
 * @brief Computes the weighted max-norm of a vector.
 *
 * This function computes the weighted max-norm of the vector of length n contained in the array v,
 * with weights contained in the array w of length n.
 *
 * The result is:
 *   vmnorm = max_{i=1,...,n} (abs(v[i]) * w[i])
 *
 * @param n Length of the vector.
 * @param v Input vector of length n.
 * @param w Weights vector of length n.
 * @return The weighted max-norm.
 */
static double
vmnorm(const int n, double* restrict v, double* restrict w)
{
    double vm = 0.0;
    for (int i = 0; i < n; i++)
    {
        vm = fmax(vm, fabs(v[i]) * w[i]);
    }
    return vm;
}


/**
 * @brief Sets method coefficients and test constants for the ODE integrator.
 *
 * This function initializes the arrays elco and tesco with coefficients
 * required by the integrator, depending on the method selected by meth.
 * cfode is called by the integrator routine to set coefficients needed there. The
 * coefficients for the current method, as given by the value of meth, are set for
 * all orders and saved.
 *
 * cfode is called once at the beginning of the problem, and is not called again
 * unless and until meth is changed.
 *
 * For meth = 1 (implicit Adams), the maximum order is 12.
 * For meth = 2 (BDF), the maximum order is 5.
 * (a smaller value of the maximum order is also allowed.).
 *
 * The elco array contains the basic method coefficients. The coefficients el[i],
 * 0 <= i <= nq, for the method of order nq are stored in elco(i,nq). They are
 * given by a generating polynomial, i.e.,
 *
 *    l(x) = el[0] + el[1]*x + ... + el(nq)*x**nq.
 *
 * for the implicit adams methods, l(x) is given by
 *
 *    dl/dx = (x+1)*(x+2)*...*(x+nq-1)/(nq-1)!, l(-1) = 0.
 *
 * for the BDF methods, l(x) is given by
 *
 *    l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
 * where
 *     k = (nq!)*(1 + 1/2 + ... + 1/nq).
 *
 * the tesco array contains test constants used for the local error test and the
 * selection of step size and/or order. At order nq, tesco[k,nq] is used for the
 * selection of step size at order nq - 1 if k = 1, at order nq if k = 2, and at
 * order nq + 1 if k = 3.
 *
 * @param meth   Method selector (1 for Adams, 2 for BDF).
 * @param elco   Output array for method coefficients, shape (13,12).
 * @param tesco  Output array for test constants, shape (3,12).
 *
 */
static void
cfode(const int meth, double* elco, double* tesco)
{
    double pc[12] = {0.0};

    if (meth == 1)
    {
        elco[0]   = 1.0;
        elco[1]   = 1.0;
        tesco[0]  = 0.0;
        tesco[1]  = 2.0;
        tesco[3]  = 1.0;
        tesco[35] = 0.0;
        pc[0] = 1.0;
        double rqfac = 1.0;
        for (int nq = 1; nq < 12; nq++)
        {
            double rq1fac = rqfac;
            rqfac = rqfac / nq;
            pc[nq] = 0.0;
            for (int ib = 0; ib <= nq - 1; ib++)
            {
                pc[nq - ib] = pc[nq - ib - 1] + nq*pc[nq - ib];
            }
            // 110
            pc[0] = nq*pc[0];
            // compute integral, -1 to 0, of p(x) and x*p(x)
            double pint = pc[0];
            double xpin = pc[0]*0.5;
            double tsign = 1.0;
            for (int i = 1; i <= nq; i++)
            {
                tsign = -tsign;
                pint = pint + tsign * pc[i] / (i + 1);
                xpin = xpin + tsign * pc[i] / (i + 2);
            }
            // 120
            // Store coefficients in elco and tesco
            elco[0 + nq*13] = pint * rq1fac;
            elco[1 + nq*13] = 1.0;

            for (int i = 1; i <= nq; i++)
            {
                elco[(i + 1) + nq*13] = rq1fac * pc[i] / (i + 1);
            }
            // 130
            double agamq = rqfac * xpin;
            double ragq = 1.0 / agamq;
            tesco[1 + 3*nq] = ragq;
            if (nq < 11)
            {
                tesco[3*(nq + 1)] = ragq * rqfac / (nq + 1);
            }
            tesco[2 + 3*(nq - 1)] = ragq;
        }
        // 140
        return;
    }

    // meth == 2
    // the pc array will contain the coefficients of the polynomial
    //     p(x) = (x+1)*(x+2)*...*(x+nq).
    // initially, p(x) = 1.

    pc[0] = 1.0;
    double rq1fac = 1.0;
    for (int nq = 0; nq < 5; nq++)
    {
        // Form coefficients of p(x)*(x+nq)
        pc[nq + 1] = 0.0;
        for (int ib = 0; ib <= nq; ib++)
        {
            pc[nq - ib + 1] = pc[nq - ib] + (nq + 1)*pc[nq - ib + 1];
        }
        // 210
        pc[0] = (nq + 1)*pc[0];

        // Store coefficients in elco and tesco
        for (int i = 0; i <= nq + 1; i++)
        {
            elco[i + nq*13] = pc[i] / pc[1];
        }
        // 220
        elco[1 + nq*13] = 1.0;
        tesco[0 + 3*nq] = rq1fac;
        tesco[1 + 3*nq] = (nq + 2) / elco[nq*13];
        tesco[2 + 3*nq] = (nq + 3) / elco[nq*13];
        rq1fac = rq1fac / (nq + 1);
    }
    // 230

    return;

}


/**
 * @brief Computes and processes the matrix P = I - h*el[0]*J, where J is an approximation to the Jacobian.
 *
 * prja is called by stoda to compute and process the matrix P.
 * - J is computed by the user-supplied routine jac if miter == 0 or 3, or by finite differencing if miter == 1 or 4.
 * - J, scaled by -h*el[0], is stored in wm.
 * - The norm of J (matrix norm consistent with the weighted max-norm on vectors given by vmnorm) is computed, and J is overwritten by P.
 * - P is then subjected to LU decomposition in preparation for later solution of linear systems with P as coefficient matrix.
 *   This is done by dgetrf if miter == 0 or 1, and by dgbtrf if miter == 3 or 4.
 *
 * In addition to variables described previously, communication with prja uses the following:
 * @param y      Array containing predicted values on entry.
 * @param ftem   Work array of length n (acor in stoda).
 * @param savf   Array containing f evaluated at predicted y.
 * @param wm     Real work space for matrices. On output, it contains the LU decomposition of P.
 *               Storage of matrix elements starts at wm[2].
 *               wm also contains the following matrix-related data:
 *               - wm[0] = sqrt(uround), used in numerical Jacobian increments.
 * @param iwm    Integer work space containing pivot information, starting at iwm[20].
 *               iwm also contains the band parameters ml = iwm[0] and mu = iwm[1] if miter == 4 or 5.
 * @param el0    el[0] (input).
 * @param pdnorm Norm of Jacobian matrix (output).
 * @param ierpj  Output error flag: 0 if no trouble, >0 if P matrix found to be singular.
 * @param jcur   Output flag: 1 to indicate that the Jacobian matrix (or approximation) is now current.
 *
 * This routine also uses the common variables el0, h, tn, uround, miter, n, nfe, and nje.
 */
static void
prja(
    int* neq, double* y, double* yh, int nyh, double* ewt, double* ftem, double* savf,
    double* wm, int* iwm, lsoda_func_t f, lsoda_jac_t jac, lsoda_common_struct_t* S)
{
    (void)nyh;
    int ier = 0, lenp, mba, mband, meband, ml, ml2, mu;
    S->nje += 1;
    S->ierpj = 0;
    S->jcur = 1;
    double hl0 = S->h * S->el0;

    // Return on unused miter value.
    if (S->miter == 2) { return; }

    // Unknown miter values via F77 Arithmetic GOTO should fall through the
    // miter == 0 case. Hence the negation in the conditions.
    // Check if miter is 0, 1 or unexpected value.
    if ((S->miter != 3) && (S->miter != 4))
    {
        // If miter is 0 or unexpected value
        if (S->miter != 1)
        {
            // If miter is 0, call jac and multiply by scalar
            int lenp = S->n * S->n;
            for (int i = 0; i < lenp; i++)
            {
                wm[i + 2] = 0.0;
            }
            // 110
            jac(neq, &S->tn, y, 0, 0, &wm[2], &lenp);
            if (*neq == -1) { return; }
            double con = -hl0;
            for (int i = 0; i < lenp; i++)
            {
                wm[i + 2] = wm[i + 2]*con;
            }
        } else {
            // miter is 1, make n calls to f to approximate j
            double fac = vmnorm(S->n, savf, ewt);
            double r0 = 1000.0 * fabs(S->h) * S->uround * S->n * fac;
            if (r0 == 0.0) { r0 = 1.0; }
            double srur  = wm[0];
            int j1 = 2;
            for (int j = 0; j < S->n; j++)
            {
                double yj = y[j];
                double r = fmax(srur * fabs(yj), r0 / ewt[j]);
                y[j] = y[j] + r;
                fac = -hl0 / r;
                f(neq, &S->tn, y, ftem);
                if (*neq == -1) { return; }
                for (int i = 0; i < S->n; i++)
                {
                    wm[i + j1] = (ftem[i] - savf[i]) * fac;
                }
                // 220
                y[j] = yj;
                j1 += S->n;
            }
            // 230
            S->nfe += S->n;
        }
        // 240
        // Compute the norm of the jacobian.
        S->pdnorm = fnorm(S->n, &wm[2], ewt) / fabs(hl0);
        // Add identity matrix
        int j = 2;
        for (int i = 0; i < S->n; i++)
        {
            wm[j] = wm[j] + 1.0;
            j += S->n + 1;
        }
        // 250

        // LU decomposition of P
        dgetrf_(&S->n, &S->n, &wm[2], &S->n, &iwm[20], &ier);
        if (ier != 0) { S->ierpj = 1; }
        return;
    }

    // Handle banded matrix cases (miter 3 and 4)
    ml = iwm[0];
    mu = iwm[1];
    mband = ml + mu + 1;
    meband = mband + ml;

    if (S->miter == 3)
    {
        // Call jac and multiply by scalar.
        ml2 = ml + 2;
        lenp = meband * S->n;
        for (int i = 0; i < lenp; i++)
        {
            wm[i + 2] = 0.0;
        }
        // 410
        jac(neq, &S->tn, y, &ml, &mu, &wm[ml2], &meband);
        if (*neq == -1) { return; }
        double con = -hl0;
        for (int i = 0; i < lenp; i++)
        {
            wm[i + 2] = wm[i + 2] * con;
        }
        // 420
    } else if (S->miter == 4) {
        // Make mband calls to f to approximate j
        mba = (mband < S->n) ? mband : S->n;
        int meb1 = meband - 1;
        double srur = wm[0];
        double fac = vmnorm(S->n, savf, ewt);
        double r0 = 1000.0 * fabs(S->h) * S->uround * S->n * fac;
        if (r0 == 0.0) { r0 = 1.0; }

        for (int j = 0; j < mba; j++)
        {
            // Perturb y values for this column group
            for (int i = j; i < S->n; i += mband)
            {
                double yi = y[i];
                double r = fmax(srur * fabs(yi), r0 / ewt[i]);
                y[i] = y[i] + r;
            }
            // 530

            // Call f once for this column group
            f(neq, &S->tn, y, ftem);
            if (*neq == -1) { return; }

            // Compute finite differences for all affected columns
            for (int jj = j; jj < S->n; jj += mband)
            {
                y[jj] = yh[jj]; // yh(jj,1) -> first column of yh
                double yjj = y[jj];
                double r = fmax(srur * fabs(yjj), r0 / ewt[jj]);
                fac = -hl0 / r;
                int i1 = (jj - mu < 0) ? 0 : jj - mu;
                int i2 = (jj + ml > S->n - 1) ? S->n - 1 : jj + ml;
                int ii = jj * meb1 - ml + 2;
                for (int i = i1; i <= i2; i++)
                {
                    wm[ii + i] = (ftem[i] - savf[i]) * fac;
                }
                // 540
            }
            // 550
        }
        // 560
        S->nfe += mba;
    }

    S->pdnorm = bnorm(S->n, &wm[2], meband, ml, mu, ewt) / fabs(hl0);
    // Add identity matrix
    int ii = 2;
    for (int i = 0; i < S->n; i++)
    {
        wm[ii] = wm[ii] + 1.0;
        ii += meband;
    }
    // 580
    // LU decomposition of P
    dgbtrf_(&S->n, &S->n, &ml, &mu, &wm[2], &meband, &iwm[20], &ier);
    if (ier != 0) { S->ierpj = 1; }
    return;

}


/**
 * @brief Manages the solution of the linear system arising from a chord iteration.
 *
 * This routine is called if miter != 0.
 * - If miter is 0 or 1, it calls dgetrs to solve the system.
 * - If miter == 2, it updates the coefficient h*el0 in the diagonal matrix, and then computes the solution.
 * - If miter is 3 or 4, it calls dgbtrs.
 *
 * Communication with solsy uses the following variables:
 * @param wm   Real work space containing the inverse diagonal matrix if miter == 3,
 *             and the LU decomposition of the matrix otherwise.
 *             Storage of matrix elements starts at wm[2].
 *             wm also contains the following matrix-related data:
 *             - wm[0] = sqrt(uround) (not used here)
 *             - wm[1] = hl0, the previous value of h*el0, used if miter == 3
 * @param iwm  Integer work space containing pivot information, starting at iwm[20],
 *             if miter is 1, 2, 4, or 5. iwm also contains band parameters
 *             ml = iwm[0] and mu = iwm[1] if miter is 4 or 5.
 * @param x    The right-hand side vector on input, and the solution vector on output, of length n.
 * @param tem  Vector of work space of length n, not used in this version.
 * @param iersl Output flag (in common). iersl = 0 if no trouble occurred,
 *              iersl = 1 if a singular matrix arose with miter == 3.
 *
 * This routine also uses the common variables el0, h, miter, and n.
 */
static void
solsy(double* wm, int* iwm, double* x, double* tem, lsoda_common_struct_t* S)
{
    (void)tem;  // tem is not used in this version
    int int1 = 1, ierr = 0;
    S->iersl = 0;

    if ((S->miter == 0) || (S->miter == 1))
    {
        dgetrs_("N", &S->n, &int1, &wm[2], &S->n, &iwm[20], x, &S->n, &ierr);
        return;
    } else if (S->miter == 2) {
        double phl0 = wm[1];
        double hl0 = S->h * S->el0;
        wm[1] = hl0;
        if (hl0 != phl0)
        {
            double r = hl0 / phl0;
            for (int i = 0; i < S->n; i++)
            {
                double di = 1.0 - r * (1.0 - 1.0 / wm[i + 2]);
                if (fabs(di) == 0.0)
                {
                    S->iersl = 1;
                    return;
                }
                wm[i + 2] = 1.0 / di;
            }
            // 320
        }
        for (int i = 0; i < S->n; i++)
        {
            x[i] = wm[i + 2] * x[i];
        }
        // 340
    } else if ((S->miter == 3) || (S->miter == 4)) {
        int ml = iwm[0];
        int mu = iwm[1];
        int meband = 2*ml + mu + 1;
        dgbtrs_("N", &S->n, &ml, &mu, &int1, &wm[2], &meband, &iwm[20], x, &S->n, &ierr);
    }

    return;
}


/**
 *
 * Fortran77 code stoda.f from original ODEPACK is exceedingly entangled with goto
 * flow control. Here we take out the corrector loop and control the flow with a
 * state machine.
 *
 */

/**
 *********************************************************************************
 *                              stoda.f implementation
 *********************************************************************************
 */

/**
 * @brief First call initialization for stoda, jstart = 0 case, codeblock up to label 100
 */
static void
stoda_first_call_init(lsoda_common_struct_t* S)
{
    S->lmax = S->maxord + 1;
    S->nq = 1;
    S->l = 2;
    S->ialth = 2;
    S->rmax = 10000.0;
    S->rc = 0.0;
    S->el0 = 1.0;
    S->crate = 0.7;
    S->hold = S->h;
    S->nslp = 0;
    S->ipup = S->miter;

    // Initialize switching parameters - meth = 1 is assumed initially
    S->icount = 20;
    S->irflag = 0;
    S->pdest = 0.0;
    S->pdlast = 0.0;
    S->ratio = 5.0;

    // Setup coefficients for both methods
    cfode(2, S->elco, S->tesco);
    for (int i = 0; i < 5; i++)
    {
        S->cm2[i] = S->tesco[1 + i*3] * S->elco[i + 1 + i*13];  // tesco(2,i+1)*elco(i+2,i+1) in FORTRAN
    }
    cfode(1, S->elco, S->tesco);
    for (int i = 0; i < 12; i++)
    {
        S->cm1[i] = S->tesco[1 + i*3] * S->elco[i + 1 + i*13];
    }
}


/**
 * @brief Restart initialization for stoda, label 100 - 150
 *
 * This block handles preliminaries needed when jstart = -1. ipup is set to miter to force a matrix update.
 * if an order increase is about to be considered (ialth = 1), ialth is reset to 2 to postpone consideration
 * one more step. If the caller has changed meth, cfode is called to reset the coefficients of the method.
 * If h is to be changed, yh must be rescaled. If h or meth is being changed, ialth is reset to l = nq + 1
 * to prevent further changes in h for that many steps.
 *
 * @param S Pointer to the common structure containing integrator state.
 */
static void
stoda_restart_init(lsoda_common_struct_t* S)
{
    // 100
    S->ipup = S->miter;
    S->lmax = S->maxord + 1;
    if (S->ialth == 1) { S->ialth = 2; }

    if (S->meth != S->mused)
    {
        cfode(S->meth, S->elco, S->tesco);
        S->ialth = S->l;
    }
}



/**
 * @brief Resets the state of the ODE integrator for a new step, label 150 - 160
 *
 * The el vector and related constants are reset whenever the order nq
 * is changed, or at the start of the problem.
 *
 * @param S Pointer to the common structure containing integrator state.
 */
static void
stoda_reset(lsoda_common_struct_t* S)
{
        // 150
        for (int i = 0; i < S->l; i++) { S->el[i] = S->elco[i + (S->nq - 1)*13]; }
        // 155
        S->nqnyh = S->nq * S->nyh;
        S->rc = S->rc* S->el[0] / S->el0;
        S->el0 = S->el[0];
        S->conit = 0.5 / (S->nq + 2);
}


/**
 * @brief Step size adjustment (labels 160-178)
 */
static void
stoda_adjust_step_size(lsoda_common_struct_t* S, double* yh, int nyh, const double* sm1)
{
    // If h is being changed, the h ratio rh is checked against
    // rmax, hmin, and hmxi, and the yh array rescaled. ialth is set to
    // l = nq + 1 to prevent a change of h for that many steps, unless
    // forced by a convergence or error test failure.
    if (S->h == S->hold) { return; }

    // 160-175: Step size ratio calculation and limiting
    double rh = S->h / S->hold;
    S->h = S->hold;

    // 170-175: Apply various constraints
    rh = fmax(rh, S->hmin / fabs(S->h));
    rh = fmin(rh, S->rmax);
    rh = rh / fmax(1.0, fabs(S->h) * S->hmxi * rh);

    // if meth = 1, also restrict the new step size by the stability region.
    // if this reduces h, set irflag to 1 so that if there are roundoff
    // problems later, we can assume that is the cause of the trouble.
    // 175-178: Adams stability region constraint
    if (S->meth != 2)
    {
        S->irflag = 0;
        double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
        if (rh * pdh * 1.00001 < sm1[S->nq - 1])
        {
            rh = sm1[S->nq - 1] / pdh;
            S->irflag = 1;
        }
    }

    // 178-180: Rescale yh array
    double r = 1.0;
    for (int j = 1; j < S->l; j++)
    {
        r = r * rh;
        for (int i = 0; i < S->n; i++)
        {
            yh[i + j*nyh] = yh[i + j*nyh] * r;
        }
    }
    S->h = S->h * rh;
    S->rc = S->rc * rh;
    S->ialth = S->l;
}


/**
 * @brief Enumeration of the states for the corrector loop to manage the corrector iterations.
 *
 */
typedef enum lsoda_corrector_status {
    CORRECTOR_CONVERGED,         // Corrector loop has converged successfully.
    CORRECTOR_NO_CONVERGENCE,    // Convergence not achieved within the maximum number of iterations.
    CORRECTOR_ERROR,             // Structural errors such as singularities or f, jac returned error.
    CORRECTOR_RETRY              // Convergence not achieved, but retry is possible.
} lsoda_corrector_status_t;


/**
 * @brief Main corrector loop factored out from stoda code. Implements (roughly) labels 200-410
 *
 * Since this function has too many exit points to labels 430, 410, 450, an enumeration above
 * is used to mark the exit destinations.
 *
 */
static lsoda_corrector_status_t
stoda_corrector_loop(
    int* neq, double* y, double* yh, int nyh, double* ewt,
    double* savf, double* acor, double* wm, int* iwm, lsoda_func_t f,
    lsoda_jac_t jac, double pnorm, int* m_out, double* del_out, lsoda_common_struct_t* S)
{
    // Up to maxcor corrector iterations are taken. A convergence test is
    // made on the r.m.s. norm of each correction, weighted by the error
    // weight vector ewt. The sum of the corrections is accumulated in the
    // vector acor(i). The yh array is not altered in the corrector loop.
    int m = 0;
    double rate = 0.0;
    double del = 0.0;
    double delp = 0.0;

    for (int i = 0; i < S->n; i++) { y[i] = yh[i]; }
    f(neq, &S->tn, y, savf);
    if (*neq == -1) { return CORRECTOR_ERROR; }
    S->nfe += 1;

    if (S->ipup > 0) {
        // If indicated, the matrix p = i - h*el(1)*j is reevaluated and
        // preprocessed before starting the corrector iteration. ipup is set
        // to 0 as an indicator that this has been done.
        prja(neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac, S);
        if (*neq == -1) { return CORRECTOR_ERROR; }
        S->ipup = 0;
        S->rc = 1.0;
        S->nslp = S->nst;
        S->crate = 0.7;
        if (S->ierpj != 0) { return CORRECTOR_NO_CONVERGENCE; }
    }

    for (int i = 0; i < S->n; i++) { acor[i] = 0.0; }

    // Main corrector iteration loop
    while (1) {
        if (S->miter == 0)
        {
            // In the case of functional iteration, update y directly from
            // the result of the last function evaluation.
            for (int i = 0; i < S->n; i++)
            {
                savf[i] = S->h * savf[i] - yh[i + S->nyh];
                y[i] = savf[i] - acor[i];
            }
            // 290
            del = vmnorm(S->n, y, ewt);
            for (int i = 0; i < S->n; i++)
            {
                y[i] = yh[i] + S->el[0] * savf[i];
                acor[i] = savf[i];
            }
            // 300
        } else {
            // in the case of the chord method, compute the corrector error,
            // and solve the linear system with that as right-hand side and
            // p as coefficient matrix.
            for (int i = 0; i < S->n; i++)
            {
                y[i] = S->h * savf[i] - (yh[i + S->nyh] + acor[i]);
            }
            // 360
            solsy(wm, iwm, y, savf, S);
            if (S->iersl < 0) { return CORRECTOR_NO_CONVERGENCE; }
            if (S->iersl > 0) {
                // goto 410 - convergence failure decision point
                if ((S->miter != 0) && (S->jcur != 1))
                {
                    // Can retry with fresh Jacobian (go to 220)
                    S->icf = 1;
                    S->ipup = S->miter;
                    return CORRECTOR_RETRY;
                }
                // Cannot retry - go to 430 (failure processing)
                return CORRECTOR_NO_CONVERGENCE;
            }
            del = vmnorm(S->n, y, ewt);
            for (int i = 0; i < S->n; i++)
            {
                acor[i] = acor[i] + y[i];
                y[i] = yh[i] + S->el[0] * acor[i];
            }
            // 380
        }

        // Test for convergence. If m > 0, an estimate of the convergence
        // rate constant is stored in crate, and this is used in the test.
        //
        // We first check for a change of iterates that is the size of
        // roundoff error. If this occurs, the iteration has converged, and a
        // new rate estimate is not formed.
        // In all other cases, force at least two iterations to estimate a
        // local lipschitz constant estimate for adams methods.
        // On convergence, form pdest = local maximum lipschitz constant
        // estimate.  pdlast is the most recent nonzero estimate.

        // 400
        if (del <= 100.0 * pnorm * S->uround) { break; } // goto 450

        if ((m != 0) || (S->meth != 1))
        {
            if (m != 0)
            {
                double rm = 1024.0;
                if (del <= 1024.0 * delp) { rm = del / delp; }
                rate = fmax(rate, rm);
                S->crate = fmax(0.2 * S->crate, rm);
            }
            // 402
            double dcon = del * fmin(1.0, 1.5 * S->crate) / (S->tesco[1 + (S->nq - 1)*3] * S->conit);  // tesco(2,nq) in FORTRAN
            if (dcon <= 1.0)
            {
                S->pdest = fmax(S->pdest, rate / fabs(S->h * S->el[0]));
                if (S->pdest != 0.0) S->pdlast = S->pdest;
                break;
            }
        }
        // 405
        m = m + 1;
        if (m == S->maxcor) {
            // Too many corrector iterations
            if (S->miter != 0 && S->jcur != 1)
            {
                S->icf = 1;
                S->ipup = S->miter;
                return CORRECTOR_RETRY;
            }

            return CORRECTOR_NO_CONVERGENCE;
        }
        if (m >= 2 && del > 2.0 * delp) {
            // Diverging
            if ((S->miter != 0) && (S->jcur != 1)) {
                S->icf = 1;
                S->ipup = S->miter;
                return CORRECTOR_RETRY;
            }
            return CORRECTOR_NO_CONVERGENCE;
        }
        delp = del;
        f(neq, &S->tn, y, savf);
        if (*neq == -1) { return CORRECTOR_ERROR; }
        S->nfe += 1;
        continue;

    }

    // The corrector has converged. jcur is set to 0  to signal that the
    // jacobian involved may need updating later. The local error test is
    // made and control passes to statement 500 if it fails.

    // 450 but without the checks. Those are handled at the exit by the caller.
    S->jcur = 0;
    *m_out = m;
    *del_out = del;
    return CORRECTOR_CONVERGED;
}


/**
 * Handle corrector convergence failure
 * Returns 0 for successful recovery (retry), -1 for fatal error
 */
static int
stoda_handle_corrector_failure(double* yh, double* yh1, int nyh, int* ncf, double told, const double* sm1, lsoda_common_struct_t* S)
{
    S->icf = 2;
    (*ncf)++;
    S->rmax = 2.0;
    S->tn = told;  // Restore to value before step attempt

    // Retract yh array to values before prediction (430-445)
    int i1 = S->nqnyh;  // Start at end of column nq (0-based indexing)
    for (int jb = 1; jb <= S->nq; jb++) {
        i1 = i1 - nyh;  // Move to start of previous column
        for (int i = i1; i < S->nqnyh; i++) {  // 0-based: from i1 to nqnyh-1
            yh1[i] = yh1[i] - yh1[i + nyh];
        }
    }

    if (S->ierpj < 0 || S->iersl < 0)
    {
        // goto 680
        S->kflag = -3;
        return -1;
    }
    if (fabs(S->h) <= S->hmin * 1.00001)
    {
        // goto 670
        S->kflag = -2;  // Step size below minimum
        return -1;
    }
    if (*ncf == S->mxncf)
    {
        // goto 670
        S->kflag = -2;  // Too many convergence failures
        return -1;
    }

    // Reduce step size and prepare for retry
    double rh = 0.25;
    S->ipup = S->miter;  // Force Jacobian update on retry

    // Apply step size reduction - use existing step adjustment logic
    rh = fmax(rh, S->hmin / fabs(S->h));
    rh = fmin(rh, S->rmax);
    rh = rh / fmax(1.0, fabs(S->h) * S->hmxi * rh);

    // Goto 170
    // For stability region constraint (meth=1 only)
    if (S->meth == 1) {
        S->irflag = 0;
        double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
        if (rh * pdh * 1.00001 < sm1[S->nq - 1]) {
            rh = sm1[S->nq - 1] / pdh;
            S->irflag = 1;
        }
    }

    // Rescale yh array with new step size
    double r = 1.0;
    for (int j = 2; j <= S->l; j++) {
        r = r * rh;
        for (int i = 0; i < S->n; i++) {
            yh[i + (j-1)*nyh] = yh[i + (j-1)*nyh] * r;  // yh(i,j) in FORTRAN
        }
    }
    S->h = S->h * rh;
    S->rc = S->rc * rh;
    S->ialth = S->l;

    return 0;  // Successful recovery - caller should retry
}


static void
stoda(
    int* neq, double* y, double* yh, int nyh, double* yh1, double* ewt,
    double* svaf, double* acor, double* wm, int* iwm, lsoda_func_t f,
    lsoda_jac_t jac, lsoda_common_struct_t* S)
{
    const double sm1[12] = {0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.20, 0.15, 0.10, 0.075, 0.050, 0.025};
    lsoda_corrector_status_t corrector_status;
    int should_reset_rmax = 0;  // Flag to track if we should go to label 690 (rmax reset + exit)

    S->kflag = 0;
    double told = S->tn;
    int ncf = 0;
    S->ierpj = 0;
    S->iersl = 0;
    S->jcur = 0;
    S->icf = 0;

    // Step 1: Initialize based on jstart
    switch (S->jstart) {
        case 0:   // First call
            stoda_first_call_init(S);
            stoda_reset(S);
            break; // Go directly to prediction

        case -1:  // Restart with possible method/step change
            stoda_restart_init(S);
            stoda_reset(S);
            stoda_adjust_step_size(S, yh, nyh, sm1);
            break;

        case -2:  // Step change only
            stoda_adjust_step_size(S, yh, nyh, sm1);
            break;

        default:  // jstart > 0, normal continuation
            break; // Go directly to prediction
    }

    // Main integration loop - implements state machine to eliminate gotos
    while (1) {
        // Step 2: Prediction phase (200-220)
        // Force Jacobian update if needed
        if (fabs(S->rc - 1.0) > S->ccmax) { S->ipup = S->miter; }
        if (S->nst >= S->nslp + S->msbp) { S->ipup = S->miter; }

        S->tn = S->tn + S->h;

        // Pascal triangle multiplication on yh (200-215)
        int i1 = S->nqnyh;  // Start at end of column nq (0-based indexing)
        for (int jb = 1; jb <= S->nq; jb++)
        {
            i1 = i1 - nyh;  // Move to start of previous column
            for (int i = i1; i < S->nqnyh; i++)  // 0-based: from i1 to nqnyh-1
            {
                yh1[i] = yh1[i] + yh1[i + nyh];
            }
        }
        double pnorm = vmnorm(S->n, yh1, ewt);

        // Corrector loop (220-430)
        int m_corrector = 0;
        double del_corrector = 0.0;

        // Corrector loop 220
        corrector_status = stoda_corrector_loop(neq, y, yh, nyh, ewt, svaf, acor, wm, iwm, f, jac, pnorm, &m_corrector, &del_corrector, S);
        // 430

        // Step 4: Handle corrector results
        switch (corrector_status) {
            case CORRECTOR_CONVERGED:
                // Continue with local error testing, exit while loop
                break;

            case CORRECTOR_RETRY:
                // Retry with fresh Jacobian - re-spin loop
                continue;

            case CORRECTOR_NO_CONVERGENCE:

                if (stoda_handle_corrector_failure(yh, yh1, nyh, &ncf, told, sm1, S) != 0)
                {
                    // 670/680 -> 720
                    S->hold = S->h;
                    S->jstart = 1;
                    return;
                }
                // A modified step size value found - respin loop
                continue;

            case CORRECTOR_ERROR:
                // 680 -> 720
                S->kflag = -3;
                S->hold = S->h;
                S->jstart = 1;
                return;
        }

        // Local error test and step/order selection (labels 450-488)
        // The corrector has converged. jcur is set to 0 to signal that the
        // jacobian involved may need updating later. The local error test is
        // made and control passes to statement 500 if it fails.
        S->jcur = 0;

        // Compute local error estimate (label 450)
        double dsm;
        if (m_corrector == 0)
        {
            dsm = del_corrector / S->tesco[1 + (S->nq - 1)*3];
        } else {
            dsm = vmnorm(S->n, acor, ewt) / S->tesco[1 + (S->nq - 1)*3];
        }

        // Local error test
        if (dsm > 1.0) {
            // The error test failed. kflag keeps track of multiple failures.
            // restore tn and the yh array to their previous values, and prepare
            // to try the step again. Compute the optimum step size for this or
            // one lower order. After 2 or more failures, h is forced to decrease
            // by a factor of 0.2 or less.

            // 500
            S->kflag -= 1;
            S->tn = told;

            int i1 = S->nqnyh;
            for (int jb = 1; jb <= S->nq; jb++) {
                i1 = i1 - nyh;
                for (int i = i1; i < S->nqnyh; i++) { yh1[i] = yh1[i] - yh1[i + nyh]; }  // 510
            }
            // 515

            S->rmax = 2.0;
            if (fabs(S->h) <= S->hmin * 1.00001)
            {
                // 660
                S->kflag = -1;
                S->hold = S->h;
                S->jstart = 1;
                return;
            }
            if (S->kflag <= -3)
            {
                // 640
                if (S->kflag == -10)
                {
                    // 660
                    S->kflag = -1;
                    S->hold = S->h;
                    S->jstart = 1;
                    return;
                }

                double rh = 0.1;
                rh = fmax(S->hmin / fabs(S->h), rh);
                S->h = S->h * rh;
                for (int i = 0; i < S->n; i++) { y[i] = yh[i]; } // 645

                f(neq, &S->tn, y, svaf);
                if (*neq == -1) { return; }

                S->nfe++;
                for (int i = 0; i < S->n; i++) { yh[i + nyh] = S->h * svaf[i]; } // 650
                S->ipup = S->miter;
                S->ialth = 5;
                if (S->nq == 1) { continue; }
                S->nq = 1;
                S->l = 2;
                stoda_reset(S);
                continue;
            }

            double rhup = 0.0;
            // 540
            double rhsm = 1.0 / (1.2 * pow(dsm, 1.0 / S->l) + 0.0000012);
            double rhdn = 0.0;

            if (S->nq != 1) {
                double ddn = vmnorm(S->n, &yh[(S->l - 1) * S->nyh], ewt) / S->tesco[(S->nq - 1) * 3];
                double exdn = 1.0 / S->nq;
                rhdn = 1.0 / (1.3 * pow(ddn, exdn) + 0.0000013);
            }

            // 550
            // If meth = 1, limit rh according to the stability region also.
            if (S->meth == 1) {
                double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
                if (S->l < S->lmax) { rhup = fmin(rhup, sm1[S->l - 1] / pdh); }
                rhsm = fmin(rhsm, sm1[S->nq - 1] / pdh);
                if (S->nq > 1) { rhdn = fmin(rhdn, sm1[S->nq - 2] / pdh); }
                S->pdest = 0.0;
            }

            // 560
            int newq;
            double rh;
            if (rhsm >= rhup) {
                // 570
                if (rhsm >= rhdn) {
                    newq = S->nq;
                    rh = rhsm;
                } else {
                    // 580
                    newq = S->nq - 1;
                    rh = rhdn;
                    if (S->kflag < 0 && rh > 1.0) {
                        rh = 1.0;
                    }
                }
            } else {
                if (rhup > rhdn) {
                    // 590
                    newq = S->l;
                    rh = rhup;
                    if (rh < 1.1) {
                        // 610 -> go to 700 (no rmax reset)
                        S->ialth = 3;
                        break;
                    }
                    // Handle acor scaling for order increase (590-600)
                    double r = S->el[S->l - 1] / (double)S->l;
                    for (int i = 0; i < S->n; i++) {
                        yh[i + newq * nyh] = acor[i] * r;
                    }
                    // 630 - update nq, l and reset coefficients
                    S->nq = newq;
                    S->l = S->nq + 1;
                    stoda_reset(S);
                    continue;
                } else {
                    // 580
                    newq = S->nq - 1;
                    rh = rhdn;
                    if (S->kflag < 0 && rh > 1.0) {
                        rh = 1.0;
                    }
                }
            }

            // 620
            if (S->meth == 1) {
                double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
                if (rh * pdh * 1.00001 < sm1[newq - 1]) {
                    // Step restricted by stability - skip 10% test (go to 625)
                } else if (S->kflag == 0 && rh < 1.1) {
                    // 622: go to 610 -> go to 700 (no rmax reset)
                    S->ialth = 3;
                    break;
                }
            } else {
                // 622: meth == 2
                if (S->kflag == 0 && rh < 1.1) {
                    // go to 610 -> go to 700 (no rmax reset)
                    S->ialth = 3;
                    break;
                }
            }

            // 625
            if (S->kflag <= -2) { rh = fmin(rh, 0.2); }

            // Apply changes and retry
            if (newq != S->nq) {
                // 630: go to 150 and then 170
                S->nq = newq;
                S->l = S->nq + 1;
                stoda_reset(S);
            }
            S->h = S->h * rh;
            S->rc = S->rc * rh;
            S->ialth = S->l;
            stoda_adjust_step_size(S, yh, nyh, sm1);
            // Check if we should go to label 690 (equivalent to FORTRAN: if (iredo .eq. 0) go to 690)
            if (should_reset_rmax) {
                S->rmax = 10.0;  // Label 690: reset rmax
                break;  // goto 700 (exit function)
            }
            continue;
        }

        // Back to happy branch after 450
        // After a successful step, update the yh array.
        // Decrease icount by 1, and if it is -1, consider switching methods.
        // If a method switch is made, reset various parameters,
        // rescale the yh array, and exit. If there is no switch,
        // consider changing h if ialth = 1. Otherwise decrease ialth by 1.
        // If ialth is then 1 and nq < maxord, then acor is saved for
        // use in a possible order increase on the next step.
        // if a change in h is considered, an increase or decrease in order
        // by one is considered also. A change in h is made only if it is by a
        // factor of at least 1.1. If not, ialth is set to 3 to prevent
        // testing for that many steps.
        S->kflag = 0;
        should_reset_rmax = 1;  // Equivalent to FORTRAN iredo = 0 at line 393
        S->nst += 1;
        S->hu = S->h;
        S->nqu = S->nq;
        S->mused = S->meth;

        for (int j = 1; j <= S->l; j++)
        {
            for (int i = 0; i < S->n; i++)
            {
                yh[i + (j-1)*nyh] = yh[i + (j-1)*nyh] + S->el[j-1] * acor[i];
            }
        }
        // 460

        // Decrease icount and consider method switching
        S->icount -= 1;
        if (S->icount >= 0)
        {
            ; // No switch
        } else {
            // Consider method switching based on current method
            if (S->meth == 2) {
                // We are currently using a bdf method. Consider switching to adams.
                // Compute the step size we could have (ideally) used on this step,
                // With the current (bdf) method, and also that for the adams.
                // if nq > mxordn, we consider changing to order mxordn on switching.
                // compare the two step sizes to decide whether to switch.
                // The step size advantage must be at least 5/ratio = 1 to switch.
                // If the step size for adams would be so small as to cause roundoff
                // pollution, we stay with bdf.
                // 480
                double exsm = 1.0 / (double)S->l;
                int nqm1;
                double rh1, exm1;

                if (S->mxordn >= S->nq) {
                    double dm1 = dsm * (S->cm2[S->nq - 1] / S->cm1[S->nq - 1]);  // C 0-based
                    rh1 = 1.0 / (1.2 * pow(dm1, exsm) + 0.0000012);
                    nqm1 = S->nq;
                    exm1 = exsm;
                } else {
                    nqm1 = S->mxordn;
                    int lm1 = S->mxordn + 1;
                    exm1 = 1.0 / (double)lm1;
                    int lm1p1 = lm1 + 1;
                    double dm1 = vmnorm(S->n, &yh[(lm1p1 - 1) * S->nyh], ewt) / S->cm1[S->mxordn - 1];  // yh(1,lm1p1) in FORTRAN, C 0-based
                    rh1 = 1.0 / (1.2 * pow(dm1, exm1) + 0.0000012);
                }

                double rh1it = 2.0 * rh1;
                double pdh = S->pdnorm * fabs(S->h);
                if (pdh * rh1 > 0.00001) {
                    rh1it = sm1[nqm1 - 1] / pdh;  // C 0-based indexing
                }
                rh1 = fmin(rh1, rh1it);
                double rh2 = 1.0 / (1.2 * pow(dsm, exsm) + 0.0000012);

                if (rh1 * S->ratio >= 5.0 * rh2) {
                    double alpha = fmax(0.001, rh1);
                    double dm1_test = pow(alpha, exm1) * dsm * (S->cm2[S->nq - 1] / S->cm1[S->nq - 1]);
                    if (dm1_test > 1000.0 * S->uround * pnorm) {
                        // Switch to Adams
                        S->h = S->h * rh1;
                        S->icount = 20;
                        S->meth = 1;
                        S->miter = 0;
                        S->pdlast = 0.0;
                        S->nq = nqm1;
                        S->l = S->nq + 1;
                        // 150 -> 170
                        stoda_reset(S);
                        should_reset_rmax = 0;  // Method switch resets iredo context
                        stoda_adjust_step_size(S, yh, nyh, sm1);
                        continue;  // Restart integration with new method
                    }
                }
                // No switch

            } else {
                // 480
                // We are currently using an adams method, consider switching to bdf.
                // If the current order is greater than 5, assume the problem is
                // not stiff, and skip this section.
                // If the lipschitz constant and error estimate are not polluted
                // by roundoff, go to 470 and perform the usual test.
                // otherwise, switch to the bdf methods if the last step was
                // restricted to insure stability (irflag = 1), and stay with adams
                // method if not.  when switching to bdf with polluted error estimates,
                // in the absence of other information, double the step size.
                //
                // when the estimates are ok, we make the usual test by computing
                // the step size we could have (ideally) used on this step,
                // with the current (adams) method, and also that for the bdf.
                // if nq .gt. mxords, we consider changing to order mxords on switching.
                // compare the two step sizes to decide whether to switch.
                // the step size advantage must be at least ratio = 5 to switch.
                // Currently using BDF - consider switching to Adams (label 480)
                if (S->nq <= 5) {
                    // Check if estimates are polluted by roundoff
                    if (dsm > 100.0 * pnorm * S->uround && S->pdest != 0.0)
                    {
                        // Estimates are clean - do full switching test (label 470)
                        double exsm = 1.0 / (double)S->l;
                        double rh1 = 1.0 / (1.2 * pow(dsm, exsm) + 0.0000012);
                        double rh1it = 2.0 * rh1;
                        double pdh = S->pdlast * fabs(S->h);
                        if (pdh * rh1 > 0.00001) { rh1it = sm1[S->nq - 1] / pdh; }
                        rh1 = fmin(rh1, rh1it);

                        int nqm2;
                        double rh2;
                        if (S->nq <= S->mxords)
                        {
                            double dm2 = dsm * (S->cm1[S->nq - 1] / S->cm2[S->nq - 1]);
                            rh2 = 1.0 / (1.2 * pow(dm2, exsm) + 0.0000012);
                            nqm2 = S->nq;
                        } else {
                            nqm2 = S->mxords;
                            int lm2 = S->mxords + 1;
                            double exm2 = 1.0 / (double)lm2;
                            int lm2p1 = lm2 + 1;
                            double dm2 = vmnorm(S->n, &yh[(lm2p1 - 1) * S->nyh], ewt) / S->cm2[S->mxords - 1];
                            rh2 = 1.0 / (1.2 * pow(dm2, exm2) + 0.0000012);
                        }

                        if (rh2 >= S->ratio * rh1)
                        {
                            // 478 switch to BDF
                            S->h = S->h * rh2;
                            S->icount = 20;
                            S->meth = 2;
                            S->miter = S->jtyp;
                            S->pdlast = 0.0;
                            S->nq = nqm2;
                            S->l = S->nq + 1;
                            // 150 -> 170
                            stoda_reset(S);
                            should_reset_rmax = 0;  // Method switch resets iredo context
                            stoda_adjust_step_size(S, yh, nyh, sm1);
                            continue;
                        }
                    } else {
                        // 470
                        if (S->irflag != 0) {
                            // Switch to BDF with doubled step size
                            S->h = S->h * 2.0;
                            int nqm2 = (S->nq < S->mxords) ? S->nq : S->mxords;
                            S->icount = 20;
                            S->meth = 2;
                            S->miter = S->jtyp;
                            S->pdlast = 0.0;
                            S->nq = nqm2;
                            S->l = S->nq + 1;
                            // 150 -> 170
                            stoda_reset(S);
                            should_reset_rmax = 0;  // Method switch resets iredo context
                            stoda_adjust_step_size(S, yh, nyh, sm1);
                            continue;
                        }
                    }
                }
            }
        }

        // 488 Start of goto olympics.
        S->ialth -= 1;
        if (S->ialth == 0) {
            // 520

            // Regardless of the success or failure of the step, factors
            // rhdn, rhsm, and rhup are computed, by which h could be multiplied
            // at order nq - 1, order nq, or order nq + 1, respectively.
            // In the case of failure, rhup = 0.0 to avoid an order increase.
            // The largest of these is determined and the new order chosen
            // accordingly. If the order is to be increased, we compute one
            // additional scaled derivative.
            double rhup = 0.0;
            if (S->l != S->lmax) {
                for (int i = 0; i < S->n; i++)
                {
                    svaf[i] = acor[i] - yh[i + (S->lmax - 1) * nyh];
                }
                // 530
                double dup = vmnorm(S->n, svaf, ewt) / S->tesco[2 + (S->nq - 1) * 3];
                double exup = 1.0 / (double)(S->l + 1);
                rhup = 1.0 / (1.4 * pow(dup, exup) + 0.0000014);
            }

            // 540
            double exsm = 1.0 / (double)S->l;
            double rhsm = 1.0 / (1.2 * pow(dsm, exsm) + 0.0000012);

            double rhdn = 0.0;
            if (S->nq != 1)
            {
                double ddn = vmnorm(S->n, &yh[(S->l - 1) * S->nyh], ewt) / S->tesco[0 + (S->nq - 1) * 3];
                double exdn = 1.0 / (double)S->nq;
                rhdn = 1.0 / (1.3 * pow(ddn, exdn) + 0.0000013);
            }

            // 550
            if (S->meth == 1) {
                double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
                if (S->l < S->lmax) { rhup = fmin(rhup, sm1[S->l - 1] / pdh); }
                rhsm = fmin(rhsm, sm1[S->nq - 1] / pdh);
                if (S->nq > 1) { rhdn = fmin(rhdn, sm1[S->nq - 2] / pdh); }
                S->pdest = 0.0;
            }

            int newq;
            double rh;
            // 560
            if (rhsm >= rhup)
            {
                // 570
                if (rhsm >= rhdn)
                {
                    newq = S->nq;
                    rh = rhsm;
                } else {
                    // 580
                    newq = S->nq - 1;
                    rh = rhdn;
                    if (S->kflag < 0 && rh > 1.0) { rh = 1.0; }
                }
            } else {
                if (rhup > rhdn) {
                    // 590
                    newq = S->l;
                    rh = rhup;
                    if (rh < 1.1) {
                        // 610 -> go to 700 (no rmax reset)
                        S->ialth = 3;
                        break;
                    }

                    double r = S->el[S->l - 1] / (double)S->l;
                    for (int i = 0; i < S->n; i++) {
                        yh[i + newq * nyh] = acor[i] * r;
                    }
                    // Apply step size change and continue
                    if (newq != S->nq) {
                        S->nq = newq;
                        S->l = S->nq + 1;
                        stoda_reset(S);
                    }
                    S->h = S->h * rh;
                    S->rc = S->rc * rh;
                    S->ialth = S->l;
                    continue;
                } else {
                    // 580 - Default case: rhsm < rhup AND rhup <= rhdn
                    newq = S->nq - 1;
                    rh = rhdn;
                    if (S->kflag < 0 && rh > 1.0) { rh = 1.0; }
                }
            }

            // 620
            // If meth = 1 and h is restricted by stability, bypass 10 percent test.
            if (S->meth == 1) {
                double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
                if (rh * pdh * 1.00001 >= sm1[newq - 1])
                {
                    // goto 625
                    ;
                } else if (S->kflag == 0 && rh < 1.1) {
                    // 610 -> go to 700 (no rmax reset)
                    S->ialth = 3;
                    break;
                }
            } else {
                if (S->kflag == 0 && rh < 1.1)
                {
                    // 610 -> go to 700 (no rmax reset)
                    S->ialth = 3;
                    break;
                }
            }

            // 625
            if (S->kflag <= -2) { rh = fmin(rh, 0.2); }

            // If there is a change of order, reset nq, l, and the coefficients.
            // in any case h is reset according to rh and the yh array is rescaled.
            // then exit from 690 if the step was ok, or redo the step otherwise.
            if (newq != S->nq) {
                // 630: nq = newq, l = nq + 1, iret = 2, go to 150 -> 170
                S->nq = newq;
                S->l = S->nq + 1;
                stoda_reset(S);
            }

            // 170 is inlined here to skip over the 160-170 block in stoda_adjust_step_size.
            rh = fmax(rh, S->hmin / fabs(S->h));
            rh = fmin(rh, S->rmax);
            rh = rh / fmax(1.0, fabs(S->h) * S->hmxi * rh);

            // Adams stability region constraint (meth = 1 only)
            if (S->meth == 1) {
                S->irflag = 0;
                double pdh = fmax(fabs(S->h) * S->pdlast, 0.000001);
                if (rh * pdh * 1.00001 < sm1[S->nq - 1]) {
                    rh = sm1[S->nq - 1] / pdh;
                    S->irflag = 1;
                }
            }

            // 178-180: Rescale yh array and apply step size
            double r = 1.0;
            for (int j = 2; j <= S->l; j++) {
                r = r * rh;
                for (int i = 0; i < S->n; i++) {
                    yh[i + (j-1)*nyh] = yh[i + (j-1)*nyh] * r;
                }
            }
            S->h = S->h * rh;
            S->rc = S->rc * rh;
            S->ialth = S->l;
            continue;

        } else if ((S->ialth > 1) || (S->l == S->lmax)) {
            break;
        } else {
            for (int i = 0; i < S->n; i++)
            {
                yh[i + (S->lmax - 1)*nyh] = acor[i];
            }
            break;
        }
    }

    // 700
    double r = 1.0 / S->tesco[1 + (S->nqu - 1)*3];
    for (int i = 0; i < S->n; i++)
    {
        acor[i] = acor[i] * r;
    }
    S->hold = S->h;
    S->jstart = 1;
    return;
}

/**
 *********************************************************************************
 *                         end of stoda.f implementation
 *********************************************************************************
 */


static void
intdy(const double t, const int k, double* yh, const int nyh, double* dky, int* iflag, lsoda_common_struct_t* S)
{
    double hu = 0.0;
    *iflag = 0;
    if ((k < 0) || (k > S->nq)) { *iflag = -1; return; }
    double tp = S->tn - hu - 100.0 * S->uround * (S->tn + hu);
    if ((t-tp)*(t- S->tn) > 0.0) { *iflag = -2; return; }

    double s = (t - S->tn) / S->h;
    int ic = 1;
    if (k != 0)
    {
        for (int jj = S->l - k; jj <= S->nq; jj++) { ic = ic * jj; }
    }
    double c = (double)ic;
    for (int i = 0; i < S->n; i++) { dky[i] = c*yh[i + (S->l - 1)*nyh]; }

    if (k != S->nq)
    {
        for (int jb = 1; jb <= S->nq - k; jb++)
        {
            if (k != 0)
            {
                for (int jj = S->nq - jb + 1 - k; jj <= S->nq - jb; jj++) { ic = ic * jj; }
            }
            // 35
            c = (double)ic;
            for (int i = 0; i < S->n; i++)
            {
                dky[i] = c * yh[i + (S->nq - jb)*nyh] + s * dky[i];
            }
            // 40
        }
        // 50
        if (k == 0) { return; }
    }
    // 55
    double r = pow(S->h, -k);
    for (int i = 0; i < S->n; i++) { dky[i] = r*dky[i]; }

    return;

}


/**
 *
 * This subroutine sets the error weight vector ewt according to
 *      ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
 *  with the subscript on rtol and/or atol possibly replaced by 1 above,
 *  depending on the value of itol.
 *
 */
static void
ewset(const int n, const int itol, double* rtol, double* atol, double* ycur, double* ewt)
{
    switch (itol)
    {
        case 1:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = *rtol * fabs(ycur[i]) + *atol;
            }
            break;

        case 2:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = *rtol * fabs(ycur[i]) + atol[i];
            }
            break;

        case 3:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = rtol[i] * fabs(ycur[i]) + *atol;
            }
            break;
        case 4:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = rtol[i] * fabs(ycur[i]) + atol[i];
            }
            break;
    }

    return;
}


void lsoda_mark_error(int* istate, int* illin)
{
    *istate = -3;  // Set istate to indicate an error
    (*illin)++;
    if (*illin == 5) {
        *istate = -8;  // Set istate to a more severe error after multiple issues
    }
}


void lsoda(
    lsoda_func_t f, int neq, double* restrict y, double* t, double* tout, int itol, double* rtol, double* atol,
    int* itask, int* istate, int* iopt, double* restrict rwork, int lrw, int* restrict iwork, int liw, lsoda_jac_t jac, const int jt)
{

    const int mord[2] = {12, 5};  // Maximum orders for Adams and BDF methods
    const int mxstp0 = 500;       // Default maximum steps
    const int mxhnl0 = 10;        // Default maximum nil step warnings

    lsoda_common_struct_t S = {0};  // Initialize state structure
    int n = 0, jtyp, ml = 0, mu = 0;
    int ixpr, mxstep, mxhnil, mxordn = 12, mxords = 5;
    double hmxi, hmin, h0 = 0.0, hmax, tsw = 0.0;
    int len1n, len1s, lwm, lenwm, len1c, len1, len2, lenrw, lenrwn, lenrws, lenrwc;
    int liwm, leniw, leniwc;
    int lewt = 0, lsavf, lacor, insufr, insufi;
    double rtoli, atoli;
    int i;
    int ntrep = 0;

    // Phase 1: Input Validation and Initial Dispatch
    // Block A: Test istate and itask for legality and branch appropriately

    if (*istate < 1 || *istate > 3) { lsoda_mark_error(istate, &S.illin); return; }
    if (*itask < 1 || *itask > 5) { lsoda_mark_error(istate, &S.illin); return; }

    if (*istate != 2)
    {
        // Initial state dispatch
        if (*istate == 1) {
            S.init = 0;
            if (*tout == *t) {
                ntrep++;
                if (ntrep < 5) return;
                *istate = -8;
                printf("Repeated tout equal to t: %f\n", *tout);
                return;
            }
            ntrep = 0;
        } else {  // *istate == 3
            if (S.init == 0) { lsoda_mark_error(istate, &S.illin); return; }
            ntrep = 0;
            // Continue to Block B (parameter changes)
        }

        // Check legality of non-optional inputs neq, itol, iopt, jt, ml, and mu
        if (neq <= 0) { lsoda_mark_error(istate, &S.illin); return; }

        if (*istate == 1) {
            n = neq;
        } else {
            if (neq > S.n) { lsoda_mark_error(istate, &S.illin); return; }
            n = neq;
        }

        if (itol < 1 || itol > 4) { lsoda_mark_error(istate, &S.illin); return; }
        if (*iopt < 0 || *iopt > 1) { lsoda_mark_error(istate, &S.illin); return; }

        if (jt == 3 || jt < 1 || jt > 5) { lsoda_mark_error(istate, &S.illin); return; }
        jtyp = jt;
        if (jt > 2) {
            ml = iwork[0];
            mu = iwork[1];
            if (ml < 0 || ml >= n) { lsoda_mark_error(istate, &S.illin); return; }
            if (mu < 0 || mu >= n) { lsoda_mark_error(istate, &S.illin); return; }
        }

        // Process and check optional inputs
        if (*iopt == 1) {
            ixpr = iwork[4];
            if (ixpr < 0 || ixpr > 1) { lsoda_mark_error(istate, &S.illin); return; }

            mxstep = iwork[5];
            if (mxstep < 0) { lsoda_mark_error(istate, &S.illin); return; }
            if (mxstep == 0) mxstep = mxstp0;

            mxhnil = iwork[6];
            if (mxhnil < 0) { lsoda_mark_error(istate, &S.illin); return; }
            if (mxhnil == 0) mxhnil = mxhnl0;

            if (*istate == 1) {
                h0 = rwork[4];
                mxordn = iwork[7];
                if (mxordn < 0)  { lsoda_mark_error(istate, &S.illin); return; }
                if (mxordn == 0) mxordn = 100;
                mxordn = (mxordn < mord[0]) ? mxordn : mord[0];
                mxords = iwork[8];
                if (mxords < 0) { lsoda_mark_error(istate, &S.illin); return; }
                if (mxords == 0) mxords = 100;
                mxords = (mxords < mord[1]) ? mxords : mord[1];
                if ((*tout - *t) * h0 < 0.0) { lsoda_mark_error(istate, &S.illin); return; }
            }

            hmax = rwork[5];
            if (hmax < 0.0) { lsoda_mark_error(istate, &S.illin); return; }
            hmxi = 0.0;
            if (hmax > 0.0) hmxi = 1.0/hmax;

            hmin = rwork[6];
            if (hmin < 0.0) { lsoda_mark_error(istate, &S.illin); return; }
        } else {
            // Default optional inputs
            ixpr = 0;
            mxstep = mxstp0;
            mxhnil = mxhnl0;
            hmxi = 0.0;
            hmin = 0.0;
            if (*istate == 1) {
                h0 = 0.0;
                mxordn = mord[0];
                mxords = mord[1];
            }
        }

        // Workspace setup and length checking
        if (*istate == 1) S.meth = 1;
        if (*istate == 1) S.nyh = n;
        S.lyh = 20;
        len1n = 20 + (mxordn + 1) * S.nyh;
        len1s = 20 + (mxords + 1) * S.nyh;
        lwm = len1s;
        if (jt <= 2) { lenwm = n * n + 2; }
        if (jt >= 4) { lenwm = (2 * ml + mu + 1) * n + 2; }
        len1s = len1s + lenwm;
        len1c = len1n;
        if (S.meth == 2) len1c = len1s;
        len1 = (len1n > len1s) ? len1n : len1s;
        len2 = 3 * n;
        lenrw = len1 + len2;
        lenrwn = len1n + len2;
        lenrws = len1s + len2;
        lenrwc = len1c + len2;
        iwork[16] = lenrw;  // Store required rwork length
        liwm = 1;
        leniw = 20 + n;
        leniwc = 20;
        if (S.meth == 2) leniwc = leniw;
        iwork[16] = leniw;  // Store required iwork length

        // Check workspace lengths
        if (*istate == 1 && lrw < lenrwc) { lsoda_mark_error(istate, &S.illin); return; }
        if (*istate == 1 && liw < leniwc) { lsoda_mark_error(istate, &S.illin); return; }
        if (*istate == 3 && lrw < lenrwc) {
            *istate = -7;
            return;
        }
        if (*istate == 3 && liw < leniwc) {
            *istate = -7;
            return;
        }

        lewt = len1;
        insufr = 0;
        if (lrw >= lenrw) {
            // Sufficient workspace
            lsavf = lewt + n;
            lacor = lsavf + n;
        } else {
            // Insufficient but will proceed with warning
            insufr = 2;
            lewt = len1c;
            lsavf = lewt + n;
            lacor = lsavf + n;
            // Note: xerrwv warning calls omitted as requested
        }

        insufi = 0;
        if (liw < leniw) {
            insufi = 2;
            // Note: xerrwv warning calls omitted as requested
        }

        // Check rtol and atol for legality
        rtoli = rtol[0];
        atoli = atol[0];
        for (i = 0; i < n; i++) {
            if (itol >= 3) rtoli = rtol[i];
            if (itol == 2 || itol == 4) atoli = atol[i];
            if (rtoli < 0.0) { lsoda_mark_error(istate, &S.illin); return; }
            if (atoli < 0.0) { lsoda_mark_error(istate, &S.illin); return; }
        }

        // Handle different initialization paths based on istate
        if (*istate == 1) {
            // Block C: First call initialization
            // Store required workspace lengths
            iwork[16] = lenrw;
            iwork[17] = leniw;

            // Initialize state variables
            S.uround = 2.22e-16;  // Machine precision (d1mach(4) equivalent)
            S.tn = *t;
            tsw = *t;
            S.maxord = mxordn;

            // Handle critical time for itask 4 and 5
            double tcrit = 0.0;
            if (*itask == 4 || *itask == 5) {
                tcrit = rwork[0];  // 0-based indexing
                if ((tcrit - *tout) * (*tout - *t) < 0.0)  { lsoda_mark_error(istate, &S.illin); return; }
                if (h0 != 0.0 && (*t + h0 - tcrit) * h0 > 0.0) {
                    h0 = tcrit - *t;
                }
            }

            // Initialize counters and flags
            S.jstart = 0;
            S.nhnil = 0;
            S.nst = 0;
            S.nje = 0;
            S.nslast = 0;
            S.hu = 0.0;
            S.nqu = 0;
            S.mused = 0;
            S.miter = 0;
            S.ccmax = 0.3;
            S.maxcor = 3;
            S.msbp = 20;
            S.mxncf = 10;

            // Initial call to f (lf0 points to yh(*,2))
            int lf0 = S.lyh + S.nyh;
            f(&neq, t, y, &rwork[lf0]);
            if (neq == -1) { return; }

            S.nfe = 1;

            // Load initial value vector in yh
            for (i = 0; i < n; i++) { rwork[i + S.lyh - 1] = y[i]; }

            // Load and invert the ewt array (h is temporarily set to 1.0)
            S.nq = 1;
            S.h = 1.0;

            // Call ewset (this function needs to be available)
            ewset(n, itol, rtol, atol, &rwork[S.lyh], &rwork[lewt]);

            for (i = 0; i < n; i++) {
                if (rwork[i + lewt] <= 0.0)  { lsoda_mark_error(istate, &S.illin); return; }
                rwork[i + lewt] = 1.0 / rwork[i + lewt];
            }

            // Compute initial step size h0 if not provided
            if (h0 == 0.0) {
                double tdist = fabs(*tout - *t);
                double w0 = fmax(fabs(*t), fabs(*tout));
                if (tdist < 2.0 * S.uround * w0)  { lsoda_mark_error(istate, &S.illin); return; }

                double tol = rtol[0];
                if (itol <= 2) {
                    // Find maximum rtol
                    for (i = 0; i < n; i++) {
                        if (itol >= 3) tol = fmax(tol, rtol[i]);
                    }
                }

                if (tol <= 0.0) {
                    atoli = atol[0];
                    for (i = 0; i < n; i++) {
                        if (itol == 2 || itol == 4) atoli = atol[i];
                        double ayi = fabs(y[i]);
                        if (ayi != 0.0) tol = fmax(tol, atoli / ayi);
                    }
                }

                tol = fmax(tol, 100.0 * S.uround);
                tol = fmin(tol, 0.001);

                // Call vmnorm to compute sum
                double sum = vmnorm(n, &rwork[lf0], &rwork[lewt]);
                sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
                h0 = 1.0 / sqrt(sum);
                h0 = fmin(h0, tdist);
                h0 = copysign(h0, *tout - *t);
            }

            // Adjust h0 if necessary to meet hmax bound
            double rh = fabs(h0) * hmxi;
            if (rh > 1.0) { h0 = h0 / rh; }

            // Load h with h0 and scale yh(*,2) by h0
            S.h = h0;
            for (i = 0; i < n; i++) { rwork[i + lf0 - 1] = h0 * rwork[i + lf0 - 1]; }
        }
    }

    // This is where istate == 2 jumps to after the initial setup and ready for label 200.
    // istate == 1 is now ready to continue to label 270
    // istate == 3 is ready to continue to label 200
    // Hence only istate == 2 and istate == 3 will enter the following label 200 logic.

    if (*istate != 1) {
        // Label 200: Block D - Task-specific pre-checks for continuation calls (istate == 2,3)
        S.nslast = S.nst;  // Save step count for later comparison

        // Handle Block D task-specific pre-checks for continuation calls
        switch (*itask) {
            case 1:
                if ((S.tn - *tout) * S.h < 0.0) {
                    // Continue to main loop
                    break;
                }
                // Call intdy for interpolation
                int iflag;
                intdy(*tout, 0, &rwork[S.lyh-1], S.nyh, y, &iflag, &S);
                if (iflag != 0)  { lsoda_mark_error(istate, &S.illin); return; }
                *t = *tout;
                goto success_exit;

            case 2:
                // One step mode - continue to main loop
                break;

            case 3: {
                double tp = S.tn - S.hu * (1.0 + 100.0 * S.uround);
                if ((tp - *tout) * S.h > 0.0)  { lsoda_mark_error(istate, &S.illin); return; }
                if ((S.tn - *tout) * S.h < 0.0) {
                    // Continue to main loop
                    break;
                }
                *t = S.tn;
                goto success_exit;
            }

            case 4: {
                double tcrit = rwork[0];
                if ((S.tn - tcrit) * S.h > 0.0)  { lsoda_mark_error(istate, &S.illin); return; }
                if ((tcrit - *tout) * S.h < 0.0)  { lsoda_mark_error(istate, &S.illin); return; }
                if ((S.tn - *tout) * S.h < 0.0) {
                    // Continue to main loop after setting tcrit handling
                    double hmx = fabs(S.tn) + fabs(S.h);
                    int ihit = fabs(S.tn - tcrit) <= 100.0 * S.uround * hmx;
                    if (!ihit) {
                        double tnext = S.tn + S.h * (1.0 + 4.0 * S.uround);
                        if ((tnext - tcrit) * S.h <= 0.0) {
                            // Continue to main loop
                            break;
                        }
                        S.h = (tcrit - S.tn) * (1.0 - 4.0 * S.uround);
                        if (*istate == 2 && S.jstart >= 0) S.jstart = -2;
                    }
                    break;
                }
                // Call intdy for interpolation
                int iflag;
                intdy(*tout, 0, &rwork[S.lyh-1], S.nyh, y, &iflag, &S);
                if (iflag != 0)  { lsoda_mark_error(istate, &S.illin); return; }
                *t = *tout;
                goto success_exit;
            }

            case 5: {
                double tcrit = rwork[0];
                if ((S.tn - tcrit) * S.h > 0.0)  { lsoda_mark_error(istate, &S.illin); return; }
                double hmx = fabs(S.tn) + fabs(S.h);
                int ihit = fabs(S.tn - tcrit) <= 100.0 * S.uround * hmx;
                if (ihit) {
                    *t = tcrit;
                    goto success_exit;
                }
                double tnext = S.tn + S.h * (1.0 + 4.0 * S.uround);
                if ((tnext - tcrit) * S.h <= 0.0) {
                    // Continue to main loop
                    break;
                }
                S.h = (tcrit - S.tn) * (1.0 - 4.0 * S.uround);
                if (*istate == 2 && S.jstart >= 0) S.jstart = -2;
                break;
            }
        }
    }

    // Label 250: Main integration loop - This is the convergence point for all cases
    while (1) {
        if (*istate != 1) {
            // Check if method has switched since last iteration
            if (S.meth != S.mused) {
                if (S.insufr == 1) {
                    *istate = -5;  // Insufficient rwork length
                    goto error_exit;
                }
                if (S.insufi == 1) {
                    *istate = -6;  // Insufficient iwork length
                    goto error_exit;
                }
            }

            // Check for too many steps
            if ((S.nst - S.nslast) >= S.mxstep) {
                *istate = -1;  // Too many steps
                goto error_exit;
            }

            // Set error weights (skip for istate=1 as already done in Block C)
            ewset(n, itol, rtol, atol, &rwork[S.lyh], &rwork[lewt]);
            for (i = 0; i < n; i++) {
                if (rwork[i + lewt] <= 0.0) {
                    *istate = -6;
                    goto error_exit;
                }
                rwork[i + lewt] = 1.0 / rwork[i + lewt];
            }
        }

        // Label 270: Tolerance check (this is where istate=1 enters after Block C)
        double tolsf = S.uround * vmnorm(n, &rwork[S.lyh], &rwork[lewt]);
        if (tolsf > 0.01) {
            tolsf = tolsf * 200.0;
            if (S.nst == 0) {
                // At start of problem, too much accuracy requested
                rwork[13] = tolsf;  // Store tolsf for user information
                lsoda_mark_error(istate, &S.illin);
                goto error_exit;
            }
            // Too much accuracy requested at current t
            *istate = -2;
            rwork[13] = tolsf;  // Store tolsf for user information
            *t = S.tn;
            goto error_exit;
        }

        // Check for h too small
        if ((S.tn + S.h) == S.tn) {
            S.nhnil++;
            if (S.nhnil <= S.mxhnil) {
                // Warning message would go here if ixpr == 1
            } else {
                *istate = -4;  // h is effectively 0
                goto error_exit;
            }
        }

        // Call the core integrator stoda
        stoda(&neq, y, &rwork[S.lyh], S.nyh, &rwork[S.lyh], &rwork[lewt],
              &rwork[S.lsavf], &rwork[S.lacor], &rwork[S.lwm], &iwork[S.liwm],
              f, jac, &S);

        int kgo = 1 - S.kflag;  // Convert FORTRAN logic

        // Process return from stoda
        switch (kgo) {
            case 1:  // kflag = 0, successful step
                // Continue with task-specific post-step processing
                break;

            case 2:  // kflag = -1, step size too small
                *istate = -4;
                goto error_exit;

            case 3:  // kflag = -2, too much work
                *istate = -1;
                goto error_exit;

            default:
                // This shouldn't happen
                *istate = -7;
                goto error_exit;
        }

        // Task-specific post-step processing
        switch (*itask) {
            case 1:
                if ((S.tn - *tout) * S.h < 0.0) {
                    continue;  // Go back to 250 (main loop)
                }
                // Call intdy for interpolation
                int iflag;
                intdy(*tout, 0, &rwork[S.lyh-1], S.nyh, y, &iflag, &S);
                if (iflag != 0) {
                    lsoda_mark_error(istate, &S.illin);
                    goto error_exit;
                }
                *t = *tout;
                goto success_exit;

            case 2:
                // One step mode - always exit after one step
                goto success_exit;

            case 3: {
                double tp = S.tn - S.hu * (1.0 + 100.0 * S.uround);
                if ((tp - *tout) * S.h > 0.0) {
                    lsoda_mark_error(istate, &S.illin);
                    goto error_exit;
                }
                if ((S.tn - *tout) * S.h < 0.0) {
                    continue;  // Go back to 250
                }
                *t = S.tn;
                goto success_exit;
            }

            case 4: {
                double tcrit = rwork[0];
                if ((S.tn - tcrit) * S.h > 0.0) {
                    lsoda_mark_error(istate, &S.illin);
                    goto error_exit;
                }
                double hmx = fabs(S.tn) + fabs(S.h);
                int ihit = fabs(S.tn - tcrit) <= 100.0 * S.uround * hmx;
                if (ihit) {
                    *t = tcrit;
                    goto success_exit;
                }
                double tnext = S.tn + S.h * (1.0 + 4.0 * S.uround);
                if ((tnext - tcrit) * S.h <= 0.0) {
                    continue;  // Go back to 250
                }
                S.h = (tcrit - S.tn) * (1.0 - 4.0 * S.uround);
                if (S.jstart >= 0) S.jstart = -2;
                continue;  // Go back to 250
            }

            case 5: {
                double tcrit = rwork[0];
                double hmx = fabs(S.tn) + fabs(S.h);
                int ihit = fabs(S.tn - tcrit) <= 100.0 * S.uround * hmx;
                if (ihit) {
                    *t = tcrit;
                    goto success_exit;
                }
                // Call intdy for interpolation
                int iflag;
                intdy(*tout, 0, &rwork[S.lyh-1], S.nyh, y, &iflag, &S);
                if (iflag != 0) {
                    lsoda_mark_error(istate, &S.illin);
                    goto error_exit;
                }
                *t = *tout;
                goto success_exit;
            }
        }
    }

success_exit:
    // Block G: Successful return preparation
    for (i = 0; i < n; i++) {
        y[i] = rwork[i + S.lyh - 1];
    }
    *t = S.tn;
    if (*itask == 4 || *itask == 5) {
        // Handle tcrit hit for tasks 4 and 5
        double tcrit = rwork[0];
        double hmx = fabs(S.tn) + fabs(S.h);
        int ihit = fabs(S.tn - tcrit) <= 100.0 * S.uround * hmx;
        if (ihit) *t = tcrit;
    }

    *istate = 2;
    S.illin = 0;

    // Set optional outputs
    rwork[10] = S.hu;
    rwork[11] = S.h;
    rwork[12] = S.tn;
    rwork[14] = tsw;
    iwork[10] = S.nst;
    iwork[11] = S.nfe;
    iwork[12] = S.nje;
    iwork[13] = S.nqu;
    iwork[14] = S.nq;
    iwork[18] = S.mused;
    iwork[19] = S.meth;
    return;

error_exit:;  // Semi-colon to silence compiler warning
    // Handle error returns
    // Compute imxer (index of component with largest weighted local error) if relevant
    int imxer = 1;  // Default to first component (1-based for compatibility)
    if (S.lacor > 0 && lewt > 0) {  // Only if error vectors are available
        double big = 0.0;
        for (i = 0; i < n; i++) {
            double size = fabs(rwork[i + S.lacor] * rwork[i + lewt]);
            if (size > big) {
                big = size;
                imxer = i + 1;  // Convert to 1-based indexing for compatibility
            }
        }
    }
    iwork[15] = imxer;  // Store in iwork(16) using 0-based indexing

    // Set y vector and outputs
    for (i = 0; i < n; i++) {
        y[i] = rwork[i + S.lyh - 1];
    }
    *t = S.tn;
    S.illin = 0;

    // Set optional outputs
    rwork[10] = S.hu;
    rwork[11] = S.h;
    rwork[12] = S.tn;
    rwork[14] = tsw;
    iwork[10] = S.nst;
    iwork[11] = S.nfe;
    iwork[12] = S.nje;
    iwork[13] = S.nqu;
    iwork[14] = S.nq;
    iwork[18] = S.mused;
    iwork[19] = S.meth;
    return;

}
