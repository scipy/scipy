#include "vode.h"
#include <stdio.h>

static inline int int_min(const int a, const int b) { return a < b ? a : b; }
static inline int int_max(const int a, const int b) { return a > b ? a : b; }


/**
 * @brief Computes the weighted root-mean-square norm of a vector.
 *
 * This function computes the weighted RMS norm of the vector of length n
 * contained in the array v, with weights contained in the array w.
 *
 * Formula: dvnorm = sqrt((1/n) * sum((v[i]*w[i])^2))
 *
 * @param n  Length of the vectors
 * @param v  Input vector
 * @param w  Weight vector
 * @return   The weighted RMS norm
 */
static double
dvnorm(const int n, double* restrict v, double* restrict w)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        double vw = v[i] * w[i];
        sum += vw * vw;
    }
    return sqrt(sum / (double)n);
}


/**
 * @brief Sets method coefficients for Adams or BDF integration formulas.
 *
 * DVSET is called by DVSTEP and sets coefficients for use there. For each order NQ,
 * the coefficients in EL are calculated by use of the generating polynomial lambda(x),
 * with coefficients EL(i): lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
 *
 * For BDF methods:
 *   lambda(x) = (1 + x/xi(NQ)) * product_{i=1}^{NQ-1} (1 + x/xi(i))
 *
 * For Adams methods:
 *   (d/dx) lambda(x) = c * product_{i=1}^{NQ-1} (1 + x/xi(i))
 *   lambda(-1) = 0, lambda(0) = 1, where c is a normalization constant.
 *
 * In both cases, xi(i) is defined by: H*xi(i) = t_n - t_{n-i} = H + TAU(1) + ... + TAU(i-1).
 *
 * The coefficients are stored in:
 *   - EL: Vector of length 13 storing coefficients for the corrector formula
 *   - TQ: Vector of length 5 storing constants for convergence test, error test,
 *         and selection of H at a new order
 *
 * Communication uses:
 *   - TAU: Vector of length 13 containing past NQ values of H
 *   - METH: Basic method indicator (1=Adams, 2=BDF)
 *   - NQ: Current order
 *   - L: NQ + 1, length of vector stored in EL
 *   - NQWAIT: Counter controlling frequency of order changes (order change considered if NQWAIT=1)
 *
 * @param S  Pointer to VODE common struct containing method parameters
 */
static void
dvset(vode_common_struct_t* S)
{

    double em[13], cortes = 0.1;
    double flotl = (double)(S->l);
    switch (S->meth) {
        case 1: {
            // Set coefficients for Adams methods.
            if (S->nq == 1) {
                S->el[0] = 1.0;
                S->el[1] = 1.0;
                S->tq[0] = 1.0;
                S->tq[1] = 2.0;
                S->tq[2] = 6.0 * S->tq[1];
                S->tq[4] = 1.0;
                break;
            }

            double hsum = S->h;
            em[0] = 1.0;
            double flotnq = flotl - 1.0;
            for (int i = 1; i < S->l; i++) { em[i] = 0.0; }
            for (int j = 1; j < S->nq; j++) {
                if ((j == S->nq - 1) && (S->nqwait == 1)) {
                    double s = 1.0;
                    double csum = 0.0;
                    for (int i = 1; i < S->nq; i++) {
                        csum += s * em[i-1] / (double)(i + 1);
                        s = -s;
                    }
                    S->tq[0] = em[S->nq - 2] / (flotnq * csum);
                }
                double rxi = S->h / hsum;
                for (int iback = 1; iback <= j; iback++) {
                    int i = (j + 1) - iback;
                    em[i] += em[i - 1] * rxi;
                }
                hsum += S->tau[j - 1];
            }
            // Compute integral from -1 to 0 of polynomial and of x times it.
            double s = 1.0;
            double em0 = 0.0;
            double csum = 0.0;
            for (int i = 1; i <= S->nq; i++) {
                double floti = (double)(i);
                em0 += s * em[i - 1] / floti;
                csum += s * em[i - 1] / (floti + 1.0);
                s = -s;
            }
            // In EL, form coefficients of normalized integrated polynomial.
            s = 1.0 / em0;
            S->el[0] = 1.0;
            for (int i = 1; i <= S->nq; i++) {
                S->el[i] = s * em[i - 1] / (double)(i);
            }
            double xi = hsum / S->h;
            S->tq[1] = xi * em0 / csum;
            S->tq[4] = xi / S->el[S->l - 1];
            if (S->nqwait != 1) {
                break;
            }
            // For higher order control constant, multiply polynomial by 1+x/xi(q).
            double rxi = 1.0 / xi;
            for (int iback = 1; iback <= S->nq; iback++) {
                int i = (S->l) - iback;
                em[i] += em[i - 1] * rxi;
            }
            // Compute integral of polynomial.
            s = 1.0;
            csum = 0.0;
            for (int i = 1; i <= S->l; i++) {
                csum += s * em[i - 1] / (double)(i + 1);
                s = -s;
            }
            S->tq[2] = flotl * em0 / csum;
            break;
        }
        case 2: {
            // Set coefficients for BDF methods.
            for (int j = 2; j < S->l; j++) {
                S->el[j] = 0.0;
            }
            S->el[0] = 1.0;
            S->el[1] = 1.0;
            double alph0 = -1.0;
            double ahatn0 = -1.0;
            double hsum = S->h;
            double rxi = 1.0;
            double rxis = 1.0;

            // For NQ=1, skip coefficient computation (Fortran line 2537: IF (NQ .EQ. 1) GO TO 240)
            if (S->nq != 1) {
                for (int j = 1; j < S->nq - 1; j++) {
                    // In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)).
                    hsum += S->tau[j - 1];
                    rxi = S->h / hsum;
                    alph0 -= 1.0 / (double)(j + 1);
                    for (int iback = 1; iback <= j + 1; iback++) {
                        int i = (j + 2) - iback;
                        S->el[i] += S->el[i - 1] * rxi;
                    }
                }
                // 230

                alph0 -= 1.0 / (double)(S->nq);
                rxis = -S->el[1] - alph0;
                hsum += S->tau[S->nq - 2];
                rxi = S->h / hsum;
                ahatn0 = -S->el[1] - rxi;
                for (int iback = 1; iback <= S->nq; iback++) {
                    int i = (S->nq + 1) - iback;
                    S->el[i] += S->el[i - 1] * rxis;
                }
            }
            double t1 = 1.0 - ahatn0 + alph0;
            double t2 = 1.0 + (double)(S->nq) * t1;
            S->tq[1] = fabs(alph0 * t2 / t1);
            S->tq[4] = fabs(t2 / (S->el[S->l - 1] * rxi / rxis));
            if (S->nqwait != 1) { break; }
            double cnqm1 = rxis / S->el[S->l - 1];
            double t3 = alph0 + 1.0 / (double)(S->nq);
            double t4 = ahatn0 + rxi;
            double elp = t3 / (1.0 - t4 + t3);
            S->tq[0] = fabs(elp / cnqm1);
            hsum += S->tau[S->nq - 1];
            rxi = S->h / hsum;
            double t5 = alph0 - 1.0 / (double)(S->nq + 1);
            double t6 = ahatn0 - rxi;
            elp = t2 / (1.0 - t6 + t5);
            S->tq[2] = fabs(elp * rxi * (flotl + 1.0) * t5);
            break;
        }
    }

    S->tq[3] = cortes * S->tq[1];

    return;
}


/**
 * @brief Adjusts the Nordsieck array YH when changing integration order.
 *
 * DVJUST adjusts the YH array on reduction of order, and also when the order is
 * increased for the stiff option (METH = 2). The adjustment involves updating the
 * Nordsieck history array to reflect the new polynomial degree.
 *
 * For nonstiff methods (METH=1):
 *   - Order increase: Zeros out the new column in YH
 *   - Order decrease: Reconstructs coefficients and subtracts correction terms
 *
 * For stiff methods (METH=2):
 *   - Order increase: Constructs new column using BDF formula and adds correction terms
 *   - Order decrease: Reconstructs BDF coefficients and subtracts correction terms
 *
 * Communication with DVJUST uses:
 *   - IORD: Integer flag for order change direction (+1 = increase, -1 = decrease)
 *   - HSCAL: Step size H used in scaling of Nordsieck array YH
 *           (If IORD = +1, DVJUST assumes that HSCAL = TAU(1))
 *
 * Special case: If NQ = 2 and IORD != 1, the routine returns immediately without
 * making any adjustments.
 *
 * @param yh    Nordsieck history array (n x lmax)
 * @param ldyh  Leading dimension of yh array
 * @param iord  Order change indicator: +1 for increase, -1 for decrease
 * @param S     Pointer to VODE common struct
 */
static void
dvjust(double* restrict yh, const int ldyh, const int iord, vode_common_struct_t* S)
{
    if ((S->nq == 2) && (iord != 1)) {
        return;
    }

    switch (S->meth) {
        case 1: {
            // Nonstiff option
            // Check to see if the order is being increased or decreased.
            if (iord == 1) {
                // Order increase
                // Zero out next column in YH array.
                for (int i = 0; i < S->n; i++) {
                    yh[i + (S->l) * ldyh] = 0.0;
                }
                return;
            } else {
                // Order decrease
                for (int j = 0; j < S->lmax; j++) {
                    S->el[j] = 0.0;
                }
                S->el[1] = 1.0;
                double hsum = 0.0;
                for (int j = 0; j < S->nq - 2; j++) {
                    // Construct coefficients of x*(x+xi(1))*...*(x+xi(j)).
                    hsum += S->tau[j];
                    double xi = hsum / S->hscal;
                    for (int iback = 0; iback < j + 2; iback++) {
                        int i = (j + 2) - iback;
                        S->el[i] = S->el[i] * xi + S->el[i - 1];
                    }
                }
                // Construct coefficients of integrated polynomial.
                for (int j = 1; j < S->nq - 1; j++) {
                    S->el[j + 1] = (double)(S->nq) * S->el[j] / (double)(j + 1);
                }
                // Subtract correction terms from YH array.
                for (int j = 2; j < S->nq; j++) {
                    for (int i = 0; i < S->n; i++) {
                        yh[i + j * ldyh] -= yh[i + (S->l - 1) * ldyh] * S->el[j];
                    }
                }
            }
            break;
        }
        case 2: {
            // Stiff option
            // Check to see if the order is being increased or decreased.
            if (iord == 1) {
                // Order increase
                for (int j = 0; j < S->lmax; j++) {
                    S->el[j] = 0.0;
                }
                S->el[2] = 1.0;
                double alph0 = -1.0;
                double alph1 = 1.0;
                double prod = 1.0;
                double xiold = 1.0;
                double hsum = S->hscal;
                if (S->nq != 1) {
                    for (int j = 0; j < S->nq - 1; j++) {
                        // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
                        int jp1 = j + 1;
                        hsum += S->tau[jp1];
                        double xi = hsum / S->hscal;
                        prod *= xi;
                        alph0 -= 1.0 / (double)(j + 2);
                        alph1 += 1.0 / xi;
                        for (int iback = 0; iback < j + 2; iback++) {
                            int i = (j + 3) - iback;
                            S->el[i] = S->el[i] * xiold + S->el[i - 1];
                        }
                        xiold = xi;
                    }
                }
                double t1 = (-alph0 - alph1) / prod;
                // Load column L + 1 in YH array.
                int lp1 = S->l + 1;
                for (int i = 0; i < S->n; i++) {
                    yh[i + (lp1 - 1) * ldyh] = t1 * yh[i + (S->lmax - 1) * ldyh];
                }
                // Add correction terms to YH array.
                int nqp1 = S->nq + 1;
                for (int j = 2; j < nqp1; j++) {
                    daxpy_(&S->n, &S->el[j], &yh[0 + (lp1 - 1) * ldyh], &(int){1}, &yh[0 + j * ldyh], &(int){1});
                }
            } else {
                // Order decrease
                for (int j = 0; j < S->lmax; j++) {
                    S->el[j] = 0.0;
                }
                S->el[2] = 1.0;
                double hsum = 0.0;
                for (int j = 0; j < S->nq - 2; j++) {
                    // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
                    hsum += S->tau[j];
                    double xi = hsum / S->hscal;
                    for (int iback = 0; iback < j + 2; iback++) {
                        int i = (j + 3) - iback;
                        S->el[i] = S->el[i] * xi + S->el[i - 1];
                    }
                }
                // Subtract correction terms from YH array.
                for (int j = 2; j < S->nq; j++) {
                    for (int i = 0; i < S->n; i++) {
                        yh[i + j * ldyh] -= yh[i + (S->l - 1) * ldyh] * S->el[j];
                    }
                }
            }
            break;
        }
    }

    return;
}


/**
 * @brief Computes an initial step size H0 for the first step in ZVODE.
 *
 * This routine computes the step size H0 to be attempted on the first step
 * when the user has not supplied a value. It uses an iterative procedure to
 * approximate the initial second derivative and defines h from the weighted
 * RMS norm condition: norm(h**2 * yddot / 2) = 1.
 * A bias factor of 1/2 is applied to the resulting h.
 *
 * @param n       Size of ODE system
 * @param t0      Initial value of independent variable
 * @param y0      Vector of initial conditions
 * @param ydot    Vector of initial first derivatives
 * @param f       User function for right-hand side f(t,y)
 * @param rpar    User real work array
 * @param ipar    User integer work array
 * @param tout    First output value of independent variable
 * @param uround  Machine unit roundoff
 * @param ewt     Error weight vector
 * @param itol    Tolerance type indicator
 * @param atol    Absolute tolerance (scalar or array)
 * @param y       Work array of length n
 * @param temp    Work array of length n
 * @param h0      Output: step size to be attempted
 * @param niter   Output: number of iterations performed
 * @param ier     Output: error flag (0=success, -1=tout too close to t0)
 */
static void
dvhin(
    int n,
    double t0,
    double* y0,
    double* ydot,
    vode_func_t f,
    double* rpar,
    int* ipar,
    double tout,
    double uround,
    double* ewt,
    int itol,
    double* atol,
    double* y,
    double* temp,
    double* h0,
    int* niter,
    int* ier
)
{
    double afi, atoli, delyi, h, half = 0.5, hg, hlb, hnew, hrat,
        hub, hun = 100.0, pt1 = 0.1, t1, tdist, tround, two = 2.0, yddnrm;
    int i, iter;

    *niter = 0;
    tdist = fabs(tout - t0);
    tround = uround * fmax(fabs(t0), fabs(tout));
    if (tdist < two * tround) {
        *ier = -1;
        return;
    }

    // Set a lower bound on h based on the roundoff level in T0 and TOUT.
    hlb = hun * tround;
    // Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT.
    hub = pt1 * tdist;
    atoli = atol[0];
    for (i = 0; i < n; i++) {
        if (itol == 2 || itol == 4) {
            atoli = atol[i];
        }
        delyi = pt1 * fabs(y0[i]) + atoli;
        afi = fabs(ydot[i]);
        if (afi * hub > delyi) {
            hub = delyi / afi;
        }
    }

    // Set initial guess for h as geometric mean of upper and lower bounds.
    iter = 0;
    hg = sqrt(hlb * hub);
    // If the bounds have crossed, exit with the mean value.
    if (hub < hlb) {
        *h0 = hg;
        *h0 = copysign(*h0, tout - t0);
        *niter = iter;
        *ier = 0;
        return;
    }

    // Looping point for iteration.
    while (1) {
        // Estimate the second derivative as a difference quotient in f.
        h = copysign(hg, tout - t0);
        t1 = t0 + h;
        for (i = 0; i < n; i++) {
            y[i] = y0[i] + h * ydot[i];
        }
        f(n, t1, y, temp, rpar, ipar);
        for (i = 0; i < n; i++) {
            temp[i] = (temp[i] - ydot[i]) / h;
        }
        yddnrm = dvnorm(n, temp, ewt);
        // Get the corresponding new value of h.
        if (yddnrm * hub * hub > two) {
            hnew = sqrt(two / yddnrm);
        } else {
            hnew = sqrt(hg * hub);
        }
        iter = iter + 1;

        // Test the stopping conditions.
        // Stop if the new and previous h values differ by a factor of .lt. 2.
        // Stop if four iterations have been done.  Also, stop with previous h
        // if HNEW/HG .gt. 2 after first iteration, as this probably means that
        // the second derivative value is bad because of cancellation error.
        if (iter >= 4) {
            break;
        }
        hrat = hnew / hg;
        if ((hrat > half) && (hrat < two)) {
            break;
        }
        if ((iter >= 2) && (hnew > two * hg)) {
            hnew = hg;
            break;
        }
        hg = hnew;
    }

    // Iteration done. Apply bounds, bias factor, and sign. Then exit.
    *h0 = hnew * half;
    *h0 = fmin(fmax(*h0, hlb), hub);
    *h0 = copysign(*h0, tout - t0);
    *niter = iter;
    *ier = 0;

    return;
}


/**
 * @brief Interpolation routine for computing derivatives of the solution.
 *
 * This routine computes interpolated values of the K-th derivative of the
 * dependent variable vector y, and stores it in DKY. The values are obtained
 * by interpolation using the Nordsieck history array YH.
 *
 * @param S      Pointer to VODE common struct
 * @param t      Value of independent variable where answers are desired
 * @param k      Order of derivative (0 <= k <= nq)
 * @param yh     Nordsieck history array
 * @param ldyh   Leading dimension of yh array
 * @param dky    Output array for k-th derivative values
 * @param iflag  Output flag: 0=success, -1=illegal k, -2=illegal t
 */
static void
dvindy(
    vode_common_struct_t* S,
    const double t,
    const int k,
    double* restrict yh,
    const int ldyh,
    double* restrict dky,
    int* iflag
)
{
    *iflag = 0;
    if ((k < 0) || (k > S->nq)) {
        *iflag = -1;
        return;
    }

    double tfuzz = 100.0 * S->uround * (S->tn + S->hu);
    double tp = S->tn - S->hu - tfuzz;
    double tn1 = S->tn + tfuzz;
    if ((t - tp) * (t - tn1) > 0.0) {
        *iflag = -2;
        return;
    }

    double s = (t - S->tn) / S->h;
    int ic = 1;
    if (k != 0) {
        for (int jj = S->l - k; jj <= S->nq; jj++) {
            ic = ic * jj;
        }
    }
    double c = (double)ic;
    for (int i = 0; i < S->n; i++) {
        dky[i] = c * yh[i + (S->l - 1) * ldyh];
    }

    if (k != S->nq) {
        int jb2 = S->nq - k;
        for (int jb = 1; jb <= jb2; jb++) {
            int j = S->nq - jb;
            int jp1 = j + 1;
            ic = 1;
            if (k != 0) {
                int jj1 = jp1 - k;
                for (int jj = jj1; jj <= j; jj++) {
                    ic = ic * jj;
                }
            }
            c = (double)ic;
            for (int i = 0; i < S->n; i++) {
                dky[i] = c * yh[i + (jp1 - 1) * ldyh] + s * dky[i];
            }
        }
        if (k == 0) { return; }
    }

    // Scale by H**(-K)
    double r = pow(S->h, -k);
    int int1 = 1;
    dscal_(&S->n, &r, dky, &int1);

    return;
}


/**
 * @brief Manages the solution of linear systems arising from chord iteration.
 *
 * This routine is called if miter != 0.
 * - If miter is 1 or 2, it calls dgetrs to solve the system
 * - If miter == 3, it updates the coefficient h*rl1 in the diagonal matrix
 * - If miter is 4 or 5, it calls dgbtrs
 *
 * @param S      Pointer to VODE common struct
 * @param wm     Work array containing LU decomposition or diagonal matrix
 * @param iwm    Integer work array with pivot information
 * @param x      Right-hand side on input, solution on output
 * @param iersl  Error flag: 0=success, 1=singular matrix (miter=3)
 */
static void
dvsol(vode_common_struct_t* S, double* restrict wm, int* restrict iwm, double* restrict x, int* iersl)
{
    *iersl = 0;
    int ier = 0;
    int int1 = 1;

    switch (S->miter)
    {
        case 1:
        case 2: {
            dgetrs_("N", &S->n, &int1, &wm[2], &S->n, &iwm[30], x, &S->n, &ier);
            break;
        }
        case 3: {
            double phrl1 = wm[1];
            double hrl1 = S->h * S->rl1;
            wm[1] = hrl1;
            if (hrl1 != phrl1) {
                double r = hrl1 / phrl1;
                for (int i = 0; i < S->n; i++) {
                    double di = 1.0 - r * (1.0 - 1.0 / wm[2 + i]);
                    if (fabs(di) == 0.0) {
                        *iersl = 1;
                        return;
                    }
                    wm[2 + i] = 1.0 / di;
                }
            }
            for (int i = 0; i < S->n; i++) {
                x[i] = wm[2 + i] * x[i];
            }
            break;
        }
        case 4:
        case 5: {
            int ml = iwm[0];
            int mu = iwm[1];
            int meband = 2 * ml + mu + 1;
            dgbtrs_("N", &S->n, &ml, &mu, &int1, &wm[2], &meband, &iwm[30], x, &S->n, &ier);
            break;
        }
    }

    return;
}


/**
 * @brief Sets the error weight vector EWT based on tolerance parameters.
 *
 * This subroutine sets the error weight vector EWT according to:
 *     EWT(i) = RTOL(i)*abs(YCUR(i)) + ATOL(i)
 * with the subscript on RTOL and/or ATOL possibly replaced by 1,
 * depending on the value of ITOL.
 *
 * @param n     Number of equations
 * @param itol  Tolerance type indicator (1-4)
 * @param rtol  Relative tolerance (scalar or array)
 * @param atol  Absolute tolerance (scalar or array)
 * @param ycur  Current solution vector
 * @param ewt   Output: error weight vector
 */
static void
dewset(const int n, const int itol, double* rtol, double* atol, double* ycur, double* ewt)
{
    switch (itol)
    {
        case 1:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = rtol[0] * fabs(ycur[i]) + atol[0];
            }
            break;

        case 2:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = rtol[0] * fabs(ycur[i]) + atol[i];
            }
            break;

        case 3:
            for (int i = 0; i < n; i++)
            {
                ewt[i] = rtol[i] * fabs(ycur[i]) + atol[0];
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


/**
 * @brief Computes and processes the Jacobian matrix for implicit ODE methods.
 *
 * DVJAC is called by DVNLSD to compute and process the matrix P = I - h*rl1*J,
 * where J is an approximation to the Jacobian. Here J is computed by the
 * user-supplied routine JAC if MITER = 1 or 4, or by finite differencing if
 * MITER = 2, 3, or 5. If MITER = 3, a diagonal approximation to J is used.
 *
 * If JSV = -1, J is computed from scratch in all cases. If JSV = 1 and
 * MITER = 1, 2, 4, or 5, and if the saved value of J is considered acceptable,
 * then P is constructed from the saved J.
 *
 * J is stored in wm and replaced by P. If MITER != 3, P is then subjected to
 * LU decomposition in preparation for later solution of linear systems with P
 * as coefficient matrix. This is done by DGETRF if MITER = 1 or 2, and by
 * DGBTRF if MITER = 4 or 5.
 *
 * @param y      Vector containing predicted values on entry
 * @param yh     Nordsieck history array
 * @param ldyh   Leading dimension of yh array (>= n)
 * @param ewt    Error weight vector of length n
 * @param ftem   Temporary array for function evaluations
 * @param savf   Array containing f evaluated at predicted y
 * @param wm     Real work space for matrices. On output contains the inverse
 *               diagonal matrix if MITER=3, or LU decomposition of P if MITER=1,2,4,5.
 *               Matrix elements start at wm[2]. Saved Jacobian starts at wm[locjs-1].
 *               wm[0] = sqrt(uround) for numerical Jacobian step.
 *               wm[1] = h*rl1, saved for later use if MITER=3.
 * @param iwm    Integer work space containing pivot information starting at iwm[30]
 *               if MITER=1,2,4,5. Also contains band parameters ml=iwm[0] and mu=iwm[1]
 *               if MITER=4 or 5.
 * @param f      User-supplied function for computing dy/dt
 * @param jac    User-supplied Jacobian function
 * @param ierpj  Output error flag: 0=success, 1=singular matrix detected
 * @param rpar   User real work array
 * @param ipar   User integer work array
 * @param S      Pointer to VODE common struct
 */
static void
dvjac(
    double* restrict y,
    double* restrict yh,
    const int ldyh,
    double* restrict ewt,
    double* restrict ftem,
    double* restrict savf,
    double* restrict wm,
    int* restrict iwm,
    vode_func_t f,
    vode_jac_t jac,
    int* ierpj,
    double* restrict rpar,
    int* restrict ipar,
    vode_common_struct_t* S
)
{
    *ierpj = 0;
    double hrl1 = S->h * S->rl1;

    // See whether J should be evaluated (JOK = -1) or not (JOK = 1).
    int jok = S->jsv;
    if (jok == 1) {
        if ((S->nst == 0) || (S->nst > S->nslj + S->msbj)) { jok = -1; }
        if ((S->icf == 1) && (S->drc < S->ccmxj)) { jok = -1; }
        if (S->icf == 2) { jok = -1; }
    }

    // End of setting JOK.
    if ((jok == -1) && (S->miter == 1)) {
        // If JOK = -1 and MITER = 1, call JAC for a new Jacobian.
        S->nje += 1;
        S->nslj = S->nst;
        S->jcur = 1;
        int lenp = S->n * S->n;
        for (int i = 0; i < lenp; i++) {
            wm[i + 2] = 0.0;
        }
        jac(S->n, S->tn, y, 0, 0, &wm[2], S->n, rpar, ipar);
        if (S->jsv == 1) {
            dcopy_(&lenp, &wm[2], &(int){1}, &wm[S->locjs - 1], &(int){1});
        }
    }

    if ((jok == -1) && (S->miter == 2)) {
        // If MITER = 2, make N calls to F to approximate the Jacobian.
        S->nje += 1;
        S->nslj = S->nst;
        S->jcur = 1;
        double fac = dvnorm(S->n, savf, ewt);
        double r0 = 1000.0 * fabs(S->h) * S->uround * (double)(S->n) * fac;
        if (r0 == 0.0) { r0 = 1.0; }
        double srur = wm[0];
        int j1 = 2;
        for (int j = 0; j < S->n; j++) {
            double yj = y[j];
            double r = fmax(srur * fabs(yj), r0 / ewt[j]);
            y[j] += r;
            fac = 1.0 / r;
            f(S->n, S->tn, y, ftem, rpar, ipar);
            for (int i = 0; i < S->n; i++) {
                wm[i + j1] = (ftem[i] - savf[i]) * fac;
            }
            y[j] = yj;
            j1 += S->n;
        }
        S->nfe += S->n;
        int lenp = S->n * S->n;
        if (S->jsv == 1) {
            dcopy_(&lenp, &wm[2], &(int){1}, &wm[S->locjs - 1], &(int){1});
        }
    }

    if ((jok == 1) && ((S->miter == 1) || (S->miter == 2))) {
        S->jcur = 0;
        int lenp = S->n * S->n;
        dcopy_(&lenp, &wm[S->locjs - 1], &(int){1}, &wm[2], &(int){1});
    }

    if ((S->miter == 1) || (S->miter == 2)) {
        // Multiply Jacobian by scalar, add identity, and do LU factorization.
        int lenp = S->n * S->n;
        dscal_(&lenp, &(double){-hrl1}, &wm[2], &(int){1});
        int j = 2;
        for (int i = 0; i < S->n; i++) {
            wm[j] += 1.0;
            j += S->n + 1;
        }
        S->nlu += 1;
        int ier = 0;
        dgetrf_(&S->n, &S->n, &wm[2], &S->n, &iwm[30], &ier);
        if (ier != 0) {
            *ierpj = 1;
        }
        return;
    }

    if (S->miter == 3) {
        // If MITER = 3, construct a diagonal approximation to J and P.
        S->nje += 1;
        S->jcur = 1;
        wm[1] = hrl1;
        double r = S->rl1 * 0.1;
        for (int i = 0; i < S->n; i++) {
            y[i] = y[i] + r * (S->h * savf[i] - yh[i + ldyh]);
        }
        f(S->n, S->tn, y, &wm[2], rpar, ipar);
        S->nfe += 1;
        for (int i = 0; i < S->n; i++) {
            double r0 = S->h * savf[i] - yh[i + ldyh];
            double di = 0.1 * r0 - S->h * (wm[i + 2] - savf[i]);
            wm[i + 2] = 1.0;
            if (fabs(r0) < S->uround / ewt[i]) {
                continue;
            }
            if (fabs(di) == 0.0) {
                *ierpj = 1;
                return;
            }
            wm[i + 2] = 0.1 * r0 / di;
        }
        return;
    }

    // Set constants for MITER = 4 or 5.
    int ml = iwm[0];
    int mu = iwm[1];
    int ml3 = ml + 3;
    int mband = ml + mu + 1;
    int meband = mband + ml;
    int lenp = meband * S->n;

    if ((jok == -1) && (S->miter == 4)) {
        // If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian.
        S->nje += 1;
        S->nslj = S->nst;
        S->jcur = 1;
        for (int i = 0; i < lenp; i++) {
            wm[i + 2] = 0.0;
        }
        jac(S->n, S->tn, y, ml, mu, &wm[ml3 - 1], meband, rpar, ipar);
        if (S->jsv == 1) {
            // Call DACOPY equivalent - copy banded matrix
            for (int ic = 0; ic < S->n; ic++) {
                dcopy_(&mband, &wm[ml3 - 1 + ic * meband], &(int){1},
                       &wm[S->locjs - 1 + ic * mband], &(int){1});
            }
        }
    }

    if ((jok == -1) && (S->miter == 5)) {
        // If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian.
        S->nje += 1;
        S->nslj = S->nst;
        S->jcur = 1;
        int mba = int_min(mband, S->n);
        int meb1 = meband - 1;
        double srur = wm[0];
        double fac = dvnorm(S->n, savf, ewt);
        double r0 = 1000.0 * fabs(S->h) * S->uround * (double)(S->n) * fac;
        if (r0 == 0.0) { r0 = 1.0; }
        for (int j = 0; j < mba; j++) {
            for (int i = j; i < S->n; i += mband) {
                double yi = y[i];
                double r = fmax(srur * fabs(yi), r0 / ewt[i]);
                y[i] = y[i] + r;
            }
            f(S->n, S->tn, y, ftem, rpar, ipar);
            for (int jj = j; jj < S->n; jj += mband) {
                y[jj] = yh[jj];
                double yjj = y[jj];
                double r = fmax(srur * fabs(yjj), r0 / ewt[jj]);
                fac = 1.0 / r;
                int i1 = int_max(jj - mu, 0);
                int i2 = int_min(jj + ml, S->n - 1);
                int ii = (jj + 1) * meb1 - ml + 2;
                for (int i = i1; i <= i2; i++) {
                    wm[ii + i] = (ftem[i] - savf[i]) * fac;
                }
            }
        }
        S->nfe += mba;
        if (S->jsv == 1) {
            // Call DACOPY equivalent
            for (int ic = 0; ic < S->n; ic++) {
                dcopy_(&mband, &wm[ml3 - 1 + ic * meband], &(int){1}, &wm[S->locjs - 1 + ic * mband], &(int){1});
            }
        }
    }

    if (jok == 1) {
        S->jcur = 0;
        // Call DACOPY equivalent
        for (int ic = 0; ic < S->n; ic++) {
            dcopy_(&mband, &wm[S->locjs - 1 + ic * mband], &(int){1}, &wm[ml3 - 1 + ic * meband], &(int){1});
        }
    }

    // Multiply Jacobian by scalar, add identity, and do LU decomposition.
    double con = -hrl1;
    dscal_(&lenp, &con, &wm[2], &(int){1});
    int ii = mband + 1;
    for (int i = 0; i < S->n; i++) {
        wm[ii] += 1.0;
        ii += meband;
    }
    S->nlu += 1;
    int ier = 0;
    dgbtrf_(&S->n, &S->n, &ml, &mu, &wm[2], &meband, &iwm[30], &ier);
    if (ier != 0) {
        *ierpj = 1;
    }

    return;
}


static void
dvnlsd(
    vode_common_struct_t* S,
    double* restrict y,
    double* restrict yh,
    const int ldyh,
    double* restrict savf,
    double* restrict ewt,
    double* acor,
    int* restrict iwm,
    double* restrict wm,
    vode_func_t f,
    vode_jac_t jac,
    int* nflag,
    double* restrict rpar,
    int* restrict ipar
)
{
    const int msbp = 20;
    const double crdown = 0.3, rdiv = 2.0;
    const int maxcor = 3;
    int m = 0, iersl = 0, ierpj = 0;
    double del, delp = 0.0;

    if (S->jstart == 0) { S->nslp = 0; }
    if (*nflag == 0) { S->icf = 0; }
    if (*nflag == -2) { S->ipup = S->miter; }
    if ((S->jstart == 0) || (S->jstart == -1)) { S->ipup = S->miter; }
    // If this is functional iteration, set CRATE == 1
    if (S->miter != 0) {
        // RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
        // When RC differs from 1 by more than CCMAX, IPUP is set to MITER
        // to force DVJAC to be called, if a Jacobian is involved.
        // In any case, DVJAC is called at least every MSBP steps.
        if ((fabs(S->rc - 1.0) > 0.3) || (S->nst >= S->nslp + msbp)) { S->ipup = S->miter; }
    } else {
        S->crate = 1.0;
    }

    // Up to MAXCOR corrector iterations are taken.  A convergence test is
    // made on the r.m.s. norm of each correction, weighted by the error
    // weight vector EWT.  The sum of the corrections is accumulated in the
    // vector ACOR(i). The YH array is not altered in the corrector loop.

    // Corrector iteration section has two nested while loops:
    // - Outer loop implements Fortran label 220 (retry with new Jacobian if needed)
    // - Inner loop implements Fortran label 270 (up to MAXCOR corrector iterations)
    //
    // Inner loop has only TWO exit paths:
    //   1. SUCCESS (GO TO 450): return directly from function
    //   2. FAILURE (GO TO 410): break to label 410 handler
    // Therefore, if we exit the inner loop via break, we ALWAYS go to label 410.


    // Outer retry loop - Fortran label 220
    while (1) {
        m = 0;
        delp = 0.0;
        dcopy_(&S->n, yh, &(int){1}, y, &(int){1});
        f(S->n, S->tn, y, savf, rpar, ipar);
        S->nfe += 1;

        if (S->ipup > 0) {
            // If indicated, the matrix P = I - h*rl1*J is reevaluated and
            // preprocessed before starting the corrector iteration.  IPUP is set
            // to 0 as an indicator that this has been done.
            dvjac(y, yh, ldyh, ewt, acor, savf, wm, iwm, f, jac, &ierpj, rpar, ipar, S);
            S->ipup = 0;
            S->rc = 1.0;
            S->drc = 0.0;
            S->crate = 1.0;
            S->nslp = S->nst;
            // If matrix is singular, give up immediately (Fortran line 2896: GO TO 430)
            if (ierpj != 0) {
                *nflag = -1;
                S->icf = 1;
                S->ipup = S->miter;
                return;
            }
        }

        // Fortran label 250: Zero ACOR
        for (int i = 0; i < S->n; i++) {
            acor[i] = 0.0;
        }

        // Inner corrector iteration loop - Fortran label 270
        while (1) {

            if (S->miter == 0) {
                // Functional iteration: update Y directly from last function evaluation
                // Fortran lines 2901-2909
                for (int i = 0; i < S->n; i++) {
                    savf[i] = S->rl1 * (S->h * savf[i] - yh[i + ldyh]);
                }
                for (int i = 0; i < S->n; i++) {
                    y[i] = savf[i] - acor[i];
                }
                del = dvnorm(S->n, y, ewt);
                for (int i = 0; i < S->n; i++) {
                    y[i] = yh[i] + savf[i];
                }
                dcopy_(&S->n, savf, &(int){1}, acor, &(int){1});
            } else {
                // Chord iteration: solve linear system
                // Fortran lines 2917-2928
                for (int i = 0; i < S->n; i++) {
                    y[i] = (S->rl1 * S->h) * savf[i] - (S->rl1 * yh[i + ldyh] + acor[i]);
                }
                dvsol(S, wm, iwm, y, &iersl);
                S->nni += 1;
                if (iersl > 0) {
                    break;  // Exit to label 410 handler
                }
                if ((S->meth == 2) && (S->rc != 1.0)) {
                    double cscale = 2.0 / (1.0 + S->rc);
                    dscal_(&S->n, &cscale, y, &(int){1});
                }
                del = dvnorm(S->n, y, ewt);
                daxpy_(&S->n, &(double){1.0}, y, &(int){1}, acor, &(int){1});
                for (int i = 0; i < S->n; i++) {
                    y[i] = yh[i] + acor[i];
                }
            }

            // Fortran label 400: Test for convergence
            if (m != 0) {
                S->crate = fmax(crdown * S->crate, del / delp);
            }
            double dcon = del * fmin(1.0, S->crate) / S->tq[3];
            // Fortran line 2935: IF (DCON .LE. ONE) GO TO 450
            if (dcon <= 1.0) {
                // Fortran label 450: Successful convergence - RETURN
                *nflag = 0;
                S->jcur = 0;
                S->icf = 0;
                if (m == 0) {
                    S->acnrm = del;
                } else {
                    S->acnrm = dvnorm(S->n, acor, ewt);
                }
                return;
            }

            // Not converged yet. Increment iteration counter.
            m = m + 1;

            if (m == maxcor) { break; }
            if ((m >= 2) && (del > rdiv * delp)) { break; }

            // Continue corrector iteration; compute new f and GO TO 270
            delp = del;
            f(S->n, S->tn, y, savf, rpar, ipar);
            S->nfe += 1;
        }
        // End inner while loop

        // 410 Convergence failure handler
        // We ALWAYS reach here after breaking from inner loop (no success path leads here)
        // IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
        if ((S->miter == 0) || (S->jcur == 1)) {
            // Fortran label 430: Give up
            *nflag = -1;
            S->icf = 2;
            S->ipup = S->miter;
            return;
        }

        // Retry with new Jacobian. Protection against infinite retry:
        // After retry, JCUR will be 1 (we just updated Jacobian), so if we
        // fail again, the condition above will be true and we'll give up.
        S->icf = 1;
        S->ipup = S->miter;
    }
    // End outer while loop

}


/**
 * @brief Perform one step of the integration of an ODE system.
 *
 * DVSTEP performs one step of the integration of an initial value problem for a
 * system of ordinary differential equations. It calls the nonlinear system solver
 * VNLS (which in turn calls DVNLSD) for the solution of the nonlinear system
 * arising in the time step. Thus it is independent of the problem Jacobian
 * structure and the type of nonlinear system solution method.
 *
 * DVSTEP returns a completion flag KFLAG (in the state struct). A return with
 * KFLAG = -1 or -2 means either |H| = HMIN or 10 consecutive failures occurred.
 * On a return with KFLAG negative, the values of TN and the YH array are as of
 * the beginning of the last step, and H is the last step size attempted.
 *
 * Algorithm Overview:
 * -------------------
 * 1. Handle initialization (JSTART=0), normal continuation (JSTART>0), or
 *    restart (JSTART=-1)
 * 2. Adjust order if needed (via DVJUST) and rescale history array YH
 * 3. Compute predictor by multiplying YH by Pascal triangle matrix
 * 4. Call DVSET to compute integration coefficients
 * 5. Call VNLS to solve the nonlinear system (corrector phase)
 * 6. If convergence fails, retract YH, reduce H, and retry
 * 7. Test local error (DSM = ACNRM/TQ(2))
 * 8. If error test fails, retract YH and retry with smaller H/order
 * 9. If successful, update YH, TAU arrays
 * 10. Select new order and step size for next step
 *
 * Error Handling:
 * ---------------
 * - Convergence failure: Reduce H by factor ETACF (0.25), retry
 * - Error test failure: Reduce H based on error estimate, retry
 * - Multiple failures (≥3): Reduce order by 1
 * - Severe failures (≥7): Return with KFLAG = -1
 * - |H| ≤ HMIN: Return with KFLAG = -2
 *
 * Order Selection:
 * ----------------
 * After a successful step (if ETAMAX ≠ 1), compute step size ratios for:
 * - ETAQ: same order (NQ)
 * - ETAQM1: order decrease (NQ-1)
 * - ETAQP1: order increase (NQ+1)
 * Choose the order giving the largest step size ratio, subject to:
 * - Must increase H by at least factor THRESH (1.5)
 * - Constrained by ETAMAX and HMXI
 *
 * @param S      Pointer to VODE state structure
 * @param y      Solution vector (length N)
 * @param yh     History array (LDYH × LMAX), column-major
 * @param ldyh   Leading dimension of YH (≥ N)
 * @param yh1    1D view of YH for flat indexing
 * @param ewt    Error weight vector (length N)
 * @param savf   Work array (length N), also f(t,y) storage
 * @param vsav   Work array for VNLS (length N)
 * @param acor   Accumulated corrections (length N), returns local error estimate
 * @param wm     Real work array for matrix operations
 * @param iwm    Integer work array for matrix operations
 * @param f      User function for dy/dt = f(t,y)
 * @param jac    User Jacobian function (if applicable)
 * @param rpar   User real parameter array
 * @param ipar   User integer parameter array
 *
 * State Variables Used:
 * ---------------------
 * Input/Output: KFLAG, JCUR, TN, H, HSCAL, NQ, L, NQWAIT, NEWQ, NEWH, ETA,
 *               ETAMAX, RC, PRL1, TAU, ACNRM
 * Output only: NST, HU, NQU, NCFN, NETF, NFE, JSTART
 * Input only: N, MAXORD, HMIN, HMXI, KUTH, EL, TQ, UROUND
 */
static void
dvstep(
    vode_common_struct_t* S,
    double* restrict y,
    double* yh,
    const int ldyh,
    double* yh1,
    double* restrict ewt,
    double* savf,
    double* restrict acor,
    double* restrict wm,
    int* restrict iwm,
    vode_func_t f,
    vode_jac_t jac,
    double* restrict rpar,
    int* restrict ipar
)
{
    const int kfc = -3, kfh = -7, mxncf = 10;
    const double addon = 1.0e-6, bias1 = 6.0, bias2 = 6.0, bias3 = 10.0, etacf = 0.25, etamin = 0.1, etamxf = 0.2,
                 etamx1 = 1.0e4, etamx2 = 10.0, etamx3 = 10.0, onepsm = 1.00001, thresh = 1.5;
    int int1 = 1;

    // Local variables
    double told, r, dsm, flotl, ddn, dup, cnquot;
    double etaqp1;  // Step size ratio for order+1 (computed locally each time)
    int ncf, nflag, i, i1, i2, j, jb, iback;
    int need_rescale;  // Flag for whether rescaling is needed

    // Initialize for this step
    S->kflag = 0;
    told = S->tn;
    ncf = 0;
    S->jcur = 0;
    nflag = 0;

    // Three-way branch based on JSTART
    if (S->jstart > 0) {
        // ===================================================================
        // Label 20: Normal continuation step (JSTART > 0)
        // Take preliminary actions on a normal continuation step.
        // If the driver changed H, then ETA must be reset and NEWH set to 1.
        // If a change of order was dictated on the previous step, then
        // it is done here and appropriate adjustments in the history are made.
        // On an order decrease, the history array is adjusted by DVJUST.
        // On an order increase, the history array is augmented by a column.
        // On a change of step size H, the history array YH is rescaled.
        // ===================================================================

        if (S->kuth == 1) {
            S->eta = fmin(S->eta, S->h / S->hscal);
            S->newh = 1;
        }

        // Label 50: Handle order changes if NEWH != 0
        if (S->newh != 0) {
            if (S->newq == S->nq) {
                // No order change, just rescaling needed
                need_rescale = 1;
            } else if (S->newq < S->nq) {
                // Order decrease
                dvjust(yh, ldyh, -1, S);
                S->nq = S->newq;
                S->l = S->nq + 1;
                S->nqwait = S->l;
                need_rescale = 1;
            } else {  // S->newq > S->nq
                // Order increase
                dvjust(yh, ldyh, 1, S);
                S->nq = S->newq;
                S->l = S->nq + 1;
                S->nqwait = S->l;
                need_rescale = 1;
            }
        } else {
            // NEWH == 0, skip rescaling
            need_rescale = 0;
        }

    } else if (S->jstart == -1) {
        // ===================================================================
        // Label 100: Restart (JSTART == -1)
        // This block handles preliminaries needed when JSTART = -1.
        // If N was reduced, zero out part of YH to avoid undefined references.
        // If MAXORD was reduced to a value less than the tentative order NEWQ,
        // then NQ is set to MAXORD, and a new H ratio ETA is chosen.
        // Otherwise, we take the same preliminary actions as for JSTART > 0.
        // In any case, NQWAIT is reset to L = NQ + 1 to prevent further
        // changes in order for that many steps.
        // The new H ratio ETA is limited by the input H if KUTH = 1,
        // by HMIN if KUTH = 0, and by HMXI in any case.
        // Finally, the history array YH is rescaled.
        // ===================================================================

        S->lmax = S->maxord + 1;

        // Zero out part of YH if N was reduced
        if (S->n != ldyh) {
            i1 = (S->newq + 1) * ldyh;
            i2 = (S->maxord + 1) * ldyh;
            if (i1 < i2) {
                for (i = i1; i < i2; i++) {
                    yh1[i] = 0.0;
                }
            }
        }

        // Label 120: Check if NEWQ <= MAXORD
        if (S->newq > S->maxord) {
            // MAXORD was reduced to a value less than the tentative order NEWQ
            flotl = (double)S->lmax;

            if (S->maxord < S->nq - 1) {
                ddn = dvnorm(S->n, savf, ewt) / S->tq[0];
                S->eta = 1.0 / (pow(bias1 * ddn, 1.0 / flotl) + addon);
            }

            if (S->maxord == S->nq && S->newq == S->nq + 1) {
                S->eta = S->etaq;
            }

            if (S->maxord == S->nq - 1 && S->newq == S->nq + 1) {
                S->eta = S->etaqm1;
                dvjust(yh, ldyh, -1, S);
            }

            if (S->maxord == S->nq - 1 && S->newq == S->nq) {
                ddn = dvnorm(S->n, savf, ewt) / S->tq[0];
                S->eta = 1.0 / (pow(bias1 * ddn, 1.0 / flotl) + addon);
                dvjust(yh, ldyh, -1, S);
            }

            S->eta = fmin(S->eta, 1.0);
            S->nq = S->maxord;
            S->l = S->lmax;
        }

        // Label 140: Limit ETA by H constraints
        if (S->kuth == 1) {
            S->eta = fmin(S->eta, fabs(S->h / S->hscal));
        }
        if (S->kuth == 0) {
            S->eta = fmax(S->eta, S->hmin / fabs(S->hscal));
        }
        S->eta = S->eta / fmax(1.0, fabs(S->hscal) * S->hmxi * S->eta);
        S->newh = 1;
        S->nqwait = S->l;

        // If NEWQ <= MAXORD, handle order change (same as Label 50 logic)
        if (S->newq <= S->maxord) {
            if (S->newq == S->nq) {
                need_rescale = 1;
            } else if (S->newq < S->nq) {
                dvjust(yh, ldyh, -1, S);
                S->nq = S->newq;
                S->l = S->nq + 1;
                S->nqwait = S->l;
                need_rescale = 1;
            } else {  // S->newq > S->nq
                dvjust(yh, ldyh, 1, S);
                S->nq = S->newq;
                S->l = S->nq + 1;
                S->nqwait = S->l;
                need_rescale = 1;
            }
        } else {
            // NEWQ > MAXORD, fall through to rescaling
            need_rescale = 1;
        }

    } else {
        // ===================================================================
        // JSTART == 0: First call
        // On the first call, the order is set to 1, and other variables are
        // initialized. ETAMAX is the maximum ratio by which H can be increased
        // in a single step. It is normally 10, but is larger during the
        // first step to compensate for the small initial H. If a failure
        // occurs (in corrector convergence or error test), ETAMAX is set to 1
        // for the next increase.
        // ===================================================================

        S->lmax = S->maxord + 1;
        S->nq = 1;
        S->l = 2;
        S->nqnyh = S->nq * ldyh;
        S->tau[0] = S->h;
        S->prl1 = 1.0;
        S->rc = 0.0;
        S->etamax = etamx1;
        S->nqwait = 2;
        S->hscal = S->h;

        // Skip rescaling and go directly to predictor (Label 200)
        need_rescale = 0;
    }

    // ===================================================================
    // MAIN RETRY LOOP
    // Combines Label 150 (rescaling) and Label 200 (predictor) with retry logic
    // ===================================================================

    while (1) {
        // ===============================================================
        // Label 150: Rescale the history array for a change in H by ETA
        // ===============================================================

        if (need_rescale) {
            r = 1.0;
            for (j = 1; j < S->l; j++) {
                r *= S->eta;
                dscal_(&S->n, &r, &yh[j * ldyh], &int1);
            }
            S->h = S->hscal * S->eta;
            S->hscal = S->h;
            S->rc *= S->eta;
            S->nqnyh = S->nq * ldyh;
            need_rescale = 0;  // Clear flag after rescaling
        }

        // ===============================================================
        // Label 200: Compute predicted values
        // This section computes the predicted values by effectively
        // multiplying the YH array by the Pascal triangle matrix.
        // DVSET is called to calculate all integration coefficients.
        // RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
        // ===============================================================

        S->tn += S->h;
        i1 = S->nqnyh;
        for (jb = 0; jb < S->nq; jb++) {
            i1 -= ldyh;
            for (i = i1; i < S->nqnyh; i++) {
                yh1[i] += yh1[i + ldyh];
            }
        }

        dvset(S);
        S->rl1 = 1.0 / S->el[1];
        S->rc = S->rc * (S->rl1 / S->prl1);
        S->prl1 = S->rl1;

        // ===============================================================
        // Call the nonlinear system solver
        // ===============================================================

        dvnlsd(S, y, yh, ldyh, savf, ewt, acor, iwm, wm, f, jac, &nflag, rpar, ipar);

        if (nflag != 0) {
            // ===========================================================
            // Convergence failure: DVNLSD routine failed to achieve convergence
            // The YH array is retracted to its values before prediction.
            // The step size H is reduced and the step is retried, if possible.
            // Otherwise, an error exit is taken.
            // ===========================================================

            ncf++;
            S->ncfn++;
            S->etamax = 1.0;
            S->tn = told;

            // Retract YH array
            i1 = S->nqnyh;
            for (jb = 0; jb < S->nq; jb++) {
                i1 -= ldyh;
                for (i = i1; i < S->nqnyh; i++) {
                    yh1[i] -= yh1[i + ldyh];
                }
            }

            // Check for fatal errors
            if (nflag < -1) {
                // Label 680: Set KFLAG based on NFLAG
                if (nflag == -2) S->kflag = -3;
                if (nflag == -3) S->kflag = -4;
                break;  // Exit to Label 720
            }

            if (fabs(S->h) <= S->hmin * onepsm) {
                // Label 670: H too small
                S->kflag = -2;
                break;  // Exit to Label 720
            }

            if (ncf == mxncf) {
                // Label 670: Too many convergence failures
                S->kflag = -2;
                break;  // Exit to Label 720
            }

            // Retry with reduced H
            S->eta = etacf;
            S->eta = fmax(S->eta, S->hmin / fabs(S->h));
            nflag = -1;
            need_rescale = 1;
            continue;  // GO TO 150 (retry with rescaling)
        }

        // ===============================================================
        // Label 450: The corrector has converged (NFLAG = 0)
        // The local error test is made
        // ===============================================================

        dsm = S->acnrm / S->tq[1];

        if (dsm > 1.0) {
            // ===========================================================
            // Label 500: Error test failed
            // KFLAG keeps track of multiple failures.
            // Restore TN and the YH array to their previous values, and prepare
            // to try the step again. Compute the optimum step size for the
            // same order. After repeated failures, H is forced to decrease
            // more rapidly.
            // ===========================================================

            S->kflag--;
            S->netf++;
            nflag = -2;
            S->tn = told;

            // Retract YH array
            i1 = S->nqnyh;
            for (jb = 0; jb < S->nq; jb++) {
                i1 -= ldyh;
                for (i = i1; i < S->nqnyh; i++) {
                    yh1[i] -= yh1[i + ldyh];
                }
            }

            if (fabs(S->h) <= S->hmin * onepsm) {
                // Label 660: H too small
                S->kflag = -1;
                break;  // Exit to Label 720
            }

            S->etamax = 1.0;

            if (S->kflag > kfc) {
                // Regular error failure: compute new H at current order
                flotl = (double)S->l;
                S->eta = 1.0 / (pow(bias2 * dsm, 1.0 / flotl) + addon);
                S->eta = fmax(S->eta, S->hmin / fabs(S->h));
                S->eta = fmax(S->eta, etamin);
                if (S->kflag <= -2 && S->eta > etamxf) {
                    S->eta = etamxf;
                }
                need_rescale = 1;
                continue;  // GO TO 150 (retry with rescaling)
            }

            // ===========================================================
            // Label 530: Control reaches this section if 3 or more consecutive
            // failures have occurred. It is assumed that the elements of the YH
            // array have accumulated errors of the wrong order. The order is
            // reduced by one, if possible. Then H is reduced by a factor of 0.1
            // and the step is retried. After a total of 7 consecutive failures,
            // an exit is taken with KFLAG = -1.
            // ===========================================================

            if (S->kflag == kfh) {
                // Label 660: Too many failures
                S->kflag = -1;
                break;  // Exit to Label 720
            }

            if (S->nq == 1) {
                // Label 540: Special case for NQ=1
                // Cannot reduce order, so reduce H and recompute YH(,2)
                S->eta = fmax(etamin, S->hmin / fabs(S->h));
                S->h = S->h * S->eta;
                S->hscal = S->h;
                S->tau[0] = S->h;

                f(S->n, S->tn, y, savf, rpar, ipar);
                S->nfe++;

                for (i = 0; i < S->n; i++) {
                    yh[i + ldyh] = S->h * savf[i];
                }

                S->nqwait = 10;
                need_rescale = 0;  // Already updated H directly
                continue;  // GO TO 200 (retry without rescaling)
            }

            // Reduce order by one
            S->eta = fmax(etamin, S->hmin / fabs(S->h));
            dvjust(yh, ldyh, -1, S);
            S->l = S->nq;
            S->nq--;
            S->nqwait = S->l;
            need_rescale = 1;
            continue;  // GO TO 150 (retry with rescaling)
        }

        // ===============================================================
        // Successful step: Update the YH and TAU arrays and decrement NQWAIT
        // If NQWAIT is then 1 and NQ < MAXORD, then ACOR is saved
        // for use in a possible order increase on the next step.
        // If ETAMAX = 1 (a failure occurred this step), keep NQWAIT >= 2.
        // ===============================================================

        S->kflag = 0;
        S->nst++;
        S->hu = S->h;
        S->nqu = S->nq;

        // Update TAU array
        for (iback = 1; iback <= S->nq; iback++) {
            i = S->l - iback;
            S->tau[i] = S->tau[i - 1];
        }
        S->tau[0] = S->h;

        // Update YH array
        for (j = 0; j < S->l; j++) {
            daxpy_(&S->n, &S->el[j], acor, &int1, &yh[j * ldyh], &int1);
        }

        S->nqwait--;

        if ((S->l != S->lmax) && (S->nqwait == 1)) {
            dcopy_(&S->n, acor, &int1, &yh[(S->lmax - 1) * ldyh], &int1);
            S->conp = S->tq[4];
        }

        // Label 490: Check if a failure occurred this step
        if (S->etamax != 1.0) {
            // No failure - proceed to order selection (Label 560)

            // ===============================================================
            // Label 560: Order selection
            // If NQWAIT = 0, an increase or decrease in order by one is considered.
            // Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
            // be multiplied at order q, q-1, or q+1, respectively.
            // The largest of these is determined, and the new order and
            // step size set accordingly.
            // A change of H or NQ is made only if H increases by at least a
            // factor of THRESH. If an order change is considered and rejected,
            // then NQWAIT is set to 2 (reconsider it after 2 steps).
            // ===============================================================

            flotl = (double)S->l;
            S->etaq = 1.0 / (pow(bias2 * dsm, 1.0 / flotl) + addon);

            if (S->nqwait != 0) {
                // Label 600: Use current order
                S->eta = S->etaq;
                S->newq = S->nq;
            } else {
                // Consider order change
                S->nqwait = 2;
                S->etaqm1 = 0.0;

                if (S->nq > 1) {
                    // Label 570: Compute ETAQM1 for order decrease
                    ddn = dvnorm(S->n, &yh[(S->l - 1) * ldyh], ewt) / S->tq[0];
                    S->etaqm1 = 1.0 / (pow(bias1 * ddn, 1.0 / (flotl - 1.0)) + addon);
                }

                etaqp1 = 0.0;
                if (S->l != S->lmax) {
                    // Label 575: Compute ETAQP1 for order increase
                    cnquot = (S->tq[4] / S->conp) * pow(S->h / S->tau[1], S->l);
                    for (i = 0; i < S->n; i++) {
                        savf[i] = acor[i] - cnquot * yh[i + (S->lmax - 1) * ldyh];
                    }
                    dup = dvnorm(S->n, savf, ewt) / S->tq[2];
                    etaqp1 = 1.0 / (pow(bias3 * dup, 1.0 / (flotl + 1.0)) + addon);
                }

                // Label 580-590: Determine which ETA is largest
                if (S->etaq >= etaqp1) {
                    if (S->etaq < S->etaqm1) {
                        // Label 610: Choose order decrease
                        S->eta = S->etaqm1;
                        S->newq = S->nq - 1;
                    } else {
                        // Label 600: Choose current order
                        S->eta = S->etaq;
                        S->newq = S->nq;
                    }
                } else {
                    if (etaqp1 > S->etaqm1) {
                        // Label 620: Choose order increase
                        S->eta = etaqp1;
                        S->newq = S->nq + 1;
                        dcopy_(&S->n, acor, &int1, &yh[(S->lmax - 1) * ldyh], &int1);
                    } else {
                        // Label 610: Choose order decrease
                        S->eta = S->etaqm1;
                        S->newq = S->nq - 1;
                    }
                }
            }

            // Label 630: Test tentative new H against THRESH, ETAMAX, and HMXI
            if (S->eta < thresh || S->etamax == 1.0) {
                // Label 640: No change in H or order
                S->newq = S->nq;
                S->newh = 0;
                S->eta = 1.0;
                S->hnew = S->h;
            } else {
                // Accept change
                S->eta = fmin(S->eta, S->etamax);
                S->eta = S->eta / fmax(1.0, fabs(S->h) * S->hmxi * S->eta);
                S->newh = 1;
                S->hnew = S->h * S->eta;
            }

        } else {
            // Label 490: Failure occurred this step (ETAMAX == 1)
            // Keep NQWAIT >= 2 and don't change order
            if (S->nqwait < 2) {
                S->nqwait = 2;
            }
            S->newq = S->nq;
            S->newh = 0;
            S->eta = 1.0;
            S->hnew = S->h;
        }

        // Label 690: Set ETAMAX for next step (both success paths converge here)
        S->etamax = etamx3;
        if (S->nst <= 10) {
            S->etamax = etamx2;
        }

        // Label 700: Scale ACOR
        r = 1.0 / S->tq[1];
        dscal_(&S->n, &r, acor, &int1);

        // Exit the retry loop - step successful
        break;

    }

    // Label 720: Final cleanup and return
    S->jstart = 1;

    return;
}


/**
 * @brief Sets optional output values in RWORK and IWORK arrays.
 *
 * This helper function sets the standard optional output values that are
 * common to both successful and error returns from DVODE. The caller may
 * need to override specific values after calling this function:
 *
 * For successful returns (ISTATE=2):
 *   - Call this function, then set RWORK[11] = S->hnew and IWORK[14] = S->newq
 *
 * For error returns (ISTATE < 0):
 *   - Call this function, then set RWORK[11] = S->h and IWORK[14] = S->nq
 *
 * Standard outputs set by this function:
 *   RWORK[10] = HU    (step size successfully used on last step)
 *   RWORK[11] = HNEW  (step size to be attempted on next step) - may be overridden
 *   RWORK[12] = TN    (current value of independent variable)
 *   IWORK[10] = NST   (number of steps taken)
 *   IWORK[11] = NFE   (number of f evaluations)
 *   IWORK[12] = NJE   (number of Jacobian evaluations)
 *   IWORK[13] = NQU   (method order on last step)
 *   IWORK[14] = NEWQ  (order to be attempted on next step) - may be overridden
 *   IWORK[18] = NLU   (number of LU decompositions)
 *   IWORK[19] = NNI   (number of nonlinear iterations)
 *   IWORK[20] = NCFN  (number of convergence failures)
 *   IWORK[21] = NETF  (number of error test failures)
 *
 * @param S      Pointer to VODE common struct
 * @param rwork  Real work array
 * @param iwork  Integer work array
 */
static void
dvode_set_optional_output(vode_common_struct_t* S, double* rwork, int* iwork)
{
    rwork[10] = S->hu;
    rwork[11] = S->hnew;
    rwork[12] = S->tn;
    iwork[10] = S->nst;
    iwork[11] = S->nfe;
    iwork[12] = S->nje;
    iwork[13] = S->nqu;
    iwork[14] = S->newq;
    iwork[18] = S->nlu;
    iwork[19] = S->nni;
    iwork[20] = S->ncfn;
    iwork[21] = S->netf;
}


/**
 * @brief Main DVODE driver routine for ODE integration.
 *
 * DVODE solves the initial value problem for stiff or nonstiff systems
 * of first order ODEs: dy/dt = f(t,y)
 *
 * @param S       Pointer to VODE common struct (maintains state between calls)
 * @param f       User function computing dy/dt
 * @param neq     Number of first-order ODEs
 * @param y       Solution vector (input: initial values, output: computed solution)
 * @param t       Independent variable (input: initial t, output: current t)
 * @param tout    Next desired output point
 * @param itol    Tolerance type (1-4)
 * @param rtol    Relative tolerance
 * @param atol    Absolute tolerance
 * @param itask   Task indicator (1-5)
 * @param istate  State flag (input/output)
 * @param iopt    Optional input flag
 * @param rwork   Real work array
 * @param lrw     Length of rwork
 * @param iwork   Integer work array
 * @param liw     Length of iwork
 * @param jac     User Jacobian function (if applicable)
 * @param mf      Method flag
 * @param rpar    User real array
 * @param ipar    User integer array
 */
void
dvode(
    vode_common_struct_t* S,
    vode_func_t f,
    int neq,
    double* restrict y,
    double* t,
    double* tout,
    int itol,
    double* rtol,
    double* atol,
    int* itask,
    int* istate,
    int* iopt,
    double* restrict rwork,
    int lrw,
    int* restrict iwork,
    int liw,
    vode_jac_t jac,
    int mf,
    double* rpar,
    int* ipar
)
{
    // Constants
    static const int mord[2] = {12, 5};
    static const int mxstp0 = 500;
    static const int mxhnl0 = 10;
    const double zero = 0.0, one = 1.0, two = 2.0, four = 4.0;
    const double pt2 = 0.2, hun = 100.0;

    // Local variables
    double h0, hmax, hmx, tolsf, tnext, tp, atoli, rtoli, tcrit = 0.0;
    int i, iflag, jco, lenwm, lenrw, leniw, lf0, mband, mfa, ml = 0, mu = 0;
    int maxord_input, mxstep_input, mxhnil_input;
    int nslast, kgo;
    int ihit;

    // BLAS/LAPACK interface
    int int1 = 1;

    // ===================================================================
    // Block A: Initial validation on every call
    // ===================================================================
    if (*istate < 1 || *istate > 3) {
        *istate = -3;
        return;
    }
    if (*itask < 1 || *itask > 5) {
        *istate = -3;
        return;
    }

    // Check initialization for continuation calls
    if (*istate > 1) {
        if (S->init != 1) {
            *istate = -3;
            return;
        }
    }

    // Handle first call special case
    if (*istate == 1) {
        S->init = 0;
        if (*tout == *t) {
            return;
        }
    }

    // ===================================================================
    // Block B: Input validation and workspace setup (ISTATE=1 or 3)
    // ===================================================================
    if (*istate == 1 || *istate == 3) {
        // Check NEQ
        if (neq <= 0) {
            *istate = -3;
            return;
        }
        if (*istate == 3 && neq > S->n) {
            *istate = -3;
            return;
        }
        S->n = neq;

        // Check ITOL and IOPT
        if (itol < 1 || itol > 4) {
            *istate = -3;
            return;
        }
        if (*iopt < 0 || *iopt > 1) {
            *istate = -3;
            return;
        }

        // Decode MF
        S->jsv = (mf > 0) ? 1 : -1;
        mfa = abs(mf);
        S->meth = mfa / 10;
        S->miter = mfa - 10 * S->meth;

        if (S->meth < 1 || S->meth > 2) {
            *istate = -3;
            return;
        }
        if (S->miter < 0 || S->miter > 5) {
            *istate = -3;
            return;
        }

        // Check ML and MU for banded Jacobian
        if (S->miter > 3) {
            ml = iwork[0];
            mu = iwork[1];
            if (ml < 0 || ml >= neq) {
                *istate = -3;
                return;
            }
            if (mu < 0 || mu >= neq) {
                *istate = -3;
                return;
            }
        }

        // Process optional inputs
        if (*iopt == 0) {
            S->maxord = mord[S->meth - 1];
            S->mxstep = mxstp0;
            S->mxhnil = mxhnl0;
            if (*istate == 1) {
                h0 = zero;
            }
            S->hmxi = zero;
            S->hmin = zero;
        } else {
            // IOPT == 1: read optional inputs
            maxord_input = iwork[4];
            if (maxord_input < 0) {
                *istate = -3;
                return;
            }
            S->maxord = (maxord_input == 0) ? 100 : int_min(maxord_input, mord[S->meth - 1]);

            mxstep_input = iwork[5];
            if (mxstep_input < 0) {
                *istate = -3;
                return;
            }
            S->mxstep = (mxstep_input == 0) ? mxstp0 : mxstep_input;

            mxhnil_input = iwork[6];
            if (mxhnil_input < 0) {
                *istate = -3;
                return;
            }
            S->mxhnil = (mxhnil_input == 0) ? mxhnl0 : mxhnil_input;

            if (*istate == 1) {
                h0 = rwork[4];
                if ((*tout - *t) * h0 < zero) {
                    *istate = -3;
                    return;
                }
            }

            hmax = rwork[5];
            if (hmax < zero) {
                *istate = -3;
                return;
            }
            S->hmxi = zero;
            if (hmax > zero) {
                S->hmxi = one / hmax;
            }

            S->hmin = rwork[6];
            if (S->hmin < zero) {
                *istate = -3;
                return;
            }
        }

        // Set work array pointers
        S->lyh = 20;
        if (*istate == 1) {
            S->nyh = S->n;
        }
        S->lwm = S->lyh + (S->maxord + 1) * S->nyh;

        jco = (S->jsv > 0) ? 1 : 0;
        if (S->miter == 0) {
            lenwm = 0;
        } else if (S->miter == 1 || S->miter == 2) {
            lenwm = 2 + (1 + jco) * S->n * S->n;
            S->locjs = S->n * S->n + 2;
        } else if (S->miter == 3) {
            lenwm = 2 + S->n;
        } else {  // MITER == 4 or 5
            mband = ml + mu + 1;
            int lenp = (mband + ml) * S->n;
            int lenj = mband * S->n;
            lenwm = 2 + lenp + jco * lenj;
            S->locjs = lenp + 2;
        }

        S->lewt = S->lwm + lenwm;
        S->lsavf = S->lewt + S->n;
        S->lacor = S->lsavf + S->n;
        lenrw = S->lacor + S->n - 1;
        iwork[16] = lenrw;

        S->liwm = 0;
        leniw = (S->miter == 0 || S->miter == 3) ? 30 : (30 + S->n);
        iwork[17] = leniw;

        if (lenrw > lrw) {
            *istate = -3;
            return;
        }
        if (leniw > liw) {
            *istate = -3;
            return;
        }

        // Check RTOL and ATOL for legality
        rtoli = rtol[0];
        atoli = atol[0];
        for (i = 0; i < S->n; i++) {
            if (itol >= 3) {
                rtoli = rtol[i];
            }
            if (itol == 2 || itol == 4) {
                atoli = atol[i];
            }
            if (rtoli < zero) {
                *istate = -3;
                return;
            }
            if (atoli < zero) {
                *istate = -3;
                return;
            }
        }

        // Handle ISTATE=3 parameter changes
        if (*istate == 3) {
            S->jstart = -1;
            if (S->nq > S->maxord) {
                // MAXORD was reduced below NQ. Copy YH(*,MAXORD+2) into SAVF.
                dcopy_(&S->n, &rwork[S->lwm], &int1, &rwork[S->lsavf], &int1);
            }
            // Reload WM(1) = RWORK(LWM), since LWM may have changed.
            if (S->miter > 0) {
                rwork[S->lwm] = sqrt(S->uround);
            }
        }
    }

    // ===================================================================
    // Block C: Initialization for first call only (ISTATE=1)
    // ===================================================================
    if (*istate == 1) {
        S->uround = DBL_EPSILON;
        S->tn = *t;

        if (*itask == 4 || *itask == 5) {
            tcrit = rwork[0];
            if ((tcrit - *tout) * (*tout - *t) < zero) {
                *istate = -3;
                return;
            }
            if (h0 != zero && (*t + h0 - tcrit) * h0 > zero) {
                h0 = tcrit - *t;
            }
        }

        S->jstart = 0;
        if (S->miter > 0) {
            rwork[S->lwm] = sqrt(S->uround);
        }
        S->ccmxj = pt2;
        S->msbj = 50;
        S->nhnil = 0;
        S->nst = 0;
        S->nje = 0;
        S->nni = 0;
        S->ncfn = 0;
        S->netf = 0;
        S->nlu = 0;
        S->nslj = 0;
        nslast = 0;
        S->hu = zero;
        S->nqu = 0;

        // Initial call to F. (LF0 points to YH(*,2).)
        lf0 = S->lyh + S->nyh;
        f(S->n, *t, y, &rwork[lf0], rpar, ipar);
        S->nfe = 1;

        // Load the initial value vector in YH.
        dcopy_(&S->n, y, &int1, &rwork[S->lyh], &int1);

        // Load and invert the EWT array. (H is temporarily set to 1.0.)
        S->nq = 1;
        S->h = one;
        dewset(S->n, itol, rtol, atol, &rwork[S->lyh], &rwork[S->lewt]);
        for (i = 0; i < S->n; i++) {
            if (rwork[S->lewt + i] <= zero) {
                *istate = -3;
                return;
            }
            rwork[S->lewt + i] = one / rwork[S->lewt + i];
        }

        // Set or compute initial step size H0
        if (h0 == zero) {
            int niter, ier;
            dvhin(S->n, *t, &rwork[S->lyh], &rwork[lf0], f, rpar, ipar, *tout, S->uround, &rwork[S->lewt], itol, atol, y, &rwork[S->lacor], &h0, &niter, &ier);
            S->nfe += niter;
            if (ier != 0) {
                *istate = -3;
                return;
            }
        }

        // Adjust H0 if necessary to meet HMAX bound.
        {
            double rh = fabs(h0) * S->hmxi;
            if (rh > one) {
                h0 = h0 / rh;
            }
        }

        // Load H with H0 and scale YH(*,2) by H0.
        S->h = h0;
        dscal_(&S->n, &h0, &rwork[lf0], &int1);
    }

    // ===================================================================
    // Block D: Pre-step checks for continuation calls (ISTATE=2 or 3)
    // Block E: Main integration loop
    // ===================================================================
    if (*istate == 2 || *istate == 3) {
        nslast = S->nst;
        S->kuth = 0;

        // Check stop conditions based on ITASK before taking a step
        switch (*itask) {
            case 1:
                // Normal computation to TOUT
                if ((S->tn - *tout) * S->h >= zero) {
                    // Already past TOUT, interpolate
                    dvindy(S, *tout, 0, &rwork[S->lyh], S->nyh, y, &iflag);
                    if (iflag != 0) {
                        *istate = -3;
                        return;
                    }
                    *t = *tout;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;
                }
                break;

            case 2:
                // One-step mode, no pre-checks needed
                break;

            case 3:
                // Stop at internal mesh point at or beyond TOUT
                tp = S->tn - S->hu * (one + hun * S->uround);
                if ((tp - *tout) * S->h > zero) {
                    *istate = -3;
                    return;
                }
                if ((S->tn - *tout) * S->h >= zero) {
                    // Already at or past TOUT
                    dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                    *t = S->tn;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;
                }
                break;

            case 4:
                // Normal computation but don't overshoot TCRIT
                tcrit = rwork[0];
                if ((S->tn - tcrit) * S->h > zero) {
                    *istate = -3;
                    return;
                }
                if ((tcrit - *tout) * S->h < zero) {
                    *istate = -3;
                    return;
                }
                if ((S->tn - *tout) * S->h >= zero) {
                    // Already at or past TOUT
                    dvindy(S, *tout, 0, &rwork[S->lyh], S->nyh, y, &iflag);
                    if (iflag != 0) {
                        *istate = -3;
                        return;
                    }
                    *t = *tout;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;
                }
                // Check if approaching TCRIT
                hmx = fabs(S->tn) + fabs(S->h);
                ihit = (fabs(S->tn - tcrit) <= hun * S->uround * hmx);
                if (ihit) {
                    dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                    *t = S->tn;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;
                }
                tnext = S->tn + S->hnew * (one + four * S->uround);
                if ((tnext - tcrit) * S->h > zero) {
                    S->h = (tcrit - S->tn) * (one - four * S->uround);
                    S->kuth = 1;
                }
                break;

            case 5:
                // One step without passing TCRIT
                tcrit = rwork[0];
                if ((S->tn - tcrit) * S->h > zero) {
                    *istate = -3;
                    return;
                }
                hmx = fabs(S->tn) + fabs(S->h);
                ihit = (fabs(S->tn - tcrit) <= hun * S->uround * hmx);
                if (ihit) {
                    dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                    *t = S->tn;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;
                }
                tnext = S->tn + S->hnew * (one + four * S->uround);
                if ((tnext - tcrit) * S->h > zero) {
                    S->h = (tcrit - S->tn) * (one - four * S->uround);
                    S->kuth = 1;
                }
                break;
        }
    }

    // ===================================================================
    // Block E: Main integration loop
    // ===================================================================
    while (1) {
        // Check for too many steps
        if ((S->nst - nslast) >= S->mxstep) {
            // ISTATE = -1: Maximum number of steps exceeded
            dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
            *t = S->tn;
            *istate = -1;
            dvode_set_optional_output(S, rwork, iwork);
            rwork[11] = S->h;
            iwork[14] = S->nq;
            return;
        }

        // Update error weights
        dewset(S->n, itol, rtol, atol, &rwork[S->lyh], &rwork[S->lewt]);
        for (i = 0; i < S->n; i++) {
            if (rwork[S->lewt + i] <= zero) {
                // ISTATE = -6: EWT(i) became zero
                dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                *t = S->tn;
                *istate = -6;
                dvode_set_optional_output(S, rwork, iwork);
                rwork[11] = S->h;
                iwork[14] = S->nq;
                return;
            }
            rwork[S->lewt + i] = one / rwork[S->lewt + i];
        }

        // Check for too much accuracy requested
        tolsf = S->uround * dvnorm(S->n, &rwork[S->lyh], &rwork[S->lewt]);
        if (tolsf > one) {
            tolsf = tolsf * two;
            if (S->nst == 0) {
                *istate = -3;
                return;
            }
            // ISTATE = -2: Too much accuracy requested
            dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
            *t = S->tn;
            *istate = -2;
            rwork[13] = tolsf;
            dvode_set_optional_output(S, rwork, iwork);
            rwork[11] = S->h;
            iwork[14] = S->nq;
            return;
        }

        // Warning for TN+H=TN (we skip XERRWD calls per user requirement)
        // Just track S->mxhnil counter if needed

        // Call DVSTEP
        dvstep(S, y, &rwork[S->lyh], S->nyh, &rwork[S->lyh],
               &rwork[S->lewt], &rwork[S->lsavf], &rwork[S->lacor],
               &rwork[S->lwm], &iwork[S->liwm], f, jac, rpar, ipar);

        // Branch on KFLAG
        kgo = 1 - S->kflag;

        if (kgo == 1) {
            // ===================================================================
            // Block F: Successful step (KFLAG = 0)
            // ===================================================================
            S->init = 1;
            S->kuth = 0;

            // Check completion based on ITASK
            switch (*itask) {
                case 1:
                    // ITASK=1: If TOUT reached, interpolate
                    if ((S->tn - *tout) * S->h < zero) {
                        continue;  // Take another step
                    }
                    dvindy(S, *tout, 0, &rwork[S->lyh], S->nyh, y, &iflag);
                    *t = *tout;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;

                case 2:
                    // ITASK=2: One-step mode, always return
                    dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                    *t = S->tn;
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;

                case 3:
                    // ITASK=3: Stop at mesh point at or beyond TOUT
                    if ((S->tn - *tout) * S->h >= zero) {
                        dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                        *t = S->tn;
                        *istate = 2;
                        dvode_set_optional_output(S, rwork, iwork);
                        return;
                    }
                    continue;  // Take another step

                case 4:
                    // ITASK=4: Check TOUT or TCRIT
                    if ((S->tn - *tout) * S->h >= zero) {
                        dvindy(S, *tout, 0, &rwork[S->lyh], S->nyh, y, &iflag);
                        *t = *tout;
                        *istate = 2;
                        dvode_set_optional_output(S, rwork, iwork);
                        return;
                    }
                    hmx = fabs(S->tn) + fabs(S->h);
                    ihit = (fabs(S->tn - tcrit) <= hun * S->uround * hmx);
                    if (ihit) {
                        dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                        *t = tcrit;
                        *istate = 2;
                        dvode_set_optional_output(S, rwork, iwork);
                        return;
                    }
                    tnext = S->tn + S->hnew * (one + four * S->uround);
                    if ((tnext - tcrit) * S->h > zero) {
                        S->h = (tcrit - S->tn) * (one - four * S->uround);
                        S->kuth = 1;
                    }
                    continue;  // Take another step

                case 5:
                    // ITASK=5: One step without passing TCRIT
                    hmx = fabs(S->tn) + fabs(S->h);
                    ihit = (fabs(S->tn - tcrit) <= hun * S->uround * hmx);
                    dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
                    *t = S->tn;
                    if (ihit) {
                        *t = tcrit;
                    }
                    *istate = 2;
                    dvode_set_optional_output(S, rwork, iwork);
                    return;
            }

        } else if (kgo == 2) {
            // ===================================================================
            // Block H: KFLAG = -1 (error test failure)
            // ===================================================================
            // Compute IMXER (index of largest error component)
            double big = zero;
            int imxer = 1;
            for (i = 0; i < S->n; i++) {
                double size = fabs(rwork[S->lacor + i] * rwork[S->lewt + i]);
                if (size > big) {
                    big = size;
                    imxer = i + 1;  // Fortran 1-based indexing for output
                }
            }
            iwork[15] = imxer;

            // Set outputs and return
            dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
            *t = S->tn;
            *istate = -4;
            dvode_set_optional_output(S, rwork, iwork);
            rwork[11] = S->h;
            iwork[14] = S->nq;
            return;

        } else if (kgo == 3) {
            // ===================================================================
            // Block H: KFLAG = -2 (convergence failure)
            // ===================================================================
            // Compute IMXER (index of largest error component)
            double big = zero;
            int imxer = 1;
            for (i = 0; i < S->n; i++) {
                double size = fabs(rwork[S->lacor + i] * rwork[S->lewt + i]);
                if (size > big) {
                    big = size;
                    imxer = i + 1;  // Fortran 1-based indexing for output
                }
            }
            iwork[15] = imxer;

            // Set outputs and return
            dcopy_(&S->n, &rwork[S->lyh], &int1, y, &int1);
            *t = S->tn;
            *istate = -5;
            dvode_set_optional_output(S, rwork, iwork);
            rwork[11] = S->h;
            iwork[14] = S->nq;
            return;
        }
    }
}
