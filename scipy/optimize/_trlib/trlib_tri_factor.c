/* MIT License
 *
 * Copyright (c) 2016--2017 Felix Lenders
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "trlib_private.h"
#include "trlib.h"

trlib_int_t trlib_tri_factor_min(
    trlib_int_t nirblk, trlib_int_t *irblk, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t *neglin, trlib_flt_t radius,
    trlib_int_t itmax, trlib_flt_t tol_rel, trlib_flt_t tol_newton_tiny,
    trlib_int_t pos_def, trlib_int_t equality,
    trlib_int_t *warm0, trlib_flt_t *lam0, trlib_int_t *warm, trlib_flt_t *lam,
    trlib_int_t *warm_leftmost, trlib_int_t *ileftmost, trlib_flt_t *leftmost,
    trlib_int_t *warm_fac0, trlib_flt_t *diag_fac0, trlib_flt_t *offdiag_fac0,
    trlib_int_t *warm_fac, trlib_flt_t *diag_fac, trlib_flt_t *offdiag_fac,
    trlib_flt_t *sol0, trlib_flt_t *sol, trlib_flt_t *ones, trlib_flt_t *fwork,
    trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
    trlib_int_t *timing, trlib_flt_t *obj, trlib_int_t *iter_newton, trlib_int_t *sub_fail) {
    // use notation of Gould paper
    // h = h(lam) denotes solution of (T+lam I) * h = -lin

    trlib_int_t *leftmost_timing = NULL;
    trlib_int_t *eigen_timing = NULL;
    // local variables
    #if TRLIB_MEASURE_TIME
        struct timespec verystart, start, end;
        TRLIB_TIC(verystart)
        leftmost_timing = timing + 1 + TRLIB_SIZE_TIMING_LINALG;
        eigen_timing = timing + 1 + TRLIB_SIZE_TIMING_LINALG + trlib_leftmost_timing_size();
    #endif
    /* this is based on Theorem 5.8 in Gould paper,
     * the data for the first block has a 0 suffix,
     * the data for the \ell block has a l suffix */
    trlib_int_t n0 = irblk[1];                               // dimension of first block
    trlib_int_t nl;                                          // dimension of block corresponding to leftmost
    trlib_int_t nm0 = irblk[1]-1;                            // length of offdiagonal of first block
    trlib_int_t info_fac = 0;                                // factorization information
    trlib_int_t ret = 0;                                     // return code
    trlib_flt_t lam_pert = 0.0;                           // perturbation of leftmost eigenvalue as starting value for lam
    trlib_flt_t norm_sol0 = 0.0;                          // norm of h_0(lam)
    trlib_int_t jj = 0;                                      // local iteration counter
    trlib_flt_t dlam     = 0.0;                           // increment in newton iteration
    trlib_int_t inc = 1;                                     // increment in vector storage
    trlib_flt_t *w = fwork;                               // auxiliary vector to be used in newton iteration
    trlib_flt_t *diag_lam = fwork+(irblk[nirblk]);        // vector that holds diag + lam, could be saved if we would implement iterative refinement ourselves
    trlib_flt_t *work = fwork+2*(irblk[nirblk]);          // workspace for iterative refinement
    trlib_flt_t ferr = 0.0;                               // forward  error bound from iterative refinement
    trlib_flt_t berr = 0.0;                               // backward error bound from iterative refinement
    trlib_flt_t pert_low, pert_up;                        // lower and upper bound on perturbation of lambda
    trlib_flt_t dot = 0.0, dot2 = 0.0;                    // save dot products
    trlib_flt_t invD_norm_w_sq = 0.0;                     // || w ||_{D^-1}^2
    *iter_newton = 0;                                // newton iteration counter

    // FIXME: ensure diverse warmstarts work as expected

    // initialization:
    *sub_fail = 0;

    // set sol to 0 as a safeguard
    memset(sol, 0, irblk[nirblk]*sizeof(trlib_flt_t));

    // first make sure that lam0, h_0 is accurate
    TRLIB_PRINTLN_1("Solving trust region problem, radius %e; starting on first irreducible block", radius)
    // if only blocks changed that differ from the first, then there is nothing to do
    if (nirblk > 1 && *warm0) {
        TRLIB_DNRM2(norm_sol0, &n0, sol0, &inc)
        if (unicode) { TRLIB_PRINTLN_1("Solution provided via warmstart, \u03bb\u2080 = %e, \u2016h\u2080\u2016 = %e", *lam0, norm_sol0) }
        else { TRLIB_PRINTLN_1("Solution provided via warmstart, lam0 = %e, ||h0|| = %e", *lam0, norm_sol0) }
        if (norm_sol0-radius < 0.0) {
            if (unicode) { TRLIB_PRINTLN_1("  violates \u2016h\u2080\u2016 - radius \u2265 0, but is %e, switch to coldstart", norm_sol0-radius) }
            else { TRLIB_PRINTLN_1("  violates ||h0|| - radius >= 0, but is %e, switch to coldstart", norm_sol0-radius) }
            *warm0 = 0;
        }
    }
    if (nirblk == 1 || !*warm0) {
        // seek for lam0, h_0 with (T0+lam0*I) pos def and ||h_0(lam_0)|| = radius

        /* as a first step to initialize the newton iteration,
         *  find such a pair with the loosened requirement ||h_0(lam_0)|| >= radius */
        if(*warm0) {
            if(!*warm_fac0) {
                // factorize T + lam0 I
                TRLIB_DCOPY(&n0, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
                TRLIB_DAXPY(&n0, lam0, ones, &inc, diag_lam, &inc) // diag_lam <-- lam0 + diag_lam
                TRLIB_DCOPY(&n0, diag_lam, &inc, diag_fac0, &inc) // diag_fac0 <-- diag_lam
                TRLIB_DCOPY(&nm0, offdiag, &inc, offdiag_fac0, &inc) // offdiag_fac0 <-- offdiag
                TRLIB_DPTTRF(&n0, diag_fac0, offdiag_fac0, &info_fac) // compute factorization
                if (info_fac != 0) { *warm0 = 0; } // factorization failed, switch to coldstart
                else { *warm_fac0 = 1; }
            }
            if(*warm_fac0) {
                // solve for h0(lam0) and compute norm
                TRLIB_DCOPY(&n0, neglin, &inc, sol0, &inc) // sol0 <-- neglin
                TRLIB_DPTTRS(&n0, &inc, diag_fac0, offdiag_fac0, sol0, &n0, &info_fac) // sol <-- (T+lam0 I)^-1 sol
                if(info_fac!=0) {
                    if (unicode) { TRLIB_PRINTLN_2("Failure on computing h\u2080 upon initialization") }
                    else { TRLIB_PRINTLN_2("Failure on computing h0 upon initialization") }
                    TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE)
                }
                TRLIB_DNRM2(norm_sol0, &n0, sol0, &inc)
                if (norm_sol0 >= radius) { *warm0 = 1; } else { *warm0 = 0; }
            }
        }
        if(!*warm0) {
            *lam0 = 0.0;
            if (unicode) { TRLIB_PRINTLN_1("Coldstart. Seeking suitable initial \u03bb\u2080, starting with 0") }
            else { TRLIB_PRINTLN_1("Coldstart. Seeking suitable initial lam0, starting with 0") }
            TRLIB_DCOPY(&n0, diag, &inc, diag_fac0, &inc) // diag_fac0 <-- diag0
            TRLIB_DCOPY(&nm0, offdiag, &inc, offdiag_fac0, &inc) // offdiag_fac0 <-- offdiag0
            TRLIB_DCOPY(&n0, neglin, &inc, sol0, &inc) // sol0 <-- neglin0
            TRLIB_DPTTRF(&n0, diag_fac0, offdiag_fac0, &info_fac) // compute factorization
            if (info_fac == 0) {
                // test if lam0 = 0 is suitable
                TRLIB_DCOPY(&n0, neglin, &inc, sol0, &inc) // sol0 <-- neglin
                TRLIB_DPTTRS(&n0, &inc, diag_fac0, offdiag_fac0, sol0, &n0, &info_fac) // sol0 <-- T0^-1 sol0
                if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on computing h\u2080 upon initialization") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
                TRLIB_DNRM2(norm_sol0, &n0, sol0, &inc)
                if (norm_sol0<radius && equality) { info_fac = 1; } // in equality case, we have to find suitable lam
            }
            if (info_fac != 0) {
                if (unicode) { TRLIB_PRINTLN_1(" \u03bb\u2080 = 0 unsuitable \u2265 get leftmost ev of first block!") }
                else { TRLIB_PRINTLN_1(" lam0 = 0 unsuitable ==> get leftmost ev of first block!") }
                *sub_fail = trlib_leftmost_irreducible(irblk[1], diag, offdiag, *warm_leftmost, *leftmost, 1000, TRLIB_EPS_POW_75, verbose-2, unicode, " LM ", fout, leftmost_timing, leftmost, &jj); // ferr can safely be overwritten by computed leftmost for the moment as can jj with the number of rp iterations
                // if (*sub_fail != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_LM) } failure of leftmost: may lead to inefficiency, since what we are doing may be slow...
                // T - leftmost*I is singular, so do bisection to catch factorization and find suitable initial lambda
                pert_low = - TRLIB_EPS_POW_5 * fabs(*leftmost); // lower bound on allowed perturbation
                pert_up = 1.0/TRLIB_EPS; // upper bound on allowed perturbation
                jj = 0; // counter on number of tries
                lam_pert = 0.0;
                TRLIB_PRINTLN_1(" ")
                if (unicode) { TRLIB_PRINTLN_1(" perturb \u03bb\u2080 by safeguarded bisection to find suitable initial value") }
                else { TRLIB_PRINTLN_1(" perturb lam0 by safeguarded bisection to find suitable initial value") }
                while( jj < 50 ) {
                    if( jj % 20 == 0 ) {
                        if (unicode) { TRLIB_PRINTLN_2(" %2s%14s%14s%14s%3s%14s", "it", "low     ", "pert    ", "up      ", "pd", "  \u2016h\u2080\u2016 - radius") }
                        else { TRLIB_PRINTLN_2(" %2s%14s%14s%14s%3s%14s", "it", "low     ", "pert    ", "up      ", "pd", "  ||h0|| - radius") }
                    }
                    *lam0 = -(*leftmost) + lam_pert;
                    // factorize T + lam I
                    TRLIB_DCOPY(&n0, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
                    TRLIB_DAXPY(&n0, lam0, ones, &inc, diag_lam, &inc) // diag_lam <-- lam + diag_lam
                    TRLIB_DCOPY(&n0, diag_lam, &inc, diag_fac0, &inc) // diag_fac <-- diag_lam
                    TRLIB_DCOPY(&nm0, offdiag, &inc, offdiag_fac0, &inc) // offdiag_fac <-- offdiag
                    TRLIB_DPTTRF(&n0, diag_fac0, offdiag_fac0, &info_fac) // compute factorization
                    if(info_fac == 0) {
                        pert_up = lam_pert; // as this ensures a factorization, it provides a upper bound
                        TRLIB_DCOPY(&n0, neglin, &inc, sol0, &inc) // sol0 <-- neglin
                        TRLIB_DPTTRS(&n0, &inc, diag_fac0, offdiag_fac0, sol0, &n0, &info_fac) // sol0 <-- T0^-1 sol0
                        if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on computing h\u2080 in safeguarded initialization iteration") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
                        TRLIB_DNRM2(norm_sol0, &n0, sol0, &inc)
                        TRLIB_PRINTLN_2(" %2ld%14e%14e%14e%3s%14e", jj, pert_low, lam_pert, pert_up, " +", norm_sol0 - radius)
                        if(norm_sol0 >= radius) {
                            break;
                        }
                        else {
                            // lam_pert has to be decreased since we caught factorization but got solution that was too small
                            lam_pert = .5*(pert_low+pert_up);
                        }
                    }
                    else {
                        TRLIB_PRINTLN_2(" %2ld%14e%14e%14e%3s", jj, pert_low, lam_pert, pert_up, " -")
                        pert_low = lam_pert; // as factorization fails, it provides a upper bound
                        // now increase perturbation, either by bisection if there is a useful upper bound,
                        // otherwise by a small increment
                        if( pert_up == 1.0/TRLIB_EPS ) {
                            if( lam_pert == 0.0 ) { lam_pert = (1.0+fabs(*leftmost))*TRLIB_EPS_POW_75; } else { lam_pert = 2.0*lam_pert; }
                        }
                        else {
                            lam_pert = .5*(pert_low+pert_up);
                        }
                    }
                    ++jj;
                }
                // ensure that we get a factorization and compute solution with it
                if(info_fac != 0) {
                    lam_pert = pert_up;
                    *lam0 = -(*leftmost) + lam_pert;
                    TRLIB_DCOPY(&n0, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
                    TRLIB_DAXPY(&n0, lam0, ones, &inc, diag_lam, &inc) // diag_lam <-- lam + diag_lam
                    TRLIB_DCOPY(&n0, diag_lam, &inc, diag_fac0, &inc) // diag_fac <-- diag_lam
                    TRLIB_DCOPY(&nm0, offdiag, &inc, offdiag_fac0, &inc) // offdiag_fac <-- offdiag
                    TRLIB_DPTTRF(&n0, diag_fac0, offdiag_fac0, &info_fac) // compute factorization
                    if(info_fac != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_FACTOR) }
                    TRLIB_DPTTRS(&n0, &inc, diag_fac0, offdiag_fac0, sol0, &n0, &info_fac) // sol0 <-- T0^-1 sol0
                    if (info_fac != 0) {
                        if (unicode) { TRLIB_PRINTLN_2("Failure on computing h\u2080 upon after factorization was ensured") }
                        else { TRLIB_PRINTLN_2("Failure on computing h\u2080 upon after factorization was ensured") }
                        TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE)
                    }
                    TRLIB_DNRM2(norm_sol0, &n0, sol0, &inc)
                }
            }
        }
    }
    if (norm_sol0 >= radius) { // perform newton iteration if possible
        /* now a suitable pair lam0, h0 has been found.
         * perform a newton iteration on 0 = 1/||h0(lam0)|| - 1/radius */
        if (unicode) { TRLIB_PRINTLN_1("Starting Newton iteration for \u03bb\u2080 with initial choice %e", *lam0) }
        else { TRLIB_PRINTLN_1("Starting Newton iteration for lam0 with initial choice %e", *lam0) }
        while (1) {
            /* compute newton correction to lam, by
                (1) Factoring T0 + lam0 I = LDL^T
                (2) Solving (T0+lam0 I)*h0 = -lin
                (3) L*w = h0/||h0||
                (4) compute increment (||h0||-Delta)/Delta/||w||_{D^-1}^2 */

            // steps (1) and (2) have already been performed on initializaton or previous iteration

            /* step (3) L*w = h/||h||
               compute ||w||_{D^-1}^2 in same loop */
            ferr = 1.0/norm_sol0; TRLIB_DCOPY(&n0, sol0, &inc, w, &inc) TRLIB_DSCAL(&n0, &ferr, w, &inc) // w <-- sol/||sol||
            invD_norm_w_sq = w[0]*w[0]/diag_fac0[0];
            for( jj = 1; jj < n0; ++jj ) { w[jj] = w[jj] - offdiag_fac0[jj-1]*w[jj-1]; invD_norm_w_sq += w[jj]*w[jj]/diag_fac0[jj]; }

            // step (4) compute increment (||h||-Delta)/Delta/||w||_{D^-1}^2
            dlam = (norm_sol0-radius)/(radius*invD_norm_w_sq);

            // iteration completed
            *iter_newton += 1;

            // test if dlam is not tiny or newton limit exceeded, return eventually
            if ( *lam0 + dlam == *lam0 || fabs(dlam) <= tol_newton_tiny * fmax(1.0, fabs(*lam0)) || *iter_newton > itmax) {
                if (unicode) { TRLIB_PRINTLN_1("%s%e%s%e", "Newton breakdown, d\u03bb = ", dlam, " \u03bb = ", *lam0) }
                else { TRLIB_PRINTLN_1("%s%e%s%e", "Newton breakdown, dlam = ", dlam, " \u03bb = ", *lam0) }
                if(*iter_newton > itmax) { ret = TRLIB_TTR_ITMAX; break; }
                ret = TRLIB_TTR_NEWTON_BREAK; break;
            }

            // prepare next iteration

            // update lam
            *lam0 += dlam;

            // step (1) Factoring T0 + lam0 I = LDL^T
            TRLIB_DCOPY(&n0, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
            TRLIB_DAXPY(&n0, lam0, ones, &inc, diag_lam, &inc) // diag_lam <-- lam + diag_lam
            TRLIB_DCOPY(&n0, diag_lam, &inc, diag_fac0, &inc) // diag_fac <-- diag_lam
            TRLIB_DCOPY(&nm0, offdiag, &inc, offdiag_fac0, &inc) // offdiag_fac <-- offdiag
            TRLIB_DPTTRF(&n0, diag_fac0, offdiag_fac0, &info_fac) // compute factorization
            if (info_fac != 0) {
                if (unicode) { TRLIB_PRINTLN_2("Fail on factorization, \u03bb = %e, d\u03bb = %e! Exiting Newton Iteration", *lam0, dlam) }
                else {TRLIB_PRINTLN_2("Fail on factorization, lam = %e, dlam = %e! Exiting Newton Iteration", *lam0, dlam) }
                TRLIB_RETURN(TRLIB_TTR_FAIL_FACTOR)
            }

            // step (2) Solving (T+lam I)*h = -lin
            TRLIB_DCOPY(&n0, neglin, &inc, sol0, &inc) // sol <-- neglin
            TRLIB_DPTTRS(&n0, &inc, diag_fac0, offdiag_fac0, sol0, &n0, &info_fac) // sol <-- (T+lam I)^-1 sol
            if (info_fac != 0) {
                if (unicode) { TRLIB_PRINTLN_2("Failure on computing h\u2080 in newton iteration") }
                else { TRLIB_PRINTLN_2("Failure on computing h0 in newton iteration") }
                TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE)
            }
            if (refine) { TRLIB_DPTRFS(&n0, &inc, diag_lam, offdiag, diag_fac0, offdiag_fac0, neglin, &n0, sol0, &n0, &ferr, &berr, work, &info_fac) }
            if (info_fac != 0) {
                if (unicode) { TRLIB_PRINTLN_2("Failure on computing h\u2080 in refinement in newton iteration") }
                else { TRLIB_PRINTLN_2("Failure on computing h\u2080 in refinement in newton iteration") }
                TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE)
            }
            TRLIB_DNRM2(norm_sol0, &n0, sol0, &inc)

            if (*iter_newton % 20 == 1) {
                if (unicode) { TRLIB_PRINTLN_1("%6s%14s%14s%14s", " iter ", "       \u03bb      ", "      d\u03bb      ", " \u2016h\u2080(\u03bb)\u2016-radius") }
                else { TRLIB_PRINTLN_1("%6s%14s%14s%14s", " iter ", "     lam      ", "     dlam     ", "  tr resdidual") }
            }
            TRLIB_PRINTLN_1("%6ld%14e%14e%14e", *iter_newton, *lam0, dlam, norm_sol0-radius)

            // test for convergence
            if (norm_sol0 - radius <= tol_rel * radius) {
                // what if norm_sol < radius in a significant way?
                // theory tells this should not happen...

                ret = TRLIB_TTR_CONV_BOUND; break;
            }
        }
    }
    *warm0 = 1;

    // test if we trust region problem is solved on first irreducible with sufficient accuracy,
    // otherwise build linear combination with eigenvector to leftmost that solves trust region constraint
    if ( fabs(radius - norm_sol0) >= TRLIB_EPS_POW_5*radius ) {
        if(*lam0 == 0.0 && !equality) { ret = TRLIB_TTR_CONV_INTERIOR; }
        else {
            if (unicode) { TRLIB_PRINTLN_1(" Found \u03bb\u2080 with tr residual %e! Bail out with h\u2080 + \u03b1 eig", radius - norm_sol0) }
            else { TRLIB_PRINTLN_1(" Found lam0 with tr residual %e! Bail out with h0 + alpha eig", radius - norm_sol0) }
            *sub_fail = trlib_eigen_inverse(n0, diag, offdiag,
                    *leftmost, 10, TRLIB_EPS_POW_5, ones,
                    diag_fac, offdiag_fac, sol,
                    verbose-2, unicode, " EI", fout, eigen_timing, &ferr, &berr, &jj); // can safely overwrite ferr, berr, jj with results. only interesting: eigenvector
            if (*sub_fail != 0 && *sub_fail != -1) { TRLIB_PRINTLN_2("Failure in eigenvector computation: %ld", *sub_fail) TRLIB_RETURN(TRLIB_TTR_FAIL_EIG) }
            if (*sub_fail == -1) { TRLIB_PRINTLN_2("In eigenvector computation itmax reached, continue with approximate eigenvector") }
            // compute solution as linear combination of h0 and eigenvector
            // ||h0 + t*eig||^2 = ||h_0||^2 + t * <h0, eig> + t^2 = radius^2
            TRLIB_DDOT(dot, &n0, sol0, &inc, sol, &inc); // dot = <h0, eig>
            trlib_quadratic_zero( norm_sol0*norm_sol0 - radius*radius, 2.0*dot, TRLIB_EPS_POW_75, verbose - 3, unicode, prefix, fout, &ferr, &berr);
            // select solution that corresponds to smaller objective
            // quadratic as a function of t without offset
            // q(t) = 1/2 * leftmost * t^2 + (leftmost * <eig, h0> + <eig, lin>) * t
            TRLIB_DDOT(dot2, &n0, sol, &inc, neglin, &inc) // dot2 = - <eig, lin>
            if( .5*(*leftmost)*ferr*ferr + ((*leftmost)*dot - dot2)*ferr <= .5*(*leftmost)*berr*berr + ((*leftmost)*dot - dot2)*berr) {
                TRLIB_DAXPY(&n0, &ferr, sol, &inc, sol0, &inc)
            }
            else {
                TRLIB_DAXPY(&n0, &berr, sol, &inc, sol0, &inc)
            }
            ret = TRLIB_TTR_HARD_INIT_LAM;
        }
    }


    /* now in a situation were accurate lam0, h_0 exists to first irreducible block
     * invoke Theorem 5.8:
     * (i)  if lam0 >= -leftmost the pair lam0, h_0 solves the problem
     * (ii) if lam0 < -leftmost a solution has to be constructed to lam = -leftmost */

    // quick exit: only one irreducible block
    if (nirblk == 1) {
        *lam = *lam0; *warm = 1;
        TRLIB_DCOPY(&n0, sol0, &inc, sol, &inc) // sol <== sol0
        // compute objective. first store 2*gradient in w, then compute obj = .5*(sol, w)
        TRLIB_DCOPY(&n0, neglin, &inc, w, &inc) ferr = -2.0; TRLIB_DSCAL(&n0, &ferr, w, &inc) ferr = 1.0; // w <-- -2 neglin
        TRLIB_DLAGTM("N", &n0, &inc, &ferr, offdiag, diag, offdiag, sol, &n0, &ferr, w, &n0) // w <-- T*sol + w
        TRLIB_DDOT(dot, &n0, sol, &inc, w, &inc) *obj = 0.5*dot; // obj = .5*(sol, w)
        TRLIB_RETURN(ret)
    }

    // now that we have accurate lam, h_0 invoke Theorem 5.8
    // check if lam <= leftmost --> in that case the first block information describes everything
    if (unicode) { TRLIB_PRINTLN_1("\nCheck if \u03bb\u2080 provides global solution, get leftmost ev for irred blocks") }
    else { TRLIB_PRINTLN_1("\nCheck if lam0 provides global solution, get leftmost ev for irred blocks") }
    if(!*warm_leftmost) {
        *sub_fail = trlib_leftmost(nirblk, irblk, diag, offdiag, 0, leftmost[nirblk-1], 1000, TRLIB_EPS_POW_75, verbose-2, unicode, " LM ", fout, leftmost_timing, ileftmost, leftmost);
        *warm_leftmost = 1;
    }
    TRLIB_PRINTLN_1("    leftmost = %e (block %ld)", leftmost[*ileftmost], *ileftmost)
    if(*lam0 >= -leftmost[*ileftmost]) {
        if (unicode) { TRLIB_PRINTLN_1("  \u03bb\u2080 \u2265 -leftmost \u21d2 \u03bb = \u03bb\u2080, exit: h\u2080(\u03bb\u2080)") }
        else { TRLIB_PRINTLN_1("  lam0 >= -leftmost => lam = lam0, exit: h0(lam0)") }
        *lam = *lam0; *warm = 1;
        TRLIB_DCOPY(&n0, sol0, &inc, sol, &inc) // sol <== sol0
        // compute objective. first store 2*gradient in w, then compute obj = .5*(sol, w)
        TRLIB_DCOPY(&n0, neglin, &inc, w, &inc) ferr = -2.0; TRLIB_DSCAL(&n0, &ferr, w, &inc) ferr = 1.0; // w <-- -2 neglin
        TRLIB_DLAGTM("N", &n0, &inc, &ferr, offdiag, diag, offdiag, sol, &n0, &ferr, w, &n0) // w <-- T*sol + w
        TRLIB_DDOT(dot, &n0, sol, &inc, w, &inc) *obj = 0.5*dot; // obj = .5*(sol, w)
        TRLIB_RETURN(ret)
    }
    else {
        if (unicode) { TRLIB_PRINTLN_1("  -leftmost > \u03bb\u2080 \u21d2 \u03bb = -leftmost, exit: h\u2080(-leftmost) + \u03b1 u") }
        else  { TRLIB_PRINTLN_1("  -leftmost > lam0 => lam = -leftmost, exit: h0(-leftmost) + alpha u") }

        // Compute solution of (T0 - leftmost*I)*h0 = neglin
        *lam = -leftmost[*ileftmost]; *warm = 1;
        TRLIB_DCOPY(&n0, neglin, &inc, sol, &inc) // neglin <-- sol
        if(!*warm_fac){
            TRLIB_DCOPY(&n0, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
            TRLIB_DAXPY(&n0, lam, ones, &inc, diag_lam, &inc) // diag_lam <-- lam + diag_lam
            TRLIB_DCOPY(&n0, diag_lam, &inc, diag_fac, &inc) // diag_fac <-- diag_lam
            TRLIB_DCOPY(&nm0, offdiag, &inc, offdiag_fac, &inc) // offdiag_fac <-- offdiag
            TRLIB_DPTTRF(&n0, diag_fac, offdiag_fac, &info_fac) // compute factorization
            if (info_fac != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_FACTOR) }
        }
        *warm_fac = 1;
        TRLIB_DPTTRS(&n0, &inc, diag_fac, offdiag_fac, sol, &n0, &info_fac) // sol <-- (T+lam I)^-1 sol
        if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on computing h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
        if (refine) { TRLIB_DPTRFS(&n0, &inc, diag_lam, offdiag, diag_fac, offdiag_fac, neglin, &n0, sol, &n0, &ferr, &berr, work, &info_fac) }
        if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on refining h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
        TRLIB_DNRM2(norm_sol0, &n0, sol, &inc)

        // compute normalized eigenvector u corresponding to leftmost of block ileftmost
        nl = irblk[*ileftmost+1]-irblk[*ileftmost];
        *sub_fail = trlib_eigen_inverse(nl, diag+irblk[*ileftmost], offdiag+irblk[*ileftmost],
                leftmost[*ileftmost], 10, TRLIB_EPS_POW_5, ones,
                diag_fac+irblk[*ileftmost], offdiag_fac+irblk[*ileftmost],
                sol+irblk[*ileftmost],
                verbose-2, unicode, " EI", fout, eigen_timing, &ferr, &berr, &jj); // can safely overwrite ferr, berr, jj with results. only interesting: eigenvector
        if (*sub_fail != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_EIG) }

        // solution is of form [h,0,...,0,alpha*u,0,...,0]
        // alpha = sqrt( radius^2 - ||h||^2 )
        ferr = sqrt( radius*radius - norm_sol0*norm_sol0 );
        TRLIB_DSCAL(&nl, &ferr, sol+irblk[*ileftmost], &inc)

        if (unicode) { TRLIB_PRINTLN_1("    with \u2016h\u2080(-leftmost)\u2016 = %e, \u03b1 = %e", norm_sol0, ferr) }
        else { TRLIB_PRINTLN_1("    with ||h0(-leftmost)|| = %e, alpha = %e", norm_sol0, ferr) }

        ret = TRLIB_TTR_HARD;

        // compute objective. first store 2*gradient in w, then compute obj = .5*(sol, w)
        *obj = 0.5*leftmost[*ileftmost]*ferr*ferr;
        TRLIB_DCOPY(&n0, neglin, &inc, w, &inc) ferr = -2.0; TRLIB_DSCAL(&n0, &ferr, w, &inc) ferr = 1.0; // w <-- -2 neglin
        TRLIB_DLAGTM("N", &n0, &inc, &ferr, offdiag, diag, offdiag, sol, &n0, &ferr, w, &n0) // w <-- T*sol + w
        TRLIB_DDOT(dot, &n0, sol, &inc, w, &inc) *obj = *obj+0.5*dot; // obj = .5*(sol, w)
        TRLIB_RETURN(ret);
    }
}

trlib_int_t trlib_tri_factor_regularized_umin(
    trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t *neglin, trlib_flt_t lam,
    trlib_flt_t *sol,
    trlib_flt_t *ones, trlib_flt_t *fwork,
    trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
    trlib_int_t *timing, trlib_flt_t *norm_sol, trlib_int_t *sub_fail) {

    // local variables
    #if TRLIB_MEASURE_TIME
        struct timespec verystart, start, end;
        TRLIB_TIC(verystart)
    #endif

    trlib_flt_t *diag_lam = fwork;        // vector that holds diag + lam, could be saved if we would implement iterative refinement ourselves
    trlib_flt_t *diag_fac = fwork+n;      // vector that holds diagonal of factor of diag + lam
    trlib_flt_t *offdiag_fac = fwork+2*n; // vector that holds offdiagonal of factor of diag + lam
    trlib_flt_t *work = fwork+3*n;        // workspace for iterative refinement
    trlib_flt_t ferr = 0.0;               // forward  error bound from iterative refinement
    trlib_flt_t berr = 0.0;               // backward error bound from iterative refinement
    trlib_int_t inc = 1;                  // vector increment
    trlib_int_t info_fac = 0;             // LAPACK return code
    trlib_int_t nm = n-1;

    // factorize T + lam0 I
    TRLIB_DCOPY(&n, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
    TRLIB_DAXPY(&n, &lam, ones, &inc, diag_lam, &inc) // diag_lam <-- lam0 + diag_lam
    TRLIB_DCOPY(&n, diag_lam, &inc, diag_fac, &inc) // diag_fac <-- diag_lam
    TRLIB_DCOPY(&nm, offdiag, &inc, offdiag_fac, &inc) // offdiag_fac <-- offdiag
    TRLIB_DPTTRF(&n, diag_fac, offdiag_fac, &info_fac) // compute factorization
    if (info_fac != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_FACTOR); } // factorization failed, switch to coldastart

    TRLIB_DCOPY(&n, neglin, &inc, sol, &inc) // sol <-- neglin
    TRLIB_DPTTRS(&n, &inc, diag_fac, offdiag_fac, sol, &n, &info_fac) // sol <-- (T+lam I)^-1 sol
    if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on backsolve for h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
    if (refine) { TRLIB_DPTRFS(&n, &inc, diag_lam, offdiag, diag_fac, offdiag_fac, neglin, &n, sol, &n, &ferr, &berr, work, &info_fac) }
    if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on iterative refinement for h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }

    TRLIB_DNRM2(*norm_sol, &n, sol, &inc)
    TRLIB_RETURN(TRLIB_TTR_CONV_INTERIOR);
}

trlib_int_t trlib_tri_factor_get_regularization(
    trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t *neglin, trlib_flt_t *lam,
    trlib_flt_t sigma, trlib_flt_t sigma_l, trlib_flt_t sigma_u,
    trlib_flt_t *sol,
    trlib_flt_t *ones, trlib_flt_t *fwork,
    trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
    trlib_int_t *timing, trlib_flt_t *norm_sol, trlib_int_t *sub_fail) {

    // local variables
    #if TRLIB_MEASURE_TIME
        struct timespec verystart, start, end;
        TRLIB_TIC(verystart)
    #endif

    trlib_flt_t *diag_lam = fwork;        // vector that holds diag + lam, could be saved if we would implement iterative refinement ourselves
    trlib_flt_t *diag_fac = fwork+n;      // vector that holds diagonal of factor of diag + lam
    trlib_flt_t *offdiag_fac = fwork+2*n; // vector that holds offdiagonal of factor of diag + lam
    trlib_flt_t *work = fwork+3*n;        // workspace for iterative refinement
    trlib_flt_t *aux  = fwork+5*n;        // auxiliary vector ds/n
    trlib_flt_t ferr = 0.0;               // forward  error bound from iterative refinement
    trlib_flt_t berr = 0.0;               // backward error bound from iterative refinement
    trlib_int_t inc = 1;                  // vector increment
    trlib_int_t info_fac;                 // LAPACK return code
    trlib_int_t nm = n-1;
    trlib_flt_t lambda_l = 0.0;           // lower bound on lambda
    trlib_flt_t lambda_u = 1e20;          // upper bound on lambda
    trlib_int_t jj = 0;                   // local loop variable
    trlib_flt_t dlam = 0.0;               // step in lambda
    trlib_flt_t dn = 0.0;                 // derivative of norm

    // get suitable lambda for which factorization exists
    info_fac = 1;
    while(info_fac != 0 && jj < 500) {
        // factorize T + lam0 I
        TRLIB_DCOPY(&n, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
        TRLIB_DAXPY(&n, lam, ones, &inc, diag_lam, &inc) // diag_lam <-- lam0 + diag_lam
        TRLIB_DCOPY(&n, diag_lam, &inc, diag_fac, &inc) // diag_fac <-- diag_lam
        TRLIB_DCOPY(&nm, offdiag, &inc, offdiag_fac, &inc) // offdiag_fac <-- offdiag
        TRLIB_DPTTRF(&n, diag_fac, offdiag_fac, &info_fac) // compute factorization
        if(info_fac == 0) { break; }
        if(*lam > lambda_u) { break; }
        lambda_l = *lam;
        //if (*lam == 0.0) { *lam = 1.0; }
        *lam = 2.0 * (*lam); jj++;
    }
    if (info_fac != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_FACTOR); } // factorization failed
    TRLIB_PRINTLN_1("Initial Regularization Factor that allows Cholesky: %e", *lam);

    TRLIB_DCOPY(&n, neglin, &inc, sol, &inc) // sol <-- neglin
    TRLIB_DPTTRS(&n, &inc, diag_fac, offdiag_fac, sol, &n, &info_fac) // sol <-- (T+lam I)^-1 sol
    if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on backsolve for h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
    if (refine) { TRLIB_DPTRFS(&n, &inc, diag_lam, offdiag, diag_fac, offdiag_fac, neglin, &n, sol, &n, &ferr, &berr, work, &info_fac) }
    if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on iterative refinement for h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }

    TRLIB_DNRM2(*norm_sol, &n, sol, &inc)

    jj = 0;
    TRLIB_PRINTLN_2("%ld\t Reg %e\t Reg/Norm %e\t lb %e ub %e", jj, *lam, *lam/(*norm_sol), sigma_l, sigma_u);

    // check if accetable
    if( *norm_sol * sigma_l <= *lam && *lam <= *norm_sol * sigma_u ) {
        TRLIB_PRINTLN_1("Exit with Regularization Factor %e and Reg/Norm %e", *lam, *lam/(*norm_sol))
        TRLIB_RETURN(TRLIB_TTR_CONV_INTERIOR);
    }
    else {
        /* do safeguarded newton iteration on f(lam) = lam/n(lam) with n(lam) = ||s(lam)||
         * then dn = 1/n <s, ds> with - (H+lam I) ds = s
         * thus get - f/df = (n*lam - n*n*sigma) / ( lam dn - n) with dn = <s, ds/n> and (H + lam I) ds/n = - s/n
         * note that f is increasing for lam such that (H+lam I) is spd
         */

        jj = 0;

        while(jj < 500) {

            // first get the vector ds/n
            TRLIB_DCOPY(&n, sol, &inc, aux, &inc) // aux <-- sol
            dn = -1.0/(*norm_sol); // scaling for right hand side
            TRLIB_DSCAL(&n, &dn, aux, &inc) // aux <-- -sol/||norm_sol||
            TRLIB_DDOT(dn, &n, sol, &inc, aux, &inc) // dn = <s, ds/n>

            // compute step correction
            dlam = (*lam*(*norm_sol)-*norm_sol*(*norm_sol)*sigma) / (*lam*dn - *norm_sol);

            // check feasibility of step
            if (*lam + dlam <= lambda_u && lambda_l <= *lam + dlam) {
                *lam = *lam + dlam;
            }
            else { *lam = .5*( lambda_l + lambda_u); }

            // compute next function value

            // factorize T + lam0 I
            TRLIB_DCOPY(&n, diag, &inc, diag_lam, &inc) // diag_lam <-- diag
            TRLIB_DAXPY(&n, lam, ones, &inc, diag_lam, &inc) // diag_lam <-- lam0 + diag_lam
            TRLIB_DCOPY(&n, diag_lam, &inc, diag_fac, &inc) // diag_fac <-- diag_lam
            TRLIB_DCOPY(&nm, offdiag, &inc, offdiag_fac, &inc) // offdiag_fac <-- offdiag
            TRLIB_DPTTRF(&n, diag_fac, offdiag_fac, &info_fac) // compute factorization
            if (info_fac != 0) { TRLIB_RETURN(TRLIB_TTR_FAIL_FACTOR); } // factorization failed

            TRLIB_DCOPY(&n, neglin, &inc, sol, &inc) // sol <-- neglin
            TRLIB_DPTTRS(&n, &inc, diag_fac, offdiag_fac, sol, &n, &info_fac) // sol <-- (T+lam I)^-1 sol
            if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on backsolve for h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }
            if (refine) { TRLIB_DPTRFS(&n, &inc, diag_lam, offdiag, diag_fac, offdiag_fac, neglin, &n, sol, &n, &ferr, &berr, work, &info_fac) }
            if (info_fac != 0) { TRLIB_PRINTLN_2("Failure on iterative refinement for h") TRLIB_RETURN(TRLIB_TTR_FAIL_LINSOLVE) }

            TRLIB_DNRM2(*norm_sol, &n, sol, &inc)

            jj++;
            TRLIB_PRINTLN_2("%ld\t Reg %e\t Reg/Norm %e\t lb %e ub %e", jj, *lam, *lam/(*norm_sol), sigma_l, sigma_u);

            // check if accetable
            if( *norm_sol * sigma_l <= *lam && *lam <= *norm_sol * sigma_u ) {
                TRLIB_PRINTLN_1("Exit with Regularization Factor %e and Reg/Norm %e", *lam, *lam/(*norm_sol))
                TRLIB_RETURN(TRLIB_TTR_CONV_INTERIOR);
            }
            else { // contract bounds
                if(*lam > *norm_sol * sigma_u) { lambda_u = *lam; }
                if(*lam < *norm_sol * sigma_l) { lambda_l = *lam; }
            }

        }
        TRLIB_PRINTLN_1("Failure: no convergence to determine regularaization factor, last iterate %e", *lam) TRLIB_RETURN(TRLIB_TTR_NEWTON_BREAK)

    }

}

trlib_int_t trlib_tri_factor_regularize_posdef(
    trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t tol_away, trlib_flt_t security_step, trlib_flt_t *regdiag) {

    /* modify diagonal to be able to factorize
       Cholesky recurrence for diagonal is
       diag_fac[0] = diag[0]
       diag_fac[i+1] = diag[i+1] - offdiag[i]*offdiag[i] / diag_fac[i]
       we have to ensure diag_fac > 0 */

    trlib_flt_t diag_fac = 0.0;
    trlib_int_t pivot = 0;

    regdiag[0] = 0.0;
    if (diag[0] <= tol_away) { regdiag[0] = security_step*tol_away; }
    diag_fac = diag[0] + regdiag[0];

    for(pivot = 0; pivot < n-1; ++pivot) {
        regdiag[pivot+1] = 0.0;
        if ( diag[pivot+1] - offdiag[pivot]*offdiag[pivot]/diag_fac <= tol_away * diag_fac ) {
            regdiag[pivot+1] = security_step * fabs(offdiag[pivot]*offdiag[pivot]/diag_fac - diag[pivot+1]);
        }
        diag_fac = diag[pivot+1] + regdiag[pivot+1] - offdiag[pivot]*offdiag[pivot]/diag_fac;
    }

    return 0;
}


trlib_int_t trlib_tri_timing_size() {
#if TRLIB_MEASURE_TIME
    return 1+TRLIB_SIZE_TIMING_LINALG+trlib_leftmost_timing_size()+trlib_eigen_timing_size();
#endif
    return 0;
}

trlib_int_t trlib_tri_factor_memory_size(trlib_int_t n) {
    return 6*n;
}
