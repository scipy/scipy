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

#include "trlib.h"
#include "trlib_private.h"

#include "_c99compat.h"

trlib_int_t trlib_leftmost(
        trlib_int_t nirblk, trlib_int_t *irblk, trlib_flt_t *diag, trlib_flt_t *offdiag,
        trlib_int_t warm, trlib_flt_t leftmost_minor, trlib_int_t itmax, trlib_flt_t tol_abs,
        trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_int_t *timing, trlib_int_t *ileftmost, trlib_flt_t *leftmost) {
    trlib_int_t ret = 0, curit = 0;
    if(! warm) {
        trlib_int_t curret = 0;
        trlib_int_t ii = 0;
        ret = 0;
        for(ii = 0; ii < nirblk; ++ii) {
            curret = trlib_leftmost_irreducible(irblk[ii+1]-irblk[ii], diag+irblk[ii], offdiag+irblk[ii], 0, 0.0, itmax,
                tol_abs, verbose, unicode, prefix, fout, timing, leftmost+ii, &curit);
            if (curret == 0) { ret = curret; }
        }
        *ileftmost = 0;
        for(ii = 1; ii < nirblk; ++ii) {
            if (leftmost[ii] < leftmost[*ileftmost]) { *ileftmost = ii; }
        }
    }
    else { 
        ret = trlib_leftmost_irreducible(irblk[nirblk] - irblk[nirblk-1], diag+irblk[nirblk-1], offdiag+irblk[nirblk-1],
                1, leftmost_minor, itmax, tol_abs, verbose, unicode, prefix, fout, timing, leftmost+nirblk-1, &curit);
        if (leftmost[nirblk-1] < leftmost[*ileftmost]) { *ileftmost = nirblk-1; }
    }
    return ret;
}

trlib_int_t trlib_leftmost_irreducible(
        trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
        trlib_int_t warm, trlib_flt_t leftmost_minor, trlib_int_t itmax, trlib_flt_t tol_abs,
        trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_int_t *timing, trlib_flt_t *leftmost, trlib_int_t *iter_pr) {
    // Local variables
    #if TRLIB_MEASURE_TIME
        struct timespec verystart, start, end;
        TRLIB_TIC(verystart)
    #endif
    trlib_int_t jj = 0;                     // local counter variable
    trlib_flt_t low = 0.0;                  // lower bracket variable: low <= leftmost       for desired value
    trlib_flt_t up = 0.0;                   // upper bracket variable:        leftmost <= up for desired value
    trlib_flt_t leftmost_attempt = 0.0;     // trial step for leftmost eigenvalue
    trlib_flt_t dleftmost = 0.0;            // increment
    trlib_flt_t prlp = 0.0;                 // value of Parlett-Reid-Last-Pivot function
    trlib_flt_t obyprlp = 0.0;              // quotient used in Cholesky computation
    trlib_flt_t dprlp = 0.0;                // derivative of Parlett-Reid-Last-Pivot function wrt to leftmost
    trlib_flt_t ddprlp = 0.0;               // second derivative of Parlett-Reid-Last-Pivot function wrt to leftmost
    trlib_int_t n_neg_piv = 0;              // number of negative pivots in factorization
    trlib_flt_t quad_abs = 0.0;             // absolute  coefficient in quadratic model
    trlib_flt_t quad_lin = 0.0;             // linear    coefficient in quadratic model
    trlib_flt_t quad_qua = 0.0;             // quadratic coefficient in quadratic model
    trlib_flt_t zerodum = 0.0;              // dummy return variables from quadratic equation
    trlib_flt_t oabs0 = 0.0, oabs1 = 0.0;   // temporaries in Gershgorin limit computation

    trlib_int_t continue_outer_loop = 0;    // local spaghetti code control variable
    trlib_int_t model_type = 0;
    trlib_int_t ii = 0;
    *leftmost = 0.0;                        // estimation of desired leftmost eigenvalue
    *iter_pr = 0;                           // iteration counter

    // trivial case: one-dimensional. return diagonal value
    if (n == 1) { *leftmost = diag[0]; TRLIB_RETURN(TRLIB_LMR_CONV) }

    /* set bracket interval derived from Gershgorin circles
       Gershgorin:
        eigenvalues are contained in the union of balls centered at
        diag_i with radius sum of absolute values in column i, except diagonal element
       this estimation is rough and could be improved by circle component analysis
              determine if worth doing */

    oabs0 = fabs(offdiag[0]); oabs1 = fabs(offdiag[n-2]);
    low = fmin( diag[0] - oabs0, diag[n-1] - oabs1 );
    up  = fmax( diag[0] + oabs0, diag[n-1] - oabs1 );
    for(ii = 1; ii < n-1; ++ii ) {
        oabs1 = fabs(offdiag[ii]);
        low = fmin( low, diag[ii] - oabs0 - oabs1 );
        up  = fmax( up,  diag[ii] + oabs0 + oabs1 );
        oabs0 = oabs1;
    }

    /* set leftmost to sensible initialization
       on warmstart, provided leftmost is eigenvalue of principal (n-1) * (n-1) submatrix
          by eigenvalue interlacing theorem desired value <= provided leftmost
       on coldstart, start close lower bound as hopefully this is a good estimation */
    if ( warm ) {
        // provided leftmost is an upper bound and a pole of Parlett-Reid Value, thus pertub a bit
        up = fmin(up, leftmost_minor); *leftmost = leftmost_minor - .1*(up-low); //*leftmost = leftmost_minor - TRLIB_EPS_POW_4;
    }  
    else { leftmost_minor = 0.0; *leftmost = low + .1*(up-low); }; // ensure sanity on leftmost_minor and start with lower bound
    // Parlett-Reid Iteration, note we can assume n > 1
    itmax = itmax*n;

    while (1) {
        /* iterate to obtain Parlett-Reid last pivot value of -leftmost == 0.0
           this iteration uses a safeguard bracket [low, up] such that always low <= leftmost <= up
           note that T - t*I is positive definite for t <= desired leftmost
           steps of iteration:
          
           (1) compute Parlett-Reid last pivot value which is D_n in a LDL^T factorization of T
               obtain derivative d D_n / d leftmost as byproduct in the very same recursion
               track if breakdown would occur in factorization, happens if either
               (a) a pivot become zero
               (b) more than one negative pivot present
               if breakdown would occurs this means that Parlett-Reid value is infinite
                 end iteration at premature point and restart with adapted bounds and estimation:
               (a) a pivot became zero:  
                   if last pivot zero   --> goal reached, exit
                   if previous zero     --> T - leftmost I not positive definite, thus desired value <= leftmost
               (b) multiple neg privots --> T - leftmost I            indefinite, thus desired value <= leftmost
           (2) compute a trial update for leftmost. Possibilities
               (a) Newton update
               (b) zero of suitable model of analytic expression,
                   analytic expression is given by prlp(t) = det(T-t*I)/det(U-t*I) with U principal (n-1)*(n-1) submatrix
                   prlp(t) has a pole at leftmost_minor, so better use lifted
                       lprlp(t) = (leftmost_minor-t)*prlp(t)
                   Gould proposes model    m(t) = a+bt+t^2  (correct asymptotic, very good for t far away)
                   Other choices  would be m(t) = a+bt      (Newton on lifted model)
                                           m(t) = a+bt+ct^2 (Taylor, very good close to zero, possibly really bad far away)
               do (b) if warmstart where user provided leftmost(U), otherwise go route (a)
          
           (3) take trial step if inside bracket, otherwise bisect
          
           stop iteration if either bracket is sufficiently small or Parlett-Reid value is close enough to zero */

        *iter_pr += 1;
        
        // test if iteration limit exceeded
        if ( *iter_pr > itmax ) { TRLIB_RETURN(TRLIB_LMR_ITMAX) }

        // initialize: no negative pivots so far
        n_neg_piv = 0;

        // print iteration headline every 10 iterations
        if (*iter_pr % 10 == 1) {
            TRLIB_PRINTLN_1("%6s%8s%14s%14s%14s%14s%14s%6s%6s", "  it  ", " action ", "     low      ", "   leftmost   ", "      up      ", "   dleftmost  ", "      prlp    ", " nneg ", "  br  ")
        }
        TRLIB_PRINTLN_1("%6ld%8s%14e%14e%14e", *iter_pr, "  entry ", low, *leftmost, up)

        // compute pivot and derivative of LDL^T factorization of T - leftmost I
        continue_outer_loop = 0;
        for( jj = 0; jj < n; ++jj ) {
            /* compute jj-th pivot
               special case for jj == 0 since offdiagonal is missing */
            if (jj == 0) { prlp = diag[0] - *leftmost; dprlp = -1.0; ddprlp = 0.0; }
            else{
                // update pivot as pivot      = d_j - leftmost - o_{j-1}^2/pivot
                // thus dpivot/dleftmost      =     - 1.0      + o_{j-1}^2/pivot^2 * dpivot
                //      d^2 pivot/dleftmost^2 =                 (o_{j-1}^2/pivot^2)(ddpivot - 2 dpivot^2/pivot)
                obyprlp = offdiag[jj-1]/prlp;
                //dprlp  = -1.0 + obyprlp*dprlp*obyprlp;// * dprlp;
                //prlp  = diag[jj] - offdiag[jj-1]*offdiag[jj-1]/prlp - *leftmost;
                ddprlp = obyprlp*obyprlp*(ddprlp - 2.0*dprlp*dprlp/prlp);
                dprlp  = -1.0 + offdiag[jj-1]*offdiag[jj-1]*dprlp / (prlp*prlp);
                prlp   = diag[jj] - *leftmost - offdiag[jj-1]*obyprlp;
            }

            // check for breakdown
            if (prlp == 0.0) {
                // if last pivot and no negative pivots encountered --> finished
                if (n_neg_piv == 0 && jj+1 == n) { TRLIB_RETURN(TRLIB_LMR_CONV) }
                else{
                    /* if not last pivot or negative pivots encountered:
                       estimation provides a new upper bound; reset estimation */
                    up = *leftmost;
                    *leftmost = 0.5 * (low+up);
                    continue_outer_loop = 1;
                    break; // continue outer loop
                }
            }
            else if ( prlp < 0.0 ) {
                n_neg_piv += 1;
                up = *leftmost;
                if (n_neg_piv > 1) {
                    // more than one negative pivot: factorization would fail, to the right of singularity!
                    *leftmost = 0.5 * (low+up);
                    continue_outer_loop = 1;
                    break; // continue outer loop
                }
            }
        }

        if (continue_outer_loop) { 
            TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14s%14e%6ld%6ld", "", " bisecp ", low, *leftmost, up, "", prlp, n_neg_piv, jj)
            continue; 
        }

        // we have survived computing the Last-Pivot value without finding a zero pivot and at most one negative pivot

        // adapt bracket, no negative pivots encountered: leftmost provides new lower bound, otherwise upper bound
        if (n_neg_piv == 0) { low = *leftmost; }
        else { up = *leftmost; }

        // test if bracket interval is small or last pivot has converged to zero
        if (up-low <= tol_abs * fmax(1.0, fmax(fabs(low), fabs(up))) || fabs(prlp) <= tol_abs) { 
            TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14s%14e%6ld%6ld", "", "  conv  ", low, *leftmost, up, "", prlp, n_neg_piv, jj)
            TRLIB_RETURN(TRLIB_LMR_CONV)
        }

        /* compute trial step for new leftmost
           on coldstart do Newton iteration, on warmstart find zero of model of analytic expression */
        // select suitable model for analytic expression in dependence of warm
        // warm: 1 --> heuristic depending on estimation if already close to zero
        // warm: 2 --> use asymptotic quadratic model of lifted prlp
        // warm: 3 --> use taylor     quadratic model of lifted prlp
        // warm: 4 --> use linear (newton)      model of lifted prlp
        if (warm) {
            if ( warm == 2 || (warm == 1 && up-low >= .1*fmax(1.0, fabs(*leftmost))) ) {
                /* use analytic model m(t) = (t-a)(t-b)/(t-leftmost_minor) for prlp(t)
                   fit a, b such that m matches function value and derivative
                   at current estimation and compute left zero of numerator */
                quad_lin = -(2.0*(*leftmost)+prlp+((*leftmost)-leftmost_minor)*dprlp);
                quad_abs = -(((*leftmost)-leftmost_minor)*prlp+(*leftmost)*(quad_lin+(*leftmost)));
                trlib_quadratic_zero(quad_abs, quad_lin, TRLIB_EPS_POW_75, 0, 0, "", NULL, &leftmost_attempt, &zerodum);
                model_type = 2; dleftmost = leftmost_attempt - *leftmost;
            }
            if( warm > 2 || (warm == 1 && up-low <.1*fmax(1.0, fabs(*leftmost))) ) {
                /* use quadratic taylor model for pole lifted function (leftmost_minor-t)*prlp(t) */
                quad_qua = -(dprlp + .5*((*leftmost)-leftmost_minor)*ddprlp);
                quad_lin = -(prlp + ((*leftmost)-leftmost_minor)*dprlp);
                quad_abs = (leftmost_minor-(*leftmost))*prlp;
                if ( warm == 4 || fabs(quad_qua) <= TRLIB_EPS * TRLIB_EPS || quad_lin*quad_lin - 4.0*quad_qua*quad_abs <= 0.0 ) { // tiny curvature, resort to Newton step of lifted function
                    model_type = 3; dleftmost = -quad_abs/quad_lin; leftmost_attempt = *leftmost + dleftmost;
                }
                else { // scale coefficients to normalize equation
                    quad_lin = quad_lin/quad_qua;
                    quad_abs = quad_abs/quad_qua;
                    trlib_quadratic_zero(quad_abs, quad_lin, TRLIB_EPS_POW_75, 0, 0, "", NULL, &dleftmost, &zerodum);
                    model_type = 4; leftmost_attempt = *leftmost + dleftmost;
                }
            }
        }
        else { model_type = 1; dleftmost = -prlp/dprlp; leftmost_attempt = *leftmost + dleftmost; } // Newton step

        // assess if we can use trial step
        if (low <= leftmost_attempt && leftmost_attempt <= up) { 
            if( fabs(dleftmost) <= tol_abs * fmax(1.0, fmax(fabs(low), fabs(up))) ) { TRLIB_RETURN(TRLIB_LMR_NEWTON_BREAK) }
            *leftmost = leftmost_attempt;
        }
        else { 
            // if warmstart information available, lifted newton step may still be feasible
            if(warm) {
                quad_lin = -(prlp + ((*leftmost)-leftmost_minor)*dprlp);
                quad_abs = (leftmost_minor-(*leftmost))*prlp;
                model_type = 3; dleftmost = -quad_abs/quad_lin; leftmost_attempt = *leftmost + dleftmost;
            }
            // if that fails, newton step may still be feasible
            if(low > leftmost_attempt || up < leftmost_attempt) {
                model_type = 1; dleftmost = -prlp/dprlp; leftmost_attempt = *leftmost + dleftmost;
            }
            // now out of options, do bisection
            if(low > leftmost_attempt || up < leftmost_attempt) {
                model_type = 0; *leftmost = .5*(low+up);
            }
        }
        if ( verbose > 0 ) {
            if ( model_type == 0 ) { TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14e%14e%6ld%6ld", "", " bisecs ", low, *leftmost, up, .5*(up-low), prlp, n_neg_piv, jj) }
            if ( model_type == 1 ) { TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14e%14e%6ld%6ld", "", "  piv 1 ", low, *leftmost, up, dleftmost, prlp, n_neg_piv, jj) }
            if ( model_type == 2 ) { TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14e%14e%6ld%6ld", "", " lpiv q ", low, *leftmost, up, dleftmost, prlp, n_neg_piv, jj) }
            if ( model_type == 3 ) { TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14e%14e%6ld%6ld", "", " lpiv 1 ", low, *leftmost, up, dleftmost, prlp, n_neg_piv, jj) }
            if ( model_type == 4 ) { TRLIB_PRINTLN_1("%6s%8s%14e%14e%14e%14e%14e%6ld%6ld", "", " lpiv 2 ", low, *leftmost, up, dleftmost, prlp, n_neg_piv, jj) }
        }

    }
}

trlib_int_t trlib_leftmost_timing_size() {
#if TRLIB_MEASURE_TIME
    return 1;
#endif
    return 0;
}

