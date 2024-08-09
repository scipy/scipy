/*
 * Copyright (C) 2024 SciPy developers
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Names of the SciPy Developers may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *
 * This file is a C translation of the Fortran code known as QUADPACK written by
 * R. Piessens, E. de Doncker-Kapenga, C. Überhuber , and D. Kahaner with the
 * original description below.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * QUADPACK is a FORTRAN subroutine package for the numerical
 * computation of definite one-dimensional integrals. It originated
 * from a joint project of R. Piessens and E. de Doncker (Appl.
 * Math. and Progr. Div.- K.U.Leuven, Belgium), C. Ueberhuber (Inst.
 * Fuer Math.- Techn.U.Wien, Austria), and D. Kahaner (Nation. Bur.
 * of Standards- Washington D.C., U.S.A.).
 * The routine names for the DOUBLE PRECISION versions are preceded
 * by the letter D.
 *
 * - QNG  : Is a simple non-adaptive automatic integrator, based on
 *          a sequence of rules with increasing degree of algebraic
 *          precision (Patterson, 1968).
 *
 * - QAG  : Is a simple globally adaptive integrator using the
 *          strategy of Aind (Piessens, 1973). It is possible to
 *          choose between 6 pairs of Gauss-Kronrod quadrature
 *          formulae for the rule evaluation component. The pairs
 *          of high degree of precision are suitable for handling
 *          integration difficulties due to a strongly oscillating
 *          integrand.
 *
 * - QAGS : Is an integrator based on globally adaptive interval
 *          subdivision in connection with extrapolation (de Doncker,
 *          1978) by the Epsilon algorithm (Wynn, 1956).
 *
 * - QAGP : Serves the same purposes as QAGS, but also allows
 *          for eventual user-supplied information, i.e. the
 *          abscissae of internal singularities, discontinuities
 *          and other difficulties of the integrand function.
 *          The algorithm is a modification of that in QAGS.
 *
 * - QAGI : Handles integration over infinite intervals. The
 *          infinite range is mapped onto a finite interval and
 *          then the same strategy as in QAGS is applied.
 *
 * - QAWO : Is a routine for the integration of COS(OMEGA*X)*F(X)
 *          or SIN(OMEGA*X)*F(X) over a finite interval (A,B).
 *          OMEGA is specified by the user
 *          The rule evaluation component is based on the
 *          modified Clenshaw-Curtis technique.
 *          An adaptive subdivision scheme is used connected with
 *          an extrapolation procedure, which is a modification
 *          of that in QAGS and provides the possibility to deal
 *          even with singularities in F.
 *
 * - QAWF : Calculates the Fourier cosine or Fourier sine
 *          transform of F(X), for user-supplied interval (A,
 *          INFINITY), OMEGA, and F. The procedure of QAWO is
 *          used on successive finite intervals, and convergence
 *          acceleration by means of the Epsilon algorithm (Wynn,
 *          1956) is applied to the series of the integral
 *          contributions.
 *
 * - QAWS : Integrates W(X)*F(X) over (A,B) with A.LT.B finite,
 *          and   W(X) = ((X-A)**ALFA)*((B-X)**BETA)*V(X)
 *          where V(X) = 1 or LOG(X-A) or LOG(B-X)
 *                         or LOG(X-A)*LOG(B-X)
 *          and   ALFA.GT.(-1), BETA.GT.(-1).
 *          The user specifies A, B, ALFA, BETA and the type of
 *          the function V.
 *          A globally adaptive subdivision strategy is applied,
 *          with modified Clenshaw-Curtis integration on the
 *          subintervals which contain A or B.
 *
 * - QAWC : Computes the Cauchy Principal Value of F(X)/(X-C)
 *          over a finite interval (A,B) and for
 *          user-determined C.
 *          The strategy is globally adaptive, and modified
 *          Clenshaw-Curtis integration is used on the subranges
 *          which contain the point X = C.
 *
 *    Each of the routines above also has a "more detailed" version
 * with a name ending in E, as QAGE.  These provide more
 * information and control than the easier versions.
 *
 *
 *    The preceding routines are all automatic.  That is, the user
 * inputs his problem and an error tolerance.  The routine
 * attempts to perform the integration to within the requested
 * absolute or relative error.
 *    There are, in addition, a number of non-automatic integrators.
 * These are most useful when the problem is such that the
 * user knows that a fixed rule will provide the accuracy
 * required.  Typically they return an error estimate but make
 * no attempt to satisfy any particular input error request.
 *
 *   QK15 QK21 QK31 QK41 QK51 QK61
 *        Estimate the integral on [a,b] using 15, 21,..., 61
 *        point rule and return an error estimate.
 *   QK15I 15 point rule for (semi)infinite interval.
 *   QK15W 15 point rule for special singular weight functions.
 *   QC25C 25 point rule for Cauchy Principal Values
 *   QC25F 25 point rule for sin/cos integrand.
 *   QMOMO Integrates k-th degree Chebychev polynomial times
 *         function with various explicit singularities.
 *
 * Support functions from linpack, slatec, and blas have been omitted
 * by default but can be obtained by asking.  For example, suppose you
 * already have installed linpack and the blas, but not slatec.  Then
 * use a request like  "send dqag from quadpack slatec".
 *
 *
 * [see also toms/691]
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * The original Fortran code can be found at https://www.netlib.org/quadpack/
 *
 * References:
 *
 * [1]: Robert Piessens, Elise Doncker-Kapenga, Christoph Überhuber, David
 *      Kahaner. "QUADPACK: A subroutine package for automatic integration",
 *      Springer, 1983, https://doi.org/10.1007/978-3-642-61786-7
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "__quadpack.h"
#include <math.h>


typedef double quadpack_w_func(const double,const double,const double,const double,const double,const int);

// Internal functions
static void dqc25c(double(*)(double*),const double,const double,const double,double*,double*,int*,int*);
static void dqc25f(double(*)(double*),const double,const double,const double,
                   const int,const int,const int, const int, double*,double*,
                   int*,double*,double*,int*,double*);
static void dqc25s(double(*)(double*),const double,const double,const double,
                   const double,const double,const double,double*,double*,double*,
                   double*,double*,double*,double*,const int,int*);
static void dqcheb(const double*,double*,double*,double*);
static void dqelg(int*,double*,double*,double*,double*,int*);
static void dqk15(double(*)(double*),const double,const double,double*,double*,double*,double*);
static void dqk15i(double(*)(double*),const double,const int,const double,const double,double*,double*,double*,double*);
static void dqk15w(double(*)(double*),double(),const double,const double,
                   const double,const double,const int,const double,const double,
                   double*,double*,double*,double*);
static void dqk21(double(*)(double*),const double,const double,double*,double*,double*,double*);
static void dqk31(double(*)(double*),const double,const double,double*,double*,double*,double*);
static void dqk41(double(*)(double*),const double,const double,double*,double*,double*,double*);
static void dqk51(double(*)(double*),const double,const double,double*,double*,double*,double*);
static void dqk61(double(*)(double*),const double,const double,double*,double*,double*,double*);
static void dqmomo(const double,const double,double*,double*,double*,double*,const int);
static void dqng(double(*)(double*),const double,const double,const double,const double,double*,double*,int*,int*);
static void dqpsrt(const int,const int,int*,double*,const double*,int *,int*);
static double dqwgtc(const double,const double,const double,const double,const double,const int);
static double dqwgtf(const double,const double,const double,const double,const double,const int);
static double dqwgts(const double,const double,const double,const double,const double,const int);

// Constants
static const double uflow = 2.2250738585072014e-308;  /* np.finfo(np.float64).tiny */
static const double oflow = 1.7976931348623157e+308;  /* np.finfo(np.float64).max  */
static const double epmach = 2.220446049250313e-016;  /* np.finfo(np.float64).eps  */


// Exported functions DQAGIE, DQAGPE, DQAGSE, DQAWCE, DQAWFE, DQAWOE, DQAWSE
void
dqagie(double(*fcn)(double* x), const double bound, const int inf,
       const double epsabs, const double epsrel, const int limit, double* result,
       double* abserr, int* neval, int* ier, double* alist, double* blist,
       double* rlist, double* elist, int* iord, int* last)
{
    // ***begin prologue  dqagie
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a3a1,h2a4a1
    // ***keywords  automatic integrator, infinite intervals,
    //              general-purpose, transformation, extrapolation,
    //              globally adaptive
    // ***author  piessens,robert,appl. math & progr. div - k.u.leuven
    //            de doncker,elise,appl. math & progr. div - k.u.leuven
    // ***purpose  the routine calculates an approximation result to a given
    //             integral   i = integral of f over (bound,+infinity)
    //             or i = integral of f over (-infinity,bound)
    //             or i = integral of f over (-infinity,+infinity),
    //             hopefully satisfying following claim for accuracy
    //             abs(i-result).le.max(epsabs,epsrel*abs(i))
    // ***description
    //
    //  integration over infinite intervals
    //  standard fortran subroutine
    //
    //             f      - double precision
    //                      function subprogram defining the integrand
    //                      function f(x). the actual name for f needs to be
    //                      declared e x t e r n a l in the driver program.
    //
    //             bound  - double precision
    //                      finite bound of integration range
    //                      (has no meaning if interval is doubly-infinite)
    //
    //             inf    - double precision
    //                      indicating the kind of integration range involved
    //                      inf = 1 corresponds to  (bound,+infinity),
    //                      inf = -1            to  (-infinity,bound),
    //                      inf = 2             to (-infinity,+infinity).
    //
    //             epsabs - double precision
    //                      absolute accuracy requested
    //             epsrel - double precision
    //                      relative accuracy requested
    //                      if  epsabs.le.0
    //                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                      the routine will end with ier = 6.
    //
    //             limit  - integer
    //                      gives an upper bound on the number of subintervals
    //                      in the partition of (a,b), limit.ge.1
    //
    //          on return
    //             result - double precision
    //                      approximation to the integral
    //
    //             abserr - double precision
    //                      estimate of the modulus of the absolute error,
    //                      which should equal or exceed abs(i-result)
    //
    //             neval  - integer
    //                      number of integrand evaluations
    //
    //             ier    - integer
    //                      ier = 0 normal and reliable termination of the
    //                              routine. it is assumed that the requested
    //                              accuracy has been achieved.
    //                    - ier.gt.0 abnormal termination of the routine. the
    //                              estimates for result and error are less
    //                              reliable. it is assumed that the requested
    //                              accuracy has not been achieved.
    //             error messages
    //                      ier = 1 maximum number of subdivisions allowed
    //                              has been achieved. one can allow more
    //                              subdivisions by increasing the value of
    //                              limit (and taking the according dimension
    //                              adjustments into account). however,if
    //                              this yields no improvement it is advised
    //                              to analyze the integrand in order to
    //                              determine the integration difficulties.
    //                              if the position of a local difficulty can
    //                              be determined (e.g. singularity,
    //                              discontinuity within the interval) one
    //                              will probably gain from splitting up the
    //                              interval at this point and calling the
    //                              integrator on the subranges. if possible,
    //                              an appropriate special-purpose integrator
    //                              should be used, which is designed for
    //                              handling the type of difficulty involved.
    //                          = 2 the occurrence of roundoff error is
    //                              detected, which prevents the requested
    //                              tolerance from being achieved.
    //                              the error may be under-estimated.
    //                          = 3 extremely bad integrand behaviour occurs
    //                              at some points of the integration
    //                              interval.
    //                          = 4 the algorithm does not converge.
    //                              roundoff error is detected in the
    //                              extrapolation table.
    //                              it is assumed that the requested tolerance
    //                              cannot be achieved, and that the returned
    //                              result is the best which can be obtained.
    //                          = 5 the integral is probably divergent, or
    //                              slowly convergent. it must be noted that
    //                              divergence can occur with any other value
    //                              of ier.
    //                          = 6 the input is invalid, because
    //                              (epsabs.le.0 and
    //                               epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                              result, abserr, neval, last, rlist(1),
    //                              elist(1) and iord(1) are set to zero.
    //                              alist(1) and blist(1) are set to 0
    //                              and 1 respectively.
    //
    //             alist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the left
    //                      end points of the subintervals in the partition
    //                      of the transformed integration range (0,1).
    //
    //             blist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the right
    //                      end points of the subintervals in the partition
    //                      of the transformed integration range (0,1).
    //
    //             rlist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the integral
    //                      approximations on the subintervals
    //
    //             elist  - double precision
    //                      vector of dimension at least limit,  the first
    //                      last elements of which are the moduli of the
    //                      absolute error estimates on the subintervals
    //
    //             iord   - integer
    //                      vector of dimension limit, the first k
    //                      elements of which are pointers to the
    //                      error estimates over the subintervals,
    //                      such that elist(iord(1)), ..., elist(iord(k))
    //                      form a decreasing sequence, with k = last
    //                      if last.le.(limit/2+2), and k = limit+1-last
    //                      otherwise
    //
    //             last   - integer
    //                      number of subintervals actually produced
    //                      in the subdivision process
    //
    // ***references  (none)
    // ***routines called  d1mach,dqelg,dqk15i,dqpsrt
    // ***end prologue  dqagie
    //
    //             the dimension of rlist2 is determined by the value of
    //             limexp in subroutine dqelg.
    //
    //
    //             list of major variables
    //             -----------------------
    //
    //            alist     - list of left end points of all subintervals
    //                        considered up to now
    //            blist     - list of right end points of all subintervals
    //                        considered up to now
    //            rlist(i)  - approximation to the integral over
    //                        (alist(i),blist(i))
    //            rlist2    - array of dimension at least (limexp+2),
    //                        containing the part of the epsilon table
    //                        which is still needed for further computations
    //            elist(i)  - error estimate applying to rlist(i)
    //            maxerr    - pointer to the interval with largest error
    //                        estimate
    //            errmax    - elist(maxerr)
    //            erlast    - error on the interval currently subdivided
    //                        (before that subdivision has taken place)
    //            area      - sum of the integrals over the subintervals
    //            errsum    - sum of the errors over the subintervals
    //            errbnd    - requested accuracy max(epsabs,epsrel*
    //                        abs(result))
    //            *****1    - variable for the left subinterval
    //            *****2    - variable for the right subinterval
    //            last      - index for subdivision
    //            nres      - number of calls to the extrapolation routine
    //            numrl2    - number of elements currently in rlist2. if an
    //                        appropriate approximation to the compounded
    //                        integral has been obtained, it is put in
    //                        rlist2(numrl2) after numrl2 has been increased
    //                        by one.
    //            small     - length of the smallest interval considered up
    //                        to now, multiplied by 1.5
    //            erlarg    - sum of the errors over the intervals larger
    //                        than the smallest interval considered up to now
    //            extrap    - logical variable denoting that the routine
    //                        is attempting to perform extrapolation. i.e.
    //                        before subdividing the smallest interval we
    //                        try to decrease the value of erlarg.
    //            noext     - logical variable denoting that extrapolation
    //                        is no longer allowed (true-value)
    //
    //             machine dependent constants
    //             ---------------------------
    //
    //            epmach is the largest relative spacing.
    //            uflow is the smallest positive magnitude.
    //            oflow is the largest positive magnitude.
    //
    int ierror, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, L, maxerr;
    int nres, nrmax, numrl2, extrap, noext;

    double a1, a2, abseps, area, area1, area12, area2, b1, b2, boun, correc, defabs;
    double defab1, defab2, dres, erlarg, erlast, errbnd, errmax, error1, error2;
    double error12, errsum, ertest, resabs, reseps, small;
    double rlist2[52], res3la[3];
    small = 0.0;
    correc = 0.0;

    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;
    alist[0] = 0.0;
    blist[0] = 1.0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    if ((epsabs <= 0.0) && (epsrel < fmax(50.0*epmach, 0.5e-28))) { *ier = 6; }
    if (*ier == 6) { return; }

    // first approximation to the integral
    // -----------------------------------
    // determine the interval to be mapped onto (0,1).
    // if inf = 2 the integral is computed as i = i1+i2, where
    // i1 = integral of f over (-infinity,0),
    // i2 = integral of f over (0,+infinity).

    boun = bound;
    if (inf == 2) { boun = 0.0; }
    dqk15i(fcn, boun, inf, 0.0, 1.0, result, abserr, &defabs, &resabs);

    // Test on accuracy.
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 0;
    dres = fabs(*result);
    errbnd = fmax(epsabs, epsrel*dres);
    if ((*abserr <= 100.0 * epmach * defabs) && (*abserr > errbnd)) { *ier = 2; }
    if (limit == 1) { *ier = 1; }
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) || (*abserr == 0.0)) { goto LINE130; }

    // Initialization for main loop.
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 0;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    nres = 0;
    ktmin = 0;
    numrl2 = 1;
    extrap = 0;
    noext = 0;
    ierror = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = (dres > (1.0 - 50.0 * epmach) * defabs ? 1 : -1);

    // Main for-loop.
    for (L = 1; L < limit; L++)
    {
        *last = L + 1;

        // Bisect the subinterval with nrmax-th largest error estimate.
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        dqk15i(fcn, boun, inf, a1, b1, &area1, &error1, &resabs, &defab1);
        dqk15i(fcn, boun, inf, a2, b2, &area2, &error2, &resabs, &defab2);

        // Improve previous approximations to integral and error and test for accuracy.
        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum = errsum + error12 - errmax;
        area = area + area12 - rlist[maxerr];

        if ((defab1 != error1) && (defab2 != error2))
        {
            if (!((fabs(rlist[maxerr] - area12) > 1.0e-5*fabs(area12)) || (error12 < 0.99*errmax)))
            {
                if (extrap) { iroff2++; } else { iroff1++; }
            }
            // 10
            if ((L > 9) && (error12 > errmax)){ iroff3++; }
        }
        // 15

        rlist[maxerr] = area1;
        rlist[L] = area2;
        errbnd = fmax(epsabs, epsrel*fabs(area));

        // Test for roundoff error and eventually set error flag.
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20)) { *ier = 2; }
        if (iroff2 >= 5) { ierror = 3; }

        // Set error flag in the case that the number of subintervals equals limit.
        if (*last == limit) { *ier = 1; }

        // Set error flag in the case of bad integrand behavior at some points
        // in the integration range.
        if (fmax(fabs(a1), fabs(b2)) <= (1.0 + 100.0*epmach)*(fabs(a2) + 1000.0*uflow)) { *ier = 4; }

        // Append the newly-created intervals to the list.
        if (!(error2 > error1))
        {
            alist[L] = a2;
            blist[maxerr] = b1;
            blist[L] = b2;
            elist[maxerr] = error1;
            elist[L] = error2;
        } else {
            alist[maxerr] = a2;
            alist[L] = a1;
            blist[L] = b1;
            rlist[maxerr] = area2;
            rlist[L] = area1;
            elist[maxerr] = error2;
            elist[L] = error1;
        }
        // 20

        // Call dqpsrt to maintain the descending ordering in the list of error
        // estimates and select the subinterval with nrmax-th largest error
        // estimate (to be bisected next).
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);

        if (errsum <= errbnd) { goto LINE115; }
        if (*ier != 0) { break; }
        if (L == 1) { goto LINE80; }
        if (noext) { continue; }
        erlarg = erlarg - erlast;
        if (fabs(b1 - a1) > small) { erlarg = erlarg + error12; }
        if (!(extrap))
        {
            // Test whether the interval to be bisected next is the smallest interval.
            if(blist[maxerr] - alist[maxerr] > small) { continue; }
            extrap = 1;
            nrmax = 1;
        }
        // 40

        if ((ierror == 3) || (erlarg <= ertest)) { goto LINE60; }
        // The smallest interval has the largest error. Before bisecting
        // decrease the sum of the errors over the larger intervals (erlarg)
        // and perform extrapolation.
        jupbnd = (*last > 2 + (limit/2) ? limit + 3 - *last : *last);
        for (k = nrmax; k < jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small) { goto LINE90; }
            nrmax++;
        }
        // 50
LINE60:
        // Perform extrapolation.
        numrl2++;
        rlist2[numrl2] = area;
        dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ktmin += 1;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) { *ier = 5; }
        if (abseps >= *abserr) { goto LINE70; }
        ktmin = 0;
        *abserr = abseps;
        *result = reseps;
        correc = erlarg;
        ertest = fmax(epsabs, epsrel*fabs(reseps));
        if (*abserr <= ertest) { break; }
LINE70:
        // Prepare bisection of the smallest interval.
        if (numrl2 == 0) { noext = 1; }
        if (*ier == 5) { break; }
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = 0;
        small = small * 0.5;
        erlarg = errsum;
        continue;
LINE80:
        small = 0.375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[1] = area;
LINE90:
        ; // no-op
    }

    // Set final result and error estimate.
    if (*abserr == oflow) { goto LINE115; }
    if ((*ier + ierror) == 0) { goto LINE110; }
    if (ierror == 3) { *abserr = *abserr + correc; }
    if (*ier == 0) { *ier = 3; }
    if ((*result != 0.0) && (area != 0.0)) { goto LINE105; }
    if (*abserr > errsum) { goto LINE115; }
    if (area == 0.0) { goto LINE130; }
    goto LINE110;
LINE105:
    if (*abserr/fabs(*result) > errsum/fabs(area)) { goto LINE115; }
LINE110:
    // Test on divergence.
    if ((ksgn == -1) && (fmax(fabs(*result),fabs(area)) <= defabs *0.01)) { goto LINE130; }
    if ((0.01 > *result/area) || (*result/area > 100.0) || (errsum > fabs(area))) { *ier = 6; }
    goto LINE130;
LINE115:
    // Compute global integral.
    *result = 0.0;
    for (k = 0; k <= L; k++) { *result = *result + rlist[k]; }
    *abserr = errsum;
LINE130:
    *neval = 30*(*last) - 15;
    if (inf == 2) { *neval *= 2; }
    if (*ier > 2) { *ier -= 1; }

    return;
}


void
dqagpe(double(*fcn)(double* x), const double a, const double b, int npts2,
       double* points, const double epsabs, const double epsrel, const int limit,
       double* result, double* abserr, int* neval, int* ier, double* alist,
       double* blist, double* rlist, double* elist, double* pts, int* iord,
       int* level, int* ndin, int* last)
{
    // ***begin prologue  dqagpe
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a1
    // ***keywords  automatic integrator, general-purpose,
    //              singularities at user specified points,
    //              extrapolation, globally adaptive.
    // ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  the routine calculates an approximation result to a given
    //             definite integral i = integral of f over (a,b), hopefully
    //             satisfying following claim for accuracy abs(i-result).le.
    //             max(epsabs,epsrel*abs(i)). break points of the integration
    //             interval, where local difficulties of the integrand may
    //             occur(e.g. singularities,discontinuities),provided by user.
    // ***description
    //
    //         computation of a definite integral
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //          on entry
    //             f      - double precision
    //                      function subprogram defining the integrand
    //                      function f(x). the actual name for f needs to be
    //                      declared e x t e r n a l in the driver program.
    //
    //             a      - double precision
    //                      lower limit of integration
    //
    //             b      - double precision
    //                      upper limit of integration
    //
    //             npts2  - integer
    //                      number equal to two more than the number of
    //                      user-supplied break points within the integration
    //                      range, npts2.ge.2.
    //                      if npts2.lt.2, the routine will end with ier = 6.
    //
    //             points - double precision
    //                      vector of dimension npts2, the first (npts2-2)
    //                      elements of which are the user provided break
    //                      points. if these points do not constitute an
    //                      ascending sequence there will be an automatic
    //                      sorting.
    //
    //             epsabs - double precision
    //                      absolute accuracy requested
    //             epsrel - double precision
    //                      relative accuracy requested
    //                      if  epsabs.le.0
    //                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                      the routine will end with ier = 6.
    //
    //             limit  - integer
    //                      gives an upper bound on the number of subintervals
    //                      in the partition of (a,b), limit.ge.npts2
    //                      if limit.lt.npts2, the routine will end with
    //                      ier = 6.
    //
    //          on return
    //             result - double precision
    //                      approximation to the integral
    //
    //             abserr - double precision
    //                      estimate of the modulus of the absolute error,
    //                      which should equal or exceed abs(i-result)
    //
    //             neval  - integer
    //                      number of integrand evaluations
    //
    //             ier    - integer
    //                      ier = 0 normal and reliable termination of the
    //                              routine. it is assumed that the requested
    //                              accuracy has been achieved.
    //                      ier.gt.0 abnormal termination of the routine.
    //                              the estimates for integral and error are
    //                              less reliable. it is assumed that the
    //                              requested accuracy has not been achieved.
    //             error messages
    //                      ier = 1 maximum number of subdivisions allowed
    //                              has been achieved. one can allow more
    //                              subdivisions by increasing the value of
    //                              limit (and taking the according dimension
    //                              adjustments into account). however, if
    //                              this yields no improvement it is advised
    //                              to analyze the integrand in order to
    //                              determine the integration difficulties. if
    //                              the position of a local difficulty can be
    //                              determined (i.e. singularity,
    //                              discontinuity within the interval), it
    //                              should be supplied to the routine as an
    //                              element of the vector points. if necessary
    //                              an appropriate special-purpose integrator
    //                              must be used, which is designed for
    //                              handling the type of difficulty involved.
    //                          = 2 the occurrence of roundoff error is
    //                              detected, which prevents the requested
    //                              tolerance from being achieved.
    //                              the error may be under-estimated.
    //                          = 3 extremely bad integrand behaviour occurs
    //                              at some points of the integration
    //                              interval.
    //                          = 4 the algorithm does not converge.
    //                              roundoff error is detected in the
    //                              extrapolation table. it is presumed that
    //                              the requested tolerance cannot be
    //                              achieved, and that the returned result is
    //                              the best which can be obtained.
    //                          = 5 the integral is probably divergent, or
    //                              slowly convergent. it must be noted that
    //                              divergence can occur with any other value
    //                              of ier.gt.0.
    //                          = 6 the input is invalid because
    //                              npts2.lt.2 or
    //                              break points are specified outside
    //                              the integration range or
    //                              (epsabs.le.0 and
    //                               epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
    //                              or limit.lt.npts2.
    //                              result, abserr, neval, last, rlist(1),
    //                              and elist(1) are set to zero. alist(1) and
    //                              blist(1) are set to a and b respectively.
    //
    //             alist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the left end points
    //                      of the subintervals in the partition of the given
    //                      integration range (a,b)
    //
    //             blist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the right end points
    //                      of the subintervals in the partition of the given
    //                      integration range (a,b)
    //
    //             rlist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the integral
    //                      approximations on the subintervals
    //
    //             elist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the moduli of the
    //                      absolute error estimates on the subintervals
    //
    //             pts    - double precision
    //                      vector of dimension at least npts2, containing the
    //                      integration limits and the break points of the
    //                      interval in ascending sequence.
    //
    //             level  - integer
    //                      vector of dimension at least limit, containing the
    //                      subdivision levels of the subinterval, i.e. if
    //                      (aa,bb) is a subinterval of (p1,p2) where p1 as
    //                      well as p2 is a user-provided break point or
    //                      integration limit, then (aa,bb) has level l if
    //                      abs(bb-aa) = abs(p2-p1)*2**(-l).
    //
    //             ndin   - integer
    //                      vector of dimension at least npts2, after first
    //                      integration over the intervals (pts(i)),pts(i+1),
    //                      i = 0,1, ..., npts2-2, the error estimates over
    //                      some of the intervals may have been increased
    //                      artificially, in order to put their subdivision
    //                      forward. if this happens for the subinterval
    //                      numbered k, ndin(k) is put to 1, otherwise
    //                      ndin(k) = 0.
    //
    //             iord   - integer
    //                      vector of dimension at least limit, the first k
    //                      elements of which are pointers to the
    //                      error estimates over the subintervals,
    //                      such that elist(iord(1)), ..., elist(iord(k))
    //                      form a decreasing sequence, with k = last
    //                      if last.le.(limit/2+2), and k = limit+1-last
    //                      otherwise
    //
    //             last   - integer
    //                      number of subintervals actually produced in the
    //                      subdivisions process
    //
    // ***references  (none)
    // ***routines called  d1mach,dqelg,dqk21,dqpsrt
    // ***end prologue  dqagpe
    //
    //             the dimension of rlist2 is determined by the value of
    //             limexp in subroutine epsalg (rlist2 should be of dimension
    //             (limexp+2) at least).
    //
    //
    //             list of major variables
    //             -----------------------
    //
    //            alist     - list of left end points of all subintervals
    //                        considered up to now
    //            blist     - list of right end points of all subintervals
    //                        considered up to now
    //            rlist(i)  - approximation to the integral over
    //                        (alist(i),blist(i))
    //            rlist2    - array of dimension at least limexp+2
    //                        containing the part of the epsilon table which
    //                        is still needed for further computations
    //            elist(i)  - error estimate applying to rlist(i)
    //            maxerr    - pointer to the interval with largest error
    //                        estimate
    //            errmax    - elist(maxerr)
    //            erlast    - error on the interval currently subdivided
    //                        (before that subdivision has taken place)
    //            area      - sum of the integrals over the subintervals
    //            errsum    - sum of the errors over the subintervals
    //            errbnd    - requested accuracy max(epsabs,epsrel*
    //                        abs(result))
    //            *****1    - variable for the left subinterval
    //            *****2    - variable for the right subinterval
    //            last      - index for subdivision
    //            nres      - number of calls to the extrapolation routine
    //            numrl2    - number of elements in rlist2. if an appropriate
    //                        approximation to the compounded integral has
    //                        been obtained, it is put in rlist2(numrl2) after
    //                        numrl2 has been increased by one.
    //            erlarg    - sum of the errors over the intervals larger
    //                        than the smallest interval considered up to now
    //            extrap    - logical variable denoting that the routine
    //                        is attempting to perform extrapolation. i.e.
    //                        before subdividing the smallest interval we
    //                        try to decrease the value of erlarg.
    //            noext     - logical variable denoting that extrapolation is
    //                        no longer allowed (true-value)
    //
    //             machine dependent constants
    //             ---------------------------
    //
    //            epmach is the largest relative spacing.
    //            uflow is the smallest positive magnitude.
    //            oflow is the largest positive magnitude.
    //
    int i, ierror, ind1, ind2, iroff1, iroff2, iroff3, j, jlow,jupbnd;
    int k, ksgn, ktmin, L, levcur, levmax, maxerr, nint, npts, nres, nrmax;
    int numrl2, extrap, noext;
    double abseps, area, area1, area12, area2, a1, a2, b1, b2, correc;
    double defabs, defab1, defab2, dres, erlarg, erlast, errbnd, errmax, error1;
    double error12, error2, errsum, ertest, resa, resabs, reseps, sign, temp;
    double res3la[3], rlist2[52];
    k = 0;
    correc = 0.0;

    // Test validity of parameters.
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    level[0] = 0;
    npts = npts2 - 2;
    if ((npts2 < 2) || (limit <= npts) ||
        ((epsabs <= 0.0) && (epsrel < fmax(50.0*epmach, 0.5e-28)))) { return; }
    *ier = 0;

    // If any break points are provided, sort them into an ascending sequence.
    sign = (a > b ? -1.0 : 1.0);
    pts[0] = fmin(a, b);
    if (npts != 0)
    {
        for (i = 0; i < npts; i++) { pts[i+1] = points[i]; }
    }
    // 15

    pts[npts + 1] = fmax(a, b);
    nint = npts + 1;
    a1 = pts[0];
    if (npts != 0)
    {
        for (i = 0; i < nint; i++)
        {
            for (j = i+1; j < nint+1; j++)
            {
                if (pts[i] <= pts[j]) { continue; }
                temp = pts[i];
                pts[i] = pts[j];
                pts[j] = temp;
            }
        }
        // 20
        if ((pts[0] != fmin(a, b)) || (pts[nint] != fmax(a, b)))
        {
            *ier = 6;
            return;
        }
    }
    // 40

    // Compute first integral and error approximations.
    resabs = 0.0;
    for (i = 0; i < nint; i++)
    {
        b1 = pts[i+1];
        dqk21(fcn, a1, b1, &area1, &error1, &defabs, &resa);
        *abserr = *abserr + error1;
        *result = *result + area1;
        ndin[i] = 0;
        if ((error1 == resa) && (error1 != 0.0)) { ndin[i] = 1; }
        resabs = resabs + defabs;
        level[i] = 0;
        elist[i] = error1;
        alist[i] = a1;
        blist[i] = b1;
        rlist[i] = area1;
        iord[i] = i;
        a1 = b1;
    }
    // 50

    errsum = 0.0;
    for (i = 0; i < nint; i++)
    {
        if (ndin[i] == 1) { elist[i] = *abserr; }
        errsum = errsum + elist[i];
    }
    // 55

    // Test on accuracy.
    *last = nint;
    *neval = 21 * nint;
    dres = fabs(*result);
    errbnd = fmax(epsabs, epsrel*dres);
    if ((*abserr <= 100.0 * epmach * resabs) && (*abserr > errbnd)) { *ier = 2; }
    if (nint != 1)
    {
        for (i = 0; i < npts; i++)
        {
            jlow = i+1;
            ind1 = iord[i];
            for (j = jlow; j < nint; j++)
            {
                ind2 = iord[j];
                if (elist[ind1] > elist[ind2]) { continue; }
                ind1 = ind2;
                k = j;
            }
            // 60

            if (ind1 == iord[i]) { continue; }
            iord[k] = iord[i];
            iord[i] = ind1;
        }
        // 70

        if (limit < npts2) { *ier = 1; }
    }
    // 80
    if ((*ier != 0) || (*abserr <= errbnd)) { goto LINE210; }

    // Initialization
    rlist2[0] = *result;
    maxerr = iord[0];
    errmax = elist[maxerr];
    area = *result;
    nrmax = 0;
    nres = 0;
    numrl2 = 0;
    ktmin = 0;
    extrap = 0;
    noext = 0;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierror = 0;
    *abserr = oflow;
    ksgn = (dres >= (1.0 - 50*epmach)*resabs ?  1 : -1);

    // Main for-loop.
    for (L = npts2 - 1; L < limit; L++)
    {
        *last = L + 1;

        // Bisect the subinterval with the nrmax-th largest error estimate.
        levcur = level[maxerr] + 1;
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        dqk21(fcn, a1, b1, &area1, &error1, &resa, &defab1);
        dqk21(fcn, a2, b2, &area2, &error2, &resa, &defab2);

        // Improve previous approximations to integral and error and test accuracy.
        *neval += 42;
        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum = errsum + error12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 != error1) && (defab2 != error2))
        {
            if (!((fabs(rlist[maxerr] - area12) > 1.0e-5*fabs(area12)) || (error12 < 0.99*errmax)))
            {
                if (extrap) { iroff2++; } else { iroff1++; }
            }
            // 90
            if ((L > 9) && error12 > errmax) { iroff3++; }
        }
        // 95
        level[maxerr] = levcur;
        level[L] = levcur;
        rlist[maxerr] = area1;
        rlist[L] = area2;
        errbnd = fmax(epsabs, epsrel*fabs(area));

        // Test for roundoff error and eventually set error flag.
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20)) { *ier = 2; }
        if (iroff2 >= 5) { ierror = 3; }

        // Set error flag in the case that the number of subintervals equals limit.
        if (*last == limit) { *ier = 1; }

        // Set error flag in the case of bad integrand behavior at a point of the
        // integration range.
        if (fmax(fabs(a1), fabs(b2)) <= (1.0 +100.0*epmach)*(fabs(a2) + 1000.0*uflow))
        {
            *ier = 4;
        }

        // Append the newly-created intervals to the list.
        if (!(error2 > error1))
        {
            alist[L] = a2;
            blist[maxerr] = b1;
            blist[L] = b2;
            elist[maxerr] = error1;
            elist[L] = error2;
        } else {
            alist[maxerr] = a2;
            alist[L] = a1;
            blist[L] = b1;
            rlist[maxerr] = area2;
            rlist[L] = area1;
            elist[maxerr] = error2;
            elist[L] = error1;
        }
        // 110

        // Call subroutine dqpsrt to maintain the descending ordering in the list
        // of error estimates and select the subinterval with nrmax-th largest
        // error estimate (to be bisected next).
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);

        if (errsum <= errbnd) { goto LINE190; }
        if (*ier != 0) { break; }
        if (noext) { continue; }
        erlarg = erlarg - erlast;
        if (levcur+1 <= levmax) { erlarg = erlarg + error12; }
        if (!(extrap))
        {
            if(level[maxerr]+1 <= levmax) { continue; }
            extrap = 1;
            nrmax = 1;
        }
        // 120

        if ((ierror != 3) && (!(erlarg <= ertest)))
        {
            // The smallest interval has the largest error. Before bisecting
            // decrease the sum of the errors over the larger intervals (erlarg)
            // and perform extrapolation.
            jupbnd = (*last > 2 + (limit/2) ? limit + 3 - *last : *last);
            for (k = nrmax; k < jupbnd; k++) {
                maxerr = iord[nrmax];
                errmax = elist[maxerr];
                if (level[maxerr] +1 <= levmax) { goto LINE160; }  // break->continue
                nrmax++;
            }
            // 130
        }
        // 140

        // Perform extrapolation.
        numrl2++;
        rlist2[numrl2] = area;
        if (numrl2 <= 1) { goto LINE155; }
        dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) { *ier = 5; }
        if (abseps >= *abserr) { goto LINE150; }
        ktmin = 0;
        *abserr = abseps;
        *result = reseps;
        correc = erlarg;
        ertest = fmax(epsabs, epsrel*fabs(reseps));
        if (*abserr < ertest) { break; }
LINE150:
        // Prepare bisection of the smallest interval.
        if (numrl2 == 0) noext = 1;
        if (*ier >= 5) { break; }
LINE155:
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = 0;
        levmax += 1;
        erlarg = errsum;
LINE160:
        ; // no-op.
    }

    // 170
    if (*abserr == oflow) { goto LINE190; }
    if ((*ier + ierror) == 0) { goto LINE180; }
    if (ierror == 3) *abserr = *abserr + correc;
    if (*ier == 0) { *ier = 3; }
    if ((*result != 0.0) && (area != 0.0)) { goto LINE175; }
    if (*abserr > errsum) { goto LINE190; }
    if (area == 0.0) { goto LINE210; }
    goto LINE180;
LINE175:
    if (*abserr/fabs(*result) > errsum/fabs(area)) { goto LINE190; }
LINE180:
    // Test on divergence.
    if ((ksgn == -1) && (fmax(fabs(*result), fabs(area)) <= resabs*0.01)) { goto LINE210; }
    if ((0.01 > *result/area) || (*result/area > 100.0) || (errsum > fabs(area))) { *ier = 6; }
    goto LINE210;
LINE190:
    // Compute global integral sum.
    *result = 0.0;
    for (k = 0; k <= L; k++) { *result = *result + rlist[k]; }
    *abserr = errsum;
LINE210:
    if (*ier > 2) { *ier -= 1; }
    *result = (*result)*sign;

    return;
}


void
dqagse(double(*fcn)(double* x), const double a, const double b,
       const double epsabs, const double epsrel, const int limit, double* result,
       double* abserr, int* neval, int* ier, double* alist, double* blist,
       double* rlist, double* elist, int* iord, int* last)
{
    // ***begin prologue  dqagse
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a1
    // ***keywords  automatic integrator, general-purpose,
    //              (end point) singularities, extrapolation,
    //              globally adaptive
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  the routine calculates an approximation result to a given
    //             definite integral i = integral of f over (a,b),
    //             hopefully satisfying following claim for accuracy
    //             abs(i-result).le.max(epsabs,epsrel*abs(i)).
    // ***description
    //
    //         computation of a definite integral
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //          on entry
    //             f      - double precision
    //                      function subprogram defining the integrand
    //                      function f(x). the actual name for f needs to be
    //                      declared e x t e r n a l in the driver program.
    //
    //             a      - double precision
    //                      lower limit of integration
    //
    //             b      - double precision
    //                      upper limit of integration
    //
    //             epsabs - double precision
    //                      absolute accuracy requested
    //             epsrel - double precision
    //                      relative accuracy requested
    //                      if  epsabs.le.0
    //                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                      the routine will end with ier = 6.
    //
    //             limit  - integer
    //                      gives an upperbound on the number of subintervals
    //                      in the partition of (a,b)
    //
    //          on return
    //             result - double precision
    //                      approximation to the integral
    //
    //             abserr - double precision
    //                      estimate of the modulus of the absolute error,
    //                      which should equal or exceed abs(i-result)
    //
    //             neval  - integer
    //                      number of integrand evaluations
    //
    //             ier    - integer
    //                      ier = 0 normal and reliable termination of the
    //                              routine. it is assumed that the requested
    //                              accuracy has been achieved.
    //                      ier.gt.0 abnormal termination of the routine
    //                              the estimates for integral and error are
    //                              less reliable. it is assumed that the
    //                              requested accuracy has not been achieved.
    //             error messages
    //                          = 1 maximum number of subdivisions allowed
    //                              has been achieved. one can allow more sub-
    //                              divisions by increasing the value of limit
    //                              (and taking the according dimension
    //                              adjustments into account). however, if
    //                              this yields no improvement it is advised
    //                              to analyze the integrand in order to
    //                              determine the integration difficulties. if
    //                              the position of a local difficulty can be
    //                              determined (e.g. singularity,
    //                              discontinuity within the interval) one
    //                              will probably gain from splitting up the
    //                              interval at this point and calling the
    //                              integrator on the subranges. if possible,
    //                              an appropriate special-purpose integrator
    //                              should be used, which is designed for
    //                              handling the type of difficulty involved.
    //                          = 2 the occurrence of roundoff error is detec-
    //                              ted, which prevents the requested
    //                              tolerance from being achieved.
    //                              the error may be under-estimated.
    //                          = 3 extremely bad integrand behaviour
    //                              occurs at some points of the integration
    //                              interval.
    //                          = 4 the algorithm does not converge.
    //                              roundoff error is detected in the
    //                              extrapolation table.
    //                              it is presumed that the requested
    //                              tolerance cannot be achieved, and that the
    //                              returned result is the best which can be
    //                              obtained.
    //                          = 5 the integral is probably divergent, or
    //                              slowly convergent. it must be noted that
    //                              divergence can occur with any other value
    //                              of ier.
    //                          = 6 the input is invalid, because
    //                              epsabs.le.0 and
    //                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
    //                              result, abserr, neval, last, rlist(1),
    //                              iord(1) and elist(1) are set to zero.
    //                              alist(1) and blist(1) are set to a and b
    //                              respectively.
    //
    //             alist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the left end points
    //                      of the subintervals in the partition of the
    //                      given integration range (a,b)
    //
    //             blist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the right end points
    //                      of the subintervals in the partition of the given
    //                      integration range (a,b)
    //
    //             rlist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the integral
    //                      approximations on the subintervals
    //
    //             elist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the moduli of the
    //                      absolute error estimates on the subintervals
    //
    //             iord   - integer
    //                      vector of dimension at least limit, the first k
    //                      elements of which are pointers to the
    //                      error estimates over the subintervals,
    //                      such that elist(iord(1)), ..., elist(iord(k))
    //                      form a decreasing sequence, with k = last
    //                      if last.le.(limit/2+2), and k = limit+1-last
    //                      otherwise
    //
    //             last   - integer
    //                      number of subintervals actually produced in the
    //                      subdivision process
    //
    // ***references  (none)
    // ***routines called  d1mach,dqelg,dqk21,dqpsrt
    // ***end prologue  dqagse
    //
    //
    //             the dimension of rlist2 is determined by the value of
    //             limexp in subroutine dqelg (rlist2 should be of dimension
    //             (limexp+2) at least).
    //
    //             list of major variables
    //             -----------------------
    //
    //            alist     - list of left end points of all subintervals
    //                        considered up to now
    //            blist     - list of right end points of all subintervals
    //                        considered up to now
    //            rlist(i)  - approximation to the integral over
    //                        (alist(i),blist(i))
    //            rlist2    - array of dimension at least limexp+2 containing
    //                        the part of the epsilon table which is still
    //                        needed for further computations
    //            elist(i)  - error estimate applying to rlist(i)
    //            maxerr    - pointer to the interval with largest error
    //                        estimate
    //            errmax    - elist(maxerr)
    //            erlast    - error on the interval currently subdivided
    //                        (before that subdivision has taken place)
    //            area      - sum of the integrals over the subintervals
    //            errsum    - sum of the errors over the subintervals
    //            errbnd    - requested accuracy max(epsabs,epsrel*
    //                        abs(result))
    //            *****1    - variable for the left interval
    //            *****2    - variable for the right interval
    //            last      - index for subdivision
    //            nres      - number of calls to the extrapolation routine
    //            numrl2    - number of elements currently in rlist2. if an
    //                        appropriate approximation to the compounded
    //                        integral has been obtained it is put in
    //                        rlist2(numrl2) after numrl2 has been increased
    //                        by one.
    //            small     - length of the smallest interval considered up
    //                        to now, multiplied by 1.5
    //            erlarg    - sum of the errors over the intervals larger
    //                        than the smallest interval considered up to now
    //            extrap    - logical variable denoting that the routine is
    //                        attempting to perform extrapolation i.e. before
    //                        subdividing the smallest interval we try to
    //                        decrease the value of erlarg.
    //            noext     - logical variable denoting that extrapolation
    //                        is no longer allowed (true value)
    //
    //             machine dependent constants
    //             ---------------------------
    //
    //            epmach is the largest relative spacing.
    //            uflow is the smallest positive magnitude.
    //            oflow is the largest positive magnitude.
    //
    int extrap, ierror, iroff1, iroff2, iroff3, k, ksgn, ktmin, L, maxerr;
    int noext, nres, nrmax, numrl2;
    double a1, a2, area, area1, area12, area2, b1, b2, defab1, defab2, defabs;
    double dres, erlarg, erlast, errbnd, errmax, error1, error12, error2, errsum;
    double ertest, small, jupbnd, correc, resabs, reseps, abseps;
    double res3la[3], rlist2[52];
    small = 0.0;

    *ier = 0;
    *neval = 0;
    *last = 0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    *result = 0.0;
    *abserr = 0.0;
    ierror = 0;
    erlarg = 0.0;
    correc = 0.0;

    if ((epsabs <= 0.0) && (epsrel < fmax(50.0*epmach, 0.5e-28)))
    {
        *ier = 6;
        return;
    }

    // First approximation to the integral.
    dqk21(fcn, a, b, result, abserr, &defabs, &resabs);

    // Test on accuracy
    dres = fabs(*result);
    errbnd = fmax(epsabs, epsrel*dres);
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 0;

    if ((*abserr <= 100.0*epmach*defabs) && (*abserr > errbnd)) { *ier = 2; }
    if (limit == 1) { *ier = 1; }
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) || (*abserr == 0.0)) { goto LINE140; }

    // Initialization.
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 0;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    nres = 0;
    numrl2 = 1;
    ktmin = 0;
    extrap = 0;
    noext = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = (dres >= (1.0 - 50.0*epmach)*defabs ? 1 : -1);

    // Main for-loop
    for (L = 1; L < limit; L++)
    {
        *last = L + 1;

        // Bisect the subinterval with nrmax-th largest error estimate.
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        dqk21(fcn, a1, b1, &area1, &error1, &resabs, &defab1);
        dqk21(fcn, a2, b2, &area2, &error2, &resabs, &defab2);

        // Improve previous approximations to integral and error and test for accuracy.
        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum = errsum + error12 - errmax;
        area = area + area12 - rlist[maxerr];

        if ((defab1 != error1) && (defab2 != error2))
        {
            if (!((fabs(rlist[maxerr] - area12) > 1.0e-5*fabs(area12)) || (error12 < 0.99*errmax)))
            {
                if (extrap) { iroff2++; } else { iroff1++; }
            }
            // 10
            if ((L > 9) && (error12 > errmax)) { iroff3++; }
        }
        // 15

        rlist[maxerr] = area1;
        rlist[L] = area2;
        errbnd = fmax(epsabs, epsrel*fabs(area));

        // Test for roundoff error and eventually set error flag.
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20)) { *ier = 2; }
        if (iroff2 >= 5) { ierror = 3; }

        // Set error flag in the case that the number of subintervals equals limit.
        if (*last == limit) { *ier = 1; }

        // Set error flag in the case of bad integrand behavior at some points
        // in the integration range.
        if (fmax(fabs(a1), fabs(b2)) <= (1.0 + 100.0*epmach)*(fabs(a2) + 1000.0*uflow))
        {
            *ier = 4;
        }

        // Append the newly-created intervals to the list.
        if (!(error2 > error1))
        {
            alist[L] = a2;
            blist[maxerr] = b1;
            blist[L] = b2;
            elist[maxerr] = error1;
            elist[L] = error2;
        } else {
            alist[maxerr] = a2;
            alist[L] = a1;
            blist[L] = b1;
            rlist[maxerr] = area2;
            rlist[L] = area1;
            elist[maxerr] = error2;
            elist[L] = error1;
        }

        // call subroutine dqpsrt to maintain the descending ordering in the
        // list of error estimates and select the subinterval with nrmax-th
        // largest error estimate (to be bisected next).
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);

        if (errsum <= errbnd) { goto LINE115; }
        if (*ier != 0) { break; }
        if (L == 1) { goto LINE80; }
        if (noext) { continue; }
        erlarg = erlarg - erlast;
        if (fabs(b1 - a1) > small) { erlarg = erlarg + error12; }
        if (!(extrap))
        {
            // Test whether the interval to be bisected next is the smallest interval.
            if (fabs(blist[maxerr] - alist[maxerr]) > small) { continue; }
            extrap = 1;
            nrmax = 1;
        }
        // 40
        if ((ierror == 3) || (erlarg <= ertest)) { goto LINE60; }

        // The smallest interval has the largest error. Before bisecting
        // decrease the sum of the errors over the larger intervals (erlarg)
        // and perform extrapolation.
        jupbnd = (*last > 2 + (limit/2) ? limit + 3 - *last : *last);
        for (k = nrmax; k < jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small) { goto LINE90; } // break->continue
            nrmax++;
        }
        // 50
LINE60:
        // Perform extrapolation.
        numrl2++;
        rlist2[numrl2] = area;
        dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ktmin += 1;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) { *ier = 5; }
        if (abseps >= *abserr) { goto LINE70; }
        ktmin = 0;
        *abserr = abseps;
        *result = reseps;
        correc = erlarg;
        ertest = fmax(epsabs, epsrel*fabs(reseps));
        if (*abserr <= ertest) { break; }
LINE70:
        // Prepare bisection of the smallest interval.
        if (numrl2 == 0) { noext = 1; }
        if (*ier == 5) { break; }
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = 0;
        small = small * 0.5;
        erlarg = errsum;
        continue;
LINE80:
        small = fabs(b - a)*0.375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[1] = area;
LINE90:
        ; // no-op
    }

    // Set final result and error estimate
    if (*abserr == oflow) { goto LINE115; }
    if ((*ier + ierror) == 0) { goto LINE110; }
    if (ierror == 3) { *abserr = *abserr + correc; }
    if (*ier == 0) { *ier = 3; }
    if ((*result != 0.0) && (area != 0.0)) { goto LINE105; }
    if (*abserr > errsum) { goto LINE115; }
    if (area == 0.0) { goto LINE130; }
    goto LINE110;
LINE105:
    if (*abserr/fabs(*result) > errsum/fabs(area)) { goto LINE115; }
LINE110:
    // Test on divergence.
    if ((ksgn == -1) && (fmax(fabs(*result), fabs(area)) <= defabs*0.01)) { goto LINE130; }
    if ((0.01 > (*result/area)) || ((*result/area) > 100.0) || (errsum > fabs(area))) { *ier = 6; }
    goto LINE130;

LINE115:
    // Compute global integral sum.
    *result = 0.0;
    for (k = 0; k <= L; k++)
    {
        *result = *result + rlist[k];
    }
    // 120
    *abserr = errsum;
LINE130:
    if (*ier > 2) { *ier -= 1; }
LINE140:
    *neval = 42*(*last) - 21;

    return;
}


void
dqawce(double(*fcn)(double* x), const double a, const double b, const double c,
       const double epsabs, const double epsrel, const int limit, double* result,
       double* abserr, int* neval, int* ier, double* alist, double* blist,
       double* rlist, double* elist, int* iord, int* last)
{
    // ***begin prologue  dqawce
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a1,j4
    // ***keywords  automatic integrator, special-purpose,
    //              cauchy principal value, clenshaw-curtis method
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***  purpose  the routine calculates an approximation result to a
    //               cauchy principal value i = integral of f*w over (a,b)
    //               (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
    //               following claim for accuracy
    //               abs(i-result).le.max(epsabs,epsrel*abs(i))
    // ***description
    //
    //         computation of a cauchy principal value
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //          on entry
    //             f      - double precision
    //                      function subprogram defining the integrand
    //                      function f(x). the actual name for f needs to be
    //                      declared e x t e r n a l in the driver program.
    //
    //             a      - double precision
    //                      lower limit of integration
    //
    //             b      - double precision
    //                      upper limit of integration
    //
    //             c      - double precision
    //                      parameter in the weight function, c.ne.a, c.ne.b
    //                      if c = a or c = b, the routine will end with
    //                      ier = 6.
    //
    //             epsabs - double precision
    //                      absolute accuracy requested
    //             epsrel - double precision
    //                      relative accuracy requested
    //                      if  epsabs.le.0
    //                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                      the routine will end with ier = 6.
    //
    //             limit  - integer
    //                      gives an upper bound on the number of subintervals
    //                      in the partition of (a,b), limit.ge.1
    //
    //          on return
    //             result - double precision
    //                      approximation to the integral
    //
    //             abserr - double precision
    //                      estimate of the modulus of the absolute error,
    //                      which should equal or exceed abs(i-result)
    //
    //             neval  - integer
    //                      number of integrand evaluations
    //
    //             ier    - integer
    //                      ier = 0 normal and reliable termination of the
    //                              routine. it is assumed that the requested
    //                              accuracy has been achieved.
    //                      ier.gt.0 abnormal termination of the routine
    //                              the estimates for integral and error are
    //                              less reliable. it is assumed that the
    //                              requested accuracy has not been achieved.
    //             error messages
    //                      ier = 1 maximum number of subdivisions allowed
    //                              has been achieved. one can allow more sub-
    //                              divisions by increasing the value of
    //                              limit. however, if this yields no
    //                              improvement it is advised to analyze the
    //                              the integrand, in order to determine the
    //                              the integration difficulties. if the
    //                              position of a local difficulty can be
    //                              determined (e.g. singularity,
    //                              discontinuity within the interval) one
    //                              will probably gain from splitting up the
    //                              interval at this point and calling
    //                              appropriate integrators on the subranges.
    //                          = 2 the occurrence of roundoff error is detec-
    //                              ted, which prevents the requested
    //                              tolerance from being achieved.
    //                          = 3 extremely bad integrand behaviour
    //                              occurs at some interior points of
    //                              the integration interval.
    //                          = 6 the input is invalid, because
    //                              c = a or c = b or
    //                              (epsabs.le.0 and
    //                               epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
    //                              or limit.lt.1.
    //                              result, abserr, neval, rlist(1), elist(1),
    //                              iord(1) and last are set to zero. alist(1)
    //                              and blist(1) are set to a and b
    //                              respectively.
    //
    //             alist   - double precision
    //                       vector of dimension at least limit, the first
    //                        last  elements of which are the left
    //                       end points of the subintervals in the partition
    //                       of the given integration range (a,b)
    //
    //             blist   - double precision
    //                       vector of dimension at least limit, the first
    //                        last  elements of which are the right
    //                       end points of the subintervals in the partition
    //                       of the given integration range (a,b)
    //
    //             rlist   - double precision
    //                       vector of dimension at least limit, the first
    //                        last  elements of which are the integral
    //                       approximations on the subintervals
    //
    //             elist   - double precision
    //                       vector of dimension limit, the first  last
    //                       elements of which are the moduli of the absolute
    //                       error estimates on the subintervals
    //
    //             iord    - integer
    //                       vector of dimension at least limit, the first k
    //                       elements of which are pointers to the error
    //                       estimates over the subintervals, so that
    //                       elist(iord(1)), ..., elist(iord(k)) with k = last
    //                       if last.le.(limit/2+2), and k = limit+1-last
    //                       otherwise, form a decreasing sequence
    //
    //             last    - integer
    //                       number of subintervals actually produced in
    //                       the subdivision process
    //
    // ***references  (none)
    // ***routines called  d1mach,dqc25c,dqpsrt
    // ***end prologue  dqawce
    //
    //             list of major variables
    //             -----------------------
    //
    //            alist     - list of left end points of all subintervals
    //                        considered up to now
    //            blist     - list of right end points of all subintervals
    //                        considered up to now
    //            rlist(i)  - approximation to the integral over
    //                        (alist(i),blist(i))
    //            elist(i)  - error estimate applying to rlist(i)
    //            maxerr    - pointer to the interval with largest
    //                        error estimate
    //            errmax    - elist(maxerr)
    //            area      - sum of the integrals over the subintervals
    //            errsum    - sum of the errors over the subintervals
    //            errbnd    - requested accuracy max(epsabs,epsrel*
    //                        abs(result))
    //            *****1    - variable for the left subinterval
    //            *****2    - variable for the right subinterval
    //            last      - index for subdivision
    //
    //
    //             machine dependent constants
    //             ---------------------------
    //
    //            epmach is the largest relative spacing.
    //            uflow is the smallest positive magnitude.
    //
    int iroff1, iroff2, k, krule, L, maxerr, nev, nrmax;
    double a1, a2, aa, area, area1, area12, area2, b1, b2, bb, errbnd, errmax;
    double error1, error12, error2, errsum;

    *ier = 6;
    *neval = 0;
    *last = 0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    *result = 0.0;
    *abserr = 0.0;

    if ((c == a) || (c == b) ||
        ((epsabs <= 0.0) && (epsrel < fmax(50.0*epmach, 0.5e-28)))) { return; }
    *ier = 0;
    // First approximation to the integral
    aa = (a > b ? b : a);
    bb = (a > b ? a : b);

    krule = 1;
    dqc25c(fcn, aa, bb, c, result, abserr, &krule, neval);
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 0;
    alist[0] = a;
    blist[0] = b;

    // Test on accuracy.
    errbnd = fmax(epsabs, epsrel * fabs(*result));
    if (limit == 1) { *ier = 1; }
    if ((*abserr < fmin(0.01*fabs(*result), errbnd)) || (*ier == 1))
    {
        // 70
        if (aa == b) { *result = -(*result); }
        return;
    }

    // Initialization
    alist[0] = aa;
    blist[0] = bb;
    rlist[0] = *result;
    errmax = *abserr;
    maxerr = 0;
    area = *result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

    // Main for-loop.
    for (L = 1; L < limit; L++)
    {
        *last = L + 1;

        // Bisect the subinterval with nrmax-th largest error estimate.
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr] + blist[maxerr]);
        b2 = blist[maxerr];
        if ((c <= b1) && (c > a1)) { b1 = 0.5 * (c + b2); }
        if ((c >  b1) && (c < b2)) { b1 = 0.5 * (a1 + c); }
        a2 = b1;
        krule = 2;
        dqc25c(fcn, a1, b1, c, &area1, &error1, &krule, &nev);
        *neval += nev;
        dqc25c(fcn, a2, b2, c, &area2, &error2, &krule, &nev);
        *neval += nev;

        // Improve previous approximations to integral and error and test for accuracy.
        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum = errsum + error12 - errmax;
        area = area + area12 - rlist[maxerr];

        if ((fabs(rlist[maxerr]-area12) < (1.0e-5*fabs(area12))) &&
            (error12 >= 0.99 * errmax) && (krule == 0)) { iroff1++; }
        if ((L > 9) && (error12 > errmax) && (krule == 0)) { iroff2++; }
        rlist[maxerr] = area1;
        rlist[L] = area2;
        errbnd = fmax(epsabs, epsrel*fabs(area));
        if (!(errsum <= errbnd))
        {
            // Test for roundoff error and eventually set error flag.
            if ((iroff1 >= 6) && (iroff2 > 20)) { *ier = 2; }

            // Set error flag in the case that number of interval bisections exceeds limit.
            if (*last == limit) { *ier = 1; }

            // Set error flag in the case of bad integrand behavior at a point
            // of the integration range.
            if (fmax(fabs(a1), fabs(b2)) <= (1.0 + 100.0*epmach)*(fabs(a2) + 1000.0*uflow)) { *ier = 3; }
        }
        // Append the newly-created intervals to the list.
        if (!(error2 > error1))
        {
            alist[L] = a2;
            blist[maxerr] = b1;
            blist[L] = b2;
            elist[maxerr] = error1;
            elist[L] = error2;
        } else {
            alist[maxerr] = a2;
            alist[L] = a1;
            blist[L] = b1;
            rlist[maxerr] = area2;
            rlist[L] = area1;
            elist[maxerr] = error2;
            elist[L] = error1;
        }
        // 30

        // Call subroutine dqpsrt to maintain the descending ordering in the
        // list of error estimates and select the subinterval with nrmax-th
        // largest error estimate (to be bisected next).
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);
        if( (*ier != 0) || (errsum <= errbnd)) { break; }
    }
    // 50

    *result = 0.0;
    for (k = 0; k <= L; k++) {
        *result = *result + rlist[k];
    }
    // 60

    *abserr = errsum;
    // 70
    if (aa == b) { *result = -*result; }

    return;
}


void
dqawfe(double(*fcn)(double* x), const double a, const double omega, const int integr,
       const double epsabs, const int limlst, const int limit, const int maxp1,
       double* result, double* abserr, int* neval, int* ier, double* rslst,
       double* erlst, int* ierlst, int* lst, double* alist, double* blist,
       double* rlist, double* elist, int* iord, int* nnlog, double* chebmo)
{
    // ***begin prologue  dqawfe
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a3a1
    // ***keywords  automatic integrator, special-purpose,
    //              fourier integrals,
    //              integration between zeros with dqawoe,
    //              convergence acceleration with dqelg
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            dedoncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  the routine calculates an approximation result to a
    //             given fourier integal
    //             i = integral of f(x)*w(x) over (a,infinity)
    //             where w(x)=cos(omega*x) or w(x)=sin(omega*x),
    //             hopefully satisfying following claim for accuracy
    //             abs(i-result).le.epsabs.
    // ***description
    //
    //         computation of fourier integrals
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //          on entry
    //             f      - double precision
    //                      function subprogram defining the integrand
    //                      function f(x). the actual name for f needs to
    //                      be declared e x t e r n a l in the driver program.
    //
    //             a      - double precision
    //                      lower limit of integration
    //
    //             omega  - double precision
    //                      parameter in the weight function
    //
    //             integr - integer
    //                      indicates which weight function is used
    //                      integr = 1      w(x) = cos(omega*x)
    //                      integr = 2      w(x) = sin(omega*x)
    //                      if integr.ne.1.and.integr.ne.2, the routine will
    //                      end with ier = 6.
    //
    //             epsabs - double precision
    //                      absolute accuracy requested, epsabs.gt.0
    //                      if epsabs.le.0, the routine will end with ier = 6.
    //
    //             limlst - integer
    //                      limlst gives an upper bound on the number of
    //                      cycles, limlst.ge.1.
    //                      if limlst.lt.3, the routine will end with ier = 6.
    //
    //             limit  - integer
    //                      gives an upper bound on the number of subintervals
    //                      allowed in the partition of each cycle, limit.ge.1
    //                      each cycle, limit.ge.1.
    //
    //             maxp1  - integer
    //                      gives an upper bound on the number of
    //                      chebyshev moments which can be stored, i.e.
    //                      for the intervals of lengths abs(b-a)*2**(-l),
    //                      l=0,1, ..., maxp1-2, maxp1.ge.1
    //
    //          on return
    //             result - double precision
    //                      approximation to the integral x
    //
    //             abserr - double precision
    //                      estimate of the modulus of the absolute error,
    //                      which should equal or exceed abs(i-result)
    //
    //             neval  - integer
    //                      number of integrand evaluations
    //
    //             ier    - ier = 0 normal and reliable termination of
    //                              the routine. it is assumed that the
    //                              requested accuracy has been achieved.
    //                      ier.gt.0 abnormal termination of the routine. the
    //                              estimates for integral and error are less
    //                              reliable. it is assumed that the requested
    //                              accuracy has not been achieved.
    //             error messages
    //                     if omega.ne.0
    //                      ier = 1 maximum number of  cycles  allowed
    //                              has been achieved., i.e. of subintervals
    //                              (a+(k-1)c,a+kc) where
    //                              c = (2*int(abs(omega))+1)*pi/abs(omega),
    //                              for k = 1, 2, ..., lst.
    //                              one can allow more cycles by increasing
    //                              the value of limlst (and taking the
    //                              according dimension adjustments into
    //                              account).
    //                              examine the array iwork which contains
    //                              the error flags on the cycles, in order to
    //                              look for eventual local integration
    //                              difficulties. if the position of a local
    //                              difficulty can be determined (e.g.
    //                              singularity, discontinuity within the
    //                              interval) one will probably gain from
    //                              splitting up the interval at this point
    //                              and calling appropriate integrators on
    //                              the subranges.
    //                          = 4 the extrapolation table constructed for
    //                              convergence acceleration of the series
    //                              formed by the integral contributions over
    //                              the cycles, does not converge to within
    //                              the requested accuracy. as in the case of
    //                              ier = 1, it is advised to examine the
    //                              array iwork which contains the error
    //                              flags on the cycles.
    //                          = 6 the input is invalid because
    //                              (integr.ne.1 and integr.ne.2) or
    //                               epsabs.le.0 or limlst.lt.3.
    //                               result, abserr, neval, lst are set
    //                               to zero.
    //                          = 7 bad integrand behaviour occurs within one
    //                              or more of the cycles. location and type
    //                              of the difficulty involved can be
    //                              determined from the vector ierlst. here
    //                              lst is the number of cycles actually
    //                              needed (see below).
    //                              ierlst(k) = 1 the maximum number of
    //                                            subdivisions (= limit) has
    //                                            been achieved on the k th
    //                                            cycle.
    //                                        = 2 occurrence of roundoff error
    //                                            is detected and prevents the
    //                                            tolerance imposed on the
    //                                            k th cycle, from being
    //                                            achieved.
    //                                        = 3 extremely bad integrand
    //                                            behaviour occurs at some
    //                                            points of the k th cycle.
    //                                        = 4 the integration procedure
    //                                            over the k th cycle does
    //                                            not converge (to within the
    //                                            required accuracy) due to
    //                                            roundoff in the
    //                                            extrapolation procedure
    //                                            invoked on this cycle. it
    //                                            is assumed that the result
    //                                            on this interval is the
    //                                            best which can be obtained.
    //                                        = 5 the integral over the k th
    //                                            cycle is probably divergent
    //                                            or slowly convergent. it
    //                                            must be noted that
    //                                            divergence can occur with
    //                                            any other value of
    //                                            ierlst(k).
    //                     if omega = 0 and integr = 1,
    //                     the integral is calculated by means of dqagie
    //                     and ier = ierlst(1) (with meaning as described
    //                     for ierlst(k), k = 1).
    //
    //             rslst  - double precision
    //                      vector of dimension at least limlst
    //                      rslst(k) contains the integral contribution
    //                      over the interval (a+(k-1)c,a+kc) where
    //                      c = (2*int(abs(omega))+1)*pi/abs(omega),
    //                      k = 1, 2, ..., lst.
    //                      note that, if omega = 0, rslst(1) contains
    //                      the value of the integral over (a,infinity).
    //
    //             erlst  - double precision
    //                      vector of dimension at least limlst
    //                      erlst(k) contains the error estimate corresponding
    //                      with rslst(k).
    //
    //             ierlst - integer
    //                      vector of dimension at least limlst
    //                      ierlst(k) contains the error flag corresponding
    //                      with rslst(k). for the meaning of the local error
    //                      flags see description of output parameter ier.
    //
    //             lst    - integer
    //                      number of subintervals needed for the integration
    //                      if omega = 0 then lst is set to 1.
    //
    //             alist, blist, rlist, elist - double precision
    //                      vector of dimension at least limit,
    //
    //             iord, nnlog - integer
    //                      vector of dimension at least limit, providing
    //                      space for the quantities needed in the subdivision
    //                      process of each cycle
    //
    //             chebmo - double precision
    //                      array of dimension at least (maxp1,25), providing
    //                      space for the chebyshev moments needed within the
    //                      cycles
    //
    // ***references  (none)
    // ***routines called  d1mach,dqagie,dqawoe,dqelg
    // ***end prologue  dqawfe
    //
    //
    //             the dimension of  psum  is determined by the value of
    //             limexp in subroutine dqelg (psum must be of dimension
    //             (limexp+2) at least).
    //
    //            list of major variables
    //            -----------------------
    //
    //            c1, c2    - end points of subinterval (of length cycle)
    //            cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
    //            psum      - vector of dimension at least (limexp+2)
    //                        (see routine dqelg)
    //                        psum contains the part of the epsilon table
    //                        which is still needed for further computations.
    //                        each element of psum is a partial sum of the
    //                        series which should sum to the value of the
    //                        integral.
    //            errsum    - sum of error estimates over the subintervals,
    //                        calculated cumulatively
    //            epsa      - absolute tolerance requested over current
    //                        subinterval
    //            chebmo    - array containing the modified chebyshev
    //                        moments (see also routine dqc25f)
    //
    int ktmin, l, last, ll, iter, momcom, nev, nres, numrl2;
    double abseps, correc, cycle, c1, c2, dl, drl, ep, eps, epsa;
    double errsum, fact, p1, reseps;
    double psum[52], res3la[3] = {0.0};
    const double p = 0.9;
    const double pi = 3.14159265358979323846264338327950;

    // Test on validity of parameters.
    *result = 0;
    *abserr = 0;
    *neval = 0;
    *lst = 0;
    *ier = 6;
    if (((integr != 1) && (integr != 2)) || (epsabs <= 0.0) || (limlst < 3)) { return; }
    *ier = 0;

    if (omega == 0.0)
    {
        // Integration by dqagie if omega is zero.
        if (integr == 1)
        {
            dqagie(fcn, 0.0, 1.0, epsabs, 0.0, limit, result, abserr, neval, ier,
                   alist, blist, rlist, elist, iord, &last);
        }
        rslst[0] = *result;
        erlst[0] = *abserr;
        ierlst[0] = *ier;
        *lst = 1;
        return;
    }

    // Initializations.
    l = (int)(fabs(omega));
    dl = 2*l + 1.0;
    cycle = dl*pi / fabs(omega);
    *ier = 0;
    ktmin = 0;
    *neval = 0;
    numrl2 = -1;
    nres = 0;
    c1 = a;
    c2 = cycle + a;
    p1 = 1.0 - p;
    eps = (epsabs > (uflow / p1) ? epsabs*p1 : epsabs);
    ep = eps;
    fact = 1.0;
    correc = 0.0;
    *abserr = 0.0;
    errsum = 0.0;
    ll = 0;  // Initialized

    // Main for-loop
    for (iter = 0; iter < limlst; iter++)
    {
        *lst = iter + 1;

        epsa = eps*fact;
        dqawoe(fcn, c1, c2, omega, integr, epsa, 0.0, limit, *lst, maxp1,
               &rslst[iter], &erlst[iter], &nev, &ierlst[iter], &last,
               alist, blist, rlist, elist, iord, nnlog, &momcom, chebmo);

        *neval += nev;
        fact = fact*p;
        errsum = errsum + erlst[iter];
        drl = 50.0*fabs(rslst[iter]);

        // Test on accuracy with partial sum
        if (((errsum + drl) <= epsabs) && (iter >= 5)) { goto LINE80; }
        correc = fmax(correc, erlst[iter]);
        if (ierlst[iter] != 0)
        {
            eps = fmax(ep, correc*p1);
            *ier = 7;
        }
        if ((*ier == 7) && ((errsum + drl) <= (correc * 10.0)) && (iter > 4)) { goto LINE80; }
        numrl2 += 1;
        if (iter > 0) { goto LINE20; }
        psum[0] = rslst[0];
        goto LINE40;
LINE20:
        psum[numrl2] = psum[ll] + rslst[iter];
        if (iter == 1) { goto LINE40; }

        // Test on maximum number of subintervals.
        if (*lst == limlst) { *ier = 1; }

        // Perform new extrapolation.
        dqelg(&numrl2, psum, &reseps, &abseps, res3la, &nres);

        // Test whether extrapolated result is influenced by roundouff.
        ktmin += 1;
        if ((ktmin >= 15) && (*abserr <= 0.001*(errsum + drl))) { *ier = 4; }
        if ((abseps > *abserr) && (iter != 2)) { goto LINE30; }
        *abserr = abseps;
        *result = reseps;
        ktmin = 0;

        // If ier is not zero, check whether direct result (partial sum) or
        // extrapolated result yields the best integral approximation.
        if (((*abserr + 10.0*correc) <= epsabs) ||
            ((*abserr <= epsabs) && (10.0*correc >= epsabs))) { goto LINE60; }
LINE30:
        if ((*ier != 0) && (*ier != 7)) { goto LINE60; }
LINE40:
        ll = numrl2;
        c1 = c2;
        c2 = c2 + cycle;
    }
    // 50

    // Set final result and error estimate.
LINE60:
    *abserr = *abserr + 10.0*correc;
    if (*ier == 0) { return; }
    if ((*result != 0.0) && (psum[numrl2] != 0.0)) { goto LINE70; }
    if (*abserr > errsum) { goto LINE80; }
    if (psum[numrl2] == 0.0) { return; }
LINE70:
    if ((*abserr / fabs(*result) > (errsum+drl) / fabs(psum[numrl2]))) { goto LINE80; }
    if ((*ier >= 1) && (*ier != 7)) { *abserr = *abserr + drl; }
    return;
LINE80:
    *result = psum[numrl2];
    *abserr = errsum + drl;

    return;
}


void
dqawoe(double(*fcn)(double* x), const double a, const double b, const double omega,
       const int integr, const double epsabs, const double epsrel, const int limit,
       const int icall, const int maxp1, double* result, double* abserr, int* neval,
       int* ier, int* last, double* alist, double* blist, double* rlist,
       double* elist, int* iord, int* nnlog, int* momcom, double* chebmo)
{
    //***begin prologue  dqawoe
    //***date written   800101   (yymmdd)
    //***revision date  830518   (yymmdd)
    //***category no.  h2a2a1
    //***keywords  automatic integrator, special-purpose,
    //             integrand with oscillatory cos or sin factor,
    //             clenshaw-curtis method, (end point) singularities,
    //             extrapolation, globally adaptive
    //***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //           de doncker,elise,appl. math. & progr. div. - k.u.leuven
    //***purpose  the routine calculates an approximation result to a given
    //            definite integral
    //            i = integral of f(x)*w(x) over (a,b)
    //            where w(x) = cos(omega*x) or w(x)=sin(omega*x),
    //            hopefully satisfying following claim for accuracy
    //            abs(i-result).le.max(epsabs,epsrel*abs(i)).
    //***description
    //
    //        computation of oscillatory integrals
    //        standard fortran subroutine
    //        double precision version
    //
    //        parameters
    //         on entry
    //            f      - double precision
    //                     function subprogram defining the integrand
    //                     function f(x). the actual name for f needs to be
    //                     declared e x t e r n a l in the driver program.
    //
    //            a      - double precision
    //                     lower limit of integration
    //
    //            b      - double precision
    //                     upper limit of integration
    //
    //            omega  - double precision
    //                     parameter in the integrand weight function
    //
    //            integr - integer
    //                     indicates which of the weight functions is to be
    //                     used
    //                     integr = 1      w(x) = cos(omega*x)
    //                     integr = 2      w(x) = sin(omega*x)
    //                     if integr.ne.1 and integr.ne.2, the routine
    //                     will end with ier = 6.
    //
    //            epsabs - double precision
    //                     absolute accuracy requested
    //            epsrel - double precision
    //                     relative accuracy requested
    //                     if  epsabs.le.0
    //                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                     the routine will end with ier = 6.
    //
    //            limit  - integer
    //                     gives an upper bound on the number of subdivisions
    //                     in the partition of (a,b), limit.ge.1.
    //
    //            icall  - integer
    //                     if dqawoe is to be used only once, icall must
    //                     be set to 1.  assume that during this call, the
    //                     chebyshev moments (for clenshaw-curtis integration
    //                     of degree 24) have been computed for intervals of
    //                     lengths (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
    //                     if icall.gt.1 this means that dqawoe has been
    //                     called twice or more on intervals of the same
    //                     length abs(b-a). the chebyshev moments already
    //                     computed are then re-used in subsequent calls.
    //                     if icall.lt.1, the routine will end with ier = 6.
    //
    //            maxp1  - integer
    //                     gives an upper bound on the number of chebyshev
    //                     moments which can be stored, i.e. for the
    //                     intervals of lengths abs(b-a)*2**(-l),
    //                     l=0,1, ..., maxp1-2, maxp1.ge.1.
    //                     if maxp1.lt.1, the routine will end with ier = 6.
    //
    //         on return
    //            result - double precision
    //                     approximation to the integral
    //
    //            abserr - double precision
    //                     estimate of the modulus of the absolute error,
    //                     which should equal or exceed abs(i-result)
    //
    //            neval  - integer
    //                     number of integrand evaluations
    //
    //            ier    - integer
    //                     ier = 0 normal and reliable termination of the
    //                             routine. it is assumed that the
    //                             requested accuracy has been achieved.
    //                   - ier.gt.0 abnormal termination of the routine.
    //                             the estimates for integral and error are
    //                             less reliable. it is assumed that the
    //                             requested accuracy has not been achieved.
    //            error messages
    //                     ier = 1 maximum number of subdivisions allowed
    //                             has been achieved. one can allow more
    //                             subdivisions by increasing the value of
    //                             limit (and taking according dimension
    //                             adjustments into account). however, if
    //                             this yields no improvement it is advised
    //                             to analyze the integrand, in order to
    //                             determine the integration difficulties.
    //                             if the position of a local difficulty can
    //                             be determined (e.g. singularity,
    //                             discontinuity within the interval) one
    //                             will probably gain from splitting up the
    //                             interval at this point and calling the
    //                             integrator on the subranges. if possible,
    //                             an appropriate special-purpose integrator
    //                             should be used which is designed for
    //                             handling the type of difficulty involved.
    //                         = 2 the occurrence of roundoff error is
    //                             detected, which prevents the requested
    //                             tolerance from being achieved.
    //                             the error may be under-estimated.
    //                         = 3 extremely bad integrand behaviour occurs
    //                             at some points of the integration
    //                             interval.
    //                         = 4 the algorithm does not converge.
    //                             roundoff error is detected in the
    //                             extrapolation table.
    //                             it is presumed that the requested
    //                             tolerance cannot be achieved due to
    //                             roundoff in the extrapolation table,
    //                             and that the returned result is the
    //                             best which can be obtained.
    //                         = 5 the integral is probably divergent, or
    //                             slowly convergent. it must be noted that
    //                             divergence can occur with any other value
    //                             of ier.gt.0.
    //                         = 6 the input is invalid, because
    //                             (epsabs.le.0 and
    //                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
    //                             or (integr.ne.1 and integr.ne.2) or
    //                             icall.lt.1 or maxp1.lt.1.
    //                             result, abserr, neval, last, rlist(1),
    //                             elist(1), iord(1) and nnlog(1) are set
    //                             to zero. alist(1) and blist(1) are set
    //                             to a and b respectively.
    //
    //            last  -  integer
    //                     on return, last equals the number of
    //                     subintervals produces in the subdivision
    //                     process, which determines the number of
    //                     significant elements actually in the
    //                     work arrays.
    //            alist  - double precision
    //                     vector of dimension at least limit, the first
    //                      last  elements of which are the left
    //                     end points of the subintervals in the partition
    //                     of the given integration range (a,b)
    //
    //            blist  - double precision
    //                     vector of dimension at least limit, the first
    //                      last  elements of which are the right
    //                     end points of the subintervals in the partition
    //                     of the given integration range (a,b)
    //
    //            rlist  - double precision
    //                     vector of dimension at least limit, the first
    //                      last  elements of which are the integral
    //                     approximations on the subintervals
    //
    //            elist  - double precision
    //                     vector of dimension at least limit, the first
    //                      last  elements of which are the moduli of the
    //                     absolute error estimates on the subintervals
    //
    //            iord   - integer
    //                     vector of dimension at least limit, the first k
    //                     elements of which are pointers to the error
    //                     estimates over the subintervals,
    //                     such that elist(iord(1)), ...,
    //                     elist(iord(k)) form a decreasing sequence, with
    //                     k = last if last.le.(limit/2+2), and
    //                     k = limit+1-last otherwise.
    //
    //            nnlog  - integer
    //                     vector of dimension at least limit, containing the
    //                     subdivision levels of the subintervals, i.e.
    //                     iwork(i) = l means that the subinterval
    //                     numbered i is of length abs(b-a)*2**(1-l)
    //
    //         on entry and return
    //            momcom - integer
    //                     indicating that the chebyshev moments
    //                     have been computed for intervals of lengths
    //                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
    //                     momcom.lt.maxp1
    //
    //            chebmo - double precision
    //                     array of dimension (maxp1,25) containing the
    //                     chebyshev moments
    //
    //***references  (none)
    //***routines called  d1mach,dqc25f,dqelg,dqpsrt
    //***end prologue  dqawoe
    //
    //            the dimension of rlist2 is determined by  the value of
    //            limexp in subroutine dqelg (rlist2 should be of
    //            dimension (limexp+2) at least).
    //
    //            list of major variables
    //            -----------------------
    //
    //           alist     - list of left end points of all subintervals
    //                       considered up to now
    //           blist     - list of right end points of all subintervals
    //                       considered up to now
    //           rlist(i)  - approximation to the integral over
    //                       (alist(i),blist(i))
    //           rlist2    - array of dimension at least limexp+2
    //                       containing the part of the epsilon table
    //                       which is still needed for further computations
    //           elist(i)  - error estimate applying to rlist(i)
    //           maxerr    - pointer to the interval with largest
    //                       error estimate
    //           errmax    - elist(maxerr)
    //           erlast    - error on the interval currently subdivided
    //           area      - sum of the integrals over the subintervals
    //           errsum    - sum of the errors over the subintervals
    //           errbnd    - requested accuracy max(epsabs,epsrel*
    //                       abs(result))
    //           *****1    - variable for the left subinterval
    //           *****2    - variable for the right subinterval
    //           last      - index for subdivision
    //           nres      - number of calls to the extrapolation routine
    //           numrl2    - number of elements in rlist2. if an appropriate
    //                       approximation to the compounded integral has
    //                       been obtained it is put in rlist2(numrl2) after
    //                       numrl2 has been increased by one
    //           small     - length of the smallest interval considered
    //                       up to now, multiplied by 1.5
    //           erlarg    - sum of the errors over the intervals larger
    //                       than the smallest interval considered up to now
    //           extrap    - logical variable denoting that the routine is
    //                       attempting to perform extrapolation, i.e. before
    //                       subdividing the smallest interval we try to
    //                       decrease the value of erlarg
    //           noext     - logical variable denoting that extrapolation
    //                       is no longer allowed (true  value)
    //
    //            machine dependent constants
    //            ---------------------------
    //
    //           epmach is the largest relative spacing.
    //           uflow is the smallest positive magnitude.
    //           oflow is the largest positive magnitude.
    //
    int extrap, extall, ierror, iroff1, iroff2, iroff3, jupbnd, k, ksgn;
    int ktmin, L, maxerr, nev, noext, nres, nrmax, nrmom, numrl2;
    double abseps, area, area1, area12, area2, a1, a2, b1, b2, correc, defab1;
    double defab2, defabs, domega, dres, erlarg, erlast, errbnd, errmax, error1;
    double error12, error2, errsum, ertest, resabs, reseps, small, width;
    double rlist2[52], res3la[3];

    // last and limit, actual and max allowed no. of subintervals, are not
    // indices, however used as such. We use L for "last"-indexed arrays

    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    nnlog[0] = 0;

    if (((integr != 1) && (integr != 2)) ||
        ((epsabs <= 0.0) && (epsrel < fmax(50.0*epmach, 0.5e-28))) ||
        (icall < 1) || (maxp1 < 1)) { return; }
    *ier = 0;

    // First approximation to the integral

    domega = fabs(omega);
    nrmom = 0;
    if (icall < 2) { *momcom = 0; }
    dqc25f(fcn, a, b, domega, integr, nrmom, maxp1, 0, result, abserr, neval, &defabs, &resabs, momcom, chebmo);

    // Test on accuracy
    dres = fabs(*result);
    errbnd = fmax(epsabs, epsrel*dres);
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 0;

    if ((*abserr <= 100.0*epmach * defabs) && (*abserr > errbnd)) { *ier = 2; }
    if (limit == 1) { *ier = 1; }
    if ((*ier != 0) || (*abserr <= errbnd))
    {
        // 200
        if ((integr == 2) && (omega < 0.0)) { *result= -*result; }
        return;
    }

    // Initializations
    errmax = *abserr;
    maxerr = 0;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    extrap = 0;
    noext = 0;
    ierror = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin = 0;
    small = fabs(b - a)*0.75;
    nres = 0;
    numrl2 = -1;
    extall = 0;

    correc = 0.0;  // Initialized
    erlarg = 0.0;  // Initialized
    ertest = 0.0;  // Initialized

    if ((0.5*fabs(b - a)*domega) <= 2.0)
    {
        numrl2 = 0;
        extall = 1;
        rlist2[0] = *result;
    } else if ((0.25*fabs(b - a)*domega) <= 2.0) {
        extall = 1;
    }
    ksgn = (dres >= (1.0 - 50.0*epmach)*defabs ? 1 : -1);

    // Main for-loop.
    for (L = 1; L < limit; L++)
    {
        *last = L + 1;

        // Bisect the subinterval with the nrmax-th largest error estimate.
        nrmom = nnlog[maxerr] + 1;
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        dqc25f(fcn, a1, b1, domega, integr, nrmom, maxp1, 0, &area1, &error1, &nev, &resabs, &defab1, momcom, chebmo);
        *neval += nev;
        dqc25f(fcn, a2, b2, domega, integr, nrmom, maxp1, 1, &area2, &error2, &nev, &resabs, &defab2, momcom, chebmo);
        *neval += nev;

        // Improve previous approximations to integral and error and test for accuracy.
        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum = errsum + error12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 != error1) && (defab2 != error2))
        {
            if (!((fabs(rlist[maxerr] - area12) > 1.0e-5*fabs(area12)) || (error12 < 0.99*errmax)))
            {
                if (extrap) { iroff2++; } else { iroff1++; }
            }
            // 20
            if ((L > 9) && (error12 > errmax)) { iroff3++; }
        }
        // 25

        rlist[maxerr] = area1;
        rlist[L] = area2;
        nnlog[maxerr] = nrmom;
        nnlog[L] = nrmom;
        errbnd = fmax(epsabs, epsrel*fabs(area));

        // Test for roundoff eerror and eventually set error flag.
        if ((iroff1 + iroff2 >= 10) || iroff3 >= 20) { *ier = 2; }
        if (iroff2 >= 5) { ierror = 3; }

        // Set error flag in the case that the number of subintervals equals limit.
        if (*last == limit) { *ier = 1; }

        // Set error flag in the case of bad integrand behavior at a point of the integration range.
        if (fmax(fabs(a1), fabs(b2)) <= (1.0 + 100.0*epmach)*(fabs(a2) + 1.0e3*epmach)) { *ier = 4; }

        // Append the newly-created intervals to the list.
        if (!(error2 > error1))
        {
            alist[L] = a2;
            blist[maxerr] = b1;
            blist[L] = b2;
            elist[maxerr] = error1;
            elist[L] = error2;
        } else {
            alist[maxerr] = a2;
            alist[L] = a1;
            blist[L] = b1;
            rlist[maxerr] = area2;
            rlist[L] = area1;
            elist[maxerr] = error2;
            elist[L] = error1;
        }
        // 40

        // Call subroutine dqpsrt to maintain the descending ordering in the list
        // of error estimates and select the subinterval with nrmax-th error
        // estimate (to bisected next).
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);

        if (errsum <= errbnd) { goto LINE170; }
        if (*ier != 0) { break; }
        if ((L == 1) && (extall)) { goto LINE120; }
        if (noext) { continue; }
        if (!extall) { goto LINE50; }
        erlarg = erlarg - erlast;
        if (fabs(b1 - a1) > small) { erlarg = erlarg + error12; }
        if (extrap) { goto LINE70; }
LINE50:
        // Test whether the interval to be bisected next is the smallest interval.
        width = fabs(blist[maxerr] - alist[maxerr]);
        if (width > small) { continue; }
        if (!(extall))
        {
            // Test whether we can start with the extrapolation procedure (we do
            // this if we integrate over the next interval with use of a
            // Gauss-Kronrod rule - see dqc25f).
            small = small*0.5;
            if (0.25*width*domega > 2.0) { continue; }
            extall = 1;
            goto LINE130;
        }
        // 60
        extrap = 1;
        nrmax = 1;
LINE70:
        if ((ierror == 3) || (erlarg <= ertest)) { goto LINE90; }

        // The smallest interval has the largest error. Before bisecting, decrease
        // the sum of the erorrs over the larger intervals (erlarg) and perform
        // extrapolation.
        jupbnd = (*last > 2 + (limit/2) ? limit + 3 - *last : *last);

        for (k = nrmax; k < jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small) { goto LINE140; }
            nrmax++;
        }
        // 80
LINE90:
        // Perform extrapolation.
        numrl2++;
        rlist2[numrl2] = area;
        if (numrl2 < 2) { goto LINE110; }
        dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) { *ier = 5; }
        if (abseps >= *abserr) { goto LINE100; }
        ktmin = 0;
        *abserr = abseps;
        *result = reseps;
        correc = erlarg;
        ertest = fmax(epsabs, epsrel * fabs(reseps));
        if (*abserr <= ertest) { break; }
LINE100:
        // Prepare bisection of the smallest interval
        if (numrl2 == 0) { noext = 1; }
        if (*ier == 5) { break; }
LINE110:
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = 0;
        small = small*0.5;
        erlarg = errsum;
        continue;
LINE120:
        small = small*0.5;
        numrl2++;
        rlist2[numrl2] = area;
LINE130:
        ertest = errbnd;
        erlarg = errsum;
LINE140:
        ;  // no-op.
    }
    // 140

    // Set the final result.
    if ((*abserr == oflow) || (nres == 0)) { goto LINE170; }
    if ((*ier + ierror) == 0) { goto LINE165; }
    if (ierror == 3) { *abserr = *abserr + correc; }
    if (*ier == 0) { *ier = 3; }
    if ((*result != 0.0) && (area != 0.0)) { goto LINE160;}
    if (*abserr > errsum) { goto LINE170; }
    if (area == 0.0) { goto LINE190; }
    goto LINE165;
LINE160:
    if ((*abserr/fabs(*result)) > errsum/fabs(area)) { goto LINE170; }
LINE165:
    // Test on divergence.
    if ((ksgn == -1) && (fmax(fabs(*result), fabs(area)) <= defabs*0.01)) { goto LINE190; }
    if ((0.01 > *result/area) || (*result/area > 100.0) || (errsum >= fabs(area))) { *ier = 6; }
    goto LINE190;
LINE170:
    // Compute global integral sum.
    *result = 0.0;
    for (k = 0; k <= L; k++)
    {
        *result = *result + rlist[k];
    }
    // 180
    *abserr = errsum;
LINE190:
    if (*ier > 2) { *ier -= 1;}
    if ((integr == 2) && (omega < 0.0)) { *result= -*result; }

    return;
}


void
dqawse(double(*fcn)(double* x), const double a, const double b, const double alfa,
       const double beta, const int integr, const double epsabs, const double epsrel,
       const int limit, double* result, double* abserr, int* neval, int* ier,
       double* alist, double* blist, double* rlist, double* elist, int* iord,
       int* last)
{
    // ***begin prologue  dqawse
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a1
    // ***keywords  automatic integrator, special-purpose,
    //              algebraico-logarithmic end point singularities,
    //              clenshaw-curtis method
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  the routine calculates an approximation result to a given
    //             definite integral i = integral of f*w over (a,b),
    //             (where w shows a singular behaviour at the end points,
    //             see parameter integr).
    //             hopefully satisfying following claim for accuracy
    //             abs(i-result).le.max(epsabs,epsrel*abs(i)).
    // ***description
    //
    //         integration of functions having algebraico-logarithmic
    //         end point singularities
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //          on entry
    //             f      - double precision
    //                      function subprogram defining the integrand
    //                      function f(x). the actual name for f needs to be
    //                      declared e x t e r n a l in the driver program.
    //
    //             a      - double precision
    //                      lower limit of integration
    //
    //             b      - double precision
    //                      upper limit of integration, b.gt.a
    //                      if b.le.a, the routine will end with ier = 6.
    //
    //             alfa   - double precision
    //                      parameter in the weight function, alfa.gt.(-1)
    //                      if alfa.le.(-1), the routine will end with
    //                      ier = 6.
    //
    //             beta   - double precision
    //                      parameter in the weight function, beta.gt.(-1)
    //                      if beta.le.(-1), the routine will end with
    //                      ier = 6.
    //
    //             integr - integer
    //                      indicates which weight function is to be used
    //                      = 1  (x-a)**alfa*(b-x)**beta
    //                      = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
    //                      = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
    //                      = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    //                      if integr.lt.1 or integr.gt.4, the routine
    //                      will end with ier = 6.
    //
    //             epsabs - double precision
    //                      absolute accuracy requested
    //             epsrel - double precision
    //                      relative accuracy requested
    //                      if  epsabs.le.0
    //                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                      the routine will end with ier = 6.
    //
    //             limit  - integer
    //                      gives an upper bound on the number of subintervals
    //                      in the partition of (a,b), limit.ge.2
    //                      if limit.lt.2, the routine will end with ier = 6.
    //
    //          on return
    //             result - double precision
    //                      approximation to the integral
    //
    //             abserr - double precision
    //                      estimate of the modulus of the absolute error,
    //                      which should equal or exceed abs(i-result)
    //
    //             neval  - integer
    //                      number of integrand evaluations
    //
    //             ier    - integer
    //                      ier = 0 normal and reliable termination of the
    //                              routine. it is assumed that the requested
    //                              accuracy has been achieved.
    //                      ier.gt.0 abnormal termination of the routine
    //                              the estimates for the integral and error
    //                              are less reliable. it is assumed that the
    //                              requested accuracy has not been achieved.
    //             error messages
    //                          = 1 maximum number of subdivisions allowed
    //                              has been achieved. one can allow more
    //                              subdivisions by increasing the value of
    //                              limit. however, if this yields no
    //                              improvement, it is advised to analyze the
    //                              integrand in order to determine the
    //                              integration difficulties which prevent the
    //                              requested tolerance from being achieved.
    //                              in case of a jump discontinuity or a local
    //                              singularity of algebraico-logarithmic type
    //                              at one or more interior points of the
    //                              integration range, one should proceed by
    //                              splitting up the interval at these
    //                              points and calling the integrator on the
    //                              subranges.
    //                          = 2 the occurrence of roundoff error is
    //                              detected, which prevents the requested
    //                              tolerance from being achieved.
    //                          = 3 extremely bad integrand behaviour occurs
    //                              at some points of the integration
    //                              interval.
    //                          = 6 the input is invalid, because
    //                              b.le.a or alfa.le.(-1) or beta.le.(-1), or
    //                              integr.lt.1 or integr.gt.4, or
    //                              (epsabs.le.0 and
    //                               epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                              or limit.lt.2.
    //                              result, abserr, neval, rlist(1), elist(1),
    //                              iord(1) and last are set to zero. alist(1)
    //                              and blist(1) are set to a and b
    //                              respectively.
    //
    //             alist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the left
    //                      end points of the subintervals in the partition
    //                      of the given integration range (a,b)
    //
    //             blist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the right
    //                      end points of the subintervals in the partition
    //                      of the given integration range (a,b)
    //
    //             rlist  - double precision
    //                      vector of dimension at least limit,the first
    //                       last  elements of which are the integral
    //                      approximations on the subintervals
    //
    //             elist  - double precision
    //                      vector of dimension at least limit, the first
    //                       last  elements of which are the moduli of the
    //                      absolute error estimates on the subintervals
    //
    //             iord   - integer
    //                      vector of dimension at least limit, the first k
    //                      of which are pointers to the error
    //                      estimates over the subintervals, so that
    //                      elist(iord(1)), ..., elist(iord(k)) with k = last
    //                      if last.le.(limit/2+2), and k = limit+1-last
    //                      otherwise form a decreasing sequence
    //
    //             last   - integer
    //                      number of subintervals actually produced in
    //                      the subdivision process
    //
    // ***references  (none)
    // ***routines called  d1mach,dqc25s,dqmomo,dqpsrt
    // ***end prologue  dqawse
    //
    //             list of major variables
    //             -----------------------
    //
    //            alist     - list of left end points of all subintervals
    //                        considered up to now
    //            blist     - list of right end points of all subintervals
    //                        considered up to now
    //            rlist(i)  - approximation to the integral over
    //                        (alist(i),blist(i))
    //            elist(i)  - error estimate applying to rlist(i)
    //            maxerr    - pointer to the interval with largest
    //                        error estimate
    //            errmax    - elist(maxerr)
    //            area      - sum of the integrals over the subintervals
    //            errsum    - sum of the errors over the subintervals
    //            errbnd    - requested accuracy max(epsabs,epsrel*
    //                        abs(result))
    //            *****1    - variable for the left subinterval
    //            *****2    - variable for the right subinterval
    //            last      - index for subdivision
    //
    //
    //             machine dependent constants
    //             ---------------------------
    //
    //            epmach is the largest relative spacing.
    //            uflow is the smallest positive magnitude.
    //
    int iroff1, iroff2, k, L, maxerr, nev, nrmax;
    double area, area1, area12, area2, a1, a2, b1, b2, centre, errbnd, errmax;
    double error1, error12, error2, errsum, resas1, resas2;
    double ri[25], rj[25], rh[25], rg[25];

    // last and limit, actual and max allowed no. of subintervals, are not
    // indices, however used as such. We use L for "last"-indexed arrays

    *ier = 6;
    *neval = 0;
    *last = 0;
    L = 0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    *result = 0.0;
    *abserr = 0.0;

    if ((b <= a) || ((epsabs == 0.0) && (epsrel < fmax(50*epmach, 0.5e-28))) ||
        (alfa <= -1.0) || (beta <= -1.0) || (integr < 1) || (integr > 4) ||
        (limit < 2)) { return; }
    *ier = 0;

    // Compute the modified Chebyshev moments.
    dqmomo(alfa, beta, ri, rj, rg, rh, integr);

    // Integrate over the intervals (a, (a+b)/2) and ((a+b)/2, b).
    centre = 0.5*(b + a);
    dqc25s(fcn, a, b, a, centre, alfa, beta, ri, rj, rg, rh, &area1, &error1, &resas1, integr, &nev);
    *neval = nev;
    dqc25s(fcn, a, b, centre, b, alfa, beta, ri, rj, rg, rh, &area2, &error2, &resas2, integr, &nev);
    *last = 2;
    *neval += nev;
    *result = area1 + area2;
    *abserr = error1 + error2;

    // Test on Accuracy.
    errbnd = fmax(epsabs, epsrel*fabs(*result));

    // Initialization
    if (!(error2 > error1))
    {
        alist[0] = a;
        alist[1] = centre;
        blist[0] = centre;
        blist[1] = b;
        rlist[0] = area1;
        rlist[1] = area2;
        elist[0] = error1;
        elist[1] = error2;
    } else {
        alist[0] = centre;
        alist[1] = a;
        blist[0] = b;
        blist[1] = centre;
        rlist[0] = area2;
        rlist[1] = area1;
        elist[0] = error2;
        elist[1] = error1;
    }
    // 20

    iord[0] = 0;
    iord[1] = 1;
    if (limit == 2) { *ier = 1; }
    if ((*abserr <= errbnd) || (*ier == 1)) { return; }

    errmax = elist[0];
    maxerr = 0;
    nrmax = 0;
    area = *result;
    errsum = *abserr;
    iroff1 = 0;
    iroff2 = 0;

    // Main for-loop.
    for (L = 2; L < limit; L++)
    {
        *last = L + 1;

        // Bisect the subinterval with largest error estimate.
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];

        dqc25s(fcn, a, b, a1, b1, alfa, beta, ri, rj, rg, rh, &area1, &error1, &resas1, integr, &nev);
        *neval += nev;
        dqc25s(fcn, a, b, a2, b2, alfa, beta, ri, rj, rg, rh, &area2, &error2, &resas2, integr, &nev);
        *neval += nev;

        // Improve previous approximations integral and error and test for accuracy.

        area12 = area1 + area2;
        error12 = error1 + error2;
        errsum = errsum + error12 - errmax;
        area = area + area12 - rlist[maxerr];
        if (!((a == a1) || (b == b2))){
            if ((resas1 != error1) && (resas2 != error2))
            {
                // Test for roundoff error.
                if ((fabs(rlist[maxerr]-area12) < (1.0e-5 * fabs(area12))) && (error12 >= (0.99 *errmax)))
                {
                    iroff1++;
                }
                if ((L > 9) && (error12 > errmax))
                {
                    iroff2++;
                }
            }
        }
        // 30
        rlist[maxerr] = area1;
        rlist[L] = area2;

        // Test on accuracy.
        errbnd = fmax(epsabs, epsrel*fabs(area));
        if (!(errsum <= errbnd))
        {
            // Set error flag in the case that the number of interval bisections
            // exceeds limit.
            if (*last == limit) { *ier = 1; }

            // Set error flag in the case of roundoff error.
            if ((iroff1 >= 6) || (iroff2 >= 20)) { *ier = 2; }

            // Set error flag in the case of bad integrand behavior at interior
            // points of integration range.
            if (fmax(fabs(a1), fabs(b2)) <= ((1.0 + 100.0*epmach)*(fabs(a2) + 1.0e3*uflow)))
            {
                *ier = 3;
            }
        }
        // 35

        // Append the newly-created intervals to the list.
        if (!(error2 > error1))
        {
            alist[L] = a2;
            blist[maxerr] = b1;
            blist[L] = b2;
            elist[maxerr] = error1;
            elist[L] = error2;
        } else {
            alist[maxerr] = a2;
            alist[L] = a1;
            blist[L] = b1;
            rlist[maxerr] = area2;
            rlist[L] = area1;
            elist[maxerr] = error2;
            elist[L] = error1;
        }
        // 50

        // Call subroutine dqpsrt to maintain the descending ordering in the
        // list of error estimates and select the subinterval with largest error
        // estimate (to be bisected next).
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);
        if ((*ier != 0) || (errsum <= errbnd)) { break; }
        // 60
    }
    // 70

    // Compute the final result.
    *result = 0.0;
    for (k = 0; k <= L; k++)
    {
        *result = *result + rlist[k];
    }
    *abserr = errsum;

    return;
}


void
dqc25c(double(*fcn)(double* x), const double a, const double b, const double c,
       double* result, double* abserr, int* krul, int* neval)
{
    // ***begin prologue  dqc25c
    // ***date written   810101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a2,j4
    // ***keywords  25-point clenshaw-curtis integration
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f*w over (a,b) with
    //             error estimate, where w(x) = 1/(x-c)
    // ***description
    //
    //         integration rules for the computation of cauchy
    //         principal value integrals
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //            f      - double precision
    //                     function subprogram defining the integrand function
    //                     f(x). the actual name for f needs to be declared
    //                     e x t e r n a l  in the driver program.
    //
    //            a      - double precision
    //                     left end point of the integration interval
    //
    //            b      - double precision
    //                     right end point of the integration interval, b.gt.a
    //
    //            c      - double precision
    //                     parameter in the weight function
    //
    //            result - double precision
    //                     approximation to the integral
    //                     result is computed by using a generalized
    //                     clenshaw-curtis method if c lies within ten percent
    //                     of the integration interval. in the other case the
    //                     15-point kronrod rule obtained by optimal addition
    //                     of abscissae to the 7-point gauss rule, is applied.
    //
    //            abserr - double precision
    //                     estimate of the modulus of the absolute error,
    //                     which should equal or exceed abs(i-result)
    //
    //            krul   - integer
    //                     key which is decreased by 1 if the 15-point
    //                     gauss-kronrod scheme has been used
    //
    //            neval  - integer
    //                     number of integrand evaluations
    //
    // .......................................................................
    // ***references  (none)
    // ***routines called  dqcheb,dqk15w,dqwgtc
    // ***end prologue  dqc25c
    //
    int i, isym, k, kp;
    double ak22, amom0, amom1, amom2, cc, centr, hlgth, resabs, resasc, res12, res24, u;
    double fval[25], cheb12[13], cheb24[25];
    double mut_D;

    static const double x[11]= {
       0.991444861373810411144557526928563,
       0.965925826289068286749743199728897,
       0.923879532511286756128183189396788,
       0.866025403784438646763723170752936,
       0.793353340291235164579776961501299,
       0.707106781186547524400844362104849,
       0.608761429008720639416097542898164,
       0.500000000000000000000000000000000,
       0.382683432365089771728459984030399,
       0.258819045102520762348898837624048,
       0.130526192220051591548406227895489
    };

    cc = (2.0*c - b - a)/(b - a);
    kp = 0;  // Not used in dqwgtc, hence value is irrelevant, just initialized.
    if (fabs(cc) >= 1.1)
    {
        // Apply the 15-point Gauss-Kronrod scheme.
        *krul -= 1;
        dqk15w(fcn, dqwgtc, c, 0.0, 0.0, 0.0, kp, a, b, result, abserr, &resabs, &resasc);
        *neval = 15;
        if (resasc == *abserr) { *krul += 1; }
        return;
    }
    // 10

    // Use the generalized clenshaw-curtis method.
    hlgth = 0.5 * (b - a);
    centr = 0.5 * (b + a);
    *neval = 25;

    mut_D = centr + hlgth;
    fval[0]  = 0.5 * (*fcn)(&mut_D);
    fval[12] = (*fcn)(&centr);
    mut_D = centr - hlgth;
    fval[24] = 0.5 * (*fcn)(&mut_D);

    for (i = 1; i < 12; i++)
    {
        u = hlgth*x[i-1];
        isym = 24 - i;
        mut_D = centr + u;
        fval[i] = (*fcn)(&mut_D);
        mut_D = centr - u;
        fval[isym] = (*fcn)(&mut_D);
    }

    // Compute the Chebyshev series expansion.
    dqcheb(x, fval, cheb12, cheb24);

    // The modified Chebyshev moments are computed by forward recursion, using
    // amom0 and amom1 as starting values.
    amom0 = log(fabs((1.0 - cc) / (1.0 + cc)));
    amom1 = 2.0 + (cc * amom0);
    res12 = cheb12[0]*amom0 + cheb12[1]*amom1;
    res24 = cheb24[0]*amom0 + cheb24[1]*amom1;
    for (k = 2; k < 13; k++) {
        amom2 = 2.0*cc*amom1 - amom0;
        ak22 = (k-1) * (k-1);
        if ((k % 2) == 1) { amom2 = amom2 - (4.0 / (ak22 - 1.0)); }
        res12 = res12 + cheb12[k] * amom2;
        res24 = res24 + cheb24[k] * amom2;
        amom0 = amom1;
        amom1 = amom2;
    }
    // 30

    for (k = 13; k < 25; k++) {
        amom2 = 2.0*cc*amom1 - amom0;
        ak22 = (k-1)*(k-1);
        if ((k % 2) == 1) { amom2 = amom2 - (4.0 / (ak22 - 1.0)); }
        res24 = res24 + cheb24[k] * amom2;
        amom0 = amom1;
        amom1 = amom2;
    }
    *result = res24;
    *abserr = fabs(res24 - res12);

    return;
}


void
dqc25f(double(*fcn)(double* x), const double a, const double b, const double omega,
       const int integr, const int nrmom, const int maxp1, const int ksave,
       double* result, double* abserr, int* neval, double* resabs, double* resasc,
       int* momcom, double* chebmo)
{
    // ***begin prologue  dqc25f
    // ***date written   810101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a2
    // ***keywords  integration rules for functions with cos or sin
    //              factor, clenshaw-curtis, gauss-kronrod
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute the integral i=integral of f(x) over (a,b)
    //             where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
    //             compute j = integral of abs(f) over (a,b). for small value
    //             of omega or small intervals (a,b) the 15-point gauss-kronro
    //             rule is used. otherwise a generalized clenshaw-curtis
    //             method is used.
    // ***description
    //
    //         integration rules for functions with cos or sin factor
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //          on entry
    //            f      - double precision
    //                     function subprogram defining the integrand
    //                     function f(x). the actual name for f needs to
    //                     be declared e x t e r n a l in the calling program.
    //
    //            a      - double precision
    //                     lower limit of integration
    //
    //            b      - double precision
    //                     upper limit of integration
    //
    //            omega  - double precision
    //                     parameter in the weight function
    //
    //            integr - integer
    //                     indicates which weight function is to be used
    //                        integr = 1   w(x) = cos(omega*x)
    //                        integr = 2   w(x) = sin(omega*x)
    //
    //            nrmom  - integer
    //                     the length of interval (a,b) is equal to the length
    //                     of the original integration interval divided by
    //                     2**nrmom (we suppose that the routine is used in an
    //                     adaptive integration process, otherwise set
    //                     nrmom = 0). nrmom must be zero at the first call.
    //
    //            maxp1  - integer
    //                     gives an upper bound on the number of chebyshev
    //                     moments which can be stored, i.e. for the
    //                     intervals of lengths abs(bb-aa)*2**(-l),
    //                     l = 0,1,2, ..., maxp1-2.
    //
    //            ksave  - integer
    //                     key which is one when the moments for the
    //                     current interval have been computed
    //
    //          on return
    //            result - double precision
    //                     approximation to the integral i
    //
    //            abserr - double precision
    //                     estimate of the modulus of the absolute
    //                     error, which should equal or exceed abs(i-result)
    //
    //            neval  - integer
    //                     number of integrand evaluations
    //
    //            resabs - double precision
    //                     approximation to the integral j
    //
    //            resasc - double precision
    //                     approximation to the integral of abs(f-i/(b-a))
    //
    //          on entry and return
    //            momcom - integer
    //                     for each interval length we need to compute the
    //                     chebyshev moments. momcom counts the number of
    //                     intervals for which these moments have already been
    //                     computed. if nrmom.lt.momcom or ksave = 1, the
    //                     chebyshev moments for the interval (a,b) have
    //                     already been computed and stored, otherwise we
    //                     compute them and we increase momcom.
    //
    //            chebmo - double precision
    //                     array of dimension at least (maxp1,25) containing
    //                     the modified chebyshev moments for the first momcom
    //                     momcom interval lengths
    //
    //  ......................................................................
    // ***references  (none)
    // ***routines called  d1mach,dgtsl,dqcheb,dqk15w,dqwgtf
    // ***end prologue  dqc25f
    //
    int i, isym, j, k, m, noequ;
    double ac, an, an2, as, asap, ass, centr, conc, cons, cospar, estc, ests;
    double hlgth, parint, par2, par22, resc12, resc24, ress12, ress24, sinpar;
    double cheb12[13], cheb24[25], d[25], d1[25], d2[25], fval[25], v[28];
    double mut_D, temp, fact;
    m = 0;

    static const double x[11]= {
       0.991444861373810411144557526928563,
       0.965925826289068286749743199728897,
       0.923879532511286756128183189396788,
       0.866025403784438646763723170752936,
       0.793353340291235164579776961501299,
       0.707106781186547524400844362104849,
       0.608761429008720639416097542898164,
       0.500000000000000000000000000000000,
       0.382683432365089771728459984030399,
       0.258819045102520762348898837624048,
       0.130526192220051591548406227895489
    };

    centr = 0.5*(b + a);
    hlgth = 0.5*(b - a);
    parint = omega*hlgth;

    // Compute the integral using the 15-point gauss-kronrod formula if the
    // value of the parameter in the integrand is small.
    if (fabs(parint) <= 2.0) {
        dqk15w(fcn, dqwgtf, omega, 0.0, 0.0, 0.0, integr, a, b, result, abserr, resabs, resasc);
        *neval = 15;
        return;
    }

    // Compute the integral using the generalized clenshaw-curtis method.
    conc = hlgth * cos(centr * omega);
    cons = hlgth * sin(centr * omega);
    *resasc = oflow;
    *neval = 25;

    // Check whether the Chebyshev moments for this interval have already been computed.
    if ((nrmom >= *momcom) && (ksave != 1))
    {
        // Compute a new set of Chebyshev moments.
        m = *momcom;
        par2 = parint * parint;
        par22 = par2 + 2.0;
        sinpar = sin(parint);
        cospar = cos(parint);

        // Compute the Chebyshev moments with respect to cosine.
        v[0] = 2.0 * sinpar / parint;
        v[1] = (8.0*cospar + (par2 + par2 - 8.0)*sinpar / parint) / par2;
        v[2] = (32.0*(par2 - 12.0)*cospar + (2.0*((par2 - 80.0)*par2 + 192.0)*sinpar) / parint) / (par2*par2);
        ac = 8.0 * cospar;
        as = 24.0 * parint * sinpar;

        if (fabs(parint) <= 24.0)
        {
            // compute the chebyshev moments as the solutions of a boundary value
            // problem with 1 initial value (v(3)) and 1 end value (computed using
            // an asymptotic formula).
            noequ = 24;
            an = 6.0;
            for (k = 0; k < noequ; k++)
            {
                an2 = an * an;
                d[k] = -2.0*(an2 - 4.0)*(par22 - an2 - an2);
                d2[k] = (an - 1.0)*(an - 2.0)*par2;
                d1[k+1] = (an + 3.0)*(an + 4.0)*par2;
                v[k+3] = as - (an2 - 4.0)*ac;
                an = an + 2.0;
            }
            // 20

            an2 = an * an;
            d[noequ] = -2.0*(an2 - 4.0)*(par22 - an2 - an2);
            v[noequ+3] = as - (an2 - 4.0)*ac;
            v[3] = v[3] - (56.0 * par2 * v[2]);
            ass = parint * sinpar;
            asap = (((((210.0*par2 -1.0)*cospar - (105.0*par2 - 63.0)*ass)/an2
                    - (1.0 - 15.0*par2)*cospar + 15.0 * ass)/an2
                    - cospar + 3.0*ass)/an2 - cospar) / an2;
            v[noequ + 3] = v[noequ + 3] - (2.0*asap * par2 * (an - 1.0) * (an - 2.0));

            // solve the tridiagonal system by means of gaussian
            // elimination with partial pivoting.

            // SciPy Translation Note: Copied from LAPACK DGTSV
            // Assuming nonsingularity according to original quadpack code
            // d1 starts from second element, v starts from third element.
            // Equivalent FORTRAN DGTSV call:
            // call dgtsv(noequ,1,d1(2),d,d2,v(4),noequ,iers)

            for (i = 0; i < noequ-1; i++) {
                if (fabs(d[i]) >= fabs(d1[i+1]))
                {
                    // No row interchange required.
                    fact = d1[i+1] / d[i];
                    d[i+1] -= fact*d2[i];
                    v[3+i+1] -= fact*v[3+i];
                    d1[i+1] = 0.0;
                } else {
                    // Interchange rows i and i + 1
                    fact = d[i] / d1[i+1];
                    d[i] = d1[i+1];
                    temp = d[i+1];
                    d[i+1] = d2[i] - fact*temp;
                    if (i != noequ-1)
                    {
                        d1[i+1] = d2[i+1];
                        d2[i+1] = -fact*d1[i+1];
                    }
                    d2[i] = temp;
                    temp = v[3+i];
                    v[3+i] = v[3+i+1];
                    v[3+i+1] = temp - fact*v[3+i+1];
                }
            }
            // Back substitute
            v[3+noequ] /= d[noequ];
            v[3+noequ-1] = (v[3+noequ-1] - d2[noequ-1]*v[3+noequ])/d[noequ-1];

            for (i = noequ-2; i >= 0; i--) {
                v[3+i] = (v[3+i] - d2[i]*v[3+i+1] - d1[i+1]*v[3+i+2]) / d[i];
            }

        } else {
            // 30
            // Compute the chebyshev moments by means of forward recursion.

            an = 4.0;
            for (i = 3; i < 13; i++) {
                an2 = an * an;
                v[i] = ((an2 - 4.0) * (2.0*(par22 - an2 - an2)*v[i-1] - ac)
                        + as - par2*(an+1.0)*(an+2.0)*v[i-2]) / (par2*(an-1.0)*(an-2.0));
                an = an + 2.0;
            }
            // 40
        }
        // 50
        for (j = 0; j < 13; j++)
        {
            chebmo[*momcom + maxp1*(2*j)] = v[j];
        }
        // 60

        // Compute the chebyshev moments with respect to sine.
        v[0] = 2.0 * (sinpar - parint * cospar) / par2;
        v[1] = (18.0 - 48.0/par2)*sinpar/par2 + (-2.0 + 48.0/par2)*cospar/parint;
        ac = -24.0 * parint * cospar;
        as = -8.0 * sinpar;

        if (fabs(parint) <= 24.0)
        {
            // Compute the chebyshev moments as the solutions of a boundary value
            // problem with 1 initial value (v(2)) and 1 end value (computed using
            // an asymptotic formula).
            an = 5.0;
            for (k = 0; k < noequ; k++)
            {
                an2 = an * an;
                d[k] = -2.0*(an2 - 4.0)*(par22 - an2 - an2);
                d2[k] = (an - 1.0)*(an - 2.0)*par2;
                d1[k+1] = (an + 3.0)*(an + 4.0)*par2;
                v[k+2] = ac + (an2 - 4.0)*as;
                an = an + 2.0;
            }
            // 70

            an2 = an * an;
            d[noequ] = -2.0*(an2 - 4.0)*(par22 - an2 - an2);
            v[noequ+2] = ac + (an2 - 4.0)*as;
            v[2] = v[2] - 42.0*par2*v[1];
            ass = parint*cospar;
            asap = (((((105.0*par2-63.0)*ass + (210.0*par2-1.0)*sinpar)/an2
                    + (15.0*par2 - 1.0)*sinpar - 15.0*ass)/an2
                    - 3.0*ass - sinpar)/an2-sinpar)/an2;
            v[noequ+2] = v[noequ+2] - 2.0*asap*par2*(an-1.0)*(an-2.0);

            // Equivalent FORTRAN DGTSV call:
            // call dgtsv(noequ,1,d1(2),d,d2,v(3),noequ,iers)
            for (i = 0; i < noequ-1; i++) {
                if (fabs(d[i]) >= fabs(d1[i+1]))
                {
                    // No row interchange required.
                    fact = d1[i+1] / d[i];
                    d[i+1] -= fact*d2[i];
                    v[2+i+1] -= fact*v[2+i];
                    d1[i+1] = 0.0;
                } else {
                    // Interchange rows i and i + 1
                    fact = d[i] / d1[i+1];
                    d[i] = d1[i+1];
                    temp = d[i+1];
                    d[i+1] = d2[i] - fact*temp;
                    if (i != noequ-1)
                    {
                        d1[i+1] = d2[i+1];
                        d2[i+1] = -fact*d1[i+1];
                    }
                    d2[i] = temp;
                    temp = v[2+i];
                    v[2+i] = v[2+i+1];
                    v[2+i+1] = temp - fact*v[2+i+1];
                }
            }
            // Back substitute
            v[2+noequ] /= d[noequ];
            v[2+noequ-1] = (v[2+noequ-1] - d2[noequ-1]*v[2+noequ])/d[noequ-1];

            for (i = noequ-2; i >= 0; i--) {
                v[2+i] = (v[2+i] - d2[i]*v[2+i+1] - d1[i+1]*v[2+i+2]) / d[i];
            }

        } else {
            // 80
            an = 3.0;
            for (i = 2; i < 12; i++)
            {
                an2 = an*an;
                v[i] = ((an2-4.0)*(2.0*(par22-an2-an2)*v[i-1]+as)
                    + ac-par2*(an+1.0)*(an+2.0)*v[i-2])
                    / (par2*(an-1.0)*(an-2.0));
                an = an + 2.0;
            }
            // 90
        }
        // 100
        for (j = 0; j < 12; j++)
        {
            chebmo[m + maxp1*(2*j + 1)] = v[j];
        }
        // 110
    }
    // 120
    if (nrmom < *momcom)
    {
        m = nrmom;
    }
    if ((*momcom < (maxp1 - 1)) && (nrmom >= *momcom))
    {
        *momcom += 1;
    }
    // 120

    // Compute the coefficients of the chebyshev expansions of degrees 12
    // and 24 of the function f.
    mut_D = centr + hlgth;
    fval[0]  = 0.5 * (*fcn)(&mut_D);
    fval[12] = (*fcn)(&centr);
    mut_D = centr - hlgth;
    fval[24] = 0.5 * (*fcn)(&mut_D);
    for (i = 1; i < 12; i++)
    {
        isym = 24 - i;
        mut_D = centr + hlgth*x[i-1];
        fval[i] = (*fcn)(&mut_D);
        mut_D = centr - hlgth*x[i-1];
        fval[isym] = (*fcn)(&mut_D);
    }
    // 130

    dqcheb(x, fval, cheb12, cheb24);

    // Compute the integral and error estimates.
    resc12 = cheb12[12]*chebmo[m + maxp1*12];
    ress12 = 0.0;
    for (k = 10; k >= 0; k -= 2)
    {
        resc12 = resc12 + cheb12[k]*chebmo[m + maxp1*k];
        ress12 = ress12 + cheb12[k+1]*chebmo[m + maxp1*(k+1)];
    }
    // 140

    resc24 = cheb24[24]*chebmo[m + maxp1*24];
    ress24 = 0.0;
    *resabs = fabs(cheb24[24]);
    for (k = 22; k >= 0; k -= 2)
    {
        resc24 = resc24 + cheb24[k]*chebmo[m + maxp1*k];
        ress24 = ress24 + cheb24[k+1]*chebmo[m + maxp1*(k+1)];
        *resabs = *resabs + fabs(cheb24[k]) + fabs(cheb24[k+1]);
    }
    // 150

    estc = fabs(resc24 - resc12);
    ests = fabs(ress24 - ress12);
    *resabs = (*resabs)*fabs(hlgth);
    if (integr != 2)
    {
        *result = conc*resc24 - cons*ress24;
        *abserr = fabs(conc*estc) + fabs(cons*ests);
    } else {
        *result = conc*ress24 + cons*resc24;
        *abserr = fabs(conc*ests) + fabs(cons*estc);
    }

    return;
}


void
dqc25s(double(*fcn)(double* x), const double a, const double b, const double bl,
       const double br, const double alfa, const double beta, double* ri, double* rj,
       double* rg, double* rh, double* result, double* abserr, double* resasc,
       const int integr, int* nev)
{
    // ***begin prologue  dqc25s
    // ***date written   810101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a2
    // ***keywords  25-point clenshaw-curtis integration
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f*w over (bl,br), with error
    //             estimate, where the weight function w has a singular
    //             behaviour of algebraico-logarithmic type at the points
    //             a and/or b. (bl,br) is a part of (a,b).
    // ***description
    //
    //         integration rules for integrands having algebraico-logarithmic
    //         end point singularities
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //            f      - double precision
    //                     function subprogram defining the integrand
    //                     f(x). the actual name for f needs to be declared
    //                     e x t e r n a l  in the driver program.
    //
    //            a      - double precision
    //                     left end point of the original interval
    //
    //            b      - double precision
    //                     right end point of the original interval, b.gt.a
    //
    //            bl     - double precision
    //                     lower limit of integration, bl.ge.a
    //
    //            br     - double precision
    //                     upper limit of integration, br.le.b
    //
    //            alfa   - double precision
    //                     parameter in the weight function
    //
    //            beta   - double precision
    //                     parameter in the weight function
    //
    //            ri,rj,rg,rh - double precision
    //                     modified chebyshev moments for the application
    //                     of the generalized clenshaw-curtis
    //                     method (computed in subroutine dqmomo)
    //
    //            result - double precision
    //                     approximation to the integral
    //                     result is computed by using a generalized
    //                     clenshaw-curtis method if b1 = a or br = b.
    //                     in all other cases the 15-point kronrod
    //                     rule is applied, obtained by optimal addition of
    //                     abscissae to the 7-point gauss rule.
    //
    //            abserr - double precision
    //                     estimate of the modulus of the absolute error,
    //                     which should equal or exceed abs(i-result)
    //
    //            resasc - double precision
    //                     approximation to the integral of abs(f*w-i/(b-a))
    //
    //            integr - integer
    //                     which determines the weight function
    //                     = 1   w(x) = (x-a)**alfa*(b-x)**beta
    //                     = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
    //                     = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
    //                     = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
    //                                  log(b-x)
    //
    //            nev    - integer
    //                     number of integrand evaluations
    // ***references  (none)
    // ***routines called  dqcheb,dqk15w
    // ***end prologue  dqc25s
    //
    int i, isym;
    double centr, dc, factor, fix, hlgth, mut_D, resabs, res12, res24, u;
    double cheb12[13] = { 0.0 };
    double cheb24[25] = { 0.0 };
    double fval[25] = { 0.0 };

    static const double x[11]= {
       0.991444861373810411144557526928563,
       0.965925826289068286749743199728897,
       0.923879532511286756128183189396788,
       0.866025403784438646763723170752936,
       0.793353340291235164579776961501299,
       0.707106781186547524400844362104849,
       0.608761429008720639416097542898164,
       0.500000000000000000000000000000000,
       0.382683432365089771728459984030399,
       0.258819045102520762348898837624048,
       0.130526192220051591548406227895489
    };

    *nev = 25;
    if ((bl == a) && ((alfa != 0.0) || (integr == 2) || (integr == 4)))
    {
        // this part of the program is executed only if a = bl.
        // ----------------------------------------------------
        // compute the chebyshev series expansion of the
        // following function
        // f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta*f(0.5*(br-a)*x+0.5*(br+a))

        // 10
        hlgth = 0.5*(br - bl);
        centr = 0.5*(br + bl);
        fix = b - centr;

        mut_D = centr + hlgth;
        fval[0]  = 0.5 * (*fcn)(&mut_D) * pow(fix-hlgth, beta);
        fval[12] = (*fcn)(&centr) * pow(fix, beta);
        mut_D = centr - hlgth;
        fval[24] = 0.5 * (*fcn)(&mut_D) * pow(fix+hlgth, beta);

        for (i = 1; i < 12; i++)
        {
            u = hlgth * x[i-1];
            isym = 24 - i;

            mut_D = centr + u;
            fval[i] = (*fcn)(&mut_D) * pow(fix - u, beta);
            mut_D = centr - u;
            fval[isym] = (*fcn)(&mut_D) * pow(fix + u, beta);
        }
        // 20
        factor = pow(hlgth, alfa + 1.0);
        *result = 0.0;
        *abserr = 0.0;
        res12 = 0.0;
        res24 = 0.0;

        if (integr <= 2)
        {
            dqcheb(x, fval, cheb12, cheb24);

            // integr = 1 or 2
            for (i = 0; i < 13; i++)
            {
                res12 = res12 + cheb12[i] * ri[i];
            }
            // 30
            for (i = 0; i < 25; i++)
            {
                res24 = res24 + cheb24[i] * ri[i];
            }
            // 40
            if (integr == 1)
            {
                // 130 -> 270
                *result = (*result + res24) * factor;
                *abserr = (*abserr + fabs(res24-res12)) * factor;
                return;
            }

            // integr = 2
            dc = log(br - bl);
            *result = res24*dc;
            *abserr = fabs((res24 - res12)*dc);
            res12 = 0.0;
            res24 = 0.0;
            for (i = 0; i < 13; i++)
            {
                res12 = res12 + cheb12[i] * rg[i];
            }
            // 50
            for (i = 0; i < 25; i++) {
                res24 = res24 + cheb24[i] * rg[i];
            }
            // 60

            // 130 -> 270
            *result = (*result + res24) * factor;
            *abserr = (*abserr + fabs(res24-res12)) * factor;
            return;

        }
        // 70

        // compute the chebyshev series expansion of the following function
        // f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
        fval[0] = fval[0]*log(fix - hlgth);
        fval[12] = fval[12]*log(fix);
        fval[24] = fval[24]*log(fix + hlgth);
        for (i = 1; i < 12; i++)
        {
            u = hlgth * x[i-1];
            isym = 24 - i;
            fval[i] = fval[i]*log(fix - u);
            fval[isym] = fval[isym]*log(fix + u);
        }
        // 80

        dqcheb(x,fval,cheb12,cheb24);

        //  integr = 3 or 4
        for (i = 0; i < 13; i++)
        {
            res12 = res12 + cheb12[i] * ri[i];
        }
        // 90
        for (i = 0; i < 25; i++)
        {
            res24 = res24 + cheb24[i] * ri[i];
        }
        // 100

        if (integr == 3)
        {
            // 130 -> 270
            *result = (*result + res24) * factor;
            *abserr = (*abserr + fabs(res24-res12)) * factor;
            return;
        }
        // integr = 4
        dc = log(br - bl);
        *result = res24 * dc;
        *abserr = fabs((res24-res12)*dc);
        res12 = 0.0;
        res24 = 0.0;
        for (i = 0; i < 13; i++)
        {
            res12 = res12 + cheb12[i] * rg[i];
        }
        // 110
        for (i = 0; i < 25; i++)
        {
            res24 = res24 + cheb24[i] * rg[i];
        }
        // 120
        *result = (*result + res24) * factor;
        *abserr = (*abserr + fabs(res24-res12)) * factor;
        return;
    }

    if ((br == b) && ((beta != 0.0) || (integr == 3) || (integr == 4)))
    {
        // 140
        // this part of the program is executed only if b = br.
        // ----------------------------------------------------
        // compute the chebyshev series expansion of the following function
        // f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))

        hlgth = 0.5*(b - bl);
        centr = 0.5*(br + bl);
        fix = centr - a;

        mut_D = centr + hlgth;
        fval[0]  = 0.5 * (*fcn)(&mut_D) * pow(fix + hlgth, alfa);
        fval[12] = (*fcn)(&centr) * pow(fix, alfa);
        mut_D = centr - hlgth;
        fval[24] = 0.5 * (*fcn)(&mut_D) * pow(fix - hlgth, alfa);

        for (i = 1; i < 12; i++)
        {
            u = hlgth * x[i-1];
            isym = 24 - i;

            mut_D = centr + u;
            fval[i] = (*fcn)(&mut_D) * pow(fix + u, alfa);
            mut_D = centr - u;
            fval[isym] = (*fcn)(&mut_D) * pow(fix - u, alfa);
        }
        // 150
        factor = pow(hlgth, beta+1.0);
        *result = 0.0;
        *abserr = 0.0;
        res12 = 0.0;
        res24 = 0.0;
        if ((integr != 2) && (integr != 4))
        {
            // integr = 1 or 3
            dqcheb(x,fval,cheb12,cheb24);

            for (i = 0; i < 13; i++)
            {
                res12 = res12 + cheb12[i] * rj[i];
            }
            // 160
            for (i = 0; i < 25; i++)
            {
                res24 = res24 + cheb24[i] * rj[i];
            }
            // 170

            if (integr == 1)
            {
                // 260 -> 270
                *result = (*result + res24) * factor;
                *abserr = (*abserr + fabs(res24-res12)) * factor;
                return;
            }

            // integr == 3
            dc = log(br - bl);
            *result = res24 * dc;
            *abserr = fabs((res24 - res12) * dc);
            res12 = 0.0;
            res24 = 0.0;
            for ( i = 0; i < 13; i++)
            {
                res12 = res12 + cheb12[i] * rh[i];
            }
            // 180
            for (i = 0; i < 25; i++)
            {
                res24 = res24 + cheb24[i] * rh[i];
            }
            // 190
            *result = (*result + res24) * factor;
            *abserr = (*abserr + fabs(res24-res12)) * factor;
            return;
        }
        // 200

        // compute the chebyshev series expansion of the following function
        // f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
        fval[0] = fval[0]*log(hlgth + fix);
        fval[12] = fval[12]*log(fix);
        fval[24] = fval[24]*log(fix - hlgth);
        for (i = 1; i < 12; i++)
        {
            u = hlgth * x[i-1];
            isym = 24 - i;

            fval[i] = fval[i]*log(fix + u);
            fval[isym] = fval[isym]*log(fix - u);
        }
        // 210
        dqcheb(x,fval,cheb12,cheb24);

        //  integr = 2 or 4
        for (i = 0 ; i < 13; i++)
        {
            res12 = res12 + cheb12[i] * rj[i];
        }
        // 220
        for (i = 0; i < 25; i++)
        {
            res24 = res24 + cheb24[i] * rj[i];
        }
        // 230

        if (integr == 2)
        {
            // 260 -> 270
            *result = (*result + res24) * factor;
            *abserr = (*abserr + fabs(res24-res12)) * factor;
            return;
        }

        dc = log(br - bl);
        *result = res24 * dc;
        *abserr = fabs((res24-res12) * dc);
        res12 = 0.0;
        res24 = 0.0;

        // integr == 4
        for (i = 0; i < 13; i++)
        {
            res12 = res12 + cheb12[i] * rh[i];
        }
        // 240
        for (i = 0; i < 25; i++)
        {
            res24 = res24 + cheb24[i] * rh[i];
        }
        // 250

        *result = (*result + res24) * factor;
        *abserr = (*abserr + fabs(res24-res12)) * factor;
        return;
    }

    // If a > bl and b < br, apply the 15-point gauss-kronrod scheme.
    dqk15w(fcn, dqwgts, a, b, alfa, beta, integr, bl, br, result, abserr, &resabs, resasc);
    *nev = 15;
    return;

}


void
dqcheb(const double* x, double* fval, double* cheb12, double* cheb24)
{
    // ***begin prologue  dqcheb
    // ***refer to  dqc25c,dqc25f,dqc25s
    // ***routines called  (none)
    // ***revision date  830518   (yymmdd)
    // ***keywords  chebyshev series expansion, fast fourier transform
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  this routine computes the chebyshev series expansion
    //             of degrees 12 and 24 of a function using a
    //             fast fourier transform method
    //             f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
    //             f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
    //             where t(k,x) is the chebyshev polynomial of degree k.
    // ***description
    //
    //         chebyshev series expansion
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //           on entry
    //            x      - double precision
    //                     vector of dimension 11 containing the
    //                     values cos(k*pi/24), k = 1, ..., 11
    //
    //            fval   - double precision
    //                     vector of dimension 25 containing the
    //                     function values at the points
    //                     (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
    //                     where (a,b) is the approximation interval.
    //                     fval(1) and fval(25) are divided by two
    //                     (these values are destroyed at output).
    //
    //           on return
    //            cheb12 - double precision
    //                     vector of dimension 13 containing the
    //                     chebyshev coefficients for degree 12
    //
    //            cheb24 - double precision
    //                     vector of dimension 25 containing the
    //                     chebyshev coefficients for degree 24
    //
    // ***end prologue  dqcheb
    //
    int i, j;
    double alam, alam1, alam2, part1, part2, part3;
    double v[12];

    // x[11], fval[25], cheb12[13], cheb24[25]

    for (i = 0; i < 12; i++)
    {
        j = 24 - i;
        v[i] = fval[i] - fval[j];
        fval[i] = fval[i] + fval[j];
    }
    // 10

    alam1 = v[0] - v[8];
    alam2 = x[5] * (v[2] - v[6] - v[10]);
    cheb12[3] = alam1 + alam2;
    cheb12[9] = alam1 - alam2;

    alam1 = v[1] - v[7] - v[9];
    alam2 = v[3] - v[5] - v[11];
    alam = x[2]*alam1 + x[8]*alam2;
    cheb24[3] = cheb12[3] + alam;
    cheb24[21] = cheb12[3] - alam;

    alam = x[8]*alam1 - x[2]*alam2;
    cheb24[9] = cheb12[9] + alam;
    cheb24[15] = cheb12[9] - alam;

    part1 = x[3] * v[4];
    part2 = x[7] * v[8];
    part3 = x[5] * v[6];

    alam1 = v[0] + part1 + part2;
    alam2 = x[1]*v[2] + part3 + x[9]*v[10];
    cheb12[1] = alam1 + alam2;
    cheb12[11] = alam1 - alam2;

    alam = x[0]*v[1] + x[2]*v[3] + x[4]*v[5] + x[6]*v[7] + x[8]*v[9] + x[10]*v[11];
    cheb24[1] = cheb12[1] + alam;
    cheb24[23] = cheb12[1] - alam;

    alam = x[10]*v[1] - x[8]*v[3] + x[6]*v[5] - x[4]*v[7] + x[2]*v[9] - x[0]*v[11];
    cheb24[11] = cheb12[11] + alam;
    cheb24[13] = cheb12[11] - alam;

    alam1 = v[0] - part1 + part2;
    alam2 = x[9] * v[2] - part3 + x[1] * v[10];
    cheb12[5] = alam1 + alam2;
    cheb12[7] = alam1 - alam2;

    alam = x[4]*v[1] - x[8]*v[3] - x[0]*v[5] - x[10]*v[7] + x[2]*v[9] + x[6]*v[11];
    cheb24[5] = cheb12[5] + alam;
    cheb24[19] = cheb12[5] - alam;

    alam = x[6]*v[1] - x[2]*v[3] - x[10]*v[5] + x[0]*v[7] - x[8]*v[9] - x[4]*v[11];
    cheb24[7] = cheb12[7] + alam;
    cheb24[17] = cheb12[7] - alam;

    for (i = 0; i < 6; i++)
    {
        j = 12 - i;
        v[i] = fval[i] - fval[j];
        fval[i] = fval[i] + fval[j];
    }
    // 20

    alam1 = v[0] + x[7]*v[4];
    alam2 = x[3]*v[2];
    cheb12[2] = alam1 + alam2;
    cheb12[10] = alam1 - alam2;
    cheb12[6] = v[0] - v[4];

    alam = x[1]*v[1] + x[5]*v[3] + x[9]*v[5];
    cheb24[2] = cheb12[2] + alam;
    cheb24[22] = cheb12[2] - alam;

    alam = x[5]*(v[1] - v[3] - v[5]);
    cheb24[6] = cheb12[6] + alam;
    cheb24[18] = cheb12[6] - alam;

    alam = x[9]*v[1] - x[5]*v[3] + x[1]*v[5];
    cheb24[10] = cheb12[10] + alam;
    cheb24[14] = cheb12[10] - alam;

    for (i = 0; i < 3; i++) {
        j = 6 - i;
        v[i] = fval[i] -fval[j];
        fval[i] = fval[i] + fval[j];
    }
    // 30

    cheb12[4] = v[0] + x[7] * v[2];
    cheb12[8] = fval[0] - x[7] * fval[2];

    alam = x[3] * v[1];
    cheb24[4] = cheb12[4] + alam;
    cheb24[20] = cheb12[4] - alam;

    alam = x[7] * fval[1] - fval[3];
    cheb24[8] = cheb12[8] + alam;
    cheb24[16] = cheb12[8] - alam;
    cheb12[0] = fval[0] + fval[2];

    alam = fval[1] + fval[3];
    cheb24[0] = cheb12[0] + alam;
    cheb24[24] = cheb12[0] - alam;
    cheb12[12] = v[0] - v[2];
    cheb24[12] = cheb12[12];

    alam = 1.0 / 6.0;
    for (i = 1; i < 12; i++)
    {
        cheb12[i] = cheb12[i] * alam;
    }
    // 40

    alam = 0.5 * alam;
    cheb12[0] = cheb12[0] * alam ;
    cheb12[12] = cheb12[12] * alam;
    for (i = 1; i < 24; i ++)
    {
        cheb24[i] = cheb24[i]*alam;
    }
    // 50

    cheb24[0] = 0.5*alam*cheb24[0];
    cheb24[24] = 0.5*alam*cheb24[24];

    return;
}


void
dqelg(int* n, double* epstab, double* result, double* abserr, double* res3la, int* nres)
{
    // ***begin prologue  dqelg
    // ***refer to  dqagie,dqagoe,dqagpe,dqagse
    // ***routines called  d1mach
    // ***revision date  830518   (yymmdd)
    // ***keywords  epsilon algorithm, convergence acceleration,
    //              extrapolation
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math & progr. div. - k.u.leuven
    // ***purpose  the routine determines the limit of a given sequence of
    //             approximations, by means of the epsilon algorithm of
    //             p.wynn. an estimate of the absolute error is also given.
    //             the condensed epsilon table is computed. only those
    //             elements needed for the computation of the next diagonal
    //             are preserved.
    // ***description
    //
    //            epsilon algorithm
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //               n      - integer
    //                        epstab(n) contains the new element in the
    //                        first column of the epsilon table.
    //
    //               epstab - double precision
    //                        vector of dimension 52 containing the elements
    //                        of the two lower diagonals of the triangular
    //                        epsilon table. the elements are numbered
    //                        starting at the right-hand corner of the
    //                        triangle.
    //
    //               result - double precision
    //                        resulting approximation to the integral
    //
    //               abserr - double precision
    //                        estimate of the absolute error computed from
    //                        result and the 3 previous results
    //
    //               res3la - double precision
    //                        vector of dimension 3 containing the last 3
    //                        results
    //
    //               nres   - integer
    //                        number of calls to the routine
    //                        (should be zero at first call)
    //
    // ***end prologue  dqelg
    //
    //            list of major variables
    //            -----------------------
    //
    //            e0     - the 4 elements on which the computation of a new
    //            e1       element in the epsilon table is based
    //            e2
    //            e3                 e0
    //                         e3    e1    new
    //                               e2
    //            newelm - number of elements to be computed in the new
    //                     diagonal
    //            error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
    //            result - the element in the new diagonal with least value
    //                     of error
    //
    //            machine dependent constants
    //            ---------------------------
    //
    //            epmach is the largest relative spacing.
    //            oflow is the largest positive magnitude.
    //            limexp is the maximum number of elements the epsilon
    //            table can contain. if this number is reached, the upper
    //            diagonal of the epsilon table is deleted.
    //
    int i, indx, k1, k2, k3, limexp, n2_rem, newelm, starting_n;
    double delta1, delta2, delta3, epsinf, error, err1, err2, err3;
    double e0, e1, e1abs, e2, e3, res, ss, tol1, tol2, tol3;

    // n is an index of epstab hence must be 0-indexed but the algorithm here
    // is tricky in terms of jump conditions. Hence np1 = n + 1 will be used for
    // those parts that require 1-indexed logic.

    starting_n = *n;
    *nres += 1;
    *abserr = oflow;
    *result = epstab[*n];
    if (*n < 2)
    {
        *abserr = fmax(*abserr, 5.0*epmach*fabs(*result));
        return;
    }

    limexp = 49; // indices 50 and 51 for the lower diagonals
    epstab[*n + 2] = epstab[*n];
    newelm = *n / 2;
    epstab[*n] = oflow;
    k1 = *n;
    for (i = 0; i < newelm; i++)
    {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1 + 2];
        e0 = epstab[k3];
        e1 = epstab[k2];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = fmax(fabs(e2), e1abs)*epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = fmax(e1abs, fabs(e0))*epmach;
        if (!((err2 > tol2) || (err3 > tol3)))
        {
            // If e0, e1 and e2 are equal to within machine accuracy,
            // convergence is assumed.
            // result = e2
            // abserr = abs(e1-e0)+abs(e2-e1)
            *result = res;
            *abserr = err2 + err3;
            *abserr = fmax(*abserr, 5.0*epmach*fabs(*result));
            return;
        }
        // 10

        e3 = epstab[k1];
        epstab[k1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = fmax(e1abs, fabs(e3))*epmach;

        // if two elements are very close to each other, omit a part of the
        // table by adjusting the value of n.
        if ((err1 <= tol1) || (err2 <= tol2) || (err3 <= tol3))
        {
            *n = i + i;
            break;  // Goto 50
        }
        ss = 1.0/delta1 + 1.0/delta2 - 1.0/delta3;
        epsinf = fabs(ss*e1);

        // Test to detect irregular behaviour in the table, and  eventually omit
        // a part of the table adjusting the value of n.
        if (!(epsinf > 1e-4))
        {
            *n = i + i;
            break;  // Goto 50
        }
        // 30

        // Compute a new element and eventually adjust the value of result.
        res = e1 + 1.0 / ss;
        epstab[k1] = res;
        k1 -= 2;
        error = err2 + fabs(res - e2) + err3;
        if (!(error > *abserr))
        {
            *abserr = error;
            *result = res;
        }
    }
    //40

    // Shift the table
    // 50
    if (*n == limexp) { *n = 2*(limexp / 2); }

    n2_rem = (starting_n % 2);
    for (i = 0; i <= newelm; i++)
    {
        epstab[2*i + n2_rem] = epstab[2*i + 2 + n2_rem];
    }
    // 60

    if (*n != starting_n)
    {
        indx = starting_n - *n;
        for (i = 0; i <= *n; i++)
        {
            epstab[i] = epstab[indx];
            indx++;
        }
        // 70
    }
    // 80
    if (*nres < 4)
    {
        res3la[*nres - 1] = *result;
        *abserr = oflow;
    } else {
        // 90
        *abserr = fabs(*result-res3la[2]) + fabs(*result-res3la[1]) + fabs(*result-res3la[0]);
        res3la[0] = res3la[1];
        res3la[1] = res3la[2];
        res3la[2] = *result;
    }
    // 100

    *abserr = fmax(*abserr, 5.0*epmach*fabs(*result));

    return;
}


void
dqk15(double(*fcn)(double* x), const double a, const double b,
      double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk15
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a2
    // ***keywords  15-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div - k.u.leuven
    // ***purpose  to compute i = integral of f over (a,b), with error
    //                            estimate
    //                        j = integral of abs(f) over (a,b)
    // ***description
    //
    //            integration rules
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //             on entry
    //               f      - double precision
    //                        function subprogram defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the calling program.
    //
    //               a      - double precision
    //                        lower limit of integration
    //
    //               b      - double precision
    //                        upper limit of integration
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 15-point
    //                        kronrod rule (resk) obtained by optimal addition
    //                        of abscissae to the7-point gauss rule(resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the absolute error,
    //                        which should not exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral j
    //
    //               resasc - double precision
    //                        approximation to the integral of abs(f-i/(b-a))
    //                        over (a,b)
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk15
    //
    //            the abscissae and weights are given for the interval (-1,1).
    //            because of symmetry only the positive abscissae and their
    //            corresponding weights are given.
    //
    //            xgk    - abscissae of the 15-point kronrod rule
    //                     xgk(2), xgk(4), ...  abscissae of the 7-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ...  abscissae which are optimally
    //                     added to the 7-point gauss rule
    //
    //            wgk    - weights of the 15-point kronrod rule
    //
    //            wg     - weights of the 7-point gauss rule
    //
    //
    //  gauss quadrature weights and kronron quadrature abscissae and weights
    //  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    //  bell labs, nov. 1981.
    //
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk,reskh;
    double fv1[7], fv2[7];
    double mut_D;

static const double wg[4] = {
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327
};
static const double xgk[8] = {
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.000000000000000000000000000000000
};
static const double wgk[8] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714
};

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 15-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc = (*fcn)(&centr);
    resg = wg[3];
    resk = wgk[7]*fc;
    *resabs = fabs(resk);

    for (j = 0; j < 3; j++) {
        jtw = 2*j + 1;
        absc = hlgth * xgk[jtw];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10

    for (j = 0; j < 4; j++) {
        jtwm1 = 2 * j;
        absc = hlgth * xgk[jtwm1];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk*0.5;
    *resasc = wgk[7]*fabs(fc - reskh);
    for (j = 0; j < 7; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}


void
dqk15i(double(*fcn)(double* x), const double boun, const int inf, const double a, const double b,
       double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk15i
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a3a2,h2a4a2
    // ***keywords  15-point transformed gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  the original (infinite integration range is mapped
    //             onto the interval (0,1) and (a,b) is a part of (0,1).
    //             it is the purpose to compute
    //             i = integral of transformed integrand over (a,b),
    //             j = integral of abs(transformed integrand) over (a,b).
    // ***description
    //
    //            integration rule
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //             on entry
    //               f      - double precision
    //                        function subprogram defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the calling program.
    //
    //               boun   - double precision
    //                        finite bound of original integration
    //                        range (set to zero if inf = +2)
    //
    //               inf    - integer
    //                        if inf = -1, the original interval is
    //                                    (-infinity,bound),
    //                        if inf = +1, the original interval is
    //                                    (bound,+infinity),
    //                        if inf = +2, the original interval is
    //                                    (-infinity,+infinity) and
    //                        the integral is computed as the sum of two
    //                        integrals, one over (-infinity,0) and one over
    //                        (0,+infinity).
    //
    //               a      - double precision
    //                        lower limit for integration over subrange
    //                        of (0,1)
    //
    //               b      - double precision
    //                        upper limit for integration over subrange
    //                        of (0,1)
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 15-point
    //                        kronrod rule(resk) obtained by optimal addition
    //                        of abscissae to the 7-point gauss rule(resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the absolute error,
    //                        which should equal or exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral j
    //
    //               resasc - double precision
    //                        approximation to the integral of
    //                        abs((transformed integrand)-i/(b-a)) over (a,b)
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk15i
    //
    //            the abscissae and weights are supplied for the interval
    //            (-1,1).  because of symmetry only the positive abscissae and
    //            their corresponding weights are given.
    //
    //            xgk    - abscissae of the 15-point kronrod rule
    //                     xgk(2), xgk(4), ... abscissae of the 7-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ...  abscissae which are optimally
    //                     added to the 7-point gauss rule
    //
    //            wgk    - weights of the 15-point kronrod rule
    //
    //            wg     - weights of the 7-point gauss rule, corresponding
    //                     to the abscissae xgk(2), xgk(4), ...
    //                     wg(1), wg(3), ... are set to zero.
    //
    int j;
    double absc, absc1, absc2, centr, dinf, fc, fsum, fval1, fval2, hlgth, resg;
    double resk, reskh, tabsc1, tabsc2;
    double fv1[7], fv2[7];

    static const double wg[8] = {
        0.000000000000000000000000000000000,
        0.129484966168869693270611432679082,
        0.000000000000000000000000000000000,
        0.279705391489276667901467771423780,
        0.000000000000000000000000000000000,
        0.381830050505118944950369775488975,
        0.000000000000000000000000000000000,
        0.417959183673469387755102040816327
    };
    static const double xgk[8] = {
        0.991455371120812639206854697526329,
        0.949107912342758524526189684047851,
        0.864864423359769072789712788640926,
        0.741531185599394439863864773280788,
        0.586087235467691130294144838258730,
        0.405845151377397166906606412076961,
        0.207784955007898467600689403773245,
        0.000000000000000000000000000000000
    };
    static const double wgk[8] = {
        0.022935322010529224963732008058970,
        0.063092092629978553290700663189204,
        0.104790010322250183839876322541518,
        0.140653259715525918745189590510238,
        0.169004726639267902826583426598550,
        0.190350578064785409913256402421014,
        0.204432940075298892414161999234649,
        0.209482141084727828012999174891714
    };

    dinf = (inf > 1 ? 1 : inf);
    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    tabsc1 = boun + dinf*(1.0 - centr)/centr;
    fval1 = (*fcn)(&tabsc1);
    if (inf == 2)
    {
        tabsc1 = -tabsc1;
        fval1 = fval1 + (*fcn)(&tabsc1);
        tabsc1 = -tabsc1;
    }
    fc = (fval1 / centr) / centr;

    // Compute the 15-point kronrod approximation to the integral, and
    // estimate the absolute error.

    resg = fc * wg[7];
    resk = fc * wgk[7];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        absc = hlgth * xgk[j];
        absc1 = centr - absc;
        absc2 = centr + absc;
        tabsc1 = boun + dinf*(1.0 - absc1)/absc1;
        tabsc2 = boun + dinf*(1.0 - absc2)/absc2;

        fval1 = (*fcn)(&tabsc1);
        fval2 = (*fcn)(&tabsc2);

        if (inf == 2) {
            tabsc1 = -tabsc1;
            tabsc2 = -tabsc2;
            fval1 = fval1 + (*fcn)(&tabsc1);
            fval2 = fval2 + (*fcn)(&tabsc2);
            tabsc1 = -tabsc1;
            tabsc2 = -tabsc2;
        }
        fval1 = (fval1 / absc1) / absc1;
        fval2 = (fval2 / absc2) / absc2;
        fv1[j] = fval1;
        fv2[j] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j]*fsum;
        resk = resk + wgk[j]*fsum;
        *resabs = *resabs + wgk[j]*(fabs(fval1) + fabs(fval2));
    }
    // 10

    reskh = resk*0.5;
    *resasc = wgk[7]*fabs(fc - reskh);
    for (j = 0; j < 7; j++)
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*hlgth;
    *resasc = (*resasc)*hlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}

void
dqk15w(double(*fcn)(double* x), quadpack_w_func w, const double p1, const double p2,
       const double p3, const double p4, const int kp, const double a, const double b,
       double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk15w
    // ***date written   810101   (yymmdd)
    // ***revision date  830518   (mmddyy)
    // ***category no.  h2a2a2
    // ***keywords  15-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f*w over (a,b), with error
    //                            estimate
    //                        j = integral of abs(f*w) over (a,b)
    // ***description
    //
    //            integration rules
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //              on entry
    //               f      - double precision
    //                        function subprogram defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the driver program.
    //
    //               w      - double precision
    //                        function subprogram defining the integrand
    //                        weight function w(x). the actual name for w
    //                        needs to be declared e x t e r n a l in the
    //                        calling program.
    //
    //               p1, p2, p3, p4 - double precision
    //                        parameters in the weight function
    //
    //               kp     - integer
    //                        key for indicating the type of weight function
    //
    //               a      - double precision
    //                        lower limit of integration
    //
    //               b      - double precision
    //                        upper limit of integration
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 15-point
    //                        kronrod rule (resk) obtained by optimal addition
    //                        of abscissae to the 7-point gauss rule (resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the absolute error,
    //                        which should equal or exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral of abs(f)
    //
    //               resasc - double precision
    //                        approximation to the integral of abs(f-i/(b-a))
    //
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk15w
    //
    //            the abscissae and weights are given for the interval (-1,1).
    //            because of symmetry only the positive abscissae and their
    //            corresponding weights are given.
    //
    //            xgk    - abscissae of the 15-point gauss-kronrod rule
    //                     xgk(2), xgk(4), ... abscissae of the 7-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ... abscissae which are optimally
    //                     added to the 7-point gauss rule
    //
    //            wgk    - weights of the 15-point gauss-kronrod rule
    //
    //            wg     - weights of the 7-point gauss rule
    //
    int j, jtw, jtwm1;
    double absc, absc1, absc2, centr, dhlgth, fc, fsum, fval1, fval2, hlgth;
    double resg,resk,reskh;
    double fv1[7], fv2[7];

    static const double wg[4] = {
        0.1294849661688697,
        0.2797053914892767,
        0.3818300505051189,
        0.4179591836734694
    };
    static const double xgk[8] = {
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993944,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0000000000000000
    };
    static const double wgk[8] = {
        0.2293532201052922e-01,
        0.6309209262997855e-01,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278
    };

    centr = 0.5 * (b + a);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 15-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc = (*fcn)(&centr) * (*w)(centr, p1, p2, p3, p4, kp);
    resg = wg[3] * fc;
    resk = wgk[7] * fc;
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2*j + 1;
        absc = hlgth * xgk[jtw];
        absc1 = centr - absc;
        absc2 = centr + absc;

        fval1 = (*fcn)(&absc1) * (*w)(absc1, p1, p2, p3, p4, kp);
        fval2 = (*fcn)(&absc2) * (*w)(absc2, p1, p2, p3, p4, kp);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10
    for (j = 0; j < 4; j++) {
        jtwm1 = 2*j;
        absc = hlgth * xgk[jtwm1];
        absc1 = centr - absc;
        absc2 = centr + absc;

        fval1 = (*fcn)(&absc1) * (*w)(absc1, p1, p2, p3, p4, kp);
        fval2 = (*fcn)(&absc2) * (*w)(absc2, p1, p2, p3, p4, kp);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk*0.5;
    *resasc = wgk[7]*fabs(fc - reskh);
    for (j = 0; j < 7; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;

}


void
dqk21(double(*fcn)(double* x), const double a, const double b,
      double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk21
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a2
    // ***keywords  21-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f over (a,b), with error
    //                            estimate
    //                        j = integral of abs(f) over (a,b)
    // ***description
    //
    //            integration rules
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //             on entry
    //               f      - double precision
    //                        function subprogram defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the driver program.
    //
    //               a      - double precision
    //                        lower limit of integration
    //
    //               b      - double precision
    //                        upper limit of integration
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 21-point
    //                        kronrod rule (resk) obtained by optimal addition
    //                        of abscissae to the 10-point gauss rule (resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the absolute error,
    //                        which should not exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral j
    //
    //               resasc - double precision
    //                        approximation to the integral of abs(f-i/(b-a))
    //                        over (a,b)
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk21
    //
    //            the abscissae and weights are given for the interval (-1,1).
    //            because of symmetry only the positive abscissae and their
    //            corresponding weights are given.
    //
    //            xgk    - abscissae of the 21-point kronrod rule
    //                     xgk(2), xgk(4), ...  abscissae of the 10-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ...  abscissae which are optimally
    //                     added to the 10-point gauss rule
    //
    //            wgk    - weights of the 21-point kronrod rule
    //
    //            wg     - weights of the 10-point gauss rule
    //
    //
    //  gauss quadrature weights and kronron quadrature abscissae and weights
    //  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    //  bell labs, nov. 1981.
    //
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[10], fv2[10];
    double mut_D;

    static const double wg[5] = {
        0.066671344308688137593568809893332,
        0.149451349150580593145776339657697,
        0.219086362515982043995534934228163,
        0.269266719309996355091226921569469,
        0.295524224714752870173892994651338
    };
    static const double xgk[11] = {
        0.995657163025808080735527280689003,
        0.973906528517171720077964012084452,
        0.930157491355708226001207180059508,
        0.865063366688984510732096688423493,
        0.780817726586416897063717578345042,
        0.679409568299024406234327365114874,
        0.562757134668604683339000099272694,
        0.433395394129247190799265943165784,
        0.294392862701460198131126603103866,
        0.148874338981631210884826001129720,
        0.000000000000000000000000000000000
    };
    static const double wgk[11] = {
        0.011694638867371874278064396062192,
        0.032558162307964727478818972459390,
        0.054755896574351996031381300244580,
        0.075039674810919952767043140916190,
        0.093125454583697605535065465083366,
        0.109387158802297641899210590325805,
        0.123491976262065851077958109831074,
        0.134709217311473325928054001771707,
        0.142775938577060080797094273138717,
        0.147739104901338491374841515972068,
        0.149445554002916905664936468389821
    };

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 21-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc=(*fcn)(&centr);
    resg = 0.0;
    resk = wgk[10]*fc;
    *resabs = fabs(resk);

    for (j = 0; j < 5; j++) {
        jtw = 2*j + 1;
        absc = hlgth * xgk[jtw];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10

    for (j = 0; j < 5; j++) {
        jtwm1 = 2 * j;
        absc = hlgth * xgk[jtwm1];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk * 0.5;
    *resasc = wgk[10]*fabs(fc - reskh);
    for (j = 0; j < 10; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}


void
dqk31(double(*fcn)(double* x), const double a, const double b,
      double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk31
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a2
    // ***keywords  31-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f over (a,b) with error
    //                            estimate
    //                        j = integral of abs(f) over (a,b)
    // ***description
    //
    //            integration rules
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //             on entry
    //               f      - double precision
    //                        function subprogram defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the calling program.
    //
    //               a      - double precision
    //                        lower limit of integration
    //
    //               b      - double precision
    //                        upper limit of integration
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 31-point
    //                        gauss-kronrod rule (resk), obtained by optimal
    //                        addition of abscissae to the 15-point gauss
    //                        rule (resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the modulus,
    //                        which should not exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral j
    //
    //               resasc - double precision
    //                        approximation to the integral of abs(f-i/(b-a))
    //                        over (a,b)
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk31
    //
    //            the abscissae and weights are given for the interval (-1,1).
    //            because of symmetry only the positive abscissae and their
    //            corresponding weights are given.
    //
    //            xgk    - abscissae of the 31-point kronrod rule
    //                     xgk(2), xgk(4), ...  abscissae of the 15-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ...  abscissae which are optimally
    //                     added to the 15-point gauss rule
    //
    //            wgk    - weights of the 31-point kronrod rule
    //
    //            wg     - weights of the 15-point gauss rule
    //
    //
    //  gauss quadrature weights and kronron quadrature abscissae and weights
    //  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    //  bell labs, nov. 1981.
    //
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk,reskh;
    double fv1[15], fv2[15];
    double mut_D;

    static const double wg[8] = {
        0.030753241996117268354628393577204,
        0.070366047488108124709267416450667,
        0.107159220467171935011869546685869,
        0.139570677926154314447804794511028,
        0.166269205816993933553200860481209,
        0.186161000015562211026800561866423,
        0.198431485327111576456118326443839,
        0.202578241925561272880620199967519
    };
    static const double xgk[16] = {
        0.998002298693397060285172840152271,
        0.987992518020485428489565718586613,
        0.967739075679139134257347978784337,
        0.937273392400705904307758947710209,
        0.897264532344081900882509656454496,
        0.848206583410427216200648320774217,
        0.790418501442465932967649294817947,
        0.724417731360170047416186054613938,
        0.650996741297416970533735895313275,
        0.570972172608538847537226737253911,
        0.485081863640239680693655740232351,
        0.394151347077563369897207370981045,
        0.299180007153168812166780024266389,
        0.201194093997434522300628303394596,
        0.101142066918717499027074231447392,
        0.000000000000000000000000000000000
    };
    static const double wgk[16] = {
        0.005377479872923348987792051430128,
        0.015007947329316122538374763075807,
        0.025460847326715320186874001019653,
        0.035346360791375846222037948478360,
        0.044589751324764876608227299373280,
        0.053481524690928087265343147239430,
        0.062009567800670640285139230960803,
        0.069854121318728258709520077099147,
        0.076849680757720378894432777482659,
        0.083080502823133021038289247286104,
        0.088564443056211770647275443693774,
        0.093126598170825321225486872747346,
        0.096642726983623678505179907627589,
        0.099173598721791959332393173484603,
        0.100769845523875595044946662617570,
        0.101330007014791549017374792767493
    };

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 31-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc=(*fcn)(&centr);
    resg = wg[7]*fc;
    resk = wgk[15]*fc;
    *resabs = fabs(resk);

    for (j = 0; j < 7; j++) {
        jtw = 2*j + 1;
        absc = hlgth * xgk[jtw];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10

    for (j = 0; j < 8; j++) {
        jtwm1 = 2 * j;
        absc = hlgth * xgk[jtwm1];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk * 0.5;
    *resasc = wgk[15]*fabs(fc - reskh);
    for (j = 0; j < 15; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}


void
dqk41(double(*fcn)(double* x), const double a, const double b,
      double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk41
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a2
    // ***keywords  41-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f over (a,b), with error
    //                            estimate
    //                        j = integral of abs(f) over (a,b)
    // ***description
    //
    //            integration rules
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //             on entry
    //               f      - double precision
    //                        function subprogram defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the calling program.
    //
    //               a      - double precision
    //                        lower limit of integration
    //
    //               b      - double precision
    //                        upper limit of integration
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 41-point
    //                        gauss-kronrod rule (resk) obtained by optimal
    //                        addition of abscissae to the 20-point gauss
    //                        rule (resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the absolute error,
    //                        which should not exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral j
    //
    //               resasc - double precision
    //                        approximation to the integal of abs(f-i/(b-a))
    //                        over (a,b)
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk41
    //
    //            the abscissae and weights are given for the interval (-1,1).
    //            because of symmetry only the positive abscissae and their
    //            corresponding weights are given.
    //
    //            xgk    - abscissae of the 41-point gauss-kronrod rule
    //                     xgk(2), xgk(4), ...  abscissae of the 20-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ...  abscissae which are optimally
    //                     added to the 20-point gauss rule
    //
    //            wgk    - weights of the 41-point gauss-kronrod rule
    //
    //            wg     - weights of the 20-point gauss rule
    //
    //
    //  gauss quadrature weights and kronron quadrature abscissae and weights
    //  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    //  bell labs, nov. 1981.
    //
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[20], fv2[20];
    double mut_D;

    static const double wg[10] = {
        0.017614007139152118311861962351853,
        0.040601429800386941331039952274932,
        0.062672048334109063569506535187042,
        0.083276741576704748724758143222046,
        0.101930119817240435036750135480350,
        0.118194531961518417312377377711382,
        0.131688638449176626898494499748163,
        0.142096109318382051329298325067165,
        0.149172986472603746787828737001969,
        0.152753387130725850698084331955098
    };
    static const double xgk[21] = {
        0.998859031588277663838315576545863,
        0.993128599185094924786122388471320,
        0.981507877450250259193342994720217,
        0.963971927277913791267666131197277,
        0.940822633831754753519982722212443,
        0.912234428251325905867752441203298,
        0.878276811252281976077442995113078,
        0.839116971822218823394529061701521,
        0.795041428837551198350638833272788,
        0.746331906460150792614305070355642,
        0.693237656334751384805490711845932,
        0.636053680726515025452836696226286,
        0.575140446819710315342946036586425,
        0.510867001950827098004364050955251,
        0.443593175238725103199992213492640,
        0.373706088715419560672548177024927,
        0.301627868114913004320555356858592,
        0.227785851141645078080496195368575,
        0.152605465240922675505220241022678,
        0.076526521133497333754640409398838,
        0.000000000000000000000000000000000
    };
    static const double wgk[21] = {
        0.003073583718520531501218293246031,
        0.008600269855642942198661787950102,
        0.014626169256971252983787960308868,
        0.020388373461266523598010231432755,
        0.025882133604951158834505067096153,
        0.031287306777032798958543119323801,
        0.036600169758200798030557240707211,
        0.041668873327973686263788305936895,
        0.046434821867497674720231880926108,
        0.050944573923728691932707670050345,
        0.055195105348285994744832372419777,
        0.059111400880639572374967220648594,
        0.062653237554781168025870122174255,
        0.065834597133618422111563556969398,
        0.068648672928521619345623411885368,
        0.071054423553444068305790361723210,
        0.073030690332786667495189417658913,
        0.074582875400499188986581418362488,
        0.075704497684556674659542775376617,
        0.076377867672080736705502835038061,
        0.076600711917999656445049901530102
    };

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 41-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc=(*fcn)(&centr);
    resg = 0.0;
    resk = wgk[20]*fc;
    *resabs = fabs(resk);

    for (j = 0; j < 10; j++) {
        jtw = 2*j + 1;
        absc = hlgth * xgk[jtw];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10

    for (j = 0; j < 10; j++) {
        jtwm1 = 2 * j;
        absc = hlgth * xgk[jtwm1];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk * 0.5;
    *resasc = wgk[20]*fabs(fc - reskh);
    for (j = 0; j < 20; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}


void
dqk51(double(*fcn)(double* x), const double a, const double b,
      double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk51
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a2
    // ***keywords  51-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f over (a,b) with error
    //                            estimate
    //                        j = integral of abs(f) over (a,b)
    // ***description
    //
    //            integration rules
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters
    //             on entry
    //               f      - double precision
    //                        function subroutine defining the integrand
    //                        function f(x). the actual name for f needs to be
    //                        declared e x t e r n a l in the calling program.
    //
    //               a      - double precision
    //                        lower limit of integration
    //
    //               b      - double precision
    //                        upper limit of integration
    //
    //             on return
    //               result - double precision
    //                        approximation to the integral i
    //                        result is computed by applying the 51-point
    //                        kronrod rule (resk) obtained by optimal addition
    //                        of abscissae to the 25-point gauss rule (resg).
    //
    //               abserr - double precision
    //                        estimate of the modulus of the absolute error,
    //                        which should not exceed abs(i-result)
    //
    //               resabs - double precision
    //                        approximation to the integral j
    //
    //               resasc - double precision
    //                        approximation to the integral of abs(f-i/(b-a))
    //                        over (a,b)
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk51
    //
    //            the abscissae and weights are given for the interval (-1,1).
    //            because of symmetry only the positive abscissae and their
    //            corresponding weights are given.
    //
    //            xgk    - abscissae of the 51-point kronrod rule
    //                     xgk(2), xgk(4), ...  abscissae of the 25-point
    //                     gauss rule
    //                     xgk(1), xgk(3), ...  abscissae which are optimally
    //                     added to the 25-point gauss rule
    //
    //            wgk    - weights of the 51-point kronrod rule
    //
    //            wg     - weights of the 25-point gauss rule
    //
    //
    //  gauss quadrature weights and kronron quadrature abscissae and weights
    //  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    //  bell labs, nov. 1981.
    //
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[25], fv2[25];
    double mut_D;

    static const double wg[13] = {
        0.011393798501026287947902964113235,
        0.026354986615032137261901815295299,
        0.040939156701306312655623487711646,
        0.054904695975835191925936891540473,
        0.068038333812356917207187185656708,
        0.080140700335001018013234959669111,
        0.091028261982963649811497220702892,
        0.100535949067050644202206890392686,
        0.108519624474263653116093957050117,
        0.114858259145711648339325545869556,
        0.119455763535784772228178126512901,
        0.122242442990310041688959518945852,
        0.123176053726715451203902873079050
    };
    static const double xgk[26] = {
        0.999262104992609834193457486540341,
        0.995556969790498097908784946893902,
        0.988035794534077247637331014577406,
        0.976663921459517511498315386479594,
        0.961614986425842512418130033660167,
        0.942974571228974339414011169658471,
        0.920747115281701561746346084546331,
        0.894991997878275368851042006782805,
        0.865847065293275595448996969588340,
        0.833442628760834001421021108693570,
        0.797873797998500059410410904994307,
        0.759259263037357630577282865204361,
        0.717766406813084388186654079773298,
        0.673566368473468364485120633247622,
        0.626810099010317412788122681624518,
        0.577662930241222967723689841612654,
        0.526325284334719182599623778158010,
        0.473002731445714960522182115009192,
        0.417885382193037748851814394594572,
        0.361172305809387837735821730127641,
        0.303089538931107830167478909980339,
        0.243866883720988432045190362797452,
        0.183718939421048892015969888759528,
        0.122864692610710396387359818808037,
        0.061544483005685078886546392366797,
        0.000000000000000000000000000000000
    };
    static const double wgk[26] = {
        0.001987383892330315926507851882843,
        0.005561932135356713758040236901066,
        0.009473973386174151607207710523655,
        0.013236229195571674813656405846976,
        0.016847817709128298231516667536336,
        0.020435371145882835456568292235939,
        0.024009945606953216220092489164881,
        0.027475317587851737802948455517811,
        0.030792300167387488891109020215229,
        0.034002130274329337836748795229551,
        0.037116271483415543560330625367620,
        0.040083825504032382074839284467076,
        0.042872845020170049476895792439495,
        0.045502913049921788909870584752660,
        0.047982537138836713906392255756915,
        0.050277679080715671963325259433440,
        0.052362885806407475864366712137873,
        0.054251129888545490144543370459876,
        0.055950811220412317308240686382747,
        0.057437116361567832853582693939506,
        0.058689680022394207961974175856788,
        0.059720340324174059979099291932562,
        0.060539455376045862945360267517565,
        0.061128509717053048305859030416293,
        0.061471189871425316661544131965264,
        0.061580818067832935078759824240066
    };

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 51-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc = (*fcn)(&centr);
    resg = wg[12]*fc;
    resk = wgk[25]*fc;
    *resabs = fabs(resk);

    for (j = 0; j < 12; j++) {
        jtw = 2*j + 1;
        absc = hlgth * xgk[jtw];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10

    for (j = 0; j < 13; j++) {
        jtwm1 = 2 * j;
        absc = hlgth * xgk[jtwm1];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk * 0.5;
    *resasc = wgk[25]*fabs(fc - reskh);
    for (j = 0; j < 25; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}


void
dqk61(double(*fcn)(double* x), const double a, const double b,
      double* result, double* abserr, double* resabs, double* resasc)
{
    // ***begin prologue  dqk61
    // ***date written   800101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a1a2
    // ***keywords  61-point gauss-kronrod rules
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  to compute i = integral of f over (a,b) with error
    //                            estimate
    //                        j = integral of dabs(f) over (a,b)
    // ***description
    //
    //         integration rule
    //         standard fortran subroutine
    //         double precision version
    //
    //
    //         parameters
    //          on entry
    //            f      - double precision
    //                     function subprogram defining the integrand
    //                     function f(x). the actual name for f needs to be
    //                     declared e x t e r n a l in the calling program.
    //
    //            a      - double precision
    //                     lower limit of integration
    //
    //            b      - double precision
    //                     upper limit of integration
    //
    //          on return
    //            result - double precision
    //                     approximation to the integral i
    //                     result is computed by applying the 61-point
    //                     kronrod rule (resk) obtained by optimal addition of
    //                     abscissae to the 30-point gauss rule (resg).
    //
    //            abserr - double precision
    //                     estimate of the modulus of the absolute error,
    //                     which should equal or exceed dabs(i-result)
    //
    //            resabs - double precision
    //                     approximation to the integral j
    //
    //            resasc - double precision
    //                     approximation to the integral of dabs(f-i/(b-a))
    //
    //
    // ***references  (none)
    // ***routines called  d1mach
    // ***end prologue  dqk61
    //
    //            the abscissae and weights are given for the
    //            interval (-1,1). because of symmetry only the positive
    //            abscissae and their corresponding weights are given.
    //
    //            xgk   - abscissae of the 61-point kronrod rule
    //                    xgk(2), xgk(4)  ... abscissae of the 30-point
    //                    gauss rule
    //                    xgk(1), xgk(3)  ... optimally added abscissae
    //                    to the 30-point gauss rule
    //
    //            wgk   - weights of the 61-point kronrod rule
    //
    //            wg    - weigths of the 30-point gauss rule
    //
    //
    //  gauss quadrature weights and kronron quadrature abscissae and weights
    //  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    //  bell labs, nov. 1981.
    //
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[30], fv2[30];
    double mut_D;

    static const double wg[15] = {
        0.007968192496166605615465883474674,
        0.018466468311090959142302131912047,
        0.028784707883323369349719179611292,
        0.038799192569627049596801936446348,
        0.048402672830594052902938140422808,
        0.057493156217619066481721689402056,
        0.065974229882180495128128515115962,
        0.073755974737705206268243850022191,
        0.080755895229420215354694938460530,
        0.086899787201082979802387530715126,
        0.092122522237786128717632707087619,
        0.096368737174644259639468626351810,
        0.099593420586795267062780282103569,
        0.101762389748405504596428952168554,
        0.102852652893558840341285636705415
    };
    static const double xgk[31] = {
        0.999484410050490637571325895705811,
        0.996893484074649540271630050918695,
        0.991630996870404594858628366109486,
        0.983668123279747209970032581605663,
        0.973116322501126268374693868423707,
        0.960021864968307512216871025581798,
        0.944374444748559979415831324037439,
        0.926200047429274325879324277080474,
        0.905573307699907798546522558925958,
        0.882560535792052681543116462530226,
        0.857205233546061098958658510658944,
        0.829565762382768397442898119732502,
        0.799727835821839083013668942322683,
        0.767777432104826194917977340974503,
        0.733790062453226804726171131369528,
        0.697850494793315796932292388026640,
        0.660061064126626961370053668149271,
        0.620526182989242861140477556431189,
        0.579345235826361691756024932172540,
        0.536624148142019899264169793311073,
        0.492480467861778574993693061207709,
        0.447033769538089176780609900322854,
        0.400401254830394392535476211542661,
        0.352704725530878113471037207089374,
        0.304073202273625077372677107199257,
        0.254636926167889846439805129817805,
        0.204525116682309891438957671002025,
        0.153869913608583546963794672743256,
        0.102806937966737030147096751318001,
        0.051471842555317695833025213166723,
        0.000000000000000000000000000000000
    };
    static const double wgk[31] = {
        0.001389013698677007624551591226760,
        0.003890461127099884051267201844516,
        0.006630703915931292173319826369750,
        0.009273279659517763428441146892024,
        0.011823015253496341742232898853251,
        0.014369729507045804812451432443580,
        0.016920889189053272627572289420322,
        0.019414141193942381173408951050128,
        0.021828035821609192297167485738339,
        0.024191162078080601365686370725232,
        0.026509954882333101610601709335075,
        0.028754048765041292843978785354334,
        0.030907257562387762472884252943092,
        0.032981447057483726031814191016854,
        0.034979338028060024137499670731468,
        0.036882364651821229223911065617136,
        0.038678945624727592950348651532281,
        0.040374538951535959111995279752468,
        0.041969810215164246147147541285970,
        0.043452539701356069316831728117073,
        0.044814800133162663192355551616723,
        0.046059238271006988116271735559374,
        0.047185546569299153945261478181099,
        0.048185861757087129140779492298305,
        0.049055434555029778887528165367238,
        0.049795683427074206357811569379942,
        0.050405921402782346840893085653585,
        0.050881795898749606492297473049805,
        0.051221547849258772170656282604944,
        0.051426128537459025933862879215781,
        0.051494729429451567558340433647099
    };

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    // Compute the 61-point kronrod approximation to the integral, and
    // estimate the absolute error.

    fc=(*fcn)(&centr);
    resg = 0.0;
    resk = wgk[30] * fc;
    *resabs = fabs(resk);

    for (j = 0; j < 15; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg = resg + wg[j] * fsum;
        resk = resk + wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    // 10

    for (j = 0; j < 15; j++) {
        jtwm1 = 2 * j;
        absc = hlgth * xgk[jtwm1];

        mut_D = centr - absc;
        fval1 = (*fcn)(&mut_D);
        mut_D = centr + absc;
        fval2 = (*fcn)(&mut_D);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    // 15

    reskh = resk * 0.5;
    *resasc = wgk[30] * fabs(fc - reskh);
    for (j = 0; j < 30; j++ )
    {
        *resasc = *resasc + wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));
    }
    // 20

    *result = resk * hlgth;
    *resabs = (*resabs)*dhlgth;
    *resasc = (*resasc)*dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if ((*resasc != 0.0) && (*abserr != 0.0))
    {
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr)/(*resasc)), 1.5));
    }
    if (*resabs > uflow/(50.0 * epmach))
    {
        *abserr = fmax(epmach * 50.0 * (*resabs),(*abserr));
    }

    return;
}


void
dqmomo(const double alfa, const double beta, double* ri, double* rj, double* rg,
       double* rh, const int integr)
{
    // ***begin prologue  dqmomo
    // ***date written   820101   (yymmdd)
    // ***revision date  830518   (yymmdd)
    // ***category no.  h2a2a1,c3a2
    // ***keywords  modified chebyshev moments
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  this routine computes modified chebsyshev moments. the k-th
    //             modified chebyshev moment is defined as the integral over
    //             (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
    //             polynomial of degree k.
    // ***description
    //
    //         modified chebyshev moments
    //         standard fortran subroutine
    //         double precision version
    //
    //         parameters
    //            alfa   - double precision
    //                     parameter in the weight function w(x), alfa.gt.(-1)
    //
    //            beta   - double precision
    //                     parameter in the weight function w(x), beta.gt.(-1)
    //
    //            ri     - double precision
    //                     vector of dimension 25
    //                     ri(k) is the integral over (-1,1) of
    //                     (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
    //
    //            rj     - double precision
    //                     vector of dimension 25
    //                     rj(k) is the integral over (-1,1) of
    //                     (1-x)**beta*t(k-1,x), k = 1, ..., 25.
    //
    //            rg     - double precision
    //                     vector of dimension 25
    //                     rg(k) is the integral over (-1,1) of
    //                     (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
    //
    //            rh     - double precision
    //                     vector of dimension 25
    //                     rh(k) is the integral over (-1,1) of
    //                     (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
    //
    //            integr - integer
    //                     input parameter indicating the modified
    //                     moments to be computed
    //                     integr = 1 compute ri, rj
    //                            = 2 compute ri, rj, rg
    //                            = 3 compute ri, rj, rh
    //                            = 4 compute ri, rj, rg, rh
    //
    // ***references  (none)
    // ***routines called  (none)
    // ***end prologue  dqmomo
    //
    int i;
    double alfp1, an, betp1, alfp2, betp2, ralf, rbet;

    alfp1 = alfa + 1.0;
    betp1 = beta + 1.0;
    alfp2 = alfa + 2.0;
    betp2 = beta + 2.0;
    ralf = pow(2.0, alfp1);
    rbet = pow(2.0, betp1);

    // Compute ri, rj using a forward recurrence relation.
    ri[0] = ralf / alfp1;
    rj[0] = rbet / betp1;
    ri[1] = ri[0] * alfa / alfp2;
    rj[1] = rj[0] * beta / betp2;

    an = 2.0;
    for (i = 2; i < 25; i++)
    {
        ri[i] = -(ralf + an*(an-alfp2)*ri[i-1])/((an - 1.0)*(an + alfp1));
        rj[i] = -(rbet + an*(an-betp2)*rj[i-1])/((an - 1.0)*(an + betp1));
        an = an + 1.0;
    }
    if (integr == 1)
    {
        for (i = 1; i < 25; i += 2){ rj[i] = -rj[i]; }
        return;
    }

    if (integr != 3)
    {
        // integr 2 or 4
        // Compute rg using a forware recurrence relation.
        rg[0] = -ri[0] / alfp1;
        rg[1] = -(ralf + ralf) / (alfp2 * alfp2) - rg[0];
        an = 2.0;
        for (i = 2; i < 25; i++)
        {
            rg[i] = -(an*(an - alfp2)*rg[i-1] - an*ri[i-1] + (an - 1.0)*ri[i]) / ((an - 1.0)*(an + alfp1));
            an = an + 1.0;
        }
        if(integr == 2)
        {
            for (i = 1; i < 25; i += 2) { rj[i] = -rj[i]; }
            return;
        }
    }
    // 40

    // Compute rh using a forward recurrence relation.
    rh[0] = -rj[0] / betp1;
    rh[1] = -(rbet + rbet) / (betp2 * betp2) - rh[0];
    an = 2.0;
    for (i = 2; i < 25; i++)
    {
        rh[i] = -(an*(an - betp2)*rh[i-1] - an*rj[i-1] + (an - 1.0)*rj[i]) / ((an - 1.0)*(an + betp1));
        an = an + 1.0;
    }
    for (i = 1; i < 25; i += 2) { rh[i] = -rh[i]; }
    // 60
    for (i = 1; i < 25; i += 2) { rj[i] = -rj[i]; }
    // 80

    return;
}


void
dqng(double(*fcn)(double* x), const double a, const double b,
     const double epsabs, const double epsrel,
     double* result, double* abserr, int* neval, int* ier)
{
    // ***begin prologue  dqng
    // ***date written   800101   (yymmdd)
    // ***revision date  810101   (yymmdd)
    // ***category no.  h2a1a1
    // ***keywords  automatic integrator, smooth integrand,
    //              non-adaptive, gauss-kronrod(patterson)
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl math & progr. div. - k.u.leuven
    //            kahaner,david,nbs - modified (2/82)
    // ***purpose  the routine calculates an approximation result to a
    //             given definite integral i = integral of f over (a,b),
    //             hopefully satisfying following claim for accuracy
    //             abs(i-result).le.max(epsabs,epsrel*abs(i)).
    // ***description
    //
    //  non-adaptive integration
    //  standard fortran subroutine
    //  double precision version
    //
    //            f      - double precision
    //                     function subprogram defining the integrand function
    //                     f(x). the actual name for f needs to be declared
    //                     e x t e r n a l in the driver program.
    //
    //            a      - double precision
    //                     lower limit of integration
    //
    //            b      - double precision
    //                     upper limit of integration
    //
    //            epsabs - double precision
    //                     absolute accuracy requested
    //            epsrel - double precision
    //                     relative accuracy requested
    //                     if  epsabs.le.0
    //                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
    //                     the routine will end with ier = 6.
    //
    //          on return
    //            result - double precision
    //                     approximation to the integral i
    //                     result is obtained by applying the 21-point
    //                     gauss-kronrod rule (res21) obtained by optimal
    //                     addition of abscissae to the 10-point gauss rule
    //                     (res10), or by applying the 43-point rule (res43)
    //                     obtained by optimal addition of abscissae to the
    //                     21-point gauss-kronrod rule, or by applying the
    //                     87-point rule (res87) obtained by optimal addition
    //                     of abscissae to the 43-point rule.
    //
    //            abserr - double precision
    //                     estimate of the modulus of the absolute error,
    //                     which should equal or exceed abs(i-result)
    //
    //            neval  - integer
    //                     number of integrand evaluations
    //
    //            ier    - ier = 0 normal and reliable termination of the
    //                             routine. it is assumed that the requested
    //                             accuracy has been achieved.
    //                     ier.gt.0 abnormal termination of the routine. it is
    //                             assumed that the requested accuracy has
    //                             not been achieved.
    //            error messages
    //                     ier = 1 the maximum number of steps has been
    //                             executed. the integral is probably too
    //                             difficult to be calculated by dqng.
    //                         = 6 the input is invalid, because
    //                             epsabs.le.0 and
    //                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
    //                             result, abserr and neval are set to zero.
    //
    // ***references  (none)
    // ***routines called  d1mach,xerror
    // ***end prologue  dqng
    //
    //            the following data statements contain the
    //            abscissae and weights of the integration rules used.
    //
    //            x1      abscissae common to the 10-, 21-, 43- and 87-
    //                    point rule
    //            x2      abscissae common to the 21-, 43- and 87-point rule
    //            x3      abscissae common to the 43- and 87-point rule
    //            x4      abscissae of the 87-point rule
    //            w10     weights of the 10-point formula
    //            w21a    weights of the 21-point formula for abscissae x1
    //            w21b    weights of the 21-point formula for abscissae x2
    //            w43a    weights of the 43-point formula for abscissae x1, x3
    //            w43b    weights of the 43-point formula for abscissae x3
    //            w87a    weights of the 87-point formula for abscissae x1,
    //                    x2, x3
    //            w87b    weights of the 87-point formula for abscissae x4
    //
    //
    //  gauss-kronrod-patterson quadrature coefficients for use in
    //  quadpack routine qng.  these coefficients were calculated with
    //  101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
    //
    int ipx, k, l;
    double fv1[5], fv2[5], fv3[5], fv4[5], savfun[21];
    double absc, centr, dhlgth, fcentr, fval, fval1, fval2, hlgth;
    double res10,res21,res43,res87, resabs,resasc,reskh;

    static const double x1[5] = {
        0.973906528517171720077964012084452,
        0.865063366688984510732096688423493,
        0.679409568299024406234327365114874,
        0.433395394129247190799265943165784,
        0.148874338981631210884826001129720
    };
    static const double w10[5] = {
        0.066671344308688137593568809893332,
        0.149451349150580593145776339657697,
        0.219086362515982043995534934228163,
        0.269266719309996355091226921569469,
        0.295524224714752870173892994651338
    };
    static const double x2[5] = {
        0.995657163025808080735527280689003,
        0.930157491355708226001207180059508,
        0.780817726586416897063717578345042,
        0.562757134668604683339000099272694,
        0.294392862701460198131126603103866
    };
    static const double w21a[5] = {
        0.032558162307964727478818972459390,
        0.075039674810919952767043140916190,
        0.109387158802297641899210590325805,
        0.134709217311473325928054001771707,
        0.147739104901338491374841515972068
    };
    static const double w21b[6] = {
        0.011694638867371874278064396062192,
        0.054755896574351996031381300244580,
        0.093125454583697605535065465083366,
        0.123491976262065851077958109831074,
        0.142775938577060080797094273138717,
        0.149445554002916905664936468389821
    };
    static const double x3[11] = {
        0.999333360901932081394099323919911,
        0.987433402908088869795961478381209,
        0.954807934814266299257919200290473,
        0.900148695748328293625099494069092,
        0.825198314983114150847066732588520,
        0.732148388989304982612354848755461,
        0.622847970537725238641159120344323,
        0.499479574071056499952214885499755,
        0.364901661346580768043989548502644,
        0.222254919776601296498260928066212,
        0.074650617461383322043914435796506
    };
    static const double w43a[10] = {
        0.016296734289666564924281974617663,
        0.037522876120869501461613795898115,
        0.054694902058255442147212685465005,
        0.067355414609478086075553166302174,
        0.073870199632393953432140695251367,
        0.005768556059769796184184327908655,
        0.027371890593248842081276069289151,
        0.046560826910428830743339154433824,
        0.061744995201442564496240336030883,
        0.071387267268693397768559114425516
    };
    static const double w43b[12] = {
        0.001844477640212414100389106552965,
        0.010798689585891651740465406741293,
        0.021895363867795428102523123075149,
        0.032597463975345689443882222526137,
        0.042163137935191811847627924327955,
        0.050741939600184577780189020092084,
        0.058379395542619248375475369330206,
        0.064746404951445885544689259517511,
        0.069566197912356484528633315038405,
        0.072824441471833208150939535192842,
        0.074507751014175118273571813842889,
        0.074722147517403005594425168280423
    };
    static const double x4[22] = {
        0.999902977262729234490529830591582,
        0.997989895986678745427496322365960,
        0.992175497860687222808523352251425,
        0.981358163572712773571916941623894,
        0.965057623858384619128284110607926,
        0.943167613133670596816416634507426,
        0.915806414685507209591826430720050,
        0.883221657771316501372117548744163,
        0.845710748462415666605902011504855,
        0.803557658035230982788739474980964,
        0.757005730685495558328942793432020,
        0.706273209787321819824094274740840,
        0.651589466501177922534422205016736,
        0.593223374057961088875273770349144,
        0.531493605970831932285268948562671,
        0.466763623042022844871966781659270,
        0.399424847859218804732101665817923,
        0.329874877106188288265053371824597,
        0.258503559202161551802280975429025,
        0.185695396568346652015917141167606,
        0.111842213179907468172398359241362,
        0.037352123394619870814998165437704
    };
    static const double w87a[21] = {
        0.008148377384149172900002878448190,
        0.018761438201562822243935059003794,
        0.027347451050052286161582829741283,
        0.033677707311637930046581056957588,
        0.036935099820427907614589586742499,
        0.002884872430211530501334156248695,
        0.013685946022712701888950035273128,
        0.023280413502888311123409291030404,
        0.030872497611713358675466394126442,
        0.035693633639418770719351355457044,
        0.000915283345202241360843392549948,
        0.005399280219300471367738743391053,
        0.010947679601118931134327826856808,
        0.016298731696787335262665703223280,
        0.021081568889203835112433060188190,
        0.025370969769253827243467999831710,
        0.029189697756475752501446154084920,
        0.032373202467202789685788194889595,
        0.034783098950365142750781997949596,
        0.036412220731351787562801163687577,
        0.037253875503047708539592001191226
    };
    static const double w87b[23] = {
        0.000274145563762072350016527092881,
        0.001807124155057942948341311753254,
        0.004096869282759164864458070683480,
        0.006758290051847378699816577897424,
        0.009549957672201646536053581325377,
        0.012329447652244853694626639963780,
        0.015010447346388952376697286041943,
        0.017548967986243191099665352925900,
        0.019938037786440888202278192730714,
        0.022194935961012286796332102959499,
        0.024339147126000805470360647041454,
        0.026374505414839207241503786552615,
        0.028286910788771200659968002987960,
        0.030052581128092695322521110347341,
        0.031646751371439929404586051078883,
        0.033050413419978503290785944862689,
        0.034255099704226061787082821046821,
        0.035262412660156681033782717998428,
        0.036076989622888701185500318003895,
        0.036698604498456094498018047441094,
        0.037120549269832576114119958413599,
        0.037334228751935040321235449094698,
        0.037361073762679023410321241766599
    };
    double mut_D;
    *result = 0.0;
    *abserr = 0.0;
    *neval = 0;
    *ier = 6;
    if ((epsabs <= 0.0) && (epsrel < fmax(50.0*epmach, 0.5e-28))) { return; }

    hlgth = 0.5*(b - a);
    dhlgth = fabs(hlgth);
    centr = 0.5*(b + a);
    fcentr = (*fcn)(&centr);
    *neval = 21;
    *ier = 1;

    for (l = 1; l <=3; l++) {
        switch (l) {
            case 1:
                // Compute the integral using the 10- and 21-point formula.

                res10 = 0.0;
                res21 = w21b[5] * fcentr;
                resabs = w21b[5] * fabs(fcentr);
                for (k = 0;k < 5; k++) {
                    absc = hlgth * x1[k];

                    mut_D = centr + absc;
                    fval1 = (*fcn)(&mut_D);
                    mut_D = centr - absc;
                    fval2 = (*fcn)(&mut_D);

                    fval = fval1 + fval2;
                    res10 = res10 + w10[k] * fval;
                    res21 = res21 + w21a[k] * fval;
                    resabs = resabs + w21a[k] * (fabs(fval1) + fabs(fval2));
                    savfun[k] = fval;
                    fv1[k] = fval1;
                    fv2[k] = fval2;
                }
                // 10
                ipx = 4;
                for (k = 0; k < 5; k++) {
                    ipx++;
                    absc = hlgth * x2[k];

                    mut_D = centr + absc;
                    fval1 = (*fcn)(&mut_D);
                    mut_D = centr - absc;
                    fval2 = (*fcn)(&mut_D);

                    fval = fval1 + fval2;
                    res21 = res21 + w21b[k] * fval;
                    resabs = resabs + w21b[k] * (fabs(fval1) + fabs(fval2));
                    savfun[ipx] = fval;
                    fv3[k] = fval1;
                    fv4[k] = fval2;
                }
                // 15

                // Test for convergence.
                *result = res21 * hlgth;
                resabs = resabs*dhlgth;
                reskh = 0.5 * res21;
                resasc = w21b[5] * fabs(fcentr - reskh);
                for (k = 0; k < 5; k++)
                {
                 resasc = resasc + w21a[k]*(fabs(fv1[k] -reskh) + fabs(fv2[k] - reskh))
                           + w21b[k] * (fabs(fv3[k] - reskh) + fabs(fv4[k]-reskh));
                }
                *abserr = fabs((res21 - res10) * hlgth);
                resasc = resasc*dhlgth;
                break;
            case 2:
                // Compute the integral using the 43-point formula.

                res43 = w43b[11] * fcentr;
                *neval = 43;
                for (k = 0; k < 10; k++)
                {
                 res43 = res43 + savfun[k] * w43a[k];
                }
                // 30
                for (k = 0; k < 11; k++) {
                    ipx++;
                    absc = hlgth * x3[k];

                    mut_D = centr + absc;
                    fval1 = (*fcn)(&mut_D);
                    mut_D = centr - absc;
                    fval2 = (*fcn)(&mut_D);

                    fval = fval1 + fval2;
                    res43 = res43 + fval * w43b[k];
                    savfun[ipx] = fval;
                }
                //40

                // Test for convergence.
                *result = res43 * hlgth;
                *abserr = fabs((res43-res21) * hlgth);
                break;
            case 3:
                // Compute the integral using the 87-point formula.

                res87 = w87b[22] * fcentr;
                *neval = 87;
                for (k = 0; k < 21; k++)
                {
                    res87 = res87 + (savfun[k] * w87a[k]);
                }
                // 50
                for (k = 0; k < 22; k++)
                {
                    absc = hlgth * x4[k];

                    mut_D = centr + absc;
                    fval1 = (*fcn)(&mut_D);
                    mut_D = centr - absc;
                    fval2 = (*fcn)(&mut_D);

                    fval = fval1 + fval2;
                    res87 = res87 + w87b[k] * fval;
                }
                // 60

                // Test for convergence.
                *result = res87 * hlgth;
                *abserr = fabs((res87 - res43) * hlgth);
                break;
        }
        // 65

        if ((resasc != 0.0) && (*abserr != 0.0))
        {
            *abserr = resasc*fmin(1.0, pow(200.0*(*abserr)/resasc, 1.5));
        }
        if (resabs > uflow/(50.0 * epmach))
        {
            *abserr = fmax((epmach * 50.0) * resabs, *abserr);
        }
        if (*abserr <= fmax(epsabs, epsrel*fabs(*result))) {
            *ier = 0;
            break;
        }
    }
    // 70

    // 999
    return;
}


void
dqpsrt(const int limit, const int last, int* maxerr, double* ermax,
       const double* elist, int *iord, int* nrmax)
{
    // ***begin prologue  dqpsrt
    // ***refer to  dqage,dqagie,dqagpe,dqawse
    // ***routines called  (none)
    // ***revision date  810101   (yymmdd)
    // ***keywords  sequential sorting
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  this routine maintains the descending ordering in the
    //             list of the local error estimated resulting from the
    //             interval subdivision process. at each call two error
    //             estimates are inserted using the sequential search
    //             method, top-down for the largest error estimate and
    //             bottom-up for the smallest error estimate.
    // ***description
    //
    //            ordering routine
    //            standard fortran subroutine
    //            double precision version
    //
    //            parameters (meaning at output)
    //               limit  - integer
    //                        maximum number of error estimates the list
    //                        can contain
    //
    //               last   - integer
    //                        number of error estimates currently in the list
    //
    //               maxerr - integer
    //                        maxerr points to the nrmax-th largest error
    //                        estimate currently in the list
    //
    //               ermax  - double precision
    //                        nrmax-th largest error estimate
    //                        ermax = elist(maxerr)
    //
    //               elist  - double precision
    //                        vector of dimension last containing
    //                        the error estimates
    //
    //               iord   - integer
    //                        vector of dimension last, the first k elements
    //                        of which contain pointers to the error
    //                        estimates, such that
    //                        elist(iord(1)),...,  elist(iord(k))
    //                        form a decreasing sequence, with
    //                        k = last if last.le.(limit/2+2), and
    //                        k = limit+1-last otherwise
    //
    //               nrmax  - integer
    //                        maxerr = iord(nrmax)
    //
    // ***end prologue  dqpsrt
    //
    int i, ibeg, j, jbnd, jupbn, k;
    double errmin, errmax;
    // Index vars maxerr, iord, nrmax

    if (last <= 2)
    {
        iord[0] = 0;
        iord[1] = 1;
        *maxerr = iord[*nrmax];
        *ermax = elist[*maxerr];
        return;
    }

    // This part of the routine is only executed if, due to a difficult integrand,
    // subdivision increased the error estimate. In the normal case the insert
    // procedure should start after the nrmax-th largest error estimate.
    errmax = elist[*maxerr];
    while ((*nrmax > 0) && (errmax > elist[iord[*nrmax-1]]))
    {
        iord[*nrmax] = iord[*nrmax - 1];
        *nrmax -= 1;
    }

    // Compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    jupbn = (last > limit/2 + 2 ? limit - last + 2 : last - 1);
    errmin = elist[last - 1];

    // insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    jbnd = jupbn-1;
    ibeg = *nrmax + 1;

    if (ibeg <= jbnd)
    {
        for (i = ibeg; i <= jbnd; i++)
        {
            if (errmax >= elist[iord[i]]) { goto LINE60; }
            iord[i-1] = iord[i];
        }
        // 40
    }
    // 50
    iord[jbnd] = *maxerr;
    iord[jupbn] = last-1;
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return;
LINE60:
    // 60
    // Insert errmin by traversing the list bottom-up
    iord[i-1] = *maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; j++)
    {
        if (errmin < elist[iord[k]]) { break; }
        iord[k+1] = iord[k];
        k--;
    }
    //70

    if (errmin < elist[iord[k]])
    {
        iord[k+1] = last-1;
    } else {
        iord[i] = last-1;
    }

    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];

    return;
}


double
dqwgtc(const double x, const double c, const double p2, const double p3, const double p4, const int integr)
{
    // ***begin prologue  dqwgtc
    // ***refer to dqk15w
    // ***routines called  (none)
    // ***revision date  810101   (yymmdd)
    // ***keywords  weight function, cauchy principal value
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  this function subprogram is used together with the
    //             routine qawc and defines the weight function.
    // ***end prologue  dqwgtc
    //
    return 1.0 / (x - c);
}


double
dqwgtf(const double x, const double omega, const double p2, const double p3, const double p4, const int integr)
{
    // ***begin prologue  dqwgtf
    // ***refer to   dqk15w
    // ***routines called  (none)
    // ***revision date 810101   (yymmdd)
    // ***keywords  cos or sin in weight function
    // ***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. * progr. div. - k.u.leuven
    // ***end prologue  dqwgtf
    //
    double omx = omega * x;
    if (integr == 1)
        return cos(omx);
    else
        return sin(omx);
}


double
dqwgts(const double x, const double a, const double b, const double alfa, const double beta, const int integr)
{
    // ***begin prologue  dqwgts
    // ***refer to dqk15w
    // ***routines called  (none)
    // ***revision date  810101   (yymmdd)
    // ***keywords  weight function, algebraico-logarithmic
    //              end-point singularities
    // ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    //            de doncker,elise,appl. math. & progr. div. - k.u.leuven
    // ***purpose  this function subprogram is used together with the
    //             routine dqaws and defines the weight function.
    // ***end prologue  dqwgts
    //
    double xma = x - a;
    double bmx = b - x;
    double ret = pow(xma, alfa)*pow(bmx, beta);
    switch (integr)
    {
        case 1:
            break;
        case 2:
            ret = ret*log(xma);
            break;
        case 3:
            ret = ret*log(bmx);
            break;
        case 4:
            ret = ret*log(xma)*log(bmx);
        default:
            break;
    }
    return ret;
}
