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
 * This file is a C translation of the Fortran code written by Jorge J. Moré,
 * Burton S. Garbow, and Kenneth E. Hillstrom, with the original description
 * below. Original Fortran docstrings are included at the top of each function.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Minpack includes software for solving nonlinear equations and nonlinear
 * least squares problems.  Five algorithmic paths each include a core
 * subroutine and an easy-to-use driver.  The algorithms proceed either from
 * an analytic specification of the Jacobian matrix or directly from the
 * problem functions.  The paths include facilities forsystems of equations
 * with a banded Jacobian matrix, for least squares problems with a large
 * amount of data, and for checking the consistency of the Jacobian matrix
 * with the functions.
 *
 *
 * Disclaimer:
 * Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
 *
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 *
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials
 * provided with the distribution.
 *
 * 3. The end-user documentation included with the
 * redistribution, if any, must include the following
 * acknowledgment:
 *
 *    "This product includes software developed by the
 *    University of Chicago, as Operator of Argonne National
 *    Laboratory.
 *
 * Alternately, this acknowledgment may appear in the software
 * itself, if and wherever such third-party acknowledgments
 * normally appear.
 *
 * 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
 * WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
 * UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
 * THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
 * OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
 * OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
 * USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
 * THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
 * DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
 * UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
 * BE CORRECTED.
 *
 * 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
 * HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
 * ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
 * INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
 * ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
 * PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
 * SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
 * (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
 * EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 * POSSIBILITY OF SUCH LOSS OR DAMAGES.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * The original Fortran code can be found at https://www.netlib.org/minpack/
 *
 * References:
 *
 * [1]: J. J. Moré, B. S. Garbow, and K. E. Hillstrom, User Guide for
 *      MINPACK-1, Argonne National Laboratory Report ANL-80-74, Argonne,
 *      Illinois, 1980.
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "__minpack.h"
#include <stdlib.h>
#include <math.h>


// Internal routines
static void dogleg(const int,const double*,const double*,const double*,const double*,double*,double*,double*);
static double enorm(const int,const double*);
static void fdjac1(int(*)(int*,double*,double*,int*),const int,double*,const double*,double*,const int,int*,const int,const int,const double,double*,double*);
static void fdjac2(int(*)(int*,int*,double*,double*,int*),const int,const int,double*,const double*,double*,const int,int*,const double,double*);
static void lmpar(const int,double *,const int,const int*,const double*,const double*,const double,double*,double*,double*,double*,double*);
static void qform(const int,const int,double*,const int,double*);
static void qrfac(const int,const int,double*,const int,const int,int*,double*,double*,double*);
static void qrsolv(const int,double*,const int,const int*,const double*,const double*,double*,double*,double*);
static void r1mpyq(const int,const int,double*,const int,const double*,const double*);
static void r1updt(const int,const int,double*,const double*,double*,double*,int*);
static void rwupdt(const int,double*,const int,const double*,double*,double*,double*,double*);


static const double dpmpar[3] = {
    2.220446049250313e-16,    /* np.finfo(np.float64).eps  */
    2.2250738585072014e-308,  /* np.finfo(np.float64).tiny */
    1.7976931348623157e+308,  /* np.finfo(np.float64).max  */
};


// Exported functions CHKDER, HYBRD, HYBRJ, LMDIF, LMDER, LMSTR
void CHKDER(const int m, const int n, double* x, const double* fvec, const double* fjac,
            const int ldfjac, double* xp, const double* fvecp, const int mode, double* err)
{
    //     **********
    //
    //     subroutine chkder
    //
    //     this subroutine checks the gradients of m nonlinear functions
    //     in n variables, evaluated at a point x, for consistency with
    //     the functions themselves. the user must call chkder twice,
    //     first with mode = 1 and then with mode = 2.
    //
    //     mode = 1. on input, x must contain the point of evaluation.
    //               on output, xp is set to a neighboring point.
    //
    //     mode = 2. on input, fvec must contain the functions and the
    //                         rows of fjac must contain the gradients
    //                         of the respective functions each evaluated
    //                         at x, and fvecp must contain the functions
    //                         evaluated at xp.
    //               on output, err contains measures of correctness of
    //                          the respective gradients.
    //
    //     the subroutine does not perform reliably if cancellation or
    //     rounding errors cause a severe loss of significance in the
    //     evaluation of a function. therefore, none of the components
    //     of x should be unusually small (in particular, zero) or any
    //     other value which may cause loss of significance.
    //
    //     the subroutine statement is
    //
    //       subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
    //
    //     where
    //
    //       m is a positive integer input variable set to the number
    //         of functions.
    //
    //       n is a positive integer input variable set to the number
    //         of variables.
    //
    //       x is an input array of length n.
    //
    //       fvec is an array of length m. on input when mode = 2,
    //         fvec must contain the functions evaluated at x.
    //
    //       fjac is an m by n array. on input when mode = 2,
    //         the rows of fjac must contain the gradients of
    //         the respective functions evaluated at x.
    //
    //       ldfjac is a positive integer input parameter not less than m
    //         which specifies the leading dimension of the array fjac.
    //
    //       xp is an array of length n. on output when mode = 1,
    //         xp is set to a neighboring point of x.
    //
    //       fvecp is an array of length m. on input when mode = 2,
    //         fvecp must contain the functions evaluated at xp.
    //
    //       mode is an integer input variable set to 1 on the first call
    //         and 2 on the second. other values of mode are equivalent
    //         to mode = 1.
    //
    //       err is an array of length m. on output when mode = 2,
    //         err contains measures of correctness of the respective
    //         gradients. if there is no severe loss of significance,
    //         then if err(i) is 1.0 the i-th gradient is correct,
    //         while if err(i) is 0.0 the i-th gradient is incorrect.
    //         for values of err between 0.0 and 1.0, the categorization
    //         is less certain. in general, a value of err(i) greater
    //         than 0.5 indicates that the i-th gradient is probably
    //         correct, while a value of err(i) less than 0.5 indicates
    //         that the i-th gradient is probably incorrect.
    //
    //     subprograms called
    //
    //       minpack supplied ... dpmpar
    //
    //       fortran supplied ... dabs,dlog10,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    double temp;
    const double epsmch = dpmpar[0];
    const double eps = sqrt(epsmch);
    const double factor = 100.0;
    const double epsf = factor*epsmch;
    const double epslog = log10(eps);

    if (mode != 2) {
        for (int j = 0; j < n; j++)
        {
            temp = fmax(eps*fabs(x[j]), eps);
            xp[j] = x[j] + temp;
        }
        // 10
    } else {
        // mode = 2
        for (int i = 0; i < m; i++)
        {
            err[i] = 0.0;
        }
        // 30
        for (int j = 0; j < n; j++)
        {
            temp = fabs(x[j]);
            temp = (temp == 0.0 ? 1.0 : temp);
            for (int i = 0; i < m; i++)
            {
                err[i] += temp*fjac[i + ldfjac*j];
            }
            // 40
        }
        // 50

        for (int i = 0; i < m; i++)
        {
            temp = 1.0;
            if ((fvec[i] != 0.0) && (fvecp[i] != 0.0) &&
                (fabs(fvecp[i] - fvec[i]) >= epsf*fabs(fvec[i])))
            {
                temp = eps*fabs((fvecp[i] - fvec[i])/eps - err[i]) /
                       (fabs(fvec[i]) + fabs(fvecp[i]));
            }
            err[i] = 1.0;
            if ((temp > epsmch) && (temp < eps))
            {
                err[i] = (log10(temp) - epslog)/epslog;
            }
            if (temp >= eps)
            {
                err[i] = 0.0;
            }
        }
        // 60
    }
    // 70

    return;
}


void HYBRD(int(*fcn)(int* n, double* x, double* fvec, int* iflag), const int n,
           double* x, double* fvec, const double xtol, const int maxfev,
           const int ml, const int mu, const double epsfcn, double* diag,
           const int mode, const double factor, const int nprint,
           int* info, int* nfev, double* fjac, const int ldfjac, double* r,
           const int lr, double* qtf, double* wa1, double* wa2, double* wa3,
           double* wa4)
{
    //     **********
    //
    //     subroutine hybrd
    //
    //     the purpose of hybrd is to find a zero of a system of
    //     n nonlinear functions in n variables by a modification
    //     of the powell hybrid method. the user must provide a
    //     subroutine which calculates the functions. the jacobian is
    //     then calculated by a forward-difference approximation.
    //
    //     the subroutine statement is
    //
    //       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
    //                        diag,mode,factor,nprint,info,nfev,fjac,
    //                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions. fcn must be declared
    //         in an external statement in the user calling
    //         program, and should be written as follows.
    //
    //         subroutine fcn(n,x,fvec,iflag)
    //         integer n,iflag
    //         double precision x(n),fvec(n)
    //         ----------
    //         calculate the functions at x and
    //         return this vector in fvec.
    //         ---------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of hybrd.
    //         in this case set iflag to a negative integer.
    //
    //       n is a positive integer input variable set to the number
    //         of functions and variables.
    //
    //       x is an array of length n. on input x must contain
    //         an initial estimate of the solution vector. on output x
    //         contains the final estimate of the solution vector.
    //
    //       fvec is an output array of length n which contains
    //         the functions evaluated at the output x.
    //
    //       xtol is a nonnegative input variable. termination
    //         occurs when the relative error between two consecutive
    //         iterates is at most xtol.
    //
    //       maxfev is a positive integer input variable. termination
    //         occurs when the number of calls to fcn is at least maxfev
    //         by the end of an iteration.
    //
    //       ml is a nonnegative integer input variable which specifies
    //         the number of subdiagonals within the band of the
    //         jacobian matrix. if the jacobian is not banded, set
    //         ml to at least n - 1.
    //
    //       mu is a nonnegative integer input variable which specifies
    //         the number of superdiagonals within the band of the
    //         jacobian matrix. if the jacobian is not banded, set
    //         mu to at least n - 1.
    //
    //       epsfcn is an input variable used in determining a suitable
    //         step length for the forward-difference approximation. this
    //         approximation assumes that the relative errors in the
    //         functions are of the order of epsfcn. if epsfcn is less
    //         than the machine precision, it is assumed that the relative
    //         errors in the functions are of the order of the machine
    //         precision.
    //
    //       diag is an array of length n. if mode = 1 (see
    //         below), diag is internally set. if mode = 2, diag
    //         must contain positive entries that serve as
    //         multiplicative scale factors for the variables.
    //
    //       mode is an integer input variable. if mode = 1, the
    //         variables will be scaled internally. if mode = 2,
    //         the scaling is specified by the input diag. other
    //         values of mode are equivalent to mode = 1.
    //
    //       factor is a positive input variable used in determining the
    //         initial step bound. this bound is set to the product of
    //         factor and the euclidean norm of diag*x if nonzero, or else
    //         to factor itself. in most cases factor should lie in the
    //         interval (.1,100.). 100. is a generally recommended value.
    //
    //       nprint is an integer input variable that enables controlled
    //         printing of iterates if it is positive. in this case,
    //         fcn is called with iflag = 0 at the beginning of the first
    //         iteration and every nprint iterations thereafter and
    //         immediately prior to return, with x and fvec available
    //         for printing. if nprint is not positive, no special calls
    //         of fcn with iflag = 0 are made.
    //
    //       info is an integer output variable. if the user has
    //         terminated execution, info is set to the (negative)
    //         value of iflag. see description of fcn. otherwise,
    //         info is set as follows.
    //
    //         info = 0   improper input parameters.
    //
    //         info = 1   relative error between two consecutive iterates
    //                    is at most xtol.
    //
    //         info = 2   number of calls to fcn has reached or exceeded
    //                    maxfev.
    //
    //         info = 3   xtol is too small. no further improvement in
    //                    the approximate solution x is possible.
    //
    //         info = 4   iteration is not making good progress, as
    //                    measured by the improvement from the last
    //                    five jacobian evaluations.
    //
    //         info = 5   iteration is not making good progress, as
    //                    measured by the improvement from the last
    //                    ten iterations.
    //
    //       nfev is an integer output variable set to the number of
    //         calls to fcn.
    //
    //       fjac is an output n by n array which contains the
    //         orthogonal matrix q produced by the qr factorization
    //         of the final approximate jacobian.
    //
    //       ldfjac is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array fjac.
    //
    //       r is an output array of length lr which contains the
    //         upper triangular matrix produced by the qr factorization
    //         of the final approximate jacobian, stored rowwise.
    //
    //       lr is a positive integer input variable not less than
    //         (n*(n+1))/2.
    //
    //       qtf is an output array of length n which contains
    //         the vector (q transpose)*fvec.
    //
    //       wa1, wa2, wa3, and wa4 are work arrays of length n.
    //
    //     subprograms called
    //
    //       user-supplied ...... fcn
    //
    //       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
    //                            qform,qrfac,r1mpyq,r1updt
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, iflag, iter, j, jeval, sing, l, msum, mut_n, ncfail, ncsuc, nslow1, nslow2;
    double actred, delta, fnorm, fnorm1, pnorm, prered, ratio, ssum, temp, xnorm;
    int iwa[1];
    double epsmch = dpmpar[0];

    mut_n = n;
    *info = 0;
    iflag = 0;
    *nfev = 0;

    // Check the input parameters for errors
    if ((n <= 0) || (xtol < 0.0) || (maxfev <= 0) || (ml < 0) || (mu < 0) ||
	    (factor <= 0.0) || (ldfjac < n) || (lr < n * (n + 1) / 2)) { goto EXIT300; }
    if (mode == 2)
    {
        for (j = 0; j < n; j++)
        {
            if (diag[j] <= 0.0) { goto EXIT300; }
        }
        // 10
    }
    // 20

    // Evaluate the function at the starting point and calculate its norm.
    iflag = 1;
    (*fcn)(&mut_n, x, fvec, &iflag);
    *nfev = 1;
    if (iflag < 0) { goto EXIT300; }
    fnorm = enorm(n, fvec);

    // Determine the number of calls to fcn needed to compute the jacobian.
    msum = (ml+mu+1 > n ? n : ml+mu+1);

    // Initialize iteration counter and monitors.
    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

    // Beginning of the outer loop

    while (1)
    {
        jeval = 1;

        // Calculate the jacobian.
        iflag = 2;
        fdjac1(fcn, n, x, fvec, fjac, ldfjac, &iflag, ml, mu, epsfcn, wa1, wa2);
        *nfev += msum;
        if (iflag < 0) { goto EXIT300; }

        // Compute the qr factorization of the jacobian.
        qrfac(n, n, fjac, ldfjac, 0, iwa, wa1, wa2, wa3);

        // On the first iteration and if mode is 1, scale according to the norms
        // of the columns of the initial jacobian.
        if (iter == 1)
        {
            if (mode != 2)
            {
                for (j = 0; j < n; j++)
                {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0)
                    {
                        diag[j] = 1.0;
                    }
                }
                // 40
            }
            // 50

            // On the first iteration, calculate the norm of the scaled x and
            // initialize the step bound delta.
            for (j = 0; j < n; j++)
            {
                wa3[j] = diag[j]*x[j];
            }
            // 60

            xnorm = enorm(n, wa3);
            delta = factor*xnorm;
            delta = (factor*xnorm == 0.0 ? factor : factor*xnorm);
        }
        // 70

        // Form (q.T)*fvec and store in qtf
        for (i = 0; i < n; i++)
        {
            qtf[i] = fvec[i];
        }
        // 80

        for (j = 0; j < n; j++)
        {
            if (fjac[j + ldfjac*j] != 0.0)
            {
                ssum = 0.0;
                for (i = j; i < n; i++)
                {
                    ssum += fjac[i + ldfjac*j]*qtf[i];
                }
                // 90
                temp = -ssum / fjac[j + ldfjac*j];
                for (i = j; i < n; i++)
                {
                    qtf[i] += fjac[i + ldfjac*j]*temp;
                }
                // 100
            }
            // 110
        }
        // 120

        // Copy the triangular factor of the qr factorization into r.
        sing = 0;
        for (j = 0; j < n; j++)
        {
            l = j;
            if (j > 0)
            {
                for (i = 0; i < j; i++)
                {
                    r[l] = fjac[i + ldfjac*j];
                    l += n - i - 1;
                }
                // 130
            }
            // 140
            r[l] = wa1[j];
            if (wa1[j] == 0.0) { sing = 1; }
        }
        // 150

        // Accumulate the orthogonal factor in fjac.
        qform(n, n, fjac, ldfjac, wa1);

        // Rescale if necessary
        if (mode != 2)
        {
            for (j = 0; j < n; j++)
            {
                diag[j] = fmax(diag[j], wa2[j]);
            }
            // 160
        }
        // 170

        // Beginning of the inner loop.
        while (1)
        {
            if (nprint > 0) {
                iflag = 0;
                if ((iter - 1) % nprint == 0) {
                    (*fcn)(&mut_n, x, fvec, &iflag);
                }
                if (iflag < 0) { goto EXIT300; }
            }
            // 190

            // Determine the direction p.
            dogleg(n, r, diag, qtf, &delta, wa1, wa2, wa3);

            // Store the direction of p.
            for (j = 0; j < n; j++)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            // 200

            pnorm = enorm(n, wa3);

            // On the first iteration, adjust the initial step bound
            if (iter == 1) { delta = fmin(delta, pnorm); }

            // Evaluate the function at x + p and calculate its norm.
            iflag = 1;
            (*fcn)(&mut_n, wa2, wa4, &iflag);
            *nfev += 1;
            if (iflag < 0) { goto EXIT300; }
            fnorm1 = enorm(n, wa4);

            // Compute the scaled actual reduction
            actred = (fnorm1 < fnorm ? 1.0 - pow(fnorm1/fnorm, 2.0) : -1.0);

            // Compute the scaled predicted reduction.
            l = 0;
            for (i = 0; i < n; i++)
            {
                ssum = 0.0;
                for (j = i; j < n; j++)
                {
                    ssum += r[l] * wa1[j];
                    l++;
                }
                // 210
                wa3[i] = qtf[i] + ssum;
            }
            // 220

            temp = enorm(n, wa3);
            prered = (temp < fnorm ? 1.0 - pow(temp/fnorm, 2.0) : 0.0);

            // Compute the ratio of the actual to the predicted reduction.
            ratio = (prered > 0.0 ? actred/prered : 0.0);

            // Update the step bound
            if (ratio < 0.1)
            {
                ncsuc = 0;
                ncfail++;
                delta = 0.5*delta;
            } else {
                // 230
                ncfail = 0;
                ncsuc++;
                if ((ratio >= 0.5) || (ncsuc > 1)) { delta = fmax(delta, pnorm/0.5); }
                if (fabs(ratio - 1.0) <= 0.1) { delta = pnorm / 0.5; }
            }
            // 240

            // Test for successful iteration.
            if (ratio >= 0.0001)
            {
                // Successful iteration. Update x, fvec, and their norms.
                for (j = 0; j < n; j++)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                    fvec[j] = wa4[j];
                }
                // 250
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }
            // 260

            // Determine the progress of the iteration.
            nslow1++;
            if (actred >= 0.001) { nslow1 = 0; }
            if (jeval) { nslow2++; }
            if (actred >= 0.1) { nslow2 = 0; }

            // Test for convergence.
            if ((delta <= xtol*xnorm) || (fnorm == 0.0)) { *info = 1; }
            if (*info != 0) { goto EXIT300; }

            // Tests for termination and stringent tolerances
            if (*nfev >= maxfev) { *info = 2; }
            if (0.1*fmax(0.1*delta, pnorm) <= epsmch*xnorm) { *info = 3; }
            if (nslow2 == 5) { *info = 4; }
            if (nslow1 == 10) { *info = 5; }
            if (*info != 0) { goto EXIT300; }

            // Criterion for recalculating jacobian approximation by forward
            // differences.
            if (ncfail == 2) { break; }

            // Calculate the rank one modification to the jacobian and update
            // qtf if necessary.
            for (j = 0; j < n; j++)
            {
                ssum = 0.0;
                for (i = 0; i < n; i++)
                {
                    ssum += fjac[i + ldfjac*j]*wa4[i];
                }
                // 270
                wa2[j] = (ssum - wa3[j]) / pnorm;
                wa1[j] = diag[j]*((diag[j] * wa1[j]) / pnorm);
                if (ratio >= 0.0001) { qtf[j] = ssum; }
            }
            // 280

            // Compute the qr factorization of the updated jacobian.
            r1updt(n, n, r, wa1, wa2, wa3, &sing);
            r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
            r1mpyq(1, n, qtf, 1, wa2, wa3);

            // End of the inner loop.
            jeval = 0;
            // 290
        }
        // End of the outer loop.
    }
    // 300

    // Termination, either normal or user imposed.
EXIT300:
    if (iflag < 0) { *info = iflag; }
    iflag = 0;
    if (nprint > 0) { (*fcn)(&mut_n, x, fvec, &iflag);}

    return;
}


void HYBRJ(int(*fcn)(int* n, double* x, double* fvec, double* fjac, int* ldfjac,
           int* iflag), const int n, double* x, double* fvec, double* fjac,
           const int ldfjac, const double xtol, const int maxfev, double* diag,
           const int mode, const double factor, const int nprint, int* info,
           int* nfev, int* njev, double* r, const int ldr, double* qtf,
           double* wa1, double* wa2, double* wa3, double* wa4)
{
    //     **********
    //
    //     subroutine hybrj
    //
    //     the purpose of hybrj is to find a zero of a system of
    //     n nonlinear functions in n variables by a modification
    //     of the powell hybrid method. the user must provide a
    //     subroutine which calculates the functions and the jacobian.
    //
    //     the subroutine statement is
    //
    //       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
    //                        mode,factor,nprint,info,nfev,njev,r,lr,qtf,
    //                        wa1,wa2,wa3,wa4)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions and the jacobian. fcn must
    //         be declared in an external statement in the user
    //         calling program, and should be written as follows.
    //
    //         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
    //         integer n,ldfjac,iflag
    //         double precision x(n),fvec(n),fjac(ldfjac,n)
    //         ----------
    //         if iflag = 1 calculate the functions at x and
    //         return this vector in fvec. do not alter fjac.
    //         if iflag = 2 calculate the jacobian at x and
    //         return this matrix in fjac. do not alter fvec.
    //         ---------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of hybrj.
    //         in this case set iflag to a negative integer.
    //
    //       n is a positive integer input variable set to the number
    //         of functions and variables.
    //
    //       x is an array of length n. on input x must contain
    //         an initial estimate of the solution vector. on output x
    //         contains the final estimate of the solution vector.
    //
    //       fvec is an output array of length n which contains
    //         the functions evaluated at the output x.
    //
    //       fjac is an output n by n array which contains the
    //         orthogonal matrix q produced by the qr factorization
    //         of the final approximate jacobian.
    //
    //       ldfjac is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array fjac.
    //
    //       xtol is a nonnegative input variable. termination
    //         occurs when the relative error between two consecutive
    //         iterates is at most xtol.
    //
    //       maxfev is a positive integer input variable. termination
    //         occurs when the number of calls to fcn with iflag = 1
    //         has reached maxfev.
    //
    //       diag is an array of length n. if mode = 1 (see
    //         below), diag is internally set. if mode = 2, diag
    //         must contain positive entries that serve as
    //         multiplicative scale factors for the variables.
    //
    //       mode is an integer input variable. if mode = 1, the
    //         variables will be scaled internally. if mode = 2,
    //         the scaling is specified by the input diag. other
    //         values of mode are equivalent to mode = 1.
    //
    //       factor is a positive input variable used in determining the
    //         initial step bound. this bound is set to the product of
    //         factor and the euclidean norm of diag*x if nonzero, or else
    //         to factor itself. in most cases factor should lie in the
    //         interval (.1,100.). 100. is a generally recommended value.
    //
    //       nprint is an integer input variable that enables controlled
    //         printing of iterates if it is positive. in this case,
    //         fcn is called with iflag = 0 at the beginning of the first
    //         iteration and every nprint iterations thereafter and
    //         immediately prior to return, with x and fvec available
    //         for printing. fvec and fjac should not be altered.
    //         if nprint is not positive, no special calls of fcn
    //         with iflag = 0 are made.
    //
    //       info is an integer output variable. if the user has
    //         terminated execution, info is set to the (negative)
    //         value of iflag. see description of fcn. otherwise,
    //         info is set as follows.
    //
    //         info = 0   improper input parameters.
    //
    //         info = 1   relative error between two consecutive iterates
    //                    is at most xtol.
    //
    //         info = 2   number of calls to fcn with iflag = 1 has
    //                    reached maxfev.
    //
    //         info = 3   xtol is too small. no further improvement in
    //                    the approximate solution x is possible.
    //
    //         info = 4   iteration is not making good progress, as
    //                    measured by the improvement from the last
    //                    five jacobian evaluations.
    //
    //         info = 5   iteration is not making good progress, as
    //                    measured by the improvement from the last
    //                    ten iterations.
    //
    //       nfev is an integer output variable set to the number of
    //         calls to fcn with iflag = 1.
    //
    //       njev is an integer output variable set to the number of
    //         calls to fcn with iflag = 2.
    //
    //       r is an output array of length lr which contains the
    //         upper triangular matrix produced by the qr factorization
    //         of the final approximate jacobian, stored rowwise.
    //
    //       lr is a positive integer input variable not less than
    //         (n*(n+1))/2.
    //
    //       qtf is an output array of length n which contains
    //         the vector (q transpose)*fvec.
    //
    //       wa1, wa2, wa3, and wa4 are work arrays of length n.
    //
    //     subprograms called
    //
    //       user-supplied ...... fcn
    //
    //       minpack-supplied ... dogleg,dpmpar,enorm,
    //                            qform,qrfac,r1mpyq,r1updt
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,mod
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, iflag, iter, j, jeval, l, ncfail, ncsuc, nslow1, nslow2, sing;
    int mut_n, mut_ldf, mut_ldr;
    double actred, delta, fnorm, fnorm1, pnorm, prered, ratio, ssum,temp,xnorm;
    int iwa[1];
    double epsmch = dpmpar[0];

    // To avoid const poisoning use mutable version of inputs
    mut_n = n;
    mut_ldf = ldfjac;
    mut_ldr = ldr;

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

    // Check the input parameters for errors.
    if ((n <= 0) || (ldfjac < n) || (xtol < 0.0) || (maxfev <= 0) ||
        (factor <= 0.0) || (ldr < (n * (n + 1)) / 2)) { goto EXIT300; }

    if (mode == 2) {
        for (j = 0; j < n; j++) {
            if (diag[j] <= 0.0) { goto EXIT300; }
        }
        // 10
    }
    // 20

    // Evaluate the function at the starting point and calculate its norm.
    iflag = 1;
    (*fcn)(&mut_n, x, fvec, fjac, &mut_ldf, &iflag);
    *nfev = 1;
    if (iflag < 0) { goto EXIT300; }
    fnorm = enorm(n, fvec);


    // Initialize iteration counter and monitors.
    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

    // Beginning of the outer loop

    while (1)
    {
        jeval = 1;

        // Calculate the jacobian.
        iflag = 2;
        (*fcn)(&mut_n, x, fvec, fjac, &mut_ldf, &iflag);
        *njev += 1;
        if (iflag < 0) { goto EXIT300; }

        // Compute the qr factorization of the jacobian.
        qrfac(n, n, fjac, ldfjac, 0, iwa, wa1, wa2, wa3);

        // On the first iteration and if mode is 1, scale according to the norms
        // of the columns of the initial jacobian.
        if (iter == 1)
        {
            if (mode != 2)
            {
                for (j = 0; j < n; j++)
                {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0)
                    {
                        diag[j] = 1.0;
                    }
                }
                // 40
            }
            // 50

            // On the first iteration, calculate the norm of the scaled x and
            // initialize the step bound delta.
            for (j = 0; j < n; j++)
            {
                wa3[j] = diag[j]*x[j];
            }
            // 60

            xnorm = enorm(n, wa3);
            delta = factor*xnorm;
            delta = (delta == 0.0 ? factor : delta);
        }
        // 70

        // Form (q.T)*fvec and store in qtf
        for (i = 0; i < n; i++)
        {
            qtf[i] = fvec[i];
        }
        // 80

        for (j = 0; j < n; j++)
        {
            if (fjac[j + ldfjac*j] != 0.0)
            {
                ssum = 0.0;
                for (i = j; i < n; i++)
                {
                    ssum += fjac[i + ldfjac*j]*qtf[i];
                }
                // 90
                temp = -ssum / fjac[j + ldfjac*j];
                for (i = j; i < n; i++)
                {
                    qtf[i] += fjac[i + ldfjac*j]*temp;
                }
                // 100
            }
            // 110
        }
        // 120

        // Copy the triangular factor of the qr factorization into r.
        sing = 0;
        for (j = 0; j < n; j++)
        {
            l = j;
            if (j > 0)
            {
                for (i = 0; i < j; i++)
                {
                    r[l] = fjac[i + ldfjac*j];
                    l += n - i - 1;
                }
                // 130
            }
            // 140
            r[l] = wa1[j];
            if (wa1[j] == 0.0) { sing = 1; }
        }
        // 150

        // Accumulate the orthogonal factor in fjac.
        qform(n, n, fjac, ldfjac, wa1);

        // Rescale if necessary
        if (mode != 2)
        {
            for (j = 0; j < n; j++)
            {
                diag[j] = fmax(diag[j], wa2[j]);
            }
            // 160
        }
        // 170

        // Beginning of the inner loop.
        while (1)
        {
            if (nprint > 0)
            {
                iflag = 0;
                if ((iter - 1) % nprint == 0) {
                    (*fcn)(&mut_n, x, fvec, fjac, &mut_ldf, &iflag);
                }
                if (iflag < 0) { goto EXIT300; }
            }
            // 190

            // Determine the direction of p.
            dogleg(n, r, diag, qtf, &delta, wa1, wa2, wa3);

            // Store the direction of p.
            for (j = 0; j < n; j++)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            // 200

            pnorm = enorm(n, wa3);

            // On the first iteration, adjust the initial step bound
            if (iter == 1) { delta = fmin(delta, pnorm); }

            // Evaluate the function at x + p and calculate its norm.
            iflag = 1;
            (*fcn)(&mut_n, wa2, wa4, fjac, &mut_ldf, &iflag);
            *nfev += 1;
            if (iflag < 0) { goto EXIT300; }
            fnorm1 = enorm(n, wa4);

            // Compute the scaled actual reduction
            actred = (fnorm1 < fnorm ? 1.0 - pow(fnorm1/fnorm, 2.0) : -1.0);

            // Compute the scaled predicted reduction.
            l = 0;
            for (i = 0; i < n; i++)
            {
                ssum = 0.0;
                for (j = i; j < n; j++)
                {
                    ssum += r[l] * wa1[j];
                    l++;
                }
                // 210
                wa3[i] = qtf[i] + ssum;
            }
            // 220

            temp = enorm(n, wa3);
            prered = (temp < fnorm ? 1.0 - pow(temp/fnorm, 2.0) : 0.0);

            // Compute the ratio of the actual to the predicted reduction.
            ratio = (prered > 0.0 ? actred/prered : 0.0);

            // Update the step bound
            if (ratio < 0.1)
            {
                ncsuc = 0;
                ncfail++;
                delta = 0.5*delta;
            } else {
                // 230
                ncfail = 0;
                ncsuc++;
                if ((ratio >= 0.5) || (ncsuc > 1)) { delta = fmax(delta, pnorm/0.5); }
                if (fabs(ratio - 1.0) <= 0.1) { delta = pnorm / 0.5; }
            }
            // 240

            // Test for successful iteration.
            if (ratio >= 0.0001)
            {
                // Successful iteration. Update x, fvec, and their norms.
                for (j = 0; j < n; j++)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                    fvec[j] = wa4[j];
                }
                // 250
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }
            // 260

            // Determine the progress of the iteration.
            nslow1++;
            if (actred >= 0.001) { nslow1 = 0; }
            if (jeval) { nslow2++; }
            if (actred >= 0.1) { nslow2 = 0; }

            // Test for convergence.
            if ((delta <= xtol*xnorm) || (fnorm == 0.0)) { *info = 1; }
            if (*info != 0) { goto EXIT300; }

            // Tests for termination and stringent tolerances
            if (*nfev >= maxfev) { *info = 2; }
            if (0.1*fmax(0.1*delta, pnorm) <= epsmch*xnorm) { *info = 3; }
            if (nslow2 == 5) { *info = 4; }
            if (nslow1 == 10) { *info = 5; }
            if (*info != 0) { goto EXIT300; }

            // Criterion for recalculating jacobian.
            if (ncfail == 2) { break; }

            // Calculate the rank one modification to the jacobian and update
            // qtf if necessary.
            for (j = 0; j < n; j++)
            {
                ssum = 0.0;
                for (i = 0; i < n; i++)
                {
                    ssum += fjac[i + ldfjac*j]*wa4[i];
                }
                // 270
                wa2[j] = (ssum- wa3[j])/pnorm;
                wa1[j] = diag[j]*((diag[j]*wa1[j])/pnorm);
                if (ratio >= 0.0001) { qtf[j] = ssum; }
            }
            // 280

            // Compute the qr factorization of the updated jacobian.
            r1updt(n, n, r, wa1, wa2, wa3, &sing);
            r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
            r1mpyq(1, n, qtf, 1, wa2, wa3);

            // End of the inner loop.
            jeval = 0;
            // 290
        }
        // End of the outer loop.
    }
    // 300

    // Termination, either normal or user imposed.
EXIT300:
    if (iflag < 0) { *info = iflag; }
    iflag = 0;
    if (nprint > 0) { (*fcn)(&mut_n, x, fvec, fjac, &mut_ldf, &iflag); }

    return;
}


void LMDIF(int(*fcn)(int* m, int* n, double* x, double* fvec, int* iflag),
           const int m, const int n, double* x, double* fvec, const double ftol,
           const double xtol, const double gtol, const int maxfev, const double epsfcn,
           double* diag, const int mode, const double factor, const int nprint,
           int* info, int* nfev, double* fjac, const int ldfjac, int* ipvt,
           double* qtf, double* wa1, double* wa2, double* wa3, double* wa4)
{
    //     **********
    //
    //     subroutine lmdif
    //
    //     the purpose of lmdif is to minimize the sum of the squares of
    //     m nonlinear functions in n variables by a modification of
    //     the levenberg-marquardt algorithm. the user must provide a
    //     subroutine which calculates the functions. the jacobian is
    //     then calculated by a forward-difference approximation.
    //
    //     the subroutine statement is
    //
    //       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
    //                        diag,mode,factor,nprint,info,nfev,fjac,
    //                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions. fcn must be declared
    //         in an external statement in the user calling
    //         program, and should be written as follows.
    //
    //         subroutine fcn(m,n,x,fvec,iflag)
    //         integer m,n,iflag
    //         double precision x(n),fvec(m)
    //         ----------
    //         calculate the functions at x and
    //         return this vector in fvec.
    //         ----------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of lmdif.
    //         in this case set iflag to a negative integer.
    //
    //       m is a positive integer input variable set to the number
    //         of functions.
    //
    //       n is a positive integer input variable set to the number
    //         of variables. n must not exceed m.
    //
    //       x is an array of length n. on input x must contain
    //         an initial estimate of the solution vector. on output x
    //         contains the final estimate of the solution vector.
    //
    //       fvec is an output array of length m which contains
    //         the functions evaluated at the output x.
    //
    //       ftol is a nonnegative input variable. termination
    //         occurs when both the actual and predicted relative
    //         reductions in the sum of squares are at most ftol.
    //         therefore, ftol measures the relative error desired
    //         in the sum of squares.
    //
    //       xtol is a nonnegative input variable. termination
    //         occurs when the relative error between two consecutive
    //         iterates is at most xtol. therefore, xtol measures the
    //         relative error desired in the approximate solution.
    //
    //       gtol is a nonnegative input variable. termination
    //         occurs when the cosine of the angle between fvec and
    //         any column of the jacobian is at most gtol in absolute
    //         value. therefore, gtol measures the orthogonality
    //         desired between the function vector and the columns
    //         of the jacobian.
    //
    //       maxfev is a positive integer input variable. termination
    //         occurs when the number of calls to fcn is at least
    //         maxfev by the end of an iteration.
    //
    //       epsfcn is an input variable used in determining a suitable
    //         step length for the forward-difference approximation. this
    //         approximation assumes that the relative errors in the
    //         functions are of the order of epsfcn. if epsfcn is less
    //         than the machine precision, it is assumed that the relative
    //         errors in the functions are of the order of the machine
    //         precision.
    //
    //       diag is an array of length n. if mode = 1 (see
    //         below), diag is internally set. if mode = 2, diag
    //         must contain positive entries that serve as
    //         multiplicative scale factors for the variables.
    //
    //       mode is an integer input variable. if mode = 1, the
    //         variables will be scaled internally. if mode = 2,
    //         the scaling is specified by the input diag. other
    //         values of mode are equivalent to mode = 1.
    //
    //       factor is a positive input variable used in determining the
    //         initial step bound. this bound is set to the product of
    //         factor and the euclidean norm of diag*x if nonzero, or else
    //         to factor itself. in most cases factor should lie in the
    //         interval (.1,100.). 100. is a generally recommended value.
    //
    //       nprint is an integer input variable that enables controlled
    //         printing of iterates if it is positive. in this case,
    //         fcn is called with iflag = 0 at the beginning of the first
    //         iteration and every nprint iterations thereafter and
    //         immediately prior to return, with x and fvec available
    //         for printing. if nprint is not positive, no special calls
    //         of fcn with iflag = 0 are made.
    //
    //       info is an integer output variable. if the user has
    //         terminated execution, info is set to the (negative)
    //         value of iflag. see description of fcn. otherwise,
    //         info is set as follows.
    //
    //         info = 0  improper input parameters.
    //
    //         info = 1  both actual and predicted relative reductions
    //                   in the sum of squares are at most ftol.
    //
    //         info = 2  relative error between two consecutive iterates
    //                   is at most xtol.
    //
    //         info = 3  conditions for info = 1 and info = 2 both hold.
    //
    //         info = 4  the cosine of the angle between fvec and any
    //                   column of the jacobian is at most gtol in
    //                   absolute value.
    //
    //         info = 5  number of calls to fcn has reached or
    //                   exceeded maxfev.
    //
    //         info = 6  ftol is too small. no further reduction in
    //                   the sum of squares is possible.
    //
    //         info = 7  xtol is too small. no further improvement in
    //                   the approximate solution x is possible.
    //
    //         info = 8  gtol is too small. fvec is orthogonal to the
    //                   columns of the jacobian to machine precision.
    //
    //       nfev is an integer output variable set to the number of
    //         calls to fcn.
    //
    //       fjac is an output m by n array. the upper n by n submatrix
    //         of fjac contains an upper triangular matrix r with
    //         diagonal elements of nonincreasing magnitude such that
    //
    //                t     t           t
    //               p *(jac *jac)*p = r *r,
    //
    //         where p is a permutation matrix and jac is the final
    //         calculated jacobian. column j of p is column ipvt(j)
    //         (see below) of the identity matrix. the lower trapezoidal
    //         part of fjac contains information generated during
    //         the computation of r.
    //
    //       ldfjac is a positive integer input variable not less than m
    //         which specifies the leading dimension of the array fjac.
    //
    //       ipvt is an integer output array of length n. ipvt
    //         defines a permutation matrix p such that jac*p = q*r,
    //         where jac is the final calculated jacobian, q is
    //         orthogonal (not stored), and r is upper triangular
    //         with diagonal elements of nonincreasing magnitude.
    //         column j of p is column ipvt(j) of the identity matrix.
    //
    //       qtf is an output array of length n which contains
    //         the first n elements of the vector (q transpose)*fvec.
    //
    //       wa1, wa2, and wa3 are work arrays of length n.
    //
    //       wa4 is a work array of length m.
    //
    //     subprograms called
    //
    //       user-supplied ...... fcn
    //
    //       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, iflag, iter, j, l, mut_m, mut_n;
    double actred, delta, dirder, fnorm, fnorm1, gnorm, par, pnorm;
    double prered, ratio, ssum, temp, temp1, temp2, xnorm;
    const double epsmch = dpmpar[0];

    // Avoid const poisoning
    mut_m = m;
    mut_n = n;
    // Initialize
    delta = 0.0;
    xnorm = 0.0;


    *info = 0;
    iflag = 0;
    *nfev = 0;

    // Check the input parameters for errors.
    if ((n <= 0) || (m < n) || (ldfjac < m) || (ftol < 0.0) || (xtol < 0.0) ||
        (gtol < 0.0) || (maxfev <= 0) || (factor <= 0.0)) { goto EXIT300; }
    if (mode == 2)
    {
        for (j = 0; j < n; j++)
        {
            if (diag[j] <= 0.0) { goto EXIT300; }
        }
        // 10
    }
    // 20

    // Evaluate the function at the starting point and calculate its norm.
    iflag = 1;
    (*fcn)(&mut_m, &mut_n, x, fvec, &iflag);
    *nfev = 1;
    if (iflag < 0) { goto EXIT300; }

    fnorm = enorm(m, fvec);

    // Initialize Levenberg-Marquardt parameter and iteration counter.
    par = 0.0;
    iter = 1;

    // Beginning of the outer loop

    while (1)
    {
        // Calculate the jacobian matrix.
        iflag = 2;
        fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, &iflag, epsfcn, wa4);
        *nfev += n;
        if (iflag < 0) { goto EXIT300; }

        // If requested call fcn to enable printing of iterates.
        if (nprint > 0)
        {
            iflag = 0;
            if ((iter - 1) % nprint == 0) {
                (*fcn)(&mut_m, &mut_n, x, fvec, &iflag);
            }
            if (iflag < 0) { goto EXIT300; }
        }
        // 40

        // Compute the qr factorization of the jacobian.
        qrfac(m, n, fjac, ldfjac, 1, ipvt, wa1, wa2, wa3);

        // On the first iteration and if mode is 1, scale according to the norms
        // of the columns if the initial jacobian.

        if (iter == 1)
        {
            if (mode != 2)
            {
                for (j = 0; j < n; j++)
                {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) { diag[j] = 1.0; }
                }
                // 50
            }
            // 60

            // On the first iteration, calculate the norm of the scaled x and
            // initialize the step bound delta.
            for (j = 0; j < n; j++)
            {
                wa3[j] = diag[j]*x[j];
            }
            // 70

            xnorm = enorm(n, wa3);
            delta = factor*xnorm;
            delta = (delta == 0.0 ? factor : delta);
        }
        // 80

        // Form (q.T)*fvec and store the first n components in qtf.
        for (i = 0; i < m; i++)
        {
            wa4[i] = fvec[i];
        }
        // 90

        for (j = 0; j < n; j++)
        {
            if (fjac[j + ldfjac*j] != 0.0)
            {
                ssum = 0.0;
                for (i = j; i < m; i++)
                {
                    ssum += fjac[i + ldfjac*j]*wa4[i];
                }
                // 100
                temp = -ssum / fjac[j + ldfjac*j];
                for (i = j; i < m; i++)
                {
                    wa4[i] += fjac[i + ldfjac*j]*temp;
                }
                // 110
            }
            // 120

            fjac[j + ldfjac*j] = wa1[j];
            qtf[j] = wa4[j];
        }
        // 130

        // Compute the norm of the scaled gradient.

        gnorm = 0.0;
        if (fnorm != 0.0)
        {
            for (j = 0; j < n; j++) {
                l = ipvt[j];
                if (wa2[l] != 0.0)
                {
                    ssum = 0.0;
                    for (i = 0; i <= j; i++)
                    {
                        ssum += fjac[i + ldfjac*j] * (qtf[i] / fnorm);
                    }
                    // 140
                    gnorm = fmax(gnorm, fabs(ssum / wa2[l]));
                }
                // 150
            }
            // 160
        }
        // 170

        // Test for convergence of the gradient norm.
        if (gnorm <= gtol) { *info = 4; }
        if (*info != 0) { goto EXIT300; }

        // Rescale if necessary
        if (mode != 2)
        {
            for (j = 0; j < n; j++)
            {
                diag[j] = fmax(diag[j], wa2[j]);
            }
            // 180
        }
        // 190

        // Beginning of the inner loop
        do
        {
            // Determine the Levenberg-Marquardt parameter
            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4);

            // Store the direction p and x+p. Calculate the norm of p.
            for (j = 0; j < n; j++)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            // 210
            pnorm = enorm(n, wa3);

            // On the first iteration, adjust the initial step bound.
            if (iter == 1) { delta = fmin(delta, pnorm); }

            // Evaluate the function at x+p and calculate its norm.
            iflag = 1;
            (*fcn)(&mut_m, &mut_n, wa2, wa4, &iflag);
            *nfev += 1;
            if (iflag < 0){ goto EXIT300; }
            fnorm1 = enorm(m, wa4);

            // Compute the scaled actual reduction.
            actred = (0.1*fnorm1 < fnorm ? 1.0 - pow(fnorm1/fnorm, 2.0) : -1.0);

            // Compute the scaled predicted reduction and the scaled directional
            // derivative.
            for (j = 0; j < n; j++)
            {
                wa3[j] = 0.0;
                l = ipvt[j];
                temp = wa1[l];
                for (i = 0; i <= j; i++)
                {
                    wa3[i] += fjac[i + ldfjac*j]*temp;
                }
                // 220
            }
            // 230

            temp1 = enorm(n, wa3)/fnorm;
            temp2 = (sqrt(par)*pnorm)/fnorm;
            prered = temp1*temp1 + temp2*temp2/0.5;
            dirder = -(temp1*temp1 + temp2*temp2);

            // Compute the ratio of the actual to the predicted reduction.
            ratio = ( prered != 0.0 ? actred/prered : 0.0);

            // Update the step bound
            if (ratio <= 0.25)
            {
                temp = (actred < 0.0 ? (0.5*dirder)/(dirder + 0.5*actred): 0.5);
                if ((0.1*fnorm1 >= fnorm) || (temp < 0.1)) { temp = 0.1; }
                delta = temp*fmin(delta, pnorm/0.1);
                par /= temp;
            } else {
                // 240
                if ((par == 0.0) || (ratio >= 0.75))
                {
                    delta = pnorm/0.5;
                    par *= 0.5;
                }
                // 250
            }
            // 260

            // Test for succesful iteration.
            if (ratio >= 0.0001)
            {
                // Successful iteration. Update x, fvec, and their norms.
                for (j = 0; j < n; ++j)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                // 270
                for (i = 0; i < m; i++)
                {
                    fvec[i] = wa4[i];
                }
                // 280
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }
            // 290

            // Tests for convergence.
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (0.5*ratio <= 1.0)) { *info = 1; }
            if (delta <= xtol*xnorm) { *info = 2; }
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (0.5*ratio <= 1.0) && (*info == 2)) { *info = 3; }
            if (*info != 0){ goto EXIT300; }

            // Tests for termination and stringent tolerances.
            if (*nfev >= maxfev) { *info = 5; }
            if ((fabs(actred) <= epsmch) && (prered <= epsmch) && (0.5*ratio <= 1.0)) { *info = 6; }
            if (delta <= epsmch*xnorm) { *info = 7; }
            if (gnorm <= epsmch) { *info = 8; }
            if (*info != 0) { goto EXIT300; }

            // End of the inner loop. Repeat if iteration unsuccessful.
        } while (ratio < 0.0001);

        // End of the outer loop.
    }
    // 300

EXIT300:
    if (iflag < 0) { *info = iflag; }
    iflag = 0;
    if (nprint > 0) { (*fcn)(&mut_m, &mut_n, x, fvec, &iflag); }

    return;
}


void LMDER(int(*fcn)(int* m, int* n,double* x, double* fvec, double* fjac,
                     int* ldfjac, int* iflag),
           const int m, const int n, double* x, double* fvec, double* fjac,
           const int ldfjac, const double ftol, const double xtol, const double gtol,
           const int maxfev, double* diag, const int mode, const double factor,
           const int nprint, int* info, int* nfev, int* njev, int* ipvt,
           double* qtf, double* wa1, double* wa2, double* wa3, double* wa4)
{
    //     **********
    //
    //     subroutine lmder
    //
    //     the purpose of lmder is to minimize the sum of the squares of
    //     m nonlinear functions in n variables by a modification of
    //     the levenberg-marquardt algorithm. the user must provide a
    //     subroutine which calculates the functions and the jacobian.
    //
    //     the subroutine statement is
    //
    //       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
    //                        maxfev,diag,mode,factor,nprint,info,nfev,
    //                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions and the jacobian. fcn must
    //         be declared in an external statement in the user
    //         calling program, and should be written as follows.
    //
    //         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    //         integer m,n,ldfjac,iflag
    //         double precision x(n),fvec(m),fjac(ldfjac,n)
    //         ----------
    //         if iflag = 1 calculate the functions at x and
    //         return this vector in fvec. do not alter fjac.
    //         if iflag = 2 calculate the jacobian at x and
    //         return this matrix in fjac. do not alter fvec.
    //         ----------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of lmder.
    //         in this case set iflag to a negative integer.
    //
    //       m is a positive integer input variable set to the number
    //         of functions.
    //
    //       n is a positive integer input variable set to the number
    //         of variables. n must not exceed m.
    //
    //       x is an array of length n. on input x must contain
    //         an initial estimate of the solution vector. on output x
    //         contains the final estimate of the solution vector.
    //
    //       fvec is an output array of length m which contains
    //         the functions evaluated at the output x.
    //
    //       fjac is an output m by n array. the upper n by n submatrix
    //         of fjac contains an upper triangular matrix r with
    //         diagonal elements of nonincreasing magnitude such that
    //
    //                t     t           t
    //               p *(jac *jac)*p = r *r,
    //
    //         where p is a permutation matrix and jac is the final
    //         calculated jacobian. column j of p is column ipvt(j)
    //         (see below) of the identity matrix. the lower trapezoidal
    //         part of fjac contains information generated during
    //         the computation of r.
    //
    //       ldfjac is a positive integer input variable not less than m
    //         which specifies the leading dimension of the array fjac.
    //
    //       ftol is a nonnegative input variable. termination
    //         occurs when both the actual and predicted relative
    //         reductions in the sum of squares are at most ftol.
    //         therefore, ftol measures the relative error desired
    //         in the sum of squares.
    //
    //       xtol is a nonnegative input variable. termination
    //         occurs when the relative error between two consecutive
    //         iterates is at most xtol. therefore, xtol measures the
    //         relative error desired in the approximate solution.
    //
    //       gtol is a nonnegative input variable. termination
    //         occurs when the cosine of the angle between fvec and
    //         any column of the jacobian is at most gtol in absolute
    //         value. therefore, gtol measures the orthogonality
    //         desired between the function vector and the columns
    //         of the jacobian.
    //
    //       maxfev is a positive integer input variable. termination
    //         occurs when the number of calls to fcn with iflag = 1
    //         has reached maxfev.
    //
    //       diag is an array of length n. if mode = 1 (see
    //         below), diag is internally set. if mode = 2, diag
    //         must contain positive entries that serve as
    //         multiplicative scale factors for the variables.
    //
    //       mode is an integer input variable. if mode = 1, the
    //         variables will be scaled internally. if mode = 2,
    //         the scaling is specified by the input diag. other
    //         values of mode are equivalent to mode = 1.
    //
    //       factor is a positive input variable used in determining the
    //         initial step bound. this bound is set to the product of
    //         factor and the euclidean norm of diag*x if nonzero, or else
    //         to factor itself. in most cases factor should lie in the
    //         interval (.1,100.).100. is a generally recommended value.
    //
    //       nprint is an integer input variable that enables controlled
    //         printing of iterates if it is positive. in this case,
    //         fcn is called with iflag = 0 at the beginning of the first
    //         iteration and every nprint iterations thereafter and
    //         immediately prior to return, with x, fvec, and fjac
    //         available for printing. fvec and fjac should not be
    //         altered. if nprint is not positive, no special calls
    //         of fcn with iflag = 0 are made.
    //
    //       info is an integer output variable. if the user has
    //         terminated execution, info is set to the (negative)
    //         value of iflag. see description of fcn. otherwise,
    //         info is set as follows.
    //
    //         info = 0  improper input parameters.
    //
    //         info = 1  both actual and predicted relative reductions
    //                   in the sum of squares are at most ftol.
    //
    //         info = 2  relative error between two consecutive iterates
    //                   is at most xtol.
    //
    //         info = 3  conditions for info = 1 and info = 2 both hold.
    //
    //         info = 4  the cosine of the angle between fvec and any
    //                   column of the jacobian is at most gtol in
    //                   absolute value.
    //
    //         info = 5  number of calls to fcn with iflag = 1 has
    //                   reached maxfev.
    //
    //         info = 6  ftol is too small. no further reduction in
    //                   the sum of squares is possible.
    //
    //         info = 7  xtol is too small. no further improvement in
    //                   the approximate solution x is possible.
    //
    //         info = 8  gtol is too small. fvec is orthogonal to the
    //                   columns of the jacobian to machine precision.
    //
    //       nfev is an integer output variable set to the number of
    //         calls to fcn with iflag = 1.
    //
    //       njev is an integer output variable set to the number of
    //         calls to fcn with iflag = 2.
    //
    //       ipvt is an integer output array of length n. ipvt
    //         defines a permutation matrix p such that jac*p = q*r,
    //         where jac is the final calculated jacobian, q is
    //         orthogonal (not stored), and r is upper triangular
    //         with diagonal elements of nonincreasing magnitude.
    //         column j of p is column ipvt(j) of the identity matrix.
    //
    //       qtf is an output array of length n which contains
    //         the first n elements of the vector (q transpose)*fvec.
    //
    //       wa1, wa2, and wa3 are work arrays of length n.
    //
    //       wa4 is a work array of length m.
    //
    //     subprograms called
    //
    //       user-supplied ...... fcn
    //
    //       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, iflag, iter, j, l, mut_ldf, mut_m, mut_n;
    double actred, delta, dirder, fnorm, fnorm1, gnorm, par, pnorm, prered;
    double ratio, ssum, temp, temp1, temp2, xnorm;
    const double epsmch = dpmpar[0];

    // Avoid const poisoning
    mut_ldf = ldfjac;
    mut_m = m;
    mut_n = n;
    // Initialize
    delta = 0.0;
    xnorm = 0.0;

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

    if ((n <= 0) || (m < n) || (ldfjac < m) || (ftol < 0.0) || (xtol < 0.0) ||
	    (gtol < 0.0) || (maxfev <= 0) || (factor <= 0.0)) { goto EXIT300; }
    if (mode == 2) {
        for (j = 0; j < n; j++)
        {
            if (diag[j] <= 0.) { goto EXIT300; }
        }
        // 10
    }
    // 20
    fnorm = enorm(m, fvec);

    iflag = 1;
    (*fcn)(&mut_m, &mut_n, x, fvec, fjac, &mut_ldf, &iflag);
    *nfev = 1;
    if (iflag < 0) { goto EXIT300; }
    fnorm = enorm(m, fvec);

    // Initialize Levenberg-Marquardt parameter and iteration counter.
    par = 0.0;
    iter = 1;

    // Beginning the outer loop
    while (1)
    {
        // Calculate the jacobian matrix.
        iflag = 2;
        (*fcn)(&mut_m, &mut_n, x, fvec, fjac, &mut_ldf, &iflag);
        *njev += 1;
        if (iflag < 0) { goto EXIT300; }

        // If requested, call fcn to enable printing of the iterates
        if (nprint > 0)
        {
            iflag = 0;
            if ((iter - 1) % nprint == 0) {
                (*fcn)(&mut_m, &mut_n, x, fvec, fjac, &mut_ldf, &iflag);
            }
            if (iflag < 0) { goto EXIT300; }
        }
        // 40

        // Compute the qr factorization of the jacobian.
        qrfac(m, n, fjac, ldfjac, 1, ipvt, wa1, wa2, wa3);

        // On the first iteration and if mode is 1, scale according to the norms
        // of the columns of the initial jacobian.
        if (iter == 1)
        {
            if (mode != 2)
            {
                for (j = 0; j < n; j++)
                {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) { diag[j] = 1.0; }
                }
                // 50
            }
            // 60

            // On the first iteration, calculate the norm of the scaled x and
            // initialize the step bound delta.
            for (j = 0; j < n; j++)
            {
                wa3[j] = diag[j]*x[j];
            }
            // 70
            xnorm = enorm(n, wa3);
            delta = factor*xnorm;
            if (delta == 0.0) { delta = factor; }
        }
        // 80

        // Form (q.T)*fvec and store the first n components in qtf.
        for (i = 0; i < m; i++)
        {
            wa4[i] = fvec[i];
        }
        // 90
        for (j = 0; j < n; j++)
        {
            if (fjac[j + ldfjac*j] != 0.0)
            {
                ssum = 0.0;
                for (i = j; i < m; i++)
                {
                    ssum += fjac[i + ldfjac*j]*wa4[i];
                }
                // 100

                temp = -ssum/fjac[j + ldfjac*j];
                for (i = j; i < m; i++)
                {
                    wa4[i] += fjac[i + ldfjac*j]*temp;
                }
                // 110
            }
            // 120
            fjac[j + ldfjac*j] = wa1[j];
            qtf[j] = wa4[j];
        }
        // 130

        // Compute the norm of the scaled gradient
        gnorm = 0.0;
        if (fnorm != 0.0)
        {
            for (j = 0; j < n; j++)
            {
                l = ipvt[j];
                if (wa2[l] != 0.0)
                {
                    ssum = 0.0;
                    for (i = 0; i <= j; i++)
                    {
                        ssum += fjac[i + ldfjac*j]*(qtf[i]/fnorm);
                    }
                    // 140
                    gnorm = fmax(gnorm, fabs(ssum/wa2[l]));
                }
                // 150
            }
            // 160
        }
        // 170

        // Test for convergence of the gradient norm.
        if (gnorm <= gtol) { *info = 4; }
        if (*info != 0) { goto EXIT300; }

        // Rescale if necessary
        if (mode != 2)
        {
            for (j = 0; j < n; j++)
            {
                diag[j] = fmax(diag[j], wa2[j]);
            }
            // 180
        }
        // 190

        // Beginning of the inner loop
        do
        {
            // Determine the Levenberg-Marquardt parameter.
            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4);

            // Store the direction p and x+p. Calculate the norm of p.
            for (j = 0; j < n; j++)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            // 210
            pnorm = enorm(n, wa3);

            // On the first iteration, adjust the initial step bound.
            if (iter == 1) { delta = fmin(delta, pnorm); }

            // Evaluate the function at x+p and calculate its norm.
            iflag = 1;
            (*fcn)(&mut_m, &mut_n, wa2, wa4, fjac, &mut_ldf, &iflag);
            *nfev += 1;
            if (iflag < 0) { goto EXIT300; }
            fnorm1 = enorm(m, wa4);

            // Compute the scaled actual reduction.
            actred = (0.1*fnorm1 < fnorm ? 1.0 - pow(fnorm1/fnorm, 2.0) : -1.0);

            // Compute the scaled predicted reduction and the scaled directional
            // derivative.
            for (j = 0; j < n; j++)
            {
                wa3[j] = 0.0;
                l = ipvt[j];
                temp = wa1[l];
                for (i = 0; i <= j; i++)
                {
                    wa3[i] += fjac[i + ldfjac*j]*temp;
                }
                // 220
            }
            // 230

            temp1 = enorm(n, wa3)/fnorm;
            temp2 = (sqrt(par)*pnorm)/fnorm;
            prered = temp1*temp1 + temp2*temp2/0.5;
            dirder = -(temp1*temp1 + temp2*temp2);

            // Compute the ratio of the actual to the predicted reduction.
            ratio = ( prered != 0.0 ? actred/prered : 0.0);

            // Update the step bound
            if (ratio <= 0.25)
            {
                temp = (actred < 0.0 ? (0.5*dirder)/(dirder + 0.5*actred): 0.5);
                if ((0.1*fnorm1 >= fnorm) || (temp < 0.1)) { temp = 0.1; }
                delta = temp*fmin(delta, pnorm/0.1);
                par /= temp;
            } else {
                // 240
                if ((par == 0.0) || (ratio >= 0.75))
                {
                    delta = pnorm/0.5;
                    par *= 0.5;
                }
                // 250
            }
            // 260

            // Test for succesful iteration.
            if (ratio >= 0.0001)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                // 270
                for (i = 0; i < m; i++)
                {
                    fvec[i] = wa4[i];
                }
                // 280
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }
            // 290

            // Tests for convergence.
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (0.5*ratio <= 1.0)) { *info = 1; }
            if (delta <= xtol*xnorm) { *info = 2; }
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (0.5*ratio <= 1.0) && (*info == 2)) { *info = 3; }
            if (*info != 0) { goto EXIT300; }

            // Tests for termination and stringent tolerances.
            if (*nfev >= maxfev) { *info = 5; }
            if ((fabs(actred) <= epsmch) && (prered <= epsmch) && (0.5*ratio <= 1.0)) { *info = 6; }
            if (delta <= epsmch*xnorm) { *info = 7; }
            if (gnorm <= epsmch) { *info = 8; }
            if (*info != 0) { goto EXIT300; }

            // End of the inner loop. Repeat if iteration unsuccessful.
        } while (ratio < 0.0001);
    }
    // 300

EXIT300:
    if (iflag < 0) { *info = iflag; }
    iflag = 0;
    if (nprint > 0) { (*fcn)(&mut_m, &mut_n, x, fvec, fjac, &mut_ldf, &iflag); }

    return;
}


void LMSTR(int(*fcn)(int* m, int* n, double* x, double* fvec, double* wa3, int* iflag),
           const int m, const int n, double* x, double* fvec, double* fjac,
           const int ldfjac, const double ftol, const double xtol, const double gtol,
           const int maxfev, double* diag, const int mode, const double factor,
           const int nprint, int* info, int* nfev, int* njev, int* ipvt, double* qtf,
           double* wa1, double* wa2, double* wa3, double* wa4)
{
    //     **********
    //
    //     subroutine lmstr
    //
    //     the purpose of lmstr is to minimize the sum of the squares of
    //     m nonlinear functions in n variables by a modification of
    //     the levenberg-marquardt algorithm which uses minimal storage.
    //     the user must provide a subroutine which calculates the
    //     functions and the rows of the jacobian.
    //
    //     the subroutine statement is
    //
    //       subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
    //                        maxfev,diag,mode,factor,nprint,info,nfev,
    //                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions and the rows of the jacobian.
    //         fcn must be declared in an external statement in the
    //         user calling program, and should be written as follows.
    //
    //         subroutine fcn(m,n,x,fvec,fjrow,iflag)
    //         integer m,n,iflag
    //         double precision x(n),fvec(m),fjrow(n)
    //         ----------
    //         if iflag = 1 calculate the functions at x and
    //         return this vector in fvec.
    //         if iflag = i calculate the (i-1)-st row of the
    //         jacobian at x and return this vector in fjrow.
    //         ----------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of lmstr.
    //         in this case set iflag to a negative integer.
    //
    //       m is a positive integer input variable set to the number
    //         of functions.
    //
    //       n is a positive integer input variable set to the number
    //         of variables. n must not exceed m.
    //
    //       x is an array of length n. on input x must contain
    //         an initial estimate of the solution vector. on output x
    //         contains the final estimate of the solution vector.
    //
    //       fvec is an output array of length m which contains
    //         the functions evaluated at the output x.
    //
    //       fjac is an output n by n array. the upper triangle of fjac
    //         contains an upper triangular matrix r such that
    //
    //                t     t           t
    //               p *(jac *jac)*p = r *r,
    //
    //         where p is a permutation matrix and jac is the final
    //         calculated jacobian. column j of p is column ipvt(j)
    //         (see below) of the identity matrix. the lower triangular
    //         part of fjac contains information generated during
    //         the computation of r.
    //
    //       ldfjac is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array fjac.
    //
    //       ftol is a nonnegative input variable. termination
    //         occurs when both the actual and predicted relative
    //         reductions in the sum of squares are at most ftol.
    //         therefore, ftol measures the relative error desired
    //         in the sum of squares.
    //
    //       xtol is a nonnegative input variable. termination
    //         occurs when the relative error between two consecutive
    //         iterates is at most xtol. therefore, xtol measures the
    //         relative error desired in the approximate solution.
    //
    //       gtol is a nonnegative input variable. termination
    //         occurs when the cosine of the angle between fvec and
    //         any column of the jacobian is at most gtol in absolute
    //         value. therefore, gtol measures the orthogonality
    //         desired between the function vector and the columns
    //         of the jacobian.
    //
    //       maxfev is a positive integer input variable. termination
    //         occurs when the number of calls to fcn with iflag = 1
    //         has reached maxfev.
    //
    //       diag is an array of length n. if mode = 1 (see
    //         below), diag is internally set. if mode = 2, diag
    //         must contain positive entries that serve as
    //         multiplicative scale factors for the variables.
    //
    //       mode is an integer input variable. if mode = 1, the
    //         variables will be scaled internally. if mode = 2,
    //         the scaling is specified by the input diag. other
    //         values of mode are equivalent to mode = 1.
    //
    //       factor is a positive input variable used in determining the
    //         initial step bound. this bound is set to the product of
    //         factor and the euclidean norm of diag*x if nonzero, or else
    //         to factor itself. in most cases factor should lie in the
    //         interval (.1,100.). 100. is a generally recommended value.
    //
    //       nprint is an integer input variable that enables controlled
    //         printing of iterates if it is positive. in this case,
    //         fcn is called with iflag = 0 at the beginning of the first
    //         iteration and every nprint iterations thereafter and
    //         immediately prior to return, with x and fvec available
    //         for printing. if nprint is not positive, no special calls
    //         of fcn with iflag = 0 are made.
    //
    //       info is an integer output variable. if the user has
    //         terminated execution, info is set to the (negative)
    //         value of iflag. see description of fcn. otherwise,
    //         info is set as follows.
    //
    //         info = 0  improper input parameters.
    //
    //         info = 1  both actual and predicted relative reductions
    //                   in the sum of squares are at most ftol.
    //
    //         info = 2  relative error between two consecutive iterates
    //                   is at most xtol.
    //
    //         info = 3  conditions for info = 1 and info = 2 both hold.
    //
    //         info = 4  the cosine of the angle between fvec and any
    //                   column of the jacobian is at most gtol in
    //                   absolute value.
    //
    //         info = 5  number of calls to fcn with iflag = 1 has
    //                   reached maxfev.
    //
    //         info = 6  ftol is too small. no further reduction in
    //                   the sum of squares is possible.
    //
    //         info = 7  xtol is too small. no further improvement in
    //                   the approximate solution x is possible.
    //
    //         info = 8  gtol is too small. fvec is orthogonal to the
    //                   columns of the jacobian to machine precision.
    //
    //       nfev is an integer output variable set to the number of
    //         calls to fcn with iflag = 1.
    //
    //       njev is an integer output variable set to the number of
    //         calls to fcn with iflag = 2.
    //
    //       ipvt is an integer output array of length n. ipvt
    //         defines a permutation matrix p such that jac*p = q*r,
    //         where jac is the final calculated jacobian, q is
    //         orthogonal (not stored), and r is upper triangular.
    //         column j of p is column ipvt(j) of the identity matrix.
    //
    //       qtf is an output array of length n which contains
    //         the first n elements of the vector (q transpose)*fvec.
    //
    //       wa1, wa2, and wa3 are work arrays of length n.
    //
    //       wa4 is a work array of length m.
    //
    //     subprograms called
    //
    //       user-supplied ...... fcn
    //
    //       minpack-supplied ... dpmpar,enorm,lmpar,qrfac,rwupdt
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
    //     jorge j. more
    //
    //     **********
    int i, iflag, iter, j, l, mut_ldf, mut_m, mut_n, sing;
    double actred, delta, dirder, fnorm, fnorm1, gnorm, par, pnorm, prered, ratio;
    double ssum, temp, temp1, temp2, xnorm;

    const double epsmch = dpmpar[0];

    // Avoid const poisoning
    mut_ldf = ldfjac;
    mut_m = m;
    mut_n = n;

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

    // Check the input parameters for errors.

    if ((n <= 0) || (m < n) || (ldfjac < n) || (ftol < 0.0) || (xtol < 0.0) ||
	    (gtol < 0.0) || (maxfev <= 0) || (factor <= 0.0)) { goto EXIT340; }
    if (mode == 2) {
        for (j = 0; j < n; ++j) {
            if (diag[j] <= 0.) { goto EXIT340; }
        }
        // 10
    }
    // 20

    // Evaluate the function at the starting point and calculate its norm.
    iflag = 1;
    (*fcn)(&mut_m, &mut_n, x, fvec, wa3, &iflag);
    *nfev += 1;
    if (iflag < 0) { goto EXIT340; }
    fnorm = enorm(m, fvec);

    // Initialize Levenberg-Marquardt parameter and iteration counter.
    par = 0.0;
    iter = 1;

    // Beginning of the outer loop.
    while (1)
    {
        // If requested, call fcn to enable printing of the iterates
        if (nprint > 0)
        {
            iflag = 0;
            if ((iter - 1) % nprint == 0) { (*fcn)(&mut_m, &mut_n, x, fvec, wa3, &iflag); }
            if (iflag < 0) { goto EXIT340; }
        }
        // 40

        // Compute the qr factorization of the jacobian matrix calculated one
        // row at a time, while simultaneously forming (q.T)*fvec and storing
        // the first n components in qtf.
        for (j = 0; j < n; j++)
        {
            qtf[j] = 0.0;
            for (i = 0; i < n; i++)
            {
                fjac[i + ldfjac*j] = 0.0;
            }
            // 50
        }
        // 60
        iflag = 2;
        for (i = 0; i < m; i++)
        {
            (*fcn)(&mut_m, &mut_n, x, fvec, wa3, &iflag);
            if (iflag < 0) { goto EXIT340; }
            temp = fvec[i];
            rwupdt(n, fjac, ldfjac, wa3, qtf, &temp, wa1, wa2);
            iflag++;
        }
        // 70
        *njev += 1;

        // If the jacobian is rank deficient, call qrfac to reorder its columns
        // and update the components of qtf.
        sing = 0;
        for (j = 0; j < n; j++)
        {
            if (fjac[j + ldfjac*j] != 0.0) { sing = 1; }
            ipvt[j] = j;
            wa2[j] = enorm(j, &fjac[ldfjac*j]);
        }
        // 80
        if (sing)
        {
            qrfac(n, n, fjac, ldfjac, 1, ipvt, wa1, wa2, wa3);
            for (j = 0; j < n; j++)
            {
                if (fjac[j + ldfjac*j] != 0.0)
                {
                    ssum = 0.0;
                    for (i = j; i < n; i++)
                    {
                        ssum += fjac[i + ldfjac*j] * qtf[i];
                    }
                    // 90
                    temp = -ssum / fjac[j + ldfjac*j];
                    for (i = j; i < n; i++)
                    {
                        qtf[i] += fjac[i + ldfjac*j] * temp;
                    }
                    // 100
                }
                // 110
                fjac[j + ldfjac*j] = wa1[j];
            }
            // 120
        }
        // 130

        // On the first iteration and if mode is 1, scale according to the norms
        // if the columns of the initial jacobian.
        if (iter == 1)
        {
            if (mode != 2)
            {
                for (j = 0; j < n; j++)
                {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) { diag[j] = 1.0; }
                }
                // 140
            }
            // 150

            // On the first iteration, calculate the norm of the scaled x and
            // initialize the step bound delta.
            for (j = 0; j <n; j++)
            {
                wa3[j] = diag[j]*x[j];
            }
            // 160

            xnorm = enorm(n, wa3);
            delta = factor*xnorm;
            if (delta == 0.0) { delta = factor; }
        }
        // 170

        // Compute the norm of the scaled gradient.
        gnorm = 0.0;
        if (fnorm != 0.0)
        {
            for (j = 0; j < n; j++)
            {
                l = ipvt[j];
                if (wa2[l] != 0.0)
                {
                    ssum = 0.0;
                    for (i = 0; i <= j; i++)
                    {
                        ssum += fjac[i + ldfjac*j]*(qtf[i]/fnorm);
                    }
                    // 180
                    gnorm = fmax(gnorm, fabs(ssum/wa2[l]));
                }
                // 190
            }
        }
        // 210

        // Test for convergence of the gradient norm.
        if (gnorm <= gtol) { *info = 4; }
        if (*info != 0) { goto EXIT340; }

        // Rescale if necessary.
        if (mode != 2)
        {
            for (j = 0; j < n; j++)
            {
                diag[j] = fmax(diag[j], wa2[j]);
            }
            // 220
        }
        // 230

        // Beginning of the inner loop.
        do
        {
            // Determine the Levenberg-Marquardt parameter.
            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4);

            // Store the direction p and x+p. Calculate the norm of p.

            for (j = 0; j < n; j++)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            // 250

            pnorm = enorm(n, wa3);

            // On the first iteration, adjust the initial step bound.
            if (iter == 1) { delta = fmin(delta, pnorm); }

            // Evaluate the function at x+p and calculate its norm.
            iflag = 1;
            (*fcn)(&mut_m, &mut_n, x, fvec, wa3, &iflag);
            *nfev += 1;
            if (iflag < 0) { goto EXIT340; }
            fnorm1= enorm(m, wa4);

            // Compute the scaled actual reduction
            actred = (0.1*fnorm1 < fnorm ? 1.0 - pow(fnorm1/fnorm, 2.0) : -1.0);

            // Compute the scaled predicted reduction and the scaled directional
            // derivative.
            for (j = 0; j < n; j++)
            {
                wa3[j] = 0.0;
                l = ipvt[j];
                temp = wa1[l];
                for (i = j; i <= j; i++)
                {
                    wa3[i] += fjac[i + ldfjac*j]*temp;
                }
                // 260
            }
            // 270

            temp1 = enorm(n, wa3)/fnorm;
            temp2 = (sqrt(par)*pnorm)/fnorm;
            prered = temp1*temp1 + temp2*temp2/0.5;
            dirder = -(temp1*temp1 + temp2*temp2);

            // Compute the ratio of the actual to the predicted reduction.
            ratio = (prered != 0.0 ? actred/prered : 0.0);

            // Update the step bound.
            if (ratio <= 0.25)
            {
                temp = (actred < 0.0 ? (0.5*dirder)/(dirder + 0.5*actred) : 0.5);
                if ((0.1*fnorm1 >= fnorm) || (temp < 0.1)) { temp = 0.1; }
                par /= temp;
            } else {
                // 280
                if ((par == 0.0) || (ratio >= 0.75))
                {
                    delta = pnorm/0.5;
                    par *= 0.5;
                }
                // 290
            }
            // 300

            // Test for succesful iteration.
            if (ratio >= 0.0001)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                // 310
                for (i = 0; i < m; i++)
                {
                    fvec[i] = wa4[i];
                }
                // 320
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }
            // 330

            // Tests for convergence.
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (0.5*ratio <= 1.0)) { *info = 1; }
            if (delta <= xtol*xnorm) { *info = 2; }
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (0.5*ratio <= 1.0) && (*info == 2)) { *info = 3; }
            if (*info != 0) { goto EXIT340; }

            // Tests for termination and stringent tolerances.
            if (*nfev >= maxfev) { *info = 5; }
            if ((fabs(actred) <= epsmch) && (prered <= epsmch) && (0.5*ratio <= 1.0)) { *info = 6; }
            if (delta <= epsmch*xnorm) { *info = 7; }
            if (gnorm <= epsmch) { *info = 8; }
            if (info != 0) { goto EXIT340; }

            // End of the inner loop. Repeat if iteration unsuccessful.
        } while (ratio < 0.0001);
        // End of the outer loop.
    }
    // 340

EXIT340:
    if (iflag < 0) { *info = iflag; }
    iflag = 0;
    if (nprint > 0) { (*fcn)(&mut_m, &mut_n, x, fvec, wa3, &iflag); }

    return;
}


void dogleg(const int n, const double* r, const double* diag, const double* qtb,
            const double* delta, double* x, double* wa1, double* wa2)
{
    //     **********
    //
    //     subroutine dogleg
    //
    //     given an m by n matrix a, an n by n nonsingular diagonal
    //     matrix d, an m-vector b, and a positive number delta, the
    //     problem is to determine the convex combination x of the
    //     gauss-newton and scaled gradient directions that minimizes
    //     (a*x - b) in the least squares sense, subject to the
    //     restriction that the euclidean norm of d*x be at most delta.
    //
    //     this subroutine completes the solution of the problem
    //     if it is provided with the necessary information from the
    //     qr factorization of a. that is, if a = q*r, where q has
    //     orthogonal columns and r is an upper triangular matrix,
    //     then dogleg expects the full upper triangle of r and
    //     the first n components of (q transpose)*b.
    //
    //     the subroutine statement is
    //
    //       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
    //
    //     where
    //
    //       n is a positive integer input variable set to the order of r.
    //
    //       r is an input array of length lr which must contain the upper
    //         triangular matrix r stored by rows.
    //
    //       lr is a positive integer input variable not less than
    //         (n*(n+1))/2.
    //
    //       diag is an input array of length n which must contain the
    //         diagonal elements of the matrix d.
    //
    //       qtb is an input array of length n which must contain the first
    //         n elements of the vector (q transpose)*b.
    //
    //       delta is a positive input variable which specifies an upper
    //         bound on the euclidean norm of d*x.
    //
    //       x is an output array of length n which contains the desired
    //         convex combination of the gauss-newton direction and the
    //         scaled gradient direction.
    //
    //       wa1 and wa2 are work arrays of length n.
    //
    //     subprograms called
    //
    //       minpack-supplied ... dpmpar,enorm
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, jj, k, l;
    double alpha, bnorm, epsmch, gnorm, qnorm, sgnorm, ssum, temp;

    epsmch = dpmpar[0];

    // First, calculate the Gauss-Newton direction
    jj = (n * (n + 1)) / 2;

    for (k = 0; k < n; k++)
    {
        j = n - k - 1;
        jj -= k + 1;
        l = jj + 1;
        ssum = 0.0;
        if (n - 1 > j)
        {
            for (i = j + 1; i < n; i++)
            {
                ssum += r[l]*x[i];
                l++;
            }
            // 10
        }
        // 20
        temp = r[jj];
        if (temp == 0.0)
        {
            l = j;
            for (i = 0; i <= j; i++)
            {
                temp = fmax(temp, fabs(r[l]));
                l += n - i - 1;
            }
            // 30
            temp = epsmch*temp;
            if (temp == 0.0) { temp = epsmch; }
        }
        // 40
        x[j] = (qtb[j] - ssum)/temp;
    }
    // 50

    // Test whether the gauss-newton direction is acceptable.
    for (j = 0; j < n; j++)
    {
        wa1[j] = 0.0;
        wa2[j] = diag[j] * x[j];
    }
    // 60
    qnorm = enorm(n, wa2);
    if (qnorm <= *delta) { return; }

    // The Gauss-Newton direction is not acceptable.
    // Next, calculate the scaled gradient direction.
    l = 0;
    for (j = 0; j < n; j++) {
        temp = qtb[j];
        for (i = j; i < n; i++) {
            wa1[i] += r[l] * temp;
            l++;
        }
        // 70
        wa1[j] /= diag[j];
    }
    // 80

    gnorm = enorm(n, wa1);
    sgnorm = 0.0;
    alpha = *delta / qnorm;
    if (gnorm != 0.0)
    {
        // Calculate the point along the scaled gradient
        // at which the quadratic is minimized.
        for (j = 0; j < n; j++)
        {
            wa1[j] /= gnorm;
            wa1[j] /= diag[j];
        }
        // 90

        l = 0;
        for (j = 0; j < n; j++) {
            ssum = 0.0;
            for (i = j; i < n; i++) {
                ssum += r[l] * wa1[i];
                l++;
            }
            // 100
            wa2[j] = ssum;
        }
        // 110
        temp = enorm(n, wa2);
        sgnorm = (gnorm / temp) / temp;

        // Test whether the scaled gradient direction is acceptable.
        alpha = 0.0;
        if (sgnorm < *delta)
        {
            // The scaled gradient direction is not acceptable.
            // Finally, calculate the point along the dogleg
            // at which the quadratic is minimized.
            bnorm = enorm(n, qtb);
            temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm / *delta);
            temp = temp - (*delta / qnorm)*pow(sgnorm / *delta, 2.0)
                        + sqrt(pow(temp - (*delta / qnorm), 2.0)
                               + (
                                   (1.0 - pow(*delta / qnorm, 2.0) )
                                  *(1.0 - pow((sgnorm / *delta), 2.0) )
                                 )
                              );
            alpha = ( (*delta/qnorm)*(1.0 - pow(sgnorm / *delta, 2.0)) )/temp;
        }
    }
    // 120
    // Form appropriate convex combination of the Gauss-Newton direction
    // and the scaled gradient direction.
    temp = (1.0 - alpha)*fmin(sgnorm, *delta);
    for (j = 0; j < n; j++) {
        x[j] = temp*wa1[j] + alpha*x[j];
    }
    // 130
    return;
}


double enorm(const int n, const double* x)
{
    //     **********
    //
    //     function enorm
    //
    //     given an n-vector x, this function calculates the
    //     euclidean norm of x.
    //
    //     the euclidean norm is computed by accumulating the sum of
    //     squares in three different sums. the sums of squares for the
    //     small and large components are scaled so that no overflows
    //     occur. non-destructive underflows are permitted. underflows
    //     and overflows do not occur in the computation of the unscaled
    //     sum of squares for the intermediate components.
    //     the definitions of small, intermediate and large components
    //     depend on two constants, rdwarf and rgiant. the main
    //     restrictions on these constants are that rdwarf**2 not
    //     underflow and rgiant**2 not overflow. the constants
    //     given here are suitable for every known computer.
    //
    //     the function statement is
    //
    //       double precision function enorm(n,x)
    //
    //     where
    //
    //       n is a positive integer input variable.
    //
    //       x is an input array of length n.
    //
    //     subprograms called
    //
    //       fortran-supplied ... dabs,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i;
    double agiant, s1, s2, s3, xabs, x1max, x3max;

    const double rgiant = pow(2.0, 63.5);
    const double rdwarf = 0.5 / rgiant;

    s1 = 0.0;
    s2 = 0.0;
    s3 = 0.0;
    x1max = 0.0;
    x3max = 0.0;
    agiant = rgiant / n;

    for (i = 0; i < n; i++)
    {
        xabs = fabs(x[i]);
        if ((xabs <= rdwarf) || (xabs >= agiant))
        {
            if (xabs > rdwarf)
            {
                // sum for large components
                if (xabs > x1max)
                {
                    s1 = 1.0 + s1*pow(x1max / xabs, 2.0);
                    x1max = xabs;
                } else {
                    // 10
                    s1 += pow(xabs / x1max, 2.0);
                }
                // 20
            } else {
                // 30

                //sum for small components
                if (xabs > x3max)
                {
                    s3 = 1.0 + s3*pow(x3max / xabs, 2.0);
                    x3max = xabs;
                } else {
                    // 40
                    if (xabs != 0.0)
                    {
                        s3 += pow(xabs / x3max, 2.0);
                    }
                }
                // 50
            }
        } else {
            // 70

            // sum for intermediate components
            s2 += xabs*xabs;
        }
        // 80
        // 90
    }

    // calculation of norm.
    if (s1 != 0.0)
    {
        return x1max*sqrt(s1+(s2/x1max)/x1max);
    }
    // 100
    if (s2 != 0.0)
    {
        if (s2 >= x3max)
        {
            return sqrt(s2*(1.0 + (x3max/s2)*(x3max*s3)));
        } else {
            return sqrt(x3max*((s2/x3max)+(x3max*s3)));
        }
    }
    // 110

    return x3max*sqrt(s3);
}


void fdjac1(int(*fcn)(int* n, double* x, double* fvec, int* iflag), const int n, double* x,
            const double* fvec, double* fjac, const int ldfjac, int* iflag, const int ml,
            const int mu, const double epsfcn, double* wa1, double* wa2)
{
    //     **********
    //
    //     subroutine fdjac1
    //
    //     this subroutine computes a forward-difference approximation
    //     to the n by n jacobian matrix associated with a specified
    //     problem of n functions in n variables. if the jacobian has
    //     a banded form, then function evaluations are saved by only
    //     approximating the nonzero terms.
    //
    //     the subroutine statement is
    //
    //       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
    //                         wa1,wa2)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions. fcn must be declared
    //         in an external statement in the user calling
    //         program, and should be written as follows.
    //
    //         subroutine fcn(n,x,fvec,iflag)
    //         integer n,iflag
    //         double precision x(n),fvec(n)
    //         ----------
    //         calculate the functions at x and
    //         return this vector in fvec.
    //         ----------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of fdjac1.
    //         in this case set iflag to a negative integer.
    //
    //       n is a positive integer input variable set to the number
    //         of functions and variables.
    //
    //       x is an input array of length n.
    //
    //       fvec is an input array of length n which must contain the
    //         functions evaluated at x.
    //
    //       fjac is an output n by n array which contains the
    //         approximation to the jacobian matrix evaluated at x.
    //
    //       ldfjac is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array fjac.
    //
    //       iflag is an integer variable which can be used to terminate
    //         the execution of fdjac1. see description of fcn.
    //
    //       ml is a nonnegative integer input variable which specifies
    //         the number of subdiagonals within the band of the
    //         jacobian matrix. if the jacobian is not banded, set
    //         ml to at least n - 1.
    //
    //       epsfcn is an input variable used in determining a suitable
    //         step length for the forward-difference approximation. this
    //         approximation assumes that the relative errors in the
    //         functions are of the order of epsfcn. if epsfcn is less
    //         than the machine precision, it is assumed that the relative
    //         errors in the functions are of the order of the machine
    //         precision.
    //
    //       mu is a nonnegative integer input variable which specifies
    //         the number of superdiagonals within the band of the
    //         jacobian matrix. if the jacobian is not banded, set
    //         mu to at least n - 1.
    //
    //       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
    //         least n, then the jacobian is considered dense, and wa2 is
    //         not referenced.
    //
    //     subprograms called
    //
    //       minpack-supplied ... dpmpar
    //
    //       fortran-supplied ... dabs,dmax1,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, k, msum, mut_n;
    double h, temp;
    const double epsmch = dpmpar[0];
    const double eps = sqrt(fmax(epsfcn, epsmch));
    mut_n = n;
    msum = ml + mu + 1;
    if (msum >= n)
    {
        // Computation of dense approximate jacobian.
        for (j = 0; j < n; j++)
        {
            temp = x[j];
            h = eps*fabs(temp);
            if (h == 0.0) { h = eps; }
            x[j] = temp + h;
            (*fcn)(&mut_n, x, wa1, iflag);
            if (*iflag < 0) { return; }
            x[j] = temp;
            for (i = 0; i < n; i++)
            {
                fjac[i + ldfjac*j] = (wa1[i] - fvec[i])/h;
            }
            // 10
        }
        // 20
        return;
    }
    // 40

    // Computation of banded approximate jacobian.
    for (k = 1; k <= msum; k++)
    {
	    for (j = k; j <= n; j += msum)
        {
            wa2[j-1] = x[j-1];
            h = eps*fabs(wa2[j-1]);
            if (h == 0.0) { h = eps; }
            x[j-1] = wa2[j-1] + h;
        }
        // 60
        (*fcn)(&mut_n, x, wa1, iflag);
        if (*iflag < 0) { return; }

        for (j = k; j <= n; j += msum)
        {
            x[j-1] = wa2[j-1];
            h = eps * fabs(wa2[j-1]);
            if (h == 0.) { h = eps; }
            for (i = 1; i <= n; i++)
            {
                fjac[(i-1) + ldfjac*(j-1)] = 0.0;
                if ((i >= j - mu) && (i <= j + ml))
                {
                    fjac[(i-1) + ldfjac*(j-1)] = (wa1[i-1] - fvec[i-1]) / h;
                }
            }
            // 70
        }
        // 80
	}
    // 90
    return;
}


void fdjac2(int(*fcn)(int *m, int *n, double* x, double* fvec, int* iflag),
            const int m, const int n, double* x, const double* fvec,
            double* fjac, const int ldfjac, int* iflag, const double epsfcn,
            double* wa)
{
    //     **********
    //
    //     subroutine fdjac2
    //
    //     this subroutine computes a forward-difference approximation
    //     to the m by n jacobian matrix associated with a specified
    //     problem of m functions in n variables.
    //
    //     the subroutine statement is
    //
    //       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
    //
    //     where
    //
    //       fcn is the name of the user-supplied subroutine which
    //         calculates the functions. fcn must be declared
    //         in an external statement in the user calling
    //         program, and should be written as follows.
    //
    //         subroutine fcn(m,n,x,fvec,iflag)
    //         integer m,n,iflag
    //         double precision x(n),fvec(m)
    //         ----------
    //         calculate the functions at x and
    //         return this vector in fvec.
    //         ----------
    //         return
    //         end
    //
    //         the value of iflag should not be changed by fcn unless
    //         the user wants to terminate execution of fdjac2.
    //         in this case set iflag to a negative integer.
    //
    //       m is a positive integer input variable set to the number
    //         of functions.
    //
    //       n is a positive integer input variable set to the number
    //         of variables. n must not exceed m.
    //
    //       x is an input array of length n.
    //
    //       fvec is an input array of length m which must contain the
    //         functions evaluated at x.
    //
    //       fjac is an output m by n array which contains the
    //         approximation to the jacobian matrix evaluated at x.
    //
    //       ldfjac is a positive integer input variable not less than m
    //         which specifies the leading dimension of the array fjac.
    //
    //       iflag is an integer variable which can be used to terminate
    //         the execution of fdjac2. see description of fcn.
    //
    //       epsfcn is an input variable used in determining a suitable
    //         step length for the forward-difference approximation. this
    //         approximation assumes that the relative errors in the
    //         functions are of the order of epsfcn. if epsfcn is less
    //         than the machine precision, it is assumed that the relative
    //         errors in the functions are of the order of the machine
    //         precision.
    //
    //       wa is a work array of length m.
    //
    //     subprograms called
    //
    //       user-supplied ...... fcn
    //
    //       minpack-supplied ... dpmpar
    //
    //       fortran-supplied ... dabs,dmax1,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, mut_m, mut_n;
    double h, temp;
    const double epsmch = dpmpar[0];
    const double eps = sqrt(fmax(epsfcn, epsmch));

    for (j = 0; j < n; j++)
    {
        temp = x[j];
        h = eps*fabs(temp);
        if (h == 0.0) { h = eps; }
        mut_m = m;
        mut_n = n;
        x[j] = temp + h;
        (*fcn)(&mut_m, &mut_n, x, wa, iflag);
        if (*iflag < 0) { return; }
        x[j] = temp;
        for (i = 0; i < m; i++)
        {
            fjac[i + ldfjac*j] = (wa[i] - fvec[i])/h;
        }
        // 10
    }
    // 20

}


void lmpar(const int n, double *r, const int ldr, const int* ipvt, const double* diag,
           const double* qtb, const double delta, double* par, double* x, double* sdiag,
           double* wa1, double* wa2)
{
    //     **********
    //
    //     subroutine lmpar
    //
    //     given an m by n matrix a, an n by n nonsingular diagonal
    //     matrix d, an m-vector b, and a positive number delta,
    //     the problem is to determine a value for the parameter
    //     par such that if x solves the system
    //
    //           a*x = b ,     sqrt(par)*d*x = 0 ,
    //
    //     in the least squares sense, and dxnorm is the euclidean
    //     norm of d*x, then either par is zero and
    //
    //           (dxnorm-delta) .le. 0.1*delta ,
    //
    //     or par is positive and
    //
    //           abs(dxnorm-delta) .le. 0.1*delta .
    //
    //     this subroutine completes the solution of the problem
    //     if it is provided with the necessary information from the
    //     qr factorization, with column pivoting, of a. that is, if
    //     a*p = q*r, where p is a permutation matrix, q has orthogonal
    //     columns, and r is an upper triangular matrix with diagonal
    //     elements of nonincreasing magnitude, then lmpar expects
    //     the full upper triangle of r, the permutation matrix p,
    //     and the first n components of (q transpose)*b. on output
    //     lmpar also provides an upper triangular matrix s such that
    //
    //            t   t                   t
    //           p *(a *a + par*d*d)*p = s *s .
    //
    //     s is employed within lmpar and may be of separate interest.
    //
    //     only a few iterations are generally needed for convergence
    //     of the algorithm. if, however, the limit of 10 iterations
    //     is reached, then the output par will contain the best
    //     value obtained so far.
    //
    //     the subroutine statement is
    //
    //       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
    //                        wa1,wa2)
    //
    //     where
    //
    //       n is a positive integer input variable set to the order of r.
    //
    //       r is an n by n array. on input the full upper triangle
    //         must contain the full upper triangle of the matrix r.
    //         on output the full upper triangle is unaltered, and the
    //         strict lower triangle contains the strict upper triangle
    //         (transposed) of the upper triangular matrix s.
    //
    //       ldr is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array r.
    //
    //       ipvt is an integer input array of length n which defines the
    //         permutation matrix p such that a*p = q*r. column j of p
    //         is column ipvt(j) of the identity matrix.
    //
    //       diag is an input array of length n which must contain the
    //         diagonal elements of the matrix d.
    //
    //       qtb is an input array of length n which must contain the first
    //         n elements of the vector (q transpose)*b.
    //
    //       delta is a positive input variable which specifies an upper
    //         bound on the euclidean norm of d*x.
    //
    //       par is a nonnegative variable. on input par contains an
    //         initial estimate of the levenberg-marquardt parameter.
    //         on output par contains the final estimate.
    //
    //       x is an output array of length n which contains the least
    //         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
    //         for the output par.
    //
    //       sdiag is an output array of length n which contains the
    //         diagonal elements of the upper triangular matrix s.
    //
    //       wa1 and wa2 are work arrays of length n.
    //
    //     subprograms called
    //
    //       minpack-supplied ... dpmpar,enorm,qrsolv
    //
    //       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, iter, j, k, l, nsing;
    double dxnorm, dwarf, fp, gnorm, parc, parl, paru, ssum, temp;

    dwarf = dpmpar[1];

    // Compute and store in x the Gauss-Newton direction. If the jacobian
    // is rank-deficient, obtain a least squares solution.

    nsing = n;  // Used as an unattainable index, signaling nonsingularity.
    for (j = 0; j < n; j++) {
	    wa1[j] = qtb[j];
	    if ((r[j + j*ldr] == 0.0) && (nsing == n)) { nsing = j; }
	    if (nsing < n) { wa1[j] = 0.0; }
    }
    // 10
    if (nsing > 0) {
        for (k = 1; k <= nsing; k++) {
            j = nsing - k;
            wa1[j] /= r[j + j*ldr];
            temp = wa1[j];
            if (j >= 1) {
                for (i = 0; i < j; i++) {
                    wa1[i] -= r[i + j*ldr]*temp;
                }
                // 20
            }
            // 30
        }
        // 40
    }
    // 50
    for (j = 0; j < n; j++) {
	    l = ipvt[j];
	    x[l] = wa1[j];
    }
    // 60

    // Initialize the iteration counter.
    // Evaluate the function at the origin and test for acceptance
    // of the Gauss-Newton direction.
    iter = 0;
    for (j = 0; j < n; j++)
    {
	    wa2[j] = diag[j] * x[j];
    }
    // 70

    dxnorm = enorm(n, wa2);
    fp = dxnorm - delta;
    if (fp <= 0.1*delta)
    {
        // 220
        if (iter == 0) { *par = 0.0; }
        return;
    }

    // If the jacobian is not rank deficient, the newton
    // step provides a lower bound, parl, for the zero of
    // the function. otherwise set this bound to zero.
    parl = 0.0;
    if (nsing == n) {
        for (j = 0; j < n; j++)
        {
            l = ipvt[j];
            wa1[j] = diag[l]*(wa2[l]/dxnorm);
        }
        // 80
        for (j = 0; j < n; j++)
        {
            ssum = 0.0;
            if (j > 0) {
                for (i = 0; i < j; i++) {
                    ssum += r[i + j*ldr]*wa1[i];
                }
                // 90
            }
            // 100
            wa1[j] = (wa1[j] - ssum) / r[j + j*ldr];
        }
        // 110
        temp = enorm(n, wa1);
        parl = ((fp / delta) / temp) / temp;
    }
    // 120

    // Calculate an upper bound, paru, for the zero of the function.
    for (j = 0; j < n; j++) {
        ssum = 0.0;
        for (i = 0; i <= j; i++)
        {
            ssum += r[i + j*ldr] * qtb[i];
        }
        // 130
        l = ipvt[j];
        wa1[j] = ssum / diag[l];
    }
    // 140

    gnorm = enorm(n, wa1);
    paru = gnorm / delta;
    if (paru == 0.0) { paru = dwarf / fmin(delta, 0.1); }

    // If the input par lies outside of the interval (parl, paru),
    // set par to the closer endpoint.
    *par = fmin(fmax(*par, parl), paru);
    if (*par == 0.0) { *par = gnorm / dxnorm; }

    // Beginning of an iteration.
    // 150
    while (1)
    {
        iter++;

        // Evaluate the function at the current value of par.
        if (*par == 0.) { *par = fmax(dwarf, paru*0.001); }
        temp = sqrt(*par);

        for (j = 0; j < n; j++) {
            wa1[j] = temp * diag[j];
        }
        // 160

        qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);

        for (j = 0; j < n; j++) {
            wa2[j] = diag[j] * x[j];
        }
        // 170

        dxnorm = enorm(n, wa2);
        temp = fp;
        fp = dxnorm - delta;
        // If the function is small enough, accept the current value
        // of par. also test for the exceptional cases where parl
        // is zero or the number of iterations has reached 10.

        if ((fabs(fp) <= 0.1*delta) || ((parl == 0.0) && (fp <= temp) && (temp < 0.0)) || (iter == 10)) {
            // 220
            if (iter == 0) { *par = 0.0; }
            return;
        }

        // Compute the Newton correction.
        for (j = 0; j < n; j++)
        {
            l = ipvt[j];
            wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        // 180

        for (j = 0; j < n; j++)
        {
            wa1[j] /= sdiag[j];
            temp = wa1[j];
            if (j != n - 1)
            {
                for (i = j + 1; i < n; i++)
                {
                    wa1[i] -= r[i + j*ldr]*temp;
                }
            }
            // 200
        }
        // 210

        temp = enorm(n, wa1);
        parc = ((fp/delta)/temp)/temp;

        // Depending on the sign of the function, update parl or paru
        if (fp > 0.) { parl = fmax(parl, *par); }
        if (fp < 0.) { paru = fmin(paru, *par); }

        // compute an improved estimate for par.
        *par = fmax(parl, *par + parc);

        // end of an iteration.
    }

    if (iter == 0)
    {
        *par = 0.0;
    }
    return;
}


void qform(const int m, const int n, double* q, const int ldq, double* wa)
{
    //     **********
    //
    //     subroutine qform
    //
    //     this subroutine proceeds from the computed qr factorization of
    //     an m by n matrix a to accumulate the m by m orthogonal matrix
    //     q from its factored form.
    //
    //     the subroutine statement is
    //
    //       subroutine qform(m,n,q,ldq,wa)
    //
    //     where
    //
    //       m is a positive integer input variable set to the number
    //         of rows of a and the order of q.
    //
    //       n is a positive integer input variable set to the number
    //         of columns of a.
    //
    //       q is an m by m array. on input the full lower trapezoid in
    //         the first min(m,n) columns of q contains the factored form.
    //         on output q has been accumulated into a square matrix.
    //
    //       ldq is a positive integer input variable not less than m
    //         which specifies the leading dimension of the array q.
    //
    //       wa is a work array of length m.
    //
    //     subprograms called
    //
    //       fortran-supplied ... min0
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, k, l, minmn, np1;
    double ssum, temp;

    // Zero out upper triangle of q in the first min(m,n) columns.
    minmn = (m > n ? n : m);
    if (minmn > 1)
    {
        for (j = 1; j < minmn; j++)
        {
            for (i = 0; i < j; i++)
            {
                q[i + j*ldq] = 0.0;
            }
            //10
        }
        // 20
    }
    // 30

    // Initialize remaining columns to those of the identity matrix.
    np1 = n + 1;
    if (m >= np1)
    {
        for (j = n; j < m; j++)
        {
            for (i = 0; i < m; i++)
            {
                q[i + j*ldq] = 0.0;
            }
            // 40
            q[j + j*ldq] = 1.0;
        }
        // 50
    }
    // 60

    // Accumulate q from its factored form.
    for (l = 0; l < minmn; l++)
    {
        k = minmn - l - 1;
        for (i = k; i < m; i++)
        {
            wa[i] = q[i + ldq*k];
            q[i + ldq*k] = 0.0;
        }
        // 70
        q[k + ldq*k] = 1.0;
        if (wa[k] != 0.0)
        {
            for (j = k; j < m; j++)
            {
                ssum = 0.0;
                for (i = k; i < m; i++)
                {
                    ssum += q[i + ldq*j]*wa[i];
                }
                // 80
                temp = ssum / wa[k];
                for (i = k; i < m; i++)
                {
                    q[i + ldq*j] -= temp*wa[i];
                }
                // 90
            }
            // 100
        }
        // 110
    }
    // 120

    return;
}


void qrfac(const int m, const int n, double* a, const int lda, const int pivot,
           int* ipvt, double* rdiag, double* acnorm, double* wa)
{
    //     **********
    //
    //     subroutine qrfac
    //
    //     this subroutine uses householder transformations with column
    //     pivoting (optional) to compute a qr factorization of the
    //     m by n matrix a. that is, qrfac determines an orthogonal
    //     matrix q, a permutation matrix p, and an upper trapezoidal
    //     matrix r with diagonal elements of nonincreasing magnitude,
    //     such that a*p = q*r. the householder transformation for
    //     column k, k = 1,2,...,min(m,n), is of the form
    //
    //                           t
    //           i - (1/u(k))*u*u
    //
    //     where u has zeros in the first k-1 positions. the form of
    //     this transformation and the method of pivoting first
    //     appeared in the corresponding linpack subroutine.
    //
    //     the subroutine statement is
    //
    //       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
    //
    //     where
    //
    //       m is a positive integer input variable set to the number
    //         of rows of a.
    //
    //       n is a positive integer input variable set to the number
    //         of columns of a.
    //
    //       a is an m by n array. on input a contains the matrix for
    //         which the qr factorization is to be computed. on output
    //         the strict upper trapezoidal part of a contains the strict
    //         upper trapezoidal part of r, and the lower trapezoidal
    //         part of a contains a factored form of q (the non-trivial
    //         elements of the u vectors described above).
    //
    //       lda is a positive integer input variable not less than m
    //         which specifies the leading dimension of the array a.
    //
    //       pivot is a logical input variable. if pivot is set true,
    //         then column pivoting is enforced. if pivot is set false,
    //         then no column pivoting is done.
    //
    //       ipvt is an integer output array of length lipvt. ipvt
    //         defines the permutation matrix p such that a*p = q*r.
    //         column j of p is column ipvt(j) of the identity matrix.
    //         if pivot is false, ipvt is not referenced.
    //
    //       lipvt is a positive integer input variable. if pivot is false,
    //         then lipvt may be as small as 1. if pivot is true, then
    //         lipvt must be at least n.
    //
    //       rdiag is an output array of length n which contains the
    //         diagonal elements of r.
    //
    //       acnorm is an output array of length n which contains the
    //         norms of the corresponding columns of the input matrix a.
    //         if this information is not needed, then acnorm can coincide
    //         with rdiag.
    //
    //       wa is a work array of length n. if pivot is false, then wa
    //         can coincide with rdiag.
    //
    //     subprograms called
    //
    //       minpack-supplied ... dpmpar,enorm
    //
    //       fortran-supplied ... dmax1,dsqrt,min0
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, k, kmax, minmn;
    double ajnorm, epsmch, ssum, temp;

    epsmch = dpmpar[0];
    // Compute the initial column norms and initialize several arrays
    for (j = 0; j < n; j++)
    {
        acnorm[j] = enorm(m, &a[j*lda]);
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        if (pivot) { ipvt[j] = j; }
    }
    // 10

    // Reduce a to r with householder transformations.
    minmn = (m > n ? n : m);
    for (j = 0; j < minmn; j++)
    {
        if (pivot)
        {
            // Bring the column of largest norm into the pivot position.
            kmax = j;
            // Search for the largest norm.
            for (k = j; k < n; k++)
            {
                if (rdiag[k] > rdiag[kmax]) { kmax = k; }
            }
            if (kmax != j)
            {
                // Swap columns j and kmax
                for (i = 0; i < m; i++)
                {
                    temp = a[i + lda*j];
                    a[i + lda*j] = a[i + lda*kmax];
                    a[i + lda*kmax] = temp;
                }
                // Swap the column norm values of j and kmax.
                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                // Swap the pivot entries
                k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }
        // 40

        // Compute Householder transformation to reduce the j-th column
        // of "a" to a multiple of the j-th unit vector.
        ajnorm = enorm(m - j, &a[lda*j + j]);
        if (ajnorm == 0.0)
        {
            // 100
            rdiag[j] = -ajnorm;
            continue;
        }

        if (a[lda*j + j] < 0.0) { ajnorm = -ajnorm; }
        // Rescale column j.
        for (i = j; i < m; i++)
        {
            a[i + lda*j] /= ajnorm;
        }
        // 50

        a[j + lda*j] += 1.0;

        // Apply the transformation to the remaining columns and
        // update the norms

        // Any column left?
        if (j == n - 1)
        {
            // 100
            rdiag[j] = -ajnorm;
            continue;
        }

        for (k = j+1; k < n; k++)
        {
            ssum = 0.0;
            for (i = j; i < m; i++)
            {
                ssum += a[i + lda*j]*a[i + lda*k];
            }
            // 60
            temp = ssum / a[j + lda*j];
            for (i = j; i < m; i++)
            {
                a[i + lda*k] -= temp*a[i + lda*j];
            }
            // 70
            if ((pivot) && (rdiag[k] != 0.0))
            {
                temp = a[j + lda*k] / rdiag[k];
                rdiag[k] *= sqrt(fmax(0.0, 1.0 - temp*temp));
                if (0.05*pow(rdiag[k]/wa[k], 2.0) <= epsmch)
                {
                    rdiag[k] = enorm(m-j, &a[(j+1) + lda*k]);
                    wa[k] = rdiag[k];
                }
            }
            // 80
        }
        // 90
        rdiag[j] = -ajnorm;
    }
    // 110

    return;
}


void qrsolv(const int n, double* r, const int ldr, const int* ipvt,
            const double* diag, const double* qtb, double* x,
            double* sdiag, double* wa)
{
    //     **********
    //
    //     subroutine qrsolv
    //
    //     given an m by n matrix a, an n by n diagonal matrix d,
    //     and an m-vector b, the problem is to determine an x which
    //     solves the system
    //
    //           a*x = b ,     d*x = 0 ,
    //
    //     in the least squares sense.
    //
    //     this subroutine completes the solution of the problem
    //     if it is provided with the necessary information from the
    //     qr factorization, with column pivoting, of a. that is, if
    //     a*p = q*r, where p is a permutation matrix, q has orthogonal
    //     columns, and r is an upper triangular matrix with diagonal
    //     elements of nonincreasing magnitude, then qrsolv expects
    //     the full upper triangle of r, the permutation matrix p,
    //     and the first n components of (q transpose)*b. the system
    //     a*x = b, d*x = 0, is then equivalent to
    //
    //                  t       t
    //           r*z = q *b ,  p *d*p*z = 0 ,
    //
    //     where x = p*z. if this system does not have full rank,
    //     then a least squares solution is obtained. on output qrsolv
    //     also provides an upper triangular matrix s such that
    //
    //            t   t               t
    //           p *(a *a + d*d)*p = s *s .
    //
    //     s is computed within qrsolv and may be of separate interest.
    //
    //     the subroutine statement is
    //
    //       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
    //
    //     where
    //
    //       n is a positive integer input variable set to the order of r.
    //
    //       r is an n by n array. on input the full upper triangle
    //         must contain the full upper triangle of the matrix r.
    //         on output the full upper triangle is unaltered, and the
    //         strict lower triangle contains the strict upper triangle
    //         (transposed) of the upper triangular matrix s.
    //
    //       ldr is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array r.
    //
    //       ipvt is an integer input array of length n which defines the
    //         permutation matrix p such that a*p = q*r. column j of p
    //         is column ipvt(j) of the identity matrix.
    //
    //       diag is an input array of length n which must contain the
    //         diagonal elements of the matrix d.
    //
    //       qtb is an input array of length n which must contain the first
    //         n elements of the vector (q transpose)*b.
    //
    //       x is an output array of length n which contains the least
    //         squares solution of the system a*x = b, d*x = 0.
    //
    //       sdiag is an output array of length n which contains the
    //         diagonal elements of the upper triangular matrix s.
    //
    //       wa is a work array of length n.
    //
    //     subprograms called
    //
    //       fortran-supplied ... dabs,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, k, l, nsing;
    double ccos, ccotan, qtbpj, ssin, ssum, ttan, temp;

    // Copy r and (q transpose)*b to preserve input and initialize s.
    // In particular, save the diagonal elements of r in x.
    for (j = 0; j < n; j++) {
	    for (i = j; i < n; i++) {
	        r[i + j*ldr] = r[j + i*ldr];
	    }
        // 10
	    x[j] = r[j + j*ldr];
	    wa[j] = qtb[j];
    }
    // 20

    // Eliminate the diagonal matrix d using a givens rotation.
    for (j = 0; j < n; j++)
    {
        // Prepare the row of d to be eliminated, locating the
        // diagonal element using p from the qr factorization.
        l = ipvt[j];
        if (diag[l] != 0.0)
        {
            for (k = j; k < n; k++)
            {
                sdiag[k] = 0.0;
            }
            sdiag[j] = diag[l];

            // The transformations to eliminate the row of d
            // modify only a single element of (q.T b) beyond
            // the first n, which is initially 0.
            qtbpj = 0.0;
            for (k = j; k < n; k++)
            {
                // Determine a givens rotation which eliminates the
                // appropriate element in the current row of d
                if (sdiag[k] != 0.0)
                {
                    if (fabs(r[ldr*k + k]) < fabs(sdiag[k]))
                    {
                        ccotan = r[ldr*k + k] / sdiag[k];
                        ssin = 0.5 / sqrt(0.25 + 0.25*pow(ccotan, 2.0));
                        ccos = ssin*ccotan;
                    } else {
                        ttan = sdiag[k] / r[ldr*k + k];
                        ccos = 0.5 / sqrt(0.25 + 0.25*pow(ttan, 2.0));
                        ssin = ccos*ttan;
                    }

                    // Compute the modified diagonal element of r and the
                    // modified element of ((q.T b),0).
                    r[ldr*k + k] *= ccos;
                    r[ldr*k + k] += ssin*sdiag[k];
                    temp = ccos*wa[k] + ssin*qtbpj;
                    qtbpj = -ssin*wa[k] + ccos*qtbpj;
                    wa[k] = temp;

                    // Accumulate the transformation in the row of s.
                    if (k != n - 1)
                    {
                        for (i = k+1; i < n; i++)
                        {
                            temp = ccos*r[i + ldr*k] + ssin*sdiag[i];
                            sdiag[i] = -ssin*r[i + ldr*k] + ccos*sdiag[i];
                            r[i + ldr*k] = temp;
                        }
                        // 60
                    }
                }
                // 70
            }
            // 80
        }
        // 90

        // Store the diagonal element of s and restore the corresponding
        // diagonal element of r.
        sdiag[j] = r[ldr*j + j];
        r[ldr*j + j] = x[j];
    }
    // 100

    // Solve the triangular system for z. If the system is
    // singular, then obtain a least squares solution.

    nsing = n;
    for (j = 0; j < n; j++)
    {
        if ((sdiag[j] == 0.0) && (nsing == n)) { nsing = j; }
        if (nsing < n) { wa[j] = 0.; }
    }
    // 110
    if (nsing > 0) {
        for (k = 1; k <= nsing; k++) {
            j = nsing - k;
            ssum = 0.0;
            if (nsing >= j+1) {
                for (i = j+1; i < nsing; ++i) {
                    ssum += r[i + j*ldr] * wa[i];
                }
            }
            // 130
            wa[j] = (wa[j] - ssum) / sdiag[j];
        }
        // 140
    }
    // 150

    // Permute the components of z back to components of x.
    for (j = 0; j < n; j++)
    {
    	l = ipvt[j];
	    x[l] = wa[j];
    }

    return;
}


void r1mpyq(const int m, const int n, double* a, const int lda, const double* v,
            const double* w)
{
    //     **********
    //
    //     subroutine r1mpyq
    //
    //     given an m by n matrix a, this subroutine computes a*q where
    //     q is the product of 2*(n - 1) transformations
    //
    //           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    //
    //     and gv(i), gw(i) are givens rotations in the (i,n) plane which
    //     eliminate elements in the i-th and n-th planes, respectively.
    //     q itself is not given, rather the information to recover the
    //     gv, gw rotations is supplied.
    //
    //     the subroutine statement is
    //
    //       subroutine r1mpyq(m,n,a,lda,v,w)
    //
    //     where
    //
    //       m is a positive integer input variable set to the number
    //         of rows of a.
    //
    //       n is a positive integer input variable set to the number
    //         of columns of a.
    //
    //       a is an m by n array. on input a must contain the matrix
    //         to be postmultiplied by the orthogonal matrix q
    //         described above. on output a*q has replaced a.
    //
    //       lda is a positive integer input variable not less than m
    //         which specifies the leading dimension of the array a.
    //
    //       v is an input array of length n. v(i) must contain the
    //         information necessary to recover the givens rotation gv(i)
    //         described above.
    //
    //       w is an input array of length n. w(i) must contain the
    //         information necessary to recover the givens rotation gw(i)
    //         described above.
    //
    //     subroutines called
    //
    //       fortran-supplied ... dabs,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more
    //
    //     **********
    int i, j, k;
    double ccos, ssin, temp;

    if (n < 2) { return; }
    for (k = 0; k < (n-1); k++) {
	    j = n - k - 2;
        if (fabs(v[j]) > 1.0) {
            ccos = 1 / v[j];
            ssin = sqrt(1.0 - ccos*ccos);
        } else {
            ssin = v[j];
            ccos = sqrt(1.0 - ssin*ssin);
        }

        for (i = 0; i < m; i++) {
            temp = ccos * a[i + j*lda] - ssin*a[i + (n-1)*lda];
            a[i + (n-1)*lda] = ssin*a[i + j*lda] + ccos*a[i + (n-1)*lda];
            a[i + j*lda] = temp;
        }
    }
    // 20

    // Apply the second set of givens rotations to a.
    for (j = 0; j < (n-1); j++) {
        if (fabs(w[j]) > 1.) {
            ccos = 1 / w[j];
            ssin = sqrt(1 - ccos*ccos);
        } else {
            ssin = w[j];
            ccos = sqrt(1 - ssin*ssin);
        }

        for (i = 0; i < m; i++) {
            temp = ccos*a[i + lda*j] + ssin*a[i + (n-1)*lda];
            a[i + (n-1)*lda] = -ssin*a[i + j*lda] + ccos*a[i + (n-1)*lda];
            a[i + lda*j] = temp;
        }
        // 30
    }
    // 40

    return;
}


void r1updt(const int m, const int n, double* s, const double* u, double* v,
            double* w, int* sing)
{
    //     **********
    //
    //     subroutine r1updt
    //
    //     given an m by n lower trapezoidal matrix s, an m-vector u,
    //     and an n-vector v, the problem is to determine an
    //     orthogonal matrix q such that
    //
    //                   t
    //           (s + u*v )*q
    //
    //     is again lower trapezoidal.
    //
    //     this subroutine determines q as the product of 2*(n - 1)
    //     transformations
    //
    //           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    //
    //     where gv(i), gw(i) are givens rotations in the (i,n) plane
    //     which eliminate elements in the i-th and n-th planes,
    //     respectively. q itself is not accumulated, rather the
    //     information to recover the gv, gw rotations is returned.
    //
    //     the subroutine statement is
    //
    //       subroutine r1updt(m,n,s,ls,u,v,w,sing)
    //
    //     where
    //
    //       m is a positive integer input variable set to the number
    //         of rows of s.
    //
    //       n is a positive integer input variable set to the number
    //         of columns of s. n must not exceed m.
    //
    //       s is an array of length ls. on input s must contain the lower
    //         trapezoidal matrix s stored by columns. on output s contains
    //         the lower trapezoidal matrix produced as described above.
    //
    //       ls is a positive integer input variable not less than
    //         (n*(2*m-n+1))/2.
    //
    //       u is an input array of length m which must contain the
    //         vector u.
    //
    //       v is an array of length n. on input v must contain the vector
    //         v. on output v(i) contains the information necessary to
    //         recover the givens rotation gv(i) described above.
    //
    //       w is an output array of length m. w(i) contains information
    //         necessary to recover the givens rotation gw(i) described
    //         above.
    //
    //       sing is a logical output variable. sing is set true if any
    //         of the diagonal elements of the output s are zero. otherwise
    //         sing is set false.
    //
    //     subprograms called
    //
    //       minpack-supplied ... dpmpar
    //
    //       fortran-supplied ... dabs,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, kenneth e. hillstrom, jorge j. more,
    //     john l. nazareth
    //
    //     **********
    int i, j, jj, l, nmj, nm1;
    double ccos, ccotan, giant, ssin, ttan, tau, temp;
    giant = dpmpar[2];

    // Initialize the diagonal element pointer
    jj = (n*(2*m - n + 1))/2 - (m - n);

    // Move the nontrivial part of the last column of s into w
    l = jj;
    for (i = n; i <= m; i++) {
        w[i - 1] = s[l - 1];
        l++;
    }

    // Rotate the vector v into a multiple of the n-th unit vector
    // in such a way that a spike is introduced into w.
    nm1 = n - 1;
    if (nm1 >= 1)
    {
        for (nmj = 1; nmj <= nm1; nmj++) {
            j = n - nmj;
            jj -= m - j + 1;
            w[j-1] = 0.0;
            if (v[j-1] != 0.0)
            {
                // Determine a givens rotation which eliminates the
                // j-th element of v.
                if (fabs(v[n-1]) < fabs(v[j-1]))
                {
                    ccotan = v[n-1]/v[j-1];
                    ssin = 0.5 / sqrt(0.25 + 0.25*(ccotan*ccotan));
                    ccos = ssin*ccotan;
                    tau = 1.0;
                    if (fabs(ccos)*giant > 1.0) { tau = 1.0/ccos; }
                } else {
                    // 20
                    ttan = v[j-1] / v[n-1];
                    ccos = 0.5 / sqrt(0.25 + 0.25*(ttan*ttan));
                    ssin = ccos*ttan;
                    tau = ssin;
                }
                // 30
                // Apply the transformation to v and store the information
                // necessary to recover the givens rotation.
                v[n-1] = ssin*v[j-1] + ccos*v[n-1];
                v[j-1] = tau;

                // Apply the transformation to s and extend the spike in w.
                l = jj;
                for (i = j; i <= m; i++) {
                    temp = ccos * s[l-1] - ssin * w[i-1];
                    w[i-1] = ssin * s[l-1] + ccos * w[i-1];
                    s[l-1] = temp;
                    l++;
                }
                // 40
            }
            // 50
        }
        // 60
    }
    // 70

    // Add the spike from the rank-1 update to w.
    for (i = 0; i < m; i++)
    {
	    w[i] += v[n-1] * u[i];
    }
    // 80

    // Eliminate the spike
    *sing = 0;
    if (nm1 >= 1) {
        for (j = 1; j <= nm1; j++) {
            if (w[j-1] != 0.) {
                // Determine a givens rotation which eliminates the
                // j-th element of the spike.

                if (fabs(s[jj-1]) < fabs(w[j-1])) {
                    ccotan = s[jj-1] / w[j-1];
                    ssin = 0.5 / sqrt(0.25 + 0.25*pow(ccotan, 2.0));
                    ccos = ssin * ccotan;
                    tau = 1.0;
                    if (fabs(ccos)*giant > 1.0) { tau = 1 / ccos; }
                } else {
                    ttan = w[j-1] / s[jj-1];
                    ccos = 0.5 / sqrt(0.25 + 0.25*pow(ttan, 2.0));
                    ssin = ccos * ttan;
                    tau = ssin;
                }
                // 100

                // Apply the transformation to s and reduce the spike in w.
                l = jj;
                for (i = j; i <= m; i++) {
                    temp = ccos*s[l-1] + ssin*w[i-1];
                    w[i-1] = -ssin*s[l-1] + ccos*w[i-1];
                    s[l-1] = temp;
                    l++;
                }
                //110
                w[j-1] = tau;
            }
            // 120

            // Test for zero diagonal elements in the output s.
            if (s[jj-1] == 0.0)
            {
                *sing = 1;
            }
            jj += (m - j + 1);
        }
        // 130
    }
    // 140

    // Move w back into the last column of the output s.
    l = jj;
    for (i = n; i <= m; i++) {
	    s[l-1] = w[i-1];
	    l++;
    }

    if (s[jj-1] == 0.) {
	*sing = 1;
    }

    return;
}


void rwupdt(const int n, double* r, const int ldr, const double* w, double* b,
            double* alpha, double* ccos, double* ssin)
{
    //     **********
    //
    //     subroutine rwupdt
    //
    //     given an n by n upper triangular matrix r, this subroutine
    //     computes the qr decomposition of the matrix formed when a row
    //     is added to r. if the row is specified by the vector w, then
    //     rwupdt determines an orthogonal matrix q such that when the
    //     n+1 by n matrix composed of r augmented by w is premultiplied
    //     by (q transpose), the resulting matrix is upper trapezoidal.
    //     the matrix (q transpose) is the product of n transformations
    //
    //           g(n)*g(n-1)* ... *g(1)
    //
    //     where g(i) is a givens rotation in the (i,n+1) plane which
    //     eliminates elements in the (n+1)-st plane. rwupdt also
    //     computes the product (q transpose)*c where c is the
    //     (n+1)-vector (b,alpha). q itself is not accumulated, rather
    //     the information to recover the g rotations is supplied.
    //
    //     the subroutine statement is
    //
    //       subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin)
    //
    //     where
    //
    //       n is a positive integer input variable set to the order of r.
    //
    //       r is an n by n array. on input the upper triangular part of
    //         r must contain the matrix to be updated. on output r
    //         contains the updated triangular matrix.
    //
    //       ldr is a positive integer input variable not less than n
    //         which specifies the leading dimension of the array r.
    //
    //       w is an input array of length n which must contain the row
    //         vector to be added to r.
    //
    //       b is an array of length n. on input b must contain the
    //         first n elements of the vector c. on output b contains
    //         the first n elements of the vector (q transpose)*c.
    //
    //       alpha is a variable. on input alpha must contain the
    //         (n+1)-st element of the vector c. on output alpha contains
    //         the (n+1)-st element of the vector (q transpose)*c.
    //
    //       cos is an output array of length n which contains the
    //         cosines of the transforming givens rotations.
    //
    //       sin is an output array of length n which contains the
    //         sines of the transforming givens rotations.
    //
    //     subprograms called
    //
    //       fortran-supplied ... dabs,dsqrt
    //
    //     argonne national laboratory. minpack project. march 1980.
    //     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
    //     jorge j. more
    //
    //     **********
    int i, j;
    double ccotan, rowj, ttan, temp;

    for (j = 0; j < n; j++)
    {
        rowj = w[j];

        // Apply the previous transformations to r(i,j), i=1...j-1 and to w(j)
        if (j > 0)
        {
            for (i = 1; i < j; i++)
            {
                temp = ccos[i-1]*r[i-1 + j*ldr] + ssin[i-1]*rowj;
                rowj = -ssin[i-1]*r[i-1 + j*ldr] + ccos[i-1]*rowj;
                r[i-1 + j*ldr] = temp;
            }
            // 10
        }
        // 20

        // Determine a givens rotation which eliminates w(j)
        ccos[j] = 1.0;
        ssin[j] = 0.0;
        if (rowj == 0.0) { continue; }

        if (fabs(r[j*ldr + j]) < fabs(rowj))
        {
            ccotan = r[j*ldr + j] / rowj;
            ssin[j] = 0.5 / sqrt(0.25 + 0.25*pow(ccotan, 2.0));
            ccos[j] = ssin[j]*ccotan;
        } else {
            // 30
            ttan = rowj / r[j*ldr + j];
            ccos[j] = 0.5 / sqrt(0.25 + 0.25*pow(ttan, 2.0));
            ssin[j] = ccos[j]*ttan;
        }

        // Apply the current transformation to r(j, j), b(j) and alpha
        r[j*ldr + j] *= ccos[j];
        r[j*ldr + j] += ssin[j]*rowj;
        temp = ccos[j]*b[j] + ssin[j]*(*alpha);
        *alpha = -ssin[j]*b[j] + ccos[j]*(*alpha);
        b[j] = temp;
    }

    return;
}
