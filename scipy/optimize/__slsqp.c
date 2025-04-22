#include "__slsqp.h"

void __nnls(const int m, const int n, double* restrict a, double* restrict b, double* restrict x, double* restrict w, double* restrict zz, int* restrict indices, const int maxiter, double* rnorm, int* info);
static void ldp(int m, int n, double* g, double* h, double* x, double* buffer, int* indices, double* xnorm, int* mode);
static void lsi(int me, int mg, int n, double* e, double* f, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);
static void lsei(int ma, int me, int mg, int n, double* a, double* b, double* e, double* f, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);
static void lsq(int m, int meq, int n, int augment, double aug_weight, double* Lf, double* gradx, double* C, double* d, double* xl, double* xu, double* x, double* y, double* buffer, int* jw, int* mode);
static void ldl_update(int n, double* a, double* z, double sigma, double* w);

/*
 * The main SLSQP function. The function argument naming in the Fortran code is
 * exceedingly inconsistent and very difficult to follow. Hence we adopted the
 * following naming convention in SLSQP and the nested function arguments:
 *
 *  - funx: The function value at the current point. (1)
 *  - gradx: The gradient of the function at the current point. (n)
 *  - C: The equality and inequality constraint normals. (m x n)
 *  - d: The  equality and inequality constraints, (m)
 *  - xl: The lower bounds on x, (n)
 *  - xu: The upper bounds on x, (n)
 *  - sol: The solution vector, (n)
 *  - mult: The Lagrange multipliers, (m + 2*n + 2)
 *  - buffer: A buffer to hold various intermediate arrays.
 *  - indices: An array to hold the indices of the active constraints. (m + 2*n + 2)
 *
 *  The buffer size should be greater than:
 *  n*(n+1)//2 + m + 4*n + 3                                           # SLSQP
 *  (n+1)*(n+2) + (n+1)*meq + m + (mineq + 2*n + 2)*(n+1) +  3*n + 3   # LSQ
 *  mineq + 2n + 2 + 2*meq + (n+1) + (mineq + 3n + 3)*(n + 1 - meq)    # LSEI
 *  (mineq + 2n + 2 + 2)*(n + 2) + mineq + 2n + 2                      # LDP
 *  mineq + 2n + 2                                                     # NNLS
 *
 *
 *  If applicable, the following are the problem matrix naming convention:
 *
 *  - A: The coefficient matrix of cost function |Ax - b|
 *  - b: The RHS of cost function |Ax - b|
 *  - E: The (E)quality constraint matrix of Ex = f
 *  - f: The equality constraint RHS of Ex = f
 *  - G: The inequality constraint matrix of Gx >= h
 *  - h: The inequality constraint RHS of Gx >= h
 *
 */
void
__slsqp_body(
    struct SLSQP_vars* S, double* funx, double* restrict gradx,
    double* restrict C, double* restrict d, double* restrict sol,
    double* restrict mult, double* restrict xl, double* restrict xu, double* buffer,
    int* indices)
{

    int one = 1, lda = (S->m > 0 ? S->m : 1);
    int j;
    double done = 1.0, dmone = -1.0, alfmin = 0.1;
    int n = S->n;
    int m = S->m;
    int n1 = n + 1;
    int n2 = n1*n/2;

    // Chop the buffer for various array pointers.
    double* restrict bfgs       = &buffer[0];
    double* restrict x0         = &buffer[n2];
    double* restrict mu         = &buffer[n2 + n];
    double* restrict s          = &buffer[n2 + n + m];
    double* restrict u          = &buffer[n2 + n + m + n1];
    double* restrict v          = &buffer[n2 + n + m + n1 + n1];
    double* restrict lsq_buffer = &buffer[n2 + n + m + n1 + n1 + n1];

    // The badlin flag keeps track whether the SQP problem on the current
    // iteration was inconsistent or not.
    int badlin = 0;

    // Fortran code uses reverse communication for the iterations hence it
    // needs to jump back to where it left off. Thus the goto statements are
    // kept as is. Fortunately, they do not overlap too much and have a relatively
    // clean separation.
    if (S->mode ==  0) { goto MODE0; }
    if (S->mode == -1) { goto MODEM1; }
    if (S->mode == 1) { goto MODE1; }

MODE0:
    // We always use inexact line search, since exact search is broken in the
    // original Fortran code.
    S->exact = 0; // (S->acc < 0.0 ? 1 : 0);
    S->acc = fabs(S->acc);
    S->tol = 10*S->acc;
    S->iter = 0;
    S->reset = 0;
    for (int i = 0; i < n; i++) { s[i] = 0.0; }
    for (int i = 0; i < m; i++) { mu[i] = 0.0; }

RESET_BFGS:
    // Reset the BFGS matrix stored in packed format
    S->reset++;
    if (S->reset > 5) { goto LABEL255;}
    for (int i = 0; i < n2; i++) { bfgs[i] = 0.0; }
    j = 0;
    for (int i = 0; i < n; i++)
    {
        bfgs[j] = 1.0;
        j += n - i;
    }
    // 120

ITER_START:
    // Main iteration: Search direction, steplength, LDL'-update
    // 130
    S->mode = 9;
    if (S->iter >= S->itermax) { return; }
    S->iter++;

    // Search direction as solution of the QP-problem
    for (int i = 0; i < n; i++)
    {
        u[i] = -sol[i] + xl[i] ;
        v[i] = -sol[i] + xu[i] ;
    }

    S->h4 = 1.0;
    // augment and aug_weight are not used and hence 0.
    lsq(m, S->meq, n, 0, 0, bfgs, gradx, C, d, u, v, s, mult, lsq_buffer, indices, &S->mode);

    // Augmented problem for inconsistent linearization

    // If it turns out that the original SQP problem is inconsistent,
    // disallow termination with convergence on this iteration,
    // even if the augmented problem was solved.
    badlin = 0;

    // If equality constraints are not full rank and all are equality constrained
    // then the problem is inconsistent.
    if ((S->mode == 6) && (n == S->meq)) { S->mode = 4;}

    // If inconsistency detected, we augment the problem and try again.
    // Fortran code augments the problem matrices by embedding them in larger
    // buffers then calls lsq. However, these matrices are then copied into
    // another buffer inside lsq hence we can let lsq insert into the second
    // buffer without modifying the original matrices. The only change lsq needs
    // is the weightvalue of the augmented variable which starts at 100 and
    // being multiplied by 10 on each iteration. Hence we only pass that value
    // with "aug_weight".
    if (S->mode == 4)
    {
        badlin = 1;
        // Reset the RHS of the constraints to zero of the augmented system.
        for (int i = 0; i < n; i++) { s[i] = 0.0; }
        S->h3 = 0.0;
        double rho = 100.0;
        S->inconsistent = 0;
        while (1)
        {
            lsq(m, S->meq, n, 1, rho, bfgs, gradx, C, d, u, v, s, mult, lsq_buffer, indices, &S->mode);
            S->h4 = 1.0 - s[n];
            if (S->mode == 4)
            {
                rho *= 10.0;
                S->inconsistent++;
                if (S->inconsistent > 5) { return; }
                continue;
            } else if (S->mode != 1) {
                return;
            }
            break;
        }
    } else if (S->mode != 1) {
        return;
    }

    // Update multipliers for L1-test
    for (int i = 0; i < n; i++) { v[i] = gradx[i]; }
    dgemv_("T", &m, &n, &dmone, C, &lda, mult, &one, &done, v, &one);

    S->f0 = *funx;
    for (int i = 0; i < n; i++) { x0[i] = sol[i]; }
    S->gs = ddot_(&n, gradx, &one, s, &one);
    S->h1 = fabs(S->gs);
    S->h2 = 0.0;
    for (int j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h3 = d[j];
        } else {
            S->h3 = 0.0;
        }
        S->h2 = S->h2 + fmax(-d[j], S->h3);
        S->h3 = fabs(mult[j]);
        mu[j] = fmax(S->h3, (mu[j] + S->h3)/2.0);
        S->h1 = S->h1 + S->h3*fabs(d[j]);
    }

    // Check convergence
    S->mode = 0;
    if ((S->h1 < S->acc) && (S->h2 < S->acc) && (!badlin) && (*funx == *funx)) { return; }
    S->h1 = 0.0;
    for (int j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h3 = d[j];
        } else {
            S->h3 = 0.0;
        }
        S->h1 += mu[j]*fmax(-d[j], S->h3);
    }
    // 180
    S->t0 = *funx + S->h1;
    S->h3 = S->gs - S->h1*S->h4;
    S->mode = 8;
    if (S->h3 >= 0.0) { goto RESET_BFGS; }

    // Line search with an L1 test function
    S->line = 0;
    S->alpha = 1.0;

    // Inexact line search
LINE_SEARCH:

    S->line++;
    S->h3 = (S->alpha) * (S->h3);
    dscal_(&n, &S->alpha, s, &one);
    for (int i = 0; i < n; i++) { sol[i] = x0[i]; }
    daxpy_(&n, &done, s, &one, sol, &one);

    S->mode = 1;
    return;

MODE1:

    S->t = *funx;
    for (int j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h1 = d[j];
        } else {
            S->h1 = 0.0;
        }
        S->t = S->t + mu[j]*fmax(-d[j], S->h1);
    }
    S->h1 = S->t - S->t0;

    if ((S->h1 > (S->h3 / 10.0)) && (S->line <= 10))
    {
        S->alpha = fmax(S->h3/(2.0*(S->h3 - S->h1)), alfmin);
        goto LINE_SEARCH;
    }

    // Check convergence
    S->h3 = 0.0;
    for (int j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h1 = d[j];
        } else {
            S->h1 = 0.0;
        }
        S->h3 = S->h3 + fmax(-d[j], S->h1);
    }
    if (
        ((fabs(*funx - S->f0) < S->acc) || (dnrm2_(&n, s, &one) < S->acc)) &&
        (S->h3 < S->acc) &&
        (!badlin) &&
        (*funx == *funx) // To filter for finite entries
    )
    {
        S->mode = 0;
        return;
    } else {
        S->mode = -1;
    }
    return;

LABEL255:
    // Check relaxed convergence in case of positive directional derivative
    S->h3 = 0.0;
    for (int j = 0; j < m; j++)
    {
        if (j < S->meq)
        {
            S->h1 = d[j];
        } else {
            S->h1 = 0.0;
        }
        S->h3 = S->h3 + fmax(-d[j], S->h1);
    }
    if (((fabs(*funx - S->f0) < S->tol) || (dnrm2_(&n, s, &one) < S->tol)) &&
        (S->h3 < S->tol) &&
        (!badlin) &&
        (*funx == *funx)
    )
    {
        S->mode = 0;
    } else {
        S->mode = 8;
    }
    return;

MODEM1:

    // Call Jacobian at current x

    // Update Cholesky factors of Hessian matrix modified by BFGS formula
    // u[i] = gradx[i] - C.T @ mult - v[i]

    for (int i = 0; i < n; i++) { u[i] = gradx[i]; }
    dgemv_("T", &m, &n, &dmone, C, &lda, mult, &one, &done, u, &one);
    for (int i = 0; i < n; i++)
    {
        u[i] = u[i] - v[i];
    }

    // L'*S
    for (int i = 0; i < n; i++) { v[i] = s[i]; }
    dtpmv_("L", "T", "U", &n, bfgs, v, &one);

    // D*L'*S
    j = 0;
    for (int i = 0; i < n; i++) {
        v[i] = bfgs[j]*v[i];
        j += n - i;
    }

    // L*D*L'*S
    dtpmv_("L", "N", "U", &n, bfgs, v, &one);

    S->h1 = ddot_(&n, s, &one, u, &one);
    S->h2 = ddot_(&n, s, &one, v, &one);
    S->h3 = 0.2*(S->h2);
    if (S->h1 < S->h3)
    {
        S->h4 = (S->h2 - S->h3) / (S->h2 - S->h1);
        S->h1 = S->h3;
        double tmp_dbl = 1.0 - S->h4;
        dscal_(&n, &S->h4, u, &one);
        daxpy_(&n, &tmp_dbl, v, &one, u, &one);
    }

    // Test for singular update, and reset hessian if so
    if ((S->h1 == 0.0) || (S->h2 == 0.0)) { goto RESET_BFGS; }

    ldl_update(n, bfgs, u, 1.0 / S->h1, v);
    ldl_update(n, bfgs, v, -1.0 / S->h2, u);

    // End of main iteration
    goto ITER_START;

    return;
}


/*
 *          min     |A*x - b|
 *        E*x = f
 *        G*x >= h
 *      xl <= x <= xu
 *
 * Problem data is kept in Lf, gradx, C, d, xl, xu arrays in a rather tedious
 * format. C(m, n) is the constraint normals, d(n) is the constraint bounds.
 * xl(n) and xu(n) are the lower and upper bounds on x.
 *
 * Lf is the LDL' factor of the BFGS matrix also holding the diagonal entries.
 *
 * NaN entries in xl, xu, signify unconstrained variables and hence not included.
 *
 * The C matrix, for a problem with all x bounds are given and finite,
 * broken into E and G as follows:
 *
 *                      ┌────┐    ┌────┐    ┌┐
 *                  meq │    │    │ E  │  = ││ f
 *                      │   ─┼────┼>   │    ││
 *                      ┼────┼    └────┘    └┘
 *                      │    │    ┌────┐    ┌┐
 *                      │    │    │    │    ││
 *      mineq = m - meq │   ─┼────┼>   │    ││
 *                      │    │    │    │    ││
 *                      │    │    │    │    ││
 *                      └────┘    │    │ >= ││
 *                        C       ┼────┼    ┼┼
 *                              n │  I │    ││  xl
 *                                ┼────┼    ┼┼
 *                              n │ -I │    ││ -xu
 *                                └────┘    └┘
 *                                   G       h
 *
 * A and b are stored in Lf[] in LAPACK packed format where Lf holds a unit, lower
 * triangular matrix with diagonal entries are overwritten by the entries of d[]
 * and vector and gradx[].
 *
 *  Lf[] = [d[0], s[1], s[2], . , d[1], s[n + 2], d[2], ...]
 *
 *  interpreted as:
 *
 *         [d[ 0 ],                          ]
 *         [s[ 1 ], d[ 1 ], .  ,             ]
 *  Lf[] = [s[ 2 ], s[n+2], .  ,             ]
 *         [ .    ,   .   , .  , d[n-1]      ]
 *         [s[ n ], s[2*n], .  ,   .   , d[n]]
 *
 * Then, the following relations recover the A and b
 *
 *          A = sqrt(d[]) * Lf[]^T
 *          b = - inv( Lf[] * sqrt(d[]) ) * gradx[]
 *
 * The solution is returned in x() and the Lagrange multipliers are returned in y().
 *
 * For solving the problem in case of a detection of inconsistent linearization,
 * see D. Kraft, "A software package for Sequential Quadratic Programming"
 * Section 2.2.3
 *
 * In the original code, the augmented system is detected by mismatch of certain
 * integers which is making things quite unreadable. Here we explicitly pass a
 * flag.
 *
 * Inconsistent linearization augments all arrays to accomodate for the dummy
 * variable. The function is still called with the original sizes but the flag
 * allows for enlarging the problem and hence the supplied buffer should accomodate
 * for this extra space.
 *
 * The required buffer size is given by:
 * (2*(m - meq)*(n + 1)+2)*(n - meq +1) + 2*2*(m-meq)*(n + 1) + 2*(m-meq)*(n + 1)
 *  + 2*meq + ld + (ld + 2*(m-meq)*(n + 1))*(n - meq)
 *
 */
void lsq(
    int m, int meq, int n, int augment, double aug_weight, double* restrict Lf,
    double* restrict gradx, double* restrict C, double* restrict d,
    double* restrict xl, double* restrict xu, double* restrict x,
    double* restrict y, double* buffer, int* jw, int* mode)
{
    int one = 1, orign = n;
    int mineq = m - meq;
    double xnorm = 0.0;
    int cursor = 0;
    int ld = n;
    int n_wG_rows = 0;

    if (augment) {
        ld = n + 1;
        x[n]     = 1.0;
        xl[n]    = 0.0;
        xu[n]    = 1.0;
    }

    // Recover A and b from Lf and gradx
    for (int i = 0; i < (ld+2)*ld; i++) { buffer[i] = 0.0; }
    double* restrict wA = buffer;
    double* restrict wb = &buffer[ld*(ld+1)];

    // Depending on augmented, wA is either the full array or the top-left block.

    for (int j = 0; j < n; j++)
    {
        double diag = sqrt(Lf[cursor++]);      // Extract the diagonal value from Lf.
        wA[j + j * ld] = diag;                 // Place the sqrt diagonal.
        for (int i = j + 1; i < n; i++)
        {
            wA[j + i * ld] = Lf[cursor++] * diag;
        }
    }

    // Compute b = - 1/sqrt(d[]) * inv(Lf[]) * gradx[]. Lf is already in packed format.
    for (int i = 0; i < n; i++) { wb[i] = gradx[i]; }
    dtpsv_("L", "N", "U", &n, Lf, wb, &one);
    cursor = 0;
    for (int i = 0; i < n; i++)
    {
        wb[i] /= -sqrt(Lf[cursor]);
        cursor += n - i;
    }

    // If augmented, fill in the extra entry in the bottom right corner.
    if (augment) { wA[ld*ld - 1] = aug_weight; }

    // If augmented, also increase the number of variables by 1.
    if (augment) { n++; }

    // Get the equality constraints if given.
    double* restrict wE = &buffer[n*(n+1) + n];
    double* restrict wf = &buffer[n*(n+1) + n + n*meq];
    if (meq > 0)
    {
        for (int j = 0; j < n-1; j++)
        {
            for (int i = 0; i < meq; i++)
            {
                wE[i + j*meq] = C[i + j*m];
            }
        }
        if (augment)
        {
            // n is incremented hence all Ceq is now in wE. Add the extra column.
            for (int i = 0; i < meq; i++) { wE[i + (n-1)*meq] = -d[i]; }

        } else  {
            // If not augmented then handle j = n - 1 that is skipped.
            for (int i = 0; i < meq; i++) { wE[i + (n-1)*meq] = C[i + (n-1)*m]; }

        }
        for (int i = 0; i < meq; i++) { wf[i] = -d[i]; }
    }

    // Get the inequality constraints if given. First zero out wG and wh.
    double* restrict wG = &buffer[n*(n+1) + n + n*meq + meq];
    double* restrict wh = &buffer[n*(n+1) + n + n*meq + meq + (mineq + 2*n)*ld];
    // Zero out wG and wh
    for (int i = 0; i < (mineq + 2*n)*(ld + 1); i++) { wG[i] = 0.0; }

    // Convert the bounds on x to +I and -I blocks in G.
    // Augment h by xl and -xu.
    // Unbounded constraints are signified by NaN values and they do not appear
    // in G and h. Hence there is a nancount tab to keep track of them.

    // We first populate "wh" to get the number of unbounded constraints. That will
    // define the unskipped row number of wG. This is different than the original
    // Fortran code where the max allocated row number and the actual row number
    // of wG has been kept separate and it causes to be sent to every nested
    // function call. Instead we form wG and wh once with fixed size.

    int nancount = 0;
    int nrow = mineq;
    if (m > meq)
    {
        for (int i = 0; i < mineq; i++) { wh[i] = -d[meq + i]; }
    }
    for (int i = 0; i < n; i++)
    {
        if (isnan(xl[i]))
        {
            nancount++;
        } else {
            wh[nrow++] = xl[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (isnan(xu[i]))
        {
            nancount++;
        } else {
            wh[nrow++] = -xu[i];
        }
    }

    n_wG_rows = mineq + 2*n - nancount;

    // Now that we know the actual row number of wG, we can finally populate
    // the top part with C.
    if (m > meq)
    {
        for (int j = 0; j < orign; j++)
        {
            for (int i = 0; i < mineq; i++)
            {
                wG[i + j*n_wG_rows] = C[meq + i + j*m];
            }
        }
    }

    // If augmented add the extra column.
    if (augment)
    {
        for (int i = 0; i < mineq; i++)
        {
            wG[i + orign*n_wG_rows] = fmax(-d[meq + i], 0.0);
        }
    }

    // Reset counter
    nrow = mineq;
    for (int i = 0; i < n; i++)
    {
        if (!isnan(xl[i]))
        {
            wG[nrow + i*n_wG_rows] = 1.0;
            nrow++;
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (!isnan(xu[i]))
        {
            wG[nrow + i*n_wG_rows] = -1.0;
            nrow++;
        }
    }

    // Assign the remaining part of the buffer to the LSEI problem.
    double* restrict lsei_scratch = &wh[mineq + 2*n];

    lsei(ld, meq, n_wG_rows, n, wA, wb, wE, wf, wG, wh, x, lsei_scratch, jw, &xnorm, mode);

    if (*mode == 1)
    {
        // Restore the Lagrange multipliers, first equality, then inequality.
        for (int i = 0; i < meq; i++) { y[i] = lsei_scratch[i+n_wG_rows]; }
        for (int i = 0; i < mineq; i++) { y[meq + i] = lsei_scratch[i]; }

        // Set the user-defined bounds on x to NaN
        for (int i = 0; i < 2*n; i++) { y[m + i] = NAN; }
    }

    // Clamp the solution, if given, to the finite bound interval
    for (int i = 0; i < n; i++)
    {
        if ((!isnan(xl[i])) && (x[i] < xl[i])) { x[i] = xl[i]; }
        else if ((!isnan(xu[i])) && (x[i] > xu[i])) { x[i] = xu[i]; }
    }

    return;
}


/*
 * Solve equality and inequality constrained least squares problem (LSEI)
 *      min |A*x - b|, subject to E*x = f, G*x >= h.
 *
 *  ma, me, mg : number of rows in A, E, G
 *  n          : number of columns in A, x
 *  a          : matrix A (ma x n)
 *  b          : vector b (ma)
 *  e          : matrix E (me x n)
 *  f          : vector f (me)
 *  g          : matrix G (mg x n)
 *  h          : vector h (mg)
 *  x          : solution vector x (n)
 *  buffer     : work buffer (mg + 2)*(n - me +1) + 3*mg + 2*me + ma + (ma + mg)*(n - me)
 *  jw         : integer work array
 *  xnorm      : norm of the solution
 *  mode       : return code
 *
 *  The buffer pointers that will be used:
 *  buffer[0]                : Lagrange multipliers (mg + me)
 *  buffer[mg + me]          : wb, Modified b vector (ma)
 *  buffer[mg + me + ma]     : tau, Pivots for the RQ decomposition of E (me)
 *  buffer[mg + 2*me + ma]   : Scratch space
 *
 */
void
lsei(int ma, int me, int mg, int n,
     double* restrict a, double* restrict b, double* restrict e,
     double* restrict f, double* restrict g, double* restrict h,
     double* restrict x, double* restrict buffer, int* jw,
     double* xnorm, int* mode)
{
    int one = 1, nvars = 0, info = 0, lde = 0, ldg = 0;
    double done = 1.0, dmone = -1.0, dzero = 0.0, t= 0.0;
    const double epsmach = 2.220446049250313e-16;

    for (int i = 0; i < n; i++) { x[i] = 0.0; }
    // Return if the problem is over-constrained.
    if (me > n) { *mode = 2; return; }

    //    [E]         [E2 |  R]                                [x ]
    //    [A] @ Q.T = [A2 | A1]  ,and, x is partitioned as x = [--]
    //    [G]         [G2 | G1]                                [xe]

    // me = 0 skips the equality constraint related computations even though it
    // causes aliasing below. The aliased arrays are not referenced in that case.
    // Use at least 1 for the leading dimension of E even when me = 0 for LAPACK
    // calls.
    nvars = (n - me);
    double* restrict gmults      = &buffer[0];
    double* restrict emults      = &buffer[mg];
    double* restrict wb          = &buffer[me + mg];
    double* restrict tau         = &buffer[me + mg + ma];
    double* restrict a2          = &buffer[mg + 2*me + ma];
    double* restrict g2          = &buffer[mg + 2*me + ma + ma*nvars];
    double* restrict lsi_scratch = &buffer[mg + 2*me + ma + (ma + mg)*nvars];

    // RQ decomposition of equality constraint data E and application to A, G.
    // LAPACK RQ routine dgerq2 forms R on the right.
    // dgeqr2 is the unblocked versions of dgeqrf without the memory allocation.
    // Use top of the yet unutilized scratch space for throw-away work.
    lde = (me > 0 ? me : 1);
    ldg = (mg > 0 ? mg : 1);
    dgerq2_(&me, &n, e, &lde, tau, lsi_scratch, &info);

    // Right triangularize E and apply Q.T to A and G from the right.
    dormr2_("R", "T", &ma, &n, &me, e, &lde, tau, a, &ma, lsi_scratch, &info);
    dormr2_("R", "T", &mg, &n, &me, e, &lde, tau, g, &ldg, lsi_scratch, &info);

    // Check the diagonal elements of E for rank deficiency.
    for (int i = 0; i < me; i++)
    {
        if (!(fabs(e[i + (nvars + i)*me]) >= epsmach)) { *mode = 6;return; }
    }
    // Solve E*x = f and modify b.
    // Note: RQ forms R at the right of E instead of [0, 0] position.
    for (int i = 0; i < me; i++) { x[nvars + i] = f[i]; }
    dtrsv_("U", "N", "N", &me, &e[(nvars)*me], &lde, &x[nvars], &one);

    *mode = 1;
    // Zero out the inequality multiplier.
    for (int i = 0; i < mg; i++) { gmults[i] = 0.0; }

    // If the problem is fully equality-constrained, revert the basis and return.
    if (me == n) { goto ORIGINAL_BASIS; }

    // Compute the modified RHS wb = b - A1*x
    // Copy b into wb
    for (int i = 0; i < ma; i++) { wb[i] = b[i]; }
    // Compute wb -= A1*xe
    dgemv_("N", &ma, &me, &dmone, &a[ma*nvars], &ma, &x[nvars], &one, &done, wb, &one);

    // Store the transformed A2 and G2 in the buffer
    for (int j = 0; j < nvars; j++)
    {
        for (int i = 0; i < ma; i++)
        {
            a2[i + j*ma] = a[i + j*ma];
        }
        for (int i = 0; i < mg; i++)
        {
            g2[i + j*mg] = g[i + j*mg];
        }
    }

    if (mg == 0)
    {
        // No inequality constraints, solve the least squares problem directly.
        // We deliberately use the unblocked algorithm to avoid allocation.
        int lwork = ma*nvars + 3*nvars + 1;
        // Save the RHS for residual computation
        double* restrict wb_orig = &lsi_scratch[lwork];
        for (int i = 0; i < ma; i++) { wb_orig[i] = wb[i]; }

        int krank = 0;
        t = sqrt(epsmach);
        dgelsy_(&ma, &nvars, &one, a2, &ma, wb, &ma, jw, &t, &krank, lsi_scratch, &lwork, &info);

        // Copy the solution to x
        for (int i = 0; i < nvars; i++) { x[i] = wb[i]; }

        // Compute the residual and its norm, use a since a2 is overwritten.
        dgemv_("N", &ma, &nvars, &done, a, &ma, x, &one, &dmone, wb_orig, &one);
        *xnorm = dnrm2_(&ma, wb_orig, &one);

        *mode = 7;
        if (krank < nvars) { return; }
        *mode = 1;
        goto ORIGINAL_BASIS;
    }

    // Modify h, and solve the inequality constrained least squares problem.
    // h -= G1*xe
    dgemv_("N", &mg, &me, &dmone, &g[mg*nvars], &ldg, &x[nvars], &one, &done, h, &one);

    lsi(ma, mg, nvars, a2, wb, g2, h, x, lsi_scratch, jw, xnorm, mode);

    // Copy multipliers from scratch to gmults
    for (int i = 0; i < mg; i++) { gmults[i] = lsi_scratch[i]; }

    // If no equality constraints this was an LSI problem all along.
    if (me == 0) { return; }

    t = dnrm2_(&me, &x[nvars], &one);
    // Modify the norm by adding the equality solution.
    *xnorm = hypot(*xnorm, t);
    if (*mode != 1) { return; }

ORIGINAL_BASIS:
    // Convert the solution and multipliers to the original basis.
    // b = A*x - b (residuals)
    dgemv_("N", &ma, &n, &done, a, &ma, x, &one, &dmone, b, &one);
    // f = A1^T*b - G1^T*w
    dgemv_("T", &ma, &me, &done, &a[nvars*ma], &ma, b, &one, &dzero, f, &one);
    dgemv_("T", &mg, &me, &dmone, &g[nvars*mg], &ldg, gmults, &one, &done, f, &one);

    // x = Q.T*x
    dormr2_("L", "T", &n, &one, &me, e, &lde, tau, x, &n, lsi_scratch, &info);

    // Solve the triangular system for the equality multipliers, emults.
    for (int i = 0; i < me; i++) { emults[i] = f[i]; }
    dtrsv_("U", "T", "N", &me, &e[(n - me)*me], &lde, emults, &one);

    return;
}


/*
 * Solve inequality constrained least squares problem
 *      min |Ax - b|  subject to Gx >= h
 *
 * A is (ma x n), b is (ma), G is (mg x n), h is (mg), x is (n)
 * buffer is at least (mg+2)*(n+1) + 2*mg
 * jw is at least (mg)
 * xnorm is the 2-norm of the residual vector
 * mode is the integer return code
 *
 * Return codes for mode
 *  1: successful computation
 *  2: error return because of wrong dimensions
 *  3: iteration count exceeded by nnls
 *  4: inequality constraints incompatible
 *  5: matrix A is not rank n
 *
*/
void
lsi(int ma, int mg, int n, double* restrict a, double* restrict b, double* restrict g,
    double* restrict h, double* restrict x, double* restrict buffer, int* jw,
    double* xnorm, int* mode)
{
    int one = 1, tmp_int = 0, info = 0;
    double done = 1.0, dmone = -1.0, tmp_dbl = 0.0;
    const double epsmach = 2.220446049250313e-16;

    // QR decomposition of A and application to b.
    // We use the unblocked versions of the LAPACK routines to avoid
    // allocating extra "work" memory for the blocked versions.
    tmp_int = (ma < n ? ma : n);
    dgeqr2_(&ma, &n, a, &ma, buffer, &buffer[tmp_int], &info);

    // Compute Q^T b
    dorm2r_("L", "T", &ma, &one, &tmp_int, a, &ma, buffer, b, &ma, &buffer[tmp_int], &info);

    // Check the diagonal elements of R for rank deficiency.
    *mode = 5;
    *xnorm = 0.0;
    for (int i = 0; i < tmp_int; i++) {
        if (!(fabs(a[i + i*ma]) >= epsmach)) { return; }
    }
    // Transform G and h to form the LDP problem.
    // Solve XR = G where R is the upper triangular matrix from the QR.
    // The result is stored in G.
    // Note: There is an inherent assumption that ma >= n. This is a bug carried
    // over here from the original slsqp implementation.
    dtrsm_("R", "U", "N", "N", &mg, &n, &done, a, &ma, g, &mg);
    // h = h - Xf
    dgemv_("N", &mg, &n, &dmone, g, &mg, b, &one, &done, h, &one);

    // Solve the LDP problem.
    ldp(mg, n, g, h, x, buffer, jw, xnorm, mode);
    if (*mode != 1) { return; }

    // Convert to the solution of the original problem.
    daxpy_(&n, &done, b, &one, x, &one);
    dtrsv_("U", "N", "N", &n, a, &ma, x, &one);

    // If any, compute the norm of the tail of b and add to xnorm
    tmp_int = ma - n;
    tmp_dbl = dnrm2_(&tmp_int, &b[(n + 1 > ma ? ma : n + 1) - 1], &one);
    *xnorm = hypot(*xnorm, tmp_dbl);

    return;
}

/*
 * Solve least distance problem
 *  min (1/2)|x|^2  subject to  Gx >= h
 *
 * G is (m x n), h is (m)
 * buffer is at least (m+2)*(n+1) + 2*m
 * indices is int(n)
 * x is (n)
 * xnorm is the norm of the solution if succeded
 * mode is the return code integer
 *
 * Mode return values
 *  1  : solution found
 *  2  : bad input dimensions
 *  3  : iteration count exceeded by nnls
 *  4  : inequality constraints incompatible
 *
*/
void
ldp(int m, int n, double* restrict g, double* restrict h, double* restrict x,
    double* restrict buffer, int* indices, double* xnorm, int* mode)
{
    int one = 1;
    double dzero = 0.0, rnorm = 0.0;
    // Check for inputs and initialize x
    if (n <= 0) { *mode = 2; return; }
    for (int i = 0; i < n; i++) { x[i] = 0.0; }
    if (m == 0) { *mode = 1; return; }

    // Define pointers for the variables on buffer
    double* restrict a    = &buffer[0];
    double* restrict b    = &buffer[m*(n+1)];
    double* restrict zz   = &buffer[(m+1)*(n+1)];
    double* restrict y    = &buffer[(m+2)*(n+1)];
    double* restrict w    = &buffer[(m+2)*(n+1) + m];

    // Save the dual problem data into buffer
    //       dual problem [G^T] [x] = [0]
    //                    [h^T]       [1]

    // LHS, G is (m x n), h is (m). Both transposed and stacked into (n+1) x m.
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            a[i + j*(n+1)] = g[j + i*m];
        }
        // Place h in the last row.
        a[n + j*(n+1)] = h[j];
    }
    // RHS is (n+1)
    for (int i = 0; i < n; i++) { b[i] = 0.0; }
    b[n] = 1.0;

    // Solve the dual problem
    __nnls(n+1, m, a, b, y, w, zz, indices, 3*m, &rnorm, mode);
    if (*mode != 1) { return; }
    *mode = 4;
    if (rnorm <= 0.0) { return; }

    // Solve the primal problem
    double fac = 1.0 - ddot_(&m, h, &one, y, &one);
    if (!((1.0 + fac) - 1.0 > 0.0)) { return; }
    *mode = 1;
    fac = 1.0 / fac;
    dgemv_("T", &m, &n, &fac, g, &m, y, &one, &dzero, x, &one);
    *xnorm = dnrm2_(&n, x, &one);

    // Compute the lagrange multipliers for the primal problem
    for (int i = 0; i < m; i++) { buffer[i] = fac*y[i]; }
    return;
}


/*
 *
 * Updates the LDL' factors of matrix a by rank-one matrix sigma*z*z'
 * n     : order of the coefficient matrix a
 * a     : positive definite matrix of dimension n; only the lower triangle is
 *         used and is stored column by column as one dimensional array of
 *         dimension n*(n+1)/2.
 * z     : vector of dimension n of updating elements
 * sigma : scalar factor by which the modifying dyade z*z' is multiplied
 * w     : working array of dimension n
 *
 * Uses the composite-t method of fletcher and powell as described in "On the
 * modification of LDL' factorizations", DOI:10.1090/S0025-5718-1974-0359297-1
 *
 * Implemented by: Dieter Kraft, dfvlr - Institut für Dynamik der Flugsysteme
 *                 D-8031  Oberpfaffenhofen
 *
 */
static void
ldl_update(int n, double* restrict a, double* restrict z, double sigma, double* restrict w)
{
    int j, ij = 0;
    const double epsmach = 2.220446049250313e-16;
    if (sigma == 0.0) { return; }
    double alpha, beta, delta, gamma, u, v, tp, t = 1.0 / sigma;

    if (sigma <= 0.0)
    {
        // Negative update
        for (int i = 0; i < n; i++) { w[i] = z[i]; }
        for (int i = 0; i < n; i++)
        {
            v = w[i];
            t = t + v*v/a[ij];
            for (int j = i + 1; j < n; j++)
            {
                ij++;
                w[j] = w[j] - v*a[ij];
            }
            ij++;
        }
        if (t >= 0.0) { t = epsmach / sigma; }

        for (int i = 0; i < n; i++)
        {
            j = n - i - 1;
            ij -= i + 1;
            u = w[j];
            w[j] = t;
            t = t - u*u / a[ij];
        }
    }

    // Positive update
    for (int i = 0; i < n; i++)
    {
        v = z[i];
        delta = v / a[ij];
        // sigma == 0.0 is handled at the beginning.
        tp = (sigma < 0.0 ? w[i] : t + delta*v);
        alpha = tp / t;
        a[ij] = alpha*a[ij];
        if (i == n - 1) { return; }
        beta = delta / tp;
        if (alpha <= 4.0)
        {
            for (int j = i + 1; j < n; j++)
            {
                ij++;
                z[j] = z[j] - v * a[ij];
                a[ij] = a[ij] + beta * z[j];
            }
        } else {
            gamma = t / tp;
            for (int j = i + 1; j < n; j++)
            {
                ij++;
                u = a[ij];
                a[ij] = gamma * u + beta * z[j];
                z[j] = z[j] - v * u;
            }
        }
        ij++;
        t = tp;
    }

    return;
}
