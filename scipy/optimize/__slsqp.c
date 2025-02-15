#include "__slsqp.h"

void __nnls(const int m, const int n, double* restrict a, double* restrict b, double* restrict x, double* restrict w, double* restrict zz, double* restrict work, int* restrict indices, const int maxiter, double* rnorm, int* info);
static void ldp(int m, int n, double* g, double* h, double* x, double* buffer, int* indices, double* xnorm, int* mode);
static void lsi(int me, int mg, int n, double* e, double* f, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);
static void lsei(int ma, int me, int mg, int n, double* a, double* b, double* e, double* f, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);
static void lsq(int m, int meq, int n, int nl, double* S, double* t, double* C, double* d, double* xl, double* xu, double* x, double* y, double* buffer, int* jw, int* mode);


void slsqp()
{
    // Nonlinear programming by solving sequentially quadratic programming

    // TODO: Implement the function
}

void slsqpb()
{
    // Nonlinear programming by solving sequentially quadratic programming

    // TODO: Implement the function
}

/*
 *          min     |A*x - b|
 *        E*x = f
 *        G*x >= h
 *      xl <= x <= xu
 *
 * Problem data is kept in S, t, C, d, xl, xu arrays in a rather tedious format.
 * C(m, n) is the constraint matrix, d(n) is the constraint bounds.
 * xl(n) and xu(n) are the lower and upper bounds on x.
 * NaN entries signify unbounded constraints and not included in the constraints.
 * The C matrix, for a problem with all x bounds are given and finite,
 * broken into E and G as follows:
 *
 *                      ┌────┐    ┌────┐  ┌┐
 *                  meq │    │    │ E  │  ││ f
 *                      │   ─┼────┼>   │  ││
 *                      ┼────┼    └────┘  └┘
 *                      │    │    ┌────┐  ┌┐
 *                      │    │    │    │  ││
 *      mineq = m - meq │   ─┼────┼>   │  ││
 *                      │    │    │    │  ││
 *                      │    │    │    │  ││
 *                      └────┘    │    │  ││
 *                        C       ┼────┼  ┼┼
 *                              n │  I │  ││  xl
 *                                ┼────┼  ┼┼
 *                              n │ -I │  ││ -xu
 *                                └────┘  └┘
 *                                   G     h
 *
 * A and b are stored in S[] in LAPACK packed format where S holds a unit, lower
 * triangular matrix with diagonal entries are overwritten by the entries of d[]
 * and vector and t[].
 *
 *  S[] = [d[0], s[1], s[2], . , d[1], s[n + 2], d[2], ...]
 *
 *        [d[ 0 ],                          ]
 *        [s[ 1 ], d[ 1 ], .  ,             ]
 *  S[] = [s[ 2 ], s[n+2], .  ,             ]
 *        [ .    ,   .   , .  , d[n-1]      ]
 *        [s[ n ], s[2*n], .  ,   .   , d[n]]
 *
 * Then, the following relations recover the A and b
 *
 *          A = sqrt(d[]) * S[]^T
 *          b = - inv( S[] * sqrt(d[]) ) * t[]
 *
 * The solution is returned in x() and the Lagrange multipliers are returned in y().
 *
 *
*/
void lsq(int m, int meq, int n, int nl, double* S, double* t, double* C, double* d,
         double* xl, double* xu, double* x, double* y, double* buffer, int* jw,
         int* mode)
{
    int one = 1;
    int mineq = m - meq;
    int aug = 0;
    double xnorm = 0.0;

    for (int i = 0; i < (n+2)*n; i++) { buffer[i] = 0.0; }
    double* restrict wA = buffer;
    double* restrict wb = &buffer[n*(n+1)];

    // Determine whether to solve problem with inconsistent linearization
    // See Kraft, "A software package for Sequential Quadratic Programming"
    // Section 2.2.3

    // Inconsistent linearization augments an extra column to A and extra row to
    // A and b. Then sends n value increased by 1 to lsq. For that we save the
    // size in ld and decrement n if aug is set to keep the problem size consistent.

    if ((n*(n+1)/2 + 1) != nl) { aug = 1; }

    // Recover A and b from S and t
    int cursor = 0;
    int ld = n;
    if (aug) { n--; }

    // Depending on aug, wA is either full (n)x(n) or top-left block of size (n-1)x(n-1).
    for (int j = 0; j < n; j++)
    {
        double diag = sqrt(S[cursor++]);      // Extract the diagonal value from S.
        wA[j + j * ld] = diag;                // Place the sqrt diagonal.
        for (int i = j + 1; i < n; i++)
        {
            wA[j + i * ld] = S[cursor++] * diag;
        }
    }

    // Compute b = - 1/sqrt(d) * inv(S) * t(). S is already in packed format.
    for (int i = 0; i < n; i++) { wb[i] = t[i]; }
    dtpsv_("L", "N", "U", &n, wA, wb, &one);
    cursor = 0;
    for (int i = 0; i < n; i++)
    {
        wb[i] /= -sqrt(S[cursor]);
        cursor += n - i;
    }

    // Fill in the augmented system extra entries.
    if (aug) { wA[ld*ld - 1] = S[ld*(ld+1)/2 + 1]; }

    // Get the equality constraints if given.
    double* restrict wE = &buffer[n*(n+1) + n];
    double* restrict wf = &buffer[n*(n+1) + n + n*meq];
    if (meq > 0)
    {
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < meq; i++)
            {
                wE[i + j*meq] = C[i + j*m];
            }
        }
        for (int i = 0; i < meq; i++) { wf[i] = d[i]; }
    }

    // Get the inequality constraints if given.
    // Also note that there is still one more variable that can come from the
    // inconsistent linearization hence in the full case we have m + 2n
    // constraints rendering G allocated size (mineq + 2n) x (ld).
    double* restrict wG = &buffer[n*(n+1) + n + n*meq + meq];
    double* restrict wh = &buffer[n*(n+1) + n + n*meq + meq + (mineq + 2*n)*ld];
    if (m > meq)
    {
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < mineq; i++)
            {
                wG[i + j*mineq] = C[meq + i + j*m];
            }
        }
        for (int i = 0; i < mineq; i++) { wh[i] = d[meq + i]; }
    }

    // Also the bottom block of G for state bounds is zeroed out.
    for (int i = 0; i < (2*n)*ld; i++) { wG[mineq*ld + i] = 0.0; }

    // Convert the bounds on x to +I and -I blocks in G.
    // Augment h by xl and -xu.
    // Unbounded constraints are signified by NaN values and they do not appear
    // in G and h. Hence there is a nancount tab to keep track of them.

    int nancount = 0;
    // Two counters to keep track of the current state and bound entry.
    int nstate = 0;
    int nrow = mineq*ld;
    for (int i = 0; i < n; i++)
    {
        if (isnan(xl[i]))
        {
            nancount++;
        } else {
            wG[nstate + nrow*ld] = 1.0;
            wh[mineq + nstate] = xl[i];
            nstate++;
            nrow++;
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (isnan(xl[i]))
        {
            nancount++;
        } else {
            wG[nstate + nrow*ld] = -1.0;
            wh[mineq + nstate] = -xu[i];
            nstate++;
            nrow++;
        }
    }

    // Solve the problem
    lsei(m, meq, mineq + 2*n - nancount, n, wA, wb, wE, wf, wG, wh, x, buffer, jw, &xnorm, mode);

    // Restore the Lagrange multipliers
    for (int i = 0; i < m; i++) { y[i] = wh[i]; }

    // Set the user-defined bounds on x to NaN
    for (int i = 0; i < 2*n; i++) { y[m + i] = NAN; }

    // Clamp the solution, if given, to the finite bound interval
    for (int i = 0; i < n; i++)
    {
        if ((!isnan(xl[i])) && (x[i] < xl[i])) { x[i] = xl[i]; }
        if ((!isnan(xu[i])) && (x[i] > xu[i])) { x[i] = xu[i]; }
    }
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
 *  buffer     : work buffer
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
     double* a, double* b, double* e, double* f, double* g, double* h,
     double* x, double* buffer, int* jw, double* xnorm, int* mode)
{
    int one = 1, nvars = 0, info = 0, lde = 0, ldg = 0;
    double done = 1.0, dmone = -1.0, dzero = 0.0, t= 0.0;
    const double epsmach = 2.220446049250313e-16;

    // Return if the problem is over-constrained.
    for (int i = 0; i < n; i++) { x[i] = 0.0; }
    if (me > n) { *mode = 2; return; }
    if (ma < n) { *mode = 5; return; }

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
        if (!(fabs(e[i + i*me]) >= epsmach)) { *mode = 6;return; }
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
        double* wb_orig = &lsi_scratch[lwork];
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
    *xnorm = sqrt((*xnorm)*(*xnorm) + t*t);
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
 * A is (ma x n), b is (ma), G is (mg x n), h is (mg)
 * buffer is ((mg+2)*(n+1) + 3*mg)
 * jw is (mg)
 * xnorm is the 2-norm of the residual vector
 * mode is the integer return code
 *
 * The following assumptions must hold
 *  - ma >= n
 *  - A is rank n
 *  - x size is (n)
 *
 * Return codes for mode
 *  1: successful computation
 *  2: error return because of wrong dimensions
 *  3: iteration count exceeded by nnls
 *  4: inequality constraints incompatible
 *  5: matrix e is not rank n
 *
*/
void
lsi(int ma, int mg, int n, double* a, double* b, double* g, double* h,
    double* x, double* buffer, int* jw, double* xnorm, int* mode)
{
    int one = 1, tmp_int = 0, info = 0;
    double done = 1.0, dmone = -1.0, tmp_dbl = 0.0;
    const double epsmach = 2.220446049250313e-16;

    // QR decomposition of a and application to b.
    // We use the unblocked versions of the LAPACK routines to avoid
    // allocating extra "work" memory for the blocked versions.
    tmp_int = (ma < n ? ma : n);
    dgeqr2_(&ma, &n, a, &ma, buffer, &buffer[tmp_int], &info);

    // Check the diagonal elements of R for rank deficiency.
    *mode = 5;
    *xnorm = 0.0;
    // Original code has a bug that it checks random entries and returns 5.
    if (ma < n) { return; }

    for (int i = 0; i < tmp_int; i++) {
        if (!(fabs(a[i + i*ma]) >= epsmach)) { return; }
    }
    // Compute Q^T b
    dorm2r_("L", "T", &ma, &one, &tmp_int, a, &ma, buffer, b, &ma, &buffer[tmp_int], &info);

    // Transform G and h to form the LDP problem.
    // Solve XR = G where R is the upper triangular matrix from the QR.
    // The result is stored in G.

    // Note: There is an inherent assumption that ma >= n. This is a bug carried
    // over here from the original slsqp implementation.
    dtrsm_("R", "U", "N", "N", &mg, &n, &done, a, &ma, g, &mg);
    // h = h - Xf
    dgemv_("N", &mg, &n, &dmone, g, &mg, b, &one, &done, h, &one);

    ldp(mg, n, g, h, x, buffer, jw, xnorm, mode);
    if (*mode != 1) { return; }

    // Solve the original problem
    for (int i = 0; i < n; i++) { x[i] += b[i]; }
    dtrsv_("U", "N", "N", &n, a, &ma, x, &one);

    // If any, compute the norm of the tail of b and add to xnorm
    tmp_int = ma - n;
    tmp_dbl = dnrm2_(&tmp_int, &b[(n + 1 > ma ? ma : n + 1) - 1], &one);
    *xnorm = sqrt((*xnorm)*(*xnorm) + tmp_dbl*tmp_dbl);

    return;
}

/*
 * Solve least distance problem
 *  min (1/2)|x|^2  subject to  Gx >= h
 *
 * G is (m x n), h is (m)
 * buffer is at least (m+2)*(n+1) + 3*m
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
ldp(int m, int n, double* g, double* h, double* x,
    double* buffer, int* indices, double* xnorm, int* mode)
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
    double* restrict work = &buffer[(m+2)*(n+1) + 2*m];

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
    __nnls(n+1, m, a, b, y, w, zz, work, indices, 3*m, &rnorm, mode);
    if (*mode != 1) { return; }
    *mode = 4;
    if (rnorm <= 0.0) { return; }

    // Solve the primal problem
    double fac = 1.0 - ddot_(&m, h, &one, y, &one);
    if ((1.0 + fac) - 1.0 == 0.0) { return; }
    *mode = 1;
    fac = 1.0 / fac;
    dgemv_("T", &m, &n, &fac, g, &m, y, &one, &dzero, x, &one);
    *xnorm = dnrm2_(&n, x, &one);

    // Compute the lagrange multipliers for the primal problem
    for (int i = 0; i < m; i++) { buffer[i] = fac*y[i]; }
    return;
}
