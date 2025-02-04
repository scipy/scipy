#include "__slsqp.h"

void nonnegative_lsq_imp(const int m, const int n, double* restrict a, double* restrict b, double* restrict x, double* restrict w, double* restrict zz, double* restrict work, int* restrict indices, const int maxiter, double* rnorm, int* info);
static void ldp(int m, int n, double* g, double* h, double* x, double* buffer, int* indices, double* xnorm, int* mode);
static void lsi(int me, int mg, int n, double* e, double* f, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);



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
 *  The buffer pointers:
 *  buffer[0]              : Lagrange multipliers (me + mg)
 *  buffer[me + mg]        : wb, Modified b vector (ma)
 *  buffer[me + mg + ma]   : wh, Modified h vector (mg)
 *  buffer[me + 2*mg + ma] : The remaining buffer for modified arrays and LSI problem
 *
 */
void
lsei(int ma, int me, int mg, int n,
     double* a, double* b, double* e, double* f, double* g, double* h,
     double* x, double* buffer, int* jw, double* xnorm, int* mode)
{
    int one = 1, tmp_int = 0;
    double done = 1.0, dmone = -1.0;
    const double epsmach = 2.220446049250313e-16;

    if (me > n) { *mode = 2; return; }

    //    [E]         [E2   R]
    //    [A] @ Q.T = [A2  A1]
    //    [G]         [G2  G1]
    // me = 0 skips the equality constraint related computations
    // RQ decomposition of equality constraint data E and application to A, G.
    dgerq2_(&me, &n, e, &me, buffer, &buffer[tmp_int], mode);

    // Right triangularize E and apply Q.T to A and G from the right.
    dormr2_("R", "T", &ma, &n, &me, e, &me, buffer, a, &ma, &buffer[tmp_int], mode);
    dormr2_("R", "T", &ma, &n, &me, e, &me, buffer, g, &mg, &buffer[tmp_int], mode);

    // Check the diagonal elements of E for rank deficiency.
    for (int i = 0; i < me; i++) {
        if (fabs(e[i + i*me]) < epsmach) { *mode = 6;return; }
    }
    // Solve E*x = f and modify b.
    // Note: RQ forms R at the right of E instead of [0, 0] position.
    dtrsv_("U", "N", "N", &me, &e[(n - me)*me], &ma, x, &one);

    // Zero out the lagrange multipliers for the inequality constraints.
    for (int i = me; i < mg; i++) { buffer[i] = 0.0; }

    // If the problem is fully equality constrained then we are done. No free
    // variables left.
    if (me == n) { *mode = 1;return; }

    // After the triangularization of E, the free variables are reduced by
    // the number of equality constraints. Hence the resulting LSI problem has
    // n - me inequality constraints.
    int nineq = n - me;

    // Compute the temporary vector wf = f - E2*x
    double* restrict wb = &buffer[me + mg];
    for (int i = 0; i < ma; i++) { wb[i] = b[i]; }
    tmp_int = (n - me);
    dgemv_("N", &ma, &tmp_int, &dmone, a, &ma, x, &one, &done, wb, &one);

    // Store the transformed A and G corresponding to the free variables












}


/*
 * Solve inequality constrained least squares problem
 *      min |Ax - b|  subject to Gx >= h
 *
 *
 *
 *
*/
void
lsi(int ma, int mg, int n, double* a, double* b, double* g, double* h,
    double* x, double* buffer, int* jw, double* xnorm, int* mode)
{
    int one = 1, tmp_int = 0;
    double done = 1.0, dmone = -1.0, tmp_dbl = 0.0;
    const double epsmach = 2.220446049250313e-16;
    // QR decomposition of a and application to b.
    // We use the unblocked versions of the LAPACK routines to avoid
    // allocating extra "work" memory for the blocked versions.
    tmp_int = (ma < n ? ma : n);
    dgeqr2_(&ma, &n, a, &ma, buffer, &buffer[tmp_int], mode);
    // Check the diagonal elements of R for rank deficiency.
    *mode = 5;
    for (int i = 0; i < tmp_int; i++) {
        if (fabs(a[i + i*ma]) < epsmach) { return; }
    }
    // Compute Q^T b
    dorm2r_("L", "T", &ma, &one, &ma, a, &ma, buffer, b, &ma, &buffer[tmp_int], mode);
    // Transform G and h to form the LDP problem.
    // Solve XR = G where R is the upper triangular matrix from the QR.
    // The result is stored in G.
    dtrsm_("R", "U", "N", "N", &mg, &n, &done, a, &ma, g, &mg);
    // h = h - Xf
    dgemv_("N", &mg, &n, &dmone, g, &mg, b, &one, &done, h, &one);

    // Solve the LDP problem
    ldp(mg, n, g, h, x, buffer, jw, xnorm, mode);
    if (*mode != 0) { return; }

    // Solve the original problem
    for (int i = 0; i < n; i++) { x[i] += b[i]; }
    dtrsv_("U", "N", "N", &n, a, &ma, x, &one);

    // If any, compute the norm of the tail of b and add to xnorm
    if (n < ma) {
        tmp_int = ma - n;
        tmp_dbl = dnrm2_(&tmp_int, &b[n], &one);
        *xnorm = sqrt((*xnorm)*(*xnorm) + tmp_dbl*tmp_dbl);
    }
    return;
}

/*
 * Solve least distance problem
 *  min (1/2)|x|^2  subject to  Gx >= h
 *
 * G is (m x n), h is (m)
 * buffer is at least (m+2)*(n+1) + 3*m)
 * indices is int(n)
 * x is (n)
 * xnorm is the norm of the solution if succeded
 * mode is the return code integer
 *
 * Mode return values
 *  0  : solution found
 *  1  : iteration count exceeded by nnls
 *  2  : inequality constraints incompatible
 * -1  : bad input dimensions
 *
*/
void
ldp(int m, int n, double* g, double* h, double* x,
    double* buffer, int* indices, double* xnorm, int* mode)
{
    int one = 1;
    double dzero = 0.0, rnorm = 0.0;
    // Check for inputs and initialize x
    if (n <= 0) { *mode = -1; return; }
    for (int i = 0; i < n; i++) { x[i] = 0.0; }
    if (m == 0) { *mode = 0; return; }

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
            a[i + j*m] = g[j + i*m];
        }
        // Last row
        a[n + j*m] = h[j];
    }
    // RHS is (n+1)
    for (int i = 0; i < n; i++) { b[i] = 0.0; }
    b[n] = 1.0;

    // Solve the dual problem
    nonnegative_lsq_imp(n+1, m, a, b, y, w, zz, work, indices, 3*m, &rnorm, mode);
    if (*mode == 1) { return; }
    if (rnorm <= 0.0) { *mode = 2; return; }

    // Solve the primal problem
    double fac = 1.0 - ddot_(&m, h, &one, y, &one);
    if ((1.0 + fac) - 1.0 == 0.0) { *mode = 2; return; }
    *mode = 0;
    fac = 1.0 / fac;
    dgemv_("T", &m, &n, &fac, g, &m, y, &one, &dzero, x, &one);
    *xnorm = dnrm2_(&n, x, &one);

    // Compute the lagrange multipliers for the primal problem
    for (int i = 0; i < m; i++) { buffer[i] = fac*y[i]; }
    return;
}
