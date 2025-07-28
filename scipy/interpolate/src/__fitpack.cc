#include <string>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "__fitpack.h"

namespace fitpack{


/*
 * B-spline evaluation routine.
 */
void
_deBoor_D(const double *t, double x, int k, int ell, int m, double *result) {
    /*
     * On completion the result array stores
     * the k+1 non-zero values of beta^(m)_i,k(x):  for i=ell, ell-1, ell-2, ell-k.
     * Where t[ell] <= x < t[ell+1].
     */
    /*
     * Implements a recursive algorithm similar to the original algorithm of
     * deBoor.
     */
    double *hh = result + k + 1;
    double *h = result;
    double xb, xa, w;
    int ind, j, n;

    /*
     * Perform k-m "standard" deBoor iterations
     * so that h contains the k+1 non-zero values of beta_{ell,k-m}(x)
     * needed to calculate the remaining derivatives.
     */
    result[0] = 1.0;
    for (j = 1; j <= k - m; j++) {
        memcpy(hh, h, j*sizeof(double));
        h[0] = 0.0;
        for (n = 1; n <= j; n++) {
            ind = ell + n;
            xb = t[ind];
            xa = t[ind - j];
            if (xb == xa) {
                h[n] = 0.0;
                continue;
            }
            w = hh[n - 1]/(xb - xa);
            h[n - 1] += w*(xb - x);
            h[n] = w*(x - xa);
        }
    }

    /*
     * Now do m "derivative" recursions
     * to convert the values of beta into the mth derivative
     */
    for (j = k - m + 1; j <= k; j++) {
        memcpy(hh, h, j*sizeof(double));
        h[0] = 0.0;
        for (n = 1; n <= j; n++) {
            ind = ell + n;
            xb = t[ind];
            xa = t[ind - j];
            if (xb == xa) {
                h[m] = 0.0;
                continue;
            }
            w = j*hh[n - 1]/(xb - xa);
            h[n - 1] -= w;
            h[n] = w;
        }
    }
}


/*
 *  Find an interval such that t[interval] <= xval < t[interval+1].
 */
int64_t
_find_interval(const double* tptr, int64_t len_t,
               int k,
               double xval,
               int64_t prev_l,
               int extrapolate)
{
    ConstRealArray1D t = ConstRealArray1D(tptr, len_t);

    int64_t n = t.nelem - k - 1;
    double tb = t(k);
    double te = t(n);

    if (xval != xval) {
        // nan
        return -1;
    }

    if (((xval < tb) || (xval > te)) && !extrapolate) {
        return -1;
    }
    int64_t l = (k < prev_l) && (prev_l < n) ? prev_l : k;

    // xval is in support, search for interval s.t. t[interval] <= xval < t[l+1]
    while ((xval < t(l)) && (l != k)) {
        l -= 1;
    }

    l += 1;
    while ((xval >= t(l)) && (l != n)) {
        l += 1;
    }

    return l-1;
}


/*
 * Fill the (m, k+1) matrix of non-zero b-splines.
 *
 * A row gives b-splines which are non-zero at the corresponding value of `x`.
 * Also for each row store the `offset`: with full matrices, the non-zero elements
 * in row `i` start at `offset[i]`. IOW,
 * A_full[i, offset[i]: offset[i] + k + 1] = A_packed[i, :].
 *
 * What we construct here is `A_packed` and `offset` arrays.
 *
 * We also take into account possible weights for each `x` value: they
 * multiply the rows of the data matrix.
 *
 * To reconstruct the full dense matrix, `A_full`, we would need to know the
 * number of its columns, `nc`. So we return it, too.
 */
void
data_matrix( /* inputs */
            const double *xptr, int64_t m,      // x, shape (m,)
            const double *tptr, int64_t len_t,  // t, shape (len_t,)
            int k,
            const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
            int extrapolate,
            /* outputs */
            double *Aptr,                       // A, shape(m, k+1)
            int64_t *offset_ptr,                // offset, shape (m,)
            int64_t *nc,                        // the number of coefficient
            /* work array*/
            double *wrk)                        // work, shape (2k+2)
{
    auto x = ConstRealArray1D(xptr, m);
    auto t = ConstRealArray1D(tptr, len_t);
    auto w = ConstRealArray1D(wptr, m);
    auto A = RealArray2D(Aptr, m, k+1);
    auto offset = Array1D<int64_t, false>(offset_ptr, m);

    int64_t ind = k;
    for (int i=0; i < m; ++i) {
        double xval = x(i);

        // find the interval
        ind = _find_interval(t.data, len_t, k, xval, ind, extrapolate);
        if (!extrapolate && (ind < 0)){
            // should not happen here, validation is expected on the python side
            throw std::runtime_error("find_interval: out of bounds with x = " + std::to_string(xval));
        }
        offset(i) = ind - k;

        // compute non-zero b-splines
        _deBoor_D(t.data, xval, k, ind, 0, wrk);

        for (int64_t j=0; j < k+1; ++j) {
            A(i, j) = wrk[j] * w(i);
        }
    }

    *nc = len_t - k - 1;
}

/*
 * Constructs the data matrix A for periodic B-spline basis functions.
 *
 * Inputs:
 *  - xptr: array of x values, shape (m,)
 *  - m: number of data points
 *  - tptr: knot vector, shape (len_t,)
 *  - len_t: length of knot vector
 *  - k: spline degree
 *  - wptr: weights for each x value, shape (m,)
 *  - extrapolate: whether extrapolation outside knot range is allowed
 *
 * Outputs:
 *  - Aptr: output matrix A, shape (m, k+1), holds basis function values
 *  - H1ptr: periodic correction matrix H1, shape (m, k+1)
 *  - H2ptr: periodic correction matrix H2, shape (m, k)
 *  - offset_ptr: index array of shape (m,), each x lies in knot interval [offset(i) + k, offset(i) + k + 1)
 *  - nc: number of spline coefficients (len_t - k - 1)
 *
 * Work array:
 *  - wrk: temporary workspace of size (2k + 2), used by De Boor’s algorithm
 */
void
data_matrix_periodic( /* inputs */
            const double *xptr, int64_t m,
            const double *tptr, int64_t len_t,
            int k,
            const double *wptr,
            int extrapolate,
            /* outputs */
            double *Aptr,
            double *H1ptr,
            double *H2ptr,
            int64_t *offset_ptr,
            int64_t *nc,
            /* work array*/
            double *wrk)
{
    // Wrap raw input pointers into readable array types
    auto x = ConstRealArray1D(xptr, m);
    auto t = ConstRealArray1D(tptr, len_t);
    auto w = ConstRealArray1D(wptr, m);
    auto A = RealArray2D(Aptr, m, k + 1);
    auto H1 = RealArray2D(H1ptr, m, k + 1);
    auto H2 = RealArray2D(H2ptr, m, k);
    auto offset = Array1D<int64_t, false>(offset_ptr, m);

    // Initialize H1 and H2 to zeros
    for( int64_t i = 0; i < m; i++ ) {
        for( int64_t j = 0; j < k; j++ ) {
            H1(i, j) = 0.0;
            H2(i, j) = 0.0;
        }
        H1(i, k) = 0.0;
    }

    // Start search for x interval from index k (first possible non-zero region)
    int64_t ind = k;

    // Loop through each x value to construct rows of A, H1, H2
    for (int i=0; i < m; ++i) {
        double xval = x(i);

        // Find index `ind` such that t[ind] <= xval < t[ind + 1]
        ind = _find_interval(t.data, len_t, k, xval, ind, extrapolate);

        // Handle invalid case: x out of knot range and extrapolation not allowed
        if (!extrapolate && ind < 0) {
            throw std::runtime_error("find_interval: out of bounds with x = " + std::to_string(xval));
        }

        // Set offset for this row; it indicates the starting index in the full coefficient array
        // where the k+1 non-zero basis functions for x[i] should be applied
        int64_t offseti = ind - k;
        offset(i) = offseti;

        // Compute the k+1 non-zero B-spline basis values at xval
        _deBoor_D(t.data, xval, k, ind, 0, wrk);

        // Fill A matrix with weighted basis values
        for (int64_t j = 0; j < k + 1; ++j) {
            A(i, j) = wrk[j] * w(i);
        }

        // ---- PERIODIC CORRECTION HANDLING ----
        // ----------------------------------------------------------------------------
        // In periodic B-splines, basis functions near the domain boundaries "wrap around"
        // due to periodicity. This causes the data matrix to be conceptually banded
        // but with values that go past the left or right ends.
        //
        // To handle this, we split the full dense data matrix D row-wise into:
        //   - A(i, :)  → in-range basis values (no wrap)
        //   - H1(i, :) → wrapped values from right end that wrap back into domain (positive wrap)
        //   - H2(i, :) → wrapped values from left underflow (negative wrap)
        //
        // This preserves sparsity and enables efficient periodic QR reduction later.
        //
        // Example (assume degree k = 2):
        //
        // Full dense row D[i]: [ 0   0   b0   b1   b2   0   0   0   b3   b4 ]
        // Index (j):           [ 0   1    2    3    4   5   6   7    8    9 ]
        //                                      ↑    ↑    ↑
        //                                      A    A    A
        //                                                        ↑    ↑
        //                                                       H2   H2
        //
        // Suppose offset[i] = 2, so:
        //
        //   A[i]  = [ b0, b1, b2 ]       → basis at indices 2, 3, 4
        //   H1[i] = [  0,  0,  0 ]       → no positive wrap-around here
        //   H2[i] = [ b4, b3 ]           → values that wrapped around left
        //                                 and landed in indices -1, 0
        //
        // This decomposition allows us to keep the data matrix compact,
        // and still capture all contributions due to periodicity.
        // ----------------------------------------------------------------------------


        // If x lies in the last 2k intervals of the knot vector, apply correction.
        int64_t l = ind + 1;
        if (l >= len_t - 2*k) {
            int64_t j;
            // First, reset H1 and H2 for the current row to zero
            for (int64_t j = 0; j < k; ++j) {
                H1(i, j) = 0;
                H2(i, j) = 0;
            }
            H1(i, k) = 0;

            // Initialize j to the index offset for periodic correction
            j = l - len_t + 2*k;

            // Loop over all (k+1) non-zero B-spline basis functions
            for (int64_t i1 = 0; i1 < k+1; i1++) {
                j += 1;

                // l0: index after wrapping right end of the knot vector (used for periodic correction)
                int64_t l0 = j;

                // l1: corresponding "wrapped" index in original coefficient space
                int64_t l1 = l0 - k;

                // If l1 is still out of valid range, wrap it again until it fits
                while (l1 > std::max(int64_t(0), len_t - 3*k - 1)) {
                    // Apply modular arithmetic to wrap l1 to the start
                    l0 = l1 - len_t + 3*k + 1;
                    l1 = l0 - k;
                }

                // Now assign weighted basis value to H1 or H2 depending on index
                if( l1 > 0 ) {
                    // Normal wrap-around: put contribution into H1
                    H1(i, l1 - 1) = wrk[i1] * w(i);
                } else {
                    // Negative wrap-around: accumulate contribution into H2
                    // This ensures full periodic coverage even for overlapping parts
                    H2(i, l0 - 1) += wrk[i1] * w(i);
                }
            }
        }
    }

    // Set total number of coefficients = number of intervals = len(t) - k - 1
    *nc = len_t - k - 1;
}


/*
 *   Solve the LSQ problem ||y - A@c||^2 via QR factorization.
 *
    QR factorization follows FITPACK: we reduce A row-by-row by Givens rotations.
    To zero out the lower triangle, we use in the row `i` and column `j < i`,
    the diagonal element in that column. That way, the sequence is
    (here `[x]` are the pair of elements to Givens-rotate)

     [x] x x x       x  x  x x      x  x  x x      x x  x  x      x x x x
     [x] x x x  ->   0 [x] x x  ->  0 [x] x x  ->  0 x  x  x  ->  0 x x x
      0  x x x       0 [x] x x      0  0  x x      0 0 [x] x      0 0 x x
      0  x x x       0  x  x x      0 [x] x x      0 0 [x] x      0 0 0 x

    The matrix A has a special structure: each row has at most (k+1) non-zeros, so
    is encoded as a PackedMatrix instance.

    On exit, the return matrix, also of shape (m, k+1), contains
    elements of the upper triangular matrix `R[i, i: i + k + 1]`.
    When we process the element (i, j), we store the rotated row in R[i, :],
    and *shift it to the left*, so that the the diagonal element is always in the
    zero-th place. This way, the process above becomes

     [x] x x x       x  x x x       x  x x x       x  x x x      x x x x
     [x] x x x  ->  [x] x x -  ->  [x] x x -   ->  x  x x -  ->  x x x -
      x  x x -      [x] x x -       x  x - -      [x] x - -      x x - -
      x  x x -       x  x x -      [x] x x -      [x] x - -      x - - -

    The most confusing part is that when rotating the row `i` with a row `j`
    above it, the offsets differ: for the upper row  `j`, `R[j, 0]` is the diagonal
    element, while for the row `i`, `R[i, 0]` is the element being annihilated.

    NB. This row-by-row Givens reduction process follows FITPACK:
    https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L112-L161
    A possibly more efficient way could be to note that all data points which
    lie between two knots all have the same offset: if `t[i] < x_1 .... x_s < t[i+1]`,
    the `s-1` corresponding rows form an `(s-1, k+1)`-sized "block".
    Then a blocked QR implementation could look like
    https://people.sc.fsu.edu/~jburkardt/f77_src/band_qr/band_qr.f

    The `startrow` optional argument accounts for the scenatio with a two-step
    factorization. Namely, the preceding rows are assumend to be already
    processed and are skipped.
    This is to account for the scenario where we append new rows to an already
    triangularized matrix.

    This routine MODIFIES `a` & `y` in-place.
 */
void
qr_reduce(double *aptr, const int64_t m, const int64_t nz, // a(m, nz), packed
          int64_t *offset,                                 // offset(m)
          const int64_t nc,                                // dense would be a(m, nc)
          double *yptr, const int64_t ydim1,               // y(m, ydim2)
          const int64_t startrow
)
{
    auto R = RealArray2D(aptr, m, nz);
    auto y = RealArray2D(yptr, m, ydim1);

    for (int64_t i=startrow; i < m; ++i) {
        int64_t oi = offset[i];
        for (int64_t j=oi; j < nc; ++j) {

            // rotate only the lower diagonal
            if (j >= std::min(i, nc)) {
                break;
            }

            // in dense format: diag a1[j, j] vs a1[i, j]
            double c, s, r;
            DLARTG(&R(j, 0), &R(i, 0), &c, &s, &r);

            // rotate l.h.s.
            R(j, 0) = r;
            for (int64_t l=1; l < R.ncols; ++l) {
                std::tie(R(j, l), R(i, l-1)) = fprota(c, s, R(j, l), R(i, l));
            }
            R(i, R.ncols-1) = 0.0;

            // rotate r.h.s.
            for (int64_t l=0; l < y.ncols; ++l) {
                std::tie(y(j, l), y(i, l)) = fprota(c, s, y(j, l), y(i, l));
            }
        }
        if (i < nc) {
            offset[i] = i;
        }

    } // for(i = ...
}

/**
 * Performs QR reduction for periodic B-spline fitting.
 *
 * aptr: pointer to A (m, nz), holds basis function values
 * h1ptr: periodic correction matrix H1, shape (m, k+1).
 * h2ptr: periodic correction matrix H2, shape (m, k).
 * m: number of data points.
 * nz: number of non-zero B-spline basis functions.
 * offset: index array of shape (m,), each x lies in knot interval [offset(i) + k, offset(i) + k + 1).
 * nc: number of spline coefficients (len_t - k - 1).
 * yptr: pointer to y (m + 1, ydim1): observed y-values (or RHS).
 * ydim1: number of columns in y (typically 1).
 * k: spline degree.
 * len_t: length of knot vector.
 * a1ptr: pointer to A1 ((len_t - k - 1) x (k + 1)): upper-triangular QR matrix.
 * a2ptr: pointer to A2 ((len_t - 2k - 1) x k): for periodic boundary reduction.
 * zptr: pointer to z (len_t - k - 1): transformed RHS vector.
 * init_p: if true, compute smoothing term 'p'.
 * p: Reference to smoothing term.
 */
void qr_reduce_periodic(
    double *aptr, double *h1ptr, double *h2ptr,
    const int64_t m, const int64_t nz,
    int64_t *offset,
    const int64_t nc,
    double *yptr, const int64_t ydim1,
    const int k, const int64_t len_t,
    double *a1ptr,
    double *a2ptr,
    double *zptr,
    double& fp,
    bool init_p,
    double& p
) {
    // Wrap raw pointers with helper matrix/array types
    auto H  = RealArray2D(aptr, m, nz);
    auto H1 = RealArray2D(h1ptr, m, nz);
    auto H2 = RealArray2D(h2ptr, m, nz - 1);
    auto A1 = RealArray2D(a1ptr, len_t - k - 1, k + 1);      // Main upper-triangular matrix
    auto A2 = RealArray2D(a2ptr, len_t - 2*k - 1, k);        // For enforcing periodicity
    auto z  = RealArray2D(zptr, len_t - k - 1, ydim1);              // Transformed RHS
    auto y  = RealArray2D(yptr, m + 1, ydim1);               // Input/output Y vector
    fp = 0.0;

    std::vector<double> yi(ydim1, 0.0);

    // Reset A1 and z to zero
    for( int64_t i = 1; i <= len_t - k - 1; i++ ) {
        for( int64_t ydimi = 0; ydimi < ydim1; ydimi++ ) {
            z(i - 1, ydimi) = 0.0;
        }
        for( int64_t j = 1; j <= k + 1; j++ ) {
            A1(i - 1, j - 1) = 0.0;
        }
    }

    int64_t jper = 0;                    // Flag to initialize periodic block once
    int64_t nk1 = len_t - k - 1;         // Total number of B-spline coefficients
    int64_t n7 = nk1 - k;                // coefficients not used for periodicity
    int64_t n10 = n7 - k;                // coefficients strictly away from periodic boundary

    // Main loop over all data points
    for( int64_t it = 1; it <= m; it++ ) {
        for( int64_t ydimi = 0; ydimi < ydim1; ydimi++ ) {
            yi[ydimi] = y(it - 1, ydimi);
        }
        int64_t ind = offset[it - 1] + k;
        int64_t l = ind + 1;                          // l = leftmost non-zero + degree + 1
        int64_t l5 = l - k - 1;

        // Case 1: Point falls within the periodic boundary range
        if (l5 >= n10) {

            // Initialize A2 for the periodic case — run only once
            if( jper == 0 ) {
                for( int64_t i = 1; i <= n7; i++ ) {
                    for( int64_t j = 1; j <= k; j++ ) {
                        A2(i - 1, j - 1) = 0;
                    }
                }

                // Copy tail of A1 into top of A2 for enforcing periodic constraints
                int64_t jk = n10 + 1;
                for( int64_t i = 1; i <= k; i++ ) {
                    int64_t ik = jk;
                    for( int64_t j = 1; j <= k + 1; j++ ) {
                        if( ik <= 0 ) break;
                        A2(ik - 1, i - 1) = A1(ik - 1, j - 1);
                        ik--;
                    }
                    jk++;
                }
                jper = 1;
            }

            // Phase 1: Apply Givens rotations to reduce A1 and A2 to triangular form
            if (n10 > 0) {
                for( int64_t j = 1; j <= n10; j++ ) {
                    double piv = H1(it - 1, 0);
                    if( piv == 0.0 ) {
                        // Shift H1 left
                        for (int64_t h1i = 1; h1i <= k; h1i++) {
                            H1(it - 1, h1i - 1) = H1(it - 1, h1i);
                        }
                        H1(it - 1, k) = 0.0;
                    } else {
                        // Apply Givens rotation
                        double c, s, r;
                        DLARTG(&A1(j - 1, 0), &piv, &c, &s, &r);
                        A1(j - 1, 0) = r;

                        for( int64_t ydimi = 0; ydimi < ydim1; ydimi++ ) {
                            std::tie(z(j - 1, ydimi), yi[ydimi]) =
                                fprota(c, s, z(j - 1, ydimi), yi[ydimi]);
                        }
                        for( int64_t h2i = 1; h2i <= k; h2i++ ) {
                            std::tie(A2(j - 1, h2i - 1), H2(it - 1, h2i - 1)) = fprota(
                                c, s, A2(j - 1, h2i - 1), H2(it - 1, h2i - 1));
                        }

                        if( j == n10 ) {
                            break;
                        }

                        int64_t i2 = std::min(n10 - j, (int64_t)k) + 1;
                        for( int64_t h1i = 1; h1i < i2; h1i++ ) {
                            std::tie(A1(j - 1, h1i), H1(it - 1, h1i)) = fprota(c, s, A1(j - 1, h1i), H1(it - 1, h1i));
                        }

                        for( int64_t h1i = 1; h1i <= i2; h1i++ ) {
                            H1(it - 1, h1i - 1) = H1(it - 1, h1i);
                        }
                        H1(it - 1, i2 - 1) = 0.0;
                    }
                }
            }

            // Phase 2: Givens rotations on A2 portion
            for( int64_t j = 1; j <= k; j++ ) {
                int64_t ij = n10 + j;
                double piv = H2(it - 1, j - 1);
                if (ij <= 0 || piv == 0.0) {
                    continue;
                }

                double c, s, r;
                DLARTG(&A2(ij - 1, j - 1), &piv, &c, &s, &r);
                A2(ij - 1, j - 1) = r;
                for( int64_t ydimi = 0; ydimi < ydim1; ydimi++ ) {
                    std::tie(z(ij - 1, ydimi), yi[ydimi]) =
                        fprota(c, s, z(ij - 1, ydimi), yi[ydimi]);
                }

                if( j == k ) {
                    break;
                }

                for( int64_t h2i = j + 1; h2i <= k; h2i++ ) {
                    std::tie(A2(ij - 1, h2i - 1), H2(it - 1, h2i - 1)) = fprota(
                        c, s, A2(ij - 1, h2i - 1), H2(it - 1, h2i - 1));
                }
            }

        } else {
            // Case 2: Observation lies completely outside periodic range
            int64_t j = l5;
            for( int64_t i = 1; i <= k + 1; i++ ) {
                j++;
                double piv = H(it - 1, i - 1);
                if( piv == 0.0 ) {
                    continue;
                }

                double c, s, r;
                DLARTG(&A1(j - 1, 0), &piv, &c, &s, &r);
                A1(j - 1, 0) = r;
                for( int64_t ydimi = 0; ydimi < ydim1; ydimi++ ) {
                    std::tie(z(j - 1, ydimi), yi[ydimi]) = fprota(c, s, z(j - 1, ydimi), yi[ydimi]);
                }

                if( i == k + 1 ) {
                    break;
                }
                int64_t i2 = 1;
                for( int64_t i1 = i + 1; i1 <= k + 1; i1++ ) {
                    i2++;
                    std::tie(A1(j - 1, i2 - 1), H(it - 1, i1 - 1)) = fprota(
                        c, s, A1(j - 1, i2 - 1), H(it - 1, i1 - 1));
                }
            }
        }
        for( int64_t ydimi = 0; ydimi < ydim1; ydimi++ ) {
            fp += yi[ydimi] * yi[ydimi];
        }
    }

    // Compute smoothing term 'p'
    if( init_p ) {
        p = 0.0;
        int64_t l = n7;
        for( int64_t i = 1; i <= k; i++ ) {
            int64_t j = k + 1 - i;
            p += A2(l - 1, j - 1);
            l--;
            if( l == 0 ) {
                break;
            }
        }

        if( l != 0 ) {
            for( int64_t i = 0; i < n10; i++ ) {
                p += A1(i, 0);
            }
        }
        p = n7/p;
    }
}

/**
 * Applies QR reduction on augmented matrices using Givens rotations.
 *
 * Example for k=1, len_t=7 (small sizes):
 *
 * Dimensions:
 *   g1: (4 x 3), g2: (4 x 2)
 *   h1: (3 x 3), h2: (3 x 2)
 *   c : (5 x 1), offset: (3)
 *
 * Initial matrices (example values):
 * g1 = [ [3.0, 2.0, 1.0],
 *        [4.0, 3.0, 2.0],
 *        [5.0, 4.0, 3.0],
 *        [6.0, 5.0, 4.0] ]
 *
 * g2 = [ [1.0, 2.0],
 *        [2.0, 3.0],
 *        [3.0, 4.0],
 *        [4.0, 5.0] ]
 *
 * h1 = [ [2.0, 1.0, 0.0],
 *        [1.0, 2.0, 1.0],
 *        [0.0, 1.0, 2.0] ]
 *
 * h2 = [ [1.0, 0.0],
 *        [0.0, 1.0],
 *        [1.0, 1.0] ]
 *
 * c = [ [1.0],
 *       [1.5],
 *       [2.0],
 *       [2.5],
 *       [3.0] ]
 *
 * offset = [1, 1, 2]
 *
 * Iteration it=1:
 * - offset(0) = 1
 * - For j=1 to 2:
 *    * piv = h1(0,0) = 2.0
 *    * Apply Givens rotation on (g1(0,0)=3.0, piv=2.0) to zero out h1(0,0):
 *       - Compute (cos, sin, r) so that:
 *         [cos sin; -sin cos] * [3.0; 2.0] = [r; 0]
 *       - Example result: r=3.6056, cos=0.8321, sin=0.5547
 *    * Update g1(0,0) = r = 3.6056
 *    * Rotate c(0,0) and yi=0 with (cos,sin), c(0,0) ≈ 0.8321, yi ≈ -0.5547
 *    * For h2i in 0 to 1:
 *      - Rotate pairs (g2(0,h2i), h2(0,h2i)) with same rotation
 *    * For h1i in 1 to min(k+1, remaining cols):
 *      - Rotate pairs (g1(0,h1i), h1(0,h1i))
 *    * Shift h1 row left and zero last element
 *
 * Result:
 *   h1(0,0) is zeroed,
 *   corresponding elements in g1, g2, h2, c updated to maintain orthogonality.
 *
 * Subsequent iterations it=2, it=3 do similar Givens rotations to reduce
 * matrices further towards upper-triangular form for stable solving.
 *
 * Purpose:
 * - Systematically eliminate entries in h1, h2 (constraint matrices)
 * - Update g1, g2 (main matrices) and c (solution vector)
 * - Prepare system for efficient and stable spline coefficient computation
 */
void qr_reduce_augmented_matrices(
    double* g1ptr, double* g2ptr,
    double* h1ptr, double* h2ptr,
    double* cptr, double* offsetptr,
    int k, int64_t len_t, int64_t ydim2
) {
    auto g1 = RealArray2D(g1ptr, len_t - 2*k - 1, k + 2);  // Main block matrix part 1
    auto g2 = RealArray2D(g2ptr, len_t - 2*k - 1, k + 1);  // Main block matrix part 2
    auto h1 = RealArray2D(h1ptr, len_t - 2*k - 2, k + 2);  // Constraint matrix part 1 (to be zeroed out)
    auto h2 = RealArray2D(h2ptr, len_t - 2*k - 2, k + 1);  // Constraint matrix part 2 (to be zeroed out)
    auto c = RealArray2D(cptr, len_t - k - 1, ydim2);      // Right-hand side vector (updated during rotations)
    auto offset = RealArray1D(offsetptr, len_t - 2*k - 2); // Offset array, starting row for rotations

    std::vector<double> yi(ydim2, 0.0);

    // Iterate over rows in the augmented matrix system (h1, h2 rows)
    for( int64_t it = 1; it <= len_t - 2*k - 2; it++ ) {
        for( int64_t ydimi = 0; ydimi < ydim2; ydimi++ ) {
            yi[ydimi] = 0.0; // Auxiliary variable to hold rotated value of c
        }

        // Perform Givens rotations to zero out entries in h1 and update g1, g2, h2, c accordingly
        for( int64_t j = offset(it - 1); j <= len_t - 3*k - 2; j++ ) {
            // Pivot element from h1 that we want to zero out using rotation
            double piv = h1(it - 1, 0);

            double cos, sin, r;

            // Compute Givens rotation parameters (cos, sin) to zero out piv when applied to g1(j-1,0)
            // This solves: [cos sin; -sin cos] * [g1(j-1,0); piv] = [r; 0]
            DLARTG(&g1(j - 1, 0), &piv, &cos, &sin, &r);

            // Update g1 with the rotated value r (magnitude after rotation)
            g1(j - 1, 0) = r;

            // Apply the same rotation to c(j-1, 0) and yi (c vector updated to maintain system consistency)
            for( int64_t ydimi = 0; ydimi < ydim2; ydimi++ ) {
                std::tie(c(j - 1, ydimi), yi[ydimi]) =
                    fprota(cos, sin, c(j - 1, ydimi), yi[ydimi]); // Auxiliary variable to hold rotated value of c
            }


            // Rotate pairs of elements in g2 and h2 with the same rotation to maintain orthogonality
            for( int64_t h2i = 0; h2i < k + 1; h2i++ ) {
                std::tie(g2(j - 1, h2i), h2(it - 1, h2i)) = fprota(cos, sin, g2(j - 1, h2i), h2(it - 1, h2i));
            }

            // If reached the last j in the loop, break since no more elements to rotate
            if( j == (len_t - 3*k - 2) ) {
                break ;
            }

            // Calculate number of elements to rotate between g1 and h1 for this step
            int64_t i2 = std::min(len_t - 3*k - 2 - j, (int64_t) k + 1) + 1;

            // Rotate pairs of elements in g1 and h1 with the same rotation for indices > 0
            for( int64_t h1i = 1; h1i < i2; h1i++ ) {
                std::tie(g1(j - 1, h1i), h1(it - 1, h1i)) = fprota(cos, sin, g1(j - 1, h1i), h1(it - 1, h1i));
            }

            // Shift h1 row left by one position to remove the zeroed element and set last to zero
            for( int64_t h1i = 1; h1i <= (i2 - 1); h1i++ ) {
                h1(it - 1, h1i - 1) = h1(it - 1, h1i);
            }
            h1(it - 1, i2 - 1) = 0.0; // Clear last shifted position (zeroed out)
        }

        // Now handle rotations on h2 part to zero out entries similarly
        for( int64_t j = 1; j <= k + 1; j++ ) {
            int64_t ij = len_t - 3*k - 2 + j; // Index in g2 to rotate with h2

            if( ij <= 0 ) {
                continue; // Skip if index out of bounds (safety check)
            }

            // Pivot element from h2 to be zeroed
            double piv = h2(it - 1, j - 1);

            double cos, sin, r;

            // Compute Givens rotation parameters to zero out pivot against g2(ij-1, j-1)
            DLARTG(&g2(ij - 1, j - 1), &piv, &cos, &sin, &r);

            // Update g2 with rotated magnitude
            g2(ij - 1, j - 1) = r;

            // Rotate c(ij-1, 0) and yi with the same rotation to keep system consistent
            for( int64_t ydimi = 0; ydimi < ydim2; ydimi++ ) {
                std::tie(c(ij - 1, ydimi), yi[ydimi]) =
                    fprota(cos, sin, c(ij - 1, ydimi), yi[ydimi]);
            }

            // If this is the last column in h2 for this row, break early
            if( j == k + 1) {
                break;
            }

            // Rotate remaining pairs in g2 and h2 in the same row with the rotation
            for( int64_t h2i = j + 1; h2i <= k + 1; h2i++ ) {
                std::tie(g2(ij - 1, h2i - 1), h2(it - 1, h2i - 1)) = fprota(
                    cos, sin, g2(ij - 1, h2i - 1), h2(it - 1, h2i - 1));
            }
        }
    }
}

void _compute_residuals(
    /* inputs */
    const double *xptr, int64_t m,      // x, shape (m,)
    const double *yptr, int64_t ydim2,            // y(m, ydim2)
    const double *tptr, int64_t len_t,  // t, shape (len_t,)
    const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
    int k,
    int extrapolate,
    int64_t nc,
    const double *cptr,                                 // c(nc, ydim2)
    /* output */
    double *fp,
    double *residualsptr                            // residuals(m,)
) {
    auto x = ConstRealArray1D(xptr, m);
    auto y = ConstRealArray2D(yptr, m, ydim2);
    auto t = ConstRealArray1D(tptr, len_t);
    auto w = ConstRealArray1D(wptr, m);
    auto residuals = RealArray1D(residualsptr, m);
    auto c = ConstRealArray2D(cptr, nc, ydim2);

    *fp = 0.0;
    double l = k + 2;
    int64_t ind = k;
    std::vector<double> wrk(2*k + 2);
    for( int64_t it = 1; it <= m; it++ ) {
        if( !(x(it - 1) < t(l - 1) || l > len_t - k - 1) ) {
            l = l + 1;
        }

        double xval = x(it - 1);

        // find the interval
        ind = _find_interval(t.data, len_t, k, xval, ind, extrapolate);

        // compute non-zero b-splines
        _deBoor_D(t.data, xval, k, ind, 0, wrk.data());

        residuals(it - 1) = 0.0;
        for( int64_t yi = 0; yi < ydim2; yi++ ) {
            double term = 0.0;
            int64_t l0 = l - k - 2;

            for( int64_t j = 1; j <= k + 1; j++ ) {
                l0 = l0 + 1;
                term = term + c(l0 - 1, yi) * wrk[j - 1];
            }
            double delta = w(it - 1) * (term - y(it - 1, yi));
            term = delta * delta;
            residuals(it - 1) += term;
            *fp = *fp + term;
        }
    }
}


/*
 * Back substitution solve of `R @ c = y` with an upper triangular R.
 *
 * R is in the 'packed' format:
 *    1. Each row has at most `nz` non-zero elements, stored in a (m, nz) array
 *    2. All non-zeros are consecutive.
 * Since R is upper triangular, non-zeros in row `i` start at `i`. IOW,
 * the 'offset' array is `np.arange(nc)`.
 *
 * If `R` were dense, it would have had shape `(m, nc)`. Since R is upper triangular,
 * the first `nc` rows contain non-zeros; the last `m-nc` rows are all zeros
 * (in fact, they may contain whatever, and are not referenced in the routine).
 *
 * `y` array is always two-dimensional, and has shape `(m, ydim2)`.
 * IOW, if the original data is 1D, `y.shape == (m, 1)`.
 *
 * The output `c` array has the shape `(nc, ydim2)`.
 */
void
fpback( /* inputs*/
       const double *Rptr, int64_t m, int64_t nz,    // R(m, nz), packed
       int64_t nc,                                   // dense R would be (m, nc)
       const double *xptr, int64_t m_,      // x, shape (m,)
       const double *tptr, int64_t len_t,  // t, shape (len_t,)
       int k,
       const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
       int extrapolate,
       const double* ywptr,
       const double *yptr, int64_t ydim2,            // y(m, ydim2)
        /* output */
       double *cptr,                                 // c(nc, ydim2)
       double *fp,
       double *residualsptr)                            // residuals(m,)
{
    auto R = ConstRealArray2D(Rptr, m, nz);
    auto yw = ConstRealArray2D(ywptr, m, ydim2);
    auto c = RealArray2D(cptr, nc, ydim2);

    // c[nc-1, ...] = y[nc-1] / R[nc-1, 0]
    for (int64_t l=0; l < ydim2; ++l) {
        c(nc - 1, l) = yw(nc - 1, l) / R(nc - 1, 0);
    }

    //for i in range(nc-2, -1, -1):
    //    nel = min(nz, nc-i)
    //    c[i, ...] = ( y[i] - (R[i, 1:nel, None] * c[i+1:i+nel, ...]).sum(axis=0) ) / R[i, 0]
    for (int64_t i=nc-2; i >= 0; --i) {
        int64_t nel = std::min(nz, nc - i);
        for (int64_t l=0; l < ydim2; ++l){
            double ssum = yw(i, l);
            for (int64_t j=1; j < nel; ++j) {
                ssum -= R(i, j) * c(i + j, l);
            }
            ssum /= R(i, 0);
            c(i, l) = ssum;
        }
    }

    _compute_residuals(
        xptr, m_, yptr, ydim2, tptr, len_t,
        wptr, k, extrapolate, nc, cptr, fp, residualsptr
    );
}

void _fpbacp(
        ConstRealArray2D& A1,
        ConstRealArray2D& A2,
        ConstRealArray2D& Z,
        int k, int kp, int64_t len_t,
        RealArray2D& c) {

    int64_t nc = len_t - k - 1;
    int64_t n = nc - k;
    int64_t n2 = n - kp;
    int64_t l = n;
    std::vector<double> store(Z.ncols, 0.0);
    for( int64_t i = 1; i <= kp; i++ ) {
        for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
            store[ydimi] = Z(l - 1, ydimi);
        }
        int64_t j = kp + 2 - i;
        if( i != 1 ) {
            int64_t l0 = l;
            for( int64_t l1 = j; l1 <= kp; l1++ ) {
                l0 = l0 + 1;
                for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
                    store[ydimi] = store[ydimi] - c(l0 - 1, ydimi) * A2(l - 1, l1 - 1);
                }
            }
        }
        for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
            c(l - 1, ydimi) = store[ydimi]/A2(l - 1, j - 2);
        }
        l = l - 1;
        if( l == 0 ) {
            return ;
        }
    }
    for( int64_t i = 1; i <= n2; i++ ) {
        for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
            store[ydimi] = Z(i - 1, ydimi);
        }
        l = n2;
        for( int64_t j = 1; j <= kp; j++ ) {
            l = l + 1;
            for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
                store[ydimi] = store[ydimi] - c(l - 1, ydimi) * A2(i - 1, j - 1);
            }
        }
        for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
            c(i - 1, ydimi) = store[ydimi];
        }
    }
    int64_t i = n2;
    for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
        c(i - 1, ydimi) = c(i - 1, ydimi)/A1(i - 1, 0);
    }
    if( i == 1 ) {
        return ;
    }
    for( int64_t j = 2; j <= n2; j++ ) {
        i = i - 1;
        for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
            store[ydimi] = c(i - 1, ydimi);
        }
        int64_t i1 = kp;
        if( j <= kp ) {
            i1 = j - 1;
        }
        l = i;
        for( int64_t l0 = 1; l0 <= i1; l0++ ) {
            l = l + 1;
            for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
                store[ydimi] = store[ydimi] - c(l - 1, ydimi) * A1(i - 1, l0);
            }
        }
        for( int64_t ydimi = 0; ydimi < Z.ncols; ydimi++ ) {
            c(i - 1, ydimi) = store[ydimi]/A1(i - 1, 0);
        }
    }
}

void
fpbacp( /* inputs*/
       const double *A1ptr,
       int64_t a1_rows,
       const double *A2ptr,
       int64_t a2_rows,
       const double *Zptr,
       int k, int kp,
       const double *xptr, int64_t m,      // x, shape (m,)
       const double *yptr, int64_t ydim2,            // y(m, ydim2)
       const double *tptr, int64_t len_t,  // t, shape (len_t,)
       const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
       int extrapolate,
       /* output */
       double *cptr,
       double *fp,
       double *residualsptr) {

    int64_t nc = len_t - k - 1;
    auto A1 = ConstRealArray2D(A1ptr, a1_rows, kp + 1);
    auto A2 = ConstRealArray2D(A2ptr, a2_rows, kp);
    auto Z = ConstRealArray2D(Zptr, nc, ydim2);
    auto c = RealArray2D(cptr, nc, ydim2);

    _fpbacp(A1, A2, Z, k, kp, len_t, c);

    int64_t offset = len_t - 2*k - 1;
    for( int64_t i = 0; i < k; i++ ) {
        for( int64_t ydimi = 0; ydimi < ydim2; ydimi++ ) {
            c(i + offset, ydimi) = c(i, ydimi);
        }
    }

    _compute_residuals(
        xptr, m, yptr, ydim2, tptr, len_t,
        wptr, k, extrapolate, nc, cptr, fp, residualsptr
    );

}


/*
 * A helper for _fpknot:
 *
 * Split the `x` array into knot "runs" and sum the residuals per "run".
 *
 * Here a "run" is a set of `x` values which lie between consecutive knots:
 * these are `x(i)` which for a given `j` satisfy `t(j+k) <= x(i) <= t(j+k+1)`.
 *
 * The _split routine returns two vectors: a vector of indices into `x`, `ix`,
 * and a vector of partial sums of residuals, `fparts`.
 * The i-th entry is the i-th "run". IOW the pair (fparts, ix) means that
 * `fparts[i]` is the sum of residuals over the "run" x[ix[i]] <= xvalue <= x[ix[i+1]].
 *
 * This routine is a (best-effort) translation of
 * https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpcurf.f#L190-L215
 *
 */
pair_t
_split(ConstRealArray1D x, ConstRealArray1D t, int k, ConstRealArray1D residuals)
{
    /*
     * c  search for knot interval t(number+k) <= x <= t(number+k+1) where
     * c  fpint(number) is maximal on the condition that nrdata(number)
     * c  not equals zero.
     */
    int64_t interval = k+1;
    int64_t nc = t.nelem - k - 1;

    std::vector<int64_t> ix;
    ix.push_back(0);

    std::vector<double> fparts;
    double fpart = 0.0;

    for(int64_t i=0; i < x.nelem; i++) {
        double xv = x(i);
        double rv = residuals(i);
        fpart += rv;

        if ((xv >= t(interval)) && (interval < nc)) {
            // end of the current interval: split the weight at xv by 1/2
            // between two intervals
            double carry = rv / 2.0;
            fpart -= carry;
            fparts.push_back(fpart);

            fpart = carry;
            interval++;

            ix.push_back(i);
        }
    } // for i

    // the last interval
    ix.push_back(x.nelem - 1);
    fparts.push_back(fpart);

    return std::make_tuple(fparts, ix);
}


/*
 * Find a position for a new knot, a la FITPACK.
 *
 *   (Approximately) replicate FITPACK's logic:
 *     1. split the `x` array into knot intervals, ``t(j+k) <= x(i) <= t(j+k+1)``
 *     2. find the interval with the maximum sum of residuals
 *     3. insert a new knot into the middle of that interval.
 *
 *   NB: a new knot is in fact an `x` value at the middle of the interval.
 *   So *the knots are a subset of `x`*.
 *
 *   This routine is an analog of
 *   https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpcurf.f#L190-L215
 *   (cf _split function)
 *
 *   and https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpknot.f
 */
double
fpknot(const double *x_ptr, int64_t m,
       const double *t_ptr, int64_t len_t,
       int k,
       const double *residuals_ptr)
{
    auto x = ConstRealArray1D(x_ptr, m);
    auto t = ConstRealArray1D(t_ptr, len_t);
    auto residuals = ConstRealArray1D(residuals_ptr, m);

    std::vector<double> fparts;
    std::vector<int64_t> ix;
    std::tie(fparts, ix) = _split(x, t, k, residuals);

    int64_t idx_max = -101;
    double fpart_max = -1e100;
    for (size_t i=0; i < fparts.size(); i++) {
        bool is_better = (ix[i+1] - ix[i] > 1) && (fparts[i] > fpart_max);
        if(is_better) {
            idx_max = static_cast<int64_t>(i);
            fpart_max = fparts[i];
        }
    }

    if (idx_max == -101) {
        throw std::runtime_error("Internal error. Please report it to SciPy developers.");
    }

    // round up, like Dierckx does? This is really arbitrary though.
    int64_t idx_newknot = (ix[idx_max] + ix[idx_max+1] + 1) / 2;

    return x(idx_newknot);
}


/*
 * Evaluate the spline function
*/
void
_evaluate_spline(
    const double *tptr, int64_t len_t,         // t, shape (len_t,)
    const double *cptr, int64_t n, int64_t m,  // c, shape (n, m)
    int64_t k,
    const double *xp_ptr, int64_t s,           // xp, shape (s,)
    int64_t nu,
    int extrapolate,
    double *out_ptr,                           // out, shape (s, m) NOT CHECKED
    double *wrk                                // scratch, shape (2k+2,)
)
{
    auto t = ConstRealArray1D(tptr, len_t);
    auto c = ConstRealArray2D(cptr, n, m);
    auto xp = ConstRealArray1D(xp_ptr, s);
    auto out = RealArray2D(out_ptr, s, m);

    int64_t interval = k;
    for(int64_t ip=0; ip < s; ip++) {

        double xval = xp(ip);
        interval = _find_interval(t.data, len_t, k, xval, interval, extrapolate);
        if (interval < 0) {
            // xval was nan etc
            for (int64_t jp=0; jp < m; jp++) {
                out(ip, jp) = std::numeric_limits<double>::quiet_NaN();
            }
        }
        else {
            // Evaluate (k+1) b-splines which are non-zero on the interval.
            // on return, first k+1 elements of work are B_{m-k},..., B_{m}
            _deBoor_D(t.data, xval, k, interval, nu, wrk);

            // Form linear combinations
            for (int64_t jp=0; jp < m; jp++) {
                out(ip, jp) = 0.0;
                for (int64_t a=0; a < k+1 ; a++){
                    out(ip, jp) += c(interval + a -k, jp) * wrk[a];
                }
            }
        }
    }
}


/*
 * Spline colocation matrix in the LAPACK banded storage
 */
void
_coloc_matrix(const double *xptr, int64_t m,       // x, shape(m,)
              const double *tptr, int64_t len_t,   // t, shape(len_t,)
              int k,
              double *abT_ptr, int64_t nbands,     // ab(nbands, len_t - k - 1) in F order!
              int offset,
              double *wrk                          // scratch, shape (2k+2)
)
{
    auto x = ConstRealArray1D(xptr, m);
    auto t = ConstRealArray1D(tptr, len_t);
    auto abT = RealArray2D(abT_ptr, len_t - k - 1, nbands); // NB: transposed in C order

    int64_t kl = k, ku = k;  // upper and lower bands; NB: nbands == 2*ku+kl+1, not checked
    int64_t left = k;
    for(int64_t j=0; j < m; j++) {
        double xval = x(j);
        left = _find_interval(t.data, len_t, k, xval, left, 0);

        // Evaluate (k+1) b-splines which are non-zero on the interval `left`.
        // on return, first k+1 elements of work are B_{m-k},..., B_{m}
        _deBoor_D(t.data, xval, k, left, 0, wrk);

        // Fill a row. For a full matrix in the C order, it would've been
        // ``A[j+offset, left-k:left+1] = wrk``.
        // In the LAPACK banded storage, need to spread the row over:
        // https://www.netlib.org/lapack/lug/node124.html
        // Additionally, for the Fortran order, we operate on the transposed matrix
        // (by just swapping the row and column indices)
        for (int64_t a=0; a < k+1 ; a++) {
            int64_t clmn = left - k + a;
            abT(clmn, kl + ku + j + offset - clmn) = wrk[a];
        }
    }
}


void
norm_eq_lsq(const double *xptr, int64_t m,            // x, shape (m,)
              const double *tptr, int64_t len_t,        // t, shape (len_t,)
              int k,
              const double *yptr, int64_t ydim2,        // y, shape(m, ydim2)
              const double *wptr,                       // w, shape (m,)
              /* outputs */
              double *abT_ptr,                          // ab, shape (k+1, m) IN FORTRAN ORDER
              double *rhs_ptr,                          // rhs, shape (m, ydim2)
              double *wrk
)
{
    auto x = ConstRealArray1D(xptr, m);
    auto t = ConstRealArray1D(tptr, len_t);
    auto y = ConstRealArray2D(yptr, m, ydim2);
    auto w = ConstRealArray1D(wptr, m);

    auto rhs = RealArray2D(rhs_ptr, m, ydim2);   // C order
    auto abT = RealArray2D(abT_ptr, m, k+1);  // transposed for the F order

    int64_t left = k;
    for(int64_t j=0; j < m ; j++) {
        double xval = x(j);
        double wval = w(j) * w(j);

        // where in t is xval
        left = _find_interval(t.data, len_t, k, xval, k, 0);

        // non-zero b-splines at xval
        _deBoor_D(t.data, xval, k, left, 0, wrk);

        // Fill non-zero values of A.T @ A: in the banded storage w/ lower=True and F order
        // With full matrices, the colloc matrix would have been
        //    A[j, left-k:left+1] = wrk
        // Here we work out A.T @ A *in the banded storage* w/ lower=True,
        // see the docstring of `scipy.linalg.cholesky_banded` for details
        int64_t r, s, row, clmn, ci;

        for (r=0; r < k+1; r++) {
            row = left - k + r;
            for (s=0; s < r+1; s++) {
                clmn = left - k + s;
                abT(clmn, r-s) += wrk[r] * wrk[s] * wval;   // NB: rows/cols swapped for F order
            }

            // ... rhs = A.T @ y
            for (ci=0; ci < ydim2; ci++) {
                rhs(row, ci) += wrk[r] * y(j, ci) * wval;
            }
        }
    }
}


/*** NDBSpline ***/

/* Evaluate an N-dim tensor product spline or its derivative */
void
_evaluate_ndbspline(const double *xi_ptr, int64_t npts, int64_t ndim,  // xi, shape(npts, ndim)
                    const double *t_ptr, int64_t max_len_t,            // t, shape (ndim, max_len_t)
                    const int64_t *len_t_ptr,                          // len_t, shape (ndim,)
                    const int64_t *k_ptr,                              // k, shape (ndim,)
                    const int64_t *nu_ptr,                             // nu, shape (ndim,)
                    int i_extrap,
                    const double *c1_ptr,  int64_t num_c1,             // flattened coefficients
                    // pre-tabulated helpers for iterating over (k+1)**ndim subarrays
                    const int64_t *strides_c1_ptr,                           // shape (ndim,)
                    const int64_t *indices_k1d_ptr,  int64_t num_k1d,        // shape (num_k1, ndim)
                    double *out_ptr, int64_t num_c_tr    // out, shape(npts, num_c_tr)
)
{
    auto xi = ConstRealArray2D(xi_ptr, npts, ndim);
    auto t = ConstRealArray2D(t_ptr, ndim, max_len_t);
    auto len_t = ConstIndexArray1D(len_t_ptr, ndim);
    auto k = ConstIndexArray1D(k_ptr, ndim);
    auto nu = ConstIndexArray1D(nu_ptr, ndim);
    auto c1 = ConstRealArray1D(c1_ptr, num_c1);
    auto strides_c1 = ConstIndexArray1D(strides_c1_ptr, ndim);
    auto indices_k1d = ConstIndexArray2D(indices_k1d_ptr, num_k1d, ndim);
    auto out = RealArray2D(out_ptr, npts, num_c_tr);

    // allocate work arrays (small, allocations unlikely to fail)
    int64_t max_k = *std::max_element(k_ptr, k_ptr + ndim);
    std::vector<double> wrk(2*max_k + 2);
    std::vector<int64_t> i(ndim);

    std::vector<double> v_b(ndim * (max_k + 1));
    auto b = RealArray2D(v_b.data(), ndim, max_k + 1);

    // the number of non-zero terms for each point in ``xi``
    int64_t volume = 1;
    for (int d=0; d < ndim; d++) {
        volume *= k(d) + 1;
    }

    // Iterate over the data points
    for (int64_t j=0; j < npts; j++){

        // For each point, iterate over the dimensions
        bool out_of_bounds = false;
        for(int d=0; d < ndim; d++) {
            double xd = xi(j, d);
            int64_t kd = k(d);

            // knots in the dimension d
            const double *td = t.data + max_len_t*d;

            // get the location of x[d] in td
            int64_t i_d = _find_interval(td, len_t(d), kd, xd, kd, i_extrap);

            if (i_d < 0) {
                out_of_bounds = true;
                break;
            }

            // compute non-zero b-splines at this value of xd in dimension d
            _deBoor_D(td, xd, kd, i_d, nu(d), wrk.data());

            for (int s=0; s < kd + 1; s++) {
                b(d, s) = wrk[s];
            }
            i[d] = i_d;
        } // for (d=...

        if (out_of_bounds) {
            // xd was nan or extrapolate=False: Fill the output array
            // for this data point, xi(j, :), and continue to the next xv in xi.

            for (int i_c=0; i_c < num_c_tr; i_c++) {
                out(j, i_c) = std::numeric_limits<double>::quiet_NaN();
            }
            continue;
        }

        // proceed to combining non-zero terms
        for(int i_c=0; i_c < num_c_tr; i_c++) {
            out(j, i_c) = 0;
        }

        // iterate over the direct product of non-zero b-splines
        for (int64_t iflat=0; iflat < volume; iflat++) {
            /* `idx_b = indiced_k1d[iflat, :]` assignment is equivalent to
             * idx_b = np.unravel_index(iflat, (k+1,)*ndim)
             * i.e. `idx_b` would be an ndim-dimensional index corresponding to
             * `iflat`.
             *
             * From the indices in ``idx_b``, we prepare to index into
             * c1.ravel() : for each dimension d, need to shift the index
             * by ``i[d] - k[d]`` (see the docstring above).
             *
             * Since the strides of `c1` are pre-computed, and the array
             * is already raveled and is guaranteed to be C-ordered, we only
             * need to compute the base index for iterating over ``num_c_tr``
             * elements which represent the trailing dimensions of ``c``.
             *
             * This all is essentially equivalent to iterating over
             * idx_cflat = np.ravel_multi_index(tuple(idx_c) + (i_c,),
             *                                  c1.shape)
             */
            int64_t idx_cflat_base = 0;
            double factor = 1.0;
            for (int d=0; d < ndim; d++) {
                int64_t idx_d = indices_k1d(iflat, d);
                factor *= b(d, idx_d);
                int64_t idx = idx_d + i[d] - k(d);
                idx_cflat_base += idx * strides_c1(d);
            }

            // finally, collect linear combinations of coef * factor
            for (int i_c=0; i_c < num_c_tr; i_c++) {
                out(j, i_c) += c1(idx_cflat_base + i_c) * factor;
            }
        }
    } // for (j=...
}


/*
 * Construct the N-D tensor product collocation matrix as a CSR array
 * Return value is 0 on a normal return, and negative on error:
 * if the data point `j` is problematic, return `-j`.
 */
int
_coloc_nd(/* inputs */
          const double *xi_ptr, int64_t npts, int64_t ndim,  // xi, shape(npts, ndim)
          const double *t_ptr, int64_t max_len_t,            // t, shape (ndim, max_len_t)
          const int64_t *len_t_ptr,                          // len_t, shape (ndim,)
          const int64_t *k_ptr,                              // k, shape (ndim,)
          /* pre-tabulated helpers for iterating over (k+1)**ndim subarrays */
          const int64_t *indices_k1d_ptr, int64_t num_k1d,        // shape (num_k1, ndim)
          const int64_t *strides_c1_ptr,                          // shape (ndim,)
          /* outputs */
          int64_t *csr_indices_ptr, int64_t volume,               // shape (npts*volume,)
          double *csr_data_ptr
)
{
    auto xi = ConstRealArray2D(xi_ptr, npts, ndim);
    auto t = ConstRealArray2D(t_ptr, ndim, max_len_t);
    auto len_t = ConstIndexArray1D(len_t_ptr, ndim);
    auto k = ConstIndexArray1D(k_ptr, ndim);

    auto strides_c1 = ConstIndexArray1D(strides_c1_ptr, ndim);
    auto indices_k1d = ConstIndexArray2D(indices_k1d_ptr, num_k1d, ndim);

    auto csr_indices = IndexArray1D(csr_indices_ptr, npts*volume);
    auto csr_data = RealArray1D(csr_data_ptr, npts*volume);

    // allocate work arrays (small, allocations unlikely to fail)
    int64_t max_k = *std::max_element(k_ptr, k_ptr + ndim);
    std::vector<double> wrk(2*max_k + 2);
    std::vector<int64_t> i(ndim);

    std::vector<double> v_b(ndim * (max_k + 1));
    auto b = RealArray2D(v_b.data(), ndim, max_k + 1);

    // Iterate over the data points
    for (int64_t j=0; j < npts; j++){

        // For each point, iterate over the dimensions
        bool out_of_bounds = false;
        for(int d=0; d < ndim; d++) {
            double xd = xi(j, d);
            int64_t kd = k(d);

            // knots in the dimension d
            const double *td = t.data + max_len_t*d;

            // get the location of x[d] in td
            int64_t i_d = _find_interval(td, len_t(d), kd, xd, kd, 1);

            if (i_d < 0) {
                out_of_bounds = true;
                break;
            }

            // compute non-zero b-splines at this value of xd in dimension d
            _deBoor_D(td, xd, kd, i_d, 0, wrk.data());

            for (int s=0; s < kd + 1; s++) {
                b(d, s) = wrk[s];
            }
            i[d] = i_d;
        } // for (d=...

        if (out_of_bounds) {
            // bail out
            return -j;
        }

        // Iterate over the products of non-zero b-splines and place them
        // into the current row of the design matrix
        for (int64_t iflat=0; iflat < volume; iflat++) {
            // The `idx_cflat` computation is an unrolled version of
            // idx_cflat = np.ravel_multi_index(tuple(idx_c), c_shape)
            //
            // `_indiced_k1d` array is pre-tabulated such that `idx_d` is a d-th component
            // of `idx_b = np.unravel_index(iflat,  tuple(kd+1 for kd in k))`
            int64_t idx_cflat = 0;
            double factor = 1.0;
            for (int d=0; d < ndim; d++) {
                int64_t idx_d = indices_k1d(iflat, d);
                factor *= b(d, idx_d);
                int64_t idx = idx_d + i[d] - k(d);
                idx_cflat += idx * strides_c1(d);
            }

            /*
             *  Fill the row of the colocation matrix in the CSR format.
             * If it were dense, it would have been just
             * >>> matr[j, idx_cflat] = factor
             *
             * Each row of the full matrix has `volume` non-zero elements.
             * Thus the CSR format `indptr` increases in steps of `volume`
             */
            csr_indices(j*volume + iflat) = idx_cflat;
            csr_data(j*volume + iflat) = factor;
        }  // for (iflat=...
    } // for( j=...

    return 0;
}

/**
 * Initialize augmented matrices for periodic B-spline fitting system.
 *
 * This function prepares the augmented banded matrices used in FITPACK's
 * periodic spline fitting algorithm. It expands and reorganizes the original
 * banded matrices (a1, a2, b) into augmented forms (g1, g2, h1, h2) that
 * incorporate periodic boundary conditions by linking the "head" and "tail"
 * parts of the spline basis. The offset array stores index adjustments for
 * wrap-around effects in the periodic system.
 *
 * Parameters:
 * a1ptr: (len_t - k - 1) x (k + 1)
 *        Original main banded matrix `a1` data.
 *        Represents the core system matrix for spline fitting normal equations.
 * a2ptr: (len_t - 2*k - 1) x k
 *        Original secondary banded matrix `a2` data.
 *        Contains coupling coefficients related to periodic boundary conditions.
 * bptr: (len_t - 2*k - 2) x (k + 2)
 *       Periodic constraints matrix `b`.
 *       Encodes the periodic tail-head interactions in the spline basis.
 * k: Spline degree.
 * len_t: Length of the knot vector.
 * g1ptr: (len_t - 2*k - 1) x (k + 2)
 *        Output augmented matrix `g1`.
 *        Augmented version of `a1` matrix, extended for periodicity.
 * g2ptr: (len_t - 2*k - 1) x (k + 1)
 *        Output augmented matrix `g2`.
 *        Augmented version of `a2` with an extra zero column for periodic effects.
 * h1ptr: (len_t - 2*k - 2) x (k + 2)
 *        Output matrix `h1`.
 *        Holds periodic contributions related to the "head" part of spline basis.
 * h2ptr: (len_t - 2*k - 2) x (k + 1)
 *        Output matrix `h2`.
 *        Holds periodic contributions complementary to `h1`.
 * offsetptr: (len_t - 2*k - 2)
 *            Output offset array.
 *            Contains index offsets used for wrap-around indexing in periodic system.
 */
void init_augmented_matrices(
    double *a1ptr, double *a2ptr, double *bptr,
    int k, int64_t len_t,
    double *g1ptr, double *g2ptr,
    double *h1ptr, double *h2ptr, double *offsetptr) {
    // ----------------------------------------------------------------------------
    // Augmented matrix construction for periodic B-spline fitting
    //
    // This function creates new matrices (g1, g2, h1, h2) from the original
    // banded matrices (a1, a2, b). These new matrices help in handling
    // periodic boundary conditions in B-spline fitting.
    //
    // In periodic splines, the start and end of the spline curve are connected.
    // This means the basis functions near the end of the knot vector also affect
    // the beginning of the curve — like wrapping around in a circle.
    // To handle this properly in matrix form, we build special "augmented" matrices
    // that include this wrapping information.
    //
    // Here's how each matrix is formed:
    //
    // - a1 → g1:
    //     g1 is just a copy of a1, but with one extra column at the end (filled with zeros).
    //     This extra column is used to hold wrapped values if needed later.
    //
    // - a2 → g2:
    //     g2 is like a2, but with one extra column at the beginning (filled with zeros).
    //     Also, some values from the end of a1 are copied into the bottom of g2 to handle
    //     wrapping from the end back to the start.
    //
    // - b → h1 and h2:
    //     b contains information about how the end of the spline interacts with the start.
    //     This is split into:
    //       → h1: stores the first part of each row of b (normally placed values)
    //       → h2: stores the remaining values from b that would overflow (wrap around)
    //
    // - offset:
    //     For each row, we store an offset value that tells us how the wrapping was done,
    //     so that later steps in the algorithm can apply corrections at the right places.
    //
    // --------------------------
    // Example (assume degree k = 2):
    //
    // Input shapes:
    //   a1:  (len_t - k - 1)     x (k + 1)
    //   a2:  (len_t - 2k - 1)    x (k)
    //   b:   (len_t - 2k - 2)    x (k + 2)
    //
    // Output shapes:
    //   g1:  (len_t - 2k - 1)    x (k + 2)
    //   g2:  (len_t - 2k - 1)    x (k + 1)
    //   h1:  (len_t - 2k - 2)    x (k + 2)
    //   h2:  (len_t - 2k - 2)    x (k + 1)
    //   offset: (len_t - 2k - 2)
    //
    // ASCII illustration of values:
    //
    //   a1 = [ a a a ]
    //   a2 = [   a a ]
    //   b  = [ b b b b ]   ← periodic influence terms (end → start)
    //
    //   g1 = [ a a a 0 ]   ← copy of a1 with extra column
    //   g2 = [ 0 a a ]     ← copy of a2 with zero column added in front,
    //                        and bottom rows filled from a1 for wrap-around
    //
    //   h1 = [ b b b b ]   ← first part of b copied here
    //   h2 = [ b b b ]     ← overflow from b (after wrap-around) goes here
    //
    // ----------------------------------------------------------------------------
    // Why is this necessary?
    //
    // If we don't do this augmentation, the spline system will miss connections
    // between the end and the beginning of the curve. That means:
    //
    // - The periodic property will break: the curve won't join smoothly at the ends
    // - The matrix will be incomplete or incorrect: wrapping basis functions
    //   will be placed in the wrong location or lost
    //
    //   What does this mean?
    //   Suppose we are fitting a periodic spline and x lies near the end of the domain.
    //   Its non-zero basis functions are supposed to "wrap around" and influence
    //   the start of the spline (because of periodicity).
    //
    //   Without the augmented matrices (h1, h2), these wrapped basis function values
    //   will either:
    //     (a) fall outside the bounds of the main system matrix (A, a1, a2), or
    //     (b) overwrite unrelated entries, or
    //     (c) be completely dropped during matrix assembly
    //
    //   As a result, the linear system we solve will be missing those terms —
    //   meaning, we're not enforcing periodicity properly, and the spline
    //   won't be continuous at the boundary.
    //
    //   h1 and h2 act like "extension slots" that safely store these out-of-bound
    //   or wrapped values in the correct position, so the final matrix still
    //   represents all required interactions correctly.
    // - Solving the system will give wrong results or fitting errors near boundaries
    //
    // These augmented matrices ensure the periodic conditions are correctly applied
    // and the final fitted spline is smooth and continuous across the domain.
    //
    // ----------------------------------------------------------------------------
    // Why do we need g1 and g2?
    //
    // g1 and g2 are extended versions of a1 and a2.
    // They include extra columns to handle the special case of periodic splines,
    // where the end of the spline connects back to the beginning.
    //
    // - g1 is a copy of a1 with one extra column added at the end (filled with 0)
    // - g2 has one extra column at the beginning (filled with 0),
    //   and a few rows at the bottom are filled from the end of a1
    //
    // If we do not create g1 and g2 like this, here is what goes wrong:
    //
    // - Some values that are supposed to "wrap around" from the end to the start
    //   will not fit in the original a1 or a2 matrix
    //
    // - These values may be:
    //     - written in the wrong place,
    //     - go outside the matrix and get ignored,
    //     - or overwrite other important values
    //
    // - This breaks the structure of the system we are trying to solve
    // - The result will be a spline that is not truly periodic —
    //   it may have a sharp jump at the start or end
    //
    // In short, g1 and g2 give us extra space in the right place
    // to safely store these wrap-around values, so the final matrix is correct.
    // ----------------------------------------------------------------------------



    // Wrap raw pointers as 2D arrays of appropriate sizes
    // These matrices come from banded system representations in FITPACK
    auto g1 = RealArray2D(g1ptr, len_t - 2*k - 1, k + 2);
    auto g2 = RealArray2D(g2ptr, len_t - 2*k - 1, k + 1);
    auto a1 = RealArray2D(a1ptr, len_t - k - 1, k + 1);
    auto a2 = RealArray2D(a2ptr, len_t - 2*k - 1, k);
    auto h1 = RealArray2D(h1ptr, len_t - 2*k - 2, k + 2);
    auto h2 = RealArray2D(h2ptr, len_t - 2*k - 2, k + 1);
    auto b = RealArray2D(bptr, len_t - 2*k - 2, k + 2);
    auto offset = RealArray1D(offsetptr, len_t - 2*k - 2);

    int64_t l0, l, l1;

    // --------------------------------------------------------------------
    // Copy matrix a1 into g1, extending the original banded system matrix.
    //
    // a1:
    //  - Represents the main banded matrix from the spline normal equations
    //  - Size corresponds to spline basis dimension minus degree adjustments
    //
    // g1:
    //  - Augmented matrix version of a1, extended with an extra column (k+2 instead of k+1)
    //  - Purpose: to incorporate periodic boundary conditions by increasing band size
    // --------------------------------------------------------------------
    for( int64_t i = 0; i < len_t - 2*k - 1; i++ ) {
        for( int64_t j = 0; j < k + 1; j++ ) {
            g1(i, j) = a1(i, j);
        }
    }
    // Last column is zero-padded because g1 has one extra column compared to a1
    for( int64_t i = 0; i < len_t - 2*k - 1; i++ ) {
        g1(i, k + 1) = 0.0;
    }

    // --------------------------------------------------------------------
    // Initialize first column of g2 to zero and fill remaining from a2
    //
    // a2:
    //  - Represents off-diagonal banded matrix terms in the original system
    //  - Corresponds to interactions related to the periodic part of the spline
    //
    // g2:
    //  - Augmented version of a2 with an additional zero column at start
    //  - First column zero to separate periodic contribution
    // --------------------------------------------------------------------
    for( int64_t i = 0; i < len_t - 2*k - 1; i++ ) {
        g2(i, 0) = 0.0;
    }
    for( int64_t i = 0; i < len_t - 2*k - 1; i++ ) {
        for( int64_t j = 1; j < k + 1; j++ ) {
            g2(i, j) = a2(i, j - 1);
        }
    }

    // --------------------------------------------------------------------
    // Special case: Fill first column of g2 at the bottom with values from a1
    // This enforces the periodic wrap-around conditions on the system matrix
    // by linking last rows back to first columns.
    //
    // This is critical in periodic splines where basis functions "wrap around".
    // --------------------------------------------------------------------
    l = len_t - 3*k - 1;
    for( int64_t j = 0; j < k + 1; j++ ) {
        if( l <= 0 ) {
            break;
        }
        g2(l - 1, 0) = a1(l - 1, j);
        l = l - 1;
    }

    // --------------------------------------------------------------------
    // Initialize h1, h2, and offset matrices which represent
    // the "periodic" part of the augmented system.
    //
    // h1, h2:
    //  - These matrices correspond to coefficients that arise from the
    //    periodic boundary conditions and wrapping in the banded matrix system.
    //  - They hold contributions from the "tail" of the spline basis interacting
    //    with the "head" due to periodicity.
    //
    // b:
    //  - Input matrix containing the periodic constraints and interactions.
    //
    // offset:
    //  - Stores the integer shifts needed for correct indexing in
    //    the augmented periodic system (wrap-around shifts).
    // --------------------------------------------------------------------
    for( int64_t it = 1; it <= len_t - 2*k - 2; it++ ) {
        // Zero out current row for h1 and h2
        for( int64_t i = 0; i < k + 1; i++ ) {
            h1(it - 1, i) = 0;
            h2(it - 1, i) = 0;
        }
        h1(it - 1, k + 1) = 0;

        if( it <= len_t - 3*k - 2 ) {
            // Fill h1 and h2 with rows from b matrix while handling wrap-around
            l = it;
            l0 = it;
            for( int64_t j = 1; j <= k + 2; j++ ) {
                if( l0 == len_t - 3*k - 1 ) {
                    // If boundary reached, wrap and continue filling h2 from b
                    l0 = 1;
                    for( int64_t l1 = j; l1 <= k + 2; l1++ ) {
                        h2(it - 1, l0 - 1) = b(it - 1, l1 - 1);
                        l0 = l0 + 1;
                    }
                    break;
                }
                // Normal assignment to h1 from b
                h1(it - 1, j - 1) = b(it - 1, j - 1);
                l0 = l0 + 1;
            }
        } else {
            // For rows beyond boundary, distribute elements from b to h1 and h2
            // according to periodic wrap-around indexing logic
            l = 1;
            int64_t i = it - (len_t - 3*k - 1);
            for( int64_t j = 1; j <= k + 2; j++ ) {
                i = i + 1;
                l0 = i;
                l1 = l0 - (k + 1);

                // Adjust indices for wrap-around within allowed bounds
                while( l1 > std::max((int64_t) 0, len_t - 3*k - 2) ) {
                    l0 = l1 - (len_t - 3*k - 2);
                    l1 = l0 - (k + 1);
                }

                if( l1 > 0 ) {
                    h1(it - 1, l1 - 1) = b(it - 1, j - 1);
                } else {
                    // Accumulate values in h2 for remaining wrapped positions
                    h2(it - 1, l0 - 1) += b(it - 1, j - 1);
                }
            }
        }
        // Store the offset for this row for index calculations in periodic solver
        offset(it - 1) = l;
    }
}


} // namespace fitpack
