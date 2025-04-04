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
       const double *yptr, int64_t ydim2,            // y(m, ydim2)
        /* output */
       double *cptr)                                 // c(nc, ydim2)
{
    auto R = ConstRealArray2D(Rptr, m, nz);
    auto y = ConstRealArray2D(yptr, m, ydim2);
    auto c = RealArray2D(cptr, nc, ydim2);

    // c[nc-1, ...] = y[nc-1] / R[nc-1, 0]
    for (int64_t l=0; l < ydim2; ++l) {
        c(nc - 1, l) = y(nc - 1, l) / R(nc-1, 0);
    }

    //for i in range(nc-2, -1, -1):
    //    nel = min(nz, nc-i)
    //    c[i, ...] = ( y[i] - (R[i, 1:nel, None] * c[i+1:i+nel, ...]).sum(axis=0) ) / R[i, 0]
    for (int64_t i=nc-2; i >= 0; --i) {
        int64_t nel = std::min(nz, nc - i);
        for (int64_t l=0; l < ydim2; ++l){
            double ssum = y(i, l);
            for (int64_t j=1; j < nel; ++j) {
                ssum -= R(i, j) * c(i + j, l);
            }
            ssum /= R(i, 0);
            c(i, l) = ssum;
        }
    }
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


} // namespace fitpack
