#include <cstring>
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
ssize_t
_find_interval(const double* tptr, ssize_t len_t,
               int k,
               double xval,
               ssize_t prev_l,
               int extrapolate)
{
    ConstRealArray1D t = ConstRealArray1D(tptr, len_t);

    ssize_t n = t.nelem - k - 1;
    double tb = t(k);
    double te = t(n);

    if (xval != xval) {
        // nan
        return -1;
    }

    if (((xval < tb) || (xval > te)) && !extrapolate) {
        return -1;
    }
    ssize_t l = (k < prev_l) && (prev_l < n) ? prev_l : k;

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
            const double *xptr, ssize_t m,      // x, shape (m,)
            const double *tptr, ssize_t len_t,  // t, shape (len_t,)
            int k,
            const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
            /* outputs */
            double *Aptr,                       // A, shape(m, k+1)
            ssize_t *offset_ptr,                // offset, shape (m,)
            ssize_t *nc,                        // the number of coefficient
            /* work array*/
            double *wrk)                        // work, shape (2k+2)
{
    auto x = ConstRealArray1D(xptr, m);
    auto t = ConstRealArray1D(tptr, len_t);
    auto w = ConstRealArray1D(wptr, m);
    auto A = RealArray2D(Aptr, m, k+1);
    auto offset = Array1D<ssize_t, false>(offset_ptr, m);

    ssize_t ind = k;
    for (int i=0; i < m; ++i) {
        double xval = x(i);

        // find the interval
        ind = _find_interval(t.data, len_t, k, xval, ind, 0);
        if (ind < 0){
            // should not happen here, validation is expected on the python side
            throw std::runtime_error("find_interval: out of bounds with x = " + std::to_string(xval));
        }
        offset(i) = ind - k;

        // compute non-zero b-splines
        _deBoor_D(t.data, xval, k, ind, 0, wrk);

        for (ssize_t j=0; j < k+1; ++j) {
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
qr_reduce(double *aptr, const ssize_t m, const ssize_t nz, // a(m, nz), packed
          ssize_t *offset,                                 // offset(m)
          const ssize_t nc,                                // dense would be a(m, nc)
          double *yptr, const ssize_t ydim1,               // y(m, ydim2)
          const ssize_t startrow
)
{
    auto R = RealArray2D(aptr, m, nz);
    auto y = RealArray2D(yptr, m, ydim1);

    for (ssize_t i=startrow; i < m; ++i) {
        ssize_t oi = offset[i];
        for (ssize_t j=oi; j < nc; ++j) {

            // rotate only the lower diagonal
            if (j >= std::min(i, nc)) {
                break;
            }

            // in dense format: diag a1[j, j] vs a1[i, j]
            double c, s, r;
            DLARTG(&R(j, 0), &R(i, 0), &c, &s, &r);

            // rotate l.h.s.
            R(j, 0) = r;
            for (ssize_t l=1; l < R.ncols; ++l) {
                std::tie(R(j, l), R(i, l-1)) = fprota(c, s, R(j, l), R(i, l));
            }
            R(i, R.ncols-1) = 0.0;

            // rotate r.h.s.
            for (ssize_t l=0; l < y.ncols; ++l) {
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
       const double *Rptr, ssize_t m, ssize_t nz,    // R(m, nz), packed
       ssize_t nc,                                   // dense R would be (m, nc)
       const double *yptr, ssize_t ydim2,            // y(m, ydim2)
        /* output */
       double *cptr)                                 // c(nc, ydim2)
{
    auto R = ConstRealArray2D(Rptr, m, nz);
    auto y = ConstRealArray2D(yptr, m, ydim2);
    auto c = RealArray2D(cptr, nc, ydim2);

    // c[nc-1, ...] = y[nc-1] / R[nc-1, 0]
    for (ssize_t l=0; l < ydim2; ++l) {
        c(nc - 1, l) = y(nc - 1, l) / R(nc-1, 0);
    }

    //for i in range(nc-2, -1, -1):
    //    nel = min(nz, nc-i)
    //    c[i, ...] = ( y[i] - (R[i, 1:nel, None] * c[i+1:i+nel, ...]).sum(axis=0) ) / R[i, 0]
    for (ssize_t i=nc-2; i >= 0; --i) {
        ssize_t nel = std::min(nz, nc - i);
        for (ssize_t l=0; l < ydim2; ++l){
            double ssum = y(i, l);
            for (ssize_t j=1; j < nel; ++j) {
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
    ssize_t interval = k+1;
    ssize_t nc = t.nelem - k - 1;

    std::vector<ssize_t> ix;
    ix.push_back(0);

    std::vector<double> fparts;
    double fpart = 0.0;

    for(ssize_t i=0; i < x.nelem; i++) {
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
fpknot(const double *x_ptr, ssize_t m,
       const double *t_ptr, ssize_t len_t,
       int k,
       const double *residuals_ptr)
{
    auto x = ConstRealArray1D(x_ptr, m);
    auto t = ConstRealArray1D(t_ptr, len_t);
    auto residuals = ConstRealArray1D(residuals_ptr, m);

    std::vector<double> fparts;
    std::vector<ssize_t> ix;
    std::tie(fparts, ix) = _split(x, t, k, residuals);

    ssize_t idx_max = -101;
    double fpart_max = -1e100;
    for (size_t i=0; i < fparts.size(); i++) {
        bool is_better = (ix[i+1] - ix[i] > 1) && (fparts[i] > fpart_max);
        if(is_better) {
            idx_max = i;
            fpart_max = fparts[i];
        }
    }

    if (idx_max == -101) {
        throw std::runtime_error("Internal error. Please report it to SciPy developers.");
    }

    // round up, like Dierckx does? This is really arbitrary though.
    ssize_t idx_newknot = (ix[idx_max] + ix[idx_max+1] + 1) / 2;

    return x(idx_newknot);
}

} // namespace fitpack
