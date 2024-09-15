#pragma once
#include <iostream>
#include <cinttypes>
#include <tuple>
#include <vector>
#include <string>
#include <limits>

#include "../_build_utils/src/npy_cblas.h"
#include "../_build_utils/src/fortran_defs.h"

#define DLARTG BLAS_FUNC(dlartg)

namespace fitpack {

// use int64_t for indexing throughout
// Using system dependent int types (ssize_t, ptrdiff_t, Py_ssize_t and npy_intp)
// is a world of pain on all of Windows; 32-bit & 64-bit linux & numpy 1.x.
typedef int64_t d_ssize_t;


/*
 * LAPACK Givens rotation
 */
extern "C" {
void BLAS_FUNC(dlartg)(double *f, double *g, double *cs, double *sn, double *r);
}


/*
 * Apply a Givens transform.
 */
template<typename T>
inline
std::tuple<T, T>
fprota(T c, T s, T a, T b)
{
    return std::make_tuple(
        c*a + s*b,
       -s*a + c*b
    );
}


/*
 * 1D and 2D array wrappers, with and without boundschecking
 */


// Bounds checking
template<bool boundscheck> inline void _bcheck(d_ssize_t index, d_ssize_t size, d_ssize_t dim);


template<>
inline void _bcheck<true>(d_ssize_t index, d_ssize_t size, d_ssize_t dim) {
    if (!((0 <= index) && (index < size))){
        auto mesg = "Out of bounds with index = " + std::to_string(index) + " of size = ";
        mesg = mesg + std::to_string(size) + " in dimension = " + std::to_string(dim);
        throw(std::runtime_error(mesg) );
    }
}

template<>
inline void _bcheck<false>(d_ssize_t index, d_ssize_t size, d_ssize_t dim) { /* noop*/ }


// Arrays: C contiguous only
template<typename T, bool boundscheck=true>
struct Array1D
{
    T* data;
    d_ssize_t nelem;
    T& operator()(const d_ssize_t i) {
        _bcheck<boundscheck>(i, nelem, 0);
        return *(data + i);
    }
    Array1D(T *ptr, d_ssize_t num_elem) : data(ptr), nelem(num_elem) {};
};



template<typename T, bool boundscheck=true>
struct Array2D
{
    T* data;
    d_ssize_t nrows;
    d_ssize_t ncols;
    T& operator()(const d_ssize_t i, const d_ssize_t j) {
        _bcheck<boundscheck>(i, nrows, 0);
        _bcheck<boundscheck>(j, ncols, 1);
        return *(data + ncols*i + j);
    }
    Array2D(T *ptr, d_ssize_t num_rows, d_ssize_t num_columns) : data(ptr), nrows(num_rows), ncols(num_columns) {};
};


// Flip boundschecking on/off here
constexpr bool BOUNDS_CHECK = true;

typedef Array2D<double, BOUNDS_CHECK> RealArray2D;
typedef Array1D<double, BOUNDS_CHECK> RealArray1D;
typedef Array1D<const double, BOUNDS_CHECK> ConstRealArray1D;
typedef Array2D<const double, BOUNDS_CHECK> ConstRealArray2D;



/*
 * B-spline evaluation routine.
 */
void
_deBoor_D(const double *t, double x, int k, int ell, int m, double *result);


/*
 *  Find an interval such that t[interval] <= xval < t[interval+1].
 */
d_ssize_t
_find_interval(const double* tptr, d_ssize_t len_t,
               int k,
               double xval,
               d_ssize_t prev_l,
               int extrapolate);



/*
 * Fill the (m, k+1) matrix of non-zero b-splines.
 */
void
data_matrix(/* inputs */
            const double *xptr, d_ssize_t m,      // x, shape (m,)
            const double *tptr, d_ssize_t len_t,  // t, shape (len_t,)
            int k,
            const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
            /* outputs */
            double *Aptr,                       // A, shape(m, k+1)
            d_ssize_t *offset_ptr,                // offset, shape (m,)
            d_ssize_t *nc,                        // the number of coefficient
            /* work array*/
            double *wrk                         // work, shape (2k+2)
);


/*
    Solve the LSQ problem ||y - A@c||^2 via QR factorization.
    This routine MODIFIES `a` & `y` in-place.
*/
void
qr_reduce(double *aptr, const d_ssize_t m, const d_ssize_t nz, // a(m, nz), packed
          d_ssize_t *offset,                                 // offset(m)
          const d_ssize_t nc,                                // dense would be a(m, nc)
          double *yptr, const d_ssize_t ydim1,               // y(m, ydim2)
          const d_ssize_t startrow=1
);


/*
 * Back substitution solve of `R @ c = y` with an upper triangular R.
 */
void
fpback( /* inputs*/
       const double *Rptr, d_ssize_t m, d_ssize_t nz,    // R(m, nz), packed
       d_ssize_t nc,                                   // dense R would be (m, nc)
       const double *yptr, d_ssize_t ydim2,            // y(m, ydim2)
        /* output */
       double *cptr                                 // c(nc, ydim2)
);


/*
 * A helper for _fpknot:
 * Split the `x` array into knot "runs" and sum the residuals per "run".
 */
typedef std::tuple<std::vector<double>, std::vector<d_ssize_t>> pair_t;

pair_t
_split(ConstRealArray1D x, ConstRealArray1D t, int k, ConstRealArray1D residuals);


/*
 * Find a position for a new knot, a la FITPACK
 */
double
fpknot(const double *x_ptr, d_ssize_t m,
       const double *t_ptr, d_ssize_t len_t,
       int k,
       const double *residuals_ptr);


/*
 * Evaluate the spline function
 */

void
_evaluate_spline(
    const double *tptr, d_ssize_t len_t,         // t, shape (len_t,)
    const double *cptr, d_ssize_t n, d_ssize_t m,  // c, shape (n, m)
    d_ssize_t k,
    const double *xp_ptr, d_ssize_t s,           // xp, shape (s,)
    d_ssize_t nu,
    int extrapolate,
    double *out_ptr,                           // out, shape (s, m) NOT CHECKED
    double *wrk                                // scratch, shape (2k+2,)
);


/*
 * Spline collocation matrix in the LAPACK banded storage
 */
void
_coloc_matrix(const double *xptr, d_ssize_t m,        // x, shape(m,)
               const double *tptr, d_ssize_t len_t,   // t, shape(len_t,)
               int k,
               double *abT, d_ssize_t nbands,         // ab(nbands, len_t - k - 1) in F order!
               int offset,
               double *wrk                          // scratch, shape (2k+2,)
);


/*
 * Construct the l.h.s. and r.h.s of the normal equations for the LSQ spline fitting.
 */
void
norm_eq_lsq(const double *xptr, d_ssize_t m,      // x, shape (m,)
              const double *tptr, d_ssize_t len_t,  // t, shape (len_t,)
              int k,
              const double *yptr, d_ssize_t ydim2,  // y, shape(m, ydim2)
              const double *wptr,                 // w, shape (m,)
              /* outputs */
              double *abT_ptr,                    // ab, shape (k+1, m) IN FORTRAN ORDER
              double *rhs_ptr,                    // rhs, shape (m, ydim2)
              double *wrk
);

} // namespace fitpack
