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

// NB. We use int64_t for indexing throughout.
// Using system dependent int types (int64_t, ptrdiff_t, Py_int64_t and npy_intp)
// is a world of pain on all of Windows; 32-bit & 64-bit linux & numpy 1.x.


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
template<bool boundscheck> inline void _bcheck(int64_t index, int64_t size, int dim);


template<>
inline void _bcheck<true>(int64_t index, int64_t size, int dim) {
    if (!((0 <= index) && (index < size))){
        auto mesg = "Out of bounds with index = " + std::to_string(index) + " of size = ";
        mesg = mesg + std::to_string(size) + " in dimension = " + std::to_string(dim);
        throw(std::runtime_error(mesg) );
    }
}

template<>
inline void _bcheck<false>(int64_t index, int64_t size, int dim) { /* noop*/ }


// Arrays: C contiguous only
template<typename T, bool boundscheck=true>
struct Array1D
{
    T* data;
    int64_t nelem;
    T& operator()(const int64_t i) {
        _bcheck<boundscheck>(i, nelem, 0);
        return *(data + i);
    }
    Array1D(T *ptr, int64_t num_elem) : data(ptr), nelem(num_elem) {};
};



template<typename T, bool boundscheck=true>
struct Array2D
{
    T* data;
    int64_t nrows;
    int64_t ncols;
    T& operator()(const int64_t i, const int64_t j) {
        _bcheck<boundscheck>(i, nrows, 0);
        _bcheck<boundscheck>(j, ncols, 1);
        return *(data + ncols*i + j);
    }
    Array2D(T *ptr, int64_t num_rows, int64_t num_columns) : data(ptr), nrows(num_rows), ncols(num_columns) {};
};


// Flip boundschecking on/off here
constexpr bool BOUNDS_CHECK = false;

typedef Array2D<double, BOUNDS_CHECK> RealArray2D;
typedef Array1D<double, BOUNDS_CHECK> RealArray1D;
typedef Array1D<const double, BOUNDS_CHECK> ConstRealArray1D;
typedef Array2D<const double, BOUNDS_CHECK> ConstRealArray2D;

typedef Array1D<int64_t, BOUNDS_CHECK> IndexArray1D;
typedef Array1D<const int64_t, BOUNDS_CHECK> ConstIndexArray1D;
typedef Array2D<const int64_t, BOUNDS_CHECK> ConstIndexArray2D;



/*
 * B-spline evaluation routine.
 */
void
_deBoor_D(const double *t, double x, int k, int ell, int m, double *result);


/*
 *  Find an interval such that t[interval] <= xval < t[interval+1].
 */
int64_t
_find_interval(const double* tptr, int64_t len_t,
               int k,
               double xval,
               int64_t prev_l,
               int extrapolate);


/*
 * Fill the (m, k+1) matrix of non-zero b-splines.
 */
void
data_matrix(/* inputs */
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
            double *wrk                         // work, shape (2k+2)
);


/*
    Solve the LSQ problem ||y - A@c||^2 via QR factorization.
    This routine MODIFIES `a` & `y` in-place.
*/
void
qr_reduce(double *aptr, const int64_t m, const int64_t nz, // a(m, nz), packed
          int64_t *offset,                                 // offset(m)
          const int64_t nc,                                // dense would be a(m, nc)
          double *yptr, const int64_t ydim1,               // y(m, ydim2)
          const int64_t startrow=1
);


/*
 * Back substitution solve of `R @ c = y` with an upper triangular R.
 */
void
fpback( /* inputs*/
       const double *Rptr, int64_t m, int64_t nz,    // R(m, nz), packed
       int64_t nc,                                   // dense R would be (m, nc)
       const double *yptr, int64_t ydim2,            // y(m, ydim2)
        /* output */
       double *cptr                                 // c(nc, ydim2)
);


/*
 * A helper for _fpknot:
 * Split the `x` array into knot "runs" and sum the residuals per "run".
 */
typedef std::tuple<std::vector<double>, std::vector<int64_t>> pair_t;

pair_t
_split(ConstRealArray1D x, ConstRealArray1D t, int k, ConstRealArray1D residuals);


/*
 * Find a position for a new knot, a la FITPACK
 */
double
fpknot(const double *x_ptr, int64_t m,
       const double *t_ptr, int64_t len_t,
       int k,
       const double *residuals_ptr);


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
);


/*
 * Spline collocation matrix in the LAPACK banded storage
 */
void
_coloc_matrix(const double *xptr, int64_t m,        // x, shape(m,)
               const double *tptr, int64_t len_t,   // t, shape(len_t,)
               int k,
               double *abT, int64_t nbands,         // ab(nbands, len_t - k - 1) in F order!
               int offset,
               double *wrk                          // scratch, shape (2k+2,)
);


/*
 * Construct the l.h.s. and r.h.s of the normal equations for the LSQ spline fitting.
 */
void
norm_eq_lsq(const double *xptr, int64_t m,      // x, shape (m,)
              const double *tptr, int64_t len_t,  // t, shape (len_t,)
              int k,
              const double *yptr, int64_t ydim2,  // y, shape(m, ydim2)
              const double *wptr,                 // w, shape (m,)
              /* outputs */
              double *abT_ptr,                    // ab, shape (k+1, m) IN FORTRAN ORDER
              double *rhs_ptr,                    // rhs, shape (m, ydim2)
              double *wrk
);


/*
 * Evaluate an ND spline function
 */
void
_evaluate_ndbspline(/* inputs */
                    const double *xi_ptr, int64_t npts, int64_t ndim,  // xi, shape(npts, ndim) 
                    const double *t_ptr, int64_t max_len_t,            // t, shape (ndim, max_len_t)
                    const int64_t *len_t_ptr,                          // len_t, shape (ndim,)
                    const int64_t *k_ptr,                              // k, shape (ndim,)
                    const int64_t *nu_ptr,                             // nu, shape (ndim,)
                    int i_extrap,
                    /* flattened coefficients */
                    const double *c1_ptr, int64_t num_c1,
                    /* pre-tabulated helpers for iterating over (k+1)**ndim subarrays */
                    const int64_t *strides_c1_ptr,                           // shape (ndim,)
                    const int64_t *indices_k1d_ptr,  int64_t num_k1d,        // shape (num_k1, ndim)
                    /* output */
                    double *out_ptr, int64_t num_c_tr            // out, shape(npts, num_c_tr)
);


int
_coloc_nd(/* inputs */
          const double *xi_ptr, int64_t npts, int64_t ndim,  // xi, shape(npts, ndim)
          const double *t_ptr, int64_t max_len_t,            // t, shape (ndim, max_len_t)
          const int64_t *len_t_ptr,                          // len_t, shape (ndim,)
          const int64_t *k_ptr,                              // k, shape (ndim,)
          /* pre-tabulated helpers for iterating over (k+1)**ndim subarrays */
          const int64_t *indices_k1d_ptr, int64_t num_k1d,   // shape (num_k1, ndim)
          const int64_t *strides_c1_ptr,                     // shape (ndim,)
          /* outputs */
          int64_t *csr_indices_ptr, int64_t volume,          // shape (npts*volume,)
          double *csr_data_ptr
);


} // namespace fitpack
