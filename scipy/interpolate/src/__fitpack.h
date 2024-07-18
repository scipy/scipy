#pragma once
#include <iostream>
#include <tuple>
#include <vector>
#include <string>

#include "../_build_utils/src/npy_cblas.h"
#include "../_build_utils/src/fortran_defs.h"

#define DLARTG BLAS_FUNC(dlartg)

/* MSVC */
#define ssize_t ptrdiff_t

namespace fitpack {


/*
 * 1D and 2D array wrappers, with and without boundschecking
 */


// Bounds checking
template<bool boundscheck> inline void _bcheck(ssize_t index, ssize_t size, ssize_t dim);


template<>
inline void _bcheck<true>(ssize_t index, ssize_t size, ssize_t dim) {
    if (!((0 <= index) && (index < size))){
        auto mesg = "Out of bounds with index = " + std::to_string(index) + " of size = ";
        mesg = mesg + std::to_string(size) + " in dimension = " + std::to_string(dim);
        throw(std::runtime_error(mesg) );
    }
}

template<>
inline void _bcheck<false>(ssize_t index, ssize_t size, ssize_t dim) { /* noop*/ }


// Arrays: C contiguous only
template<typename T, bool boundscheck=true>
struct Array1D
{
    T* data;
    ssize_t nelem;
    T& operator()(const ssize_t i) {
        _bcheck<boundscheck>(i, nelem, 0);
        return *(data + i);
    }
    Array1D(T *ptr, ssize_t num_elem) : data(ptr), nelem(num_elem) {};
};



template<typename T, bool boundscheck=true>
struct Array2D
{
    T* data;
    ssize_t nrows;
    ssize_t ncols;
    T& operator()(const ssize_t i, const ssize_t j) {
        _bcheck<boundscheck>(i, nrows, 0);
        _bcheck<boundscheck>(j, ncols, 1);
        return *(data + ncols*i + j);
    }
    Array2D(T *ptr, ssize_t num_rows, ssize_t num_columns) : data(ptr), nrows(num_rows), ncols(num_columns) {};
};


// Flip boundschecking on/off here
constexpr bool BOUNDS_CHECK = false;

typedef Array2D<double, BOUNDS_CHECK> RealArray2D;
typedef Array1D<double, BOUNDS_CHECK> RealArray1D;
typedef Array1D<const double, BOUNDS_CHECK> ConstRealArray1D;
typedef Array2D<const double, BOUNDS_CHECK> ConstRealArray2D;


/*
 * LAPACK Givens rotation
 */
extern "C" {
void DLARTG(double *f, double *g, double *cs, double *sn, double *r);
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
 * B-spline evaluation routine.
 */
void
_deBoor_D(const double *t, double x, int k, int ell, int m, double *result);


/*
 *  Find an interval such that t[interval] <= xval < t[interval+1].
 */
ssize_t
_find_interval(const double* tptr, ssize_t len_t,
               int k,
               double xval,
               ssize_t prev_l,
               int extrapolate);



/*
 * Fill the (m, k+1) matrix of non-zero b-splines.
 */
void
data_matrix(/* inputs */
            const double *xptr, ssize_t m,      // x, shape (m,)
            const double *tptr, ssize_t len_t,  // t, shape (len_t,)
            int k,
            const double *wptr,                 // weights, shape (m,) // NB: len(w) == len(x), not checked
            /* outputs */
            double *Aptr,                       // A, shape(m, k+1)
            ssize_t *offset_ptr,                // offset, shape (m,)
            ssize_t *nc,                        // the number of coefficient
            /* work array*/
            double *wrk                         // work, shape (2k+2)
);


/*
    Solve the LSQ problem ||y - A@c||^2 via QR factorization.
    This routine MODIFIES `a` & `y` in-place.
*/
void
qr_reduce(double *aptr, const ssize_t m, const ssize_t nz, // a(m, nz), packed
          ssize_t *offset,                                 // offset(m)
          const ssize_t nc,                                // dense would be a(m, nc)
          double *yptr, const ssize_t ydim1,               // y(m, ydim2)
          const ssize_t startrow=1
);


/*
 * Back substitution solve of `R @ c = y` with an upper triangular R.
 */
void
fpback( /* inputs*/
       const double *Rptr, ssize_t m, ssize_t nz,    // R(m, nz), packed
       ssize_t nc,                                   // dense R would be (m, nc)
       const double *yptr, ssize_t ydim2,            // y(m, ydim2)
        /* output */
       double *cptr                                 // c(nc, ydim2)
);


/*
 * A helper for _fpknot:
 * Split the `x` array into knot "runs" and sum the residuals per "run".
 */
typedef std::tuple<std::vector<double>, std::vector<ssize_t>> pair_t;

pair_t
_split(ConstRealArray1D x, ConstRealArray1D t, int k, ConstRealArray1D residuals);


/*
 * Find a position for a new knot, a la FITPACK
 */
double
fpknot(const double *x_ptr, ssize_t m,
       const double *t_ptr, ssize_t len_t,
       int k,
       const double *residuals_ptr);


} // namespace fitpack
