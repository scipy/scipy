
#ifndef CKDTREE_CPP_METHODS
#define CKDTREE_CPP_METHODS

#ifdef CKDTREE_METHODS_IMPL
#define CKDTREE_EXTERN extern "C"
#else
#define CKDTREE_EXTERN extern "C" 
struct ckdtree;
#endif 

extern npy_float64 infinity;
extern int number_of_processors;

#ifndef NPY_LIKELY
#define NPY_LIKELY(x) (x)
#endif

#include <cmath>

// Utility functions
// =================

inline npy_float64 
dmax(const npy_float64 x, const npy_float64 y) 
{
    if (x > y)
        return x;
    else
        return y;
};

inline npy_float64 
dabs(const npy_float64 x)
{
    if (x > 0)
        return x;
    else
        return -x;
}


// Measuring distances
// ===================

inline npy_float64 
_distance_p(const npy_float64 *x, const npy_float64 *y,
            const npy_float64 p, const npy_intp k,
            const npy_float64 upperbound)
{    
    // Compute the distance between x and y
    //
    // Computes the Minkowski p-distance to the power p between two points.
    // If the distance**p is larger than upperbound, then any number larger
    // than upperbound may be returned (the calculation is truncated).
    
    npy_intp i;
    npy_float64 r, z;
    r = 0;
    if (NPY_LIKELY(p==2.)) {
        for (i=0; i<k; ++i) {
            z = x[i] - y[i];
            r += z*z;
            if (r>upperbound)
                return r;
        }
    } else if (p==infinity) {
        for (i=0; i<k; ++i) {
            r = dmax(r,dabs(x[i]-y[i]));
            if (r>upperbound)
                return r;
        }
    } else if (p==1.) {
        for (i=0; i<k; ++i) {
            r += dabs(x[i]-y[i]);
            if (r>upperbound)
                return r;
        }
    } else {
        for (i=0; i<k; ++i) {
            r += std::pow(dabs(x[i]-y[i]),p);
            if (r>upperbound)
                 return r;
        }
    }
    return r;
} 

// k-nearest neighbor query
          
CKDTREE_EXTERN PyObject*
query_knn(const ckdtree     *self, 
          npy_float64       *dd, 
          npy_intp          *ii, 
          const npy_float64 *xx,
          const npy_intp     n,
          const npy_intp     k, 
          const npy_float64  eps, 
          const npy_float64  p, 
          const npy_float64  distance_upper_bound);
          

// Other query methods can follow here when they are 
// implemented
          
#endif

