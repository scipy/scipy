
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

#if defined(__GNUC__)
inline void 
prefetch_datapoint(const npy_float64 *x, const npy_intp m) 
{
    const int cache_line = 64;  // x86, amd64
    char *cur = (char*)x;
    char *end = (char*)(x+m);
    while (cur < end) { 
        __builtin_prefetch((void*)cur);
        cur += cache_line;
    }
}
#else
#if defined(_WIN32)
#include <xmmintrin.h>
inline void
prefetch_datapoint(const npy_float64 *x, const npy_intp m)
{
    const int cache_line = 64;  // x86, amd64
    char *cur = (char*)x;
    char *end = (char*)(x+m);
    while (cur < end) {
        _mm_prefetch((const char*)cur,_MM_HINT_T0);
        cur += cache_line;
    }
}
#else
#define prefetch_datapoint(x,y)
#endif
#endif

/*
 * Utility functions
 * =================
 */
 
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

inline npy_float64 
dabs_b(const npy_float64 x, const npy_float64 hb, const npy_float64 fb)
{
    npy_float64 x1;
    x1 = x;
    if (x < 0) x1 = -x;
    if (hb > 0 && x1 > hb)
        x1 = fb - x1;
#if 0
    printf("dabs_b x : %g x1 %g\n", x, x1);
#endif
    return x1;
}

/*
 * Measuring distances
 * ===================
 */


inline npy_float64
sqeuclidean_distance_double(const npy_float64 *u, const npy_float64 *v,
                            npy_intp n)
{
    npy_float64 s;
    npy_intp i;
    /* manually unrolled loop, might be vectorized */
    npy_float64 acc[4] = {0., 0., 0., 0.};
    for (i = 0; i < n/4; i += 4) {
        npy_float64 _u[4] = {u[i], u[i + 1], u[i + 2], u[i + 3]};
        npy_float64 _v[4] = {v[i], v[i + 1], v[i + 2], v[i + 3]};
        npy_float64 diff[4] = {_u[0] - _v[0], 
                               _u[1] - _v[1], 
                               _u[2] - _v[2],
                               _u[3] - _v[3]};
        acc[0] += diff[0] * diff[0];
        acc[1] += diff[1] * diff[1];
        acc[2] += diff[2] * diff[2];
        acc[3] += diff[3] * diff[3];
    }
    s = acc[0] + acc[1] + acc[2] + acc[3];
    if (i < n) {
        for(; i<n; ++i) {
            npy_float64 d = u[i] - v[i];
            s += d * d;
        }
    }
    return s;
} 
 
inline npy_float64 
_distance_p(const npy_float64 *x, const npy_float64 *y,
            const npy_float64 p, const npy_intp k,
            const npy_float64 upperbound, const ckdtreebox * box)
{    
   /*
    * Compute the distance between x and y
    *
    * Computes the Minkowski p-distance to the power p between two points.
    * If the distance**p is larger than upperbound, then any number larger
    * than upperbound may be returned (the calculation is truncated).
    */
    
    npy_intp i;
    npy_float64 r, r1;
    r = 0;
    if (NPY_LIKELY(p==2.)) {
        /*
        for (i=0; i<k; ++i) {
            z = x[i] - y[i];
            r += z*z;
            if (r>upperbound)
                return r;
        }*/
        if(box->fbox[0] <= 0) {
            return sqeuclidean_distance_double(x,y,k);
        } 
        for (i=0; i<k; ++i) {
            r1 = dabs_b(x[i] - y[i], box->hbox[i], box->fbox[i]);
            r += r1 * r1;
            if (r>upperbound)
                return r;
        }
    } 
    else if (p==infinity) {
        for (i=0; i<k; ++i) {
            r1 = dabs_b(x[i] - y[i], box->hbox[i], box->fbox[i]);
            r = dmax(r,r1);
            if (r>upperbound)
                return r;
        }
    } 
    else if (p==1.) {
        for (i=0; i<k; ++i) {
            r1 = dabs_b(x[i] - y[i], box->hbox[i], box->fbox[i]);
            r += r1;
            if (r>upperbound)
                return r;
        }
    } 
    else {
        for (i=0; i<k; ++i) {
            r1 = dabs_b(x[i] - y[i], box->hbox[i], box->fbox[i]);
            r += std::pow(r1,p);
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
          const npy_float64  distance_upper_bound,
          const ckdtreebox * box);
          

// Other query methods can follow here when they are implemented
          
#endif

