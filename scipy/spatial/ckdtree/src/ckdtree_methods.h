
#ifndef CKDTREE_CPP_METHODS
#define CKDTREE_CPP_METHODS

#ifdef CKDTREE_METHODS_IMPL
#define CKDTREE_EXTERN extern "C"
#else
#define CKDTREE_EXTERN extern "C" 
struct ckdtree;
#endif 

extern int number_of_processors;

#ifndef NPY_LIKELY
#define NPY_LIKELY(x) (x)
#endif

#ifndef NPY_UNLIKELY
#define NPY_UNLIKELY(x) (x)
#endif

#include <cmath>
#include <vector>
#include "numpy/npy_math.h"
#include "ordered_pair.h"
#include "coo_entries.h"

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

#endif // _WIN32
#endif // __GNUC__

/*
 * Utility functions
 * =================
 */
 
#define ckdtree_isinf(x) (x == NPY_INFINITY)
 
inline npy_float64 
dmax(const npy_float64 x, const npy_float64 y) 
{
    if (x > y)
        return x;
    else
        return y;
};

inline npy_float64 
dmin(const npy_float64 x, const npy_float64 y) 
{
    if (x < y)
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
wrap_distance(const npy_float64 x, const npy_float64 hb, const npy_float64 fb)
{
    npy_float64 x1;
    if (NPY_UNLIKELY(x < -hb)) x1 = fb + x;
    else if (NPY_UNLIKELY(x > hb)) x1 = x - fb;
    else x1 = x;
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
sqeuclidean_distance_double(const npy_float64 *u, const npy_float64 *v, npy_intp n)
{
    npy_float64 s;
    npy_intp i;
    // manually unrolled loop, might be vectorized
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


inline const npy_float64 
_wrap(const npy_float64 x, const npy_float64 box) 
{
    const npy_float64 r = std::floor(x / box);
    return x - r * box;
}

/* FIXME: this should be replaced by a function
 * in MinMaxDistBox */
const inline npy_float64 side_distance_from_min_max(
    const npy_float64 x,
    const npy_float64 min,
    const npy_float64 max,
    const npy_float64 p,
    const npy_float64 hb,
    const npy_float64 fb) 
{
    npy_float64 s, t, tmin, tmax;

    if (NPY_LIKELY(fb <= 0)) {
        /* non-periodic */
        s = 0; 
        t = x - max;
        if (t > s) {
            s = t;
        } else {
            t = min - x;
            if (t > s) s = t;
        }
    } else {
        /* periodic */
        s = 0;
        tmax = x - max;
        tmin = x - min;
        /* is the test point in this range */
        if(NPY_UNLIKELY(tmax >= 0 || tmin <= 0)) {
            /* no */
            tmax = dabs(tmax);
            tmin = dabs(tmin);
            /* tmin will be the closer edge */
            if(tmin > tmax) {
                t = tmin;
                tmin = tmax;
                tmax = t;
            }
            if(tmax < hb) {
                /* both edges are less than half a box. */
                /* no wrapping, use the closer edge */ 
                s = tmin;
            } else if(tmin > hb) {
                /* both edge are more than half a box. */
                /* wrapping on both edge, use the 
                 * wrapped further edge */
                s = fb - tmax;
            } else {
                /* the further side is wrapped */
                tmax = fb - tmax;
                if(tmin > tmax) {
                    s = tmax;
                } else {
                    s = tmin;
                }
            }
        } else {
            /* yes. min distance is 0 */
            s = 0; 
        }
    }

    if (NPY_UNLIKELY(p == 1. || ckdtree_isinf(p))) {
        s = dabs(s);
    } else if (NPY_LIKELY(p == 2.)) {
        s = s * s;
    } else {
        s = std::pow(s,p);
    }
        
    return s;
}

/* Build methods in C++ for better speed and GIL release */

CKDTREE_EXTERN PyObject*
build_ckdtree(ckdtree *self, npy_intp start_idx, npy_intp end_idx,
              npy_float64 *maxes, npy_float64 *mins, int _median, int _compact);

extern "C" PyObject*
build_weights (ckdtree *self, npy_float64 *node_weights, npy_float64 *weights);

/* Query methods in C++ for better speed and GIL release */

CKDTREE_EXTERN PyObject*
query_knn(const ckdtree     *self, 
          npy_float64       *dd, 
          npy_intp          *ii, 
          const npy_float64 *xx,
          const npy_intp     n,
          const npy_intp     *k, 
          const npy_intp     nk, 
          const npy_intp     kmax, 
          const npy_float64  eps, 
          const npy_float64  p, 
          const npy_float64  distance_upper_bound);
          
CKDTREE_EXTERN PyObject*
query_pairs(const ckdtree *self, 
            const npy_float64 r, 
            const npy_float64 p, 
            const npy_float64 eps,
            std::vector<ordered_pair> *results);
            
CKDTREE_EXTERN PyObject*
count_neighbors_unweighted(const ckdtree *self,
                const ckdtree *other,
                npy_intp n_queries,
                npy_float64 *real_r,
                npy_intp *results,
                const npy_float64 p,
                int cumulative);

CKDTREE_EXTERN PyObject*
count_neighbors_weighted(const ckdtree *self,
                const ckdtree *other,
                npy_float64 *self_weights, 
                npy_float64 *other_weights, 
                npy_float64 *self_node_weights, 
                npy_float64 *other_node_weights, 
                npy_intp n_queries,
                npy_float64 *real_r,
                npy_float64 *results,
                const npy_float64 p,
                int cumulative);
                               
CKDTREE_EXTERN PyObject*
query_ball_point(const ckdtree *self,
                 const npy_float64 *x,
                 const npy_float64 r,
                 const npy_float64 p,
                 const npy_float64 eps,
                 const npy_intp n_queries,
                 std::vector<npy_intp> **results);
                 
CKDTREE_EXTERN PyObject*                
query_ball_tree(const ckdtree *self,
                const ckdtree *other,
                const npy_float64 r,
                const npy_float64 p,
                const npy_float64 eps,
                std::vector<npy_intp> **results);                

CKDTREE_EXTERN PyObject*                 
sparse_distance_matrix(const ckdtree *self,
                       const ckdtree *other,
                       const npy_float64 p,
                       const npy_float64 max_distance,
                       std::vector<coo_entry> *results);

                  
#endif


