
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


