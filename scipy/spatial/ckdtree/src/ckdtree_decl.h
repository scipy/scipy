#ifndef CKDTREE_CPP_DECL
#define CKDTREE_CPP_DECL

/*
 * Use numpy to provide some platform independence.
 * Define these functions for your platform
 * */
#include <numpy/npy_common.h>
#include <cmath>
#include <new>

#if defined(__cpp_lib_hardware_interference_size)
    #define CACHE_LINE std::hardware_destructive_interference_size
#else
    // 64 is OK for x86 and amd64, not for Apple silicon with 128
    // On Apple silicon: Use a C++17 library that has  std::size_t hardware_destructive_interference_size
    // defined or suffer the pain of extra prefetching.
    #define CACHE_LINE 64
#endif // __cpp_lib_hardware_interference_size

#if defined(__GNUC__) 

#define CKDTREE_LIKELY(expr)    (__builtin_expect(!!(expr), 1))
#define CKDTREE_UNLIKELY(expr)  (__builtin_expect(!!(expr), 0))

inline void 
prefetch_datapoint(const npy_float64 *x, const npy_intp m) 
{
    // The data point can live on multiple cache lines
    // so we must hit all of them
    char *cur = (char*)x;
    char *end = (char*)(x+m);    
    do { 
        __builtin_prefetch((void*)cur, 0, 3);
        cur += CACHE_LINE;
    } while (CKDTREE_UNLIKELY(cur < end));
}
  
#else 

#define CKDTREE_LIKELY(x) NPY_LIKELY(x)     // FIXME: Does not do anything as of May 2025, better replace with C++20 [[likely]]
#define CKDTREE_UNLIKELY(x) NPY_UNLIKELY(x) // FIXME: Does not do anything as of May 2025, better replace with C++20 [[unlikely]]
  
#if defined(_WIN32) 
  
#include <xmmintrin.h> 
  
inline void 
prefetch_datapoint(const npy_float64 *x, const npy_intp m) 
{
    // The data point can live on multiple cache lines
    // so we must hit all of them
    char *cur = (char*)x;
    char *end = (char*)(x+m);    
    do { 
        _mm_prefetch((const char*)cur,_MM_HINT_T0);
        cur += CACHE_LINE;
    } while (CKDTREE_UNLIKELY(cur < end));
}

#else 

#define prefetch_datapoint(x,m)

#endif // _WIN32
#endif // __GNUC__

// rw is always 0
#define CKDTREE_PREFETCH(x,rw,m) prefetch_datapoint(x,m)

#define ckdtree_intp_t npy_intp
#define ckdtree_fmin(x, y)   fmin(x, y)
#define ckdtree_fmax(x, y)   fmax(x, y)
#define ckdtree_fabs(x)   fabs(x)

#include "ordered_pair.h"
#include "coo_entries.h"

struct ckdtreenode {
    ckdtree_intp_t      split_dim;
    ckdtree_intp_t      children;
    double   split;
    ckdtree_intp_t      start_idx;
    ckdtree_intp_t      end_idx;
    ckdtreenode   *less;
    ckdtreenode   *greater;
    ckdtree_intp_t      _less;
    ckdtree_intp_t      _greater;
};

struct ckdtree {
    // tree structure
    std::vector<ckdtreenode>  *tree_buffer;
    ckdtreenode   *ctree;
    // meta data
    double   *raw_data;
    ckdtree_intp_t      n;
    ckdtree_intp_t      m;
    ckdtree_intp_t      leafsize;
    double   *raw_maxes;
    double   *raw_mins;
    ckdtree_intp_t      *raw_indices;
    double   *raw_boxsize_data;
    ckdtree_intp_t size;
};

/* Build methods in C++ for better speed and GIL release */

int
build_ckdtree(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
              double *maxes, double *mins, int _median, int _compact);

int
build_weights (ckdtree *self, double *node_weights, double *weights);

/* Query methods in C++ for better speed and GIL release */

int
query_knn(const ckdtree     *self,
          double       *dd,
          ckdtree_intp_t          *ii,
          const double *xx,
          const ckdtree_intp_t     n,
          const ckdtree_intp_t     *k,
          const ckdtree_intp_t     nk,
          const ckdtree_intp_t     kmax,
          const double  eps,
          const double  p,
          const double  distance_upper_bound);

int
query_pairs(const ckdtree *self,
            const double r,
            const double p,
            const double eps,
            std::vector<ordered_pair> *results);

int
count_neighbors_unweighted(const ckdtree *self,
                const ckdtree *other,
                ckdtree_intp_t n_queries,
                double *real_r,
                ckdtree_intp_t *results,
                const double p,
                int cumulative);

int
count_neighbors_weighted(const ckdtree *self,
                const ckdtree *other,
                double *self_weights,
                double *other_weights,
                double *self_node_weights,
                double *other_node_weights,
                ckdtree_intp_t n_queries,
                double *real_r,
                double *results,
                const double p,
                int cumulative);

int
query_ball_point(const ckdtree *self,
                 const double *x,
                 const double *r,
                 const double p,
                 const double eps,
                 const ckdtree_intp_t n_queries,
                 std::vector<ckdtree_intp_t> *results,
                 const bool return_length,
                 const bool sort_output);

int
query_ball_tree(const ckdtree *self,
                const ckdtree *other,
                const double r,
                const double p,
                const double eps,
                std::vector<ckdtree_intp_t> *results
                );

int
sparse_distance_matrix(const ckdtree *self,
                       const ckdtree *other,
                       const double p,
                       const double max_distance,
                       std::vector<coo_entry> *results);


#endif
