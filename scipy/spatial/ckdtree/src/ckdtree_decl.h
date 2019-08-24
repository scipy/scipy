#ifndef CKDTREE_CPP_DECL
#define CKDTREE_CPP_DECL

/*
 * Use numpy to provide some platform independency.
 * Define these functions for your platform
 * */
#include <cmath> /* needed for isinf / sc_inf from c99compat under CLANG */
#include "_c99compat.h"
#include <numpy/npy_common.h>
#define CKDTREE_LIKELY(x) NPY_LIKELY(x)
#define CKDTREE_UNLIKELY(x)  NPY_UNLIKELY(x)
#define CKDTREE_PREFETCH(x, rw, loc)  NPY_PREFETCH(x, rw, loc)

#define ckdtree_intp_t npy_intp
#define ckdtree_isinf(x)   sc_isinf(x)
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
                 std::vector<ckdtree_intp_t> **results,
                 const int return_length);

int
query_ball_tree(const ckdtree *self,
                const ckdtree *other,
                const double r,
                const double p,
                const double eps,
                std::vector<ckdtree_intp_t> **results
                );

int
sparse_distance_matrix(const ckdtree *self,
                       const ckdtree *other,
                       const double p,
                       const double max_distance,
                       std::vector<coo_entry> *results);


#endif
