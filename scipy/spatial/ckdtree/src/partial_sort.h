#ifndef CKDTREE_PARTIAL_SORT
#define CKDTREE_PARTIAL_SORT

/* Adapted version of the code originally
 * written by @jiefangxuanyan for scikit-learn.
 *
 */

#include "ckdtree_decl.h"
#include <algorithm>

struct IndexComparator {

    const double *data;
    ckdtree_intp_t split_dim;
    ckdtree_intp_t n_dims;

    IndexComparator(const double *data,
                    ckdtree_intp_t split_dim,
                    ckdtree_intp_t n_dims) :
                                            data(data),
                                            split_dim(split_dim),
                                            n_dims(n_dims) {};

    inline bool operator()(ckdtree_intp_t a, ckdtree_intp_t b) {
        const double point_a = data[a * n_dims + split_dim];
        const double point_b = data[b * n_dims + split_dim];

        if CKDTREE_UNLIKELY (point_a == point_b) {
            return a < b;
        } else {
            return point_a < point_b;
        }
    }
};

/*
 * Partition points in the node into two groups.
 *
 */

static int
partition_node_indices(const double   *data,
                       ckdtree_intp_t *node_indices,
                       ckdtree_intp_t split_dim,
                       ckdtree_intp_t split_index,
                       ckdtree_intp_t n_dims,
                       ckdtree_intp_t n_points) {

    /* Upon return, the values in node_indices will be rearranged such that
     * (assuming numpy-style indexing):
     *
     *   data[node_indices[0:split_index], split_dim]
     *     <= data[node_indices[split_index], split_dim]
     *
     * and
     *
     *   data[node_indices[split_index], split_dim]
     *     <= data[node_indices[split_index:n_points], split_dim]
     *
     * This is eassentially a wrapper around the standard C++ function
     * ``std::nth_element``.
     *
     *
     * Parameters
     * ----------
     * data : double pointer
     *    Pointer to a 2D array of the training data, of shape [N, n_dims].
     *    N must be greater than any of the values in node_indices.
     * node_indices : int pointer
     *    Pointer to a 1D array of length n_points.  This lists the indices of
     *    each of the points within the current node.  This will be modified
     *    in-place.
     * split_dim : int
     *    the dimension on which to split.  This will usually be computed via
     *    the routine ``find_node_split_dim``
     * split_index : int
     *    the index within node_indices around which to split the points.
     *
     * Returns
     * -------
     * status : int
     *    integer exit status.  On return, the contents of node_indices are
     *    modified as noted above.
     */

    IndexComparator index_comparator(data, split_dim, n_dims);

    std::nth_element(node_indices,
                     node_indices + split_index,
                     node_indices + n_points,
                     index_comparator);

    return 0;
}

#endif
