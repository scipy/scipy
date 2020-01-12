
#ifndef CKDTREE_PARTIAL_SORT
#define CKDTREE_PARTIAL_SORT

/* Splitting routines for a balanced kd-tree
 * Code originally written by Jake Vanderplas for scikit-learn
 *
 */

inline void
index_swap(ckdtree_intp_t *arr, intptr_t i1, intptr_t i2)
{
    /* swap the values at index i1 and i2 of arr */
    ckdtree_intp_t tmp = arr[i1];
    arr[i1] = arr[i2];
    arr[i2] = tmp;
}

static void
partition_node_indices(const double *data,
                       ckdtree_intp_t *node_indices,
                       ckdtree_intp_t split_dim,
                       ckdtree_intp_t split_index,
                       ckdtree_intp_t n_features,
                       ckdtree_intp_t n_points)
{
    /* Partition points in the node into two equal-sized groups
     * Upon return, the values in node_indices will be rearranged such that
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
     * The algorithm is essentially a partial in-place quicksort around a
     * set pivot.
     *
     * Parameters
     * ----------
     * data : double pointer
     *    Pointer to a 2D array of the training data, of shape [N, n_features].
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

    ckdtree_intp_t left, right, midindex, i;
    double d1, d2;
    left = 0;
    right = n_points - 1;
    for(;;) {
        midindex = left;
        for (i=left; i<right; ++i) {
            d1 = data[node_indices[i] * n_features + split_dim];
            d2 = data[node_indices[right] * n_features + split_dim];
            if (d1 < d2) {
                index_swap(node_indices, i, midindex);
                ++midindex;
            }
        }
        index_swap(node_indices, midindex, right);
        if (midindex == split_index)
            break;
        else if (midindex < split_index)
            left = midindex + 1;
        else
            right = midindex - 1;
    }
}

#endif
