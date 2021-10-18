# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE_NOGIL=1

# Code adapted from github.com/adrianveres/Polo, licensed:
#
# The MIT License (MIT)
# Copyright (c) 2016 Adrian Veres
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free

from scipy.spatial.distance import squareform, is_valid_y, is_valid_dm

np.import_array()

@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void dual_swap(float* darr, int* iarr,
                           int i1, int i2):
    """
    [Taken from Scikit-learn.]

    swap the values at inex i1 and i2 of both darr and iarr"""
    cdef float dtmp = darr[i1]
    darr[i1] = darr[i2]
    darr[i2] = dtmp

    cdef int itmp = iarr[i1]
    iarr[i1] = iarr[i2]
    iarr[i2] = itmp


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _simultaneous_sort(float* dist, int* idx,
                            int size) except -1:
    """
    [Taken from Scikit-learn.]


    Perform a recursive quicksort on the dist array, simultaneously
    performing the same swaps on the idx array.  The equivalent in
    numpy (though quite a bit slower) is
    def simultaneous_sort(dist, idx):
        i = np.argsort(dist)
        return dist[i], idx[i]
    """
    cdef int pivot_idx, i, store_idx
    cdef float pivot_val

    # in the small-array case, do things efficiently
    if size <= 1:
        pass
    elif size == 2:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
    elif size == 3:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
        if dist[1] > dist[2]:
            dual_swap(dist, idx, 1, 2)
            if dist[0] > dist[1]:
                dual_swap(dist, idx, 0, 1)
    else:
        # Determine the pivot using the median-of-three rule.
        # The smallest of the three is moved to the beginning of the array,
        # the middle (the pivot value) is moved to the end, and the largest
        # is moved to the pivot index.
        pivot_idx = size // 2
        if dist[0] > dist[size - 1]:
            dual_swap(dist, idx, 0, size - 1)
        if dist[size - 1] > dist[pivot_idx]:
            dual_swap(dist, idx, size - 1, pivot_idx)
            if dist[0] > dist[size - 1]:
                dual_swap(dist, idx, 0, size - 1)
        pivot_val = dist[size - 1]

        # partition indices about pivot.  At the end of this operation,
        # pivot_idx will contain the pivot value, everything to the left
        # will be smaller, and everything to the right will be larger.
        store_idx = 0
        for i in range(size - 1):
            if dist[i] < pivot_val:
                dual_swap(dist, idx, i, store_idx)
                store_idx += 1
        dual_swap(dist, idx, store_idx, size - 1)
        pivot_idx = store_idx

        # recursively sort each side of the pivot
        if pivot_idx > 1:
            _simultaneous_sort(dist, idx, pivot_idx)
        if pivot_idx + 2 < size:
            _simultaneous_sort(dist + pivot_idx + 1,
                               idx + pivot_idx + 1,
                               size - pivot_idx - 1)
    return 0


cdef inline void _sort_M_slice(float[:, ::1] M,
                               float* vals, int* idx,
                               int dim1_min, int dim1_max, int dim2_val):
    """
    Simultaneously sort indices and values of M[{m}, u] using
    `_simultaneous_sort`

    This is equivalent to :
       m_sort = M[dim1_min:dim1_max, dim2_val].argsort()
       m_iter = np.arange(dim1_min, dim1_max)[m_sort]

    but much faster because we don't have to pay the numpy overhead. This
    matters a lot for the sorting of M[{k}, w] which is executed many times.
    """
    cdef int i
    for i in range(0, dim1_max - dim1_min):
        vals[i] = M[dim1_min + i, dim2_val]
        idx[i] = dim1_min + i
    _simultaneous_sort(vals, idx, dim1_max - dim1_min)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:] identify_swaps(int[:, ::1] sorted_Z,
                           double[:, ::1] sorted_D,
                           int[:, ::1] cluster_ranges):
    """
    Implements the Optimal Leaf Ordering algorithm described in
    "Fast Optimal leaf ordering for hierarchical clustering"
        Ziv Bar-Joseph, David K. Gifford, Tommi S. Jaakkola
        Bioinformatics, 2001, :doi:`10.1093/bioinformatics/17.suppl_1.S22`

    `sorted_Z` : Linkage list, with 'height' column removed.

    """
    cdef int n_points = len(sorted_Z) + 1

    cdef:
        # (n x n) floats
        float[:, ::1] M = np.zeros((n_points, n_points), dtype=np.float32)
        # (n x n x 2) booleans
        int[:, :, :] swap_status = np.zeros((n_points, n_points, 2),
                                            dtype=np.intc)
        int[:] must_swap = np.zeros((len(sorted_Z),), dtype=np.intc)

        int i, v_l, v_r, v_size,
        int v_l_min, v_l_max, v_r_min, v_r_max

        int u_clusters[2]
        int m_clusters[2]
        int w_clusters[2]
        int k_clusters[2]
        int total_u_clusters, total_w_clusters

        int u, w, m, k
        int u_min, u_max, m_min, m_max, w_min, w_max, k_min, k_max
        int swap_L, swap_R

        float* m_vals
        int* m_idx
        float* k_vals
        int* k_idx
        int mi, ki

        float min_km_dist
        float cur_min_M, current_M
        int best_m, best_k

        int best_u = 0, best_w = 0

    for i in range(len(sorted_Z)):
        # Iterate over the linkage list instead of recursion.
        #     v_l = sorted_Z[i, 0]
        #     v_r = sorted_Z[i, 1]
        # are indices of the left and right children for node i.
        #
        # If the v_l or v_r are < n_points, then v_l or v_r are singleton
        # clusters. Otherwise, it is the node defined in the (i - n_points)
        #
        #                V
        #              /   \
        #           --       --
        #          /            \
        #        V_l            V_r
        #       /   \          /   \
        #     V_l1  V_l2     V_r1  V_r2
        #     (u)   (m)      (k)   (w)
        #
        # Briefly, for every node V, the algorithm finds left-most and
        # right-most nodes u, w that minimizes U[u, w] the sum of distances of
        # every neighboring singleton node in the linear ordering.
        #
        # This is done recursively, by finding the optimizing the ordering of
        #     v_l (bounded by nodes u, m) and v_r (bounded by k, w)
        # such that
        #     U[u, w] = U[u, m] + U[k, w] + D[m, k]
        # is then minimized.
        #
        # Part of the optimization is that at every search step,
        # if (u) ~ V_l1, then (m) ~ V_l2 (and vice-versa)
        # likewise for (w) ~ V_r1, then (w) ~ V_r2.
        #
        # This means we need to search 4 pairs of (V_li, V_rj) combinations.
        # If V_l or V_r are singletons, for example, V_l = u = m.

        v_l = sorted_Z[i, 0]
        v_r = sorted_Z[i, 1]
        v_size = sorted_Z[i, 2]

        v_l_min = cluster_ranges[v_l, 0]; v_l_max = cluster_ranges[v_l, 1]
        v_r_min = cluster_ranges[v_r, 0]; v_r_max = cluster_ranges[v_r, 1]

        if v_l < n_points:
            # V_l is a singleton, so U = M = V_L.
            total_u_clusters = 1

            # This could be handled more efficiently, but in practice the code
            # would get longer for no speed gain.
            u_clusters[0] = v_l
            m_clusters[0] = v_l

        else:
            total_u_clusters = 2

            # First look for U from V_LL and M from V_LR
            u_clusters[0] = sorted_Z[v_l - n_points, 0]
            m_clusters[0] = sorted_Z[v_l - n_points, 1]

            # Then look for U from V_LR and M from V_LL
            u_clusters[1] = sorted_Z[v_l - n_points, 1]
            m_clusters[1] = sorted_Z[v_l - n_points, 0]

        if v_r < n_points:
            total_w_clusters = 1
            # V_r is a singleton, so W = K = V_R.
            w_clusters[0] = v_r
            k_clusters[0] = v_r

        else:
            total_w_clusters = 2

            # First look for W from V_RR and L from V_RL
            w_clusters[0] = sorted_Z[v_r - n_points, 1]
            w_clusters[1] = sorted_Z[v_r - n_points, 0]

            # Next look for W from V_RL and L from V_RR
            k_clusters[0] = sorted_Z[v_r - n_points, 0]
            k_clusters[1] = sorted_Z[v_r - n_points, 1]

        for swap_L in range(total_u_clusters):
            for swap_R in range(total_w_clusters):
                # Get bounds for the clusters from which we'll sample u, m, w, k
                # (see note above for details).
                # If in the chosen ordering,
                #     U came from V_ll : Don't swap V_l.
                #     U came from V_lr : Swap V_l.
                #     W came from V_rl : Swap V_r.
                #     W came from V_ll : Don't swap V_r.

                u_min = cluster_ranges[u_clusters[swap_L], 0]
                u_max = cluster_ranges[u_clusters[swap_L], 1]
                m_min = cluster_ranges[m_clusters[swap_L], 0]
                m_max = cluster_ranges[m_clusters[swap_L], 1]
                w_min = cluster_ranges[w_clusters[swap_R], 0]
                w_max = cluster_ranges[w_clusters[swap_R], 1]
                k_min = cluster_ranges[k_clusters[swap_R], 0]
                k_max = cluster_ranges[k_clusters[swap_R], 1]

                # Find the minimum of D[m, k] for the appropriate sets {m}, {k}.
                # This is C[{m}, {k}] in the paper's notation.
                min_km_dist = 1073741824 #2^30
                for m in range(m_min, m_max):
                    for k in range(k_min, k_max):
                        if sorted_D[m, k] < min_km_dist:
                            min_km_dist = sorted_D[m, k]

                m_vals = <float*>malloc(sizeof(float) * (m_max - m_min))
                m_idx = <int*>malloc(sizeof(int) * (m_max - m_min))
                k_vals = <float*>malloc(sizeof(float) * (k_max - k_min))
                k_idx = <int*>malloc(sizeof(int) * (k_max - k_min))
                if not m_vals or not m_idx or not k_vals or not k_idx:
                    free(m_vals)
                    free(m_idx)
                    free(k_vals)
                    free(k_idx)
                    raise MemoryError("failed to allocate memory in identify_swaps().")

                for u in range(u_min, u_max):
                    # Sort the values of M[{m}, u]
                    _sort_M_slice(M, m_vals, m_idx, m_min, m_max, u)

                    for w in range(w_min, w_max):
                        # Sort the values of M[{k}, w]
                        _sort_M_slice(M, k_vals, k_idx, k_min, k_max, w)

                        # Set initial value for cur_min_M.
                        # I used a large number.
                        cur_min_M = 1073741824.0 #2^30

                        for mi in range(0, m_max - m_min):
                            m = m_idx[mi]

                            if (M[u, m] + M[w, k_idx[0]] + min_km_dist
                                    >= cur_min_M):
                                # Terminate the outer loop early, there will not
                                # be a better 'k' in the current k list.
                                break
                            for ki in range(0, k_max - k_min):
                                k = k_idx[ki]

                                if M[u, m] + M[w, k] + min_km_dist >= cur_min_M:
                                    # Terminate the inner loop early
                                    break

                                current_M = M[u, m] + M[w, k] + sorted_D[m, k]
                                if current_M < cur_min_M:
                                    # We found a better m, k than previously.
                                    cur_min_M = current_M
                                    best_m = m
                                    best_k = k

                        # For the chosen (u, w), record the resulting minimal
                        # M[u, w] = M[u, m] + M[k, w] + D[m, k]
                        M[u, w] = cur_min_M
                        M[w, u] = cur_min_M
                        # whether we need to swap V_l and V_r given the current
                        # chosen (m, k) (see note above). This saves us from
                        # storing (m, k) and doing back-tracking later.
                        swap_status[u, w, 0] = swap_L
                        swap_status[w, u, 0] = swap_L
                        swap_status[u, w, 1] = swap_R
                        swap_status[w, u, 1] = swap_R

                # We are getting a fresh `w` and `u` so need to resort
                # M[{k}, w] and M[{m}, u]
                free(m_vals)
                free(m_idx)
                free(k_vals)
                free(k_idx)

        # We are now ready to find the best minimal value for M[{u}, {w}]
        cur_min_M = 1073741824.0 #2^30
        for u in range(v_l_min, v_l_max):
            for w in range(v_r_min, v_r_max):
                if M[u, w] < cur_min_M:
                    cur_min_M = M[u, w]
                    best_u = u
                    best_w = w

        # If v_l, v_r are not singletons, record whether our choice of (u, w)
        # for V requires a swap of its children.
        if v_l >= n_points:
            must_swap[v_l - n_points] = int(swap_status[best_u, best_w, 0])
        if v_r >= n_points:
            must_swap[v_r - n_points] = int(swap_status[best_u, best_w, 1])

    return must_swap


def optimal_leaf_ordering(Z, D):
    """
    Compute the optimal leaf order for Z (according to D) and return an
    optimally sorted Z.

    We start by sorting and relabelling Z and D according to the current leaf
    order in Z.

    This is because when everything is sorted each cluster (including
    singletons) can be defined by its range over (0...n_points).

    This is used extensively to loop efficiently over the various arrays in the
    algorithm.

    """
    # Import here to avoid import cycles
    from scipy.cluster.hierarchy import leaves_list, is_valid_linkage

    is_valid_linkage(Z, throw=True, name='Z')

    if is_valid_y(D):
        sorted_D = squareform(D)
    elif is_valid_dm(D):
        sorted_D = D
    else:
        raise("Not a valid distance matrix (neither condensed nor square form)")

    n_points = Z.shape[0] + 1
    n_clusters = 2*n_points - 1

    # Get the current linear ordering
    sorted_leaves = leaves_list(Z)

    # Create map from original order to sorted order.
    original_order_to_sorted_order = dict((orig_i, sorted_i) for sorted_i,orig_i
                                          in enumerate(sorted_leaves))


    # Re-write linkage map so it refers to sorted positions, rather than input
    # positions. Remove the 'height' column so we can cast the whole thing as
    # integer and simplify passing to C function above.
    sorted_Z = []
    for (v_l, v_r, _, v_size) in Z:
        if v_l < n_points:
            v_l = original_order_to_sorted_order[int(v_l)]
        if v_r < n_points:
            v_r = original_order_to_sorted_order[int(v_r)]

        sorted_Z.append([v_l, v_r, v_size])
    sorted_Z = np.array(sorted_Z).astype(np.int32).copy(order='C')


    # Sort distance matrix D by the leaf order
    sorted_D = sorted_D[sorted_leaves, :]
    sorted_D = sorted_D[:, sorted_leaves].copy(order='C')

    # Defines the range of each cluster over (0... n_points) as explained above.
    cluster_ranges = np.zeros((n_clusters, 2))
    cluster_ranges[np.arange(n_points), 0] = np.arange(n_points)
    cluster_ranges[np.arange(n_points), 1] = np.arange(n_points) + 1
    for link_i, (v_l, v_r, v_size) in enumerate(sorted_Z):
        v = link_i + n_points
        cluster_ranges[v, 0] = cluster_ranges[v_l, 0]
        cluster_ranges[v, 1] = cluster_ranges[v_r, 1]
    cluster_ranges = cluster_ranges.astype(np.int32).copy(order='C')

    # Get Swaps
    must_swap = identify_swaps(sorted_Z, sorted_D, cluster_ranges)

    # To 'rotate' around the axis of a node, we need to consider the left-right
    # children of every descendant of this target node.
    #
    # To do so efficiently, we record how many total times a given node must be
    # swapped (once if it needs to be swapped itself, once for each parent that
    # needs to be swapped) and take modulo 2 to find whether it needs to be
    # swapped at all.
    is_descendant = np.zeros((n_clusters - n_points, n_clusters - n_points),
                             dtype=int)
    for i, (v_l, v_r, v_size) in enumerate(sorted_Z):
        is_descendant[i, i] = 1
        if v_l >= n_points:
            is_descendant[i, v_l - n_points] = 1
            is_descendant[i, :] += is_descendant[v_l - n_points, :]
        if v_r >= n_points:
            is_descendant[i, v_r - n_points] = 1
            is_descendant[i, :] += is_descendant[v_r - n_points, :]


    # To "rotate" a tree node, we need to 'swap' its left-right children,
    # and do the same to all its children.
    applied_swap = (np.array(is_descendant).astype(bool)
                    * np.array(must_swap).reshape(-1, 1))
    final_swap = applied_swap.sum(axis=0) % 2

    # Create a new linkage matrix by applying swaps where needed.
    swapped_Z = []
    for i, (in_l, in_r, h, v_size) in enumerate(Z):
        if final_swap[i]:
            out_l = in_r
            out_r = in_l
        else:
            out_r = in_r
            out_l = in_l
        swapped_Z.append((out_l, out_r, h, v_size))
    swapped_Z = np.array(swapped_Z)

    return swapped_Z
