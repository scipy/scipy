# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from libc.math cimport sqrt
from libc.string cimport memset
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef extern from "numpy/npy_math.h":
    cdef enum:
        NPY_INFINITYF

ctypedef unsigned char uchar


# _hierarchy_distance_update.pxi includes the definition of linkage_distance_update
# and the distance update functions for the supported linkage methods.
include "_hierarchy_distance_update.pxi"
cdef linkage_distance_update *linkage_methods = [
    _single, _complete, _average, _centroid, _median, _ward, _weighted]


cdef inline int condensed_index(int n, int i, int j):
    """
    Calculate the condensed index of element (i, j) in an n x n condensed
    matrix.
    """
    if i < j:
        return n * i - (i * (i + 1) / 2) + (j - i - 1)
    elif i > j:
        return n * j - (j * (j + 1) / 2) + (i - j - 1)


cdef inline int is_visited(uchar *bitset, int i):
    """
    Check if node i was visited.
    """
    return bitset[i >> 3] & (1 << (i & 7))


cdef inline void set_visited(uchar *bitset, int i):
    """
    Mark node i as visited.
    """
    bitset[i >> 3] |= 1 << (i & 7)


cpdef calculate_cluster_sizes(double[:, :] Z, double[:] cs, int n):
    """
    Calculate the size of each cluster. The result is the fourth column of
    the linkage matrix.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix. The fourth column can be empty.
    cs : ndarray
        The array to store the sizes.
    n : ndarray
        The number of observations.
    """
    cdef int i, child_l, child_r

    for i in range(n - 1):
        child_l = <int>Z[i, 0]
        child_r = <int>Z[i, 1]

        if child_l >= n:
            cs[i] += cs[child_l - n]
        else:
            cs[i] += 1

        if child_r >= n:
            cs[i] += cs[child_r - n]
        else:
            cs[i] += 1


def cluster_dist(double[:, :] Z, int[:] T, double cutoff, int n):
    """
    Form flat clusters by distance criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    cutoff : double
        Clusters are formed when distances are less than or equal to `cutoff`.
    n : int
        The number of observations.
    """
    cdef double[:] max_dists = np.ndarray(n, dtype=np.double)
    get_max_dist_for_each_cluster(Z, max_dists, n)
    cluster_monocrit(Z, max_dists, T, cutoff, n)


def cluster_in(double[:, :] Z, double[:, :] R, int[:] T, double cutoff, int n):
    """
    Form flat clusters by inconsistent criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    R : ndarray
        The inconsistent matrix.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    cutoff : double
        Clusters are formed when the inconsistent values are less than or
        or equal to `cutoff`.
    n : int
        The number of observations.
    """
    cdef double[:] max_inconsists = np.ndarray(n, dtype=np.double)
    get_max_Rfield_for_each_cluster(Z, R, max_inconsists, n, 3)
    cluster_monocrit(Z, max_inconsists, T, cutoff, n)


def cluster_maxclust_dist(double[:, :] Z, int[:] T, int n, int mc):
    """
    Form flat clusters by maxclust criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    n : int
        The number of observations.
    mc : int
        The maximum number of clusters.
    """
    cdef double[:] max_dists = np.ndarray(n, dtype=np.double)
    get_max_dist_for_each_cluster(Z, max_dists, n)
    # should use an O(n) algorithm
    cluster_maxclust_monocrit(Z, max_dists, T, n, mc)


def cluster_maxclust_monocrit(double[:, :] Z, double[:] MC, int[:] T,
                              int n, int max_nc):
    """
    Form flat clusters by maxclust_monocrit criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    MC : ndarray
        The monotonic criterion array.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    n : int
        The number of observations.
    max_nc : int
        The maximum number of clusters.
    """
    cdef int i, k, i_lc, i_rc, root, nc, lower_idx, upper_idx
    cdef double thresh
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError

    max_illegal = -1
    min_legal = MC[n - 2]

    lower_idx = 0
    upper_idx = n - 1
    while upper_idx - lower_idx > 1:
        i = (lower_idx + upper_idx) >> 1
        thresh = MC[i]

        memset(visited, 0, visited_size)
        nc = 0
        k = 0
        curr_node[0] = 2 * n - 2

        while k >= 0:
            root = curr_node[k] - n
            i_lc = <int>Z[root, 0]
            i_rc = <int>Z[root, 1]

            if MC[root] <= thresh:  # this subtree forms a cluster
                nc += 1
                if nc > max_nc:  # illegal
                    break
                k -= 1
                set_visited(visited, i_lc)
                set_visited(visited, i_rc)
                continue

            if not is_visited(visited, i_lc):
                set_visited(visited, i_lc)
                if i_lc >= n:
                    k += 1
                    curr_node[k] = i_lc
                    continue
                else:  # singleton cluster
                    nc += 1
                    if nc > max_nc:
                        break

            if not is_visited(visited, i_rc):
                set_visited(visited, i_rc)
                if i_rc >= n:
                    k += 1
                    curr_node[k] = i_rc
                    continue
                else:  # singleton cluster
                    nc += 1
                    if nc > max_nc:
                        break

            k -= 1

        if nc > max_nc:
            lower_idx = i
        else:
            upper_idx = i

    PyMem_Free(visited)
    cluster_monocrit(Z, MC, T, MC[upper_idx], n)


cpdef cluster_monocrit(double[:, :] Z, double[:] MC, int[:] T,
                       double cutoff, int n):
    """
    Form flat clusters by monocrit criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    MC : ndarray
        The monotonic criterion array.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    cutoff : double
        Clusters are formed when the MC values are less than or equal to
        `cutoff`.
    n : int
        The number of observations.
    """
    cdef int k, i_lc, i_rc, root, n_cluster = 0, cluster_leader = -1
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError
    memset(visited, 0, visited_size)

    k = 0
    curr_node[0] = 2 * n - 2
    while k >= 0:
        root = curr_node[k] - n
        i_lc = <int>Z[root, 0]
        i_rc = <int>Z[root, 1]

        if cluster_leader == -1 and MC[root] <= cutoff:  # found a cluster
            cluster_leader = root
            n_cluster += 1

        if i_lc >= n and not is_visited(visited, i_lc):
            set_visited(visited, i_lc)
            k += 1
            curr_node[k] = i_lc
            continue

        if i_rc >= n and not is_visited(visited, i_rc):
            set_visited(visited, i_rc)
            k += 1
            curr_node[k] = i_rc
            continue

        if i_lc < n:
            if cluster_leader == -1:  # singleton cluster
                n_cluster += 1
            T[i_lc] = n_cluster

        if i_rc < n:
            if cluster_leader == -1:  # singleton cluster
                n_cluster += 1
            T[i_rc] = n_cluster

        if cluster_leader == root:  # back to the leader
            cluster_leader = -1
        k -= 1

    PyMem_Free(visited)


def cophenetic_distances(double[:, :] Z, double[:] d, int n):
    """
    Calculate the cophenetic distances between each observation

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    d : ndarray
        The condensed matrix to store the cophenetic distances.
    n : int
        The number of observations.
    """
    cdef int i, j, k, root, i_lc, i_rc, n_lc, n_rc, right_start
    cdef double dist
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)
    cdef int[:] members = np.ndarray(n, dtype=np.int32)
    cdef int[:] left_start = np.ndarray(n, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError
    memset(visited, 0, visited_size)

    k = 0
    curr_node[0] = 2 * n - 2
    left_start[0] = 0
    while k >= 0:
        root = curr_node[k] - n
        i_lc = <int>Z[root, 0]
        i_rc = <int>Z[root, 1]

        if i_lc >= n:  # left child is not a leaf
            n_lc = <int>Z[i_lc - n, 3]

            if not is_visited(visited, i_lc):
                set_visited(visited, i_lc)
                k += 1
                curr_node[k] = i_lc
                left_start[k] = left_start[k - 1]
                continue  # visit the left subtree
        else:
            n_lc = 1
            members[left_start[k]] = i_lc

        if i_rc >= n:  # right child is not a leaf
            n_rc = <int>Z[i_rc - n, 3]

            if not is_visited(visited, i_rc):
                set_visited(visited, i_rc)
                k += 1
                curr_node[k] = i_rc
                left_start[k] = left_start[k - 1] + n_lc
                continue  # visit the right subtree
        else:
            n_rc = 1
            members[left_start[k] + n_lc] = i_rc

        # back to the root of current subtree
        dist = Z[root, 2]
        right_start = left_start[k] + n_lc
        for i in range(left_start[k], right_start):
            for j in range(right_start, right_start + n_rc):
                d[condensed_index(n, members[i], members[j])] = dist 

        k -= 1  # back to parent node

    PyMem_Free(visited)


cdef from_pointer_representation(double[:, :] Z, double[:] Lambda, int[:] Pi,
                                 int n):
    """
    Generate a linkage matrix from its pointer representation.

    Parameters
    ----------
    Z : ndarray
        An array to store the linkage matrix.
    Lambda : ndarray
        The :math:`\\Lambda` array of the pointer representation.
    Pi : ndarray
        The :math:`\\Pi` array of the pointer representation.
    n : int
        The number of observations.
    """
    cdef int i, current_leaf, pi
    cdef np.intp_t[:] sorted_idx = np.argsort(Lambda)
    cdef int[:] node_ids = np.ndarray(n, dtype=np.int32)

    for i in range(n):
        node_ids[i] = i

    for i in range(n - 1):
        current_leaf = sorted_idx[i]
        pi = Pi[current_leaf]
        if node_ids[current_leaf] < node_ids[pi]:
            Z[i, 0] = node_ids[current_leaf]
            Z[i, 1] = node_ids[pi]
        else:
            Z[i, 0] = node_ids[pi]
            Z[i, 1] = node_ids[current_leaf]
        Z[i, 2] = Lambda[current_leaf]
        node_ids[pi] = n + i

    Z[:, 3] = 0
    calculate_cluster_sizes(Z, Z[:, 3], n)


cpdef get_max_Rfield_for_each_cluster(double[:, :] Z, double[:, :] R,
                                      double[:] max_rfs, int n, int rf):
    """
    Get the maximum statistic for each non-singleton cluster. For the i'th
    non-singleton cluster, max_rfs[i] = max{R[j, rf] j is a descendent of i}.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    R : ndarray
        The R matrix.
    max_rfs : ndarray
        The array to store the result.
    n : int
        The number of observations.
    rf : int
        Indicate which column of `R` is used.
    """
    cdef int k, i_lc, i_rc, root
    cdef double max_rf, max_l, max_r
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError
    memset(visited, 0, visited_size)

    k = 0
    curr_node[0] = 2 * n - 2
    while k >= 0:
        root = curr_node[k] - n
        i_lc = <int>Z[root, 0]
        i_rc = <int>Z[root, 1]

        if i_lc >= n and not is_visited(visited, i_lc):
            set_visited(visited, i_lc)
            k += 1
            curr_node[k] = i_lc
            continue

        if i_rc >= n and not is_visited(visited, i_rc):
            set_visited(visited, i_rc)
            k += 1
            curr_node[k] = i_rc
            continue

        max_rf = R[root, rf]
        if i_lc >= n:
            max_l = max_rfs[i_lc - n]
            if max_l > max_rf:
                max_rf = max_l
        if i_rc >= n:
            max_r = max_rfs[i_rc - n]
            if max_r > max_rf:
                max_rf = max_r
        max_rfs[root] = max_rf

        k -= 1

    PyMem_Free(visited)


cpdef get_max_dist_for_each_cluster(double[:, :] Z, double[:] MD, int n):
    """
    Get the maximum inconsistency coefficient for each non-singleton cluster.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    MD : ndarray
        The array to store the result.
    n : int
        The number of observations.
    """
    cdef int k, i_lc, i_rc, root
    cdef double max_dist, max_l, max_r
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError
    memset(visited, 0, visited_size)

    k = 0
    curr_node[0] = 2 * n - 2
    while k >= 0:
        root = curr_node[k] - n
        i_lc = <int>Z[root, 0]
        i_rc = <int>Z[root, 1]

        if i_lc >= n and not is_visited(visited, i_lc):
            set_visited(visited, i_lc)
            k += 1
            curr_node[k] = i_lc
            continue

        if i_rc >= n and not is_visited(visited, i_rc):
            set_visited(visited, i_rc)
            k += 1
            curr_node[k] = i_rc
            continue

        max_dist = Z[root, 2]
        if i_lc >= n:
            max_l = MD[i_lc - n]
            if max_l > max_dist:
                max_dist = max_l
        if i_rc >= n:
            max_r = MD[i_rc - n]
            if max_r > max_dist:
                max_dist = max_r
        MD[root] = max_dist

        k -= 1

    PyMem_Free(visited)


def inconsistent(double[:, :] Z, double[:, :] R, int n, int d):
    """
    Calculate the inconsistency statistics.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    R : ndarray
        A (n - 1) x 5 matrix to store the result. The inconsistency statistics
        `R[i]` are calculated over `d` levels below cluster i. `R[i, 0]` is the
        mean of distances. `R[i, 1]` is the standard deviation of distances.
        `R[i, 2]` is the number of clusters included. `R[i, 3]` is the
        inconsistency coefficient.

        .. math:: \\frac{\\mathtt{Z[i,2]}-\\mathtt{R[i,0]}} {R[i,1]}

    n : int
        The number of observations.
    d : int
        The number of levels included in calculation below a node.
    """
    cdef int i, k, i_lc, i_rc, root, level_count
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)
    cdef double level_sum, level_std_sum, level_std, dist

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError

    for i in range(n - 1):
        k = 0
        level_count = 0
        level_sum = 0
        level_std_sum = 0
        memset(visited, 0, visited_size)
        curr_node[0] = i

        while k >= 0:
            root = curr_node[k]

            if k < d - 1:
                i_lc = <int>Z[root, 0]
                if i_lc >= n and not is_visited(visited, i_lc):
                    set_visited(visited, i_lc)
                    k += 1
                    curr_node[k] = i_lc - n
                    continue

                i_rc = <int>Z[root, 1]
                if i_rc >= n and not is_visited(visited, i_rc):
                    set_visited(visited, i_rc)
                    k += 1
                    curr_node[k] = i_rc - n
                    continue

            dist = Z[root, 2]
            level_count += 1
            level_sum += dist
            level_std_sum += dist * dist
            k -= 1

        R[i, 0] = level_sum / level_count
        R[i, 2] = level_count
        if level_count < 2:
            level_std = (level_std_sum - (level_sum * level_sum)) / level_count
        else:
            level_std = ((level_std_sum -
                         ((level_sum * level_sum) / level_count)) /
                         (level_count - 1))
        if level_std > 0:
            level_std = sqrt(level_std)
            R[i, 1] = level_std
            R[i, 3] = (Z[i, 2] - R[i, 0]) / level_std
        else:
            R[i, 1] = 0

    PyMem_Free(visited)


def leaders(double[:, :] Z, int[:] T, int[:] L, int[:] M, int nc, int n):
    """
    Find the leader (root) of each flat cluster.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    T : ndarray
        The flat clusters assignment returned by `fcluster` or `fclusterdata`.
    L : ndarray
        `L` and `M` store the result. The leader of flat cluster `L[i]` is
        node `M[i]`.
    M : ndarray
        `L` and `M` store the result. The leader of flat cluster `L[i]` is
        node `M[i]`.
    nc : int
        The number of flat clusters.
    n : int
        The number of observations.

    Returns
    -------
    err_node : int
        Found that `T` is invalid when examining node `err_node`.
        `-1` indicates success.
    """
    cdef int k, i_lc, i_rc, root, cid_lc, cid_rc, leader_idx, result = -1
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)
    cdef int[:] cluster_ids = np.ndarray(n * 2 - 1, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError
    memset(visited, 0, visited_size)

    cluster_ids[:n] = T[:]
    cluster_ids[n:] = -1
    k = 0
    curr_node[0] = 2 * n - 2
    leader_idx = 0
    while k >= 0:
        root = curr_node[k] - n
        i_lc = <int>Z[root, 0]
        i_rc = <int>Z[root, 1]

        if i_lc >= n and not is_visited(visited, i_lc):
            set_visited(visited, i_lc)
            k += 1
            curr_node[k] = i_lc
            continue

        if i_rc >= n and not is_visited(visited, i_rc):
            set_visited(visited, i_rc)
            k += 1
            curr_node[k] = i_rc
            continue

        cid_lc = cluster_ids[i_lc]
        cid_rc = cluster_ids[i_rc]

        if cid_lc == cid_rc:  # left and right children in the same cluster
            cluster_ids[root + n] = cid_lc
        else:  # left and right children are both leaders
            if cid_lc != -1:
                if leader_idx >= nc:
                    result = root + n
                    break
                L[leader_idx] = i_lc
                M[leader_idx] = cid_lc
                leader_idx += 1

            if cid_rc != -1:
                if leader_idx >= nc:
                    result = root + n
                    break
                L[leader_idx] = i_rc
                M[leader_idx] = cid_rc
                leader_idx += 1

            cluster_ids[root + n] = -1

        k -= 1

    if result == -1:
        i_lc = <int>Z[n - 2, 0]
        i_rc = <int>Z[n - 2, 1]
        cid_lc = cluster_ids[i_lc]
        cid_rc = cluster_ids[i_rc]
        if cid_lc == cid_rc and cid_lc != -1:
            if leader_idx >= nc:
                result = 2 * n - 2
            else:
                L[leader_idx] = 2 * n - 2
                M[leader_idx] = cid_lc

    PyMem_Free(visited)
    return result  # -1 means success here


def linkage(double[:] dists, np.ndarray[double, ndim=2] Z, int n, int method):
    """
    Perform hierarchy clustering.

    Parameters
    ----------
    dists : ndarray
        A condensed matrix stores the pairwise distances of the observations.
    Z : ndarray
        A (n - 1) x 4 matrix to store the result (i.e. the linkage matrix).
    n : int
        The number of observations.
    method : int
        The linkage method. 0: single 1: complete 2: average 3: centroid
        4: median 5: ward 6: weighted
    """
    cdef int i, j, k, x, y, i_start, nx, ny, ni, id_x, id_y, id_i
    cdef double current_min
    # inter-cluster dists
    cdef double[:] D = np.ndarray(n * (n - 1) / 2, dtype=np.double)
    # map the indices to node ids
    cdef int[:] id_map = np.ndarray(n, dtype=np.int32)
    cdef linkage_distance_update new_dist

    new_dist = linkage_methods[method]

    D[:] = dists
    for i in range(n):
        id_map[i] = i

    for k in range(n - 1):
        # find two closest clusters x, y (x < y)
        current_min = NPY_INFINITYF
        for i in range(n - 1):
            if id_map[i] == -1:
                continue

            i_start = condensed_index(n, i, i + 1)
            for j in range(n - i - 1):
                if D[i_start + j] < current_min:
                    current_min = D[i_start + j]
                    x = i
                    y = i + j + 1

        id_x = id_map[x]
        id_y = id_map[y]

        # get the original numbers of points in clusters x and y
        nx = 1 if id_x < n else <int>Z[id_x - n, 3]
        ny = 1 if id_y < n else <int>Z[id_y - n, 3]

        # record the new node
        if id_x < id_y:
            Z[k] = (id_x, id_y, current_min, nx + ny)
        else:
            Z[k] = (id_y, id_x, current_min, nx + ny)
        id_map[x] = -1  # cluster x will be dropped
        id_map[y] = n + k  # cluster y will be replaced with the new cluster

        # update the distance matrix
        for i in range(n):
            id_i = id_map[i]
            if id_i == -1 or id_i == n + k:
                continue

            ni = 1 if id_i < n else <int>Z[id_i - n, 3]
            D[condensed_index(n, i, y)] = new_dist(
                D[condensed_index(n, i, x)],
                D[condensed_index(n, i, y)],
                current_min, nx, ny, ni)
            if i < x:
                D[condensed_index(n, i, x)] = NPY_INFINITYF


def prelist(double[:, :] Z, int[:] members, int n):
    """
    Perform a pre-order traversal on the linkage tree and get a list of ids
    of the leaves.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    members : ndarray
        The array to store the result.
    n : int
        The number of observations.
    """
    cdef int k, i_lc, i_rc, root, mem_idx
    cdef int[:] curr_node = np.ndarray(n, dtype=np.int32)

    cdef int visited_size = (((n * 2) - 1) >> 3) + 1
    cdef uchar *visited = <uchar *>PyMem_Malloc(visited_size)
    if not visited:
        raise MemoryError
    memset(visited, 0, visited_size)

    mem_idx = 0
    k = 0
    curr_node[0] = 2 * n - 2
    while k >= 0:
        root = curr_node[k] - n

        i_lc = <int>Z[root, 0]
        if not is_visited(visited, i_lc):
            set_visited(visited, i_lc)
            if i_lc >= n:
                k += 1
                curr_node[k] = i_lc
                continue
            else:
                members[mem_idx] = i_lc
                mem_idx += 1

        i_rc = <int>Z[root, 1]
        if not is_visited(visited, i_rc):
            set_visited(visited, i_rc)
            if i_rc >= n:
                k += 1
                curr_node[k] = i_rc
                continue
            else:
                members[mem_idx] = i_rc
                mem_idx += 1

        k -= 1

    PyMem_Free(visited)


def slink(double[:] dists, double[:, :] Z, int n):
    """
    The SLINK algorithm. Single linkage in O(n^2) time complexity.

    Parameters
    ----------
    dists : ndarray
        A condensed matrix stores the pairwise distances of the observations.
    Z : ndarray
        A (n - 1) x 4 matrix to store the result (i.e. the linkage matrix).
    n : int
        The number of observations.

    References
    ----------
    R. Sibson, "SLINK: An optimally efficient algorithm for the single-link
    cluster method", The Computer Journal 1973 16: 30-34.
    """
    cdef int i, j
    cdef double[:] M = np.ndarray(n, dtype=np.double)
    cdef double[:] Lambda = np.ndarray(n, dtype=np.double)
    cdef int[:] Pi = np.ndarray(n, dtype=np.int32)

    Pi[0] = 0
    Lambda[0] = NPY_INFINITYF
    for i in range(1, n):
        Pi[i] = i
        Lambda[i] = NPY_INFINITYF

        for j in range(i):
            M[j] = dists[condensed_index(n, i, j)]

        for j in range(i):
            if Lambda[j] >= M[j]:
                M[Pi[j]] = min(M[Pi[j]], Lambda[j])
                Lambda[j] = M[j]
                Pi[j] = i
            else:
                M[Pi[j]] = min(M[Pi[j]], M[j])

        for j in range(i):
            if Lambda[j] >= Lambda[Pi[j]]:
                Pi[j] = i

    from_pointer_representation(Z, Lambda, Pi, n)
