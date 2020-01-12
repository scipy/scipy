# Copyright Anne M. Archibald 2008
# Additional contributions by Patrick Varilly and Sturla Molden 2012
# Revision by Sturla Molden 2015
# Balanced kd-tree construction written by Jake Vanderplas for scikit-learn
# Released under the scipy license

# distutils: language = c++

from __future__ import absolute_import

import numpy as np
import scipy.sparse

cimport numpy as np
from numpy.math cimport INFINITY

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memset, memcpy
from libcpp.vector cimport vector
from libcpp.algorithm cimport sort

cimport cython

from multiprocessing import cpu_count
import threading

cdef extern from "<limits.h>":
    long LONG_MAX

cdef int number_of_processors = cpu_count()

__all__ = ['cKDTree']

cdef extern from *:
    int NPY_LIKELY(int)
    int NPY_UNLIKELY(int)


# C++ implementations
# ===================

cdef extern from "ckdtree_decl.h":
    int ckdtree_isinf(np.float64_t x) nogil

    struct ckdtreenode:
        np.intp_t split_dim
        np.intp_t children
        np.float64_t split
        np.intp_t start_idx
        np.intp_t end_idx
        ckdtreenode *less
        ckdtreenode *greater
        np.intp_t _less
        np.intp_t _greater

    struct ckdtree:
        vector[ckdtreenode]  *tree_buffer
        ckdtreenode   *ctree
        np.float64_t   *raw_data
        np.intp_t      n
        np.intp_t      m
        np.intp_t      leafsize
        np.float64_t   *raw_maxes
        np.float64_t   *raw_mins
        np.intp_t      *raw_indices
        np.float64_t   *raw_boxsize_data
        np.intp_t size

    # External build and query methods in C++.
    
    int build_ckdtree(ckdtree *self,
                         np.intp_t start_idx,
                         np.intp_t end_idx,
                         np.float64_t *maxes,
                         np.float64_t *mins,
                         int _median,
                         int _compact) nogil except +

    int build_weights(ckdtree *self,
                         np.float64_t *node_weights,
                         np.float64_t *weights) nogil except +

    int query_knn(const ckdtree *self,
                     np.float64_t *dd,
                     np.intp_t    *ii,
                     const np.float64_t *xx,
                     const np.intp_t    n,
                     const np.intp_t    *k,
                     const np.intp_t    nk,
                     const np.intp_t    kmax,
                     const np.float64_t eps,
                     const np.float64_t p,
                     const np.float64_t distance_upper_bound) nogil except +

    int query_pairs(const ckdtree *self,
                       const np.float64_t r,
                       const np.float64_t p,
                       const np.float64_t eps,
                       vector[ordered_pair] *results) nogil except +

    int count_neighbors_unweighted(const ckdtree *self,
                           const ckdtree *other,
                           np.intp_t     n_queries,
                           np.float64_t  *real_r,
                           np.intp_t     *results,
                           const np.float64_t p,
                           int cumulative) nogil except +

    int count_neighbors_weighted(const ckdtree *self,
                           const ckdtree *other,
                           np.float64_t  *self_weights,
                           np.float64_t  *other_weights,
                           np.float64_t  *self_node_weights,
                           np.float64_t  *other_node_weights,
                           np.intp_t     n_queries,
                           np.float64_t  *real_r,
                           np.float64_t     *results,
                           const np.float64_t p,
                           int cumulative) nogil except +

    int query_ball_point(const ckdtree *self,
                            const np.float64_t *x,
                            const np.float64_t *r,
                            const np.float64_t p,
                            const np.float64_t eps,
                            const np.intp_t n_queries,
                            vector[np.intp_t] **results,
                            const int return_length) nogil except +

    int query_ball_tree(const ckdtree *self,
                           const ckdtree *other,
                           const np.float64_t r,
                           const np.float64_t p,
                           const np.float64_t eps,
                           vector[np.intp_t] **results) nogil except +

    int sparse_distance_matrix(const ckdtree *self,
                                  const ckdtree *other,
                                  const np.float64_t p,
                                  const np.float64_t max_distance,
                                  vector[coo_entry] *results) nogil except +


# C++ helper functions
# ====================

cdef extern from "coo_entries.h":

    struct coo_entry:
        np.intp_t i
        np.intp_t j
        np.float64_t v

cdef extern from "ordered_pair.h":

    struct ordered_pair:
        np.intp_t i
        np.intp_t j

# coo_entry wrapper
# =================

cdef class coo_entries:

    cdef:
        readonly object __array_interface__
        vector[coo_entry] *buf

    def __cinit__(coo_entries self):
        self.buf = NULL

    def __init__(coo_entries self):
        self.buf = new vector[coo_entry]()

    def __dealloc__(coo_entries self):
        if self.buf != NULL:
            del self.buf

    # The methods ndarray, dict, coo_matrix, and dok_matrix must only
    # be called after the buffer is filled with coo_entry data. This
    # is because std::vector can reallocate its internal buffer when
    # push_back is called.

    def ndarray(coo_entries self):
        cdef:
            coo_entry *pr
            np.uintp_t uintptr
            np.intp_t n
        _dtype = [('i',np.intp),('j',np.intp),('v',np.float64)]
        res_dtype = np.dtype(_dtype, align = True)
        n = <np.intp_t> self.buf.size()
        if NPY_LIKELY(n > 0):
            pr = self.buf.data()
            uintptr = <np.uintp_t> (<void*> pr)
            dtype = np.dtype(np.uint8)
            self.__array_interface__ = dict(
                data = (uintptr, False),
                descr = dtype.descr,
                shape = (n*sizeof(coo_entry),),
                strides = (dtype.itemsize,),
                typestr = dtype.str,
                version = 3,
            )
            return np.asarray(self).view(dtype=res_dtype)
        else:
            return np.empty(shape=(0,), dtype=res_dtype)

    def dict(coo_entries self):
        cdef:
            np.intp_t i, j, k, n
            np.float64_t v
            coo_entry *pr
            dict res_dict
        n = <np.intp_t> self.buf.size()
        if NPY_LIKELY(n > 0):
            pr = self.buf.data()
            res_dict = dict()
            for k in range(n):
                i = pr[k].i
                j = pr[k].j
                v = pr[k].v
                res_dict[(i,j)] = v
            return res_dict
        else:
            return {}

    def coo_matrix(coo_entries self, m, n):
        res_arr = self.ndarray()
        return scipy.sparse.coo_matrix(
                       (res_arr['v'], (res_arr['i'], res_arr['j'])),
                                       shape=(m, n))

    def dok_matrix(coo_entries self, m, n):
        return self.coo_matrix(m,n).todok()


# ordered_pair wrapper
# ====================

cdef class ordered_pairs:

    cdef:
        readonly object __array_interface__
        vector[ordered_pair] *buf

    def __cinit__(ordered_pairs self):
        self.buf = NULL

    def __init__(ordered_pairs self):
        self.buf = new vector[ordered_pair]()

    def __dealloc__(ordered_pairs self):
        if self.buf != NULL:
            del self.buf

    # The methods ndarray and set must only be called after the buffer
    # is filled with ordered_pair data.

    def ndarray(ordered_pairs self):
        cdef:
            ordered_pair *pr
            np.uintp_t uintptr
            np.intp_t n
        n = <np.intp_t> self.buf.size()
        if NPY_LIKELY(n > 0):
            pr = self.buf.data()
            uintptr = <np.uintp_t> (<void*> pr)
            dtype = np.dtype(np.intp)
            self.__array_interface__ = dict(
                data = (uintptr, False),
                descr = dtype.descr,
                shape = (n,2),
                strides = (2*dtype.itemsize,dtype.itemsize),
                typestr = dtype.str,
                version = 3,
            )
            return np.asarray(self)
        else:
            return np.empty(shape=(0,2), dtype=np.intp)

    def set(ordered_pairs self):
        cdef:
            ordered_pair *pair
            np.intp_t i, n
            set results
        results = set()
        pair = self.buf.data()
        n = <np.intp_t> self.buf.size()
        if sizeof(long) < sizeof(np.intp_t):
            # Needed for Python 2.x on Win64
            for i in range(n):
                results.add((int(pair.i), int(pair.j)))
                pair += 1
        else:
            # other platforms
            for i in range(n):
                results.add((pair.i, pair.j))
                pair += 1
        return results



# Tree structure exposed to Python
# ================================

cdef class cKDTreeNode:
    """
    class cKDTreeNode

    This class exposes a Python view of a node in the cKDTree object.

    All attributes are read-only.

    Attributes
    ----------
    level : int
        The depth of the node. 0 is the level of the root node.
    split_dim : int
        The dimension along which this node is split. If this value is -1
        the node is a leafnode in the kd-tree. Leafnodes are not split further
        and scanned by brute force.
    split : float
        The value used to separate split this node. Points with value >= split
        in the split_dim dimension are sorted to the 'greater' subnode
        whereas those with value < split are sorted to the 'lesser' subnode.
    children : int
        The number of data points sorted to this node.
    data_points : ndarray of float64
        An array with the data points sorted to this node.
    indices : ndarray of intp
        An array with the indices of the data points sorted to this node. The
        indices refer to the position in the data set used to construct the
        kd-tree.
    lesser : cKDTreeNode or None
        Subnode with the 'lesser' data points. This attribute is None for
        leafnodes.
    greater : cKDTreeNode or None
        Subnode with the 'greater' data points. This attribute is None for
        leafnodes.

    """
    cdef:
        readonly np.intp_t    level
        readonly np.intp_t    split_dim
        readonly np.intp_t    children
        readonly np.float64_t split
        ckdtreenode           *_node
        np.ndarray            _data
        np.ndarray            _indices

    cdef void _setup(cKDTreeNode self):
        self.split_dim = self._node.split_dim
        self.children = self._node.children
        self.split = self._node.split

    property data_points:
        def __get__(cKDTreeNode self):
            return self._data[self.indices,:]

    property indices:
        def __get__(cKDTreeNode self):
            cdef np.intp_t i, start, stop
            if self.split_dim == -1:
                start = self._node.start_idx
                stop = self._node.end_idx
                return self._indices[start:stop]
            else:
                return np.hstack([self.lesser.indices,
                           self.greater.indices])

    property lesser:
        def __get__(cKDTreeNode self):
            if self.split_dim == -1:
                return None
            else:
                n = cKDTreeNode()
                n._node = self._node.less
                n._data = self._data
                n._indices = self._indices
                n.level = self.level + 1
                n._setup()
                return n

    property greater:
        def __get__(cKDTreeNode self):
            if self.split_dim == -1:
                return None
            else:
                n = cKDTreeNode()
                n._node = self._node.greater
                n._data = self._data
                n._indices = self._indices
                n.level = self.level + 1
                n._setup()
                return n


# Main cKDTree class
# ==================

cdef class cKDTree:
    """
    cKDTree(data, leafsize=16, compact_nodes=True, copy_data=False,
            balanced_tree=True, boxsize=None)

    kd-tree for quick nearest-neighbor lookup

    This class provides an index into a set of k-dimensional points
    which can be used to rapidly look up the nearest neighbors of any
    point.

    The algorithm used is described in Maneewongvatana and Mount 1999.
    The general idea is that the kd-tree is a binary trie, each of whose
    nodes represents an axis-aligned hyperrectangle. Each node specifies
    an axis and splits the set of points based on whether their coordinate
    along that axis is greater than or less than a particular value.

    During construction, the axis and splitting point are chosen by the
    "sliding midpoint" rule, which ensures that the cells do not all
    become long and thin.

    The tree can be queried for the r closest neighbors of any given point
    (optionally returning only those within some maximum distance of the
    point). It can also be queried, with a substantial gain in efficiency,
    for the r approximate closest neighbors.

    For large dimensions (20 is already large) do not expect this to run
    significantly faster than brute force. High-dimensional nearest-neighbor
    queries are a substantial open problem in computer science.

    Parameters
    ----------
    data : array_like, shape (n,m)
        The n data points of dimension m to be indexed. This array is
        not copied unless this is necessary to produce a contiguous
        array of doubles, and so modifying this data will result in
        bogus results. The data are also copied if the kd-tree is built
        with copy_data=True.
    leafsize : positive int, optional
        The number of points at which the algorithm switches over to
        brute-force. Default: 16.
    compact_nodes : bool, optional
        If True, the kd-tree is built to shrink the hyperrectangles to
        the actual data range. This usually gives a more compact tree that
        is robust against degenerated input data and gives faster queries
        at the expense of longer build time. Default: True.
    copy_data : bool, optional
        If True the data is always copied to protect the kd-tree against
        data corruption. Default: False.
    balanced_tree : bool, optional
        If True, the median is used to split the hyperrectangles instead of
        the midpoint. This usually gives a more compact tree and
        faster queries at the expense of longer build time. Default: True.
    boxsize : array_like or scalar, optional
        Apply a m-d toroidal topology to the KDTree.. The topology is generated
        by :math:`x_i + n_i L_i` where :math:`n_i` are integers and :math:`L_i`
        is the boxsize along i-th dimension. The input data shall be wrapped
        into :math:`[0, L_i)`. A ValueError is raised if any of the data is
        outside of this bound.

    Attributes
    ----------
    data : ndarray, shape (n,m)
        The n data points of dimension m to be indexed. This array is
        not copied unless this is necessary to produce a contiguous
        array of doubles. The data are also copied if the kd-tree is built
        with `copy_data=True`.
    leafsize : positive int
        The number of points at which the algorithm switches over to
        brute-force.
    m : int
        The dimension of a single data-point.
    n : int
        The number of data points.
    maxes : ndarray, shape (m,)
        The maximum value in each dimension of the n data points.
    mins : ndarray, shape (m,)
        The minimum value in each dimension of the n data points.
    tree : object, class cKDTreeNode
        This class exposes a Python view of the root node in the cKDTree object.
    size : int
        The number of nodes in the tree.

    See Also
    --------
    KDTree : Implementation of `cKDTree` in pure Python

    """
    cdef:
        ckdtree * cself
        readonly cKDTreeNode     tree
        readonly np.ndarray      data
        readonly np.ndarray      maxes
        readonly np.ndarray      mins
        readonly np.ndarray      indices
        readonly object          boxsize
        np.ndarray               boxsize_data

    property n:
        def __get__(self): return self.cself.n
    property m:
        def __get__(self): return self.cself.m
    property leafsize:
        def __get__(self): return self.cself.leafsize
    property size:
        def __get__(self): return self.cself.size

    def __cinit__(cKDTree self):
        self.cself = <ckdtree * > PyMem_Malloc(sizeof(ckdtree))
        self.cself.tree_buffer = NULL

    def __init__(cKDTree self, data, np.intp_t leafsize=16, compact_nodes=True,
            copy_data=False, balanced_tree=True, boxsize=None):
        
        cdef: 
            np.float64_t [::1] tmpmaxes, tmpmins
            np.float64_t *ptmpmaxes
            np.float64_t *ptmpmins
            ckdtree *cself = self.cself
            int compact, median

        data = np.array(data, order='C', copy=copy_data, dtype=np.float64)

        if data.ndim != 2:
            raise ValueError("data must be 2 dimensions")

        self.data = data
        cself.n = data.shape[0]
        cself.m = data.shape[1]
        cself.leafsize = leafsize

        if leafsize<1:
            raise ValueError("leafsize must be at least 1")

        if boxsize is None:
            self.boxsize = None
            self.boxsize_data = None
        else:
            self.boxsize_data = np.empty(2 * self.m, dtype=np.float64)
            boxsize = np.float64(np.broadcast_to(boxsize, self.m))
            self.boxsize_data[:self.m] = boxsize
            self.boxsize_data[self.m:] = 0.5 * boxsize

            self.boxsize = boxsize
            periodic_mask = self.boxsize > 0
            if ((self.data >= self.boxsize[None, :])[:, periodic_mask]).any():
                raise ValueError("Some input data are greater than the size of the periodic box.")
            if ((self.data < 0)[:, periodic_mask]).any():
                raise ValueError("Negative input data are outside of the periodic box.")

        self.maxes = np.ascontiguousarray(
            np.amax(self.data, axis=0) if self.n > 0 else np.zeros(self.m),
            dtype=np.float64)
        self.mins = np.ascontiguousarray(
            np.amin(self.data,axis=0) if self.n > 0 else np.zeros(self.m),
            dtype=np.float64)
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=np.intp))

        self._pre_init()

        compact = 1 if compact_nodes else 0
        median = 1 if balanced_tree else 0

        cself.tree_buffer = new vector[ckdtreenode]()

        tmpmaxes = np.copy(self.maxes)
        tmpmins = np.copy(self.mins)
        
        ptmpmaxes = &tmpmaxes[0]
        ptmpmins = &tmpmins[0]
        with nogil: 
            build_ckdtree(cself, 0, cself.n, ptmpmaxes, ptmpmins, median, compact)

        # set up the tree structure pointers
        self._post_init()

    cdef _pre_init(cKDTree self):
        cself = self.cself

        # finalize the pointers from array attributes

        cself.raw_data = <np.float64_t*> np.PyArray_DATA(self.data)
        cself.raw_maxes = <np.float64_t*> np.PyArray_DATA(self.maxes)
        cself.raw_mins = <np.float64_t*> np.PyArray_DATA(self.mins)
        cself.raw_indices = <np.intp_t*> np.PyArray_DATA(self.indices)

        if self.boxsize_data is not None:
            cself.raw_boxsize_data = <np.float64_t*>np.PyArray_DATA(self.boxsize_data)
        else:
            cself.raw_boxsize_data = NULL

    cdef _post_init(cKDTree self):
        cself = self.cself
        # finalize the tree points, this calls _post_init_traverse

        cself.ctree = cself.tree_buffer.data()

        # set the size attribute after tree_buffer is built
        cself.size = cself.tree_buffer.size()

        self._post_init_traverse(cself.ctree)

        # make the tree viewable from Python
        self.tree = cKDTreeNode()
        self.tree._node = cself.ctree
        self.tree._data = self.data
        self.tree._indices = self.indices
        self.tree.level = 0
        self.tree._setup()

    cdef _post_init_traverse(cKDTree self, ckdtreenode *node):
        cself = self.cself
        # recurse the tree and re-initialize
        # "less" and "greater" fields
        if node.split_dim == -1:
            # leafnode
            node.less = NULL
            node.greater = NULL
        else:
            node.less = cself.ctree + node._less
            node.greater = cself.ctree + node._greater
            self._post_init_traverse(node.less)
            self._post_init_traverse(node.greater)

    def __dealloc__(cKDTree self):
        cself = self.cself
        if cself.tree_buffer != NULL:
            del cself.tree_buffer
        PyMem_Free(cself)

    # -----
    # query
    # -----

    @cython.boundscheck(False)
    def query(cKDTree self, object x, object k=1, np.float64_t eps=0,
              np.float64_t p=2, np.float64_t distance_upper_bound=INFINITY,
              np.intp_t n_jobs=1):
        """
        query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf, n_jobs=1)

        Query the kd-tree for nearest neighbors

        Parameters
        ----------
        x : array_like, last dimension self.m
            An array of points to query.
        k : list of integer or integer
            The list of k-th nearest neighbors to return. If k is an
            integer it is treated as a list of [1, ... k] (range(1, k+1)).
            Note that the counting starts from 1.
        eps : non-negative float
            Return approximate nearest neighbors; the k-th returned value
            is guaranteed to be no further than (1+eps) times the
            distance to the real k-th nearest neighbor.
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
            A finite large p may cause a ValueError if overflow can occur.
        distance_upper_bound : nonnegative float
            Return only neighbors within this distance.  This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.
        n_jobs : int, optional
            Number of jobs to schedule for parallel processing. If -1 is given
            all processors are used. Default: 1.

        Returns
        -------
        d : array of floats
            The distances to the nearest neighbors.
            If ``x`` has shape ``tuple+(self.m,)``, then ``d`` has shape ``tuple+(k,)``.
            When k == 1, the last dimension of the output is squeezed.
            Missing neighbors are indicated with infinite distances.
        i : ndarray of ints
            The locations of the neighbors in ``self.data``.
            If ``x`` has shape ``tuple+(self.m,)``, then ``i`` has shape ``tuple+(k,)``.
            When k == 1, the last dimension of the output is squeezed.
            Missing neighbors are indicated with ``self.n``.

        Notes
        -----
        If the KD-Tree is periodic, the position ``x`` is wrapped into the
        box.

        When the input k is a list, a query for arange(max(k)) is performed, but
        only columns that store the requested values of k are preserved. This is
        implemented in a manner that reduces memory usage.

        Examples
        --------

        >>> import numpy as np
        >>> from scipy.spatial import cKDTree
        >>> x, y = np.mgrid[0:5, 2:8]
        >>> tree = cKDTree(np.c_[x.ravel(), y.ravel()])

        To query the nearest neighbours and return squeezed result, use

        >>> dd, ii = tree.query([[0, 0], [2.1, 2.9]], k=1)
        >>> print(dd, ii)
        [2.         0.14142136] [ 0 13]

        To query the nearest neighbours and return unsqueezed result, use

        >>> dd, ii = tree.query([[0, 0], [2.1, 2.9]], k=[1])
        >>> print(dd, ii)
        [[2.        ]
         [0.14142136]] [[ 0]
         [13]]

        To query the second nearest neighbours and return unsqueezed result, use

        >>> dd, ii = tree.query([[0, 0], [2.1, 2.9]], k=[2])
        >>> print(dd, ii)
        [[2.23606798]
         [0.90553851]] [[ 6]
         [12]]

        To query the first and second nearest neighbours, use

        >>> dd, ii = tree.query([[0, 0], [2.1, 2.9]], k=2)
        >>> print(dd, ii)
        [[2.         2.23606798]
         [0.14142136 0.90553851]] [[ 0  6]
         [13 12]]

        or, be more specific

        >>> dd, ii = tree.query([[0, 0], [2.1, 2.9]], k=[1, 2])
        >>> print(dd, ii)
        [[2.         2.23606798]
         [0.14142136 0.90553851]] [[ 0  6]
         [13 12]]

        """

        cdef:
            np.intp_t n, i, j
            int overflown
            const np.float64_t [:, ::1] xx

        xshape = np.shape(x)

        if len(xshape) == 0 or xshape[-1] != self.m:
            raise ValueError("x must consist of vectors of length %d but "
                             "has shape %s" % (int(self.m), xshape))

        n = <np.intp_t> np.prod(xshape[:-1])
        xx = np.ascontiguousarray(x, dtype=np.float64).reshape(n, self.m)

        if p < 1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        if len(xshape) == 1:
            single = True
        else:
            single = False

        nearest = False
        if np.isscalar(k):
            if k == 1:
                nearest = True
            k = np.arange(1, k + 1)

        retshape = xshape[:-1]

        # The C++ function touches all dd and ii entries,
        # setting the missing values.

        cdef: 
            np.float64_t [:, ::1] dd = np.empty((n,len(k)),dtype=np.float64)
            np.intp_t [:, ::1] ii = np.empty((n,len(k)),dtype=np.intp)
            np.intp_t [::1] kk = np.array(k, dtype=np.intp)
            np.intp_t kmax = np.max(k)

        # Do the query in an external C++ function.
        def _thread_func(np.intp_t start, np.intp_t stop):
            cdef:
                np.float64_t *pdd = &dd[start,0]
                np.intp_t *pii = &ii[start,0]
                const np.float64_t *pxx = &xx[start,0]
                np.intp_t *pkk = &kk[0]
            with nogil:
                query_knn(self.cself, pdd, pii,
                    pxx, stop-start, pkk, kk.shape[0], kmax, eps, p, distance_upper_bound)

        if (n_jobs == -1):
            n_jobs = number_of_processors

        _run_threads(_thread_func, n, n_jobs)

        # massage the output in conformabity to the documented behavior

        if sizeof(long) < sizeof(np.intp_t):
            # ... e.g. Windows 64
            overflown = False
            for i in range(n):
                for j in range(len(k)):
                    if ii[i,j] > <np.intp_t>LONG_MAX:
                        # C long overlow, return array of dtype=np.int_p
                        overflown = True
                        break
                if overflown:
                    break

            if overflown:
                ddret = np.reshape(dd,retshape+(len(k),))
                iiret = np.reshape(ii,retshape+(len(k),))
            else:
                ddret = np.reshape(dd,retshape+(len(k),))
                iiret = np.reshape(ii,retshape+(len(k),)).astype(int)

        else:
            # ... most other platforms
            ddret = np.reshape(dd,retshape+(len(k),))
            iiret = np.reshape(ii,retshape+(len(k),))

        if nearest:
            ddret = ddret[..., 0]
            iiret = iiret[..., 0]
            # the only case where we return a python scalar
            if single:
                ddret = float(ddret)
                iiret = int(iiret)

        return ddret, iiret

    # ----------------
    # query_ball_point
    # ----------------

    def query_ball_point(cKDTree self, object x, object r,
                         np.float64_t p=2., np.float64_t eps=0, n_jobs=1,
                         return_sorted=None,
                         return_length=False):
        """
        query_ball_point(self, x, r, p=2., eps=0)

        Find all points within distance r of point(s) x.

        Parameters
        ----------
        x : array_like, shape tuple + (self.m,)
            The point or points to search for neighbors of.
        r : array_like, float
            The radius of points to return, shall broadcast to the length of x.
        p : float, optional
            Which Minkowski p-norm to use.  Should be in the range [1, inf].
            A finite large p may cause a ValueError if overflow can occur.
        eps : nonnegative float, optional
            Approximate search. Branches of the tree are not explored if their
            nearest points are further than ``r / (1 + eps)``, and branches are
            added in bulk if their furthest points are nearer than
            ``r * (1 + eps)``.
        n_jobs : int, optional
            Number of jobs to schedule for parallel processing. If -1 is given
            all processors are used. Default: 1.
        return_sorted : bool, optional
            Sorts returned indicies if True and does not sort them if False. If
            None, does not sort single point queries, but does sort
            multi-point queries which was the behavior before this option
            was added.

            .. versionadded:: 1.2.0
        return_length: bool, optional
            Return the number of points inside the radius instead of a list
            of the indices.
            .. versionadded:: 1.3.0

        Returns
        -------
        results : list or array of lists
            If `x` is a single point, returns a list of the indices of the
            neighbors of `x`. If `x` is an array of points, returns an object
            array of shape tuple containing lists of neighbors.

        Notes
        -----
        If you have many points whose neighbors you want to find, you may save
        substantial amounts of time by putting them in a cKDTree and using
        query_ball_tree.

        Examples
        --------
        >>> from scipy import spatial
        >>> x, y = np.mgrid[0:4, 0:4]
        >>> points = np.c_[x.ravel(), y.ravel()]
        >>> tree = spatial.cKDTree(points)
        >>> tree.query_ball_point([2, 0], 1)
        [4, 8, 9, 12]

        """

        cdef:
            const np.float64_t[::1] vrr
            const np.float64_t[:, ::1] vxx
            object[::1] vout
            np.intp_t[::1] vlen
            list tmp
            np.intp_t i, j, n, m
            np.intp_t xndim

        xshape = np.shape(x)
        if xshape[-1] != self.m:
            raise ValueError("Searching for a %d-dimensional point in a "
                             "%d-dimensional KDTree" %
                                 (int(xshape[-1]), int(self.m)))

        vxx = np.ascontiguousarray(x, dtype=np.float64).reshape(-1, self.m)
        vrr = np.ascontiguousarray(np.broadcast_to(r, xshape[:-1]), dtype=np.float64).reshape(-1)

        retshape = xshape[:-1]

        # scalar query if xndim == 1
        xndim = len(xshape)

        # allocate an array of std::vector<npy_intp>
        n = np.prod(retshape)

        if return_length:
            result = np.empty(retshape, dtype=np.intp)
            vlen = result.reshape(-1)
        else:
            result = np.empty(retshape, dtype=object)
            vout = result.reshape(-1)

        def _thread_func(np.intp_t start, np.intp_t stop):
            cdef:
                vector[np.intp_t] **vvres
                np.intp_t i
                np.intp_t *cur
                int rlen
                const np.float64_t *pvxx
                const np.float64_t *pvrr

            rlen = <int> return_length

            try:
                vvres = (<vector[np.intp_t] **>
                    PyMem_Malloc((stop-start) * sizeof(void*)))
                if vvres == NULL:
                    raise MemoryError()

                memset(<void*> vvres, 0, (stop-start) * sizeof(void*))

                for i in range(stop - start):
                    vvres[i] = new vector[np.intp_t]()

                pvxx = &vxx[start, 0]
                pvrr = &vrr[start + 0]
                
                with nogil:
                    query_ball_point(self.cself, pvxx,
                        pvrr, p, eps, stop - start, vvres, rlen)

                for i in range(stop - start):
                    if return_length:
                        vlen[start + i] = vvres[i].front()
                        continue

                    if return_sorted:
                        with nogil:
                            sort(vvres[i].begin(), vvres[i].end())
                    elif return_sorted is None and xndim > 1:
                        # compatibility with the old bug not sorting scalar queries.
                        with nogil:
                            sort(vvres[i].begin(), vvres[i].end())

                    m = <np.intp_t> (vvres[i].size())
                    tmp = m * [None]

                    cur = vvres[i].data()
                    for j in range(m):
                        tmp[j] = cur[0]
                        cur += 1
                    vout[start + i] = tmp
            finally:
                if vvres != NULL:
                    for i in range(stop-start):
                        if vvres[i] != NULL:
                            del vvres[i]
                    PyMem_Free(vvres)

        # multithreading logic is similar to cKDTree.query
        if n_jobs == -1:
            n_jobs = number_of_processors

        _run_threads(_thread_func, n, n_jobs)

        if xndim == 1: # scalar query, unpack result.
            result = result[()]
        return result

    # ---------------
    # query_ball_tree
    # ---------------

    def query_ball_tree(cKDTree self, cKDTree other,
                        np.float64_t r, np.float64_t p=2., np.float64_t eps=0):
        """
        query_ball_tree(self, other, r, p=2., eps=0)

        Find all pairs of points whose distance is at most r

        Parameters
        ----------
        other : cKDTree instance
            The tree containing points to search against.
        r : float
            The maximum distance, has to be positive.
        p : float, optional
            Which Minkowski norm to use.  `p` has to meet the condition
            ``1 <= p <= infinity``.
            A finite large p may cause a ValueError if overflow can occur.
        eps : float, optional
            Approximate search.  Branches of the tree are not explored
            if their nearest points are further than ``r/(1+eps)``, and
            branches are added in bulk if their furthest points are nearer
            than ``r * (1+eps)``.  `eps` has to be non-negative.

        Returns
        -------
        results : list of lists
            For each element ``self.data[i]`` of this tree, ``results[i]`` is a
            list of the indices of its neighbors in ``other.data``.

        """

        cdef:
            vector[np.intp_t] **vvres
            np.intp_t i, j, n, m
            np.intp_t *cur
            list results
            list tmp

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to query_ball_tree have different "
                             "dimensionality")

        n = self.n

        try:

            # allocate an array of std::vector<npy_intp>
            vvres = (<vector[np.intp_t] **>
                PyMem_Malloc(n * sizeof(void*)))
            if vvres == NULL:
                raise MemoryError()

            memset(<void*> vvres, 0, n * sizeof(void*))

            for i in range(n):
                vvres[i] = new vector[np.intp_t]()

            # query in C++
            with nogil:
                query_ball_tree(self.cself, other.cself, r, p, eps, vvres)

            # store the results in a list of lists
            results = n * [None]
            for i in range(n):
                m = <np.intp_t> (vvres[i].size())
                if NPY_LIKELY(m > 0):
                    tmp = m * [None]
                    with nogil:
                        sort(vvres[i].begin(), vvres[i].end())
                    cur = vvres[i].data()
                    for j in range(m):
                        tmp[j] = cur[0]
                        cur += 1
                    results[i] = tmp
                else:
                    results[i] = []

        finally:
            if vvres != NULL:
                for i in range(n):
                    if vvres[i] != NULL:
                        del vvres[i]
                PyMem_Free(vvres)

        return results

    # -----------
    # query_pairs
    # -----------

    def query_pairs(cKDTree self, np.float64_t r, np.float64_t p=2.,
                    np.float64_t eps=0, output_type='set'):
        """
        query_pairs(self, r, p=2., eps=0)

        Find all pairs of points whose distance is at most r.

        Parameters
        ----------
        r : positive float
            The maximum distance.
        p : float, optional
            Which Minkowski norm to use.  ``p`` has to meet the condition
            ``1 <= p <= infinity``.
            A finite large p may cause a ValueError if overflow can occur.
        eps : float, optional
            Approximate search.  Branches of the tree are not explored
            if their nearest points are further than ``r/(1+eps)``, and
            branches are added in bulk if their furthest points are nearer
            than ``r * (1+eps)``.  `eps` has to be non-negative.
        output_type : string, optional
            Choose the output container, 'set' or 'ndarray'. Default: 'set'

        Returns
        -------
        results : set or ndarray
            Set of pairs ``(i,j)``, with ``i < j``, for which the corresponding
            positions are close. If output_type is 'ndarray', an ndarry is
            returned instead of a set.

        """

        cdef ordered_pairs results

        results = ordered_pairs()
        
        with nogil:
            query_pairs(self.cself, r, p, eps, results.buf)

        if output_type == 'set':
            return results.set()
        elif output_type == 'ndarray':
            return results.ndarray()
        else:
            raise ValueError("Invalid output type")

    def _build_weights(cKDTree self, object weights):
        """
        _build_weights(weights)

        Compute weights of nodes from weights of data points. This will sum
        up the total weight per node. This function is used internally.

        Parameters
        ----------
        weights : array_like
            weights of data points; must be the same length as the data points.
            currently only scalar weights are supported. Therefore the weights
            array must be 1 dimensional.

        Returns
        -------
        node_weights : array_like
            total weight for each KD-Tree node.

        """
        cdef: 
            np.intp_t num_of_nodes
            np.float64_t [::1] node_weights
            np.float64_t [::1] proper_weights
            np.float64_t *pnw
            np.float64_t *ppw

        num_of_nodes = self.cself.tree_buffer.size();
        node_weights = np.empty(num_of_nodes, dtype=np.float64)

        # FIXME: use templates to avoid the type conversion
        proper_weights = np.ascontiguousarray(weights, dtype=np.float64)

        if len(proper_weights) != self.n:
            raise ValueError('Number of weights differ from the number of data points')

        pnw = &node_weights[0]
        ppw = &proper_weights[0]

        with nogil:
            build_weights(self.cself, pnw, ppw)

        return node_weights

    # ---------------
    # count_neighbors
    # ---------------

    @cython.boundscheck(False)
    def count_neighbors(cKDTree self, cKDTree other, object r, np.float64_t p=2.,
                        object weights=None, int cumulative=True):
        """
        count_neighbors(self, other, r, p=2., weights=None, cumulative=True)

        Count how many nearby pairs can be formed. (pair-counting)

        Count the number of pairs (x1,x2) can be formed, with x1 drawn
        from self and x2 drawn from ``other``, and where
        ``distance(x1, x2, p) <= r``.

        Data points on self and other are optionally weighted by the ``weights``
        argument. (See below)

        The algorithm we implement here is based on [1]_. See notes for further discussion.

        Parameters
        ----------
        other : cKDTree instance
            The other tree to draw points from, can be the same tree as self.
        r : float or one-dimensional array of floats
            The radius to produce a count for. Multiple radii are searched with
            a single tree traversal.
            If the count is non-cumulative(``cumulative=False``), ``r`` defines
            the edges of the bins, and must be non-decreasing.
        p : float, optional
            1<=p<=infinity.
            Which Minkowski p-norm to use.
            Default 2.0.
            A finite large p may cause a ValueError if overflow can occur.
        weights : tuple, array_like, or None, optional
            If None, the pair-counting is unweighted.
            If given as a tuple, weights[0] is the weights of points in ``self``, and
            weights[1] is the weights of points in ``other``; either can be None to
            indicate the points are unweighted.
            If given as an array_like, weights is the weights of points in ``self``
            and ``other``. For this to make sense, ``self`` and ``other`` must be the
            same tree. If ``self`` and ``other`` are two different trees, a ``ValueError``
            is raised.
            Default: None
        cumulative : bool, optional
            Whether the returned counts are cumulative. When cumulative is set to ``False``
            the algorithm is optimized to work with a large number of bins (>10) specified
            by ``r``. When ``cumulative`` is set to True, the algorithm is optimized to work
            with a small number of ``r``. Default: True

        Returns
        -------
        result : scalar or 1-D array
            The number of pairs. For unweighted counts, the result is integer.
            For weighted counts, the result is float.
            If cumulative is False, ``result[i]`` contains the counts with
            ``(-inf if i == 0 else r[i-1]) < R <= r[i]``

        Notes
        -----
        Pair-counting is the basic operation used to calculate the two point
        correlation functions from a data set composed of position of objects.

        Two point correlation function measures the clustering of objects and
        is widely used in cosmology to quantify the large scale structure
        in our Universe, but it may be useful for data analysis in other fields
        where self-similar assembly of objects also occur.

        The Landy-Szalay estimator for the two point correlation function of
        ``D`` measures the clustering signal in ``D``. [2]_

        For example, given the position of two sets of objects,

        - objects ``D`` (data) contains the clustering signal, and

        - objects ``R`` (random) that contains no signal,

        .. math::

             \\xi(r) = \\frac{<D, D> - 2 f <D, R> + f^2<R, R>}{f^2<R, R>},

        where the brackets represents counting pairs between two data sets
        in a finite bin around ``r`` (distance), corresponding to setting
        `cumulative=False`, and ``f = float(len(D)) / float(len(R))`` is the
        ratio between number of objects from data and random.

        The algorithm implemented here is loosely based on the dual-tree
        algorithm described in [1]_. We switch between two different
        pair-cumulation scheme depending on the setting of ``cumulative``.
        The computing time of the method we use when for
        ``cumulative == False`` does not scale with the total number of bins.
        The algorithm for ``cumulative == True`` scales linearly with the
        number of bins, though it is slightly faster when only
        1 or 2 bins are used. [5]_.

        As an extension to the naive pair-counting,
        weighted pair-counting counts the product of weights instead
        of number of pairs.
        Weighted pair-counting is used to estimate marked correlation functions
        ([3]_, section 2.2),
        or to properly calculate the average of data per distance bin
        (e.g. [4]_, section 2.1 on redshift).

        .. [1] Gray and Moore,
               "N-body problems in statistical learning",
               Mining the sky, 2000,
               https://arxiv.org/abs/astro-ph/0012333

        .. [2] Landy and Szalay,
               "Bias and variance of angular correlation functions",
               The Astrophysical Journal, 1993,
               http://adsabs.harvard.edu/abs/1993ApJ...412...64L

        .. [3] Sheth, Connolly and Skibba,
               "Marked correlations in galaxy formation models",
               Arxiv e-print, 2005,
               https://arxiv.org/abs/astro-ph/0511773

        .. [4] Hawkins, et al.,
               "The 2dF Galaxy Redshift Survey: correlation functions,
               peculiar velocities and the matter density of the Universe",
               Monthly Notices of the Royal Astronomical Society, 2002,
               http://adsabs.harvard.edu/abs/2003MNRAS.346...78H

        .. [5] https://github.com/scipy/scipy/pull/5647#issuecomment-168474926

        """
        cdef:
            int r_ndim
            np.intp_t n_queries, i
            np.float64_t[::1] real_r
            np.float64_t[::1] fresults
            np.intp_t[::1] iresults
            np.float64_t[::1] w1, w1n
            np.float64_t[::1] w2, w2n
            np.float64_t *w1p = NULL
            np.float64_t *w1np = NULL
            np.float64_t *w2p = NULL
            np.float64_t *w2np = NULL
            np.float64_t *prr
            np.intp_t *pir
            np.float64_t *pfr
            int cum

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to count_neighbors have different "
                             "dimensionality")

        # Make a copy of r array to ensure it's contiguous and to modify it
        # below
        r_ndim = len(np.shape(r))
        if r_ndim > 1:
            raise ValueError("r must be either a single value or a "
                             "one-dimensional array of values")
        real_r = np.array(r, ndmin=1, dtype=np.float64, copy=True)
        if not cumulative:
            for i in range(real_r.shape[0] - 1):
                if real_r[i] > real_r[i + 1]:
                    raise ValueError("r must be non-decreasing for non-cumulative counting.");
        real_r, uind, inverse = np.unique(real_r, return_inverse=True, return_index=True)
        n_queries = real_r.shape[0]

        # Internally, we represent all distances as distance ** p
        if not ckdtree_isinf(p):
            for i in range(n_queries):
                if not ckdtree_isinf(real_r[i]):
                    real_r[i] = real_r[i] ** p

        if weights is None:
            self_weights = other_weights = None
        elif isinstance(weights, tuple):
            self_weights, other_weights = weights
        else:
            self_weights = other_weights = weights
            if other is not self:
                raise ValueError("Two different trees are used. Specify weights for both in a tuple.")

        cum = <int> cumulative

        if self_weights is None and other_weights is None:
            int_result = True
            # unweighted, use the integer arithmetics
            results = np.zeros(n_queries + 1, dtype=np.intp)

            iresults = results
            
            prr = &real_r[0]
            pir = &iresults[0]
            
            with nogil:
                count_neighbors_unweighted(self.cself, other.cself, n_queries,
                            prr, pir, p, cum)

        else:
            int_result = False

            # weighted / half weighted, use the floating point arithmetics
            if self_weights is not None:
                w1 = np.ascontiguousarray(self_weights, dtype=np.float64)
                w1n = self._build_weights(w1)
                w1p = &w1[0]
                w1np = &w1n[0]
            if other_weights is not None:
                w2 = np.ascontiguousarray(other_weights, dtype=np.float64)
                w2n = other._build_weights(w2)
                w2p = &w2[0]
                w2np = &w2n[0]

            results = np.zeros(n_queries + 1, dtype=np.float64)
            fresults = results
            
            prr = &real_r[0]
            pfr = &fresults[0]
            
            with nogil:
                count_neighbors_weighted(self.cself, other.cself,
                                    w1p, w2p, w1np, w2np,
                                    n_queries,
                                    prr, pfr, p, cum)

        results2 = np.zeros(inverse.shape, results.dtype)
        if cumulative:
            # copy out the results (taking care of duplication and sorting)
            results2[...] = results[inverse]
        else:
            # keep the identical ones zero
            # this could have been done in a more readable way.
            results2[uind] = results[inverse][uind]
        results = results2

        if r_ndim == 0:
            if int_result and results[0] <= <np.intp_t> LONG_MAX:
                return int(results[0])
            else:
                return results[0]
        else:
            return results

    # ----------------------
    # sparse_distance_matrix
    # ----------------------

    def sparse_distance_matrix(cKDTree self, cKDTree other,
                               np.float64_t max_distance,
                               np.float64_t p=2.,
                               output_type='dok_matrix'):
        """
        sparse_distance_matrix(self, other, max_distance, p=2.)

        Compute a sparse distance matrix

        Computes a distance matrix between two cKDTrees, leaving as zero
        any distance greater than max_distance.

        Parameters
        ----------
        other : cKDTree

        max_distance : positive float

        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            A finite large p may cause a ValueError if overflow can occur.

        output_type : string, optional
            Which container to use for output data. Options: 'dok_matrix',
            'coo_matrix', 'dict', or 'ndarray'. Default: 'dok_matrix'.

        Returns
        -------
        result : dok_matrix, coo_matrix, dict or ndarray
            Sparse matrix representing the results in "dictionary of keys"
            format. If a dict is returned the keys are (i,j) tuples of indices.
            If output_type is 'ndarray' a record array with fields 'i', 'j',
            and 'v' is returned,
        """

        cdef coo_entries res

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to sparse_distance_matrix have "
                             "different dimensionality")
        # do the query
        res = coo_entries()
        
        with nogil:
            sparse_distance_matrix(
                self.cself, other.cself, p, max_distance, res.buf)

        if output_type == 'dict':
            return res.dict()
        elif output_type == 'ndarray':
            return res.ndarray()
        elif output_type == 'coo_matrix':
            return res.coo_matrix(self.n, other.n)
        elif output_type == 'dok_matrix':
            return res.dok_matrix(self.n, other.n)
        else:
            raise ValueError('Invalid output type')


    # ----------------------
    # pickle
    # ----------------------

    def __getstate__(cKDTree self):
        cdef object state
        cdef np.intp_t size
        cdef ckdtree * cself = self.cself
        size = cself.tree_buffer.size() * sizeof(ckdtreenode)

        cdef np.ndarray tree = np.asarray(<char[:size]> <char*> cself.tree_buffer.data())

        state = (tree.copy(), self.data.copy(), self.n, self.m, self.leafsize,
                      self.maxes, self.mins, self.indices.copy(),
                      self.boxsize, self.boxsize_data)
        return state

    def __setstate__(cKDTree self, state):
        cdef np.ndarray tree
        cdef ckdtree * cself = self.cself
        cdef np.ndarray mytree

        # unpack the state
        (tree, self.data, self.cself.n, self.cself.m, self.cself.leafsize,
            self.maxes, self.mins, self.indices, self.boxsize, self.boxsize_data) = state

        cself.tree_buffer = new vector[ckdtreenode]()
        cself.tree_buffer.resize(tree.size // sizeof(ckdtreenode))

        mytree = np.asarray(<char[:tree.size]> <char*> cself.tree_buffer.data())

        # set raw pointers
        self._pre_init()

        # copy the tree data
        mytree[:] = tree


        # set up the tree structure pointers
        self._post_init()

def _run_threads(_thread_func, n, n_jobs):
    if n_jobs > 1:
        ranges = [(j * n // n_jobs, (j + 1) * n // n_jobs)
                        for j in range(n_jobs)]

        threads = [threading.Thread(target=_thread_func,
                   args=(start, end))
                   for start, end in ranges]
        for t in threads:
            t.daemon = True
            t.start()
        for t in threads:
            t.join()

    else:
        _thread_func(0, n)
