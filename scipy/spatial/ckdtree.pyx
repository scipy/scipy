# Copyright Anne M. Archibald 2008
# Additional contributions by Patrick Varilly and Sturla Molden 2012
# Revision by Sturla Molden 2015
# Balanced kd-tree construction written by Jake Vanderplas for scikit-learn
# Released under the scipy license

# distutils: language = c++

import numpy as np
import scipy.sparse

cimport numpy as np
from numpy.math cimport INFINITY
    
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memset, memcpy

cimport cython

from multiprocessing import cpu_count
import threading

cdef extern from "limits.h":
    long LONG_MAX

cdef extern from "ckdtree_methods.h":
    int number_of_processors
    
number_of_processors = cpu_count()

from libcpp.vector cimport vector
from libc cimport string

__all__ = ['cKDTree']

    
# Borrowed references
# ===================

cdef extern from *:
    
    struct ckdtree:
        pass
        
    int NPY_LIKELY(int)
    int NPY_UNLIKELY(int)

       
cdef extern from "ckdtree_decl.h":
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

def new_object(obj):
    return obj.__new__(obj)
 
cdef extern from "cpp_utils.h": 
    object pickle_tree_buffer(vector[ckdtreenode] *buf)    
    object unpickle_tree_buffer(vector[ckdtreenode] *buf, object src)
    ckdtreenode *tree_buffer_root(vector[ckdtreenode] *buf)
    ordered_pair *ordered_pair_vector_buf(vector[ordered_pair] *buf)
    coo_entry *coo_entry_vector_buf(vector[coo_entry] *buf)
    void *tree_buffer_pointer(vector[ckdtreenode] *buf)
    np.intp_t *npy_intp_vector_buf(vector[np.intp_t] *buf)
    np.float64_t *npy_float64_vector_buf(vector[np.float64_t] *buf)
    ctypedef void *intvector_ptr_t 



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
            pr = coo_entry_vector_buf(self.buf) 
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
            pr = coo_entry_vector_buf(self.buf)        
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
            pr = ordered_pair_vector_buf(self.buf) 
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
        pair = ordered_pair_vector_buf(self.buf)
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

cdef extern from "ckdtree_methods.h":

    # External build and query methods in C++. These will internally
    # release the GIL to avoid locking up the interpreter.
    
    int ckdtree_isinf(np.float64_t x)
    
    object build_ckdtree(ckdtree *self, 
                         np.intp_t start_idx, 
                         np.intp_t end_idx,
                         np.float64_t *maxes, 
                         np.float64_t *mins, 
                         int _median, 
                         int _compact)

    object build_weights(ckdtree *self, 
                         np.float64_t *node_weights,
                         np.float64_t *weights)
       
    object query_knn(const ckdtree *self, 
                     np.float64_t *dd, 
                     np.intp_t    *ii, 
                     const np.float64_t *xx,
                     const np.intp_t    n,
                     const np.intp_t    *k, 
                     const np.intp_t    nk, 
                     const np.intp_t    kmax, 
                     const np.float64_t eps, 
                     const np.float64_t p, 
                     const np.float64_t distance_upper_bound) 
                     
    object query_pairs(const ckdtree *self, 
                       const np.float64_t r, 
                       const np.float64_t p, 
                       const np.float64_t eps,
                       vector[ordered_pair] *results)
                       
    object count_neighbors_unweighted(const ckdtree *self,
                           const ckdtree *other,
                           np.intp_t     n_queries,
                           np.float64_t  *real_r,
                           np.intp_t     *results,
                           const np.float64_t p)

    object count_neighbors_weighted(const ckdtree *self,
                           const ckdtree *other,
                           np.float64_t  *self_weights,
                           np.float64_t  *other_weights,
                           np.float64_t  *self_node_weights,
                           np.float64_t  *other_node_weights,
                           np.intp_t     n_queries,
                           np.float64_t  *real_r,
                           np.float64_t     *results,
                           const np.float64_t p)
                           
    object query_ball_point(const ckdtree *self,
                            const np.float64_t *x,
                            const np.float64_t r,
                            const np.float64_t p,
                            const np.float64_t eps,
                            const np.intp_t n_queries,
                            vector[np.intp_t] **results)

    object query_ball_tree(const ckdtree *self,
                           const ckdtree *other,
                           const np.float64_t r,
                           const np.float64_t p,
                           const np.float64_t eps,
                           vector[np.intp_t] **results)                     
     
    object sparse_distance_matrix(const ckdtree *self,
                                  const ckdtree *other,
                                  const np.float64_t p,
                                  const np.float64_t max_distance,
                                  vector[coo_entry] *results)                    
                      
                      
cdef public class cKDTree [object ckdtree, type ckdtree_type]:
    """
    cKDTree(data, leafsize=16, compact_nodes=True, copy_data=False,
            balanced_tree=True)

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
    size : int
        The number of nodes in the tree.

    See Also
    --------
    KDTree : Implementation of `cKDTree` in pure Python

    """
    cdef:
        vector[ckdtreenode]      *tree_buffer
        ckdtreenode              *ctree 
        readonly cKDTreeNode     tree 
        readonly np.ndarray      data
        np.float64_t             *raw_data
        readonly np.intp_t       n, m
        readonly np.intp_t       leafsize
        readonly np.ndarray      maxes
        np.float64_t             *raw_maxes
        readonly np.ndarray      mins
        np.float64_t             *raw_mins
        readonly np.ndarray      indices
        np.intp_t                *raw_indices
        np.ndarray               _median_workspace
        readonly object          boxsize
        np.ndarray               boxsize_data
        np.float64_t             *raw_boxsize_data

    property size:
        def __get__(self):
            return self.tree_buffer.size()

    def __cinit__(cKDTree self):
        self.tree_buffer = NULL        
            
    def __init__(cKDTree self, data, np.intp_t leafsize=16, compact_nodes=True, 
            copy_data=False, balanced_tree=True, boxsize=None):
        cdef np.ndarray[np.float64_t, ndim=2] data_arr
        cdef np.float64_t *tmp
        cdef int _median, _compact
        cdef np.ndarray[np.float64_t, ndim=1] boxsize_arr
        data_arr = np.ascontiguousarray(data, dtype=np.float64)
        if copy_data and (data_arr is data):
            data_arr = data_arr.copy()
        self.data = data_arr
        self.n = data_arr.shape[0]
        self.m = data_arr.shape[1]
        self.leafsize = leafsize
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")

        if boxsize is None:
            self.boxsize = None
            self.raw_boxsize_data = NULL
            self.boxsize_data = None
        else:
            boxsize_arr = np.empty(2 * self.m, dtype=np.float64)
            boxsize_arr[:] = boxsize
            boxsize_arr[self.m:] = 0.5 * boxsize_arr[:self.m]
            # FIXME: how to use a matching del if new is used?
            self.boxsize_data = boxsize_arr
            self.raw_boxsize_data = <np.float64_t*> np.PyArray_DATA(boxsize_arr)
            self.boxsize = boxsize_arr[:self.m].copy()
            if (self.data >= self.boxsize[None, :]).any():
                raise ValueError("Some input data are greater than the size of the periodic box.")
            if (self.data < 0).any():
                raise ValueError("Negative input data are outside of the periodic box.")

        self.maxes = np.ascontiguousarray(np.amax(self.data,axis=0), dtype=np.float64)
        self.mins = np.ascontiguousarray(np.amin(self.data,axis=0), dtype=np.float64)
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=np.intp))

        self.raw_data = <np.float64_t*> np.PyArray_DATA(self.data)
        self.raw_maxes = <np.float64_t*> np.PyArray_DATA(self.maxes)
        self.raw_mins = <np.float64_t*> np.PyArray_DATA(self.mins)
        self.raw_indices = <np.intp_t*> np.PyArray_DATA(self.indices)

        _compact = 1 if compact_nodes else 0
        _median = 1 if balanced_tree else 0
        if _median:
            self._median_workspace = np.zeros(self.n)

        self.tree_buffer = new vector[ckdtreenode]()
        
        try:
            tmp = <np.float64_t*> PyMem_Malloc(self.m*2*sizeof(np.float64_t))
            if tmp == NULL: raise MemoryError()            
            memcpy(tmp, self.raw_maxes, self.m*sizeof(np.float64_t))
            memcpy(tmp + self.m, self.raw_mins, self.m*sizeof(np.float64_t))
            build_ckdtree(<ckdtree*> self, 0, self.n, tmp, tmp + self.m, 
                _median, _compact)
        finally:
            PyMem_Free(tmp)
                
        self._median_workspace = None
        
        # set up the tree structure pointers
        self.ctree = tree_buffer_root(self.tree_buffer)
        self._post_init(self.ctree)
        
        # make the tree viewable from Python
        self.tree = cKDTreeNode()
        self.tree._node = self.ctree
        self.tree._data = self.data
        self.tree._indices = self.indices
        self.tree.level = 0
        self.tree._setup()
        
                
    cdef int _post_init(cKDTree self, ckdtreenode *node) except -1:
        # recurse the tree and re-initialize
        # "less" and "greater" fields
        if node.split_dim == -1:
            # leafnode
            node.less = NULL
            node.greater = NULL
        else:
            node.less = self.ctree + node._less
            node.greater = self.ctree + node._greater
            self._post_init(node.less)
            self._post_init(node.greater)
        return 0
        

    def __dealloc__(cKDTree self):
        if self.tree_buffer != NULL:
            del self.tree_buffer

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
        
        >>> tree = cKDTree(data)

        To query the nearest neighbours and return squeezed result, use

        >>> dd, ii = tree.query(x, k=1)

        To query the nearest neighbours and return unsqueezed result, use

        >>> dd, ii = tree.query(x, k=[1])

        To query the second nearest neighbours and return unsqueezed result, use

        >>> dd, ii = tree.query(x, k=[2])

        To query the first and second nearest neighbours, use

        >>> dd, ii = tree.query(x, k=2)

        or, be more specific

        >>> dd, ii = tree.query(x, k=[1, 2])


        """
        
        cdef:
            np.intp_t n, i, j
            int overflown
        
        x_arr = np.asarray(x, dtype=np.float64)
        if x_arr.ndim == 0 or x_arr.shape[x_arr.ndim - 1] != self.m:
            raise ValueError("x must consist of vectors of length %d but "
                             "has shape %s" % (int(self.m), np.shape(x)))
        if p < 1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        if x_arr.ndim == 1:
            single = True
            x_arr = x_arr[np.newaxis,:]
        else:
            single = False

        nearest = False
        if np.isscalar(k):
            if k == 1:
                nearest = True
            k = np.arange(1, k + 1)
    
        retshape = np.shape(x)[:-1]
        n = <np.intp_t> np.prod(retshape)
        xx = np.ascontiguousarray(x_arr).reshape(n, self.m)
        dd = np.empty((n,len(k)),dtype=np.float64)
        dd.fill(INFINITY)
        ii = np.empty((n,len(k)),dtype=np.intp)
        ii.fill(self.n)

        # Do the query in an external C++ function. 
        # The GIL will be released in the external query function.
        def _thread_func(self, np.intp_t start, np.intp_t stop):
            cdef: 
                np.ndarray[np.intp_t,ndim=2] _ii = ii
                np.ndarray[np.float64_t,ndim=2] _dd = dd
                np.ndarray[np.float64_t,ndim=2] _xx = xx
                np.ndarray[np.intp_t,ndim=1] _k = np.array(k, dtype=np.intp)
            
            kmax = np.max(k)

            query_knn(<ckdtree*>self, &_dd[start,0], &_ii[start,0], 
                &_xx[start,0], stop-start, &_k[0], len(k), kmax, eps, p, distance_upper_bound)
        
        if (n_jobs == -1): 
            n_jobs = number_of_processors
        
        if n_jobs > 1:
            # static scheduling without load balancing is good enough
                             
            ranges = [(j * n // n_jobs, (j + 1) * n // n_jobs)
                            for j in range(n_jobs)]
            
            # There might be n_jobs+1 threads spawned here, but only n_jobs of 
            # them will do significant work.
            threads = [threading.Thread(target=_thread_func,
                                args=(self, start, stop)) for start, stop in ranges]

            # Set the daemon flag so the process can be aborted, 
            # start all threads and wait for completion.
            for t in threads:
                t.daemon = True
                t.start()
            for t in threads: 
                t.join()
        else:
            _thread_func(self, 0, n)
                
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

    def query_ball_point(cKDTree self, object x, np.float64_t r,
                         np.float64_t p=2., np.float64_t eps=0, n_jobs=1):
        """
        query_ball_point(self, x, r, p=2., eps=0)
        
        Find all points within distance r of point(s) x.

        Parameters
        ----------
        x : array_like, shape tuple + (self.m,)
            The point or points to search for neighbors of.
        r : positive float
            The radius of points to return.
        p : float, optional
            Which Minkowski p-norm to use.  Should be in the range [1, inf].
        eps : nonnegative float, optional
            Approximate search. Branches of the tree are not explored if their
            nearest points are further than ``r / (1 + eps)``, and branches are
            added in bulk if their furthest points are nearer than
            ``r * (1 + eps)``.
        n_jobs : int, optional
            Number of jobs to schedule for parallel processing. If -1 is given
            all processors are used. Default: 1.

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
        >>> points = zip(x.ravel(), y.ravel())
        >>> tree = spatial.cKDTree(points)
        >>> tree.query_ball_point([2, 0], 1)
        [4, 8, 9, 12]

        """
        
        cdef:
            np.ndarray[np.float64_t, ndim=1, mode="c"] xx
            np.ndarray[np.float64_t, ndim=2, mode="c"] vxx
            vector[np.intp_t] *vres
            vector[np.intp_t] **vvres
            np.uintp_t vvres_uintp
            np.intp_t *cur
            list tmp
            np.intp_t i, j, n, m
        
        vres = NULL
        vvres = NULL
        
        try:
               
            x = np.asarray(x, dtype=np.float64)
            if x.shape[-1] != self.m:
                raise ValueError("Searching for a %d-dimensional point in a "
                                 "%d-dimensional KDTree" % 
                                     (int(x.shape[-1]), int(self.m)))
            if len(x.shape) == 1:
                vres = new vector[np.intp_t]()
                xx = np.ascontiguousarray(x, dtype=np.float64)
                query_ball_point(<ckdtree*> self, &xx[0], r, p, eps, 1, &vres)
                n = <np.intp_t> vres.size()
                tmp = n * [None]
                if NPY_LIKELY(n > 0):
                    cur = npy_intp_vector_buf(vres)
                    for i in range(n):
                        tmp[i] = cur[0]
                        cur += 1
                result = tmp
            
            else:
                retshape = x.shape[:-1]
                
                # allocate an array of std::vector<npy_intp>
                n = np.prod(retshape)
                vvres = (<vector[np.intp_t] **> 
                    PyMem_Malloc(n * sizeof(intvector_ptr_t)))
                if vvres == NULL:
                    raise MemoryError()
                
                memset(<void*> vvres, 0, n * sizeof(intvector_ptr_t))      
            
                for i in range(n):
                    vvres[i] = new vector[np.intp_t]()
                
                result = np.empty(retshape, dtype=object)
                
                vxx = np.zeros((n,self.m), dtype=np.float64)
                i = 0
                for c in np.ndindex(retshape):
                    vxx[i,:] = x[c]
                    i += 1
                    
                # multithreading logic is similar to cKDTree.query
                                        
                if (n_jobs == -1): 
                    n_jobs = number_of_processors
        
                if n_jobs > 1:
                
                    CHUNK = n//n_jobs if n//n_jobs else n
                             
                    def _thread_func(self, _j, _vxx, r, p, eps, _vvres, CHUNK): 
                        cdef: 
                            np.intp_t j = _j
                            np.ndarray[np.float64_t,ndim=2] vxx = _vxx
                            vector[np.intp_t] **vvres                   
                            np.intp_t start = j*CHUNK
                            np.intp_t stop = start + CHUNK
                        stop = n if stop > n else stop
                        vvres = (<vector[np.intp_t] **> 
                                  (<void*> (<np.uintp_t> _vvres)))                                    
                        if start < n:
                            query_ball_point(<ckdtree*>self, &vxx[start,0], 
                                r, p, eps, stop-start, vvres+start)
                                
                    vvres_uintp = <np.uintp_t> (<void*> vvres)
                    threads = [threading.Thread(target=_thread_func,
                               args=(self, j, vxx, r, p, eps,vvres_uintp,CHUNK))
                                  for j in range(1+(n//CHUNK))]
                    for t in threads:
                        t.daemon = True
                        t.start()
                    for t in threads: 
                        t.join()
                                                                
                else:
                
                    query_ball_point(<ckdtree*>self, &vxx[0,0], r, p, eps, 
                        n, vvres)
                
                i = 0
                for c in np.ndindex(retshape):
                    m = <np.intp_t> (vvres[i].size())
                    if NPY_LIKELY(m > 0):
                        tmp = m * [None]
                        cur = npy_intp_vector_buf(vvres[i])
                        for j in range(m):
                            tmp[j] = cur[0]
                            cur += 1
                        result[c] = sorted(tmp)
                    else:
                        result[c] = []
                    i += 1
        
        finally:
            if vres != NULL: 
                del vres
                
            if vvres != NULL:
                for i in range(n):
                    if vvres[i] != NULL:
                        del vvres[i]     
                PyMem_Free(vvres)
                
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
                PyMem_Malloc(n * sizeof(intvector_ptr_t)))
            if vvres == NULL:
                raise MemoryError()
                
            memset(<void*> vvres, 0, n * sizeof(intvector_ptr_t))      
            
            for i in range(n):
                vvres[i] = new vector[np.intp_t]()
        
            # query in C++
            # the GIL will be released in the C++ code
            query_ball_tree(
                <ckdtree*> self, <ckdtree*> other, r, p, eps, vvres)
                          
            # store the results in a list of lists                                        
            results = n * [None]
            for i in range(n):
                m = <np.intp_t> (vvres[i].size())
                if NPY_LIKELY(m > 0):
                    tmp = m * [None]
                    cur = npy_intp_vector_buf(vvres[i]) 
                    for j in range(m):
                        tmp[j] = cur[0]
                        cur += 1
                    results[i] = sorted(tmp)
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
        query_pairs(<ckdtree*> self, r, p, eps, results.buf)
        
        if output_type == 'set':
            return results.set()
        elif output_type == 'ndarray':
            return results.ndarray()
        else:
            raise ValueError("Invalid output type") 

    @cython.boundscheck(False)
    def build_weights(cKDTree self, object weights):
        cdef np.intp_t num_of_nodes
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] node_weights;
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] _weights

        num_of_nodes = self.tree_buffer.size();
        node_weights = np.empty(num_of_nodes, dtype=np.float64)
        _weights = np.ascontiguousarray(weights, dtype=np.float64)

        if len(_weights) != self.n:
            raise ValueError('Number of weights differ from the number of data points')

        build_weights(<ckdtree*> self, <np.float64_t*>np.PyArray_DATA(node_weights),
                            <np.float64_t*> np.PyArray_DATA(weights))
                    
        return node_weights

    # ---------------
    # count_neighbors
    # ---------------

    @cython.boundscheck(False)
    def count_neighbors(cKDTree self, cKDTree other, object r, np.float64_t p=2., 
                        object weights=None, int cumulative=True):
        """
        count_neighbors(self, other, r, p=2., self_weights=None, other_weights=None, cumulative=True)

        Count how many nearby pairs can be formed. (pair-counting)

        Count the number of pairs (x1,x2) can be formed, with x1 drawn
        from self and x2 drawn from ``other``, and where
        ``distance(x1, x2, p) <= r``.
        This is the "two-point correlation" described in Gray and Moore 2000,
        "N-body problems in statistical learning", and the code here is based
        on their algorithm. 

        Data points on self and other are optionally weighted by 
        ``self_weights`` and ``other_weights``. 

        Parameters
        ----------
        other : cKDTree instance
            The other tree to draw points from.
        r : float or one-dimensional array of floats
            The radius to produce a count for. Multiple radii are searched with
            a single tree traversal. 
            If the count is non-cumulative(``cumulative=False``), ``r`` defines 
            the edges of the bins, and must be non-decreasing.
        p : float, optional 
            1<=p<=infinity, default 2.0
            Which Minkowski p-norm to use
        weights : tuple, array_like, or None, optional
            If given as a tuple, weights[0] is the weights of points in self, and
            weights[1] is the weights of points in other; either can be None to 
            indicate the points are unweighted.
            If given as an array_like, weights is the weights of points in self 
            and other. If self and other are two different trees, an Error is raised.
            If None, the pair-counting is unweighted.
        cumulative : bool, optional
            Whether the returned counts are cumulative. Default: True

        Returns
        -------
        result : scalar or 1-D array
            The number of pairs. For unweighted counts, the result is integer.
            For weighted counts, the result is float.
            If cumulative is False, ``result[i]`` contains the counts with
            ``r[i-1] or -inf < R <= r[i]``
        """
        cdef: 
            int r_ndim
            np.intp_t n_queries, i
            np.ndarray[np.float64_t, ndim=1, mode="c"] real_r
            np.ndarray[np.float64_t, ndim=1, mode="c"] fresults
            np.ndarray[np.intp_t, ndim=1, mode="c"] iresults
            np.ndarray[np.float64_t, ndim=1, mode="c"] w1, w1n
            np.ndarray[np.float64_t, ndim=1, mode="c"] w2, w2n
            np.float64_t *w1p
            np.float64_t *w1np
            np.float64_t *w2p
            np.float64_t *w2np

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
            if (real_r[:-1] > real_r[1:]).any():
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

        if self_weights is None and other_weights is None:
            int_result = True
            # unweighted, use the integer arithmetics
            results = np.zeros(n_queries + 1, dtype=np.intp)

            iresults = results
            count_neighbors_unweighted(<ckdtree*> self, <ckdtree*> other, n_queries,
                            &real_r[0], &iresults[0], p)

        else:
            int_result = False
            # weighted / half weighted, use the floating point arithmetics
            if self_weights is not None:
                w1 = np.ascontiguousarray(self_weights, dtype=np.float64)
                w1n = self.build_weights(w1)
                w1p = <np.float64_t*> np.PyArray_DATA(w1)
                w1np = <np.float64_t*> np.PyArray_DATA(w1n)
            else:
                w1p = NULL
                w1np = NULL
            if other_weights is not None:
                w2 = np.ascontiguousarray(other_weights, dtype=np.float64)
                w2n = other.build_weights(w2)
                w2p = <np.float64_t*> np.PyArray_DATA(w2)
                w2np = <np.float64_t*> np.PyArray_DATA(w2n)
            else:
                w2p = NULL
                w2np = NULL

            results = np.zeros(n_queries + 1, dtype=np.float64)
            fresults = results
            count_neighbors_weighted(<ckdtree*> self, <ckdtree*> other,
                                    w1p, w2p, w1np, w2np,
                                    n_queries,
                                    &real_r[0], &fresults[0], p)

        if cumulative:
            results = results.cumsum()
            results = results[inverse]
        else:
            results2 = np.zeros_like(results)
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
        
        output_type : string, optional
            Which container to use for output data. Options: 'dok_matrix',
            'coo_matrix', 'dict', or 'ndarray'. Default: 'dok_matrix'.

        Returns
        -------
        result : dok_matrix, coo_matrix, dict or ndarray
            Sparse matrix representing the results in "dictionary of keys" 
            format. If a dict is returned the keys are (i,j) tuples of indices.
            If output_type is 'ndarray' a record array with fields 'i', 'j',
            and 'k' is returned,
        """
        
        cdef coo_entries res

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to sparse_distance_matrix have "
                             "different dimensionality")                                      
        # do the query
        res = coo_entries()
        sparse_distance_matrix(
                <ckdtree*> self, <ckdtree*> other, p, max_distance, res.buf)
                
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

        
    def __reduce__(self):
        return (new_object, (cKDTree,), self.__getstate__())

    def __getstate__(cKDTree self):
        cdef object state
        cdef object tree = pickle_tree_buffer(self.tree_buffer)
        state = (tree, self.data.copy(), self.n, self.m, self.leafsize,
                      self.maxes, self.mins, self.indices.copy(), 
                      self.boxsize, self.boxsize_data)
        return state
            
    def __setstate__(cKDTree self, state):
        cdef object tree
        self.tree_buffer = new vector[ckdtreenode]()
        
        # unpack the state
        (tree, self.data, self.n, self.m, self.leafsize, 
            self.maxes, self.mins, self.indices, self.boxsize, self.boxsize_data) = state
        
        # copy kd-tree buffer 
        unpickle_tree_buffer(self.tree_buffer, tree)    
        
        # set raw pointers
        self.raw_data = <np.float64_t*>np.PyArray_DATA(self.data)
        self.raw_maxes = <np.float64_t*>np.PyArray_DATA(self.maxes)
        self.raw_mins = <np.float64_t*>np.PyArray_DATA(self.mins)
        self.raw_indices = <np.intp_t*>np.PyArray_DATA(self.indices)
        if self.boxsize_data is not None:
            self.raw_boxsize_data = <np.float64_t*>np.PyArray_DATA(self.boxsize_data)
 
        # set up the tree structure pointers
        self.ctree = tree_buffer_root(self.tree_buffer)
        self._post_init(self.ctree)
        
        # make the tree viewable from Python
        self.tree = cKDTreeNode()
        self.tree._node = self.ctree
        self.tree._data = self.data
        self.tree._indices = self.indices
        self.tree.level = 0
        self.tree._setup() 

