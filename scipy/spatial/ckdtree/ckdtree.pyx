# Copyright Anne M. Archibald 2008
# Additional contributions by Patrick Varilly and Sturla Molden 2012
# Revision by Sturla Molden 2015
# Balanced kd-tree construction written by Jake Vanderplas for scikit-learn
# Released under the scipy license

# distutils: language = c++

import numpy as np
import scipy.sparse

cimport numpy as np
    
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memset

cimport cython

from multiprocessing import cpu_count
import threading

cdef extern from "limits.h":
    long LONG_MAX
    
cdef extern from "ckdtree_cpp_methods.h":
    int number_of_processors
    np.float64_t infinity
    np.float64_t dmax(np.float64_t x, np.float64_t y)
    np.float64_t dabs(np.float64_t x)
    np.float64_t _distance_p(np.float64_t *x, np.float64_t *y,
                       np.float64_t p, np.intp_t k, np.float64_t upperbound)
infinity = np.inf
number_of_processors = cpu_count()

from libcpp.vector cimport vector
from libc cimport string

__all__ = ['cKDTree']

    
# Borrowed references
# ===================

cdef extern from *:
    
    struct ckdtree:
        pass
       
cdef extern from "ckdtree_cpp_decl.h":
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

cdef extern from "ckdtree_cpp_ordered_pair.h":

    struct ordered_pair:
        np.intp_t i
        np.intp_t j

def new_object(obj):
    return obj.__new__(obj)
 
cdef extern from "ckdtree_cpp_utils.h": 
    object pickle_tree_buffer(vector[ckdtreenode] *buf)    
    object unpickle_tree_buffer(vector[ckdtreenode] *buf, object src)
    ckdtreenode *tree_buffer_root(vector[ckdtreenode] *buf)
    ordered_pair *ordered_pair_vector_buf(vector[ordered_pair] *buf)
    void *tree_buffer_pointer(vector[ckdtreenode] *buf)
    np.intp_t *npy_intp_vector_buf(vector[np.intp_t] *buf)
    np.float64_t *npy_float64_vector_buf(vector[np.float64_t] *buf)
    ctypedef void *intvector_ptr_t 


# Splitting routines for a balanced kd-tree
# Code originally written by Jake Vanderplas for scikit-learn

cdef inline void index_swap(np.intp_t *arr, np.intp_t i1, np.intp_t i2):
    """swap the values at index i1 and i2 of arr"""
    cdef np.intp_t tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp

cdef void partition_node_indices(np.float64_t *data,
                                np.intp_t *node_indices,
                                np.intp_t split_dim,
                                np.intp_t split_index,
                                np.intp_t n_features,
                                np.intp_t n_points):
    """Partition points in the node into two equal-sized groups.
    
    Upon return, the values in node_indices will be rearranged such that
    (assuming numpy-style indexing):

        data[node_indices[0:split_index], split_dim]
          <= data[node_indices[split_index], split_dim]

    and

        data[node_indices[split_index], split_dim]
          <= data[node_indices[split_index:n_points], split_dim]

    The algorithm is essentially a partial in-place quicksort around a
    set pivot.
    
    Parameters
    ----------
    data : double pointer
        Pointer to a 2D array of the training data, of shape [N, n_features].
        N must be greater than any of the values in node_indices.
    node_indices : int pointer
        Pointer to a 1D array of length n_points.  This lists the indices of
        each of the points within the current node.  This will be modified
        in-place.
    split_dim : int
        the dimension on which to split.  This will usually be computed via
        the routine ``find_node_split_dim``
    split_index : int
        the index within node_indices around which to split the points.

    Returns
    -------
    status : int
        integer exit status.  On return, the contents of node_indices are
        modified as noted above.

    """
    cdef np.intp_t left, right, midindex, i
    cdef np.float64_t d1, d2
    left = 0
    right = n_points - 1
    
    while True:
        midindex = left
        for i in range(left, right):
            d1 = data[node_indices[i] * n_features + split_dim]
            d2 = data[node_indices[right] * n_features + split_dim]
            if d1 < d2:
                index_swap(node_indices, i, midindex)
                midindex += 1
        index_swap(node_indices, midindex, right)
        if midindex == split_index:
            break
        elif midindex < split_index:
            left = midindex + 1
        else:
            right = midindex - 1
    
    

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

cdef extern from "ckdtree_cpp_methods.h":

    # External query methods in C++. These will internally
    # release the GIL to avoid locking up the interpreter.
       
    object query_knn(const ckdtree *self, 
                     np.float64_t *dd, 
                     np.intp_t    *ii, 
                     const np.float64_t *xx,
                     const np.intp_t    n,
                     const np.intp_t    k, 
                     const np.float64_t eps, 
                     const np.float64_t p, 
                     const np.float64_t distance_upper_bound) 
                     
    object query_pairs(const ckdtree *self, 
                       const np.float64_t r, 
                       const np.float64_t p, 
                       const np.float64_t eps,
                       vector[ordered_pair] *results)
                       
    object count_neighbors(const ckdtree *self,
                           const ckdtree *other,
                           np.intp_t     n_queries,
                           np.float64_t  *real_r,
                           np.intp_t     *results,
                           np.intp_t     *idx,
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
                                  vector[np.intp_t] *results_i,
                                  vector[np.intp_t] *results_j,
                                  vector[np.float64_t] *results_v)                     
                      
                      
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
        the actual data range. This usually gives a more compact tree and 
        faster queries at the expense of longer build time. Default: True.
    copy_data : bool, optional
        If True the data is always copied to protect the kd-tree against 
        data corruption. Default: False.
    balanced_tree : bool, optional    
        If True, the median is used to split the hyperrectangles instead of 
        the midpoint. This usually gives a more compact tree and 
        faster queries at the expense of longer build time. Default: True.

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
    
    
    def __init__(cKDTree self, data, np.intp_t leafsize=16, compact_nodes=True, 
            copy_data=False, balanced_tree=True):
        cdef np.ndarray[np.float64_t, ndim=2] data_arr
        cdef np.float64_t *tmp
        cdef int _median
        data_arr = np.ascontiguousarray(data, dtype=np.float64)
        if copy_data and (data_arr is data):
            data_arr = data_arr.copy()
        self.data = data_arr
        self.n = data_arr.shape[0]
        self.m = data_arr.shape[1]
        self.leafsize = leafsize
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")

        self.maxes = np.ascontiguousarray(np.amax(self.data,axis=0), dtype=np.float64)
        self.mins = np.ascontiguousarray(np.amin(self.data,axis=0), dtype=np.float64)
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=np.intp))

        self.raw_data = <np.float64_t*> np.PyArray_DATA(self.data)
        self.raw_maxes = <np.float64_t*> np.PyArray_DATA(self.maxes)
        self.raw_mins = <np.float64_t*> np.PyArray_DATA(self.mins)
        self.raw_indices = <np.intp_t*> np.PyArray_DATA(self.indices)

        _median = 1 if balanced_tree else 0
        if _median:
            self._median_workspace = np.zeros(self.n)
        
        self.tree_buffer = NULL
        self.tree_buffer = new vector[ckdtreenode]()
        
        if not compact_nodes:
             self.__build(0, self.n, self.raw_maxes, self.raw_mins, _median)
        else:
            try:
                tmp = <np.float64_t*> PyMem_Malloc(self.m*2*sizeof(np.float64_t))
                if tmp == NULL: raise MemoryError()
                self.__build_compact(0, self.n, tmp, tmp+self.m, _median)
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
        

    def __deallocate__(cKDTree self):
        if self.tree_buffer != NULL:
            del self.tree_buffer


    @cython.cdivision(True)
    cdef np.intp_t __build(cKDTree self, np.intp_t start_idx, np.intp_t end_idx,
                       np.float64_t *maxes, np.float64_t *mins, int _median)\
                       except -1:
        cdef:
            ckdtreenode new_node
            np.intp_t   node_index
            np.intp_t   _less, _greater
            ckdtreenode *n
            ckdtreenode *root
            
            np.intp_t i, j, t, p, q, d
            np.float64_t size, split, minval, maxval
            np.float64_t *mids
            np.float64_t *tmp_data_point
        
        # put a new node into the node stack
        self.tree_buffer.push_back(new_node)
        node_index = self.tree_buffer.size() - 1
        root = tree_buffer_root(self.tree_buffer)
        n = root + node_index
                    
        if end_idx-start_idx <= self.leafsize:
            # below brute force limit
            # return leafnode
            n.split_dim = -1
            n.children = end_idx - start_idx
            n.start_idx = start_idx
            n.end_idx = end_idx
            return node_index
        else:
            d = 0 
            size = 0
            for i in range(self.m):
                if maxes[i]-mins[i] > size:
                    d = i
                    size =  maxes[i]-mins[i]
            maxval = maxes[d]
            minval = mins[d]
            if maxval==minval:
                # all points are identical; warn user?
                # return leafnode
                n.split_dim = -1
                n.children = end_idx - start_idx
                n.start_idx = start_idx
                n.end_idx = end_idx
                return node_index
                
            # construct new inner node
            if _median:
                # split on median to create a balanced tree
                # adopted from scikit-learn
                i = (end_idx-start_idx) // 2
                partition_node_indices(self.raw_data,
                                self.raw_indices + start_idx,
                                d,
                                i,
                                self.m,
                                end_idx-start_idx)               
                p = start_idx + i
                split = self.raw_data[self.raw_indices[p]*self.m+d]

            else:
                # split with the sliding midpoint rule
                # this is the default
                split = (maxval+minval)/2
    
            p = start_idx
            q = end_idx - 1
            while p <= q:
                if self.raw_data[self.raw_indices[p]*self.m+d] < split:
                    p += 1
                elif self.raw_data[self.raw_indices[q]*self.m+d] >= split:
                    q -= 1
                else:
                    t = self.raw_indices[p]
                    self.raw_indices[p] = self.raw_indices[q]
                    self.raw_indices[q] = t
                    p += 1
                    q -= 1

            # slide midpoint if necessary
            if p == start_idx:
                # no points less than split
                j = start_idx
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx+1, end_idx):
                    if self.raw_data[self.raw_indices[i]*self.m+d] < split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[start_idx]
                self.raw_indices[start_idx] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = start_idx + 1
                q = start_idx
            elif p == end_idx:
                # no points greater than split
                j = end_idx - 1
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx, end_idx-1):
                    if self.raw_data[self.raw_indices[i]*self.m+d] > split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[end_idx-1]
                self.raw_indices[end_idx-1] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = end_idx - 1
                q = end_idx - 2

            try:
                mids = <np.float64_t*> PyMem_Malloc(sizeof(np.float64_t)*self.m)
                if mids == <np.float64_t*> NULL:
                    raise MemoryError
                        
                for i in range(self.m):
                    mids[i] = maxes[i]
                mids[d] = split
                _less = self.__build(start_idx, p, mids, mins, _median)
                
                for i in range(self.m):
                    mids[i] = mins[i]
                mids[d] = split
                _greater = self.__build(p, end_idx, maxes, mids, _median)
                
                root = tree_buffer_root(self.tree_buffer)
                # recompute n because std::vector can
                # reallocate its internal buffer
                n = root + node_index
                # fill in entries
                n._less = _less 
                n._greater = _greater
                n.less = root + _less
                n.greater = root + _greater
                n.children = n.less.children + n.greater.children                
                n.split_dim = d
                n.split = split
            
            finally:
                PyMem_Free(mids)

            return node_index


    @cython.cdivision(True)
    cdef np.intp_t __build_compact(cKDTree self, np.intp_t start_idx, 
            np.intp_t end_idx, np.float64_t *mins, np.float64_t *maxes, 
              int _median) except -1:

        cdef: 
            ckdtreenode new_node
            np.intp_t   node_index
            np.intp_t   _less, _greater
            ckdtreenode *n
            ckdtreenode *root
            
            np.intp_t i, j, t, p, q, d
            np.float64_t size, split, minval, maxval
            np.float64_t tmp
            np.float64_t *tmp_data_point
                
        # put a new node into the node stack
        self.tree_buffer.push_back(new_node)
        node_index = self.tree_buffer.size() - 1
        root = tree_buffer_root(self.tree_buffer)
        n = root + node_index        
        
        if end_idx-start_idx <= self.leafsize:
            # below brute force limit
            # return leafnode
            n.split_dim = -1
            n.children = end_idx - start_idx
            n.start_idx = start_idx
            n.end_idx = end_idx
            return node_index
        else:
            d = 0 
            size = 0
            # Recompute hyperrectangle bounds. This should lead to a more 
            # compact kd-tree but comes at the expense of larger construction
            # time. However, construction time is usually dwarfed by the
            # query time by orders of magnitude.
            tmp_data_point = self.raw_data + self.raw_indices[start_idx]*self.m
            for i in range(self.m):
                maxes[i] = tmp_data_point[i]
                mins[i] = tmp_data_point[i]
            for j in range(start_idx+1, end_idx):
                tmp_data_point = self.raw_data + self.raw_indices[j]*self.m
                for i in range(self.m):
                    tmp = tmp_data_point[i]
                    maxes[i] = maxes[i] if (maxes[i] > tmp) else tmp
                    mins[i] = mins[i] if (mins[i] < tmp) else tmp
            # split on the dimension with largest spread        
            for i in range(self.m):
                if maxes[i]-mins[i] > size:
                    d = i
                    size = maxes[i] - mins[i]
            maxval = maxes[d]
            minval = mins[d]
            if maxval == minval:
                # all points are identical; warn user?
                # return leafnode
                n.split_dim = -1
                n.children = end_idx - start_idx
                n.start_idx = start_idx
                n.end_idx = end_idx
                return node_index
                
            # construct new inner node
            
            if _median:  
                # split on median to create a balanced tree
                # adopted from scikit-learn
                i = (end_idx-start_idx) // 2
                partition_node_indices(self.raw_data,
                                self.raw_indices + start_idx,
                                d,
                                i,
                                self.m,
                                end_idx-start_idx)               
                p = start_idx + i
                split = self.raw_data[self.raw_indices[p]*self.m+d]
             
            else:
                # split with sliding midpoint rule
                split = (maxval+minval) / 2
    
            p = start_idx
            q = end_idx - 1
            while p<=q:
                if self.raw_data[self.raw_indices[p]*self.m+d] < split:
                    p += 1
                elif self.raw_data[self.raw_indices[q]*self.m+d] >= split:
                    q -= 1
                else:
                    t = self.raw_indices[p]
                    self.raw_indices[p] = self.raw_indices[q]
                    self.raw_indices[q] = t
                    p += 1
                    q -= 1
    
            # slide midpoint if necessary
            if p == start_idx:
                # no points less than split
                j = start_idx
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx+1, end_idx):
                    if self.raw_data[self.raw_indices[i]*self.m+d] < split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[start_idx]
                self.raw_indices[start_idx] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = start_idx + 1
                q = start_idx
            elif p == end_idx:
                # no points greater than split
                j = end_idx - 1
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx, end_idx-1):
                    if self.raw_data[self.raw_indices[i]*self.m+d] > split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[end_idx-1]
                self.raw_indices[end_idx-1] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = end_idx - 1
                q = end_idx - 2

            _less = self.__build_compact(start_idx, p, mins, maxes, _median)
            _greater = self.__build_compact(p, end_idx, mins, maxes, _median)
                        
            root = tree_buffer_root(self.tree_buffer)
            # recompute n because std::vector can reallocate
            # its internal buffer
            n = root + node_index
            # fill in entries
            n._less = _less
            n._greater = _greater
            n.less = root + _less
            n.greater = root + _greater
            n.children = n.less.children + n.greater.children                
            n.split_dim = d
            n.split = split
            
            return node_index

    # -----
    # query
    # -----
    
    @cython.boundscheck(False)
    def query(cKDTree self, object x, np.intp_t k=1, np.float64_t eps=0,
              np.float64_t p=2, np.float64_t distance_upper_bound=infinity,
              np.intp_t n_jobs=1):
        """
        query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf, n_jobs=1)

        Query the kd-tree for nearest neighbors

        Parameters
        ----------
        x : array_like, last dimension self.m
            An array of points to query.
        k : integer
            The number of nearest neighbors to return.
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
            If x has shape tuple+(self.m,), then d has shape tuple+(k,).
            Missing neighbors are indicated with infinite distances.
        i : ndarray of ints
            The locations of the neighbors in self.data.
            If `x` has shape tuple+(self.m,), then `i` has shape tuple+(k,).
            Missing neighbors are indicated with self.n.

        """
        
        cdef:
            np.ndarray[np.intp_t, ndim=2] ii
            np.ndarray[np.float64_t, ndim=2] dd
            np.ndarray[np.float64_t, ndim=2] xx      
            np.intp_t c, n, i, j, CHUNK
        
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
        retshape = np.shape(x)[:-1]
        n = <np.intp_t> np.prod(retshape)
        xx = np.ascontiguousarray(x_arr).reshape(n, self.m)
        dd = np.empty((n,k),dtype=np.float64)
        dd.fill(infinity)
        ii = np.empty((n,k),dtype=np.intp)
        ii.fill(self.n)
        
        # Do the query in an external C++ function. 
        # The GIL will be released in the external query function.
        
        if (n_jobs == -1): 
            n_jobs = number_of_processors
        
        if n_jobs > 1:
            # static scheduling without load balancing is good enough
            CHUNK = n//n_jobs if n//n_jobs else n
                             
            def _thread_func(self, _dd, _ii, _xx, _j, n, CHUNK, p, k, eps, dub):
                cdef: 
                    np.intp_t j = _j
                    np.ndarray[np.intp_t,ndim=2] ii = _ii
                    np.ndarray[np.float64_t,ndim=2] dd = _dd
                    np.ndarray[np.float64_t,ndim=2] xx = _xx
                    np.intp_t start = j*CHUNK
                    np.intp_t stop = start + CHUNK
                stop = n if stop > n else stop
                if start < n:
                    query_knn(<ckdtree*>self, &dd[start,0], &ii[start,0], 
                        &xx[start,0], stop-start, k, eps, p, dub)
            
            # There might be n_jobs+1 threads spawned here, but only n_jobs of 
            # them will do significant work.
            threads = [threading.Thread(target=_thread_func,
                                        args=(self, dd, ii, xx, j, n, CHUNK, p,
                                              k, eps, distance_upper_bound))
                             for j in range(1+(n//CHUNK))]
            # Set the daemon flag so the process can be aborted, 
            # start all threads and wait for completion.
            for t in threads:
                t.daemon = True
                t.start()
            for t in threads: 
                t.join()
        else:
            query_knn(<ckdtree*>self, &dd[0,0], &ii[0,0], &xx[0,0], 
                n, k, eps, p, distance_upper_bound)
                
        if single:
            if k == 1:
                if sizeof(long) < sizeof(np.intp_t):
                    # ... e.g. Windows 64
                    if ii[0,0] <= <np.intp_t>LONG_MAX:
                        return dd[0,0], int(ii[0,0])
                    else:
                        return dd[0,0], ii[0,0]
                else:
                    # ... most other platforms
                    return dd[0,0], ii[0,0]
            else:
                return dd[0], ii[0]
        else:
            if sizeof(long) < sizeof(np.intp_t):
                # ... e.g. Windows 64
                for i in range(n):
                    for j in range(k):
                        if ii[i,j] > <np.intp_t>LONG_MAX:
                            # C long overlow, return array of dtype=np.int_p
                            if k==1:
                                return np.reshape(dd[...,0],retshape), np.reshape(ii[...,0],retshape)
                            else:
                                return np.reshape(dd,retshape+(k,)), np.reshape(ii,retshape+(k,))

                # no C long overlow, return array of dtype=int
                if k==1:
                    return np.reshape(dd[...,0],retshape), np.reshape(ii[...,0],retshape).astype(int)
                else:
                    return np.reshape(dd,retshape+(k,)), np.reshape(ii,retshape+(k,)).astype(int)     

            else:
                # ... most other platforms
                if k==1:
                    return np.reshape(dd[...,0],retshape), np.reshape(ii[...,0],retshape)
                else:
                    return np.reshape(dd,retshape+(k,)), np.reshape(ii,retshape+(k,))

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
                    tmp = m * [None]
                    cur = npy_intp_vector_buf(vvres[i])
                    for j in range(m):
                        tmp[j] = cur[0]
                        cur += 1
                    result[c] = sorted(tmp)
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
                tmp = m * [None]
                cur = npy_intp_vector_buf(vvres[i]) 
                for j in range(m):
                    tmp[j] = cur[0]
                    cur += 1
                results[i] = sorted(tmp)       
                                  
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
                    np.float64_t eps=0, object output_type=set):
        """
        query_pairs(self, r, p=2., eps=0)

        Find all pairs of points whose distance is at most r.

        Parameters
        ----------
        r : positive float
            The maximum distance.
        p : float, optional
            Which Minkowski norm to use.  `p` has to meet the condition
            ``1 <= p <= infinity``.
        eps : float, optional
            Approximate search.  Branches of the tree are not explored
            if their nearest points are further than ``r/(1+eps)``, and
            branches are added in bulk if their furthest points are nearer
            than ``r * (1+eps)``.  `eps` has to be non-negative.
        output_type : type, optional
            Choose the output container, set or ndarray. Default: set

        Returns
        -------
        results : set or ndarray
            Set of pairs ``(i,j)``, with ``i < j``, for which the corresponding
            positions are close. If output_type is ndarray, an ndarry is 
            returned instead of a set.

        """
                 
        cdef: 
            vector[ordered_pair] *vres
            set results
            np.intp_t i, n
            ordered_pair *pair
            np.ndarray array_res
        
        if output_type not in (set, np.ndarray):
            raise ValueError("output type must be set or ndarray")
        
        vres = NULL
        results = set()
        try:
            vres = new vector[ordered_pair]()
            query_pairs(<ckdtree*> self, r, p, eps, vres)
            n = vres.size()
            pair = ordered_pair_vector_buf(vres)
            
            if output_type == set:
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
            else:
                array_res = np.zeros((n,2), dtype=np.intp)
                string.memcpy(np.PyArray_DATA(array_res), <void*> pair, 
                                  n*2*sizeof(np.intp_t))
        finally:
            if vres != NULL:
                del vres
        
        if output_type == set:
            return results
        else:           
            return array_res


    # ---------------
    # count_neighbors
    # ---------------

    @cython.boundscheck(False)
    def count_neighbors(cKDTree self, cKDTree other, object r, np.float64_t p=2.):
        """
        count_neighbors(self, other, r, p=2.)

        Count how many nearby pairs can be formed.

        Count the number of pairs (x1,x2) can be formed, with x1 drawn
        from self and x2 drawn from `other`, and where
        ``distance(x1, x2, p) <= r``.
        This is the "two-point correlation" described in Gray and Moore 2000,
        "N-body problems in statistical learning", and the code here is based
        on their algorithm.

        Parameters
        ----------
        other : KDTree instance
            The other tree to draw points from.
        r : float or one-dimensional array of floats
            The radius to produce a count for. Multiple radii are searched with
            a single tree traversal.
        p : float, optional 
            1<=p<=infinity, default 2.0
            Which Minkowski p-norm to use

        Returns
        -------
        result : int or 1-D array of ints
            The number of pairs. Note that this is internally stored in a numpy int,
            and so may overflow if very large (2e9).

        """
        cdef: 
            int r_ndim
            np.intp_t n_queries, i
            np.ndarray[np.float64_t, ndim=1, mode="c"] real_r
            np.ndarray[np.intp_t, ndim=1, mode="c"] results, idx

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
        n_queries = real_r.shape[0]

        # Internally, we represent all distances as distance ** p
        if p != infinity:
            for i in range(n_queries):
                if real_r[i] != infinity:
                    real_r[i] = real_r[i] ** p

        results = np.zeros(n_queries, dtype=np.intp)
        idx = np.arange(n_queries, dtype=np.intp)
        
        count_neighbors(<ckdtree*> self, <ckdtree*> other, n_queries,
                        &real_r[0], &results[0], &idx[0], p)
        
        if r_ndim == 0:
            if results[0] <= <np.intp_t> LONG_MAX:
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
            'coo_matrix', 'dict', or 'recarray'. Default: 'dok_matrix'.

        Returns
        -------
        result : dok_matrix, coo_matrix, dict or ndarray
            Sparse matrix representing the results in "dictionary of keys" 
            format. If a dict is returned the keys are (i,j) tuples of indices.
            If output_type is 'recarray' a record array with fields 'i', 'j',
            and 'k' is returned,
        """
        
        cdef:
            vector[np.intp_t] *res_i
            vector[np.intp_t] *res_j
            vector[np.float64_t] *res_v
            np.ndarray res_arr
            dict res_dict
            np.intp_t i, j, k, n, s
            np.intp_t *pi 
            np.intp_t *pj
            np.intp_t *pvi 
            np.intp_t *pvj
            char *cur
            np.float64_t v
            np.float64_t *pv
            np.float64_t *pvv
            
        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to sparse_distance_matrix have "
                             "different dimensionality")
                             
        # check output type:
        if not output_type in ('dok_matrix', 'coo_matrix', 'dict', 'recarray'):
            raise ValueError('Invalid output_type')
                                     
        res_i = NULL
        res_j = NULL
        res_v = NULL
        
        try:
        
            res_i = new vector[np.intp_t]()
            res_j = new vector[np.intp_t]()
            res_v = new vector[np.float64_t]()
            
            sparse_distance_matrix(
                <ckdtree*> self, <ckdtree*> other, p, max_distance, 
                res_i, res_j, res_v)
                
            if output_type == 'dict':
                res_dict = dict()
                n = <np.intp_t> (res_i.size())
                pvi = npy_intp_vector_buf(res_i)
                pvj = npy_intp_vector_buf(res_j)
                pvv = npy_float64_vector_buf(res_v)
                for k in range(n):
                    i = pvi[k]
                    j = pvj[k]
                    v = pvv[k]
                    res_dict[(i,j)] = v
                
            else:
                n = <np.intp_t> (res_i.size())
                res_dtype = np.dtype(
                   [('i', np.intp), 
                    ('j', np.intp), 
                    ('v',np.float64)])
                res_arr = np.empty(n, dtype=res_dtype)
                s = res_arr.strides[0]
                cur = <char*> np.PyArray_DATA(res_arr)
                pvi = npy_intp_vector_buf(res_i)
                pvj = npy_intp_vector_buf(res_j)
                pvv = npy_float64_vector_buf(res_v)
                for k in range(n):
                    pi = <np.intp_t*> cur
                    pj = <np.intp_t*> (cur + sizeof(np.intp_t))
                    pv = <np.float64_t*> (cur + 2*sizeof(np.intp_t))
                    i = pvi[k]
                    j = pvj[k]
                    v = pvv[k]
                    pi[0] = i
                    pj[0] = j
                    pv[0] = v 
                    cur += s
                result = res_arr
                
                if output_type != 'recarray':
                    result = scipy.sparse.coo_matrix(
                       (res_arr['v'], (res_arr['i'], res_arr['j'])),
                                       shape=(self.n,other.n))
                else:
                    result = res_arr
        
        finally:
        
            if res_i != NULL: 
                del res_i
            
            if res_j != NULL: 
                del res_j
                
            if res_v != NULL: 
                del res_v
                
                
        if output_type != 'dok_matrix':
            return result
        else:
            return result.todok()
                             
         
    # ----------------------
    # pickle
    # ----------------------    

        
    def __reduce__(self):
        return (new_object, (cKDTree,), self.__getstate__())

    def __getstate__(cKDTree self):
        cdef object state
        cdef object tree = pickle_tree_buffer(self.tree_buffer)
        state = (tree, self.data.copy(), self.n, self.m, self.leafsize,
                      self.maxes, self.mins, self.indices.copy())
        return state
            
    def __setstate__(cKDTree self, state):
        cdef object tree
        self.tree_buffer = new vector[ckdtreenode]()
        
        # unpack the state
        (tree, self.data, self.n, self.m, self.leafsize, 
            self.maxes, self.mins, self.indices) = state
        
        # copy kd-tree buffer 
        unpickle_tree_buffer(self.tree_buffer, tree)    
        
        # set raw pointers
        self.raw_data = <np.float64_t*>np.PyArray_DATA(self.data)
        self.raw_maxes = <np.float64_t*>np.PyArray_DATA(self.maxes)
        self.raw_mins = <np.float64_t*>np.PyArray_DATA(self.mins)
        self.raw_indices = <np.intp_t*>np.PyArray_DATA(self.indices)
        
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

