# Copyright Anne M. Archibald 2008
# Released under the scipy license
import numpy as np
cimport numpy as np
cimport libc.stdlib as stdlib
cimport cython

cdef extern from "limits.h":
    long LONG_MAX

import scipy.sparse

import kdtree

cdef np.float64_t infinity = np.inf


__all__ = ['cKDTree']


# This is needed because a Python int and np.int is a C long on Python 2.x,
# which can be 32 bits on 64-bit systems (e.g. Windows).

cdef object npy_intp_dtype

if sizeof(np.npy_intp) == sizeof(np.int32_t):
    npy_intp_dtype = np.int32
elif sizeof(np.npy_intp) == sizeof(np.int64_t):    
    npy_intp_dtype = np.int64
else:
    raise ImportError, 'Unexpected length of npy_intp'    


# priority queue

cdef union heapcontents:    # FIXME: Unions are not always portable, verify this 
    np.npy_intp intdata     # union is never used in an ABI dependent way.
    char* ptrdata

cdef struct heapitem:
    np.float64_t priority
    heapcontents contents

cdef struct heap:
    np.npy_intp n
    heapitem* heap
    np.npy_intp space

cdef inline int heapcreate(heap* self, np.npy_intp initial_size) except -1:
    cdef void *tmp
    self.space = initial_size
    self.heap = <heapitem*> NULL
    tmp = stdlib.malloc(sizeof(heapitem)*self.space)
    if tmp == NULL:
        raise MemoryError
    self.heap = <heapitem*> tmp  
    self.n = 0
    return 0

cdef inline int heapdestroy(heap* self) except -1:
    if self.heap != <heapitem*> NULL:
        stdlib.free(self.heap)
    return 0

cdef inline int heapresize(heap* self, np.npy_intp new_space) except -1:
    cdef void *tmp
    if new_space<self.n:
        raise ValueError("Heap containing %d items cannot be resized to %d" % (int(self.n), int(new_space)))
    self.space = new_space
    self.heap = <heapitem*> NULL
    tmp = stdlib.realloc(<void*>self.heap, new_space*sizeof(heapitem))
    if tmp == NULL:
        raise MemoryError
    self.heap = <heapitem*> tmp
    return 0


@cython.cdivision(True)
cdef inline int heappush(heap* self, heapitem item) except -1:
    cdef np.npy_intp i
    cdef heapitem t
    self.n += 1
    if self.n>self.space:
        heapresize(self,2*self.space+1)
    i = self.n-1
    self.heap[i] = item
    while i>0 and self.heap[i].priority<self.heap[(i-1)//2].priority:
        t = self.heap[(i-1)//2]
        self.heap[(i-1)//2] = self.heap[i]
        self.heap[i] = t
        i = (i-1)//2
    return 0


cdef heapitem heappeek(heap* self):
    return self.heap[0]


@cython.cdivision(True)
cdef int heapremove(heap* self) except -1:
    cdef heapitem t
    cdef np.npy_intp i, j, k, l

    self.heap[0] = self.heap[self.n-1]
    self.n -= 1
    if self.n < self.space//4 and self.space>40: #FIXME: magic number
        heapresize(self,self.space//2+1)
    i=0
    j=1
    k=2
    while ((j<self.n and 
                self.heap[i].priority > self.heap[j].priority or
            k<self.n and 
                self.heap[i].priority > self.heap[k].priority)):
        if k<self.n and self.heap[j].priority>self.heap[k].priority:
            l = k
        else:
            l = j
        t = self.heap[l]
        self.heap[l] = self.heap[i]
        self.heap[i] = t
        i = l
        j = 2*i+1
        k = 2*i+2
    return 0


cdef int heappop(heap* self, heapitem *it) except -1:
    it[0] = heappeek(self)
    heapremove(self)
    return 0


# utility functions
cdef inline np.float64_t dmax(np.float64_t x, np.float64_t y):
    if x>y:
        return x
    else:
        return y
        
cdef inline np.float64_t dabs(np.float64_t x):
    if x>0:
        return x
    else:
        return -x
        
cdef inline np.float64_t _distance_p(np.float64_t *x, np.float64_t *y,
                                     np.float64_t p, np.npy_intp k,
                                     np.float64_t upperbound):
    """Compute the distance between x and y

    Computes the Minkowski p-distance to the power p between two points.
    If the distance**p is larger than upperbound, then any number larger
    than upperbound may be returned (the calculation is truncated).
    """
    cdef np.npy_intp i
    cdef np.float64_t r, z
    r = 0
    if p==2:
        for i in range(k):
            z = x[i] - y[i]
            r += z*z
            if r>upperbound:
                return r 
    elif p==infinity:
        for i in range(k):
            r = dmax(r,dabs(x[i]-y[i]))
            if r>upperbound:
                return r
    elif p==1:
        for i in range(k):
            r += dabs(x[i]-y[i])
            if r>upperbound:
                return r
    else:
        for i in range(k):
            r += dabs(x[i]-y[i])**p
            if r>upperbound:
                return r
    return r


# Interval arithmetic
cdef struct Rectangle:
    np.npy_intp m
    np.float64_t *mins
    np.float64_t *maxes

# 1-d pieces
# These should only be used if p != infinity
cdef inline np.float64_t min_dist_point_interval_p(np.float64_t* x,
                                                   Rectangle rect,
                                                   np.npy_intp k,
                                                   np.float64_t p):    
    """Compute the minimum distance along dimension k between x and
    a point in the hyperrectangle.
    """
    return dmax(0, dmax(rect.mins[k] - x[k], x[k] - rect.maxes[k])) ** p

cdef inline np.float64_t max_dist_point_interval_p(np.float64_t* x,
                                                   Rectangle rect,
                                                   np.npy_intp k,
                                                   np.float64_t p):
    """Compute the maximum distance along dimension k between x and
    a point in the hyperrectangle.
    """
    return dmax(rect.maxes[k] - x[k], x[k] - rect.mins[k]) ** p

cdef inline np.float64_t min_dist_interval_interval_p(Rectangle rect1,
                                                      Rectangle rect2,
                                                      np.npy_intp k,
                                                      np.float64_t p):
    """Compute the minimum distance along dimension k between points in
    two hyperrectangles.
    """
    return dmax(0, dmax(rect1.mins[k] - rect2.maxes[k],
                        rect2.mins[k] - rect1.maxes[k])) ** p

cdef inline np.float64_t max_dist_interval_interval_p(Rectangle rect1,
                                                      Rectangle rect2,
                                                      np.npy_intp k,
                                                      np.float64_t p):
    """Compute the maximum distance along dimension k between points in
    two hyperrectangles.
    """
    return dmax(rect1.maxes[k] - rect2.mins[k], rect2.maxes[k] - rect1.mins[k]) ** p

# Interval arithmetic in m-D

# These should be used only for p == infinity
cdef inline np.float64_t min_dist_point_rect_p_inf(np.float64_t* x,
                                                   Rectangle rect):
    """Compute the minimum distance between x and the given hyperrectangle."""
    cdef np.npy_intp i
    cdef np.float64_t min_dist = 0.
    for i in range(rect.m):
        min_dist = dmax(min_dist, dmax(rect.mins[i]-x[i], x[i]-rect.maxes[i]))
    return min_dist

cdef inline np.float64_t max_dist_point_rect_p_inf(np.float64_t* x,
                                                   Rectangle rect):
    """Compute the maximum distance between x and the given hyperrectangle."""
    cdef np.npy_intp i
    cdef np.float64_t max_dist = 0.
    for i in range(rect.m):
        max_dist = dmax(max_dist, dmax(rect.maxes[i]-x[i], x[i]-rect.mins[i]))
    return max_dist

cdef inline np.float64_t min_dist_rect_rect_p_inf(Rectangle rect1,
                                                  Rectangle rect2):
    """Compute the minimum distance between points in two hyperrectangles."""
    cdef np.npy_intp i
    cdef np.float64_t min_dist = 0.
    for i in range(rect1.m):
        min_dist = dmax(min_dist, dmax(rect1.mins[i] - rect2.maxes[i],
                                       rect2.mins[i] - rect1.maxes[i]))
    return min_dist

cdef inline np.float64_t max_dist_rect_rect_p_inf(Rectangle rect1,
                                                  Rectangle rect2):
    """Compute the maximum distance between points in two hyperrectangles."""
    cdef np.npy_intp i
    cdef np.float64_t max_dist = 0.
    for i in range(rect1.m):
        max_dist = dmax(max_dist, dmax(rect1.maxes[i] - rect2.mins[i],
                                       rect2.maxes[i] - rect1.mins[i]))
    return max_dist

# A pair of functions to do incremental updates of min and max distances
# between two hyperrectangles
cdef inline void __rect_preupdate(Rectangle rect1, Rectangle rect2,
                                  np.npy_intp k, np.float64_t p,
                                  np.float64_t min_distance,
                                  np.float64_t max_distance,
                                  np.float64_t *part_min,
                                  np.float64_t *part_max):
    if p != infinity:
        part_min[0] = min_distance - min_dist_interval_interval_p(rect1, rect2, k, p)
        part_max[0] = max_distance - max_dist_interval_interval_p(rect1, rect2, k, p)

cdef inline void __rect_postupdate(Rectangle rect1, Rectangle rect2,
                                   np.npy_intp k, np.float64_t p,
                                   np.float64_t *min_distance,
                                   np.float64_t *max_distance,
                                   np.float64_t part_min,
                                   np.float64_t part_max):
    if p != infinity:
        min_distance[0] = part_min + min_dist_interval_interval_p(rect1, rect2, k, p)
        max_distance[0] = part_max + max_dist_interval_interval_p(rect1, rect2, k, p)
    else:
        min_distance[0] = min_dist_rect_rect_p_inf(rect1, rect2)
        max_distance[0] = max_dist_rect_rect_p_inf(rect1, rect2)


# Tree structure

cdef struct innernode:
    np.npy_intp split_dim
    np.npy_intp children
    np.float64_t split
    innernode* less
    innernode* greater
    
cdef struct leafnode:
    np.npy_intp split_dim
    np.npy_intp children
    np.npy_intp start_idx
    np.npy_intp end_idx


# this is the standard trick for variable-size arrays:
# malloc sizeof(nodeinfo)+self.m*sizeof(np.float64_t) bytes.

cdef struct nodeinfo:
    innernode* node
    np.float64_t side_distances[0]  # FIXME: Only valid in C99, invalid C++ and C89

# Utility for building a coo matrix incrementally
cdef class coo_entries:
    cdef:
        np.npy_intp n, n_max
        np.ndarray i, j
        np.ndarray v
        np.npy_intp  *i_data, *j_data
        np.float64_t *v_data
    
    def __init__(self):
        self.n = 0
        self.n_max = 10
        self.i = np.empty(self.n_max, dtype=npy_intp_dtype)
        self.j = np.empty(self.n_max, dtype=npy_intp_dtype)
        self.v = np.empty(self.n_max, dtype=np.float64)
        self.i_data = <np.npy_intp *>np.PyArray_DATA(self.i)
        self.j_data = <np.npy_intp *>np.PyArray_DATA(self.j)
        self.v_data = <np.float64_t*>np.PyArray_DATA(self.v)

    cdef void add(coo_entries self, np.npy_intp i, np.npy_intp j, np.float64_t v):
        cdef np.npy_intp k
        if self.n == self.n_max:
            self.n_max *= 2
            self.i.resize(self.n_max)
            self.j.resize(self.n_max)
            self.v.resize(self.n_max)
            self.i_data = <np.npy_intp *>np.PyArray_DATA(self.i)
            self.j_data = <np.npy_intp *>np.PyArray_DATA(self.j)
            self.v_data = <np.float64_t*>np.PyArray_DATA(self.v)
        k = self.n
        self.i_data[k] = i
        self.j_data[k] = j
        self.v_data[k] = v
        self.n += 1

    def to_matrix(coo_entries self, shape=None):
        # Shrink arrays to size
        self.i.resize(self.n)
        self.j.resize(self.n)
        self.v.resize(self.n)
        self.i_data = <np.npy_intp *>np.PyArray_DATA(self.i)
        self.j_data = <np.npy_intp *>np.PyArray_DATA(self.j)
        self.v_data = <np.float64_t*>np.PyArray_DATA(self.v)
        self.n_max = self.n
        return scipy.sparse.coo_matrix((self.v, (self.i, self.j)),
                                       shape=shape)


cdef class cKDTree:
    """
    cKDTree(data, int leafsize=10)

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
    data : array-like, shape (n,m)
        The n data points of dimension mto be indexed. This array is 
        not copied unless this is necessary to produce a contiguous 
        array of doubles, and so modifying this data will result in 
        bogus results.
    leafsize : positive integer
        The number of points at which the algorithm switches over to
        brute-force.

    """

    cdef innernode* tree 
    cdef readonly object data
    cdef np.float64_t* raw_data
    cdef readonly np.npy_intp n, m
    cdef readonly np.npy_intp leafsize
    cdef readonly object maxes
    cdef np.float64_t* raw_maxes
    cdef readonly object mins
    cdef np.float64_t* raw_mins
    cdef object indices
    cdef np.npy_intp* raw_indices

    def __init__(cKDTree self, data, np.npy_intp leafsize=10):
        cdef np.ndarray[np.float64_t, ndim=2] inner_data
        cdef np.ndarray[np.float64_t, ndim=1] inner_maxes
        cdef np.ndarray[np.float64_t, ndim=1] inner_mins
        cdef np.ndarray[np.npy_intp, ndim=1] inner_indices
        self.data = np.ascontiguousarray(data,dtype=np.float64)
        self.n, self.m = np.shape(self.data)
        self.leafsize = leafsize
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.maxes = np.ascontiguousarray(np.amax(self.data,axis=0), dtype=np.float64)
        self.mins = np.ascontiguousarray(np.amin(self.data,axis=0), dtype=np.float64)
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=npy_intp_dtype))

        inner_data = self.data
        self.raw_data = <np.float64_t*>np.PyArray_DATA(inner_data)
        inner_maxes = self.maxes
        self.raw_maxes = <np.float64_t*>np.PyArray_DATA(inner_maxes)
        inner_mins = self.mins
        self.raw_mins = <np.float64_t*>np.PyArray_DATA(inner_mins)
        inner_indices = self.indices
        self.raw_indices = <np.npy_intp*>np.PyArray_DATA(inner_indices)

        self.tree = self.__build(0, self.n, self.raw_maxes, self.raw_mins)

    cdef innernode* __build(cKDTree self, np.npy_intp start_idx, np.npy_intp end_idx,
                            np.float64_t* maxes, np.float64_t* mins) except? <innernode*> NULL:
        cdef leafnode* n
        cdef innernode* ni
        cdef np.npy_intp i, j, t, p, q, d
        cdef np.float64_t size, split, minval, maxval
        cdef np.float64_t*mids
        if end_idx-start_idx<=self.leafsize:
            n = <leafnode*>stdlib.malloc(sizeof(leafnode))
            if n == <leafnode*> NULL: 
                raise MemoryError
            n.split_dim = -1
            n.children = end_idx - start_idx
            n.start_idx = start_idx
            n.end_idx = end_idx
            return <innernode*>n
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
                n = <leafnode*>stdlib.malloc(sizeof(leafnode))
                if n == <leafnode*> NULL: 
                    raise MemoryError
                n.split_dim = -1
                n.children = end_idx - start_idx
                n.start_idx = start_idx
                n.end_idx = end_idx
                return <innernode*>n

            split = (maxval+minval)/2

            p = start_idx
            q = end_idx-1
            while p<=q:
                if self.raw_data[self.raw_indices[p]*self.m+d]<split:
                    p+=1
                elif self.raw_data[self.raw_indices[q]*self.m+d]>=split:
                    q-=1
                else:
                    t = self.raw_indices[p]
                    self.raw_indices[p] = self.raw_indices[q]
                    self.raw_indices[q] = t
                    p+=1
                    q-=1

            # slide midpoint if necessary
            if p==start_idx:
                # no points less than split
                j = start_idx
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx+1, end_idx):
                    if self.raw_data[self.raw_indices[i]*self.m+d]<split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[start_idx]
                self.raw_indices[start_idx] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = start_idx+1
                q = start_idx
            elif p==end_idx:
                # no points greater than split
                j = end_idx-1
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx, end_idx-1):
                    if self.raw_data[self.raw_indices[i]*self.m+d]>split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[end_idx-1]
                self.raw_indices[end_idx-1] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = end_idx-1
                q = end_idx-2

            # construct new node representation
            ni = <innernode*>stdlib.malloc(sizeof(innernode))
            if ni ==  <innernode*> NULL:
                raise MemoryError

            try:
                mids = <np.float64_t*>stdlib.malloc(sizeof(np.float64_t)*self.m)
                if mids == <np.float64_t*> NULL:
                    raise MemoryError
                        
                for i in range(self.m):
                    mids[i] = maxes[i]
                mids[d] = split
                ni.less = self.__build(start_idx,p,mids,mins)

                for i in range(self.m):
                    mids[i] = mins[i]
                mids[d] = split
                ni.greater = self.__build(p,end_idx,maxes,mids)

                ni.children = ni.less.children + ni.greater.children
            
            except:
                # free ni if it cannot be returned
                if ni !=  <innernode*> NULL:
                    stdlib.free(mids)
                if mids != <np.float64_t*> NULL:
                    stdlib.free(mids)
                raise
            else:
                if mids != <np.float64_t*> NULL:
                    stdlib.free(mids)

            ni.split_dim = d
            ni.split = split
            return ni
                    
    cdef __free_tree(cKDTree self, innernode* node):
        if node.split_dim!=-1:
            self.__free_tree(node.less)
            self.__free_tree(node.greater)
        stdlib.free(node)

    def __dealloc__(cKDTree self):
        if <np.npy_intp>(self.tree) == 0:
            # should happen only if __init__ was never called
            return
        self.__free_tree(self.tree)

    cdef int __query(cKDTree self, 
            np.float64_t*result_distances, 
            np.npy_intp*result_indices, 
            np.float64_t*x, 
            np.npy_intp k, 
            np.float64_t eps, 
            np.float64_t p, 
            np.float64_t distance_upper_bound) except -1:

        cdef heap q
        cdef heap neighbors

        cdef np.npy_intp i, j
        cdef np.float64_t t
        cdef nodeinfo* inf
        cdef nodeinfo* inf2
        cdef np.float64_t d
        cdef np.float64_t epsfac
        cdef np.float64_t min_distance
        cdef np.float64_t far_min_distance
        cdef heapitem it, it2, neighbor
        cdef leafnode* node
        cdef innernode* inode
        cdef innernode* near
        cdef innernode* far
        cdef np.float64_t* side_distances


        inf = inf2 = <nodeinfo*> NULL    

        try:
            # priority queue for chasing nodes
            # entries are:
            #  minimum distance between the cell and the target
            #  distances between the nearest side of the cell and the target
            #  the head node of the cell
            heapcreate(&q,12)

            # priority queue for the nearest neighbors
            # furthest known neighbor first
            # entries are (-distance**p, i)
            heapcreate(&neighbors,k)

            # set up first nodeinfo
            inf = <nodeinfo*>stdlib.malloc(sizeof(nodeinfo)+self.m*sizeof(np.float64_t))
            if inf == <nodeinfo*> NULL:
                raise MemoryError
            inf.node = self.tree
            for i in range(self.m):
                inf.side_distances[i] = 0
                t = x[i]-self.raw_maxes[i]
                if t>inf.side_distances[i]:
                    inf.side_distances[i] = t
                else:
                    t = self.raw_mins[i]-x[i]
                    if t>inf.side_distances[i]:
                        inf.side_distances[i] = t
                if p!=1 and p!=infinity:
                    inf.side_distances[i]=inf.side_distances[i]**p

            # compute first distance
            min_distance = 0.
            for i in range(self.m):
                if p==infinity:
                    min_distance = dmax(min_distance,inf.side_distances[i])
                else:
                    min_distance += inf.side_distances[i]

            # fiddle approximation factor
            if eps==0:
                epsfac=1
            elif p==infinity:
                epsfac = 1/(1+eps)
            else:
                epsfac = 1/(1+eps)**p

            # internally we represent all distances as distance**p
            if p!=infinity and distance_upper_bound!=infinity:
                distance_upper_bound = distance_upper_bound**p

            while True:
                if inf.node.split_dim==-1:
                    node = <leafnode*>inf.node

                    # brute-force
                    for i in range(node.start_idx,node.end_idx):
                        d = _distance_p(
                                self.raw_data+self.raw_indices[i]*self.m,
                                x,p,self.m,distance_upper_bound)
                            
                        if d<distance_upper_bound:
                            # replace furthest neighbor
                            if neighbors.n==k:
                                heapremove(&neighbors)
                            neighbor.priority = -d
                            neighbor.contents.intdata = self.raw_indices[i]
                            heappush(&neighbors,neighbor)

                            # adjust upper bound for efficiency
                            if neighbors.n==k:
                                distance_upper_bound = -heappeek(&neighbors).priority
                    
                    # done with this node, get another
                    stdlib.free(inf)
                    inf = <nodeinfo*> NULL

                    if q.n==0:
                        # no more nodes to visit
                        break
                    else:
                        heappop(&q, &it)
                        inf = <nodeinfo*>it.contents.ptrdata
                        min_distance = it.priority
                else:
                    inode = <innernode*>inf.node

                    # we don't push cells that are too far onto the queue at all,
                    # but since the distance_upper_bound decreases, we might get 
                    # here even if the cell's too far
                    if min_distance>distance_upper_bound*epsfac:

                        # since this is the nearest cell, we're done, bail out
                        stdlib.free(inf)
                        inf = <nodeinfo*> NULL

                        # free all the nodes still on the heap
                        for i in range(q.n):
                            stdlib.free(q.heap[i].contents.ptrdata)
                            q.heap[i].contents.ptrdata = <char*> NULL
                        break

                    # set up children for searching
                    if x[inode.split_dim]<inode.split:
                        near = inode.less
                        far = inode.greater
                    else:
                        near = inode.greater
                        far = inode.less

                    # near child is at the same distance as the current node
                    # we're going here next, so no point pushing it on the queue
                    # no need to recompute the distance or the side_distances
                    inf.node = near

                    # far child is further by an amount depending only
                    # on the split value; compute its distance and side_distances
                    # and push it on the queue if it's near enough
                    inf2 = <nodeinfo*>stdlib.malloc(sizeof(nodeinfo)+self.m*sizeof(np.float64_t))
                    if inf2 == <nodeinfo*> NULL:
                        raise MemoryError
            
                    it2.contents.ptrdata = <char*> inf2
                    inf2.node = far
                    # most side distances unchanged
                    for i in range(self.m):
                        inf2.side_distances[i] = inf.side_distances[i]

                    # one side distance changes
                    # we can adjust the minimum distance without recomputing
                    if p == infinity:
                        # we never use side_distances in the l_infinity case
                        # inf2.side_distances[inode.split_dim] = dabs(inode.split-x[inode.split_dim])
                        far_min_distance = dmax(min_distance, dabs(inode.split-x[inode.split_dim]))
                    elif p == 1:
                        inf2.side_distances[inode.split_dim] = dabs(inode.split-x[inode.split_dim])
                        far_min_distance = min_distance - \
                                           inf.side_distances[inode.split_dim] + \
                                           inf2.side_distances[inode.split_dim]
                    else:
                        inf2.side_distances[inode.split_dim] = dabs(inode.split - 
                                                                    x[inode.split_dim])**p
                        far_min_distance = min_distance - \
                                           inf.side_distances[inode.split_dim] + \
                                           inf2.side_distances[inode.split_dim]

                    it2.priority = far_min_distance


                    # far child might be too far, if so, don't bother pushing it
                    if far_min_distance<=distance_upper_bound*epsfac:
                        heappush(&q,it2)
                    else:
                        stdlib.free(inf2)
                        inf2 = <nodeinfo*> NULL
                        # just in case
                        it2.contents.ptrdata = <char*> NULL

            # fill output arrays with sorted neighbors 
            for i in range(neighbors.n-1,-1,-1):
                heappop(&neighbors, &neighbor) # FIXME: neighbors may be realloced
                result_indices[i] = neighbor.contents.intdata
                if p==1 or p==infinity:
                    result_distances[i] = -neighbor.priority
                else:
                    result_distances[i] = (-neighbor.priority)**(1./p)

            inf = inf2 = <nodeinfo*> NULL

        finally:
            if inf2 != <nodeinfo*> NULL:
                stdlib.free(inf2)

            if inf != <nodeinfo*> NULL:
                stdlib.free(inf)
            try:
                heapdestroy(&q)
            finally:
                heapdestroy(&neighbors)

        return 0


    def query(cKDTree self, object x, np.npy_intp k=1, np.float64_t eps=0,
              np.float64_t p=2, np.float64_t distance_upper_bound=infinity):
        """query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf)
        
        Query the kd-tree for nearest neighbors

        Parameters
        ----------
        x : array_like, last dimension self.m
            An array of points to query.
        k : integer
            The number of nearest neighbors to return.
        eps : non-negative float
            Return approximate nearest neighbors; the kth returned value 
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
        cdef np.ndarray[np.npy_intp, ndim=2] ii
        cdef np.ndarray[np.float64_t, ndim=2] dd
        cdef np.ndarray[np.float64_t, ndim=2] xx
        cdef np.npy_intp c
        x = np.asarray(x).astype(np.float64)
        if np.shape(x)[-1] != self.m:
            raise ValueError("x must consist of vectors of length %d but has"
                             "shape %s" % (int(self.m), np.shape(x)))
        if p < 1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        if len(x.shape)==1:
            single = True
            x = x[np.newaxis,:]
        else:
            single = False
        retshape = np.shape(x)[:-1]
        n = np.prod(retshape)
        xx = np.reshape(x,(n,self.m))
        xx = np.ascontiguousarray(xx,dtype=np.float64)
        dd = np.empty((n,k),dtype=np.float64)
        dd.fill(infinity)
        ii = np.empty((n,k),dtype=npy_intp_dtype)
        ii.fill(self.n)
        for c in range(n):
            self.__query(
                    (<np.float64_t*>np.PyArray_DATA(dd))+c*k,
                    (<np.npy_intp*>np.PyArray_DATA(ii))+c*k,
                    (<np.float64_t*>np.PyArray_DATA(xx))+c*self.m, 
                    k, 
                    eps,
                    p, 
                    distance_upper_bound)
        if single:
            if k==1:
                if sizeof(long) < sizeof(np.npy_intp):
                    # ... e.g. Windows 64
                    if ii[0,0] <= <np.npy_intp>LONG_MAX:
                        return dd[0,0], int(ii[0,0])
                    else:
                        return dd[0,0], ii[0,0]
                else:
                    # ... most other platforms
                    return dd[0,0], ii[0,0]
            else:
                return dd[0], ii[0]
        else:
            if k==1:
                return np.reshape(dd[...,0],retshape), np.reshape(ii[...,0],retshape)
            else:
                return np.reshape(dd,retshape+(k,)), np.reshape(ii,retshape+(k,))
            
    # ----------------
    # query_ball_point
    # ----------------
    cdef int __query_ball_point_traverse_no_checking(cKDTree self,
                                                      list results,
                                                      innernode* node) except -1:
        cdef leafnode* lnode
        cdef innernode* inode
        cdef np.npy_intp i

        if node.split_dim == -1:  # leaf node
            lnode = <leafnode*> node
            for i in range(lnode.start_idx, lnode.end_idx):
                results.append(self.raw_indices[i])
        else:
            inode = <innernode*>node
            self.__query_ball_point_traverse_no_checking(results, inode.less)
            self.__query_ball_point_traverse_no_checking(results, inode.greater)

        return 0


    cdef int __query_ball_point_traverse_checking(cKDTree self,
                                                   list results,
                                                   innernode* node,
                                                   np.float64_t* x,
                                                   np.float64_t r,
                                                   np.float64_t p,
                                                   np.float64_t epsfac,
                                                   np.float64_t invepsfac,
                                                   Rectangle rect,
                                                   np.float64_t min_distance,
                                                   np.float64_t max_distance) except -1:
        cdef leafnode* lnode
        cdef innernode* inode
        
        cdef np.float64_t save_min, save_max
        cdef np.float64_t part_min_distance = 0., part_max_distance = 0.
        cdef np.float64_t d
        cdef np.npy_intp k, i, j

        if min_distance > r*epsfac:
            return 0
        elif max_distance < r*invepsfac:
            self.__query_ball_point_traverse_no_checking(results, node)
        elif node.split_dim == -1:  # leaf node
            lnode = <leafnode*>node
            # brute-force
            for i in range(lnode.start_idx, lnode.end_idx):
                d = _distance_p(
                    self.raw_data + self.raw_indices[i] * self.m,
                    x, p, self.m, r)
                if d <= r:
                    results.append(self.raw_indices[i])
        else:
            inode = <innernode*>node
                
            k = inode.split_dim
            if p != infinity:
                part_min_distance = min_distance - min_dist_point_interval_p(x, rect, k, p)
                part_max_distance = max_distance - max_dist_point_interval_p(x, rect, k, p)
            
            # Go to box with lesser component along k
            # less.maxes[k] goes from rect.maxes[k] to inode.split
            save_max = rect.maxes[k]
            rect.maxes[k] = inode.split
            if p != infinity:
                min_distance = part_min_distance + min_dist_point_interval_p(x, rect, k, p)
                max_distance = part_max_distance + max_dist_point_interval_p(x, rect, k, p)
            else:
                min_distance = min_dist_point_rect_p_inf(x, rect)
                max_distance = max_dist_point_rect_p_inf(x, rect)
                    
            self.__query_ball_point_traverse_checking(results, inode.less,
                                                      x, r, p, epsfac, invepsfac,
                                                      rect,
                                                      min_distance, max_distance)
            rect.maxes[k] = save_max
            
            # Go to box with greater component along k
            # greater.mins[k] goes from rect.mins[k] to inode.split
            save_min = rect.mins[k]
            rect.mins[k] = inode.split
            if p != infinity:
                min_distance = part_min_distance + min_dist_point_interval_p(x, rect, k, p)
                max_distance = part_max_distance + max_dist_point_interval_p(x, rect, k, p)
            else:
                min_distance = min_dist_point_rect_p_inf(x, rect)
                max_distance = max_dist_point_rect_p_inf(x, rect)
                    
            self.__query_ball_point_traverse_checking(results, inode.greater,
                                                      x, r, p, epsfac, invepsfac,
                                                      rect,
                                                      min_distance, max_distance)
            rect.mins[k] = save_min
            return 0


    cdef list __query_ball_point(cKDTree self,
                                 np.float64_t* x,
                                 np.float64_t r,
                                 np.float64_t p,
                                 np.float64_t eps):

        cdef list results
        cdef Rectangle rect
        cdef np.float64_t epsfac, invepsfac
        cdef np.float64_t min_distance, max_distance
        cdef np.npy_intp i


        # internally we represent all distances as distance**p
        if p != infinity and r != infinity:
            r = r ** p

        # fiddle approximation factor
        if eps == 0:
            epsfac = 1
        elif p == infinity:
            epsfac = 1/(1+eps)
        else:
            epsfac = 1/(1+eps)**p
        invepsfac = 1/epsfac

        # Calculate mins and maxes to outer box
        rect.m = self.m
        rect.mins = rect.maxes = <np.float64_t*> NULL
        try:
            
            rect.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect.mins == <np.float64_t*> NULL:
                raise MemoryError

            rect.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect.maxes == <np.float64_t*> NULL:
                raise MemoryError

            for i in range(self.m):
                rect.mins[i] = self.raw_mins[i]
                rect.maxes[i] = self.raw_maxes[i]

            # Computer first min and max distances
            if p == infinity:
                min_distance = min_dist_point_rect_p_inf(x, rect)
                max_distance = max_dist_point_rect_p_inf(x, rect)
            else:
                min_distance = 0.
                max_distance = 0.
                for i in range(self.m):
                    min_distance += min_dist_point_interval_p(x, rect, i, p)
                    max_distance += max_dist_point_interval_p(x, rect, i, p)
                    
            results = []
            self.__query_ball_point_traverse_checking(results, self.tree,
                                                    x, r, p, epsfac, invepsfac,
                                                    rect,
                                                    min_distance, max_distance)

        finally:
            if rect.mins != <np.float64_t*> NULL:
                stdlib.free(rect.mins)

            if rect.maxes != <np.float64_t*> NULL:
                stdlib.free(rect.maxes)

        return results


    def query_ball_point(cKDTree self, object x, np.float64_t r,
                         np.float64_t p=2., np.float64_t eps=0):
        """query_ball_point(self, x, r, p, eps)
        
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
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] xx
        
        x = np.asarray(x).astype(np.float64)
        if x.shape[-1] != self.m:
            raise ValueError("Searching for a %d-dimensional point in a " \
                             "%d-dimensional KDTree" % (int(x.shape[-1]), int(self.m)))
        if len(x.shape) == 1:
            xx = np.ascontiguousarray(x, dtype=np.float64)
            return self.__query_ball_point(<np.float64_t*>np.PyArray_DATA(xx), r, p, eps)
        else:
            retshape = x.shape[:-1]
            result = np.empty(retshape, dtype=np.object)
            for c in np.ndindex(retshape):
                xx = np.ascontiguousarray(x[c], dtype=np.float64)
                result[c] = self.__query_ball_point(
                    <np.float64_t*>np.PyArray_DATA(xx), r, p, eps)
            return result

    # ---------------
    # query_ball_tree
    # ---------------
    cdef int __query_ball_tree_traverse_no_checking(cKDTree self,
                                                     cKDTree other,
                                                     list results,
                                                     innernode* node1,
                                                     innernode* node2) except -1:
        cdef leafnode *lnode1, *lnode2
        cdef list results_i
        cdef np.npy_intp i, j
        
        if node1.split_dim == -1:  # leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # leaf node
                lnode2 = <leafnode*>node2
                
                for i in range(lnode1.start_idx, lnode1.end_idx):
                    results_i = results[self.raw_indices[i]]
                    for j in range(lnode2.start_idx, lnode2.end_idx):
                        results_i.append(other.raw_indices[j])
            else:
                
                self.__query_ball_tree_traverse_no_checking(other, results, node1, node2.less)
                self.__query_ball_tree_traverse_no_checking(other, results, node1, node2.greater)
        else:
            
            self.__query_ball_tree_traverse_no_checking(other, results, node1.less, node2)
            self.__query_ball_tree_traverse_no_checking(other, results, node1.greater, node2)

        return 0


    cdef int __query_ball_tree_traverse_checking(cKDTree self,
                                                 cKDTree other,
                                                 list results,
                                                 innernode* node1,
                                                 innernode* node2,
                                                 np.float64_t r,
                                                 np.float64_t p,
                                                 np.float64_t epsfac,
                                                 np.float64_t invepsfac,
                                                 Rectangle rect1,
                                                 Rectangle rect2,
                                                 np.float64_t min_distance,
                                                 np.float64_t max_distance) except -1:
        cdef leafnode *lnode1, *lnode2
        cdef innernode *inode1, *inode2        
        cdef list results_i
        cdef np.float64_t save_min1, save_max1
        cdef np.float64_t save_min2, save_max2
        cdef np.float64_t part_min_distance1 = 0., part_max_distance1 = 0.
        cdef np.float64_t part_min_distance2 = 0., part_max_distance2 = 0.
        cdef np.float64_t d
        cdef np.npy_intp k1, k2, i, j

        
        if min_distance > r*epsfac:
            return 0
        elif max_distance < r*invepsfac:
            self.__query_ball_tree_traverse_no_checking(other, results, node1, node2)
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = <leafnode*>node2
                
                # brute-force
                for i in range(lnode1.start_idx, lnode1.end_idx):
                    results_i = results[self.raw_indices[i]]
                    for j in range(lnode2.start_idx, lnode2.end_idx):
                        d = _distance_p(
                            self.raw_data + self.raw_indices[i] * self.m,
                            other.raw_data + other.raw_indices[j] * other.m,
                            p, self.m, r)
                        if d <= r:
                            results_i.append(other.raw_indices[j])
                            
            else:  # 1 is a leaf node, 2 is inner node
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1, node2.less,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1, node2.greater,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                rect2.mins[k2] = save_min2
            
                
        else:  # 1 is an inner node
            k1 = node1.split_dim
            __rect_preupdate(rect1, rect2, k1, p, min_distance,
                             max_distance, &part_min_distance1,
                             &part_max_distance1)
                
            # node1 goes to box with lesser component along k1
            # node1.less.maxes[k1] changes from rect1.maxes[k1] to node1.split
            save_max1 = rect1.maxes[k1]
            rect1.maxes[k1] = node1.split
            __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                              &max_distance, part_min_distance1,
                              part_max_distance1)

            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1.less, node2,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
            else: # 1 and 2 are inner nodes
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1.less, node2.less,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1.less, node2.greater,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
            rect1.maxes[k1] = save_max1
                    
            # node1 goes to box with greater component along k1
            # node1.greater.mins[k1] changes from rect1.mins[k1] to node1.split
            save_min1 = rect1.mins[k1]
            rect1.mins[k1] = node1.split
            __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                              &max_distance, part_min_distance1,
                              part_max_distance1)
                
            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1.greater, node2,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
            else: # 1 and 2 are inner nodes
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1.greater, node2.less,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from rect2.mins[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_ball_tree_traverse_checking(other, results,
                                                         node1.greater, node2.greater,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
            rect1.mins[k1] = save_min1
            return 0
            

    def query_ball_tree(cKDTree self, cKDTree other,
                        np.float64_t r, np.float64_t p=2., np.float64_t eps=0):
        """query_ball_tree(self, other, r, p, eps)

        Find all pairs of points whose distance is at most r

        Parameters
        ----------
        other : KDTree instance
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
        cdef list results
        cdef Rectangle rect1, rect2
        cdef np.float64_t epsfac, invepsfac
        cdef np.float64_t min_distance, max_distance
        cdef np.npy_intp i

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to query_ball_trees have different dimensionality")

        # internally we represent all distances as distance**p
        if p != infinity and r != infinity:
            r = r ** p

        # fiddle approximation factor
        if eps == 0:
            epsfac = 1
        elif p == infinity:
            epsfac = 1/(1+eps)
        else:
            epsfac = 1/(1+eps)**p
        invepsfac = 1/epsfac

        # Calculate mins and maxes to outer box
        rect1.m = rect2.m = self.m
        rect1.mins = rect1.maxes = rect2.mins = rect2.maxes = <np.float64_t*> NULL
        try:
            rect1.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect1.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.maxes == <np.float64_t*> NULL: 
                raise MemoryError    

            rect2.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect2.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.maxes == <np.float64_t*> NULL: 
                raise MemoryError

            for i in range(self.m):
                rect1.mins[i] = self.raw_mins[i]
                rect1.maxes[i] = self.raw_maxes[i]
                rect2.mins[i] = other.raw_mins[i]
                rect2.maxes[i] = other.raw_maxes[i]

            # Compute first min and max distances
            if p == infinity:
                min_distance = min_dist_rect_rect_p_inf(rect1, rect2)
                max_distance = max_dist_rect_rect_p_inf(rect1, rect2)
            else:
                min_distance = 0.
                max_distance = 0.
                for i in range(self.m):
                    min_distance += min_dist_interval_interval_p(rect1, rect2, i, p)
                    max_distance += max_dist_interval_interval_p(rect1, rect2, i, p)
                    
            results = [[] for i in range(self.n)]
            self.__query_ball_tree_traverse_checking(other, results,
                                                    self.tree, other.tree,
                                                    r, p, epsfac, invepsfac,
                                                    rect1, rect2,
                                                    min_distance, max_distance)

        finally:
            if rect1.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect1.mins)

            if rect1.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect1.maxes)

            if rect2.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect2.mins)

            if rect2.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect2.maxes)

        return results


    # -----------
    # query_pairs
    # -----------
    cdef int __query_pairs_traverse_no_checking(cKDTree self,
                                                 set results,
                                                 innernode* node1,
                                                 innernode* node2) except -1:
        cdef leafnode *lnode1, *lnode2
        cdef list results_i
        cdef np.npy_intp i, j
        
        if node1.split_dim == -1:  # leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # leaf node
                lnode2 = <leafnode*>node2

                # Special care here to avoid duplicate pairs
                if node1 == node2:
                    for i in range(lnode1.start_idx, lnode1.end_idx):
                        for j in range(i+1, lnode2.end_idx):
                            if self.raw_indices[i] < self.raw_indices[j]:
                                results.add((self.raw_indices[i], self.raw_indices[j]))
                            else:
                                results.add((self.raw_indices[j], self.raw_indices[i]))
                else:
                    for i in range(lnode1.start_idx, lnode1.end_idx):
                        for j in range(lnode2.start_idx, lnode2.end_idx):
                            if self.raw_indices[i] < self.raw_indices[j]:
                                results.add((self.raw_indices[i], self.raw_indices[j]))
                            else:
                                results.add((self.raw_indices[j], self.raw_indices[i]))
                            
            else:
                self.__query_pairs_traverse_no_checking(results, node1, node2.less)
                self.__query_pairs_traverse_no_checking(results, node1, node2.greater)
        else:
            if node1 == node2:
                # Avoid traversing (node1.less, node2.greater) and
                # (node1.greater, node2.less) (it's the same node pair twice
                # over, which is the source of the complication in the
                # original KDTree.query_pairs)
                self.__query_pairs_traverse_no_checking(results, node1.less, node2.less)
                self.__query_pairs_traverse_no_checking(results, node1.less, node2.greater)
                self.__query_pairs_traverse_no_checking(results, node1.greater, node2.greater)
            else:
                self.__query_pairs_traverse_no_checking(results, node1.less, node2)
                self.__query_pairs_traverse_no_checking(results, node1.greater, node2)

        return 0


    cdef int __query_pairs_traverse_checking(cKDTree self,
                                              set results,
                                              innernode* node1,
                                              innernode* node2,
                                              np.float64_t r,
                                              np.float64_t p,
                                              np.float64_t epsfac,
                                              np.float64_t invepsfac,
                                              Rectangle rect1,
                                              Rectangle rect2,
                                              np.float64_t min_distance,
                                              np.float64_t max_distance) except -1:
        cdef leafnode *lnode1, *lnode2
        cdef innernode *inode1, *inode2
        cdef list results_i
        cdef np.float64_t save_min1, save_max1
        cdef np.float64_t save_min2, save_max2
        cdef np.float64_t part_min_distance1 = 0., part_max_distance1 = 0.
        cdef np.float64_t part_min_distance2 = 0., part_max_distance2 = 0.
        cdef np.float64_t d
        cdef np.npy_intp k1, k2, i, j

        if min_distance > r*epsfac:
            return 0
        elif max_distance < r*invepsfac:
            self.__query_pairs_traverse_no_checking(results, node1, node2)
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = <leafnode*>node2
                
                # brute-force
                # Special care here to avoid duplicate pairs
                if node1 == node2:
                    for i in range(lnode1.start_idx, lnode1.end_idx):
                        for j in range(i+1, lnode2.end_idx):
                            d = _distance_p(
                                self.raw_data + self.raw_indices[i] * self.m,
                                self.raw_data + self.raw_indices[j] * self.m,
                                p, self.m, r)
                            if d <= r:
                                if self.raw_indices[i] < self.raw_indices[j]:
                                    results.add((self.raw_indices[i], self.raw_indices[j]))
                                else:
                                    results.add((self.raw_indices[j], self.raw_indices[i]))
                else:
                    for i in range(lnode1.start_idx, lnode1.end_idx):
                        for j in range(lnode2.start_idx, lnode2.end_idx):
                            d = _distance_p(
                                self.raw_data + self.raw_indices[i] * self.m,
                                self.raw_data + self.raw_indices[j] * self.m,
                                p, self.m, r)
                            if d <= r:
                                if self.raw_indices[i] < self.raw_indices[j]:
                                    results.add((self.raw_indices[i], self.raw_indices[j]))
                                else:
                                    results.add((self.raw_indices[j], self.raw_indices[i]))
                            
            else:  # 1 is a leaf node, 2 is inner node
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_pairs_traverse_checking(results,
                                                     node1, node2.less,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_pairs_traverse_checking(results,
                                                     node1, node2.greater,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
        else:  # 1 is an inner node
            k1 = node1.split_dim
            __rect_preupdate(rect1, rect2, k1, p, min_distance,
                             max_distance, &part_min_distance1,
                             &part_max_distance1)
                
            # node1 goes to box with lesser component along k1
            # node1.less.maxes[k1] changes from rect1.maxes[k1] to node1.split
            save_max1 = rect1.maxes[k1]
            rect1.maxes[k1] = node1.split
            __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                              &max_distance, part_min_distance1,
                              part_max_distance1)

            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                self.__query_pairs_traverse_checking(results,
                                                     node1.less, node2,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
            else: # 1 and 2 are inner nodes
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_pairs_traverse_checking(results,
                                                     node1.less, node2.less,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_pairs_traverse_checking(results,
                                                     node1.less, node2.greater,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
            rect1.maxes[k1] = save_max1
                    
            # node1 goes to box with greater component along k1
            # node1.greater.mins[k1] changes from rect1.mins[k1] to node1.split
            save_min1 = rect1.mins[k1]
            rect1.mins[k1] = node1.split
            __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                              &max_distance, part_min_distance1,
                              part_max_distance1)
                
            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                self.__query_pairs_traverse_checking(results,
                                                     node1.greater, node2,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
            else: # 1 and 2 are inner nodes
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)

                if node1 != node2:
                    # Avoid traversing (node1.less, node2.greater) and
                    # (node1.greater, node2.less) (it's the same node pair twice
                    # over, which is the source of the complication in the
                    # original KDTree.query_pairs)
                    
                    # node2 goes to box with lesser component along k2
                    # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                    save_max2 = rect2.maxes[k2]
                    rect2.maxes[k2] = node2.split
                    __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                      &max_distance, part_min_distance2,
                                      part_max_distance2)
                    self.__query_pairs_traverse_checking(results,
                                                         node1.greater, node2.less,
                                                         r, p, epsfac, invepsfac,
                                                         rect1, rect2,
                                                         min_distance, max_distance)
                    rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from rect2.mins[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__query_pairs_traverse_checking(results,
                                                     node1.greater, node2.greater,
                                                     r, p, epsfac, invepsfac,
                                                     rect1, rect2,
                                                     min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
            rect1.mins[k1] = save_min1
            return 0
            

    def query_pairs(cKDTree self, np.float64_t r, np.float64_t p=2.,
                    np.float64_t eps=0):
        """query_pairs(self, r, p, eps)

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

        Returns
        -------
        results : set
            Set of pairs ``(i,j)``, with ``i < j`, for which the corresponding
            positions are close.

        """
        
        cdef set results, _results
        cdef Rectangle rect1, rect2
        cdef np.float64_t epsfac, invepsfac
        cdef np.float64_t min_distance, max_distance
        cdef np.npy_intp i,j

        # internally we represent all distances as distance**p
        if p != infinity and r != infinity:
            r = r ** p

        # fiddle approximation factor
        if eps == 0:
            epsfac = 1
        elif p == infinity:
            epsfac = 1/(1+eps)
        else:
            epsfac = 1/(1+eps)**p
        invepsfac = 1/epsfac

        # Calculate mins and maxes to outer box
        rect1.m = rect2.m = self.m
        rect1.mins = rect1.maxes = rect2.mins = rect2.maxes = <np.float64_t*> NULL
        try:
            
            rect1.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect1.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.maxes == <np.float64_t*> NULL: 
                raise MemoryError    

            rect2.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect2.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.maxes == <np.float64_t*> NULL: 
                raise MemoryError

            for i in range(self.m):
                rect1.mins[i] = self.raw_mins[i]
                rect1.maxes[i] = self.raw_maxes[i]
                rect2.mins[i] = rect1.mins[i]
                rect2.maxes[i] = rect2.maxes[i]

            # Compute first min and max distances
            if p == infinity:
                min_distance = min_dist_rect_rect_p_inf(rect1, rect2)
                max_distance = max_dist_rect_rect_p_inf(rect1, rect2)
            else:
                min_distance = 0.
                max_distance = 0.
                for i in range(self.m):
                    min_distance += min_dist_interval_interval_p(rect1, rect2, i, p)
                    max_distance += max_dist_interval_interval_p(rect1, rect2, i, p)
                    
            _results = set()
            self.__query_pairs_traverse_checking(_results,
                                                self.tree, self.tree,
                                                r, p, epsfac, invepsfac,
                                                rect1, rect2,
                                                min_distance, max_distance)

        finally:
            if rect1.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect1.mins)

            if rect1.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect1.maxes)

            if rect2.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect2.mins)

            if rect2.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect2.maxes)

        ## On platforms where a C long is shorter than npy_intp,
        ## e.g. Windows 64, we convert to Python int if they don't
        ## overflow.    
        
        if sizeof(long) < sizeof(np.npy_intp):
            results = set()
            for i,j in _results:    
                if (i <= <np.npy_intp> LONG_MAX) and (j <= <np.npy_intp> LONG_MAX):
                    results.add((int(i),int(j)))
                else:
                    results.add((i,j))
                
        ## other platforms 
        else:
            results = _results

        return results


    # ---------------
    # count_neighbors
    # ---------------
    cdef int __count_neighbors_traverse(cKDTree self,
                                        cKDTree other,
                                        np.npy_intp n_queries,
                                        np.float64_t* r,
                                        np.npy_intp * results,
                                        np.npy_intp * idx,
                                        innernode* node1,
                                        innernode* node2,
                                        np.float64_t p,
                                        Rectangle rect1,
                                        Rectangle rect2,
                                        np.float64_t min_distance,
                                        np.float64_t max_distance) except -1:
        cdef leafnode *lnode1, *lnode2
        cdef innernode *inode1, *inode2
        
        cdef np.float64_t save_min1, save_max1
        cdef np.float64_t save_min2, save_max2
        cdef np.float64_t part_min_distance1 = 0., part_max_distance1 = 0.
        cdef np.float64_t part_min_distance2 = 0., part_max_distance2 = 0.
        cdef np.float64_t d
        cdef np.npy_intp  *old_idx
        cdef np.npy_intp old_n_queries, k1, k2, l, i, j

        # Speed through pairs of nodes all of whose children are close
        # and see if any work remains to be done
        old_idx = idx

        try:
            idx = <np.npy_intp *> stdlib.malloc(n_queries * sizeof(np.npy_intp ))
            if idx == <np.npy_intp *> NULL:
                raise MemoryError

            old_n_queries = n_queries
            n_queries = 0
            for i in range(old_n_queries):
                if max_distance < r[old_idx[i]]:
                    results[old_idx[i]] += node1.children * node2.children
                elif min_distance <= r[old_idx[i]]:
                    idx[n_queries] = old_idx[i]
                    n_queries += 1

            if n_queries >= 0:
                # OK, need to probe a bit deeper
                if node1.split_dim == -1:  # 1 is leaf node
                    lnode1 = <leafnode*>node1
                    if node2.split_dim == -1:  # 1 & 2 are leaves
                        lnode2 = <leafnode*>node2
                        
                        # brute-force
                        for i in range(lnode1.start_idx, lnode1.end_idx):
                            for j in range(lnode2.start_idx, lnode2.end_idx):
                                d = _distance_p(
                                    self.raw_data + self.raw_indices[i] * self.m,
                                    other.raw_data + other.raw_indices[j] * other.m,
                                    p, self.m, max_distance)

                                # I think it's usually cheaper to test d against all r's
                                # than to generate a distance array, sort it, then
                                # search for all r's via binary search
                                for l in range(n_queries):
                                    if d <= r[idx[l]]:
                                        results[idx[l]] += 1
                                    
                    else:  # 1 is a leaf node, 2 is inner node
                        k2 = node2.split_dim
                        __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                         max_distance, &part_min_distance2,
                                         &part_max_distance2)
                            
                        # node2 goes to box with lesser component along k2
                        # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                        save_max2 = rect2.maxes[k2]
                        rect2.maxes[k2] = node2.split
                        __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                          &max_distance, part_min_distance2,
                                          part_max_distance2)
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1, node2.less,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                        rect2.maxes[k2] = save_max2
                            
                        # node2 goes to box with greater component along k2
                        # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                        save_min2 = rect2.mins[k2]
                        rect2.mins[k2] = node2.split
                        __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                          &max_distance, part_min_distance2,
                                          part_max_distance2)
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1, node2.greater,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                        rect2.mins[k2] = save_min2
                    
                else:  # 1 is an inner node
                    k1 = node1.split_dim
                    __rect_preupdate(rect1, rect2, k1, p, min_distance,
                                     max_distance, &part_min_distance1,
                                     &part_max_distance1)
                        
                    # node1 goes to box with lesser component along k1
                    # node1.less.maxes[k1] changes from rect1.maxes[k1] to node1.split
                    save_max1 = rect1.maxes[k1]
                    rect1.maxes[k1] = node1.split
                    __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                                      &max_distance, part_min_distance1,
                                      part_max_distance1)
        
                    if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1.less, node2,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                    else: # 1 and 2 are inner nodes
                        k2 = node2.split_dim
                        __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                         max_distance, &part_min_distance2,
                                         &part_max_distance2)
                            
                        # node2 goes to box with lesser component along k2
                        # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                        save_max2 = rect2.maxes[k2]
                        rect2.maxes[k2] = node2.split
                        __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                          &max_distance, part_min_distance2,
                                          part_max_distance2)
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1.less, node2.less,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                        rect2.maxes[k2] = save_max2
                            
                        # node2 goes to box with greater component along k2
                        # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                        save_min2 = rect2.mins[k2]
                        rect2.mins[k2] = node2.split
                        __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                          &max_distance, part_min_distance2,
                                          part_max_distance2)
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1.less, node2.greater,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                        rect2.mins[k2] = save_min2
                        
                    rect1.maxes[k1] = save_max1
                            
                    # node1 goes to box with greater component along k1
                    # node1.greater.mins[k1] changes from rect1.mins[k1] to node1.split
                    save_min1 = rect1.mins[k1]
                    rect1.mins[k1] = node1.split
                    __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                                      &max_distance, part_min_distance1,
                                      part_max_distance1)
                        
                    if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1.greater, node2,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                    else: # 1 and 2 are inner nodes
                        k2 = node2.split_dim
                        __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                         max_distance, &part_min_distance2,
                                         &part_max_distance2)
                            
                        # node2 goes to box with lesser component along k2
                        # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                        save_max2 = rect2.maxes[k2]
                        rect2.maxes[k2] = node2.split
                        __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                          &max_distance, part_min_distance2,
                                          part_max_distance2)
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1.greater, node2.less,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                        rect2.maxes[k2] = save_max2
                            
                        # node2 goes to box with greater component along k2
                        # node2.greater.mins[k2] changes from rect2.mins[k2] to node2.split
                        save_min2 = rect2.mins[k2]
                        rect2.mins[k2] = node2.split
                        __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                          &max_distance, part_min_distance2,
                                          part_max_distance2)
                        self.__count_neighbors_traverse(other, n_queries, r, results, idx,
                                                        node1.greater, node2.greater,
                                                        p, rect1, rect2,
                                                        min_distance, max_distance)
                        rect2.mins[k2] = save_min2
                        
                    rect1.mins[k1] = save_min1
        finally:
            # Free memory
            if idx != <np.npy_intp *> NULL:
                stdlib.free(idx)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def count_neighbors(cKDTree self, cKDTree other, object r, np.float64_t p=2.):
        """count_neighbors(self, other, r, p)

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
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use

        Returns
        -------
        result : int or 1-D array of ints
            The number of pairs. Note that this is internally stored in a numpy int,
            and so may overflow if very large (2e9).

        """
        cdef np.npy_intp i, n_queries
        cdef Rectangle rect1, rect2
        cdef np.float64_t epsfac, invepsfac
        cdef np.float64_t min_distance, max_distance
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] real_r
        cdef np.ndarray[np.npy_intp, ndim=1, mode="c"] results, idx

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to query_ball_trees have different dimensionality")

        # Make a copy of r array to ensure it's contiguous and to modify it
        # below
        if np.shape(r) == ():
            real_r = np.array([r], dtype=np.float64)
            n_queries = 1
        elif len(np.shape(r))==1:
            real_r = np.array(r, dtype=np.float64)
            n_queries = r.shape[0]
        else:
            raise ValueError("r must be either a single value or a one-dimensional array of values")

                
        # internally we represent all distances as distance**p
        if p != infinity:
            for i in range(n_queries):
                if real_r[i] != infinity:
                    real_r[i] = real_r[i] ** p

        # Calculate mins and maxes to outer box
        rect1.m = rect2.m = self.m
        rect1.mins = rect1.maxes = rect2.mins = rect2.maxes = <np.float64_t*> NULL
        try:
            
            rect1.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect1.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.maxes == <np.float64_t*> NULL: 
                raise MemoryError    

            rect2.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect2.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.maxes == <np.float64_t*> NULL: 
                raise MemoryError

            for i in range(self.m):
                rect1.mins[i] = self.raw_mins[i]
                rect1.maxes[i] = self.raw_maxes[i]
                rect2.mins[i] = other.raw_mins[i]
                rect2.maxes[i] = other.raw_maxes[i]

            # Compute first min and max distances
            if p == infinity:
                min_distance = min_dist_rect_rect_p_inf(rect1, rect2)
                max_distance = max_dist_rect_rect_p_inf(rect1, rect2)
            else:
                min_distance = 0.
                max_distance = 0.
                for i in range(self.m):
                    min_distance += min_dist_interval_interval_p(rect1, rect2, i, p)
                    max_distance += max_dist_interval_interval_p(rect1, rect2, i, p)
                    
            # Go!
            results = np.zeros((n_queries,), dtype=npy_intp_dtype)
            idx = np.arange(n_queries, dtype=npy_intp_dtype)
            self.__count_neighbors_traverse(other, n_queries,
                                            <np.float64_t*>np.PyArray_DATA(real_r),
                                            <np.npy_intp *>np.PyArray_DATA(results),
                                            <np.npy_intp *>np.PyArray_DATA(idx),
                                            self.tree, other.tree,
                                            p, rect1, rect2,
                                            min_distance, max_distance)

        finally:
            if rect1.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect1.mins)

            if rect1.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect1.maxes)

            if rect2.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect2.mins)

            if rect2.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect2.maxes)

        if np.shape(r) == ():
            if results[0] <= <np.npy_intp> LONG_MAX:
                return int(results[0])
            else:
                return results[0]
        elif len(np.shape(r))==1:
            return results

    # ----------------------
    # sparse_distance_matrix
    # ----------------------
    cdef int __sparse_distance_matrix_traverse(cKDTree self, cKDTree other,
                                                coo_entries results,
                                                innernode* node1, innernode* node2,
                                                np.float64_t r,
                                                np.float64_t p,
                                                Rectangle rect1,
                                                Rectangle rect2,
                                                np.float64_t min_distance,
                                                np.float64_t max_distance) except -1:
        cdef leafnode *lnode1, *lnode2
        cdef innernode *inode1, *inode2
        cdef list results_i

        cdef np.float64_t save_min1, save_max1
        cdef np.float64_t save_min2, save_max2
        cdef np.float64_t part_min_distance1 = 0., part_max_distance1 = 0.
        cdef np.float64_t part_min_distance2 = 0., part_max_distance2 = 0.
        cdef np.float64_t d
        cdef np.npy_intp k1, k2, i, j
                
        if min_distance > r:
            return 0
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = <leafnode*>node2
                
                # brute-force
                # Special care here to avoid duplicate pairs
                if node1 == node2:
                    # self == other if we get here
                    for i in range(lnode1.start_idx, lnode1.end_idx):
                        for j in range(i+1, lnode2.end_idx):
                            d = _distance_p(
                                self.raw_data + self.raw_indices[i] * self.m,
                                self.raw_data + self.raw_indices[j] * self.m,
                                p, self.m, r)
                            if d <= r:
                                results.add(self.raw_indices[i],
                                            self.raw_indices[j], d)
                                results.add(self.raw_indices[j],
                                            self.raw_indices[i], d)
                else:
                    for i in range(lnode1.start_idx, lnode1.end_idx):
                        for j in range(lnode2.start_idx, lnode2.end_idx):
                            d = _distance_p(
                                self.raw_data + self.raw_indices[i] * self.m,
                                other.raw_data + other.raw_indices[j] * self.m,
                                p, self.m, r)
                            if d <= r:
                                results.add(self.raw_indices[i],
                                            other.raw_indices[j], d)
                            
            else:  # 1 is a leaf node, 2 is inner node
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1, node2.less,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                  &max_distance, part_min_distance2,
                                  part_max_distance2)
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1, node2.greater,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
        else:  # 1 is an inner node
            k1 = node1.split_dim
            __rect_preupdate(rect1, rect2, k1, p, min_distance,
                             max_distance, &part_min_distance1,
                             &part_max_distance1)
                
            # node1 goes to box with lesser component along k1
            # node1.less.maxes[k1] changes from rect1.maxes[k1] to node1.split
            save_max1 = rect1.maxes[k1]
            rect1.maxes[k1] = node1.split
            __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                             &max_distance, part_min_distance1,
                             part_max_distance1)

            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1.less, node2,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
            else: # 1 and 2 are inner nodes
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)
                    
                # node2 goes to box with lesser component along k2
                # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                save_max2 = rect2.maxes[k2]
                rect2.maxes[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                 &max_distance, part_min_distance2,
                                 part_max_distance2)
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1.less, node2.less,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
                rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from mins2[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                 &max_distance, part_min_distance2,
                                 part_max_distance2)
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1.less, node2.greater,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
            rect1.maxes[k1] = save_max1
                    
            # node1 goes to box with greater component along k1
            # node1.greater.mins[k1] changes from rect1.mins[k1] to node1.split
            save_min1 = rect1.mins[k1]
            rect1.mins[k1] = node1.split
            __rect_postupdate(rect1, rect2, k1, p, &min_distance,
                             &max_distance, part_min_distance1,
                             part_max_distance1)
                
            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1.greater, node2,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
            else: # 1 and 2 are inner nodes
                k2 = node2.split_dim
                __rect_preupdate(rect1, rect2, k2, p, min_distance,
                                 max_distance, &part_min_distance2,
                                 &part_max_distance2)

                if node1 != node2:
                    # Avoid traversing (node1.less, node2.greater) and
                    # (node1.greater, node2.less) (it's the same node pair twice
                    # over, which is the source of the complication in the
                    # original KDTree.sparse_distance_matrix)
                    
                    # node2 goes to box with lesser component along k2
                    # node2.less.maxes[k2] changes from rect2.maxes[k2] to node2.split
                    save_max2 = rect2.maxes[k2]
                    rect2.maxes[k2] = node2.split
                    __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                     &max_distance, part_min_distance2,
                                     part_max_distance2)
                    self.__sparse_distance_matrix_traverse(other, results,
                                                           node1.greater, node2.less,
                                                           r, p,
                                                           rect1, rect2,
                                                           min_distance, max_distance)
                    rect2.maxes[k2] = save_max2
                    
                # node2 goes to box with greater component along k2
                # node2.greater.mins[k2] changes from rect2.mins[k2] to node2.split
                save_min2 = rect2.mins[k2]
                rect2.mins[k2] = node2.split
                __rect_postupdate(rect1, rect2, k2, p, &min_distance,
                                 &max_distance, part_min_distance2,
                                 part_max_distance2)
                self.__sparse_distance_matrix_traverse(other, results,
                                                       node1.greater, node2.greater,
                                                       r, p,
                                                       rect1, rect2,
                                                       min_distance, max_distance)
                rect2.mins[k2] = save_min2
                
            rect1.mins[k1] = save_min1
            return 0
            
    def sparse_distance_matrix(cKDTree self, cKDTree other, np.float64_t r,
                               np.float64_t p=2.):
        """sparse_distance_matrix(self, r, p)

        Compute a sparse distance matrix

        Computes a distance matrix between two KDTrees, leaving as zero
        any distance greater than r.

        Parameters
        ----------
        other : cKDTree

        r : positive float
            FIXME: KDTree calls this parameter max_distance

        Returns
        -------
        result : dok_matrix
            Sparse matrix representing the results in "dictionary of keys" format.
            FIXME: Internally, built as a COO matrix, it would be more
            efficient to return this COO matrix.

        """
        cdef np.npy_intp i
        cdef coo_entries results
        cdef Rectangle rect1, rect2
        cdef np.float64_t min_distance, max_distance

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to query_ball_trees have different dimensionality")

        # internally we represent all distances as distance**p
        if p != infinity and r != infinity:
            r = r ** p

        # Calculate mins and maxes to outer box
        rect1.m = rect2.m = self.m

        rect1.mins = rect1.maxes = rect2.mins = rect2.maxes = <np.float64_t*> NULL
        try:
            rect1.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect1.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect1.maxes == <np.float64_t*> NULL: 
                raise MemoryError    

            rect2.mins = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.mins == <np.float64_t*> NULL: 
                raise MemoryError

            rect2.maxes = <np.float64_t*>stdlib.malloc(self.m * sizeof(np.float64_t))
            if rect2.maxes == <np.float64_t*> NULL: 
                raise MemoryError

            for i in range(self.m):
                rect1.mins[i] = self.raw_mins[i]
                rect1.maxes[i] = self.raw_maxes[i]
                rect2.mins[i] = other.raw_mins[i]
                rect2.maxes[i] = other.raw_maxes[i]

            # Compute first min and max distances
            if p == infinity:
                min_distance = min_dist_rect_rect_p_inf(rect1, rect2)
                max_distance = max_dist_rect_rect_p_inf(rect1, rect2)
            else:
                min_distance = 0.
                max_distance = 0.
                for i in range(self.m):
                    min_distance += min_dist_interval_interval_p(rect1, rect2, i, p)
                    max_distance += max_dist_interval_interval_p(rect1, rect2, i, p)
                    
            results = coo_entries()
            self.__sparse_distance_matrix_traverse(other, results,
                                                self.tree, other.tree,
                                                r, p, rect1, rect2,
                                                min_distance, max_distance)
        finally:
            if rect1.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect1.mins)

            if rect1.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect1.maxes)

            if rect2.mins  != <np.float64_t*> NULL: 
                stdlib.free(rect2.mins)

            if rect2.maxes != <np.float64_t*> NULL: 
                stdlib.free(rect2.maxes)

        return results.to_matrix(shape=(self.n, other.n)).todok()

