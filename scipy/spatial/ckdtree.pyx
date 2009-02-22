# Copyright Anne M. Archibald 2008
# Released under the scipy license
import numpy as np
cimport numpy as np
cimport stdlib

import kdtree

cdef double infinity = np.inf




# priority queue
cdef union heapcontents:
    int intdata
    char* ptrdata

cdef struct heapitem:
    double priority
    heapcontents contents

cdef struct heap:
    int n
    heapitem* heap
    int space

cdef inline heapcreate(heap* self,int initial_size):
    self.space = initial_size
    self.heap = <heapitem*>stdlib.malloc(sizeof(heapitem)*self.space)
    self.n=0

cdef inline heapdestroy(heap* self):
    stdlib.free(self.heap)

cdef inline heapresize(heap* self, int new_space):
    if new_space<self.n:
        raise ValueError("Heap containing %d items cannot be resized to %d" % (self.n, new_space))
    self.space = new_space
    self.heap = <heapitem*>stdlib.realloc(<void*>self.heap,new_space*sizeof(heapitem))

cdef inline heappush(heap* self, heapitem item):
    cdef int i
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

cdef heapitem heappeek(heap* self):
    return self.heap[0]

cdef heapremove(heap* self):
    cdef heapitem t
    cdef int i, j, k, l

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

cdef heapitem heappop(heap* self):
    cdef heapitem it
    it = heappeek(self)
    heapremove(self)
    return it





# utility functions
cdef inline double dmax(double x, double y):
    if x>y:
        return x
    else:
        return y
cdef inline double dabs(double x):
    if x>0:
        return x
    else:
        return -x
cdef inline double _distance_p(double*x,double*y,double p,int k,double upperbound):
    """Compute the distance between x and y

    Computes the Minkowski p-distance to the power p between two points.
    If the distance**p is larger than upperbound, then any number larger
    than upperbound may be returned (the calculation is truncated).
    """
    cdef int i
    cdef double r
    r = 0
    if p==infinity:
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



# Tree structure
cdef struct innernode:
    int split_dim
    int n_points
    double split
    innernode* less
    innernode* greater
cdef struct leafnode:
    int split_dim
    int n_points
    int start_idx
    int end_idx

# this is the standard trick for variable-size arrays:
# malloc sizeof(nodeinfo)+self.m*sizeof(double) bytes.
cdef struct nodeinfo:
    innernode* node
    double side_distances[0]

cdef class cKDTree:
    """kd-tree for quick nearest-neighbor lookup

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
    """

    cdef innernode* tree 
    cdef readonly object data
    cdef double* raw_data
    cdef readonly int n, m
    cdef readonly int leafsize
    cdef readonly object maxes
    cdef double* raw_maxes
    cdef readonly object mins
    cdef double* raw_mins
    cdef object indices
    cdef np.int32_t* raw_indices
    def __init__(cKDTree self, data, int leafsize=10):
        """Construct a kd-tree.

        Parameters:
        ===========

        data : array-like, shape (n,m)
            The n data points of dimension mto be indexed. This array is 
            not copied unless this is necessary to produce a contiguous 
            array of doubles, and so modifying this data will result in 
            bogus results.
        leafsize : positive integer
            The number of points at which the algorithm switches over to
            brute-force.
        """
        cdef np.ndarray[double, ndim=2] inner_data
        cdef np.ndarray[double, ndim=1] inner_maxes
        cdef np.ndarray[double, ndim=1] inner_mins
        cdef np.ndarray[np.int32_t, ndim=1] inner_indices
        self.data = np.ascontiguousarray(data,dtype=np.float)
        self.n, self.m = np.shape(self.data)
        self.leafsize = leafsize
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.maxes = np.ascontiguousarray(np.amax(self.data,axis=0))
        self.mins = np.ascontiguousarray(np.amin(self.data,axis=0))
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=np.int32))

        inner_data = self.data
        self.raw_data = <double*>inner_data.data
        inner_maxes = self.maxes
        self.raw_maxes = <double*>inner_maxes.data
        inner_mins = self.mins
        self.raw_mins = <double*>inner_mins.data
        inner_indices = self.indices
        self.raw_indices = <np.int32_t*>inner_indices.data

        self.tree = self.__build(0, self.n, self.raw_maxes, self.raw_mins)

    cdef innernode* __build(cKDTree self, int start_idx, int end_idx, double* maxes, double* mins):
        cdef leafnode* n
        cdef innernode* ni
        cdef int i, j, t, p, q, d
        cdef double size, split, minval, maxval
        cdef double*mids
        if end_idx-start_idx<=self.leafsize:
            n = <leafnode*>stdlib.malloc(sizeof(leafnode))
            n.split_dim = -1
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
                n.split_dim = -1
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

            mids = <double*>stdlib.malloc(sizeof(double)*self.m)
            for i in range(self.m):
                mids[i] = maxes[i]
            mids[d] = split
            ni.less = self.__build(start_idx,p,mids,mins)

            for i in range(self.m):
                mids[i] = mins[i]
            mids[d] = split
            ni.greater = self.__build(p,end_idx,maxes,mids)

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
        if <int>(self.tree) == 0:
            # should happen only if __init__ was never called
            return
        self.__free_tree(self.tree)

    cdef void __query(cKDTree self, 
            double*result_distances, 
            int*result_indices, 
            double*x, 
            int k, 
            double eps, 
            double p, 
            double distance_upper_bound):
        cdef heap q
        cdef heap neighbors

        cdef int i, j
        cdef double t
        cdef nodeinfo* inf
        cdef nodeinfo* inf2
        cdef double d
        cdef double epsfac
        cdef double min_distance
        cdef double far_min_distance
        cdef heapitem it, it2, neighbor
        cdef leafnode* node
        cdef innernode* inode
        cdef innernode* near
        cdef innernode* far
        cdef double* side_distances

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
        inf = <nodeinfo*>stdlib.malloc(sizeof(nodeinfo)+self.m*sizeof(double)) 
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
                if q.n==0:
                    # no more nodes to visit
                    break
                else:
                    it = heappop(&q)
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
                    # free all the nodes still on the heap
                    for i in range(q.n):
                        stdlib.free(q.heap[i].contents.ptrdata)
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
                inf2 = <nodeinfo*>stdlib.malloc(sizeof(nodeinfo)+self.m*sizeof(double)) 
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
                    far_min_distance = min_distance - inf.side_distances[inode.split_dim] + inf2.side_distances[inode.split_dim]
                else:
                    inf2.side_distances[inode.split_dim] = dabs(inode.split-x[inode.split_dim])**p
                    far_min_distance = min_distance - inf.side_distances[inode.split_dim] + inf2.side_distances[inode.split_dim]

                it2.priority = far_min_distance


                # far child might be too far, if so, don't bother pushing it
                if far_min_distance<=distance_upper_bound*epsfac:
                    heappush(&q,it2)
                else:
                    stdlib.free(inf2)
                    # just in case
                    it2.contents.ptrdata = <char*> 0

        # fill output arrays with sorted neighbors 
        for i in range(neighbors.n-1,-1,-1):
            neighbor = heappop(&neighbors) # FIXME: neighbors may be realloced
            result_indices[i] = neighbor.contents.intdata
            if p==1 or p==infinity:
                result_distances[i] = -neighbor.priority
            else:
                result_distances[i] = (-neighbor.priority)**(1./p)

        heapdestroy(&q)
        heapdestroy(&neighbors)

    def query(cKDTree self, object x, int k=1, double eps=0, double p=2, 
            double distance_upper_bound=infinity):
        """query the kd-tree for nearest neighbors

        Parameters:
        ===========

        x : array-like, last dimension self.m
            An array of points to query.
        k : integer
            The number of nearest neighbors to return.
        eps : nonnegative float
            Return approximate nearest neighbors; the kth returned value 
            is guaranteed to be no further than (1+eps) times the 
            distance to the real kth nearest neighbor.
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use. 
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        distance_upper_bound : nonnegative float
            Return only neighbors within this distance. This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.

        Returns:
        ========
        
        d : array of floats
            The distances to the nearest neighbors. 
            If x has shape tuple+(self.m,), then d has shape tuple+(k,).
            Missing neighbors are indicated with infinite distances.
        i : array of integers
            The locations of the neighbors in self.data.
            If x has shape tuple+(self.m,), then i has shape tuple+(k,).
            Missing neighbors are indicated with self.n.
        """
        cdef np.ndarray[int, ndim=2] ii
        cdef np.ndarray[double, ndim=2] dd
        cdef np.ndarray[double, ndim=2] xx
        cdef int c
        x = np.asarray(x).astype(np.float)
        if np.shape(x)[-1] != self.m:
            raise ValueError("x must consist of vectors of length %d but has shape %s" % (self.m, np.shape(x)))
        if p<1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        if len(x.shape)==1:
            single = True
            x = x[np.newaxis,:]
        else:
            single = False
        retshape = np.shape(x)[:-1]
        n = np.prod(retshape)
        xx = np.reshape(x,(n,self.m))
        xx = np.ascontiguousarray(xx)
        dd = np.empty((n,k),dtype=np.float)
        dd.fill(infinity)
        ii = np.empty((n,k),dtype='i')
        ii.fill(self.n)
        for c in range(n):
            self.__query(
                    (<double*>dd.data)+c*k,
                    (<int*>ii.data)+c*k,
                    (<double*>xx.data)+c*self.m, 
                    k, 
                    eps,
                    p, 
                    distance_upper_bound)
        if single:
            if k==1:
                return dd[0,0], ii[0,0]
            else:
                return dd[0], ii[0]
        else:
            if k==1:
                return np.reshape(dd[...,0],retshape), np.reshape(ii[...,0],retshape)
            else:
                return np.reshape(dd,retshape+(k,)), np.reshape(ii,retshape+(k,))

