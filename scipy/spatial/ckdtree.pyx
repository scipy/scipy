# Copyright Anne M. Archibald 2008
# Additional contributions by Patrick Varilly and Sturla Molden
# Released under the scipy license
import numpy as np
import scipy.sparse

cimport numpy as np
cimport libc.stdlib as stdlib
cimport cython

cdef extern from "limits.h":
    long LONG_MAX
cdef np.float64_t infinity = np.inf

__all__ = ['cKDTree']


# Notes on int and 64-bit cleanliness
# ===================================
#
# Never use a bare "int" for array indices; use np.intp_t instead.  A Python
# int and np.int is a C long on Python 2.x, which can be 32 bits on 64-bit
# systems (e.g. Windows).
#
# The exception is as the return type of a nominally void function, which
# instead returns 0 and signals a Python exception by returning -1.
#
# Also, when converting np.intp_t's to Python, you should explicitly cast to
# a Python "int" object if sizeof(long) < sizeof(np.intp_t).  From the
# mailing list (Sturla Molden): "On Win 64 we should use Python long instead
# of Python int if a C long (i.e. Python int) overflows, which the function
# int() will ensure.  Cython automatically converts np.npy_intp [==
# np.intp_t] to Python long on Win 64, which we want to convert to a Python
# int if it is possible.  On other platforms we don't want this extra
# overhead."

# The following utility functions help properly add int tuples to sets and
# ints to lists.  The results of the if is known at compile time, so the
# test is optimized away.

cdef inline int set_add_ordered_pair(set results,
                                     np.intp_t i, np.intp_t j) except -1:
    if i > j:
        i, j = j, i
    if sizeof(long) < sizeof(np.intp_t):
        # Win 64
        results.add((int(i), int(j)))
    else:
        # Other platforms
        results.add((i, j))
    return 0

cdef inline int list_append(list results, np.intp_t i) except -1:
    if sizeof(long) < sizeof(np.intp_t):
        # Win 64
        if i <= <np.intp_t>LONG_MAX:  # CHECK COMPARISON DIRECTION
            results.append(int(i))
        else:
            results.append(i)
    else:
        # Other platforms
        results.append(i)
    return 0
    


# Priority queue
# ==============
cdef union heapcontents:    # FIXME: Unions are not always portable, verify this 
    np.intp_t intdata     # union is never used in an ABI dependent way.
    char* ptrdata

cdef struct heapitem:
    np.float64_t priority
    heapcontents contents

cdef class heap(object):
    cdef np.intp_t n
    cdef heapitem* heap
    cdef np.intp_t space
    
    def __init__(heap self, np.intp_t initial_size):
        cdef void *tmp
        self.space = initial_size
        self.heap = <heapitem*> NULL
        tmp = stdlib.malloc(sizeof(heapitem)*self.space)
        if tmp == NULL:
            raise MemoryError
        self.heap = <heapitem*> tmp  
        self.n = 0

    def __dealloc__(heap self):
        if self.heap != <heapitem*> NULL:
            stdlib.free(self.heap)

    cdef inline int _resize(heap self, np.intp_t new_space) except -1:
        cdef void *tmp
        if new_space < self.n:
            raise ValueError("Heap containing %d items cannot be resized to %d" % (int(self.n), int(new_space)))
        self.space = new_space
        tmp = stdlib.realloc(<void*>self.heap, new_space*sizeof(heapitem))
        if tmp == NULL:
            raise MemoryError
        self.heap = <heapitem*> tmp
        return 0

    @cython.cdivision(True)
    cdef inline int push(heap self, heapitem item) except -1:
        cdef np.intp_t i
        cdef heapitem t

        self.n += 1
        if self.n > self.space:
            self._resize(2 * self.space + 1)
            
        i = self.n - 1
        self.heap[i] = item
        
        while i > 0 and self.heap[i].priority < self.heap[(i - 1) // 2].priority:
            t = self.heap[(i - 1) // 2]
            self.heap[(i - 1) // 2] = self.heap[i]
            self.heap[i] = t
            i = (i - 1) // 2
        return 0
    
    
    cdef heapitem peek(heap self):
        return self.heap[0]
    
    
    @cython.cdivision(True)
    cdef int remove(heap self) except -1:
        cdef heapitem t
        cdef np.intp_t i, j, k, l
    
        self.heap[0] = self.heap[self.n-1]
        self.n -= 1
        # No point in freeing up space as the heap empties.
        # The whole heap gets deallocated at the end of any query below
        #if self.n < self.space//4 and self.space>40: #FIXME: magic number
        #    self._resize(self.space // 2 + 1)
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
    
    cdef int pop(heap self, heapitem *it) except -1:
        it[0] = self.peek()
        self.remove()
        return 0


# Utility functions
# =================
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
        
# Utility for building a coo matrix incrementally
cdef class coo_entries:
    cdef:
        np.intp_t n, n_max
        np.ndarray i, j
        np.ndarray v
        np.intp_t *i_data
        np.intp_t *j_data
        np.float64_t *v_data
    
    def __init__(self):
        self.n = 0
        self.n_max = 10
        self.i = np.empty(self.n_max, dtype=np.intp)
        self.j = np.empty(self.n_max, dtype=np.intp)
        self.v = np.empty(self.n_max, dtype=np.float64)
        self.i_data = <np.intp_t *>np.PyArray_DATA(self.i)
        self.j_data = <np.intp_t *>np.PyArray_DATA(self.j)
        self.v_data = <np.float64_t*>np.PyArray_DATA(self.v)

    cdef int add(coo_entries self, np.intp_t i, np.intp_t j,
                 np.float64_t v) except -1:
        cdef np.intp_t k
        if self.n == self.n_max:
            self.n_max *= 2
            self.i.resize(self.n_max)
            self.j.resize(self.n_max)
            self.v.resize(self.n_max)
            self.i_data = <np.intp_t *>np.PyArray_DATA(self.i)
            self.j_data = <np.intp_t *>np.PyArray_DATA(self.j)
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
        self.i_data = <np.intp_t *>np.PyArray_DATA(self.i)
        self.j_data = <np.intp_t *>np.PyArray_DATA(self.j)
        self.v_data = <np.float64_t*>np.PyArray_DATA(self.v)
        self.n_max = self.n
        return scipy.sparse.coo_matrix((self.v, (self.i, self.j)), shape=shape)


# Classes derived from this base class are intended to efficiently accumulate
# or summarize np.intp_t pairs, for example point indices in the kdtree.
cdef class pair_accumulator:
    cdef void add_pair(self, np.intp_t i, np.intp_t j):
        pass


cdef class pair_counter(pair_accumulator):
    cdef np.intp_t _n

    def __init__(self):
        self._n = 0

    cdef void add_pair(self, np.intp_t i, np.intp_t j):
        self._n += 1

    property n:
        def __get__(self):
            return self._n


cdef class pair_array_accumulator(pair_accumulator):
    cdef:
        np.intp_t n, capacity
        np.ndarray i, j
        np.intp_t *i_data
        np.intp_t *j_data

    def __init__(self):
        self.n = 0
        self.capacity = 10
        self.i = np.empty(self.capacity, dtype=np.intp)
        self.j = np.empty(self.capacity, dtype=np.intp)
        self.i_data = <np.intp_t *>np.PyArray_DATA(self.i)
        self.j_data = <np.intp_t *>np.PyArray_DATA(self.j)

    cdef void add_pair(pair_array_accumulator self, np.intp_t i, np.intp_t j):
        cdef np.intp_t k
        cdef np.intp_t tmp
        if self.n == self.capacity:
            self.capacity *= 2
            self.i.resize(self.capacity)
            self.j.resize(self.capacity)
            self.i_data = <np.intp_t *>np.PyArray_DATA(self.i)
            self.j_data = <np.intp_t *>np.PyArray_DATA(self.j)
        k = self.n
        if i < j:
            self.i_data[k] = i
            self.j_data[k] = j
        else:
            self.i_data[k] = j
            self.j_data[k] = i
        self.n += 1

    def toarray(pair_array_accumulator self):
        # Shrink arrays to size
        self.i.resize(self.n)
        self.j.resize(self.n)
        self.i_data = <np.intp_t *>np.PyArray_DATA(self.i)
        self.j_data = <np.intp_t *>np.PyArray_DATA(self.j)
        self.capacity = self.n
        return np.vstack((self.i, self.j)).T


cdef class pair_set_accumulator(pair_accumulator):
    cdef set _pairs

    def __cinit__(self):
        self._pairs = set()

    cdef void add_pair(pair_set_accumulator self, np.intp_t i, np.intp_t j):
        set_add_ordered_pair(self._pairs, i, j)

    property pairs:
        def __get__(self):
            return self._pairs



# Measuring distances
# ===================
cdef inline np.float64_t _distance_p(np.float64_t *x, np.float64_t *y,
                                     np.float64_t p, np.intp_t k,
                                     np.float64_t upperbound):
    """Compute the distance between x and y

    Computes the Minkowski p-distance to the power p between two points.
    If the distance**p is larger than upperbound, then any number larger
    than upperbound may be returned (the calculation is truncated).
    """
    cdef np.intp_t i
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
# ===================

cdef class Rectangle:
    cdef np.intp_t m
    cdef np.float64_t *mins
    cdef np.float64_t *maxes
    cdef np.ndarray mins_arr, maxes_arr

    def __init__(self, mins_arr, maxes_arr):
        # Copy array data
        self.mins_arr = np.array(mins_arr, dtype=np.float64, order='C')
        self.maxes_arr = np.array(maxes_arr, dtype=np.float64, order='C')
        self.mins = <np.float64_t*>np.PyArray_DATA(self.mins_arr)
        self.maxes = <np.float64_t*>np.PyArray_DATA(self.maxes_arr)
        self.m = self.mins_arr.shape[0]

# 1-d pieces
# These should only be used if p != infinity
cdef inline np.float64_t min_dist_point_interval_p(np.float64_t* x,
                                                   Rectangle rect,
                                                   np.intp_t k,
                                                   np.float64_t p):    
    """Compute the minimum distance along dimension k between x and
    a point in the hyperrectangle.
    """
    return dmax(0, dmax(rect.mins[k] - x[k], x[k] - rect.maxes[k])) ** p

cdef inline np.float64_t max_dist_point_interval_p(np.float64_t* x,
                                                   Rectangle rect,
                                                   np.intp_t k,
                                                   np.float64_t p):
    """Compute the maximum distance along dimension k between x and
    a point in the hyperrectangle.
    """
    return dmax(rect.maxes[k] - x[k], x[k] - rect.mins[k]) ** p

cdef inline np.float64_t min_dist_interval_interval_p(Rectangle rect1,
                                                      Rectangle rect2,
                                                      np.intp_t k,
                                                      np.float64_t p):
    """Compute the minimum distance along dimension k between points in
    two hyperrectangles.
    """
    return dmax(0, dmax(rect1.mins[k] - rect2.maxes[k],
                        rect2.mins[k] - rect1.maxes[k])) ** p

cdef inline np.float64_t max_dist_interval_interval_p(Rectangle rect1,
                                                      Rectangle rect2,
                                                      np.intp_t k,
                                                      np.float64_t p):
    """Compute the maximum distance along dimension k between points in
    two hyperrectangles.
    """
    return dmax(rect1.maxes[k] - rect2.mins[k], rect2.maxes[k] - rect1.mins[k]) ** p

# Interval arithmetic in m-D
# ==========================

# These should be used only for p == infinity
cdef inline np.float64_t min_dist_point_rect_p_inf(np.float64_t* x,
                                                   Rectangle rect):
    """Compute the minimum distance between x and the given hyperrectangle."""
    cdef np.intp_t i
    cdef np.float64_t min_dist = 0.
    for i in range(rect.m):
        min_dist = dmax(min_dist, dmax(rect.mins[i]-x[i], x[i]-rect.maxes[i]))
    return min_dist

cdef inline np.float64_t max_dist_point_rect_p_inf(np.float64_t* x,
                                                   Rectangle rect):
    """Compute the maximum distance between x and the given hyperrectangle."""
    cdef np.intp_t i
    cdef np.float64_t max_dist = 0.
    for i in range(rect.m):
        max_dist = dmax(max_dist, dmax(rect.maxes[i]-x[i], x[i]-rect.mins[i]))
    return max_dist

cdef inline np.float64_t min_dist_rect_rect_p_inf(Rectangle rect1,
                                                  Rectangle rect2):
    """Compute the minimum distance between points in two hyperrectangles."""
    cdef np.intp_t i
    cdef np.float64_t min_dist = 0.
    for i in range(rect1.m):
        min_dist = dmax(min_dist, dmax(rect1.mins[i] - rect2.maxes[i],
                                       rect2.mins[i] - rect1.maxes[i]))
    return min_dist

cdef inline np.float64_t max_dist_rect_rect_p_inf(Rectangle rect1,
                                                  Rectangle rect2):
    """Compute the maximum distance between points in two hyperrectangles."""
    cdef np.intp_t i
    cdef np.float64_t max_dist = 0.
    for i in range(rect1.m):
        max_dist = dmax(max_dist, dmax(rect1.maxes[i] - rect2.mins[i],
                                       rect2.maxes[i] - rect1.mins[i]))
    return max_dist

# Rectangle-to-rectangle distance tracker
# =======================================
#
# The logical unit that repeats over and over is to keep track of the
# maximum and minimum distances between points in two hyperrectangles
# as these rectangles are successively split.
#
# Example
# -------
# # node1 encloses points in rect1, node2 encloses those in rect2
#
# cdef RectRectDistanceTracker dist_tracker
# dist_tracker = RectRectDistanceTracker(rect1, rect2, p)
#
# ...
#
# if dist_tracker.min_distance < ...:
#     ...
#
# dist_tracker.push_less_of(1, node1)
# do_something(node1.less, dist_tracker)
# dist_tracker.pop()
#
# dist_tracker.push_greater_of(1, node1)
# do_something(node1.greater, dist_tracker)
# dist_tracker.pop()

cdef struct RR_stack_item:
    np.intp_t which
    np.intp_t split_dim
    double min_along_dim, max_along_dim
    np.float64_t min_distance, max_distance

cdef np.intp_t LESS = 1
cdef np.intp_t GREATER = 2

cdef class RectRectDistanceTracker(object):
    cdef Rectangle rect1, rect2
    cdef np.float64_t p, epsfac, upper_bound
    cdef np.float64_t min_distance, max_distance

    cdef np.intp_t stack_size, stack_max_size
    cdef RR_stack_item *stack

    # Stack handling
    cdef int _init_stack(self) except -1:
        cdef void *tmp
        self.stack_max_size = 10
        tmp = stdlib.malloc(sizeof(RR_stack_item) *
                            self.stack_max_size)
        if tmp == NULL:
            raise MemoryError
        self.stack = <RR_stack_item*> tmp
        self.stack_size = 0
        return 0

    cdef int _resize_stack(self, np.intp_t new_max_size) except -1:
        cdef void *tmp
        self.stack_max_size = new_max_size
        tmp = stdlib.realloc(<RR_stack_item*> self.stack,
                             new_max_size * sizeof(RR_stack_item))
        if tmp == NULL:
            raise MemoryError
        self.stack = <RR_stack_item*> tmp
        return 0
    
    cdef int _free_stack(self) except -1:
        if self.stack != <RR_stack_item*> NULL:
            stdlib.free(self.stack)
        return 0
    

    def __init__(self, Rectangle rect1, Rectangle rect2,
                 np.float64_t p, np.float64_t eps, np.float64_t upper_bound):
        
        if rect1.m != rect2.m:
            raise ValueError("rect1 and rect2 have different dimensions")

        self.rect1 = rect1
        self.rect2 = rect2
        self.p = p
        
        # internally we represent all distances as distance ** p
        if p != infinity and upper_bound != infinity:
            self.upper_bound = upper_bound ** p
        else:
            self.upper_bound = upper_bound

        # fiddle approximation factor
        if eps == 0:
            self.epsfac = 1
        elif p == infinity:
            self.epsfac = 1 / (1 + eps)
        else:
            self.epsfac = 1 / (1 + eps) ** p

        self._init_stack()

        # Compute initial min and max distances
        if self.p == infinity:
            self.min_distance = min_dist_rect_rect_p_inf(rect1, rect2)
            self.max_distance = max_dist_rect_rect_p_inf(rect1, rect2)
        else:
            self.min_distance = 0.
            self.max_distance = 0.
            for i in range(rect1.m):
                self.min_distance += min_dist_interval_interval_p(rect1, rect2, i, p)
                self.max_distance += max_dist_interval_interval_p(rect1, rect2, i, p)

    def __dealloc__(self):
        self._free_stack()

    cdef int push(self, np.intp_t which, np.intp_t direction,
                  np.intp_t split_dim,
                  np.float64_t split_val) except -1:

        cdef Rectangle rect
        if which == 1:
            rect = self.rect1
        else:
            rect = self.rect2

        # Push onto stack
        if self.stack_size == self.stack_max_size:
            self._resize_stack(self.stack_max_size * 2)
            
        cdef RR_stack_item *item = &self.stack[self.stack_size]
        self.stack_size += 1
        item.which = which
        item.split_dim = split_dim
        item.min_distance = self.min_distance
        item.max_distance = self.max_distance
        item.min_along_dim = rect.mins[split_dim]
        item.max_along_dim = rect.maxes[split_dim]

        # Update min/max distances
        if self.p != infinity:
            self.min_distance -= min_dist_interval_interval_p(self.rect1, self.rect2, split_dim, self.p)
            self.max_distance -= max_dist_interval_interval_p(self.rect1, self.rect2, split_dim, self.p)

        if direction == LESS:
            rect.maxes[split_dim] = split_val
        else:
            rect.mins[split_dim] = split_val

        if self.p != infinity:
            self.min_distance += min_dist_interval_interval_p(self.rect1, self.rect2, split_dim, self.p)
            self.max_distance += max_dist_interval_interval_p(self.rect1, self.rect2, split_dim, self.p)
        else:
            self.min_distance = min_dist_rect_rect_p_inf(self.rect1, self.rect2)
            self.max_distance = max_dist_rect_rect_p_inf(self.rect1, self.rect2)
            
        return 0

    
    cdef inline int push_less_of(self, np.intp_t which,
                                 innernode *node) except -1:
        return self.push(which, LESS, node.split_dim, node.split)

    
    cdef inline int push_greater_of(self, np.intp_t which,
                                    innernode *node) except -1:
        return self.push(which, GREATER, node.split_dim, node.split)

    
    cdef inline int pop(self) except -1:
        # Pop from stack
        self.stack_size -= 1
        assert self.stack_size >= 0
        
        cdef RR_stack_item* item = &self.stack[self.stack_size]
        self.min_distance = item.min_distance
        self.max_distance = item.max_distance

        if item.which == 1:
            self.rect1.mins[item.split_dim] = item.min_along_dim
            self.rect1.maxes[item.split_dim] = item.max_along_dim
        else:
            self.rect2.mins[item.split_dim] = item.min_along_dim
            self.rect2.maxes[item.split_dim] = item.max_along_dim
        
        return 0

# Point-to-rectangle distance tracker
# ===================================
#
# The other logical unit that is used in query_ball_point is to keep track
# of the maximum and minimum distances between points in a hyperrectangle
# and another fixed point as the rectangle is successively split.
#
# Example
# -------
# # node encloses points in rect
#
# cdef PointRectDistanceTracker dist_tracker
# dist_tracker = PointRectDistanceTracker(pt, rect, p)
#
# ...
#
# if dist_tracker.min_distance < ...:
#     ...
#
# dist_tracker.push_less_of(node)
# do_something(node.less, dist_tracker)
# dist_tracker.pop()
#
# dist_tracker.push_greater_of(node)
# do_something(node.greater, dist_tracker)
# dist_tracker.pop()

cdef struct RP_stack_item:
    np.intp_t split_dim
    double min_along_dim, max_along_dim
    np.float64_t min_distance, max_distance

cdef class PointRectDistanceTracker(object):
    cdef Rectangle rect
    cdef np.float64_t *pt
    cdef np.float64_t p, epsfac, upper_bound
    cdef np.float64_t min_distance, max_distance

    cdef np.intp_t stack_size, stack_max_size
    cdef RP_stack_item *stack

    # Stack handling
    cdef int _init_stack(self) except -1:
        cdef void *tmp
        self.stack_max_size = 10
        tmp = stdlib.malloc(sizeof(RP_stack_item) *
                            self.stack_max_size)
        if tmp == NULL:
            raise MemoryError
        self.stack = <RP_stack_item*> tmp
        self.stack_size = 0
        return 0

    cdef int _resize_stack(self, np.intp_t new_max_size) except -1:
        cdef void *tmp
        self.stack_max_size = new_max_size
        tmp = stdlib.realloc(<RP_stack_item*> self.stack,
                              new_max_size * sizeof(RP_stack_item))
        if tmp == NULL:
            raise MemoryError
        self.stack = <RP_stack_item*> tmp
        return 0
    
    cdef int _free_stack(self) except -1:
        if self.stack != <RP_stack_item*> NULL:
            stdlib.free(self.stack)
        return 0

    cdef init(self, np.float64_t *pt, Rectangle rect,
              np.float64_t p, np.float64_t eps, np.float64_t upper_bound):

        self.pt = pt
        self.rect = rect
        self.p = p
        
        # internally we represent all distances as distance ** p
        if p != infinity and upper_bound != infinity:
            self.upper_bound = upper_bound ** p
        else:
            self.upper_bound = upper_bound

        # fiddle approximation factor
        if eps == 0:
            self.epsfac = 1
        elif p == infinity:
            self.epsfac = 1 / (1 + eps)
        else:
            self.epsfac = 1 / (1 + eps) ** p

        self._init_stack()

        # Compute initial min and max distances
        if self.p == infinity:
            self.min_distance = min_dist_point_rect_p_inf(pt, rect)
            self.max_distance = max_dist_point_rect_p_inf(pt, rect)
        else:
            self.min_distance = 0.
            self.max_distance = 0.
            for i in range(rect.m):
                self.min_distance += min_dist_point_interval_p(pt, rect, i, p)
                self.max_distance += max_dist_point_interval_p(pt, rect, i, p)

    def __dealloc__(self):
        self._free_stack()

    cdef int push(self, np.intp_t direction,
                  np.intp_t split_dim,
                  np.float64_t split_val) except -1:

        # Push onto stack
        if self.stack_size == self.stack_max_size:
            self._resize_stack(self.stack_max_size * 2)
            
        cdef RP_stack_item *item = &self.stack[self.stack_size]
        self.stack_size += 1
        
        item.split_dim = split_dim
        item.min_distance = self.min_distance
        item.max_distance = self.max_distance
        item.min_along_dim = self.rect.mins[split_dim]
        item.max_along_dim = self.rect.maxes[split_dim]
            
        if self.p != infinity:
            self.min_distance -= min_dist_point_interval_p(self.pt, self.rect, split_dim, self.p)
            self.max_distance -= max_dist_point_interval_p(self.pt, self.rect, split_dim, self.p)

        if direction == LESS:
            self.rect.maxes[split_dim] = split_val
        else:
            self.rect.mins[split_dim] = split_val

        if self.p != infinity:
            self.min_distance += min_dist_point_interval_p(self.pt, self.rect, split_dim, self.p)
            self.max_distance += max_dist_point_interval_p(self.pt, self.rect, split_dim, self.p)
        else:
            self.min_distance = min_dist_point_rect_p_inf(self.pt, self.rect)
            self.max_distance = max_dist_point_rect_p_inf(self.pt, self.rect)
            
        return 0

    
    cdef inline int push_less_of(self, innernode* node) except -1:
        return self.push(LESS, node.split_dim, node.split)

    
    cdef inline int push_greater_of(self, innernode* node) except -1:
        return self.push(GREATER, node.split_dim, node.split)

    
    cdef inline int pop(self) except -1:
        self.stack_size -= 1
        assert self.stack_size >= 0
        
        cdef RP_stack_item* item = &self.stack[self.stack_size]
        self.min_distance = item.min_distance
        self.max_distance = item.max_distance
        self.rect.mins[item.split_dim] = item.min_along_dim
        self.rect.maxes[item.split_dim] = item.max_along_dim
        
        return 0

# Tree structure
# ==============
cdef struct innernode:
    np.intp_t split_dim
    np.intp_t children
    np.float64_t split
    innernode* less
    innernode* greater
    
cdef struct leafnode:
    np.intp_t split_dim
    np.intp_t children
    np.intp_t start_idx
    np.intp_t end_idx


# this is the standard trick for variable-size arrays:
# malloc sizeof(nodeinfo)+self.m*sizeof(np.float64_t) bytes.

cdef struct nodeinfo:
    innernode* node
    np.float64_t side_distances[0]  # FIXME: Only valid in C99, invalid C++ and C89


# Main class
# ==========
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
    data : array_like, shape (n,m)
        The n data points of dimension m to be indexed. This array is 
        not copied unless this is necessary to produce a contiguous 
        array of doubles, and so modifying this data will result in 
        bogus results.
    leafsize : positive integer
        The number of points at which the algorithm switches over to
        brute-force.

    """

    cdef innernode* tree 
    cdef readonly np.ndarray data
    cdef np.float64_t* raw_data
    cdef readonly np.intp_t n, m
    cdef readonly np.intp_t leafsize
    cdef readonly np.ndarray maxes
    cdef np.float64_t* raw_maxes
    cdef readonly np.ndarray mins
    cdef np.float64_t* raw_mins
    cdef np.ndarray indices
    cdef np.intp_t* raw_indices

    def __init__(cKDTree self, data, np.intp_t leafsize=10):
        cdef np.ndarray[np.float64_t, ndim=2] data_arr = \
            np.ascontiguousarray(data, dtype=np.float64)
        self.data = data_arr
        self.n = data_arr.shape[0]
        self.m = data_arr.shape[1]
        self.leafsize = leafsize
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.maxes = np.ascontiguousarray(np.amax(self.data,axis=0), dtype=np.float64)
        self.mins = np.ascontiguousarray(np.amin(self.data,axis=0), dtype=np.float64)
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=np.intp))

        self.raw_data = <np.float64_t*>np.PyArray_DATA(self.data)
        self.raw_maxes = <np.float64_t*>np.PyArray_DATA(self.maxes)
        self.raw_mins = <np.float64_t*>np.PyArray_DATA(self.mins)
        self.raw_indices = <np.intp_t*>np.PyArray_DATA(self.indices)

        self.tree = self.__build(0, self.n, self.raw_maxes, self.raw_mins)

    cdef innernode* __build(cKDTree self, np.intp_t start_idx, np.intp_t end_idx,
                            np.float64_t* maxes, np.float64_t* mins) except? <innernode*> NULL:
        cdef leafnode* n
        cdef innernode* ni
        cdef np.intp_t i, j, t, p, q, d
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
        if <np.intp_t>(self.tree) == 0:
            # should happen only if __init__ was never called
            return
        self.__free_tree(self.tree)

    # -----
    # query
    # -----

    cdef int __query(cKDTree self, 
            np.float64_t*result_distances, 
            np.intp_t*result_indices, 
            np.float64_t*x, 
            np.intp_t k, 
            np.float64_t eps, 
            np.float64_t p, 
            np.float64_t distance_upper_bound) except -1:

        cdef heap q
        cdef heap neighbors

        cdef np.intp_t i, j
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


        # priority queue for chasing nodes
        # entries are:
        #  minimum distance between the cell and the target
        #  distances between the nearest side of the cell and the target
        #  the head node of the cell
        q = heap(12)

        # priority queue for the nearest neighbors
        # furthest known neighbor first
        # entries are (-distance**p, i)
        neighbors = heap(k)

        inf = inf2 = <nodeinfo*> NULL    

        try:
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
                                neighbors.remove()
                            neighbor.priority = -d
                            neighbor.contents.intdata = self.raw_indices[i]
                            neighbors.push(neighbor)

                            # adjust upper bound for efficiency
                            if neighbors.n==k:
                                distance_upper_bound = -neighbors.peek().priority
                    
                    # done with this node, get another
                    stdlib.free(inf)
                    inf = <nodeinfo*> NULL

                    if q.n==0:
                        # no more nodes to visit
                        break
                    else:
                        q.pop(&it)
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
                        q.push(it2)
                    else:
                        stdlib.free(inf2)
                        inf2 = <nodeinfo*> NULL
                        # just in case
                        it2.contents.ptrdata = <char*> NULL

            # fill output arrays with sorted neighbors 
            for i in range(neighbors.n-1,-1,-1):
                neighbors.pop(&neighbor)
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

        return 0


    @cython.boundscheck(False)
    def query(cKDTree self, object x, np.intp_t k=1, np.float64_t eps=0,
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
        cdef np.ndarray[np.intp_t, ndim=2] ii
        cdef np.ndarray[np.float64_t, ndim=2] dd
        cdef np.ndarray[np.float64_t, ndim=2] xx
        cdef np.intp_t c, n, i, j
        cdef np.ndarray x_arr = np.asarray(x, dtype=np.float64)
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
        for c in range(n):
            self.__query(&dd[c, 0], &ii[c, 0], &xx[c, 0],
                         k, eps, p, distance_upper_bound)

        if single:
            if k==1:
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
    cdef int __query_ball_point_traverse_no_checking(cKDTree self,
                                                     list results,
                                                     innernode* node) except -1:
        cdef leafnode* lnode
        cdef np.intp_t i

        if node.split_dim == -1:  # leaf node
            lnode = <leafnode*> node
            for i in range(lnode.start_idx, lnode.end_idx):
                list_append(results, self.raw_indices[i])
        else:
            self.__query_ball_point_traverse_no_checking(results, node.less)
            self.__query_ball_point_traverse_no_checking(results, node.greater)

        return 0


    @cython.cdivision(True)
    cdef int __query_ball_point_traverse_checking(cKDTree self,
                                                  list results,
                                                  innernode* node,
                                                  PointRectDistanceTracker tracker) except -1:
        cdef leafnode* lnode
        cdef np.float64_t d
        cdef np.intp_t i

        if tracker.min_distance > tracker.upper_bound * tracker.epsfac:
            return 0
        elif tracker.max_distance < tracker.upper_bound / tracker.epsfac:
            self.__query_ball_point_traverse_no_checking(results, node)
        elif node.split_dim == -1:  # leaf node
            lnode = <leafnode*>node
            # brute-force
            for i in range(lnode.start_idx, lnode.end_idx):
                d = _distance_p(
                    self.raw_data + self.raw_indices[i] * self.m,
                    tracker.pt, tracker.p, self.m, tracker.upper_bound)
                if d <= tracker.upper_bound:
                    list_append(results, self.raw_indices[i])
        else:
            tracker.push_less_of(node)
            self.__query_ball_point_traverse_checking(
                results, node.less, tracker)
            tracker.pop()
            
            tracker.push_greater_of(node)
            self.__query_ball_point_traverse_checking(
                results, node.greater, tracker)
            tracker.pop()
            
        return 0


    cdef list __query_ball_point(cKDTree self,
                                 np.float64_t* x,
                                 np.float64_t r,
                                 np.float64_t p,
                                 np.float64_t eps):

        tracker = PointRectDistanceTracker()
        tracker.init(x, Rectangle(self.mins, self.maxes),
                     p, eps, r)
        
        results = []
        self.__query_ball_point_traverse_checking(
            results, self.tree, tracker)
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
        
        x = np.asarray(x, dtype=np.float64)
        if x.shape[-1] != self.m:
            raise ValueError("Searching for a %d-dimensional point in a " \
                             "%d-dimensional KDTree" % (int(x.shape[-1]), int(self.m)))
        if len(x.shape) == 1:
            xx = np.ascontiguousarray(x, dtype=np.float64)
            return self.__query_ball_point(&xx[0], r, p, eps)
        else:
            retshape = x.shape[:-1]
            result = np.empty(retshape, dtype=object)
            for c in np.ndindex(retshape):
                xx = np.ascontiguousarray(x[c], dtype=np.float64)
                result[c] = self.__query_ball_point(&xx[0], r, p, eps)
            return result

    # ---------------
    # query_ball_tree
    # ---------------
    cdef int __query_ball_tree_traverse_no_checking(cKDTree self,
                                                    cKDTree other,
                                                    list results,
                                                    innernode* node1,
                                                    innernode* node2) except -1:
        cdef leafnode *lnode1
        cdef leafnode *lnode2
        cdef list results_i
        cdef np.intp_t i, j
        
        if node1.split_dim == -1:  # leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # leaf node
                lnode2 = <leafnode*>node2
                
                for i in range(lnode1.start_idx, lnode1.end_idx):
                    results_i = results[self.raw_indices[i]]
                    for j in range(lnode2.start_idx, lnode2.end_idx):
                        list_append(results_i, other.raw_indices[j])
            else:
                
                self.__query_ball_tree_traverse_no_checking(other, results, node1, node2.less)
                self.__query_ball_tree_traverse_no_checking(other, results, node1, node2.greater)
        else:
            
            self.__query_ball_tree_traverse_no_checking(other, results, node1.less, node2)
            self.__query_ball_tree_traverse_no_checking(other, results, node1.greater, node2)

        return 0


    @cython.cdivision(True)
    cdef int __query_ball_tree_traverse_checking(cKDTree self,
                                                 cKDTree other,
                                                 list results,
                                                 innernode* node1,
                                                 innernode* node2,
                                                 RectRectDistanceTracker tracker) except -1:
        cdef leafnode *lnode1
        cdef leafnode *lnode2
        cdef list results_i
        cdef np.float64_t d
        cdef np.intp_t i, j

        if tracker.min_distance > tracker.upper_bound * tracker.epsfac:
            return 0
        elif tracker.max_distance < tracker.upper_bound / tracker.epsfac:
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
                            tracker.p, self.m, tracker.upper_bound)
                        if d <= tracker.upper_bound:
                            list_append(results_i, other.raw_indices[j])
                            
            else:  # 1 is a leaf node, 2 is inner node

                tracker.push_less_of(2, node2)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1, node2.greater, tracker)
                tracker.pop()
            
                
        else:  # 1 is an inner node
            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                tracker.push_less_of(1, node1)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1.less, node2, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(1, node1)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1.greater, node2, tracker)
                tracker.pop()
                
            else: # 1 & 2 are inner nodes
                
                tracker.push_less_of(1, node1)
                tracker.push_less_of(2, node2)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1.less, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1.less, node2.greater, tracker)
                tracker.pop()
                tracker.pop()

                
                tracker.push_greater_of(1, node1)
                tracker.push_less_of(2, node2)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1.greater, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__query_ball_tree_traverse_checking(
                    other, results, node1.greater, node2.greater, tracker)
                tracker.pop()
                tracker.pop()
            
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

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to query_ball_tree have different dimensionality")

        # Track node-to-node min/max distances
        tracker = RectRectDistanceTracker(
            Rectangle(self.mins, self.maxes),
            Rectangle(other.mins, other.maxes),
            p, eps, r)
        
        results = [[] for i in range(self.n)]
        self.__query_ball_tree_traverse_checking(
            other, results, self.tree, other.tree, tracker)

        return results


    # -----------
    # query_pairs
    # -----------
    cdef int __query_pairs_traverse_no_checking(cKDTree self,
            pair_accumulator results,
            innernode* node1,
            innernode* node2) except -1:
        cdef leafnode *lnode1
        cdef leafnode *lnode2
        cdef np.intp_t i, j, min_j
        
        if node1.split_dim == -1:  # leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # leaf node
                lnode2 = <leafnode*>node2

                for i in range(lnode1.start_idx, lnode1.end_idx):
                    # Special care here to avoid duplicate pairs
                    if node1 == node2:
                        min_j = i + 1
                    else:
                        min_j = lnode2.start_idx
                        
                    for j in range(min_j, lnode2.end_idx):
                        results.add_pair(
                                self.raw_indices[i], self.raw_indices[j])
                            
            else:
                self.__query_pairs_traverse_no_checking(results,
                        node1, node2.less)
                self.__query_pairs_traverse_no_checking(results,
                        node1, node2.greater)
        else:
            if node1 == node2:
                # Avoid traversing (node1.less, node2.greater) and
                # (node1.greater, node2.less) (it's the same node pair twice
                # over, which is the source of the complication in the
                # original KDTree.query_pairs)
                self.__query_pairs_traverse_no_checking(results,
                        node1.less, node2.less)
                self.__query_pairs_traverse_no_checking(results,
                        node1.less, node2.greater)
                self.__query_pairs_traverse_no_checking(results,
                        node1.greater, node2.greater)
            else:
                self.__query_pairs_traverse_no_checking(results,
                        node1.less, node2)
                self.__query_pairs_traverse_no_checking(results,
                        node1.greater, node2)

        return 0


    @cython.cdivision(True)
    cdef int __query_pairs_traverse_checking(cKDTree self,
            pair_accumulator results,
            innernode* node1,
            innernode* node2,
            RectRectDistanceTracker tracker) except -1:
        cdef leafnode *lnode1
        cdef leafnode *lnode2
        cdef np.float64_t d
        cdef np.intp_t i, j, min_j

        if tracker.min_distance > tracker.upper_bound * tracker.epsfac:
            return 0
        elif tracker.max_distance < tracker.upper_bound / tracker.epsfac:
            self.__query_pairs_traverse_no_checking(results, node1, node2)
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = <leafnode*>node2
                
                # brute-force
                for i in range(lnode1.start_idx, lnode1.end_idx):
                    
                    # Special care here to avoid duplicate pairs
                    if node1 == node2:
                        min_j = i + 1
                    else:
                        min_j = lnode2.start_idx
                        
                    for j in range(min_j, lnode2.end_idx):
                        d = _distance_p(
                            self.raw_data + self.raw_indices[i] * self.m,
                            self.raw_data + self.raw_indices[j] * self.m,
                            tracker.p, self.m, tracker.upper_bound)
                        if d <= tracker.upper_bound:
                            results.add_pair(
                                    self.raw_indices[i], self.raw_indices[j])
                            
            else:  # 1 is a leaf node, 2 is inner node
                tracker.push_less_of(2, node2)
                self.__query_pairs_traverse_checking(
                    results, node1, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__query_pairs_traverse_checking(
                    results, node1, node2.greater, tracker)
                tracker.pop()
                
        else:  # 1 is an inner node
            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                tracker.push_less_of(1, node1)
                self.__query_pairs_traverse_checking(
                    results, node1.less, node2, tracker)
                tracker.pop()
                
                tracker.push_greater_of(1, node1)
                self.__query_pairs_traverse_checking(
                    results, node1.greater, node2, tracker)
                tracker.pop()
                
            else: # 1 and 2 are inner nodes
                tracker.push_less_of(1, node1)
                tracker.push_less_of(2, node2)
                self.__query_pairs_traverse_checking(
                    results, node1.less, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__query_pairs_traverse_checking(
                    results, node1.less, node2.greater, tracker)
                tracker.pop()
                tracker.pop()
                    
                tracker.push_greater_of(1, node1)
                if node1 != node2:
                    # Avoid traversing (node1.less, node2.greater) and
                    # (node1.greater, node2.less) (it's the same node pair
                    # twice over, which is the source of the complication in
                    # the original KDTree.query_pairs)
                    tracker.push_less_of(2, node2)
                    self.__query_pairs_traverse_checking(
                        results, node1.greater, node2.less, tracker)
                    tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__query_pairs_traverse_checking(
                    results, node1.greater, node2.greater, tracker)
                tracker.pop()
                tracker.pop()
                
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
            Set of pairs ``(i,j)``, with ``i < j``, for which the corresponding
            positions are close.

        """
        return self.query_pairs_set(r, p, eps)

    def _query_pairs(cKDTree self, pair_accumulator accum,
            np.float64_t r, np.float64_t p=2., np.float64_t eps=0):
        tracker = RectRectDistanceTracker(
            Rectangle(self.mins, self.maxes),
            Rectangle(self.mins, self.maxes),
            p, eps, r)
        self.__query_pairs_traverse_checking(
            accum, self.tree, self.tree, tracker)

    # The following three methods could be simplified?

    def query_pairs_count(cKDTree self, np.float64_t r, np.float64_t p=2.,
                    np.float64_t eps=0):
        accum = pair_counter()
        self._query_pairs(accum, r, p, eps)
        return accum.n

    def query_pairs_array(cKDTree self, np.float64_t r, np.float64_t p=2.,
                    np.float64_t eps=0):
        accum = pair_array_accumulator()
        self._query_pairs(accum, r, p, eps)
        return accum.toarray()

    def query_pairs_set(cKDTree self, np.float64_t r, np.float64_t p=2.,
                    np.float64_t eps=0):
        accum = pair_set_accumulator()
        self._query_pairs(accum, r, p, eps)
        return accum.pairs


    # ---------------
    # count_neighbors
    # ---------------
    cdef int __count_neighbors_traverse(cKDTree self,
                                        cKDTree other,
                                        np.intp_t n_queries,
                                        np.float64_t* r,
                                        np.intp_t * results,
                                        np.intp_t * idx,
                                        innernode* node1,
                                        innernode* node2,
                                        RectRectDistanceTracker tracker) except -1:
        cdef leafnode *lnode1
        cdef leafnode *lnode2
        cdef np.float64_t d
        cdef np.intp_t *old_idx
        cdef np.intp_t old_n_queries, l, i, j

        # Speed through pairs of nodes all of whose children are close
        # and see if any work remains to be done
        old_idx = idx
        cdef np.ndarray[np.intp_t, ndim=1] inner_idx
        inner_idx = np.empty((n_queries,), dtype=np.intp)
        idx = &inner_idx[0]

        old_n_queries = n_queries
        n_queries = 0
        for i in range(old_n_queries):
            if tracker.max_distance < r[old_idx[i]]:
                results[old_idx[i]] += node1.children * node2.children
            elif tracker.min_distance <= r[old_idx[i]]:
                idx[n_queries] = old_idx[i]
                n_queries += 1

        if n_queries > 0:
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
                                tracker.p, self.m, tracker.max_distance)
                            # I think it's usually cheaper to test d against all r's
                            # than to generate a distance array, sort it, then
                            # search for all r's via binary search
                            for l in range(n_queries):
                                if d <= r[idx[l]]:
                                    results[idx[l]] += 1
                                
                else:  # 1 is a leaf node, 2 is inner node
                    tracker.push_less_of(2, node2)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1, node2.less, tracker)
                    tracker.pop()

                    tracker.push_greater_of(2, node2)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1, node2.greater, tracker)
                    tracker.pop()
                
            else:  # 1 is an inner node
                if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                    tracker.push_less_of(1, node1)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1.less, node2, tracker)
                    tracker.pop()
                    
                    tracker.push_greater_of(1, node1)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1.greater, node2, tracker)
                    tracker.pop()
                    
                else: # 1 and 2 are inner nodes
                    tracker.push_less_of(1, node1)
                    tracker.push_less_of(2, node2)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1.less, node2.less, tracker)
                    tracker.pop()
                        
                    tracker.push_greater_of(2, node2)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1.less, node2.greater, tracker)
                    tracker.pop()
                    tracker.pop()
                        
                    tracker.push_greater_of(1, node1)
                    tracker.push_less_of(2, node2)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1.greater, node2.less, tracker)
                    tracker.pop()
                        
                    tracker.push_greater_of(2, node2)
                    self.__count_neighbors_traverse(
                        other, n_queries, r, results, idx,
                        node1.greater, node2.greater, tracker)
                    tracker.pop()
                    tracker.pop()
                    
        return 0

    @cython.boundscheck(False)
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
        cdef int r_ndim
        cdef np.intp_t n_queries, i
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] real_r
        cdef np.ndarray[np.intp_t, ndim=1, mode="c"] results, idx

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to count_neighbors have different dimensionality")

        # Make a copy of r array to ensure it's contiguous and to modify it
        # below
        r_ndim = len(np.shape(r))
        if r_ndim > 1:
            raise ValueError("r must be either a single value or a one-dimensional array of values")
        real_r = np.array(r, ndmin=1, dtype=np.float64, copy=True)
        n_queries = real_r.shape[0]

        # Internally, we represent all distances as distance ** p
        if p != infinity:
            for i in range(n_queries):
                if real_r[i] != infinity:
                    real_r[i] = real_r[i] ** p

        # Track node-to-node min/max distances
        tracker = RectRectDistanceTracker(
            Rectangle(self.mins, self.maxes),
            Rectangle(other.mins, other.maxes),
            p, 0.0, 0.0)
        
        # Go!
        results = np.zeros(n_queries, dtype=np.intp)
        idx = np.arange(n_queries, dtype=np.intp)
        self.__count_neighbors_traverse(other, n_queries,
                                        &real_r[0], &results[0], &idx[0],
                                        self.tree, other.tree,
                                        tracker)
        
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
    cdef int __sparse_distance_matrix_traverse(cKDTree self, cKDTree other,
                                               coo_entries results,
                                               innernode* node1, innernode* node2,
                                               RectRectDistanceTracker tracker) except -1:
        cdef leafnode *lnode1
        cdef leafnode *lnode2
        cdef list results_i
        cdef np.float64_t d
        cdef np.intp_t i, j, min_j
                
        if tracker.min_distance > tracker.upper_bound:
            return 0
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = <leafnode*>node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = <leafnode*>node2
                
                # brute-force
                for i in range(lnode1.start_idx, lnode1.end_idx):
                    # Special care here to avoid duplicate pairs
                    if node1 == node2:
                        min_j = i+1
                    else:
                        min_j = lnode2.start_idx
                        
                    for j in range(min_j, lnode2.end_idx):
                        d = _distance_p(
                            self.raw_data + self.raw_indices[i] * self.m,
                            other.raw_data + other.raw_indices[j] * self.m,
                            tracker.p, self.m, tracker.upper_bound)
                        if d <= tracker.upper_bound:
                            if tracker.p != 1 and tracker.p != infinity:
                                d = d**(1. / tracker.p)
                            results.add(self.raw_indices[i],
                                        other.raw_indices[j], d)
                            if node1 == node2:
                                results.add(self.raw_indices[j],
                                            other.raw_indices[i], d)

            else:  # 1 is a leaf node, 2 is inner node
                tracker.push_less_of(2, node2)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1, node2.greater, tracker)
                tracker.pop()
                
        else:  # 1 is an inner node
            if node2.split_dim == -1:  # 1 is an inner node, 2 is a leaf node
                tracker.push_less_of(1, node1)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1.less, node2, tracker)
                tracker.pop()
                
                tracker.push_greater_of(1, node1)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1.greater, node2, tracker)
                tracker.pop()
                
            else: # 1 and 2 are inner nodes
                tracker.push_less_of(1, node1)
                tracker.push_less_of(2, node2)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1.less, node2.less, tracker)
                tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1.less, node2.greater, tracker)
                tracker.pop()
                tracker.pop()
                    
                tracker.push_greater_of(1, node1)
                if node1 != node2:
                    # Avoid traversing (node1.less, node2.greater) and
                    # (node1.greater, node2.less) (it's the same node pair
                    # twice over, which is the source of the complication in
                    # the original KDTree.sparse_distance_matrix)
                    tracker.push_less_of(2, node2)
                    self.__sparse_distance_matrix_traverse(
                        other, results, node1.greater, node2.less, tracker)
                    tracker.pop()
                    
                tracker.push_greater_of(2, node2)
                self.__sparse_distance_matrix_traverse(
                    other, results, node1.greater, node2.greater, tracker)
                tracker.pop()
                tracker.pop()
                
        return 0
            
    def sparse_distance_matrix(cKDTree self, cKDTree other,
                               np.float64_t max_distance,
                               np.float64_t p=2.):
        """sparse_distance_matrix(self, other, max_distance, p=2.0)

        Compute a sparse distance matrix

        Computes a distance matrix between two KDTrees, leaving as zero
        any distance greater than max_distance.

        Parameters
        ----------
        other : cKDTree

        max_distance : positive float
        
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use. 

        Returns
        -------
        result : dok_matrix
            Sparse matrix representing the results in "dictionary of keys" format.
            FIXME: Internally, built as a COO matrix, it would be more
            efficient to return this COO matrix.

        """

        # Make sure trees are compatible
        if self.m != other.m:
            raise ValueError("Trees passed to sparse_distance_matrix have different dimensionality")

        # Calculate mins and maxes to outer box
        tracker = RectRectDistanceTracker(
            Rectangle(self.mins, self.maxes),
            Rectangle(other.mins, other.maxes),
            p, 0, max_distance)
        
        results = coo_entries()
        self.__sparse_distance_matrix_traverse(other, results,
                                               self.tree, other.tree,
                                               tracker)
        
        return results.to_matrix(shape=(self.n, other.n)).todok()
