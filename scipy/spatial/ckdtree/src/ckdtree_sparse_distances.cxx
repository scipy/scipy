
#include <Python.h>
#include "numpy/arrayobject.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <vector>
#include <string>
#include <sstream>
#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>

#define CKDTREE_METHODS_IMPL
#include "ckdtree_cpp_decl.h"
#include "ckdtree_cpp_methods.h"
#include "ckdtree_cpp_exc.h"
#include "ckdtree_cpp_rectangle"


// TODO: port me over to C++

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
        return scipy.sparse.coo_matrix((self.v, (self.i, self.j)),
                                       shape=shape)



cdef int __sparse_distance_matrix_traverse(cKDTree self, cKDTree other,
                                               coo_entries results,
                                               ckdtreenode *node1, ckdtreenode *node2,
                                               RectRectDistanceTracker tracker) except -1:
        cdef ckdtreenode *lnode1
        cdef ckdtreenode *lnode2
        cdef list results_i
        cdef np.float64_t d
        cdef np.intp_t i, j, min_j
                
        if tracker.min_distance > tracker.upper_bound:
            return 0
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = node2
                
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
        
        
extern "C" PyObject*
sparse_distance_matrix(const ckdtree *self ...)
{

    // release the GIL
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
             
             // TODO: Add C++ code to start the query
             
        } 
        catch(...) {
            translate_cpp_exception_with_gil();
        }
    }  
    // reacquire the GIL
    NPY_END_ALLOW_THREADS

    if (PyErr_Occurred()) 
        // true if a C++ exception was translated
        return NULL;
    else {
        // return None if there were no errors
        Py_RETURN_NONE;
    }
}