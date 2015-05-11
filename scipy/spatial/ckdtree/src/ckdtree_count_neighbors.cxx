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

// TODO: port this over to C++

cdef int __count_neighbors_traverse(cKDTree self,
                                    cKDTree other,
                                    np.intp_t     n_queries,
                                    np.float64_t  *r,
                                    np.intp_t     *results,
                                    np.intp_t     *idx,
                                    ckdtreenode   *node1,
                                    ckdtreenode   *node2,
                                    RectRectDistanceTracker tracker):

        cdef ckdtreenode *lnode1
        cdef ckdtreenode *lnode2
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
                lnode1 = node1
                if node2.split_dim == -1:  # 1 & 2 are leaves
                    lnode2 = node2
                    
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
        


extern "C" PyObject*
count_neighbors(const ckdtree *self ...)
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