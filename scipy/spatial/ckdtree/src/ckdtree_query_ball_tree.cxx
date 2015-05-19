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


// TODO: port this to C++

cdef int __query_ball_tree_traverse_no_checking(cKDTree self,
                                                    cKDTree other,
                                                    list results,
                                                    ckdtreenode *node1,
                                                    ckdtreenode *node2) except -1:
        cdef ckdtreenode *lnode1
        cdef ckdtreenode *lnode2
        cdef list results_i
        cdef np.intp_t i, j
        
        if node1.split_dim == -1:  # leaf node
            lnode1 = node1
            
            if node2.split_dim == -1:  # leaf node
                lnode2 = node2
                
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
                                                 ckdtreenode *node1,
                                                 ckdtreenode *node2,
                                                 RectRectDistanceTracker tracker) except -1:
        cdef ckdtreenode *lnode1
        cdef ckdtreenode *lnode2
        cdef list results_i
        cdef np.float64_t d
        cdef np.intp_t i, j

        if tracker.min_distance > tracker.upper_bound * tracker.epsfac:
            return 0
        elif tracker.max_distance < tracker.upper_bound / tracker.epsfac:
            self.__query_ball_tree_traverse_no_checking(other, results, node1, node2)
        elif node1.split_dim == -1:  # 1 is leaf node
            lnode1 = node1
            
            if node2.split_dim == -1:  # 1 & 2 are leaves
                lnode2 = node2
                
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
        
        
extern "C" PyObject*
query_ball_tree(const ckdtree *self ...)
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
        