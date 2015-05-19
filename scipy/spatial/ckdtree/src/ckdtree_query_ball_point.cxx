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

cdef int __query_ball_point_traverse_no_checking(cKDTree self,
                                                     list results,
                                                     ckdtreenode *node) except -1:
        cdef ckdtreenode *lnode
        cdef np.intp_t i

        if node.split_dim == -1:  # leaf node
            lnode = node
            for i in range(lnode.start_idx, lnode.end_idx):
                list_append(results, self.raw_indices[i])
        else:
            self.__query_ball_point_traverse_no_checking(results, node.less)
            self.__query_ball_point_traverse_no_checking(results, node.greater)

        return 0


    @cython.cdivision(True)
    cdef int __query_ball_point_traverse_checking(cKDTree self,
                                                  list results,
                                                  ckdtreenode *node,
                                                  PointRectDistanceTracker tracker) except -1:
        cdef ckdtreenode *lnode
        cdef np.float64_t d
        cdef np.intp_t i

        if tracker.min_distance > tracker.upper_bound * tracker.epsfac:
            return 0
        elif tracker.max_distance < tracker.upper_bound / tracker.epsfac:
            self.__query_ball_point_traverse_no_checking(results, node)
        elif node.split_dim == -1:  # leaf node
            lnode = <ckdtreenode*>node
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
                                 np.float64_t *x,
                                 np.float64_t r,
                                 np.float64_t p,
                                 np.float64_t eps):

        tracker = PointRectDistanceTracker()
        tracker.init(x, Rectangle(self.mins, self.maxes),
                     p, eps, r)
        
        results = []
        self.__query_ball_point_traverse_checking(
            results, self.ctree, tracker)
        return results
        
        
        
extern "C" PyObject*
query_ball_point(const ckdtree *self ...)
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
        