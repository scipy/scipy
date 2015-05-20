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
#include "ckdtree_decl.h"
#include "query_methods.h"
#include "cpp_exc.h"
#include "rectangle.h"


static void
query_ball_point_traverse_no_checking(const ckdtree *self,
                                        std::vector<npy_intp> *results,
                                        const ckdtreenode *node)
{                                                    
    const ckdtreenode *lnode;
    npy_intp i;
    
    if (node->split_dim == -1) {  /* leaf node */
        lnode = node;
        for (i = lnode->start_idx; i < lnode->end_idx; ++i)
            results->push_back(self->raw_indices[i]);
    }
    else {
        query_ball_point_traverse_no_checking(self, results, node->less);
        query_ball_point_traverse_no_checking(self, results, node->greater);
    }
}


static void 
query_ball_point_traverse_checking(const ckdtree *self,
                                     std::vector<npy_intp> *results,
                                     const ckdtreenode *node,
                                     PointRectDistanceTracker *tracker)
{
    const ckdtreenode *lnode;
    npy_float64 d;
    npy_intp i;

    if (tracker->min_distance > tracker->upper_bound * tracker->epsfac)
        return;
    else if (tracker->max_distance < tracker->upper_bound / tracker->epsfac)
        query_ball_point_traverse_no_checking(self, results, node);
    else if (node->split_dim == -1)  { /* leaf node */
        
        /* brute-force */
                
        const npy_float64 *raw_data = self->raw_data;
        const npy_intp *raw_indices = self->raw_indices;
        const npy_intp m = self->m;
        
        lnode = node;
            
        prefetch_datapoint(raw_data+raw_indices[lnode->start_idx]*m, m);
        if (lnode->start_idx < lnode->end_idx)
            prefetch_datapoint(raw_data+raw_indices[lnode->start_idx+1]*m, m);
                
        for (i = lnode->start_idx; i < lnode->end_idx; ++i) {
            
            if (i < lnode->end_idx-2)
                prefetch_datapoint(raw_data+raw_indices[i+2]*m, m);
           
            d = _distance_p(
                raw_data + raw_indices[i] * m,
                tracker->pt, tracker->p, m, tracker->upper_bound);
                
            if (d <= tracker->upper_bound) {
                results->push_back((npy_intp) raw_indices[i]);
            }
        }
    }
    else {
    
        tracker->push_less_of(node);
        query_ball_point_traverse_checking(
            self, results, node->less, tracker);
        tracker->pop();
        
        tracker->push_greater_of(node);
        query_ball_point_traverse_checking(
            self, results, node->greater, tracker);
        tracker->pop();
    }    
}

        
extern "C" PyObject*
query_ball_point(const ckdtree *self,
                 const npy_float64 *x,
                 const npy_float64 r,
                 const npy_float64 p,
                 const npy_float64 eps,
                 const npy_intp n_queries,
                 std::vector<npy_intp> **results)
{
    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
            for (npy_intp i=0; i < n_queries; ++i) {
                const npy_intp m = self->m;
                Rectangle rect(m, self->raw_mins, self->raw_maxes);             
                PointRectDistanceTracker tracker(x + i*m, rect, p, eps, r);
                query_ball_point_traverse_checking(
                    self, results[i], self->ctree, &tracker);
            }
        } 
        catch(...) {
            translate_cpp_exception_with_gil();
        }
    }  
    /* reacquire the GIL */
    NPY_END_ALLOW_THREADS

    if (PyErr_Occurred()) 
        /* true if a C++ exception was translated */
        return NULL;
    else {
        /* return None if there were no errors */
        Py_RETURN_NONE;
    }
}        