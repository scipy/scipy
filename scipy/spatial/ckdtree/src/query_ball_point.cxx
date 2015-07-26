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
#include "ckdtree_methods.h"
#include "cpp_exc.h"
#include "rectangle.h"


static void
traverse_no_checking(const ckdtree *self,
                     std::vector<npy_intp> *results,
                     const ckdtreenode *node)
{                                                    
    const npy_intp *indices = self->raw_indices;
    const ckdtreenode *lnode;
    npy_intp i;
    
    if (node->split_dim == -1) {  /* leaf node */
        lnode = node;
        const npy_intp start = lnode->start_idx;
        const npy_intp end = lnode->end_idx;
        for (i = start; i < end; ++i)
            results->push_back(indices[i]);
    }
    else {
        traverse_no_checking(self, results, node->less);
        traverse_no_checking(self, results, node->greater);
    }
}


static void 
traverse_checking(const ckdtree *self,
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
        traverse_no_checking(self, results, node);
    else if (node->split_dim == -1)  { /* leaf node */
        
        /* brute-force */
        lnode = node;
        const npy_float64 p = tracker->p;
        const npy_float64 tub = tracker->upper_bound;
        const npy_float64 *tpt = tracker->pt;        
        const npy_float64 *data = self->raw_data;
        const npy_intp *indices = self->raw_indices;
        const npy_intp m = self->m;
        const npy_intp start = lnode->start_idx;
        const npy_intp end = lnode->end_idx;
            
        prefetch_datapoint(data + indices[start] * m, m);
        if (start < end)
            prefetch_datapoint(data + indices[start+1] * m, m);
                
        for (i = start; i < end; ++i) {
            
            if (i < end-2)
                prefetch_datapoint(data + indices[i+2] * m, m);
           
            d = _distance_p(data + indices[i] * m, tpt, p, m, tub);
                
            if (d <= tub) {
                results->push_back((npy_intp) indices[i]);
            }
        }
    }
    else {
        tracker->push_less_of(node);
        traverse_checking(self, results, node->less, tracker);
        tracker->pop();
        
        tracker->push_greater_of(node);
        traverse_checking(self, results, node->greater, tracker);
        tracker->pop();
    }    
}

        
extern "C" PyObject*
query_ball_point(const ckdtree *self, const npy_float64 *x,
                 const npy_float64 r, const npy_float64 p, const npy_float64 eps,
                 const npy_intp n_queries, std::vector<npy_intp> **results)
{
    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
            for (npy_intp i=0; i < n_queries; ++i) {
                const npy_intp m = self->m;
                Rectangle rect(m, self->raw_mins, self->raw_maxes);             
                PointRectDistanceTracker tracker(x + i*m, rect, p, eps, r);
                traverse_checking(
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