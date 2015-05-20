
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


struct coo_entries {

    std::vector<npy_intp> *pi;
    std::vector<npy_intp> *pj;
    std::vector<npy_float64> *pv;
        
    coo_entries(std::vector<npy_intp> *_pi, 
                std::vector<npy_intp> *_pj, 
                std::vector<npy_float64> *_pv) {                
        pi = _pi;
        pj = _pj;
        pv = _pv;
    };
        
    inline void
    add(npy_intp i, npy_intp j, npy_float64 v) {
        pi->push_back(i);
        pj->push_back(j);
        pv->push_back(v);
    };
};


static void
sparse_distance_matrix_traverse(const ckdtree *self, 
                                const ckdtree *other,
                                coo_entries *results,
                                const ckdtreenode *node1, 
                                const ckdtreenode *node2,
                                RectRectDistanceTracker *tracker)
{
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    
    npy_float64 d;
    npy_intp i, j, min_j;
            
    if (tracker->min_distance > tracker->upper_bound)
        return;
    else if (node1->split_dim == -1) {  /* 1 is leaf node */
        lnode1 = node1;
        
        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
            lnode2 = node2;
            
            /* brute-force */
            for (i = lnode1->start_idx; i < lnode1->end_idx; ++i) {
                /* Special care here to avoid duplicate pairs */
                if (node1 == node2)
                    min_j = i+1;
                else
                    min_j = lnode2->start_idx;
                    
                for (j = min_j; j < lnode2->end_idx; ++j) {
                    d = _distance_p(
                        self->raw_data + self->raw_indices[i] * self->m,
                        other->raw_data + other->raw_indices[j] * self->m,
                        tracker->p, self->m, tracker->upper_bound);
                    if (d <= tracker->upper_bound) {
                        if ((tracker->p != 1) && (tracker->p != infinity))
                            d = std::pow(d, 1. / tracker->p);
                        results->add(self->raw_indices[i],
                                     other->raw_indices[j], d);
                        if (node1 == node2)
                            results->add(self->raw_indices[j],
                                        other->raw_indices[i], d);
                    }
                }
            }
        }
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            sparse_distance_matrix_traverse(
                self, other, results, node1, node2->less, tracker);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            sparse_distance_matrix_traverse(
                self, other, results, node1, node2->greater, tracker);
            tracker->pop();
        }
    }        
    else {  /* 1 is an inner node */
        if (node2->split_dim == -1) {  
            /* 1 is an inner node, 2 is a leaf node*/
            tracker->push_less_of(1, node1);
            sparse_distance_matrix_traverse(
                self, other, results, node1->less, node2, tracker);
            tracker->pop();
            
            tracker->push_greater_of(1, node1);
            sparse_distance_matrix_traverse(
                self, other, results, node1->greater, node2, tracker);
            tracker->pop();
        }    
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            sparse_distance_matrix_traverse(
                self, other, results, node1->less, node2->less, tracker);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            sparse_distance_matrix_traverse(
                self, other, results, node1->less, node2->greater, tracker);
            tracker->pop();
            tracker->pop();
                
            tracker->push_greater_of(1, node1);
            if (node1 != node2) {
                /*
                 * Avoid traversing (node1->less, node2->greater) and
                 * (node1->greater, node2->less) (it's the same node pair
                 * twice over, which is the source of the complication in
                 * the original KDTree.sparse_distance_matrix)
                 */
                tracker->push_less_of(2, node2);
                sparse_distance_matrix_traverse(
                    self, other, results, node1->greater, node2->less, tracker);
                tracker->pop();
            }    
            tracker->push_greater_of(2, node2);
            sparse_distance_matrix_traverse(
                self, other, results, node1->greater, node2->greater, tracker);
            tracker->pop();
            tracker->pop();
        }    
    }
}



        
extern "C" PyObject*
sparse_distance_matrix(const ckdtree *self,
                       const ckdtree *other,
                       const npy_float64 p,
                       const npy_float64 max_distance,
                       std::vector<npy_intp> *results_i,
                       std::vector<npy_intp> *results_j,
                       std::vector<npy_float64> *results_v)
{

    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
        
            Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
            Rectangle r2(other->m, other->raw_mins, other->raw_maxes); 
            
            RectRectDistanceTracker tracker(r1, r2, p, 0, max_distance);
                    
            coo_entries results(results_i, results_j, results_v);
            
            sparse_distance_matrix_traverse(
                self, other, &results, self->ctree, other->ctree, &tracker);
                                               
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