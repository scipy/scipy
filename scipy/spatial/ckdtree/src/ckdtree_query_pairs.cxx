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
#include "ckdtree_cpp_ordered_pair.h"
#include "ckdtree_cpp_methods.h"
#include "ckdtree_cpp_exc.h"
#include "ckdtree_cpp_rectangle.h"


static void
__query_pairs_traverse_no_checking(const ckdtree *self,
                                   std::vector<ordered_pair> *results,
                                   const ckdtreenode *node1,
                                   const ckdtreenode *node2)
{                                            
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    npy_intp i, j, min_j;
    
    if (node1->split_dim == -1) { /* leaf node */
        lnode1 = node1;
        
        if (node2->split_dim == -1) { /* leaf node */
            lnode2 = node2;

            for (i = lnode1->start_idx; i < lnode1->end_idx; ++i) {
            
                /* Special care here to avoid duplicate pairs */
                if (node1 == node2)
                    min_j = i + 1;
                else
                    min_j = lnode2->start_idx;
                    
                for (j = min_j; j <  lnode2->end_idx; ++j)
                    add_ordered_pair(results,
                                         self->raw_indices[i],
                                           self->raw_indices[j]);
            }
        }                
        else {
            __query_pairs_traverse_no_checking(self, results, node1, 
                node2->less);
            __query_pairs_traverse_no_checking(self, results, node1, 
                node2->greater);
        }
    }
    else {
        if (node1 == node2) {
            /*
             * Avoid traversing (node1->less, node2->greater) and
             * (node1->greater, node2->less) (it's the same node pair twice
             * over, which is the source of the complication in the
             * original KDTree.query_pairs)
             */
            __query_pairs_traverse_no_checking(self, results, node1->less, 
                node2->less);
            __query_pairs_traverse_no_checking(self, results, node1->less, 
                node2->greater);
            __query_pairs_traverse_no_checking(self, results, node1->greater, 
                node2->greater);
        }
        else {
            __query_pairs_traverse_no_checking(self, results, node1->less, 
                node2);
            __query_pairs_traverse_no_checking(self, results, node1->greater, 
                node2);
        }
    }
}    


static void
__query_pairs_traverse_checking(const ckdtree *self,
                                std::vector<ordered_pair> *results,
                                const ckdtreenode *node1,
                                const ckdtreenode *node2,
                                RectRectDistanceTracker *tracker)
{
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    npy_float64 d;
    npy_intp i, j, min_j;
    
    if (tracker->min_distance > tracker->upper_bound * tracker->epsfac)
        return;
    else if (tracker->max_distance < tracker->upper_bound / tracker->epsfac)
        __query_pairs_traverse_no_checking(self, results, node1, node2);
    else if (node1->split_dim == -1) { /* 1 is leaf node */
        lnode1 = node1;
        
        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
            lnode2 = node2;
            
            /* brute-force */
            
            const npy_float64 *raw_data = self->raw_data;
            const npy_intp *raw_indices = self->raw_indices;
            const npy_intp m = self->m;
            
            prefetch_datapoint(raw_data+raw_indices[lnode1->start_idx]*m, m);
            if (lnode1->start_idx < lnode1->end_idx)
               prefetch_datapoint(raw_data+raw_indices[lnode1->start_idx+1]*m, m);
            
            for(i = lnode1->start_idx; i < lnode1->end_idx; ++i) {
            
                if (i < lnode1->end_idx-2)
                     prefetch_datapoint(raw_data+raw_indices[i+2]*m, m);
                               
                /* Special care here to avoid duplicate pairs */
                if (node1 == node2)
                    min_j = i + 1;
                else
                    min_j = lnode2->start_idx;
                               
                prefetch_datapoint(raw_data+raw_indices[min_j]*m, m);
                if (min_j < lnode2->end_idx)
                    prefetch_datapoint(raw_data+raw_indices[min_j+1]*m, m);
                            
                for (j = min_j; j < lnode2->end_idx; ++j) {
                                        
                    if (j < lnode2->end_idx-2)
                        prefetch_datapoint(raw_data+raw_indices[j+2]*m, m);
                                        
                    d = _distance_p(
                        raw_data + raw_indices[i] * m,
                        raw_data + raw_indices[j] * m,
                        tracker->p, m, tracker->upper_bound);
                
                    if (d <= tracker->upper_bound)
                        add_ordered_pair(results,
                                             raw_indices[i],
                                             raw_indices[j]);
                }
            }
        }                      
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            __query_pairs_traverse_checking(
                self, results, node1, node2->less, tracker);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            __query_pairs_traverse_checking(
                self, results, node1, node2->greater, tracker);
            tracker->pop();
        }
    }        
    else {  /* 1 is an inner node */
        if (node2->split_dim == -1) {  /* 1 is an inner node, 2 is a leaf node */
            tracker->push_less_of(1, node1);
            __query_pairs_traverse_checking(
                self, results, node1->less, node2, tracker);
            tracker->pop();
            
            tracker->push_greater_of(1, node1);
            __query_pairs_traverse_checking(
                self, results, node1->greater, node2, tracker);
            tracker->pop();
        }    
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            __query_pairs_traverse_checking(
                self, results, node1->less, node2->less, tracker);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            __query_pairs_traverse_checking(
                self, results, node1->less, node2->greater, tracker);
            tracker->pop();
            tracker->pop();
                
            tracker->push_greater_of(1, node1);
            if (node1 != node2) {
                /*
                 * Avoid traversing (node1->less, node2->greater) and
                 * (node1->greater, node2->less) (it's the same node pair
                 * twice over, which is the source of the complication in
                 * the original KDTree.query_pairs)
                 */
                tracker->push_less_of(2, node2);
                __query_pairs_traverse_checking(
                    self, results, node1->greater, node2->less, tracker);
                tracker->pop();
            }    
            tracker->push_greater_of(2, node2);
            __query_pairs_traverse_checking(
                self, results, node1->greater, node2->greater, tracker);
            tracker->pop();
            tracker->pop();
        }
    }
}


#include <iostream>

extern "C" PyObject*
query_pairs(const ckdtree *self, 
            const npy_float64 r, 
            const npy_float64 p, 
            const npy_float64 eps,
            std::vector<ordered_pair> *results)
{

    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {    
                        
            Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
            Rectangle r2(self->m, self->raw_mins, self->raw_maxes);
                                    
            RectRectDistanceTracker tracker(r1, r2, p, eps, r);
            
            __query_pairs_traverse_checking(
                self, results, self->ctree, self->ctree, &tracker);
             
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

