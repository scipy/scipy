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

template <typename MinMaxDist> static void
traverse(const ckdtree *self, const ckdtree *other,
         npy_intp n_queries, npy_float64 *r, npy_intp *results, npy_intp *idx,
         const ckdtreenode *node1, const ckdtreenode *node2,
         RectRectDistanceTracker<MinMaxDist> *tracker)
{

    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    npy_float64 d;
    npy_intp *old_idx;
    npy_intp old_n_queries, l, i, j;
    
    /* 
     * Speed through pairs of nodes all of whose children are close
     * and see if any work remains to be done
     */
    
    old_idx = idx;
    
    std::vector<npy_intp> inner_idx(n_queries);
    idx = &inner_idx[0];

    old_n_queries = n_queries;
    n_queries = 0;
    
    for (i=0; i<old_n_queries; ++i) {
        if (tracker->max_distance < r[old_idx[i]])
            results[old_idx[i]] += node1->children * node2->children;
        else if (tracker->min_distance <= r[old_idx[i]]) {
            idx[n_queries] = old_idx[i];
            ++n_queries;
        }
    }
    
    if (n_queries > 0) {
        /* OK, need to probe a bit deeper */
        if (node1->split_dim == -1) {  /* 1 is leaf node */
            lnode1 = node1;
            if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
                lnode2 = node2;
                const npy_float64 p = tracker->p;
                const npy_float64 tmd = tracker->max_distance;                
                const npy_float64 *sdata = self->raw_data;
                const npy_intp *sindices = self->raw_indices;
                const npy_float64 *odata = other->raw_data;
                const npy_intp *oindices = other->raw_indices;
                const npy_intp m = self->m;
                const npy_intp start1 = lnode1->start_idx;
                const npy_intp start2 = lnode2->start_idx;
                const npy_intp end1 = lnode1->end_idx;
                const npy_intp end2 = lnode2->end_idx;
                
                prefetch_datapoint(sdata + sindices[start1] * m, m);
                
                if (start1 < end1)
                    prefetch_datapoint(sdata + sindices[start1+1] * m, m);
                                        
                /* brute-force */
                for (i = start1; i < end1; ++i) {
                    
                    if (i < end1-2)
                        prefetch_datapoint(sdata + sindices[i+2] * m, m);
                                      
                    prefetch_datapoint(odata + oindices[start2] * m, m);
                        
                    if (start2 < end2)
                        prefetch_datapoint(odata + oindices[start2+1] * m, m);
                  
                    for (j = start2; j < end2; ++j) {
                     
                        if (j < end2-2)
                            prefetch_datapoint(odata + oindices[j+2] * m, m);
                 
                        d = MinMaxDist::distance_p(self,
                                sdata + sindices[i] * m,
                                odata + oindices[j] * m,
                                p, m, tmd);
                        /*
                         * I think it's usually cheaper to test d against all 
                         * r's than to generate a distance array, sort it, then
                         * search for all r's via binary search
                         */
                        for (l=0; l<n_queries; ++l) {
                            if (d <= r[idx[l]]) results[idx[l]] += 1;
                        }
                    }
                }
            }
            else {  /* 1 is a leaf node, 2 is inner node */
                tracker->push_less_of(2, node2);
                traverse(self, other, n_queries, r, results, idx,
                    node1, node2->less, tracker);
                tracker->pop();

                tracker->push_greater_of(2, node2);
                traverse(self, other, n_queries, r, results, idx,
                    node1, node2->greater, tracker);
                tracker->pop();
            }
        }
        else { /* 1 is an inner node */
            if (node2->split_dim == -1) {
                /* 1 is an inner node, 2 is a leaf node */
                tracker->push_less_of(1, node1);
                traverse(self, other, n_queries, r, results, idx,
                    node1->less, node2, tracker);
                tracker->pop();
                
                tracker->push_greater_of(1, node1);
                traverse(self, other, n_queries, r, results, idx,
                    node1->greater, node2, tracker);
                tracker->pop();
            }
            else { /* 1 and 2 are inner nodes */
                tracker->push_less_of(1, node1);
                tracker->push_less_of(2, node2);
                traverse(self, other, n_queries, r, results, idx,
                    node1->less, node2->less, tracker);
                tracker->pop();
                    
                tracker->push_greater_of(2, node2);
                traverse(self, other, n_queries, r, results, idx,
                    node1->less, node2->greater, tracker);
                tracker->pop();
                tracker->pop();
                    
                tracker->push_greater_of(1, node1);
                tracker->push_less_of(2, node2);
                traverse(self, other, n_queries, r, results, idx,
                    node1->greater, node2->less, tracker);
                tracker->pop();
                    
                tracker->push_greater_of(2, node2);
                traverse(self, other, n_queries, r, results, idx,
                    node1->greater, node2->greater, tracker);
                tracker->pop();
                tracker->pop();
            }
        }
    }
}


extern "C" PyObject*
count_neighbors(const ckdtree *self, const ckdtree *other,
                npy_intp n_queries, npy_float64 *real_r, npy_intp *results,
                npy_intp *idx, const npy_float64 p)
{

    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
            
            Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
            Rectangle r2(other->m, other->raw_mins, other->raw_maxes);
            
            RectRectDistanceTracker<MinMaxDist> tracker(self, r1, r2, p, 0.0, 0.0);
            
            traverse(self, other, n_queries, real_r, results, idx,
                     self->ctree, other->ctree, &tracker);
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


