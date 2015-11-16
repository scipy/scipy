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

struct traverse_weights
{
    const ckdtreenode *self_root; /* to translate the pointer to node_index */
    npy_float64 *self_weights; 
    npy_float64 *self_node_weights; 
    const ckdtreenode *other_root; /* to translate the pointer to node_index */
    npy_float64 *other_weights; 
    npy_float64 *other_node_weights; 
};

struct Weighted {
    static inline npy_float64
    get_node_weight(const traverse_weights * w, 
                    const ckdtreenode * node1, 
                    const ckdtreenode * node2) {
        npy_float64 w1, w2;

        if(w->self_root != NULL)
            w1 = w->self_node_weights[node1 - w->self_root];
        else
            w1 = node1->children;
        
        if(w->other_root != NULL)
            w2 = w->other_node_weights[node2 - w->other_root];
        else
            w2 = node2->children;
        return w1 * w2;
    }
    static inline npy_float64
    get_weight(const traverse_weights * w, 
               const npy_intp i, 
               const npy_intp j) {
        npy_float64 w1, w2;

        if(w->self_root != NULL)
            w1 = w->self_weights[i];
        else
            w1 = 1;

        if(w->other_root != NULL)
            w2 = w->other_weights[j];
        else
            w2 = 1;
        
        return w1 * w2;
    }
};

struct Unweighted {
    static inline npy_intp
    get_node_weight(const traverse_weights * w, 
                    const ckdtreenode * node1, 
                    const ckdtreenode * node2) {
        return node1->children * node2->children;
    }
    static inline npy_intp
    get_weight(const traverse_weights * w, 
               const npy_intp i, 
               const npy_intp j) {
        return 1;
    }
};

template <typename MinMaxDist, typename WeightType, typename ResultType> static void
traverse(const ckdtree *self, const ckdtree *other,
         const traverse_weights *w,       
         npy_intp start, npy_intp end, npy_float64 *r, ResultType *results,
         const ckdtreenode *node1, const ckdtreenode *node2,
         RectRectDistanceTracker<MinMaxDist> *tracker)
{

    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    npy_float64 d;
    npy_intp l, i, j;
    
    /* 
     * Speed through pairs of nodes all of whose children are close
     * and see if any work remains to be done
     */
    
    for (i=end - 1; i >= start; --i) {
        if (tracker->max_distance <= r[i]) {
            results[i] += WeightType::get_node_weight(w, node1, node2);
        } else {
            break;
        }
    }
    end = i + 1;

    for (i=start; i<end; ++i) {
        if(tracker->min_distance <= r[i]) {
            break;
        }
    }
    start = i;
    
    if (end - start > 0) {
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
                        for (l=start; l<end; ++l) {
                            if (d <= r[l]) results[l] += WeightType::get_weight(w, sindices[i], sindices[j]);
                        }
                    }
                }
            }
            else {  /* 1 is a leaf node, 2 is inner node */
                tracker->push_less_of(2, node2);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results,
                    node1, node2->less, tracker);
                tracker->pop();

                tracker->push_greater_of(2, node2);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results,
                    node1, node2->greater, tracker);
                tracker->pop();
            }
        }
        else { /* 1 is an inner node */
            if (node2->split_dim == -1) {
                /* 1 is an inner node, 2 is a leaf node */
                tracker->push_less_of(1, node1);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results,
                    node1->less, node2, tracker);
                tracker->pop();
                
                tracker->push_greater_of(1, node1);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results, 
                    node1->greater, node2, tracker);
                tracker->pop();
            }
            else { /* 1 and 2 are inner nodes */
                tracker->push_less_of(1, node1);
                tracker->push_less_of(2, node2);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results, 
                    node1->less, node2->less, tracker);
                tracker->pop();
                    
                tracker->push_greater_of(2, node2);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results,
                    node1->less, node2->greater, tracker);
                tracker->pop();
                tracker->pop();
                    
                tracker->push_greater_of(1, node1);
                tracker->push_less_of(2, node2);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results, 
                    node1->greater, node2->less, tracker);
                tracker->pop();
                    
                tracker->push_greater_of(2, node2);
                traverse<MinMaxDist, WeightType, ResultType>(self, other, w, start, end, r, results,
                    node1->greater, node2->greater, tracker);
                tracker->pop();
                tracker->pop();
            }
        }
    }
}

template <typename WeightType, typename ResultType> PyObject*
count_neighbors(const ckdtree *self, const ckdtree *other,
                struct traverse_weights * w,
                npy_intp n_queries, npy_float64 *real_r, ResultType *results,
                const npy_float64 p)
{

#define HANDLE(cond, kls) \
    if(cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, 0.0, 0.0);\
        traverse<kls, WeightType, ResultType>(self, other, w, 0, n_queries, real_r, results, \
                 self->ctree, other->ctree, &tracker); \
    } else

    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
            
            Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
            Rectangle r2(other->m, other->raw_mins, other->raw_maxes);
            
            if(NPY_LIKELY(self->raw_boxsize_data == NULL)) {
                HANDLE(NPY_LIKELY(p == 2), MinkowskiDistP2)
                HANDLE(p == 1, MinkowskiDistP1)
                HANDLE(ckdtree_isinf(p), MinkowskiDistPinf)
                HANDLE(1, MinkowskiDistPp) 
                {}
            } else {
                HANDLE(NPY_LIKELY(p == 2), BoxMinkowskiDistP2)
                HANDLE(p == 1, BoxMinkowskiDistP1)
                HANDLE(ckdtree_isinf(p), BoxMinkowskiDistPinf)
                HANDLE(1, BoxMinkowskiDistPp) 
                {}
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

extern "C" PyObject*
count_neighbors_unweighted(const ckdtree *self, const ckdtree *other,
                npy_intp n_queries, npy_float64 *real_r, npy_intp *results,
                const npy_float64 p) {
    return count_neighbors<Unweighted, npy_intp>(self, other, NULL, n_queries, real_r, results, p);
}

extern "C" PyObject*
count_neighbors_weighted(const ckdtree *self, const ckdtree *other,
                npy_float64 *self_weights, npy_float64 *other_weights, 
                npy_float64 *self_node_weights, npy_float64 *other_node_weights, 
                npy_intp n_queries, npy_float64 *real_r, npy_float64 *results,
                const npy_float64 p) {

    struct traverse_weights w = {0};
    if(self_weights) {
        w.self_root = self->ctree;
        w.self_weights = self_weights;
        w.self_node_weights = self_node_weights;
    }
    if(other_weights) {
        w.other_root = other->ctree;
        w.other_weights = other_weights;
        w.other_node_weights = other_node_weights;
    }
    return count_neighbors<Weighted, npy_float64>(self, other, &w, n_queries, real_r, results, p);
}

