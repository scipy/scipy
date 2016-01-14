#include <Python.h>
#include "numpy/arrayobject.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <vector>
#include <algorithm>
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

struct WeightedTree {
    const ckdtree *tree;
    npy_float64 *weights; 
    npy_float64 *node_weights; 
};
struct CNBParams
{
    npy_float64 *r; 
    void * results; /* will be casted inside */
    WeightedTree self, other;
};

template <typename MinMaxDist, typename WeightType, typename ResultType> static void
traverse(
    RectRectDistanceTracker<MinMaxDist> *tracker,
    const CNBParams *params, 
    npy_float64 *start, npy_float64 *end, 
    const ckdtreenode *node1,
    const ckdtreenode *node2)
{
    static void (* const next)(RectRectDistanceTracker<MinMaxDist> *tracker,
            const CNBParams *params, 
            npy_float64 *start, npy_float64 *end, 
            const ckdtreenode *node1,
            const ckdtreenode *node2) = traverse<MinMaxDist, WeightType, ResultType>;

    ResultType *results = (ResultType*) params->results;
    
    /* 
     * Speed through pairs of nodes all of whose children are close
     * and see if any work remains to be done
     */
    
    start = std::lower_bound(start, end, tracker->min_distance);
    end = std::lower_bound(start, end, tracker->max_distance);

    /* since max_distance >= min_distance, end < start never happens */
    if (end == start) {
        results[start - params->r] += WeightType::get_weight(&params->self, node1)
                                    * WeightType::get_weight(&params->other, node2);
        return;
    }

    /* OK, need to probe a bit deeper */
    if (node1->split_dim == -1) {  /* 1 is leaf node */
        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
            npy_intp i, j;
            const npy_float64 p = tracker->p;
            const npy_float64 tmd = tracker->max_distance;                
            const npy_float64 *sdata = params->self.tree->raw_data;
            const npy_intp *sindices = params->self.tree->raw_indices;
            const npy_float64 *odata = params->other.tree->raw_data;
            const npy_intp *oindices = params->other.tree->raw_indices;
            const npy_intp m = params->self.tree->m;
            const npy_intp start1 = node1->start_idx;
            const npy_intp start2 = node2->start_idx;
            const npy_intp end1 = node1->end_idx;
            const npy_intp end2 = node2->end_idx;
            
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
             
                    npy_float64 d = MinMaxDist::distance_p(params->self.tree,
                            sdata + sindices[i] * m,
                            odata + oindices[j] * m,
                            p, m, tmd);

                    const npy_float64 *l = std::lower_bound(start, end, d);
                    results[l - params->r] += WeightType::get_weight(&params->self, sindices[i])
                                            * WeightType::get_weight(&params->other, sindices[j]);   
                }
            }
        }
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            next(tracker, params, start, end, node1, node2->less);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            next(tracker, params, start, end, node1, node2->greater);
            tracker->pop();
        }
    }
    else { /* 1 is an inner node */
        if (node2->split_dim == -1) {
            /* 1 is an inner node, 2 is a leaf node */
            tracker->push_less_of(1, node1);
            next(tracker, params, start, end, node1->less, node2);
            tracker->pop();
            
            tracker->push_greater_of(1, node1);
            next(tracker, params, start, end, node1->greater, node2);
            tracker->pop();
        }
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            next(tracker, params, start, end, node1->less, node2->less);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            next(tracker, params, start, end, node1->less, node2->greater);
            tracker->pop();
            tracker->pop();
                
            tracker->push_greater_of(1, node1);
            tracker->push_less_of(2, node2);
            next(tracker, params, start, end, node1->greater, node2->less);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            next(tracker, params, start, end, node1->greater, node2->greater);
            tracker->pop();
            tracker->pop();
        }
    }
}

template <typename WeightType, typename ResultType> void
count_neighbors(struct CNBParams *params,
                npy_intp n_queries, const npy_float64 p)
{

    const ckdtree *self = params->self.tree;
    const ckdtree *other = params->other.tree;

#define HANDLE(cond, kls) \
    if (cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, 0.0, 0.0);\
        traverse<kls, WeightType, ResultType>(&tracker, params, params->r, params->r+n_queries, \
                 self->ctree, other->ctree); \
    } else

    Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
    Rectangle r2(other->m, other->raw_mins, other->raw_maxes);
    
    if (NPY_LIKELY(self->raw_boxsize_data == NULL)) {
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

struct Unweighted {
    /* the interface for accessing weights of unweighted data. */
    static inline npy_intp
    get_weight(const WeightedTree *wt, const ckdtreenode * node)
    {
        return node->children;
    }
    static inline npy_intp
    get_weight(const WeightedTree *wt, const npy_intp i)
    {
        return 1;
    }
};

extern "C" PyObject*
count_neighbors_unweighted(const ckdtree *self, const ckdtree *other,
                npy_intp n_queries, npy_float64 *real_r, npy_intp *results,
                const npy_float64 p) {

    CNBParams params = {0};

    params.r = real_r;
    params.results = (void*) results;
    params.self.tree = self;
    params.other.tree = other;

    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
            count_neighbors<Unweighted, npy_intp>(&params, n_queries, p);
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

struct Weighted {
    /* the interface for accessing weights of weighted data. */
    static inline npy_float64
    get_weight(const WeightedTree *wt, const ckdtreenode * node)
    {
        return (wt->weights != NULL)
           ? wt->node_weights[node - wt->tree->ctree]
           : node->children;
    }
    static inline npy_float64
    get_weight(const WeightedTree *wt, const npy_intp i)
    {
        return (wt->weights != NULL)?wt->weights[i]:1;
    }
};


extern "C" PyObject*
count_neighbors_weighted(const ckdtree *self, const ckdtree *other,
                npy_float64 *self_weights, npy_float64 *other_weights, 
                npy_float64 *self_node_weights, npy_float64 *other_node_weights, 
                npy_intp n_queries, npy_float64 *real_r, npy_float64 *results,
                const npy_float64 p) 
{

    CNBParams params = {0};

    params.r = real_r;
    params.results = (void*) results;

    if (self_weights) {
        params.self.tree = self;
        params.self.weights = self_weights;
        params.self.node_weights = self_node_weights;
    }
    if (other_weights) {
        params.other.tree = other;
        params.other.weights = other_weights;
        params.other.node_weights = other_node_weights;
    }
    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
            count_neighbors<Weighted, npy_float64>(&params, n_queries, p);
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

