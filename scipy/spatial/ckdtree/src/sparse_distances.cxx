
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
#include "coo_entries.h"

template <typename MinMaxDist> static void
traverse(const ckdtree *self, const ckdtree *other, 
         std::vector<coo_entry> *results,
         const ckdtreenode *node1, const ckdtreenode *node2,
         RectRectDistanceTracker<MinMaxDist> *tracker)
{
            
    if (tracker->min_distance > tracker->upper_bound)
        return;
    else if (node1->split_dim == -1) {  /* 1 is leaf node */
    
        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */            
            /* brute-force */
            const npy_float64 p = tracker->p;
            const npy_float64 tub = tracker->upper_bound;
            const npy_float64 *sdata = self->raw_data;
            const npy_intp *sindices = self->raw_indices;
            const npy_float64 *odata = other->raw_data;
            const npy_intp *oindices = other->raw_indices;
            const npy_intp m = self->m;
            const npy_intp start1 = node1->start_idx;
            const npy_intp start2 = node2->start_idx;
            const npy_intp end1 = node1->end_idx;
            const npy_intp end2 = node2->end_idx;
                        
            prefetch_datapoint(sdata + sindices[start1] * m, m);
            if (start1 < end1)
               prefetch_datapoint(sdata + sindices[start1+1] * m, m);                         
                        
            for (npy_intp i = start1; i < end1; ++i) {
            
                if (i < end1-2)
                     prefetch_datapoint(sdata + sindices[i+2] * m, m);
            
                prefetch_datapoint(odata + oindices[start2] * m, m);
                if (start2 < end2)
                    prefetch_datapoint(sdata + oindices[start2+1] * m, m);
                    
                for (npy_intp j = start2; j < end2; ++j) {
                
                    if (j < end2-2)
                        prefetch_datapoint(odata + oindices[j+2] * m, m);
                
                    npy_float64 d = MinMaxDist::distance_p(
                            self,
                            sdata + sindices[i] * m,
                            odata + oindices[j] * m,
                            p, m, tub);
                        
                    if (d <= tub) {
                        if (NPY_LIKELY(p == 2.0))
                            d = std::sqrt(d);
                        else if ((p != 1) && (p != infinity))
                            d = std::pow(d, 1. / p);
                         
                        coo_entry e = {sindices[i], oindices[j], d};
                        results->push_back(e);
                    }
                }
            }
        }
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            traverse(self, other, results, node1, node2->less, tracker);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            traverse(self, other, results, node1, node2->greater, tracker);
            tracker->pop();
        }
    }        
    else {  /* 1 is an inner node */
        if (node2->split_dim == -1) {  
            /* 1 is an inner node, 2 is a leaf node*/
            tracker->push_less_of(1, node1);
            traverse(self, other, results, node1->less, node2, tracker);
            tracker->pop();
            
            tracker->push_greater_of(1, node1);
            traverse(self, other, results, node1->greater, node2, tracker);
            tracker->pop();
        }    
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse(self, other, results, node1->less, node2->less, tracker);
            tracker->pop();
                
            tracker->push_greater_of(2, node2);
            traverse(self, other, results, node1->less, node2->greater, tracker);
            tracker->pop();
            tracker->pop();
                
            tracker->push_greater_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse(self, other, results, node1->greater, node2->less, tracker);
            tracker->pop();
               
            tracker->push_greater_of(2, node2);
            traverse(self, other, results, node1->greater, node2->greater, 
                tracker);
            tracker->pop();
            tracker->pop();
        }    
    }
}

        
extern "C" PyObject*
sparse_distance_matrix(const ckdtree *self, const ckdtree *other,
                       const npy_float64 p,
                       const npy_float64 max_distance,
                       std::vector<coo_entry> *results)
{
#define HANDLE(cond, kls) \
    if(cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, 0, max_distance);\
        traverse(self, other, results, self->ctree, other->ctree, &tracker); \
    } else

    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS   
    {
        try {
        
            Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
            Rectangle r2(other->m, other->raw_mins, other->raw_maxes);             
            if(NPY_LIKELY(self->raw_boxsize_data == NULL)) {
                HANDLE(NPY_LIKELY(p == 2), MinMaxDistP2)
                HANDLE(p == 1, MinMaxDistP1)
                HANDLE(p == infinity, MinMaxDistPinf)
                HANDLE(1, MinMaxDistPp) 
                {}
            } else {
                HANDLE(NPY_LIKELY(p == 2), BoxMinMaxDistP2)
                HANDLE(p == 1, BoxMinMaxDistP1)
                HANDLE(p == infinity, BoxMinMaxDistPinf)
                HANDLE(1, BoxMinMaxDistPp) 
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
