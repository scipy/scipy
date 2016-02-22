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
#include "ordered_pair.h"
#include "ckdtree_methods.h"
#include "cpp_exc.h"
#include "cpp_utils.h"
#include "partial_sort.h"



static npy_intp 
build(ckdtree *self, npy_intp start_idx, npy_intp end_idx,
      npy_float64 *maxes, npy_float64 *mins, 
      const int _median, const int _compact)
{
      
    const npy_intp m = self->m;
    const npy_float64 *data = self->raw_data;
    npy_intp *indices = (npy_intp *)(self->raw_indices);
    
    ckdtreenode new_node, *n, *root;
    npy_intp node_index, _less, _greater;    
    npy_intp i, j, p, q, d;
    npy_float64 size, split, minval, maxval;
    
    /* put a new node into the node stack */
    self->tree_buffer->push_back(new_node);
    node_index = self->tree_buffer->size() - 1;
    root = tree_buffer_root(self->tree_buffer);
    n = root + node_index;
    memset(n, 0, sizeof(n[0]));

    n->start_idx = start_idx;
    n->end_idx = end_idx;
    n->children = end_idx - start_idx;
    
    if (end_idx-start_idx <= self->leafsize) {
        /* below brute force limit, return leafnode */
        n->split_dim = -1;
        return node_index;
    }
    else {
    
        if (NPY_LIKELY(_compact)) {
            /* Recompute hyperrectangle bounds. This should lead to a more 
             * compact kd-tree but comes at the expense of larger construction
             * time. However, construction time is usually dwarfed by the
             * query time by orders of magnitude.
             */
            const npy_float64 *tmp_data_point;
            tmp_data_point = data + indices[start_idx] * m;
            for(i=0; i<m; ++i) {
                maxes[i] = tmp_data_point[i];
                mins[i] = tmp_data_point[i];
            }
            for (j = start_idx + 1; j < end_idx; ++j) {
                tmp_data_point = data + indices[j] * m;
                for(i=0; i<m; ++i) {
                    npy_float64 tmp = tmp_data_point[i];
                    maxes[i] = maxes[i] > tmp ? maxes[i] : tmp;
                    mins[i] = mins[i] < tmp ? mins[i] : tmp;
                }
            }
        }
    
        /* split on the dimension with largest spread */ 
        d = 0; 
        size = 0;
        for (i=0; i<m; ++i) {
            if (maxes[i] - mins[i] > size) {
                d = i;
                size = maxes[i] - mins[i];
            }
        }
        maxval = maxes[d];
        minval = mins[d];
        if (maxval == minval) {
            /* all points are identical; warn user?
             * return leafnode
             */
            n->split_dim = -1;
            return node_index;
        }   
        
        /* construct new inner node */
        
        if (NPY_LIKELY(_median)) {
            /* split on median to create a balanced tree
             * adopted from scikit-learn
             */
            i = (end_idx - start_idx) / 2;
            partition_node_indices(data, indices + start_idx, d, i, m,
                end_idx - start_idx);               
            p = start_idx + i;
            split = data[indices[p]*m+d];
        }
        else {
            /* split with the sliding midpoint rule */
            split = (maxval + minval) / 2;
        }

        p = start_idx;
        q = end_idx - 1;
        while (p <= q) {
            if (data[indices[p] * m + d] < split)
                ++p;
            else if (data[indices[q] * m + d] >= split)
                --q;
            else {
                npy_intp t = indices[p];
                indices[p] = indices[q];
                indices[q] = t;
                ++p;
                --q;
            }
        }
        /* slide midpoint if necessary */
        if (p == start_idx) {
            /* no points less than split */
            j = start_idx;
            split = data[indices[j] * m + d];
            for (i = start_idx+1; i < end_idx; ++i) {
                if (data[indices[i] * m + d] < split) {
                    j = i;
                    split = data[indices[j] * m + d];
                }
            }
            npy_intp t = indices[start_idx];
            indices[start_idx] = indices[j];
            indices[j] = t;
            p = start_idx + 1;
            q = start_idx;
        }
        else if (p == end_idx) {
            /* no points greater than split */
            j = end_idx - 1;
            split = data[indices[j] * m + d];
            for (i = start_idx; i < end_idx-1; ++i) {
                if (data[indices[i] * m + d] > split) {
                    j = i;
                    split = data[indices[j] * m + d];
                }
            }
            npy_intp t = indices[end_idx-1];
            indices[end_idx-1] = indices[j];
            indices[j] = t;
            p = end_idx - 1;
            q = end_idx - 2;
        }
        
        if (NPY_LIKELY(_compact)) {
            _less = build(self, start_idx, p, maxes, mins, _median, _compact);
            _greater = build(self, p, end_idx, maxes, mins, _median, _compact);    
        }
        else
        {
            std::vector<npy_float64> tmp(m);
            npy_float64 *mids = &tmp[0];
            
            for (i=0; i<m; ++i) mids[i] = maxes[i];
            mids[d] = split;
            _less = build(self, start_idx, p, mids, mins, _median, _compact);
            
            for (i=0; i<m; ++i) mids[i] = mins[i];
            mids[d] = split;
            _greater = build(self, p, end_idx, maxes, mids, _median, _compact);
        }   
        
        /* recompute n because std::vector can
         * reallocate its internal buffer
         */
        root = tree_buffer_root(self->tree_buffer); 
        n = root + node_index;
        /* fill in entries */
        n->_less = _less; 
        n->_greater = _greater;
        n->less = root + _less;
        n->greater = root + _greater;
        n->split_dim = d;
        n->split = split;
        
        return node_index;
    }
}
        
        

extern "C" PyObject*
build_ckdtree(ckdtree *self, npy_intp start_idx, npy_intp end_idx,
              npy_float64 *maxes, npy_float64 *mins, int _median, int _compact)
                       
{
    
    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS
    {
        try {
            build(self, start_idx, end_idx, maxes, mins, _median, _compact);        
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

static npy_float64
add_weights(ckdtree *self, 
           npy_float64 *node_weights, 
           npy_intp node_index, 
           npy_float64 *weights)
{
      
    npy_intp *indices = (npy_intp *)(self->raw_indices);

    ckdtreenode *n, *root;

    root = tree_buffer_root(self->tree_buffer); 

    n = root + node_index;
    
    npy_float64 sum = 0;

    if (n->split_dim != -1) {
        /* internal nodes; recursively calculate the total weight */
        npy_float64 left, right;
        left = add_weights(self, node_weights, n->_less, weights);
        right = add_weights(self, node_weights, n->_greater, weights);
        sum = left + right;
    } else {
        npy_intp i;

        /* Leaf nodes */
        for (i = n->start_idx; i < n->end_idx; ++i) {
            sum += weights[indices[i]];
        }
    }

    node_weights[node_index] = sum;
    return sum;
}
        
        
extern "C" PyObject*
build_weights (ckdtree *self, npy_float64 *node_weights, npy_float64 *weights)
{
    
    /* release the GIL */
    NPY_BEGIN_ALLOW_THREADS
    {
        try {
            add_weights(self, node_weights, 0, weights);        
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

