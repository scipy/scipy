/*
 * Nearest neighbor query in C++ for cKDTree
 * Written by Sturla Molden 2015
 * SciPy license
 */

/*
 * This would break SciPy with NumPy 1.6 so just accept the compiler
 * warning for now.
 * #define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION 
 *
 */
 
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

/*
 * Priority queue
 * ==============
 */
 
union heapcontents {
    npy_intp intdata;
    void     *ptrdata;
};

struct heapitem {
    npy_float64 priority;
    heapcontents contents;
}; 

struct heap {

    std::vector<heapitem> _heap;
    npy_intp n;
    npy_intp space;
    
    heap (npy_intp initial_size) : _heap(initial_size) {
        space = initial_size;
        n = 0;
    }
        
    inline void push(heapitem &item) {
        npy_intp i;
        heapitem t;
        n++;
        
        if (n > space) _heap.resize(2*space+1);
        space = _heap.size();
        
        i = n-1;
        _heap[i] = item;
        while ((i > 0) && (_heap[i].priority < _heap[(i-1)/2].priority)) {
            t = _heap[(i-1)/2];
            _heap[(i-1)/2] = _heap[i];
            _heap[i] = t;
            i = (i-1)/2;
        }
    }
    
    inline heapitem peek() {return _heap[0];}
        
    inline void remove() {
        heapitem t;
        npy_intp i, j, k, l, nn;    
        _heap[0] = _heap[n-1];
        n--;
        /*
         * No point in freeing up space as the heap empties.
         * The whole heap gets deallocated at the end of any 
         * query below. Just keep the space to avoid unneccesary
         * reallocs.
         */
        nn = n;        
        i=0;
        j=1; 
        k=2;
        while (((j<nn) && (_heap[i].priority > _heap[j].priority)) ||
               ((k<nn) && (_heap[i].priority > _heap[k].priority))) {
            if ((k<nn) && (_heap[j].priority >_heap[k].priority))
                l = k;
            else
                l = j;
            t = _heap[l];
            _heap[l] = _heap[i];
            _heap[i] = t;
            i = l;
            j = 2*i+1;
            k = 2*i+2;
        }
    }
    
    inline heapitem pop() {
        heapitem it = _heap[0];
        remove();
        return it;
    }
};


/*
 * nodeinfo
 * ========
 */
 
struct nodeinfo {
    nodeinfo           *next;
    nodeinfo           *prev;
    const ckdtreenode  *node;
    npy_float64        *mins; // the good old struct hack       
    npy_float64        *maxes; // the good old struct hack       
    npy_float64        *side_distances; // the good old struct hack       
    npy_float64        buf[1]; // the good old struct hack       
};        

/*
 * Memory pool for nodeinfo structs
 * ================================
 */
 
struct nodeinfo_pool {

    std::vector<char*> pool;
   
    npy_intp alloc_size;
    npy_intp arena_size;
    npy_intp m;
    char *arena;
    char *arena_ptr;
    
    nodeinfo_pool(npy_intp m) {
        alloc_size = sizeof(nodeinfo) + (3 * m -1)*sizeof(npy_float64);
        alloc_size = 64*(alloc_size/64)+64;
        arena_size = 4096*((64*alloc_size)/4096)+4096;
        arena = new char[arena_size];
        arena_ptr = arena;
        pool.push_back(arena);
        this->m = m;
    }
    
    ~nodeinfo_pool() {
        for (npy_intp i = pool.size()-1; i >= 0; --i)
            delete [] pool[i];
    }
    
    inline nodeinfo *allocate() {
        nodeinfo *ni;
        npy_uintp m1 = (npy_uintp)arena_ptr;
        npy_uintp m0 = (npy_uintp)arena;
        if ((arena_size-(npy_intp)(m1-m0))<alloc_size) {
            arena = new char[arena_size];
            arena_ptr = arena;
            pool.push_back(arena);
        }
        ni = (nodeinfo*)arena_ptr;
        ni->maxes = &ni->buf[0];
        ni->mins = &ni->buf[this->m];
        ni->side_distances = &ni->buf[2 * this->m];
        arena_ptr += alloc_size;
        return ni;
    }
    inline nodeinfo *dup(nodeinfo * ni) {
        nodeinfo * ni_new = this->allocate();
        memcpy(ni_new, ni, this->alloc_size);
        ni_new->maxes = &ni_new->buf[0];
        ni_new->mins = &ni_new->buf[this->m];
        ni_new->side_distances = &ni_new->buf[2 * this->m];
        return ni_new;
    }
};

static inline npy_float64 side_distance_from_min_max(
    const npy_float64 x,
    const npy_float64 min,
    const npy_float64 max,
    const npy_float64 p,
    const npy_float64 infinity) {
    npy_float64 s, t;

    s = 0; 
    t = x - max;
    if (t > s) {
        s = t;
    } else {
        t = min - x;
        if (t > s) s = t;
    }

#if 0
    printf("sdfmm: s %g t %g x %g min %g max %g\n", s, t, x, min, max);
#endif 

    if (NPY_UNLIKELY(p == 1 || p == infinity)) {
        s = dabs(s);
    } else if (NPY_LIKELY(p == 2)) {
        s = s * s;
    } else {
        s = std::pow(s,p);
    }
        
    return s;
}

// k-nearest neighbor search for a single point x
static void 
__query_single_point(const ckdtree *self, 
                     npy_float64   *result_distances, 
                     npy_intp      *result_indices, 
                     const npy_float64  *x, 
                     const npy_intp     k, 
                     const npy_float64  eps, 
                     const npy_float64  p, 
                     npy_float64  distance_upper_bound,
                     const npy_float64 infinity)
{                
    // memory pool to allocate and automatically reclaim nodeinfo structs
    nodeinfo_pool nipool(self->m);
    
    /*
     * priority queue for chasing nodes
     * entries are:
     *  - minimum distance between the cell and the target
     *  - distances between the nearest side of the cell and the target
     *    the head node of the cell 
     */
    heap q(12);
    
    /*
     *  priority queue for chasing nodes
     *  entries are:
     *   - minimum distance between the cell and the target
     *   - distances between the nearest side of the cell and the target
     *     the head node of the cell
     */
    heap neighbors(k);
    
    npy_intp      i, m = self->m;
    nodeinfo      *inf;
    nodeinfo      *inf1, *inf2;
    npy_float64   d;
    npy_float64   epsfac;
    npy_float64   min_distance;
    npy_float64   inf1_min_distance;
    npy_float64   inf2_min_distance;
    heapitem      it, it2, neighbor;
    const ckdtreenode   *node;
    const ckdtreenode   *inode;
    
    // set up first nodeinfo
    inf = nipool.allocate();
    inf->node = self->ctree;

    // buffer for temparory nodeinfo
    inf1 = nipool.allocate();
    inf2 = nipool.allocate();
 
    for (i=0; i<m; ++i) {
        inf->mins[i] = self->raw_mins[i];
        inf->maxes[i] = self->raw_maxes[i];

        inf->side_distances[i] = side_distance_from_min_max(
            x[i], inf->mins[i], inf->maxes[i], p, infinity);
    }
    
    // compute first distance
    min_distance = 0.;
    for (i=0; i<m; ++i) {
        if (p == infinity)
            min_distance = dmax(min_distance,inf->side_distances[i]);
        else
            min_distance += inf->side_distances[i];
    }
    
    // fiddle approximation factor
    if (eps == 0.)
        epsfac = 1.;
    else if (p == infinity)
        epsfac = 1./(1+eps);
    else
        epsfac = 1./std::pow((1+eps),p);

    // internally we represent all distances as distance**p
    if ((p != infinity) && (distance_upper_bound != infinity))
        distance_upper_bound = std::pow(distance_upper_bound,p);

    for(;;) {
        if (inf->node->split_dim == -1) {
        
            node = inf->node;

            // brute-force
            {
                npy_intp start_idx = node->start_idx;
                npy_intp end_idx = node->end_idx;
                const npy_float64 *raw_data = self->raw_data;
                const npy_intp *raw_indices = self->raw_indices;
                
                prefetch_datapoint(raw_data+raw_indices[start_idx]*m, m);
                if (start_idx < end_idx)
                    prefetch_datapoint(raw_data+raw_indices[start_idx+1]*m, m);
                
                for (i=start_idx; i<end_idx; ++i) {
                       
                    if (i < end_idx-2)
                        prefetch_datapoint(raw_data+raw_indices[i+2]*m, m);
                
                    d = _distance_p(raw_data+raw_indices[i]*m, x, p, m, 
                            distance_upper_bound);
                        
                    if (d < distance_upper_bound) {
                        // replace furthest neighbor
                        if (neighbors.n == k)
                              neighbors.remove();
                        neighbor.priority = -d;
                        neighbor.contents.intdata = self->raw_indices[i];
                        neighbors.push(neighbor);
    
                        // adjust upper bound for efficiency
                        if (neighbors.n == k)
                            distance_upper_bound = -neighbors.peek().priority;
                    }
                }
            }
            // done with this node, get another                
            if (q.n == 0) {
                // no more nodes to visit
                break;
            } 
            else {
                it = q.pop();
                inf = (nodeinfo*)(it.contents.ptrdata);
                min_distance = it.priority;
            }
            
        } 
        else {
            inode = inf->node;

            /*
             * we don't push cells that are too far onto the queue at all,
             * but since the distance_upper_bound decreases, we might get 
             * here even if the cell's too far
             */
            if (min_distance > distance_upper_bound*epsfac) {
                // since this is the nearest cell, we're done, bail out 
                break;
            }
            
            for(i = 0; i < self->m; i ++) {
                inf1->mins[i] = inf->mins[i];
                inf1->maxes[i] = inf->maxes[i];
                inf1->side_distances[i] = inf->side_distances[i];
                inf2->mins[i] = inf->mins[i];
                inf2->maxes[i] = inf->maxes[i];
                inf2->side_distances[i] = inf->side_distances[i];
            }
            
            // set up children for searching
            // inf2 will be pushed to the queue
            inf1->maxes[inode->split_dim] = inode->split;
            inf1->node = inode->less;
            inf2->mins[inode->split_dim] = inode->split;
            inf2->node = inode->greater;
            
            inf1->side_distances[inode->split_dim] = 
                side_distance_from_min_max(
                    x[inode->split_dim],
                    inf1->mins[inode->split_dim],
                    inf1->maxes[inode->split_dim],
                    p, infinity);

            inf2->side_distances[inode->split_dim] = 
                side_distance_from_min_max(
                    x[inode->split_dim],
                    inf2->mins[inode->split_dim],
                    inf2->maxes[inode->split_dim],
                    p, infinity);

            /*
             * one side distance changes
             * we can adjust the minimum distance without recomputing
             */
            if (NPY_UNLIKELY(p == infinity)) {
                /*
                 * we never use side_distances in the l_infinity case
                 * so skip filling it in
                 * inf2->side_distances[inode->split_dim] = dabs(inode->split-x[inode->split_dim])
                 */
                inf1_min_distance = dmax(min_distance, dabs(inf1->side_distances[inode->split_dim]));
                inf2_min_distance = dmax(min_distance, dabs(inf2->side_distances[inode->split_dim]));
            } else {
                inf1_min_distance = min_distance -
                        inf->side_distances[inode->split_dim];
                        + inf1->side_distances[inode->split_dim];
                inf2_min_distance = min_distance -
                        inf->side_distances[inode->split_dim];
                        + inf2->side_distances[inode->split_dim];
            }

#define swap_ni(a, b) {nodeinfo * tmp; tmp = a; a = b; b = tmp;}

            /* Ensure inf1 is closer than inf2 */
            if (inf1_min_distance > inf2_min_distance) {
                npy_float64 tmp = inf2_min_distance;
                inf2_min_distance = inf1_min_distance;
                inf1_min_distance = tmp;
                swap_ni(inf1, inf2);
            }
            
#if 0
            printf("sd inf1: %g inf2: %g inf: %g\n", inf1->side_distances[inode->split_dim], 
                                                     inf2->side_distances[inode->split_dim],
                                                     inf->side_distances[inode->split_dim]
                        );
            printf("inf1: %g inf2: %g inf: %g\n", inf1_min_distance, inf2_min_distance, min_distance);
#endif
            /*
             * near child is at the same or closer than the distance as the current node
             * we're going here next, so no point pushing it on the queue
             * no need to recompute the distance or the side_distances
             */
            swap_ni(inf, inf1);
            min_distance = inf1_min_distance;

            /*
             * far child can be further
             * push it on the queue if it's near enough
             */

            if (inf2_min_distance<=distance_upper_bound*epsfac) {
                it2.priority = inf2_min_distance;
                it2.contents.ptrdata = (void*) nipool.dup(inf2);
                q.push(it2);
            }
        }
    }
    // fill output arrays with sorted neighbors 
    for (i=neighbors.n-1; i>=0; --i) {
        neighbor = neighbors.pop();
        result_indices[i] = neighbor.contents.intdata;
        if ((p==1.) || (p==infinity))
            result_distances[i] = -neighbor.priority;
        else
            result_distances[i] = std::pow((-neighbor.priority),(1./p));
    }

}

// Query n points for their k nearest neighbors

extern "C" PyObject*
query_knn(const ckdtree      *self, 
          npy_float64        *dd, 
          npy_intp           *ii, 
          const npy_float64  *xx,
          const npy_intp     n, 
          const npy_intp     k, 
          const npy_float64  eps, 
          const npy_float64  p, 
          const npy_float64  distance_upper_bound)
{
    npy_intp m = self->m;
    npy_intp i;
    
    // release the GIL
    NPY_BEGIN_ALLOW_THREADS
    {
        try {
            for (i=0; i<n; ++i) {
                npy_float64 *dd_row = dd + (i*k);
                npy_intp *ii_row = ii + (i*k);
                const npy_float64 *xx_row = xx + (i*m);                
                __query_single_point(self, dd_row, ii_row, xx_row, k, eps, p, distance_upper_bound, ::infinity);
            }    
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

