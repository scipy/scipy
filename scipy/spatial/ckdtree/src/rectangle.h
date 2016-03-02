
#ifndef CKDTREE_CPP_RECTANGLE
#define CKDTREE_CPP_RECTANGLE

#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>
#include <cmath>
#include <cstring>


#ifndef NPY_UNLIKELY
#define NPY_UNLIKELY(x) (x)
#endif

#ifndef NPY_LIKELY
#define NPY_LIKELY(x) (x)
#endif


/* Interval arithmetic
 * ===================
 */
 
struct Rectangle {
    
    npy_intp m;
    npy_float64 *mins;
    npy_float64 *maxes;
    
    std::vector<npy_float64> mins_arr;
    std::vector<npy_float64> maxes_arr;

    Rectangle(const npy_intp _m, 
              const npy_float64 *_mins, 
              const npy_float64 *_maxes) : mins_arr(_m), maxes_arr(_m) {

        /* copy array data */
        m = _m;
        mins = &mins_arr[0];
        maxes = &maxes_arr[0];        
        std::memcpy((void*)mins, (void*)_mins, m*sizeof(npy_float64));
        std::memcpy((void*)maxes, (void*)_maxes, m*sizeof(npy_float64));
    };    
         
    Rectangle(const Rectangle& rect) : mins_arr(rect.m), maxes_arr(rect.m) {
        m = rect.m;
        mins = &mins_arr[0];
        maxes = &maxes_arr[0];        
        std::memcpy((void*)mins, (void*)rect.mins, m*sizeof(npy_float64));
        std::memcpy((void*)maxes, (void*)rect.maxes, m*sizeof(npy_float64));         
    };    
    
    Rectangle() : mins_arr(0), maxes_arr(0) {
        m = 0;
        mins = NULL;
        maxes = NULL;
    };
    
};

#include "ckdtree_methods.h"
#include "distance.h"
#include "distance_box.h"

/*
 * Rectangle-to-rectangle distance tracker
 * =======================================
 *
 * The logical unit that repeats over and over is to keep track of the
 * maximum and minimum distances between points in two hyperrectangles
 * as these rectangles are successively split.
 *
 * Example
 * -------
 * node1 encloses points in rect1, node2 encloses those in rect2
 *
 * cdef RectRectDistanceTracker dist_tracker
 * dist_tracker = RectRectDistanceTracker(rect1, rect2, p)
 *
 * ...
 *
 * if dist_tracker.min_distance < ...:
 *     ...
 *
 * dist_tracker.push_less_of(1, node1)
 * do_something(node1.less, dist_tracker)
 * dist_tracker.pop()
 *
 * dist_tracker.push_greater_of(1, node1)
 * do_something(node1.greater, dist_tracker)
 * dist_tracker.pop()
 *
 * Notice that Point is just a reduced case of Rectangle where
 * mins == maxes. 
 *
 */
 
struct RR_stack_item {
    npy_intp    which;
    npy_intp    split_dim;
    npy_float64 min_along_dim;
    npy_float64 max_along_dim;
    npy_float64 min_distance;
    npy_float64 max_distance;
};    

const npy_intp LESS = 1;
const npy_intp GREATER = 2;

template<typename MinMaxDist> 
    struct RectRectDistanceTracker {
    
    const ckdtree * tree;
    Rectangle rect1; 
    Rectangle rect2;
    npy_float64 p; 
    npy_float64 epsfac;
    npy_float64 upper_bound;
    npy_float64 min_distance;
    npy_float64 max_distance;
    
    npy_intp stack_size;
    npy_intp stack_max_size;
    std::vector<RR_stack_item> stack_arr;
    RR_stack_item *stack;

    void _resize_stack(const npy_intp new_max_size) {
        stack_arr.resize(new_max_size);
        stack = &stack_arr[0];
        stack_max_size = new_max_size;
    };
    
    RectRectDistanceTracker(const ckdtree *_tree, 
                 const Rectangle& _rect1, const Rectangle& _rect2,
                 const npy_float64 _p, const npy_float64 eps, 
                 const npy_float64 _upper_bound)
        : tree(_tree), rect1(_rect1), rect2(_rect2), stack_arr(8) {
    
        if (rect1.m != rect2.m) {
            const char *msg = "rect1 and rect2 have different dimensions";
            throw std::invalid_argument(msg); // raises ValueError
        }
        
        p = _p;
        
        /* internally we represent all distances as distance ** p */
        if (NPY_LIKELY(p == 2.0))
            upper_bound = _upper_bound * _upper_bound;
        else if ((!ckdtree_isinf(p)) && (!ckdtree_isinf(_upper_bound)))
            upper_bound = std::pow(_upper_bound,p);
        else
            upper_bound = _upper_bound;
        
        /* fiddle approximation factor */
        if (NPY_LIKELY(p == 2.0)) {
            npy_float64 tmp = 1. + eps;
            epsfac = 1. / (tmp*tmp);
        }
        else if (eps == 0.)
            epsfac = 1.;
        else if (ckdtree_isinf(p)) 
            epsfac = 1. / (1. + eps);
        else
            epsfac = 1. / std::pow((1. + eps), p);

        stack = &stack_arr[0];
        stack_max_size = 8;
        stack_size = 0;

        /* Compute initial min and max distances */
        MinMaxDist::rect_rect_p(tree, rect1, rect2, p, &min_distance, &max_distance);

    };
    

    void push(const npy_intp which, const npy_intp direction,
              const npy_intp split_dim, const npy_float64 split_val) {
        
        const npy_float64 p = this->p;
        
        Rectangle *rect;
        if (which == 1)
            rect = &rect1;
        else
            rect = &rect2;

        /* push onto stack */
        if (stack_size == stack_max_size)
            _resize_stack(stack_max_size * 2);
            
        RR_stack_item *item = &stack[stack_size];
        ++stack_size;
        item->which = which;
        item->split_dim = split_dim;
        item->min_distance = min_distance;
        item->max_distance = max_distance;
        item->min_along_dim = rect->mins[split_dim];
        item->max_along_dim = rect->maxes[split_dim];

        /* update min/max distances */
        npy_float64 min, max;

        MinMaxDist::interval_interval_p(tree, rect1, rect2, split_dim, p, &min, &max);

        min_distance -= min;
        max_distance -= max;
        
        if (direction == LESS)
            rect->maxes[split_dim] = split_val;
        else
            rect->mins[split_dim] = split_val;

        MinMaxDist::interval_interval_p(tree, rect1, rect2, split_dim, p, &min, &max);

        min_distance += min;
        max_distance += max;

    };

    inline void push_less_of(const npy_intp which,
                                 const ckdtreenode *node) {
        push(which, LESS, node->split_dim, node->split);
    };
            
    inline void push_greater_of(const npy_intp which,
                                    const ckdtreenode *node) {
        push(which, GREATER, node->split_dim, node->split);
    };
    
    inline void pop() {
        /* pop from stack */
        --stack_size;
        
        /* assert stack_size >= 0 */
        if (NPY_UNLIKELY(stack_size < 0)) {
            const char *msg = "Bad stack size. This error should never occur.";
            throw std::logic_error(msg);
        }
        
        RR_stack_item* item = &stack[stack_size];
        min_distance = item->min_distance;
        max_distance = item->max_distance;

        if (item->which == 1) {
            rect1.mins[item->split_dim] = item->min_along_dim;
            rect1.maxes[item->split_dim] = item->max_along_dim;
        }
        else {
            rect2.mins[item->split_dim] = item->min_along_dim;
            rect2.maxes[item->split_dim] = item->max_along_dim;
        }
    };

};


#endif
