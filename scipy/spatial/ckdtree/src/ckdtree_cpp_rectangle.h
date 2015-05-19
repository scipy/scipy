
#ifndef CKDTREE_CPP_RECTANGLE
#define CKDTREE_CPP_RECTANGLE

#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>
#include <cmath>
#include <cstring>



extern npy_float64 infinity;

// Interval arithmetic
// ===================

struct Rectangle {
    
    npy_intp m;
    npy_float64 *mins;
    npy_float64 *maxes;
    
    std::vector<npy_float64> mins_arr;
    std::vector<npy_float64> maxes_arr;

    Rectangle(const npy_intp _m, 
              const npy_float64 *_mins, 
              const npy_float64 *_maxes) : mins_arr(_m), maxes_arr(_m) {

        // Copy array data
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


// 1-d pieces
// These should only be used if p != infinity
inline npy_float64 
min_dist_point_interval_p(const npy_float64 *x, const Rectangle& rect,
                          const npy_intp k, const npy_float64 p)
{
    // Compute the minimum distance along dimension k between x and
    // a point in the hyperrectangle.
    return std::pow(dmax(0, dmax(rect.mins[k] - x[k], x[k] - rect.maxes[k])),p);
}

inline npy_float64 
max_dist_point_interval_p(const npy_float64 *x, const Rectangle& rect,
                          const npy_intp k, const npy_float64 p)
{
    // Compute the maximum distance along dimension k between x and
    // a point in the hyperrectangle.
    return std::pow(dmax(rect.maxes[k] - x[k], x[k] - rect.mins[k]),p);
}


inline npy_float64 
min_dist_interval_interval_p(const Rectangle& rect1, const Rectangle& rect2,
                             const npy_intp k, const npy_float64 p)
{
    // Compute the minimum distance along dimension k between points in
    // two hyperrectangles.
    return std::pow(dmax(0, dmax(rect1.mins[k] - rect2.maxes[k],
                          rect2.mins[k] - rect1.maxes[k])),p);
}


inline npy_float64 
max_dist_interval_interval_p(const Rectangle& rect1, const Rectangle& rect2,
                             const npy_intp k, const npy_float64 p)
{                             
    // Compute the maximum distance along dimension k between points in
    // two hyperrectangles.
    return std::pow(dmax(rect1.maxes[k] - rect2.mins[k], 
                          rect2.maxes[k] - rect1.mins[k]),p);
}


// Interval arithmetic in m-D
// ==========================

// These should be used only for p == infinity

inline npy_float64 
min_dist_point_rect_p_inf(const npy_float64 *x, const Rectangle& rect) 
{
    // Compute the minimum distance between x and the given hyperrectangle.
    npy_intp i;
    npy_float64 min_dist = 0.;
    for (i=0; i<rect.m; ++i) {
        min_dist = dmax(min_dist, dmax(rect.mins[i]-x[i], x[i]-rect.maxes[i]));
    }
    return min_dist;
}

inline npy_float64 
max_dist_point_rect_p_inf(const npy_float64 *x, const Rectangle& rect)
{
    // Compute the maximum distance between x and the given hyperrectangle.
    npy_intp i;
    npy_float64 max_dist = 0.;
    for (i=0; i<rect.m; ++i) {
        max_dist = dmax(max_dist, dmax(rect.maxes[i]-x[i], x[i]-rect.mins[i]));
    }
    return max_dist;
}

inline npy_float64 
min_dist_rect_rect_p_inf(const Rectangle& rect1, const Rectangle& rect2)
{
    // Compute the minimum distance between points in two hyperrectangles.
    npy_intp i;
    npy_float64 min_dist = 0.;
    for (i=0; i<rect1.m; ++i) {
        min_dist = dmax(min_dist, dmax(rect1.mins[i] - rect2.maxes[i],
                                       rect2.mins[i] - rect1.maxes[i]));
    }                                   
    return min_dist;
}


inline npy_float64 
max_dist_rect_rect_p_inf(const Rectangle& rect1, const Rectangle rect2)
{
    // Compute the maximum distance between points in two hyperrectangles.
    npy_intp i;
    npy_float64 max_dist = 0.;
    for (i=0; i<rect1.m; ++i) {
        max_dist = dmax(max_dist, dmax(rect1.maxes[i] - rect2.mins[i],
                                       rect2.maxes[i] - rect1.mins[i]));
    }
    return max_dist;
}


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

struct RectRectDistanceTracker {
    
    Rectangle rect1; 
    Rectangle rect2;
    npy_float64 p; 
    npy_float64 epsfac;
    npy_float64 upper_bound;
    npy_float64 min_distance;
    npy_float64 max_distance;
    
    npy_float64 infinity;

    npy_intp stack_size;
    npy_intp stack_max_size;
    std::vector<RR_stack_item> stack_arr;
    RR_stack_item *stack;

    void _resize_stack(const npy_intp new_max_size) {
        stack_arr.resize(new_max_size);
        stack = &stack_arr[0];
        stack_max_size = new_max_size;
    };
    
    RectRectDistanceTracker(const Rectangle& _rect1, const Rectangle& _rect2,
                 const npy_float64 _p, const npy_float64 eps, 
                 const npy_float64 _upper_bound)
        : rect1(_rect1), rect2(_rect2), stack_arr(8) {
    
        infinity = ::infinity;
    
        if (rect1.m != rect2.m) {
            const char *msg = "rect1 and rect2 have different dimensions";
            throw std::invalid_argument(msg); // raises ValueError
        }
        
        p = _p;
        
        // internally we represent all distances as distance ** p
        if ((p != infinity) && (_upper_bound != infinity)) 
            upper_bound = std::pow(_upper_bound,p);
        else
            upper_bound = _upper_bound;
        
        // fiddle approximation factor
        if (eps == 0.)
            epsfac = 1.;
        else if (p == infinity) 
            epsfac = 1. / (1. + eps);
        else
            epsfac = 1. / std::pow((1. + eps), p);

        stack = &stack_arr[0];
        stack_max_size = 8;
        stack_size = 0;

        // Compute initial min and max distances
        if (p == infinity) {
            min_distance = min_dist_rect_rect_p_inf(rect1, rect2);
            max_distance = max_dist_rect_rect_p_inf(rect1, rect2);
        }
        else {
            min_distance = 0.;
            max_distance = 0.;
            for(npy_intp i=0; i<rect1.m; ++i) {
                min_distance 
                    += min_dist_interval_interval_p(rect1, rect2, i, p);
                max_distance 
                    += max_dist_interval_interval_p(rect1, rect2, i, p);
            }
        }
    };
    

    void push(const npy_intp which, const npy_intp direction,
              const npy_intp split_dim, const npy_float64 split_val) {
              
        Rectangle *rect;
        if (which == 1)
            rect = &rect1;
        else
            rect = &rect2;

        // Push onto stack
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

        // Update min/max distances
        if (p != infinity) {
            min_distance 
                -= min_dist_interval_interval_p(rect1, rect2, split_dim, p);
            max_distance 
                -= max_dist_interval_interval_p(rect1, rect2, split_dim, p);
        }
        
        if (direction == LESS)
            rect->maxes[split_dim] = split_val;
        else
            rect->mins[split_dim] = split_val;

        if (p != infinity) {
            min_distance 
                += min_dist_interval_interval_p(rect1, rect2, split_dim, p);
            max_distance 
                += max_dist_interval_interval_p(rect1, rect2, split_dim, p);
        }
        else {
            min_distance = min_dist_rect_rect_p_inf(rect1, rect2);
            max_distance = max_dist_rect_rect_p_inf(rect1, rect2);
        }     
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
        // Pop from stack
        --stack_size;
        
        // assert stack_size >= 0
        if (stack_size < 0) {
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


/*
 * Point-to-rectangle distance tracker
 * ===================================
 *
 * The other logical unit that is used in query_ball_point is to keep track
 * of the maximum and minimum distances between points in a hyperrectangle
 * and another fixed point as the rectangle is successively split.
 *
 * Example
 * -------
 * # node encloses points in rect
 *
 * cdef PointRectDistanceTracker dist_tracker
 * dist_tracker = PointRectDistanceTracker(pt, rect, p)
 *
 * ...
 *
 * if dist_tracker.min_distance < ...:
 *     ...
 *
 * dist_tracker.push_less_of(node)
 * do_something(node.less, dist_tracker)
 * dist_tracker.pop()
 *
 * dist_tracker.push_greater_of(node)
 * do_something(node.greater, dist_tracker)
 * dist_tracker.pop()
 */
 

struct RP_stack_item {
    npy_intp split_dim;
    npy_float64 min_along_dim;
    npy_float64 max_along_dim;
    npy_float64 min_distance;
    npy_float64 max_distance;
};


struct PointRectDistanceTracker {

    Rectangle   rect;
    const npy_float64 *pt;
    npy_float64 p; 
    npy_float64 epsfac; 
    npy_float64 upper_bound;
    npy_float64 min_distance;
    npy_float64 max_distance;
    
    npy_float64 infinity;

    npy_intp stack_size;
    npy_intp stack_max_size;
    std::vector<RP_stack_item> stack_arr;
    RP_stack_item *stack;
    
    void _resize_stack(const npy_intp new_max_size) {
        stack_arr.resize(new_max_size);
        stack = &stack_arr[0];
        stack_max_size = new_max_size;
    };
    
    PointRectDistanceTracker(const npy_float64 *_pt, const Rectangle& _rect,
              const npy_float64 _p, const npy_float64 eps, 
              const npy_float64 _upper_bound)
        : rect(_rect), stack_arr(8) {

        infinity = ::infinity;

        pt = _pt;
        rect = _rect;
        p = _p;
        
        // internally we represent all distances as distance ** p
        if ((p != infinity) && (_upper_bound != infinity))
            upper_bound = std::pow(_upper_bound,p);
        else
            upper_bound = _upper_bound;

        // fiddle approximation factor
        if (eps == 0.)
            epsfac = 1;
        else if (p == infinity)
            epsfac = 1. / (1. + eps);
        else
            epsfac = 1. / std::pow((1. + eps),p);

        // init stack
        stack = &stack_arr[0];
        stack_max_size = 8;
        stack_size = 0;

        // Compute initial min and max distances
        if (p == infinity) {
            min_distance = min_dist_point_rect_p_inf(pt, rect);
            max_distance = max_dist_point_rect_p_inf(pt, rect);
        }
        else {
            min_distance = 0.;
            max_distance = 0.;
            for(npy_intp i=0; i<rect.m; ++i) {
                min_distance += min_dist_point_interval_p(pt, rect, i, p);
                max_distance += max_dist_point_interval_p(pt, rect, i, p);
            }
        }
    };
    
    void push(const npy_intp direction, const npy_intp split_dim,
              const npy_float64 split_val) {

        // Push onto stack
        if (stack_size == stack_max_size)
            _resize_stack(stack_max_size * 2);
                
        RP_stack_item *item = &stack[stack_size];
        ++stack_size;
        item->split_dim = split_dim;
        item->min_distance = min_distance;
        item->max_distance = max_distance;
        item->min_along_dim = rect.mins[split_dim];
        item->max_along_dim = rect.maxes[split_dim];
            
        if (p != infinity) {
            min_distance -= min_dist_point_interval_p(pt, rect, split_dim, p);
            max_distance -= max_dist_point_interval_p(pt, rect, split_dim, p);
        }
        
        if (direction == LESS)
            rect.maxes[split_dim] = split_val;
        else
            rect.mins[split_dim] = split_val;

        if (p != infinity) {
            min_distance += min_dist_point_interval_p(pt, rect, split_dim, p);
            max_distance += max_dist_point_interval_p(pt, rect, split_dim, p);
        }
        else {
            min_distance = min_dist_point_rect_p_inf(pt, rect);
            max_distance = max_dist_point_rect_p_inf(pt, rect);
        }    
    };
    
    inline void push_less_of(const ckdtreenode* node) {
        push(LESS, node->split_dim, node->split);
    }
    
    inline void push_greater_of(const ckdtreenode* node) {
        push(GREATER, node->split_dim, node->split);
    }

    inline void pop() {
        // pop
        --stack_size;
        
        // assert stack_size >= 0
        if (stack_size < 0) {
            const char *msg = "Bad stack size. This error should never occur.";
            throw std::logic_error(msg);
        }
        
        RP_stack_item* item = &stack[stack_size];
        min_distance = item->min_distance;
        max_distance = item->max_distance;
        rect.mins[item->split_dim] = item->min_along_dim;
        rect.maxes[item->split_dim] = item->max_along_dim;
    };

};


#endif
