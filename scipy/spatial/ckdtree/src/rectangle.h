
#ifndef CKDTREE_CPP_RECTANGLE
#define CKDTREE_CPP_RECTANGLE

#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>
#include <cmath>
#include <cstring>


/* Interval arithmetic
 * ===================
 */

struct Rectangle {

    const ckdtree_intp_t m;

    /* the last const is to allow const Rectangle to use these functions;
     * also notice we had to mark buf mutable to avoid writing non const version
     * of the same accessors. */
    double * const maxes() const { return &buf[0]; }
    double * const mins() const { return &buf[0] + m; }

    Rectangle(const ckdtree_intp_t _m,
              const double *_mins,
              const double *_maxes) : m(_m), buf(2 * m) {

        /* copy array data */
        /* FIXME: use std::vector ? */
        std::memcpy((void*)mins(), (void*)_mins, m*sizeof(double));
        std::memcpy((void*)maxes(), (void*)_maxes, m*sizeof(double));
    };

    Rectangle(const Rectangle& rect) : m(rect.m), buf(rect.buf) {};

    private:
        mutable std::vector<double> buf;
};

#include "distance.h"

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
    ckdtree_intp_t    which;
    ckdtree_intp_t    split_dim;
    double min_along_dim;
    double max_along_dim;
    double min_distance;
    double max_distance;
};

const ckdtree_intp_t LESS = 1;
const ckdtree_intp_t GREATER = 2;

template<typename MinMaxDist>
    struct RectRectDistanceTracker {

    const ckdtree * tree;
    Rectangle rect1;
    Rectangle rect2;
    double p;
    double epsfac;
    double upper_bound;
    double min_distance;
    double max_distance;

    ckdtree_intp_t stack_size;
    ckdtree_intp_t stack_max_size;
    std::vector<RR_stack_item> stack_arr;
    RR_stack_item *stack;

    /* if min/max distance / adjustment is less than this,
     * we believe the incremental tracking is inaccurate */
    double inaccurate_distance_limit;

    void _resize_stack(const ckdtree_intp_t new_max_size) {
        stack_arr.resize(new_max_size);
        stack = &stack_arr[0];
        stack_max_size = new_max_size;
    };

    RectRectDistanceTracker(const ckdtree *_tree,
                 const Rectangle& _rect1, const Rectangle& _rect2,
                 const double _p, const double eps,
                 const double _upper_bound)
        : tree(_tree), rect1(_rect1), rect2(_rect2), stack_arr(8) {

        if (rect1.m != rect2.m) {
            const char *msg = "rect1 and rect2 have different dimensions";
            throw std::invalid_argument(msg); // raises ValueError
        }

        p = _p;

        /* internally we represent all distances as distance ** p */
        if (CKDTREE_LIKELY(p == 2.0))
            upper_bound = _upper_bound * _upper_bound;
        else if ((!ckdtree_isinf(p)) && (!isinf(_upper_bound)))
            upper_bound = std::pow(_upper_bound,p);
        else
            upper_bound = _upper_bound;

        /* fiddle approximation factor */
        if (CKDTREE_LIKELY(p == 2.0)) {
            double tmp = 1. + eps;
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
        if(ckdtree_isinf(max_distance)) {
            const char *msg = "Encountering floating point overflow. "
                              "The value of p too large for this dataset; "
                              "For such large p, consider using the special case p=np.inf . ";
            throw std::invalid_argument(msg); // raises ValueError
        }
        inaccurate_distance_limit = max_distance;
    };


    void push(const ckdtree_intp_t which, const intptr_t direction,
              const ckdtree_intp_t split_dim, const double split_val) {

        const double p = this->p;
        /* subnomial is 1 if round-off is expected to taint the incremental distance tracking.
         * in that case we always recompute the distances.
         * Recomputing costs more calls to pow, thus if the round-off error does not seem
         * to wipe out the value, then we still do the incremental update.
         * */
        int subnomial = 0;

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
        item->min_along_dim = rect->mins()[split_dim];
        item->max_along_dim = rect->maxes()[split_dim];

        /* update min/max distances */
        double min1, max1;
        double min2, max2;

        MinMaxDist::interval_interval_p(tree, rect1, rect2, split_dim, p, &min1, &max1);

        if (direction == LESS)
            rect->maxes()[split_dim] = split_val;
        else
            rect->mins()[split_dim] = split_val;

        MinMaxDist::interval_interval_p(tree, rect1, rect2, split_dim, p, &min2, &max2);

        subnomial = subnomial || (min_distance < inaccurate_distance_limit || max_distance < inaccurate_distance_limit);

        subnomial = subnomial || ((min1 != 0 && min1 < inaccurate_distance_limit) || max1 < inaccurate_distance_limit);
        subnomial = subnomial || ((min2 != 0 && min2 < inaccurate_distance_limit) || max2 < inaccurate_distance_limit);
        subnomial = subnomial || (min_distance < inaccurate_distance_limit || max_distance < inaccurate_distance_limit);

        if (CKDTREE_UNLIKELY(subnomial)) {
            MinMaxDist::rect_rect_p(tree, rect1, rect2, p, &min_distance, &max_distance);
        } else {
            min_distance += (min2 - min1);
            max_distance += (max2 - max1);
        }
    };

    inline void push_less_of(const ckdtree_intp_t which,
                                 const ckdtreenode *node) {
        push(which, LESS, node->split_dim, node->split);
    };

    inline void push_greater_of(const ckdtree_intp_t which,
                                    const ckdtreenode *node) {
        push(which, GREATER, node->split_dim, node->split);
    };

    inline void pop() {
        /* pop from stack */
        --stack_size;

        /* assert stack_size >= 0 */
        if (CKDTREE_UNLIKELY(stack_size < 0)) {
            const char *msg = "Bad stack size. This error should never occur.";
            throw std::logic_error(msg);
        }

        RR_stack_item* item = &stack[stack_size];
        min_distance = item->min_distance;
        max_distance = item->max_distance;

        if (item->which == 1) {
            rect1.mins()[item->split_dim] = item->min_along_dim;
            rect1.maxes()[item->split_dim] = item->max_along_dim;
        }
        else {
            rect2.mins()[item->split_dim] = item->min_along_dim;
            rect2.maxes()[item->split_dim] = item->max_along_dim;
        }
    };

};


#endif
