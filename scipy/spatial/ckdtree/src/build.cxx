#include "ckdtree_decl.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>

#define tree_buffer_root(buf) (&(buf)[0][0])

static ckdtree_intp_t
build(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
      double *maxes, double *mins,
      const int _median, const int _compact)
{

    const ckdtree_intp_t m = self->m;
    const double *data = self->raw_data;
    ckdtree_intp_t *indices = (intptr_t *)(self->raw_indices);

    ckdtreenode new_node, *n, *root;
    ckdtree_intp_t node_index, _less, _greater;
    ckdtree_intp_t i, j, p, d;
    double size, split, minval, maxval;

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

        if (CKDTREE_LIKELY(_compact)) {
            /* Recompute hyperrectangle bounds. This should lead to a more
             * compact kd-tree but comes at the expense of larger construction
             * time. However, construction time is usually dwarfed by the
             * query time by orders of magnitude.
             */
            const double *tmp_data_point;
            tmp_data_point = data + indices[start_idx] * m;
            for(i=0; i<m; ++i) {
                maxes[i] = tmp_data_point[i];
                mins[i] = tmp_data_point[i];
            }
            for (j = start_idx + 1; j < end_idx; ++j) {
                tmp_data_point = data + indices[j] * m;
                for(i=0; i<m; ++i) {
                    double tmp = tmp_data_point[i];
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
        auto index_compare = [=](ckdtree_intp_t a, ckdtree_intp_t b) {
            const double point_a = data[a * m + d];
            const double point_b = data[b * m + d];
            return point_a < point_b;
        };

        auto partition_pivot = [=](ckdtree_intp_t* first, ckdtree_intp_t* last, double pivot) {
            const auto partition_ptr = std::partition(
                first, last,
                [&](ckdtree_intp_t a) { return data[a * m + d] < pivot; });
            return partition_ptr - indices;
        };

        if (CKDTREE_LIKELY(_median)) {
            /* split on median to create a balanced tree
             * adopted from scikit-learn
             */
            const auto n_points = end_idx - start_idx;
            auto* node_indices = indices + start_idx;
            auto mid = node_indices + n_points / 2;
            std::nth_element(
                node_indices, mid, node_indices + n_points, index_compare);

            split = data[*mid * m + d];
            p = partition_pivot(node_indices, mid, split);
        }
        else {
            /* split with the sliding midpoint rule */
            split = (maxval + minval) / 2;

            p = partition_pivot(indices + start_idx, indices + end_idx, split);
        }

        /* slide midpoint if necessary */
        if (p == start_idx) {
            /* no points less than split */
            auto min_idx = *std::min_element(
                indices + start_idx, indices + end_idx, index_compare);
            split = std::nextafter(data[min_idx * m + d], HUGE_VAL);
            p = partition_pivot(indices + start_idx, indices + end_idx, split);
        }
        else if (p == end_idx) {
            /* no points greater than split */
            auto max_idx = *std::max_element(
                indices + start_idx, indices + end_idx, index_compare);
            split = data[max_idx * m + d];
            p = partition_pivot(indices + start_idx, indices + end_idx, split);
        }

        if (CKDTREE_UNLIKELY(p == start_idx || p == end_idx)) {
            // All children are equal in this dimenion, try again with new bounds
            assert(!_compact);
            self->tree_buffer->pop_back();
            std::vector<double> tmp_bounds(2 * m);
            double* tmp_mins = &tmp_bounds[0];
            std::copy_n(mins, m, tmp_mins);
            double* tmp_maxes = &tmp_bounds[m];
            std::copy_n(maxes, m, tmp_maxes);

            const auto fixed_val = data[indices[start_idx]*m + d];
            tmp_mins[d] = fixed_val;
            tmp_maxes[d] = fixed_val;

            return build(self, start_idx, end_idx, tmp_maxes, tmp_mins, _median, _compact);
        }

        if (CKDTREE_LIKELY(_compact)) {
            _less = build(self, start_idx, p, maxes, mins, _median, _compact);
            _greater = build(self, p, end_idx, maxes, mins, _median, _compact);
        }
        else
        {
            std::vector<double> tmp(m);
            double *mids = &tmp[0];

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



int build_ckdtree(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
              double *maxes, double *mins, int _median, int _compact)

{
    build(self, start_idx, end_idx, maxes, mins, _median, _compact);
    return 0;
}

static double
add_weights(ckdtree *self,
           double *node_weights,
           ckdtree_intp_t node_index,
           double *weights)
{

    ckdtree_intp_t *indices = (intptr_t *)(self->raw_indices);

    ckdtreenode *n, *root;

    root = tree_buffer_root(self->tree_buffer);

    n = root + node_index;

    double sum = 0;

    if (n->split_dim != -1) {
        /* internal nodes; recursively calculate the total weight */
        double left, right;
        left = add_weights(self, node_weights, n->_less, weights);
        right = add_weights(self, node_weights, n->_greater, weights);
        sum = left + right;
    } else {
        ckdtree_intp_t i;

        /* Leaf nodes */
        for (i = n->start_idx; i < n->end_idx; ++i) {
            sum += weights[indices[i]];
        }
    }

    node_weights[node_index] = sum;
    return sum;
}

int
build_weights (ckdtree *self, double *node_weights, double *weights)
{

    add_weights(self, node_weights, 0, weights);
    return 0;
}
