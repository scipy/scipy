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
#include <future>
#include <list>

#define tree_buffer_root(buf) (&(buf)[0][0])

static ckdtree_intp_t split_and_partition(ckdtree *self,
    ckdtree_intp_t start_idx, ckdtree_intp_t end_idx, double *maxes, double *mins,
    const int _median, const int _compact, ckdtreenode* node)
{
    const ckdtree_intp_t m = self->m;
    const double *data = self->raw_data;
    ckdtree_intp_t *indices = (intptr_t *)(self->raw_indices);

    ckdtree_intp_t p;
    double split;

    std::vector<double> tmpminmax;

    while (true) {
        if (CKDTREE_LIKELY(_compact)) {
            /* Recompute hyperrectangle bounds. This should lead to a more
             * compact kd-tree but comes at the expense of larger construction
             * time. However, construction time is usually dwarfed by the
             * query time by orders of magnitude.
             */
            const double *tmp_data_point;
            tmp_data_point = data + indices[start_idx] * m;
            for(ckdtree_intp_t i=0; i<m; ++i) {
                maxes[i] = tmp_data_point[i];
                mins[i] = tmp_data_point[i];
            }
            for (ckdtree_intp_t j = start_idx + 1; j < end_idx; ++j) {
                tmp_data_point = data + indices[j] * m;
                for(ckdtree_intp_t i=0; i<m; ++i) {
                    double tmp = tmp_data_point[i];
                    maxes[i] = maxes[i] > tmp ? maxes[i] : tmp;
                    mins[i] = mins[i] < tmp ? mins[i] : tmp;
                }
            }
        }

        /* split on the dimension with largest spread */
        ckdtree_intp_t d = 0;
        double size = 0;
        for (ckdtree_intp_t i=0; i<m; ++i) {
            if (maxes[i] - mins[i] > size) {
                d = i;
                size = maxes[i] - mins[i];
            }
        }
        const double maxval = maxes[d];
        const double minval = mins[d];
        if (maxval == minval) {
            /* all points are identical; warn user?
             * return leafnode
             */
            node->split_dim = -1;
            return -1;
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

            if (tmpminmax.empty())
            {
                tmpminmax.resize(2*m);
                std::copy_n(mins, m, tmpminmax.begin());
                std::copy_n(maxes, m, tmpminmax.begin()+m);
                mins = tmpminmax.data();
                maxes = tmpminmax.data()+m;
            }

            const auto fixed_val = data[indices[start_idx]*m + d];
            mins[d] = fixed_val;
            maxes[d] = fixed_val;
            continue;
        }

        node->split = split;
        node->split_dim = d;
        return p;
    }
}

inline void init_ckdtreenode(ckdtreenode& node, ckdtree_intp_t start_idx, intptr_t end_idx)
{
    node.split_dim = -1;
    node.children = end_idx - start_idx;
    node.split = 0.;
    node.start_idx = start_idx;
    node.end_idx = end_idx;
    node.less = node.greater = nullptr;
    node._less = node._greater = 0;
}

static ckdtree_intp_t
build(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
      double *maxes, double *mins,
      const int _median, const int _compact,
      std::vector<ckdtreenode>& buf)
{
    /* put a new node into the node stack */
    const ckdtree_intp_t node_index = buf.size();
    ckdtree_intp_t _less, _greater = 0;
    buf.emplace_back();
    ckdtreenode *node = &buf.back();
    init_ckdtreenode(*node, start_idx, end_idx);

    if (end_idx-start_idx <= self->leafsize) {
        /* below brute force limit, return leafnode */
        node->split_dim = -1;
        return node_index;
    }
    else {

        const ckdtree_intp_t p = split_and_partition(self, start_idx, end_idx, maxes, mins, _median, _compact, node);
        if(p < 0)
            return node_index; // split_and_partition decided to have a leaf node

        const double split = node->split;
        const ckdtree_intp_t d = node->split_dim;

        if (CKDTREE_LIKELY(_compact)) {
            _less = build(self, start_idx, p, maxes, mins, _median, _compact, buf);
            _greater = build(self, p, end_idx, maxes, mins, _median, _compact, buf);
        }
        else
        {
            std::vector<double> mids(self->m);

            std::copy_n(maxes, mids.size(), mids.begin());
            mids[d] = split;
            _less = build(self, start_idx, p, mids.data(), mins, _median, _compact, buf);

            std::copy_n(mins, mids.size(), mids.begin());
            mids[d] = split;
            _greater = build(self, p, end_idx, maxes, mids.data(), _median, _compact, buf);
        }

        node = &buf[node_index];
        /* fill in entries */
        node->_less = _less;
        node->less = buf.data() + node->_less;
        node->_greater = _greater;
        node->greater = buf.data() + node->_greater;

        return node_index;
    }
}

static ckdtree_intp_t
build(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
    double *maxes, double *mins,
    const int _median, const int _compact,
    std::list<std::vector<ckdtreenode>>* buflist, int workers);

struct build_worker {

    std::list<std::vector<ckdtreenode>> buflist_;
    std::vector<double> mins_, maxes_;
    std::future<ckdtree_intp_t> future_;

    ckdtree_intp_t run(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
          const double *maxes, const double *mins,
          const int _median, const int _compact,
          int workers)
    {
        const int subworkers = workers/2;
        assert(subworkers > 0);

        mins_.resize(self->m);
        std::copy_n(mins, mins_.size(), mins_.begin());
        maxes_.resize(self->m);
        std::copy_n(maxes, maxes_.size(), maxes_.begin());

        future_ = std::async(std::launch::async, [=]() {
                return build(
                    self, start_idx, end_idx, maxes_.data(), mins_.data(), _median, _compact,
                    &buflist_, subworkers);
            });

        return subworkers;
    }

    operator bool() const { return future_.valid(); }
    void wait() { future_.wait(); }
    std::list<std::vector<ckdtreenode>>& buflist() { return buflist_; }
};

static ckdtree_intp_t
build(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
      double *maxes, double *mins,
      const int _median, const int _compact,
      std::list<std::vector<ckdtreenode>>* buflist, int workers)
{
    buflist->emplace_back();
    std::vector<ckdtreenode>& buf = buflist->front();
    buf.reserve((end_idx-start_idx)/self->leafsize/workers);

    if (workers == 1)
        return build(self, start_idx, end_idx, maxes, mins, _median, _compact, buf);

    /* put a new node into the node stack */
    const ckdtree_intp_t node_index = buf.size();
    ckdtree_intp_t _less;
    buf.emplace_back();
    ckdtreenode *node = &buf.back();
    init_ckdtreenode(*node, start_idx, end_idx);

    if (end_idx-start_idx <= self->leafsize) {
        /* below brute force limit, return leafnode */
        node->split_dim = -1;
        return node_index;
    }
    else {

        const ckdtree_intp_t p = split_and_partition(self, start_idx, end_idx, maxes, mins,
            _median, _compact, node);
        if(p < 0)
            return node_index; // split_and_partition decided to have a leaf node

        const double split = node->split;
        const ckdtree_intp_t d = node->split_dim;

        // try to put this on its own cache line to avoid false sharing with the child thread
        alignas(64) build_worker worker;

        if (CKDTREE_LIKELY(_compact)) {
            workers -= worker.run(self, p, end_idx, maxes, mins, _median, _compact, workers);
            _less = build(self, start_idx, p, maxes, mins, _median, _compact, buflist, workers);
        }
        else
        {
            std::vector<double> mids(self->m);

            std::copy_n(mins, mids.size(), mids.begin());
            mids[d] = split;
            workers -= worker.run(self, p, end_idx, maxes, mids.data(), _median, _compact, workers);

            std::copy_n(maxes, mids.size(), mids.begin());
            mids[d] = split;
            _less = build(self, start_idx, p, mids.data(), mins, _median, _compact, buflist, workers);
        }

        node = &buf[node_index];
        /* fill in entries */
        node->_less = _less;

        ckdtree_intp_t _greater = 0;
        for(const auto& b : *buflist)
            _greater += b.size();
        worker.wait();
        buflist->splice(buflist->end(), worker.buflist());
        node->_greater = _greater;

        return node_index;
    }
}

int build_ckdtree(ckdtree *self, ckdtree_intp_t start_idx, intptr_t end_idx,
              double *maxes, double *mins, int _median, int _compact, int _workers)

{
    if(_workers == 1) {
        self->tree_buffer->reserve((end_idx-start_idx)/self->leafsize);
        // build tree sequentially
        build(self, start_idx, end_idx, maxes, mins, _median, _compact, *self->tree_buffer);
    }
    else {
        std::list<std::vector<ckdtreenode>> buflist;
        // build tree with sections collected in buflist
        build(self, start_idx, end_idx, maxes, mins, _median, _compact, &buflist, _workers);

        {
            // count the number of nodes we got
            ckdtree_intp_t tree_size = 0;
            for(const auto& b : buflist) {
                tree_size += b.size();
            }
            self->tree_buffer->resize(tree_size);
        }

        // merge sections from buflist in self->tree_buffer and update pointers and indices
        std::vector<std::future<void>> futures;
        futures.reserve(_workers-1);
        auto buflistLast = buflist.end();
        buflistLast--;

        auto insertBuf = [self](int offset, std::vector<ckdtreenode>* b) -> void {
                auto s = self->tree_buffer->begin()+offset;
                for(const auto& a : *b) {
                    *s = a;
                    s->_less += offset;
                    s->_greater += offset;
                    s->less = self->tree_buffer->data() + s->_less;
                    s->greater = self->tree_buffer->data() + s->_greater;
                    ++s;
                }
                // std::copy(b->cbegin(), b->cend(), self->tree_buffer->begin()+offset);
            };
        ckdtree_intp_t offset = 0;
        for(auto b = buflist.begin(); b != buflistLast; ++b) {
            futures.push_back(std::async(std::launch::async, insertBuf, offset, &*b));
            offset += b->size();
        }
        insertBuf(offset, &*buflistLast);

        // wait for threads for clarity, the destructor would take care of this anyway...
        for(auto& f : futures) {
            f.wait();
        }
    }

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
