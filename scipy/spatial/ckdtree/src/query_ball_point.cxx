#include "ckdtree_decl.h"
#include "rectangle.h"

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


static void
traverse_no_checking(const ckdtree *self,
                     const int return_length,
                     std::vector<ckdtree_intp_t> &results,
                     const ckdtreenode *node)
{
    const ckdtree_intp_t *indices = self->raw_indices;
    const ckdtreenode *lnode;
    ckdtree_intp_t i;

    if (node->split_dim == -1) {  /* leaf node */
        lnode = node;
        const ckdtree_intp_t start = lnode->start_idx;
        const ckdtree_intp_t end = lnode->end_idx;
        for (i = start; i < end; ++i) {
            if (return_length) {
                results[0] ++;
            } else {
                results.push_back(indices[i]);
            }
        }
    }
    else {
        traverse_no_checking(self, return_length, results, node->less);
        traverse_no_checking(self, return_length, results, node->greater);
    }
}


template <typename MinMaxDist> static void
traverse_checking(const ckdtree *self,
                  const int return_length,
                  std::vector<ckdtree_intp_t> &results,
                  const ckdtreenode *node,
                  RectRectDistanceTracker<MinMaxDist> *tracker
)
{
    const ckdtreenode *lnode;
    double d;
    ckdtree_intp_t i;

    if (tracker->min_distance > tracker->upper_bound * tracker->epsfac) {
        return;
    }
    else if (tracker->max_distance < tracker->upper_bound / tracker->epsfac) {
        traverse_no_checking(self, return_length, results, node);
    }
    else if (node->split_dim == -1)  { /* leaf node */

        /* brute-force */
        lnode = node;
        const double p = tracker->p;
        const double tub = tracker->upper_bound;
        const double *tpt = tracker->rect1.mins();
        const double *data = self->raw_data;
        const ckdtree_intp_t *indices = self->raw_indices;
        const ckdtree_intp_t m = self->m;
        const ckdtree_intp_t start = lnode->start_idx;
        const ckdtree_intp_t end = lnode->end_idx;

        CKDTREE_PREFETCH(data + indices[start] * m, 0, m);
        if (start < end - 1)
            CKDTREE_PREFETCH(data + indices[start+1] * m, 0, m);

        for (i = start; i < end; ++i) {

            if (i < end -2 )
                CKDTREE_PREFETCH(data + indices[i+2] * m, 0, m);

            d = MinMaxDist::point_point_p(self, data + indices[i] * m, tpt, p, m, tub);

            if (d <= tub) {
                if(return_length) {
                    results[0] ++;
                } else {
                    results.push_back((ckdtree_intp_t) indices[i]);
                }
            }
        }
    }
    else {
        tracker->push_less_of(2, node);
        traverse_checking(self, return_length, results, node->less, tracker);
        tracker->pop();

        tracker->push_greater_of(2, node);
        traverse_checking(self, return_length, results, node->greater, tracker);
        tracker->pop();
    }
}

int
query_ball_point(const ckdtree *self, const double *x,
                 const double *r, const double p, const double eps,
                 const ckdtree_intp_t n_queries,
                 std::vector<ckdtree_intp_t> *results,
                 const bool return_length,
                 const bool sort_output)
{
#define HANDLE(cond, kls) \
    if(cond) { \
        if(return_length) results[i].push_back(0); \
        RectRectDistanceTracker<kls> tracker(self, point, rect, p, eps, r[i]); \
        traverse_checking(self, return_length, results[i], self->ctree, &tracker); \
    } else

    for (ckdtree_intp_t i=0; i < n_queries; ++i) {
        const ckdtree_intp_t m = self->m;
        Rectangle rect(m, self->raw_mins, self->raw_maxes);
        if (CKDTREE_LIKELY(self->raw_boxsize_data == NULL)) {
            Rectangle point(m, x + i * m, x + i * m);
            HANDLE(CKDTREE_LIKELY(p == 2), MinkowskiDistP2)
            HANDLE(p == 1, MinkowskiDistP1)
            HANDLE(std::isinf(p), MinkowskiDistPinf)
            HANDLE(1, MinkowskiDistPp)
            {}
        } else {
            Rectangle point(m, x + i * m, x + i * m);
            int j;
            for(j=0; j<m; ++j) {
                point.maxes()[j] = point.mins()[j] = BoxDist1D::wrap_position(point.mins()[j], self->raw_boxsize_data[j]);
            }
            HANDLE(CKDTREE_LIKELY(p == 2), BoxMinkowskiDistP2)
            HANDLE(p == 1, BoxMinkowskiDistP1)
            HANDLE(std::isinf(p), BoxMinkowskiDistPinf)
            HANDLE(1, BoxMinkowskiDistPp)
            {}
        }

        if (!return_length && sort_output) {
            std::sort(results[i].begin(), results[i].end());
        }
    }
    return 0;
}
