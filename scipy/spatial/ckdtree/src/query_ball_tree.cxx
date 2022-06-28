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

#include "ckdtree_decl.h"
#include "rectangle.h"


static void
traverse_no_checking(const ckdtree *self, const ckdtree *other,
                     std::vector<ckdtree_intp_t> *results,
                     const ckdtreenode *node1, const ckdtreenode *node2)
{
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    const ckdtree_intp_t *sindices = self->raw_indices;
    const ckdtree_intp_t *oindices = other->raw_indices;
    ckdtree_intp_t i, j;

    if (node1->split_dim == -1) {   /* leaf node */
        lnode1 = node1;

        if (node2->split_dim == -1) {  /* leaf node */
            lnode2 = node2;

            const ckdtree_intp_t start1 = lnode1->start_idx;
            const ckdtree_intp_t start2 = lnode2->start_idx;
            const ckdtree_intp_t end1 = lnode1->end_idx;
            const ckdtree_intp_t end2 = lnode2->end_idx;

            for (i = start1; i < end1; ++i) {
                auto &results_i = results[sindices[i]];
                for (j = start2; j < end2; ++j)
                    results_i.push_back(oindices[j]);
            }
        }
        else {
            traverse_no_checking(self, other, results, node1, node2->less);
            traverse_no_checking(self, other, results, node1, node2->greater);
        }
    }
    else {
        traverse_no_checking(self, other, results, node1->less, node2);
        traverse_no_checking(self, other, results, node1->greater, node2);
    }
}


template <typename MinMaxDist> static void
traverse_checking(const ckdtree *self, const ckdtree *other,
                  std::vector<ckdtree_intp_t> *results,
                  const ckdtreenode *node1, const ckdtreenode *node2,
                  RectRectDistanceTracker<MinMaxDist> *tracker)
{
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    double d;
    ckdtree_intp_t i, j;

    if (tracker->min_distance > tracker->upper_bound * tracker->epsfac)
        return;
    else if (tracker->max_distance < tracker->upper_bound / tracker->epsfac)
        traverse_no_checking(self, other, results, node1, node2);
    else if (node1->split_dim == -1) { /* 1 is leaf node */
        lnode1 = node1;

        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */

            /* brute-force */
            lnode2 = node2;
            const double p = tracker->p;
            const double tub = tracker->upper_bound;
            const double tmd = tracker->max_distance;
            const double *sdata = self->raw_data;
            const ckdtree_intp_t *sindices = self->raw_indices;
            const double *odata = other->raw_data;
            const ckdtree_intp_t *oindices = other->raw_indices;
            const ckdtree_intp_t m = self->m;
            const ckdtree_intp_t start1 = lnode1->start_idx;
            const ckdtree_intp_t start2 = lnode2->start_idx;
            const ckdtree_intp_t end1 = lnode1->end_idx;
            const ckdtree_intp_t end2 = lnode2->end_idx;

            CKDTREE_PREFETCH(sdata + sindices[start1] * m, 0, m);

            if (start1 < end1 - 1)
                CKDTREE_PREFETCH(sdata + sindices[start1+1] * m, 0, m);

            for (i = start1; i < end1; ++i) {

                if (i < end1 - 2)
                    CKDTREE_PREFETCH(sdata + sindices[i+2] * m, 0, m);

                CKDTREE_PREFETCH(odata + oindices[start2] * m, 0, m);

                if (start2 < end2 - 1)
                    CKDTREE_PREFETCH(odata + oindices[start2+1] * m, 0, m);

                auto &results_i = results[sindices[i]];

                for (j = start2; j < end2; ++j) {

                    if (j < end2 - 2)
                        CKDTREE_PREFETCH(odata + oindices[j+2] * m, 0, m);

                    d = MinMaxDist::point_point_p(
                            self,
                            sdata + sindices[i] * m,
                            odata + oindices[j] * m,
                            p, m, tmd);

                    if (d <= tub)
                        results_i.push_back(other->raw_indices[j]);
                }
            }

        }
        else { /* 1 is a leaf node, 2 is inner node */

            tracker->push_less_of(2, node2);
            traverse_checking(
                self, other, results, node1, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse_checking(
                self, other, results, node1, node2->greater, tracker);
            tracker->pop();
        }
    }
    else {  /* 1 is an inner node */
        if (node2->split_dim == -1) { /* 1 is an inner node, 2 is a leaf node */
            tracker->push_less_of(1, node1);
            traverse_checking(
                self, other, results, node1->less, node2, tracker);
            tracker->pop();

            tracker->push_greater_of(1, node1);
            traverse_checking(
                self, other, results, node1->greater, node2, tracker);
            tracker->pop();
        }
        else { /* 1 & 2 are inner nodes */

            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse_checking(
                self, other, results, node1->less, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse_checking(
                self, other, results, node1->less, node2->greater, tracker);
            tracker->pop();
            tracker->pop();


            tracker->push_greater_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse_checking(
                self, other, results, node1->greater, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse_checking(
                self, other, results, node1->greater, node2->greater, tracker);
            tracker->pop();
            tracker->pop();
        }
    }
}

int
query_ball_tree(const ckdtree *self, const ckdtree *other,
                const double r, const double p, const double eps,
                std::vector<ckdtree_intp_t> *results)
{

#define HANDLE(cond, kls) \
    if(cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, eps, r); \
        traverse_checking(self, other, results, self->ctree, other->ctree, \
            &tracker); \
    } else

    Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
    Rectangle r2(other->m, other->raw_mins, other->raw_maxes);

    if(CKDTREE_LIKELY(self->raw_boxsize_data == NULL)) {
        HANDLE(CKDTREE_LIKELY(p == 2), MinkowskiDistP2)
        HANDLE(p == 1, MinkowskiDistP1)
        HANDLE(ckdtree_isinf(p), MinkowskiDistPinf)
        HANDLE(1, MinkowskiDistPp)
        {}
    } else {
        HANDLE(CKDTREE_LIKELY(p == 2), BoxMinkowskiDistP2)
        HANDLE(p == 1, BoxMinkowskiDistP1)
        HANDLE(ckdtree_isinf(p), BoxMinkowskiDistPinf)
        HANDLE(1, BoxMinkowskiDistPp)
        {}
    }

    for (ckdtree_intp_t i = 0; i < self->n; ++i) {
        std::sort(results[i].begin(), results[i].end());
    }

    return 0;
}
