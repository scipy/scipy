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

#include "ckdtree_decl.h"
#include "ordered_pair.h"
#include "rectangle.h"


static void
traverse_no_checking(const ckdtree *self,
                     std::vector<ordered_pair> *results,
                     const ckdtreenode *node1, const ckdtreenode *node2)
{
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    ckdtree_intp_t i, j, min_j;
    const ckdtree_intp_t *indices = self->raw_indices;

    if (node1->split_dim == -1) { /* leaf node */
        lnode1 = node1;

        if (node2->split_dim == -1) { /* leaf node */
            lnode2 = node2;

            const ckdtree_intp_t start1 = lnode1->start_idx;
            const ckdtree_intp_t start2 = lnode2->start_idx;
            const ckdtree_intp_t end1 = lnode1->end_idx;
            const ckdtree_intp_t end2 = lnode2->end_idx;

            for (i = start1; i < end1; ++i) {

                /* Special care here to avoid duplicate pairs */
                if (node1 == node2)
                    min_j = i + 1;
                else
                    min_j = start2;

                for (j = min_j; j < end2; ++j)
                    add_ordered_pair(results, indices[i], indices[j]);
            }
        }
        else {
            traverse_no_checking(self, results, node1, node2->less);
            traverse_no_checking(self, results, node1, node2->greater);
        }
    }
    else {
        if (node1 == node2) {
            /*
             * Avoid traversing (node1->less, node2->greater) and
             * (node1->greater, node2->less) (it's the same node pair twice
             * over, which is the source of the complication in the
             * original KDTree.query_pairs)
             */
            traverse_no_checking(self, results, node1->less, node2->less);
            traverse_no_checking(self, results, node1->less, node2->greater);
            traverse_no_checking(self, results, node1->greater, node2->greater);
        }
        else {
            traverse_no_checking(self, results, node1->less, node2);
            traverse_no_checking(self, results, node1->greater, node2);
        }
    }
}


template <typename MinMaxDist> static void
traverse_checking(const ckdtree *self,
                  std::vector<ordered_pair> *results,
                  const ckdtreenode *node1, const ckdtreenode *node2,
                  RectRectDistanceTracker<MinMaxDist> *tracker)
{
    const ckdtreenode *lnode1;
    const ckdtreenode *lnode2;
    double d;
    ckdtree_intp_t i, j, min_j;

    if (tracker->min_distance > tracker->upper_bound * tracker->epsfac)
        return;
    else if (tracker->max_distance < tracker->upper_bound / tracker->epsfac)
        traverse_no_checking(self, results, node1, node2);
    else if (node1->split_dim == -1) { /* 1 is leaf node */
        lnode1 = node1;

        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
            lnode2 = node2;

            /* brute-force */
            const double p = tracker->p;
            const double tub = tracker->upper_bound;
            const double *data = self->raw_data;
            const double epsfac = tracker->epsfac;
            const ckdtree_intp_t *indices = self->raw_indices;
            const ckdtree_intp_t m = self->m;
            const ckdtree_intp_t start1 = lnode1->start_idx;
            const ckdtree_intp_t start2 = lnode2->start_idx;
            const ckdtree_intp_t end1 = lnode1->end_idx;
            const ckdtree_intp_t end2 = lnode2->end_idx;

            CKDTREE_PREFETCH(data+indices[start1]*m, 0, m);
            if (start1 < end1 - 1)
               CKDTREE_PREFETCH(data+indices[start1+1]*m, 0, m);

            for(i = start1; i < end1; ++i) {

                if (i < end1 - 2)
                     CKDTREE_PREFETCH(data+indices[i+2]*m, 0, m);

                /* Special care here to avoid duplicate pairs */
                if (node1 == node2)
                    min_j = i + 1;
                else
                    min_j = start2;

                if (min_j < end2)
                    CKDTREE_PREFETCH(data+indices[min_j]*m, 0, m);
                if (min_j < end2 - 1)
                    CKDTREE_PREFETCH(data+indices[min_j+1]*m, 0, m);

                for (j = min_j; j < end2; ++j) {

                    if (j < end2 - 2)
                        CKDTREE_PREFETCH(data+indices[j+2]*m, 0, m);

                    d = MinMaxDist::point_point_p(
                            self,
                            data + indices[i] * m,
                            data + indices[j] * m,
                            p, m, tub);

                    if (d <= tub/epsfac)
                        add_ordered_pair(results, indices[i], indices[j]);
                }
            }
        }
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            traverse_checking(self, results, node1, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse_checking(self, results, node1, node2->greater, tracker);
            tracker->pop();
        }
    }
    else {  /* 1 is an inner node */
        if (node2->split_dim == -1) { /* 1 is an inner node, 2 is a leaf node */
            tracker->push_less_of(1, node1);
            traverse_checking(self, results, node1->less, node2, tracker);
            tracker->pop();

            tracker->push_greater_of(1, node1);
            traverse_checking(self, results, node1->greater, node2, tracker);
            tracker->pop();
        }
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse_checking(self, results, node1->less, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse_checking(self, results, node1->less, node2->greater,
                tracker);
            tracker->pop();
            tracker->pop();

            tracker->push_greater_of(1, node1);
            if (node1 != node2) {
                /*
                 * Avoid traversing (node1->less, node2->greater) and
                 * (node1->greater, node2->less) (it's the same node pair
                 * twice over, which is the source of the complication in
                 * the original KDTree.query_pairs)
                 */
                tracker->push_less_of(2, node2);
                traverse_checking(self, results, node1->greater, node2->less,
                    tracker);
                tracker->pop();
            }
            tracker->push_greater_of(2, node2);
            traverse_checking(self, results, node1->greater, node2->greater,
                tracker);
            tracker->pop();
            tracker->pop();
        }
    }
}


#include <iostream>

int
query_pairs(const ckdtree *self,
            const double r, const double p, const double eps,
            std::vector<ordered_pair> *results)
{

#define HANDLE(cond, kls) \
    if(cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, eps, r);\
        traverse_checking(self, results, self->ctree, self->ctree, \
            &tracker); \
    } else

    Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
    Rectangle r2(self->m, self->raw_mins, self->raw_maxes);

    if(CKDTREE_LIKELY(self->raw_boxsize_data == NULL)) {
        HANDLE(CKDTREE_LIKELY(p == 2), MinkowskiDistP2)
        HANDLE(p == 1, MinkowskiDistP1)
        HANDLE(std::isinf(p), MinkowskiDistPinf)
        HANDLE(1, MinkowskiDistPp)
        {}
    } else {
        HANDLE(CKDTREE_LIKELY(p == 2), BoxMinkowskiDistP2)
        HANDLE(p == 1, BoxMinkowskiDistP1)
        HANDLE(std::isinf(p), BoxMinkowskiDistPinf)
        HANDLE(1, BoxMinkowskiDistPp)
        {}
    }

    return 0;
}

