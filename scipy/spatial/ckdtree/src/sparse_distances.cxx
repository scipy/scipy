#include "ckdtree_decl.h"
#include "rectangle.h"
#include "coo_entries.h"

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

template <typename MinMaxDist> static void
traverse(const ckdtree *self, const ckdtree *other,
         std::vector<coo_entry> *results,
         const ckdtreenode *node1, const ckdtreenode *node2,
         RectRectDistanceTracker<MinMaxDist> *tracker)
{

    if (tracker->min_distance > tracker->upper_bound)
        return;
    
    else if (node1->split_dim == -1) {  /* 1 is leaf node */

        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
            /* brute-force */
            const double p = tracker->p;
            const double tub = tracker->upper_bound;
            const double *sdata = self->raw_data;
            const ckdtree_intp_t *sindices = self->raw_indices;
            const double *odata = other->raw_data;
            const ckdtree_intp_t *oindices = other->raw_indices;
            const ckdtree_intp_t m = self->m;
            const ckdtree_intp_t start1 = node1->start_idx;
            const ckdtree_intp_t start2 = node2->start_idx;
            const ckdtree_intp_t end1 = node1->end_idx;
            const ckdtree_intp_t end2 = node2->end_idx;


            for (ckdtree_intp_t i = start1; i < end1; ++i) {

                for (ckdtree_intp_t j = start2; j < end2; ++j) {

                    double d = MinMaxDist::point_point_p(
                            self,
                            sdata + sindices[i] * m,
                            odata + oindices[j] * m,
                            p, m, tub);

                    if (d <= tub) {
                        if (p == 2.0)
                            d = std::sqrt(d);
                        else if ((p != 1) && (!std::isinf(p)))
                            d = std::pow(d, 1. / p);

                        coo_entry e = {sindices[i], oindices[j], d};
                        results->push_back(e);
                    }
                }
            }
        }
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            traverse(self, other, results, node1, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse(self, other, results, node1, node2->greater, tracker);
            tracker->pop();
        }
    }
    else {  /* 1 is an inner node */
        if (node2->split_dim == -1) {
            /* 1 is an inner node, 2 is a leaf node*/
            tracker->push_less_of(1, node1);
            traverse(self, other, results, node1->less, node2, tracker);
            tracker->pop();

            tracker->push_greater_of(1, node1);
            traverse(self, other, results, node1->greater, node2, tracker);
            tracker->pop();
        }
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse(self, other, results, node1->less, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse(self, other, results, node1->less, node2->greater, tracker);
            tracker->pop();
            tracker->pop();

            tracker->push_greater_of(1, node1);
            tracker->push_less_of(2, node2);
            traverse(self, other, results, node1->greater, node2->less, tracker);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            traverse(self, other, results, node1->greater, node2->greater,
                tracker);
            tracker->pop();
            tracker->pop();
        }
    }
}


int
sparse_distance_matrix(const ckdtree *self, const ckdtree *other,
                       const double p,
                       const double max_distance,
                       std::vector<coo_entry> *results)
{
#define HANDLE(cond, kls) \
    if(cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, 0, max_distance);\
        traverse(self, other, results, self->ctree, other->ctree, &tracker); \
    } else

    Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
    Rectangle r2(other->m, other->raw_mins, other->raw_maxes);
    if(self->raw_boxsize_data == NULL) {
        HANDLE(p == 2, MinkowskiDistP2)
        HANDLE(p == 1, MinkowskiDistP1)
        HANDLE(std::isinf(p), MinkowskiDistPinf)
        HANDLE(1, MinkowskiDistPp)
        {}
    } else {
        HANDLE(p == 2, BoxMinkowskiDistP2)
        HANDLE(p == 1, BoxMinkowskiDistP1)
        HANDLE(std::isinf(p), BoxMinkowskiDistPinf)
        HANDLE(1, BoxMinkowskiDistPp)
        {}
    }

    return 0;
}
