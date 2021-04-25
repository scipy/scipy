#include <cmath>
#include <cstdlib>
#include <cstring>

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>

#include "ckdtree_decl.h"
#include "rectangle.h"

struct WeightedTree {
    const ckdtree *tree;
    double *weights;
    double *node_weights;
};

struct CNBParams
{
    double *r;
    void * results; /* will be casted inside */
    WeightedTree self, other;
    int cumulative;
};

template <typename MinMaxDist, typename WeightType, typename ResultType> static void
traverse(
    RectRectDistanceTracker<MinMaxDist> *tracker,
    const CNBParams *params,
    double *start, double *end,
    const ckdtreenode *node1,
    const ckdtreenode *node2)
{
    static void (* const next)(RectRectDistanceTracker<MinMaxDist> *tracker,
            const CNBParams *params,
            double *start, double *end,
            const ckdtreenode *node1,
            const ckdtreenode *node2) = traverse<MinMaxDist, WeightType, ResultType>;

    ResultType *results = (ResultType*) params->results;

    /*
     * Speed through pairs of nodes all of whose children are close
     * and see if any work remains to be done
     */

    double * new_start = std::lower_bound(start, end, tracker->min_distance);
    double * new_end = std::lower_bound(start, end, tracker->max_distance);


    /* since max_distance >= min_distance, end < start never happens */
    if (params->cumulative) {
        double * i;
        if (new_end != end) {
            ResultType nn = WeightType::get_weight(&params->self, node1)
                          * WeightType::get_weight(&params->other, node2);

            for (i = new_end; i < end; ++i) {
                results[i - params->r] += nn;
            }
        }
        /* any bins larger than end have been correctly counted, thus
         * thus we can truncate the queries in future of this branch of the traversal*/
        start = new_start;
        end = new_end;
    } else {
        start = new_start;
        end = new_end;

        if (end == start) {
            ResultType nn = WeightType::get_weight(&params->self, node1)
                          * WeightType::get_weight(&params->other, node2);
            results[start - params->r] += nn;
        }
    }

    if (end == start) {
        /* this pair falls into exactly one bin, no need to probe deeper. */
        return;
    }

    /* OK, need to probe a bit deeper */
    if (node1->split_dim == -1) {  /* 1 is leaf node */
        if (node2->split_dim == -1) {  /* 1 & 2 are leaves */
            ckdtree_intp_t i, j;
            const double p = tracker->p;
            const double tmd = tracker->max_distance;
            const double *sdata = params->self.tree->raw_data;
            const ckdtree_intp_t *sindices = params->self.tree->raw_indices;
            const double *odata = params->other.tree->raw_data;
            const ckdtree_intp_t *oindices = params->other.tree->raw_indices;
            const ckdtree_intp_t m = params->self.tree->m;
            const ckdtree_intp_t start1 = node1->start_idx;
            const ckdtree_intp_t start2 = node2->start_idx;
            const ckdtree_intp_t end1 = node1->end_idx;
            const ckdtree_intp_t end2 = node2->end_idx;

            CKDTREE_PREFETCH(sdata + sindices[start1] * m, 0, m);

            if (start1 < end1 - 1)
                CKDTREE_PREFETCH(sdata + sindices[start1+1] * m, 0, m);

            /* brute-force */
            for (i = start1; i < end1; ++i) {

                if (i < end1 - 2)
                    CKDTREE_PREFETCH(sdata + sindices[i+2] * m, 0, m);

                CKDTREE_PREFETCH(odata + oindices[start2] * m, 0, m);

                if (start2 < end2 - 1)
                    CKDTREE_PREFETCH(odata + oindices[start2+1] * m, 0, m);

                for (j = start2; j < end2; ++j) {

                    if (j < end2 - 2)
                        CKDTREE_PREFETCH(odata + oindices[j+2] * m, 0, m);

                    double d = MinMaxDist::point_point_p(params->self.tree,
                            sdata + sindices[i] * m,
                            odata + oindices[j] * m,
                            p, m, tmd);

                    if (params->cumulative) {
                        /*
                         * I think it's usually cheaper to test d against all
                         * r's than to generate a distance array, sort it, then
                         * search for all r's via binary search
                         */
                        double * l;
                        for (l = start; l < end; ++l) {
                            if (d <= *l) {
                                results[l - params->r] += WeightType::get_weight(&params->self, sindices[i])
                                                        * WeightType::get_weight(&params->other, oindices[j]);
                            }
                        }
                    } else {
                        const double *l = std::lower_bound(start, end, d);
                        results[l - params->r] += WeightType::get_weight(&params->self, sindices[i])
                                                * WeightType::get_weight(&params->other, oindices[j]);
                    }
                }
            }
        }
        else {  /* 1 is a leaf node, 2 is inner node */
            tracker->push_less_of(2, node2);
            next(tracker, params, start, end, node1, node2->less);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            next(tracker, params, start, end, node1, node2->greater);
            tracker->pop();
        }
    }
    else { /* 1 is an inner node */
        if (node2->split_dim == -1) {
            /* 1 is an inner node, 2 is a leaf node */
            tracker->push_less_of(1, node1);
            next(tracker, params, start, end, node1->less, node2);
            tracker->pop();

            tracker->push_greater_of(1, node1);
            next(tracker, params, start, end, node1->greater, node2);
            tracker->pop();
        }
        else { /* 1 and 2 are inner nodes */
            tracker->push_less_of(1, node1);
            tracker->push_less_of(2, node2);
            next(tracker, params, start, end, node1->less, node2->less);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            next(tracker, params, start, end, node1->less, node2->greater);
            tracker->pop();
            tracker->pop();

            tracker->push_greater_of(1, node1);
            tracker->push_less_of(2, node2);
            next(tracker, params, start, end, node1->greater, node2->less);
            tracker->pop();

            tracker->push_greater_of(2, node2);
            next(tracker, params, start, end, node1->greater, node2->greater);
            tracker->pop();
            tracker->pop();
        }
    }
}

template <typename WeightType, typename ResultType> void
count_neighbors(struct CNBParams *params,
                ckdtree_intp_t n_queries, const double p)
{

    const ckdtree *self = params->self.tree;
    const ckdtree *other = params->other.tree;

#define HANDLE(cond, kls) \
    if (cond) { \
        RectRectDistanceTracker<kls> tracker(self, r1, r2, p, 0.0, 0.0);\
        traverse<kls, WeightType, ResultType>(&tracker, params, params->r, params->r+n_queries, \
                 self->ctree, other->ctree); \
    } else

    Rectangle r1(self->m, self->raw_mins, self->raw_maxes);
    Rectangle r2(other->m, other->raw_mins, other->raw_maxes);

    if (CKDTREE_LIKELY(self->raw_boxsize_data == NULL)) {
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
}

struct Unweighted {
    /* the interface for accessing weights of unweighted data. */
    static inline ckdtree_intp_t
    get_weight(const WeightedTree *wt, const ckdtreenode * node)
    {
        return node->children;
    }
    static inline ckdtree_intp_t
    get_weight(const WeightedTree *wt, const ckdtree_intp_t i)
    {
        return 1;
    }
};


int
count_neighbors_unweighted(const ckdtree *self, const ckdtree *other,
                ckdtree_intp_t n_queries, double *real_r, intptr_t *results,
                const double p, int cumulative) {

    CNBParams params = {0};

    params.r = real_r;
    params.results = (void*) results;
    params.self.tree = self;
    params.other.tree = other;
    params.cumulative = cumulative;

    count_neighbors<Unweighted, ckdtree_intp_t>(&params, n_queries, p);

    return 0;
}

struct Weighted {
    /* the interface for accessing weights of weighted data. */
    static inline double
    get_weight(const WeightedTree *wt, const ckdtreenode * node)
    {
        return (wt->weights != NULL)
           ? wt->node_weights[node - wt->tree->ctree]
           : node->children;
    }
    static inline double
    get_weight(const WeightedTree *wt, const ckdtree_intp_t i)
    {
        return (wt->weights != NULL)?wt->weights[i]:1;
    }
};

int
count_neighbors_weighted(const ckdtree *self, const ckdtree *other,
                double *self_weights, double *other_weights,
                double *self_node_weights, double *other_node_weights,
                ckdtree_intp_t n_queries, double *real_r, double *results,
                const double p, int cumulative)
{

    CNBParams params = {0};

    params.r = real_r;
    params.results = (void*) results;
    params.cumulative = cumulative;

    params.self.tree = self;
    params.other.tree = other;
    if (self_weights) {
        params.self.weights = self_weights;
        params.self.node_weights = self_node_weights;
    }
    if (other_weights) {
        params.other.weights = other_weights;
        params.other.node_weights = other_node_weights;
    }

    count_neighbors<Weighted, double>(&params, n_queries, p);

    return 0;
}

