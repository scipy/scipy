
#ifndef CKDTREE_ORDERED_PAIR
#define CKDTREE_ORDERED_PAIR

#include <vector>

struct ordered_pair {
    ckdtree_intp_t i;
    ckdtree_intp_t j;
};

inline void
add_ordered_pair(std::vector<ordered_pair> *results,
                       const ckdtree_intp_t i, const intptr_t j)
{
    if (i > j) {
        ordered_pair p = {j,i};
        results->push_back(p);
    }
    else {
        ordered_pair p = {i,j};
        results->push_back(p);;
    }
}

#endif
