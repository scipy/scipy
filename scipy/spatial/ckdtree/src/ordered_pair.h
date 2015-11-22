
#ifndef CKDTREE_ORDERED_PAIR
#define CKDTREE_ORDERED_PAIR

struct ordered_pair {
    npy_intp i;
    npy_intp j;
};

inline void 
add_ordered_pair(std::vector<ordered_pair> *results,
                       const npy_intp i, const npy_intp j)
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


