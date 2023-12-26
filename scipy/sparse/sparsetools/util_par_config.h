#ifndef __SPTOOLS_UTIL_PAR_H__
#define __SPTOOLS_UTIL_PAR_H__

#include "util_par_algorithms.h"

using util_par::for_each_chunk;
using util_par::iota_iter;

class worker_config {
private:
    // Number of worker threads. 1 means sequential, 0 means number of CPUs.
    int num_workers = 1;
//    int num_workers = std::min(16, std::thread::hardware_concurrency());

    // Small operations are executed sequentially as they do not amortize
    // the cost of parallelization. If the estimated number of operations
    // is above this value then the parallel version is used.
    int par_threshold = 65536;
public:

    /**
     * Return the number of worker threads.
     */
    int get_workers() {
        return num_workers;
    }

    /**
     * Set the number of worker threads.
     */
    void set_workers(int workers) {
        num_workers = workers;
    }

    /**
     * Set the threshold that justifies parallel execution.
     */
    void set_par_threshold(int desired_threshold) {
        par_threshold = desired_threshold;
    }

    /**
     * Decide if the estimated number of operations justifies running in parallel.
     */
    template <class I>
    bool decide_if_par(I estimated_flops) {
        return par_threshold > 0
               && (estimated_flops >= par_threshold)
               && num_workers != 1;
    }

    util_par::threaded par_if(bool call_par) {
        return util_par::threaded(call_par ? num_workers : 1);
    }

    template <class I>
    util_par::threaded par_if_flops(I estimated_flops) {
        return par_if(decide_if_par(estimated_flops));
    }
};

extern worker_config workers;

#endif
