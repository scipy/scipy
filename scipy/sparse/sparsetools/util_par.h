#ifndef __SPTOOLS_UTIL_PAR_H__
#define __SPTOOLS_UTIL_PAR_H__

#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <numeric>

#include "poolstl.hpp"

using poolstl::iota_iter;
using poolstl::for_each_chunk;

/**
 * EXPERIMENT: If USE_POOL=1 then a single shared thread pool is used for all
 * parallel calls. It is started when the first operation large enough to justify
 * parallelism is executed.
 * If USE_POOL=0 then threads are started and joined for every operation (that is large
 * enough to justify parallelism).
 * The thread pool version requires maintaining that thread pool, but also compiles
 * to a larger size because each operation needs to be packaged. The pure-threads
 * version results in a much smaller library size and does not leave idle threads
 * when done, but it will churn through threads if many operations are being done.
 */
#define USE_POOL 0


class worker_config {
private:
    // The thread pool is started on first use with this many worker threads.
    int num_workers = 1;
//    int num_workers = std::min(16, std::thread::hardware_concurrency());
    std::atomic<bool> pool_started{false};

    std::shared_ptr<task_thread_pool::task_thread_pool> pool;
    std::once_flag flag;

    // Small operations do not amortize the cost of parallelization,
    // so they are executed sequentially. If the estimated number of operations
    // is above this value then the parallel version is used.
    int par_threshold = 65536;
public:

    /**
     * Return the worker pool. The first call starts the pool.
     */
    std::shared_ptr<task_thread_pool::task_thread_pool> get_pool() {
        std::call_once(flag, [&](){
            pool = std::make_shared<task_thread_pool::task_thread_pool>(num_workers);
            pool_started = true;
        });
        return pool;
    }

    /**
     * Return the number of worker threads.
     * If the thread pool has not been started yet then this returns the
     * number of threads that it would be started with.
     */
    int get_workers() {
        if (pool_started) {
            return get_pool()->get_num_threads();
        }
        return num_workers;
    }

    /**
     * Set the number of worker threads.
     * If the thead pool has not been started yet then when it is started it will
     * use this number of threads.
     */
    void set_workers(int workers) {
        if (pool_started) {
            get_pool()->set_num_threads((unsigned int)workers);
        }
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
     * A light operation is one that performs very little work at each memory location,
     * making it dominated by memory bandwidth rather than processing power. Such
     * operations need larger operands for parallelism to pay off.
     */
    template <class I>
    bool decide_if_par(I estimated_flops, bool light) {
        return par_threshold > 0
               && ((estimated_flops / (light ? 20 : 1)) >= par_threshold)
               && num_workers != 1;
    }

#if USE_POOL
    poolstl::execution::parallel_policy par_if(bool call_par) {
        if (call_par) {
            return poolstl::par_if(call_par, *get_pool());
        } else {
            return poolstl::par_if(false);
        }
    }

    template <class I>
    poolstl::execution::parallel_policy par_if_flops(I estimated_flops) {
        return par_if(decide_if_par(estimated_flops, false));
    }

    template <class I>
    poolstl::execution::parallel_policy par_if_mem(I estimated_flops) {
        return par_if(decide_if_par(estimated_flops, true));
    }
#else
    poolstl::execution::pure_threads_policy par_if(bool call_par) {
        return poolstl::execution::par_if_threads(call_par, num_workers);
    }

    template <class I>
    poolstl::execution::pure_threads_policy par_if_flops(I estimated_flops) {
        return par_if(decide_if_par(estimated_flops, false));
    }

    template <class I>
    poolstl::execution::pure_threads_policy par_if_mem(I estimated_flops) {
        return par_if(decide_if_par(estimated_flops, true));
    }
#endif
};

extern worker_config workers;


/**
 * Zips a key and value iterator.
 * Enables sorting both index and value arrays in-place with std::sort.
 */
template<class RandItA, class RandItB>
class kv_pair_iters {
public:
    using value_type = std::pair<typename std::iterator_traits<RandItA>::value_type,
                                 typename std::iterator_traits<RandItB>::value_type>;
    using difference_type = typename std::iterator_traits<RandItA>::difference_type;
    using pointer = value_type *;
    using iterator_category = std::random_access_iterator_tag;

    struct reference {
        reference(RandItA a, RandItB b) : a(a), b(b) {}

        reference &operator=(reference &&rhs) { *a = std::move(*rhs.a); *b = std::move(*rhs.b); return *this; }
        reference &operator=(const value_type &&rhs) { *a = std::move(rhs.first); *b = std::move(rhs.second); return *this; }
//        operator value_type() && { return std::make_pair(std::move(*a), std::move(*b)); }
        operator value_type() { return std::make_pair(*a, *b); }

        RandItA a;
        RandItB b;

        friend void swap(reference lhs, reference rhs) {
            std::iter_swap(lhs.a, rhs.a);
            std::iter_swap(lhs.b, rhs.b);
        }
    };

    kv_pair_iters() : a{}, b{} {}
    kv_pair_iters(RandItA a, RandItB b) : a(a), b(b) {}
    kv_pair_iters(const kv_pair_iters<RandItA, RandItB> &rhs) : a(rhs.a), b(rhs.b) {}

    kv_pair_iters<RandItA, RandItB> &operator=(const kv_pair_iters &rhs) { a = rhs.a; b = rhs.b; return *this; }

    reference operator*() const { return reference(a, b); }
    reference operator[](difference_type rhs) const { return ref_proxy(a + rhs, b + rhs); }
    // operator-> has no meaning in this application

    bool operator==(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a == rhs.a; }
    bool operator!=(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a != rhs.a; }
    bool operator<(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a < rhs.a; }
    bool operator>(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a > rhs.a; }
    bool operator<=(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a <= rhs.a; }
    bool operator>=(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a >= rhs.a; }

    kv_pair_iters<RandItA, RandItB> &operator+=(difference_type rhs) { a += rhs; b += rhs; return *this; }
    kv_pair_iters<RandItA, RandItB> &operator-=(difference_type rhs) { a -= rhs; b -= rhs; return *this; }

    kv_pair_iters<RandItA, RandItB> &operator++() { ++a; ++b; return *this; }
    kv_pair_iters<RandItA, RandItB> &operator--() { --a; --b; return *this; }
    kv_pair_iters<RandItA, RandItB> operator++(int) { kv_pair_iters<RandItA, RandItB> ret(*this); ++a; ++b; return ret; }
    kv_pair_iters<RandItA, RandItB> operator--(int) { kv_pair_iters<RandItA, RandItB> ret(*this); --a; --b; return ret; }

    difference_type operator-(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a - rhs.a; }
    kv_pair_iters<RandItA, RandItB> operator-(difference_type rhs) const { return kv_pair_iters(a - rhs, b - rhs); }
    kv_pair_iters<RandItA, RandItB> operator+(difference_type rhs) const { return kv_pair_iters(a + rhs, b + rhs); }

    friend inline kv_pair_iters<RandItA, RandItB> operator+(difference_type lhs, const kv_pair_iters<RandItA, RandItB> &rhs) {
        return kv_pair_iters(lhs + rhs.a, lhs + rhs.b);
    }

protected:
    RandItA a;
    RandItB b;
};

namespace std {
    /**
     * Specialize std::iterator_traits for kv_pair_iters.
     */
    template<class A, class B>
    struct iterator_traits<kv_pair_iters<A, B>> {
        using value_type =        typename kv_pair_iters<A, B>::value_type;
        using difference_type =   typename kv_pair_iters<A, B>::difference_type;
        using pointer =           typename kv_pair_iters<A, B>::pointer;
        using reference =         typename kv_pair_iters<A, B>::reference;
        using iterator_category = typename kv_pair_iters<A, B>::iterator_category;
    };
}

// GCC 8 does not have std::exclusive_scan. Reimplement if necessary.
template<typename InputIt, typename OutputIt, typename T>
OutputIt exclusive_scan(InputIt first, InputIt last, OutputIt result, T init)
{
#if (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 9)
    return std::exclusive_scan(first, last, result, init);
#else
    for (; first != last; ++first) {
        auto v = init;
        init = init + *first;
        *result++ = v;
    }
    return result;
#endif
}

// GCC 8 does not have std::inclusive_scan. Reimplement if necessary.
template<typename InputIt, typename OutputIt>
OutputIt inclusive_scan(InputIt first, InputIt last, OutputIt result)
{
#if (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 9)
    return std::inclusive_scan(first, last, result);
#else
    typename std::iterator_traits<OutputIt>::value_type init = 0;

    for (; first != last; ++first) {
        init = init + *first;
        *result++ = init;
    }
    return result;
#endif
}

#endif
