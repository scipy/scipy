#include <algorithm>
#include <cstddef>
#include <iterator>
#include <thread>

/*
 * Implement some C++17 parallel algorithms using threads.
 *
 * Instead of std::execution::par, use the util_par::threaded execution policy with
 * the specified thread count.
 *
 * For example, iterate over a vector `vec` in parallel:
 * std::for_each(util_par::threaded(4), vec.begin(), vec.end(), [](auto& element) {
 *    // loop body
 * });
 *
 * For additional algorithms see https://github.com/alugowski/poolSTL
 */
namespace util_par {
    /*
     * An execution policy that allows executing algorithms with pure threads,
     * that is starting then joining threads for each invocation.
     *
     * Modeled after std::execution::parallel_policy.
     */
    struct threaded {
        explicit threaded(unsigned int p): p(p) {}

        [[nodiscard]] unsigned int get_num_threads() const {
            if (p == 0) {
                return std::thread::hardware_concurrency();
            }
            return p;
        }

    protected:
        const unsigned int p;
    };

    namespace internal {
        /**
         * To enable/disable overload resolution of threaded policy.
         */
        template <class ExecPolicy, class Tp>
        using enable_if_threaded =
            typename std::enable_if<
                std::is_same<util_par::threaded,
                    typename std::remove_cv<typename std::remove_reference<ExecPolicy>::type>::type>::value,
                Tp>::type;

        /*
         * Chunk the range [first, last) according to the execution policy's
         * thread count, and spawn threads for each chunk. The last chunk
         * (or only chunk if single-threaded) is executed in the calling thread.
         */
        template <class RandIt, class Chunk, typename... A>
        void parallel_chunk_for(util_par::threaded &&policy,
                                RandIt first, RandIt last,
                                Chunk chunk, A&&... chunk_args) {
            std::vector<std::thread> threads;

            // Calculate number of iterations per chunk
            std::size_t num_steps = (std::size_t)std::distance(first, last);
            unsigned int p = policy.get_num_threads();
            std::size_t chunk_size = (num_steps / p) + ((num_steps % p) > 0 ? 1 : 0);

            while (first != last) {
                RandIt loop_end = first + std::min(chunk_size, (std::size_t)std::distance(first, last));

                if (loop_end == last) {
                    // execute the last (or only) chunk in the calling thread
                    chunk(first, loop_end, std::forward<A>(chunk_args)...);
                } else {
                    threads.emplace_back(std::thread(std::forward<Chunk>(chunk),
                                                     first, loop_end,
                                                     std::forward<A>(chunk_args)...));
                }

                first = loop_end;
            }

            for (auto& thread : threads) {
                if (thread.joinable()) {
                    thread.join();
                }
            }
        }
    }

    template <class RandIt, class ChunkConstructor, class UnaryFunc>
    void for_each_chunk(RandIt first, RandIt last, ChunkConstructor construct, UnaryFunc f) {
        if (first == last) {
            return;
        }

        auto chunk_data = construct();
        for (; first != last; ++first) {
            f(*first, chunk_data);
        }
    }

    /**
     * Like parallel `std::for_each`, but exposes the chunking.
     * The `construct` method is called once per parallel chunk and
     * its output is passed to `f`.
     *
     * Useful for cases where an expensive workspace can be shared between loop
     * iterations but cannot be shared by all parallel iterations.
     *
     * NOTE: Iterators are expected to be random access.
     */
    template <class ExecPolicy, class RandIt, class ChunkConstructor, class UnaryFunc>
    util_par::internal::enable_if_threaded<ExecPolicy, void>
    for_each_chunk(ExecPolicy&& policy, RandIt first, RandIt last, ChunkConstructor construct, UnaryFunc f) {
        util_par::internal::parallel_chunk_for(
                std::forward<ExecPolicy>(policy),
                first, last,
                for_each_chunk <RandIt, ChunkConstructor, UnaryFunc>,
                construct, f);
    }
}

namespace std {
    /**
     * A threaded version of std::for_each.
     * NOTE: Iterators are expected to be random access.
     */
    template <class ExecPolicy, class RandIt, class UnaryFunc>
    util_par::internal::enable_if_threaded<ExecPolicy, void>
    for_each(ExecPolicy &&policy, RandIt first, RandIt last, UnaryFunc f) {
        util_par::internal::parallel_chunk_for(
                std::forward<ExecPolicy>(policy),
                first, last,
                std::for_each<RandIt, UnaryFunc>, f);
    }
}

namespace util_par {
    /**
     * An iterator over the integers.
     *
     * Effectively a view on a fictional vector populated by std::iota, but without
     * materializing anything. Use to parallelize loops that are not over a container.
     *
     * @tparam T A type that acts as an integer.
     */
    template<typename T>
    class iota_iter {
    public:
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T *;
        using reference = T;
        using iterator_category = std::random_access_iterator_tag;

        iota_iter() : value{} {}
        explicit iota_iter(T rhs) : value(rhs) {}
        iota_iter(const iota_iter<T> &rhs) : value(rhs.value) {}

        iota_iter<T> &operator=(T rhs) { value = rhs; return *this; }
        iota_iter<T> &operator=(const iota_iter &rhs) { value = rhs.value; return *this; }

        reference operator*() const { return value; }
        reference operator[](difference_type rhs) const { return value + rhs; }
        // operator-> has no meaning in this application

        bool operator==(const iota_iter<T> &rhs) const { return value == rhs.value; }
        bool operator!=(const iota_iter<T> &rhs) const { return value != rhs.value; }
        bool operator<(const iota_iter<T> &rhs) const { return value < rhs.value; }
        bool operator>(const iota_iter<T> &rhs) const { return value > rhs.value; }
        bool operator<=(const iota_iter<T> &rhs) const { return value <= rhs.value; }
        bool operator>=(const iota_iter<T> &rhs) const { return value >= rhs.value; }

        iota_iter<T> &operator+=(difference_type rhs) { value += rhs; return *this; }
        iota_iter<T> &operator-=(difference_type rhs) { value -= rhs; return *this; }

        iota_iter<T> &operator++() { ++value; return *this; }
        iota_iter<T> &operator--() { --value; return *this; }
        iota_iter<T> operator++(int) { iota_iter<T> ret(value); ++value; return ret; }
        iota_iter<T> operator--(int) { iota_iter<T> ret(value); --value; return ret; }

        difference_type operator-(const iota_iter<T> &rhs) const { return value - rhs.value; }
        iota_iter<T> operator-(difference_type rhs) const { return iota_iter(value - rhs); }
        iota_iter<T> operator+(difference_type rhs) const { return iota_iter(value + rhs); }

        friend inline iota_iter<T> operator+(difference_type lhs, const iota_iter<T> &rhs) {
            return iota_iter(lhs + rhs.value);
        }

    protected:
        T value;
    };
}

namespace std {
    /**
     * Specialize std::iterator_traits for util_par::iota_iter.
     */
    template <typename T>
    struct iterator_traits<util_par::iota_iter<T>> {
        using value_type =        typename util_par::iota_iter<T>::value_type;
        using difference_type =   typename util_par::iota_iter<T>::difference_type;
        using pointer =           typename util_par::iota_iter<T>::pointer;
        using reference =         typename util_par::iota_iter<T>::reference;
        using iterator_category = typename util_par::iota_iter<T>::iterator_category;
    };
}
