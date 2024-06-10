#pragma once

#include "config.h"
#include "tuple_utility.h"

namespace special {

template <typename T, size_t K>
void forward_recur_shift_left(T (&res)[K]) {
    for (size_t k = 1; k < K; ++k) {
        res[k - 1] = res[k];
    }
}

template <typename... T, size_t K>
void forward_recur_shift_left(std::tuple<T (&)[K]...> &res) {
    tuple_for_each(res, [](auto &elem) { forward_recur_shift_left(elem); });
}

template <typename T, size_t K>
void forward_recur_rotate_left(T (&res)[K]) {
    T tmp = res[0];
    forward_recur_shift_left(res);
    res[K - 1] = tmp;
}

template <typename... T, size_t K>
void forward_recur_rotate_left(std::tuple<T (&)[K]...> &res) {
    tuple_for_each(res, [](auto &elem) { forward_recur_rotate_left(elem); });
}

/**
 * Compute a forward recurrence that depends on K previous values.
 *
 * @param first begin iterator
 * @param last end iterator
 * @param r recurrence
 * @param res values, initialised to the leading K values
 * @param f a function to be called as f(it, res)
 */
template <typename InputIt, typename Recurrence, typename... T, ssize_t K, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, std::tuple<T (&)[K]...> res, Func f) {
    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            std::tuple<T[K]...> coef;
            r(it, tuple_ref_each(coef));

            std::tuple<T...> tmp;
            dot(tuple_ref_each(coef), res, tuple_ref_each(tmp));

            forward_recur_shift_left(res);
            tuple_access_each(res, K - 1) = tmp;

            f(it, res);
            ++it;
        }
    }
}

} // namespace special
