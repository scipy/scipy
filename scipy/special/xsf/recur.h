#pragma once

#include "config.h"

namespace xsf {

template <typename T, size_t N>
T dot(const T (&x)[N], const T (&y)[N]) {
    T res = T(0);
    for (size_t i = 0; i < N; ++i) {
        res += x[i] * y[i];
    }

    return res;
}

template <typename T, size_t K>
void forward_recur_shift_left(T (&res)[K]) {
    for (size_t k = 1; k < K; ++k) {
        res[k - 1] = res[k];
    }
}

template <typename T, size_t K>
void forward_recur_rotate_left(T (&res)[K]) {
    T tmp = res[0];
    forward_recur_shift_left(res);
    res[K - 1] = tmp;
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
template <typename InputIt, typename Recurrence, typename T, ptrdiff_t K, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, T (&res)[K], Func f) {
    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            T coef[K];
            r(it, coef);

            T tmp = dot(coef, res);
            forward_recur_shift_left(res);
            res[K - 1] = tmp;

            f(it, res);
            ++it;
        }
    }
}

template <typename InputIt, typename Recurrence, typename T, ptrdiff_t K, typename Func>
void backward_recur(InputIt first, InputIt last, Recurrence r, T (&res)[K], Func f) {
    InputIt it = first;
    while (std::abs(it - first) != K && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        --it;
    }

    if (std::abs(last - first) > K) {
        while (it != last) {
            T coef[K];
            r(it, coef);

            T tmp = dot(coef, res);
            forward_recur_shift_left(res);
            res[K - 1] = tmp;

            f(it, res);
            --it;
        }
    }
}

} // namespace xsf
