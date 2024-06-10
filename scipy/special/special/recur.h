#pragma once

#include "config.h"
#include "tuple_utility.h"

namespace special {

template <typename T, size_t K>
void dot(std::tuple<T (&)[K]> x, const std::tuple<T (&)[K]> y, std::tuple<T &> res) {
    const auto &[x0] = x;
    const auto &[y0] = y;
    auto &[res0] = res;

    res0 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(std::tuple<T (&)[K], T (&)[K]> x, std::tuple<T (&)[K], T (&)[K]> y, std::tuple<T &, T &> res) {
    const auto &[x0, x1] = x;
    const auto &[y0, y1] = y;
    auto &[res0, res1] = res;

    res0 = 0;
    res1 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
        res1 += x0[k] * y1[k] + x1[k] * y0[k];
    }
}

template <typename T, size_t K>
void dot(
    std::tuple<T (&)[K], T (&)[K], T (&)[K]> x, std::tuple<T (&)[K], T (&)[K], T (&)[K]> y,
    std::tuple<T &, T &, T &> res
) {
    const auto &[x0, x1, x2] = x;
    const auto &[y0, y1, y2] = y;
    auto &[res0, res1, res2] = res;

    res0 = 0;
    res1 = 0;
    res2 = 0;
    for (size_t k = 0; k < K; ++k) {
        res0 += x0[k] * y0[k];
        res1 += x0[k] * y1[k] + x1[k] * y0[k];
        res2 += x0[k] * y2[k] + T(2) * x1[k] * y1[k] + x2[k] * y0[k];
    }
}

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
