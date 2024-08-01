#pragma once

#include "tuples.h"

namespace xsf {

template <typename T, size_t N>
T dot(const T (&x)[N], const T (&y)[N]) {
    T res = 0;
    for (size_t i = 0; i < N; ++i) {
        res += x[i] * y[i];
    }

    return res;
}

template <typename T, size_t K>
void dot(std::tuple<T (&)[K]> x, const std::tuple<T (&)[K]> y, std::tuple<T &> res) {
    std::get<0>(res) = dot(std::get<0>(x), std::get<0>(y));
}

template <typename T, size_t N>
void dot(std::tuple<T (&)[N], T (&)[N]> x, std::tuple<T (&)[N], T (&)[N]> y, std::tuple<T &, T &> res) {
    const auto &[x0, x1] = x;
    const auto &[y0, y1] = y;
    auto &[res0, res1] = res;

    res0 = 0;
    res1 = 0;
    for (size_t i = 0; i < N; ++i) {
        res0 += x0[i] * y0[i];
        res1 += x0[i] * y1[i] + x1[i] * y0[i];
    }
}

template <typename T, size_t N>
void dot(
    std::tuple<T (&)[N], T (&)[N], T (&)[N]> x, std::tuple<T (&)[N], T (&)[N], T (&)[N]> y,
    std::tuple<T &, T &, T &> res
) {
    const auto &[x0, x1, x2] = x;
    const auto &[y0, y1, y2] = y;
    auto &[res0, res1, res2] = res;

    res0 = 0;
    res1 = 0;
    res2 = 0;
    for (size_t i = 0; i < N; ++i) {
        res0 += x0[i] * y0[i];
        res1 += x0[i] * y1[i] + x1[i] * y0[i];
        res2 += x0[i] * y2[i] + T(2) * x1[i] * y1[i] + x2[i] * y0[i];
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
    tuples::for_each(res, [](auto &elem) { forward_recur_shift_left(elem); });
}

template <typename T, size_t K>
void forward_recur_rotate_left(T (&res)[K]) {
    T tmp = res[0];
    forward_recur_shift_left(res);
    res[K - 1] = tmp;
}

template <typename... T, size_t K>
void forward_recur_rotate_left(std::tuple<T (&)[K]...> &res) {
    tuples::for_each(res, [](auto &elem) { forward_recur_rotate_left(elem); });
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
template <typename InputIt, typename Recurrence, typename... T, ptrdiff_t K, typename Func>
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
            r(it, tuples::ref(coef));

            std::tuple<T...> tmp;
            dot(tuples::ref(coef), res, tuples::ref(tmp));

            forward_recur_shift_left(res);
            tuples::access(res, K - 1) = tmp;

            f(it, res);
            ++it;
        }
    }
}

template <typename InputIt, typename Recurrence, typename... T, ptrdiff_t K, typename Func>
void backward_recur(InputIt first, InputIt last, Recurrence r, std::tuple<T (&)[K]...> res, Func f) {
    InputIt it = first;
    while (std::abs(it - first) != K && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        --it;
    }

    if (std::abs(last - first) > K) {
        while (it != last) {
            std::tuple<T[K]...> coef;
            r(it, tuples::ref(coef));

            std::tuple<T...> tmp;
            dot(tuples::ref(coef), res, tuples::ref(tmp));

            forward_recur_shift_left(res);
            tuples::access(res, K - 1) = tmp;

            f(it, res);
            --it;
        }
    }
}

} // namespace xsf
