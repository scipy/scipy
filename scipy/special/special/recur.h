#pragma once

#include "config.h"

namespace special {

template <typename T, size_t N>
void forward_recur_rotate(T (&res)[N]) {
    T tmp = res[N - 1];
    res[N - 1] = res[0];
    for (size_t j = 0; j < N - 1; ++j) {
        res[j] = res[j + 1];
    }
    res[N - 2] = tmp;
}

template <typename T, size_t N>
void forward_recur_shift(T (&res)[N]) {
    for (size_t j = 0; j < N - 1; ++j) {
        res[j] = res[j + 1];
    }
}

template <typename Recurrence, typename InputIt, typename T, size_t N>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[N]) {
    T coef[N - 1];
    r(it, coef);

    res[N - 1] = 0;
    for (ssize_t j = N - 2; j >= 0; --j) {
        res[N - 1] += coef[j] * res[j];
    }
}

template <typename Recurrence, typename InputIt, typename T, size_t N>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[N], T (&res_jac)[N]) {
    T coef[N - 1];
    T coef_jac[N - 1];
    r(it, coef, coef_jac);

    res[N - 1] = 0;
    res_jac[N - 1] = 0;
    for (ssize_t j = N - 2; j >= 0; --j) {
        res[N - 1] += coef[j] * res[j];
        res_jac[N - 1] += coef[j] * res_jac[j] + coef_jac[j] * res[j];
    }
}

template <typename Recurrence, typename InputIt, typename T, size_t N>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[N], T (&res_jac)[N], T (&res_hess)[N]) {
    T coef[N - 1];
    T coef_jac[N - 1];
    T coef_hess[N - 1];
    r(it, coef, coef_jac, coef_hess);

    res[N - 1] = 0;
    res_jac[N - 1] = 0;
    res_hess[N - 1] = 0;
    for (ssize_t j = N - 2; j >= 0; --j) {
        res[N - 1] += coef[j] * res[j];
        res_jac[N - 1] += coef[j] * res_jac[j] + coef_jac[j] * res[j];
        res_hess[N - 1] += coef[j] * res_hess[j] + T(2) * coef_jac[j] * res_jac[j] + coef_hess[j] * res[j];
    }
}

/**
 * Compute a forward recurrence that depends on K previous values.
 *
 * @param r recurrence
 * @param init initial value and its derivatives
 * @param res value and its derivatives
 * @param callback a function to be called as callback(i, r, args...) for 0 <= i <= n
 * @param args arguments to forward to the callback
 */
template <typename Recurrence, typename T, size_t N, typename InputIt, typename Callback>
void forward_recur(Recurrence r, T (&res)[N], InputIt first, InputIt last, Callback callback) {
    constexpr ptrdiff_t K = N - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate(res);

        callback(it, r, res);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            forward_recur_shift(res);
            forward_recur_next(r, it, res);

            callback(it, r, res);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t N, typename InputIt>
void forward_recur(Recurrence r, T (&res)[N], InputIt first, InputIt last) {
    forward_recur(r, res, first, last, [](InputIt it, Recurrence r, const T(&res)[N]) {});
}

/**
 * Compute a forward recurrence that depends on K previous values.
 *
 * @param r recurrence
 * @param init initial value and its derivatives
 * @param res value and its derivatives
 * @param callback a function to be called as callback(i, r, args...) for 0 <= i <= n
 * @param args arguments to forward to the callback
 */
template <typename Recurrence, typename T, size_t N, typename InputIt, typename Callback>
void forward_recur(Recurrence r, T (&res)[N], T (&res_jac)[N], InputIt first, InputIt last, Callback callback) {
    constexpr ptrdiff_t K = N - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate(res);
        forward_recur_rotate(res_jac);

        callback(it, r, res, res_jac);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            forward_recur_shift(res);
            forward_recur_shift(res_jac);
            forward_recur_next(r, it, res, res_jac);

            callback(it, r, res, res_jac);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t N, typename InputIt>
void forward_recur(Recurrence r, T (&res)[N], T (&res_jac)[N], InputIt first, InputIt last) {
    forward_recur(r, res, res_jac, first, last, [](InputIt it, Recurrence r, const T(&res)[N], const T(&res_jac)[N]) {
    });
}

/**
 * Compute a forward recurrence that depends on K previous values.
 *
 * @param r recurrence
 * @param init initial value and its derivatives
 * @param res value and its derivatives
 * @param callback a function to be called as callback(i, r, args...) for 0 <= i <= n
 */
template <typename Recurrence, typename T, size_t N, typename InputIt, typename Callback>
void forward_recur(
    Recurrence r, T (&res)[N], T (&res_jac)[N], T (&res_hess)[N], InputIt first, InputIt last, Callback callback
) {
    constexpr ptrdiff_t K = N - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate(res);
        forward_recur_rotate(res_jac);
        forward_recur_rotate(res_hess);

        callback(it, r, res, res_jac, res_hess);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            forward_recur_shift(res);
            forward_recur_shift(res_jac);
            forward_recur_shift(res_hess);
            forward_recur_next(r, it, res, res_jac, res_hess);

            callback(it, r, res, res_jac, res_hess);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t N, typename InputIt>
void forward_recur(Recurrence r, T (&res)[N], T (&res_jac)[N], T (&res_hess)[N], InputIt first, InputIt last) {
    forward_recur(
        r, res, res_jac, res_hess, first, last,
        [](InputIt it, Recurrence r, const T(&res)[N], const T(&res_jac)[N], const T(&res_hess)[N]) {}
    );
}

} // namespace special
