#pragma once

#include "config.h"

namespace special {

template <typename T, ssize_t N>
void forward_recur_rotate(T (&res)[N]) {
    T tmp = res[N - 1];
    res[N - 1] = res[0];
    for (ssize_t j = 0; j < N - 1; ++j) {
        res[j] = res[j + 1];
    }
    res[N - 2] = tmp;
}

template <typename T, ssize_t N>
void forward_recur_shift(T (&res)[N]) {
    for (ssize_t j = 0; j < N - 1; ++j) {
        res[j] = res[j + 1];
    }
}

template <typename Recurrence, typename InputIt, typename T, ssize_t N>
void forward_recur_next2(Recurrence r, InputIt it, T (&res)[N]) {
    T coef[N];
    r(it, coef);

    T tmp = 0;
    for (ssize_t j = N - 1; j >= 0; --j) {
        tmp += coef[j] * res[j];
    }

    forward_recur_shift(res);

    res[N - 1] = tmp;
}

template <typename Recurrence, typename InputIt, typename T, ssize_t N>
void forward_recur_next2(Recurrence r, InputIt it, T (&res)[N], T (&res_jac)[N]) {
    T coef[N];
    T coef_jac[N];
    r(it, coef, coef_jac);

    T tmp = 0;
    T tmp_jac = 0;
    for (ssize_t j = N - 1; j >= 0; --j) {
        tmp += coef[j] * res[j];
        tmp_jac += coef[j] * res_jac[j] + coef_jac[j] * res[j];
    }

    forward_recur_shift(res);
    forward_recur_shift(res_jac);

    res[N - 1] = tmp;
    res_jac[N - 1] = tmp_jac;
}

template <typename Recurrence, typename InputIt, typename T, ssize_t N>
void forward_recur_next2(Recurrence r, InputIt it, T (&res)[N], T (&res_jac)[N], T (&res_hess)[N]) {
    T coef[N];
    T coef_jac[N];
    T coef_hess[N];
    r(it, coef, coef_jac, coef_hess);

    T tmp = 0;
    T tmp_jac = 0;
    T tmp_hess = 0;
    for (ssize_t j = N - 1; j >= 0; --j) {
        tmp += coef[j] * res[j];
        tmp_jac += coef[j] * res_jac[j] + coef_jac[j] * res[j];
        tmp_hess += coef[j] * res_hess[j] + T(2) * coef_jac[j] * res_jac[j] + coef_hess[j] * res[j];
    }

    forward_recur_shift(res);
    forward_recur_shift(res_jac);
    forward_recur_shift(res_hess);

    res[N - 1] = tmp;
    res_jac[N - 1] = tmp_jac;
    res_hess[N - 1] = tmp_hess;
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
template <typename InputIt, typename Recurrence, typename T, ssize_t N, typename Callback>
void forward_recur2(InputIt first, InputIt last, Recurrence r, T (&res)[N], Callback callback) {
    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate(res);

        callback(it, res);
        ++it;
    }

    if (last - first > N) {
        while (it != last) {
            forward_recur_next2(r, it, res);

            callback(it, res);
            ++it;
        }
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
template <typename InputIt, typename Recurrence, typename T, ssize_t N, typename Callback>
void forward_recur2(InputIt first, InputIt last, Recurrence r, T (&res)[N], T (&res_jac)[N], Callback callback) {
    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate(res);
        forward_recur_rotate(res_jac);

        callback(it, res, res_jac);
        ++it;
    }

    if (last - first > N) {
        while (it != last) {
            forward_recur_next2(r, it, res, res_jac);

            callback(it, res, res_jac);
            ++it;
        }
    }
}

/**
 * Compute a forward recurrence that depends on K previous values.
 *
 * @param r recurrence
 * @param init initial value and its derivatives
 * @param res value and its derivatives
 * @param callback a function to be called as callback(i, r, args...) for 0 <= i <= n
 */
template <typename InputIt, typename Recurrence, typename T, ssize_t N, typename Callback>
void forward_recur2(
    InputIt first, InputIt last, Recurrence r, T (&res)[N], T (&res_jac)[N], T (&res_hess)[N], Callback callback
) {
    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate(res);
        forward_recur_rotate(res_jac);
        forward_recur_rotate(res_hess);

        callback(it, res, res_jac, res_hess);
        ++it;
    }

    if (last - first > N) {
        while (it != last) {
            forward_recur_next2(r, it, res, res_jac, res_hess);

            callback(it, res, res_jac, res_hess);
            ++it;
        }
    }
}

} // namespace special
