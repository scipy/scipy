#pragma once

#include "config.h"

namespace special {

template <typename Recurrence, typename InputIt, typename T, ptrdiff_t KP1>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[KP1]) {
    constexpr ptrdiff_t K = KP1 - 1;

    T coef[K];
    r(it, coef);

    res[K] = 0;
    for (ptrdiff_t j = 0; j < K; ++j) {
        res[K] += coef[j] * res[j];
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
template <typename Recurrence, typename T, ptrdiff_t P, typename InputIt, typename Callback, typename... Args>
void forward_recur(Recurrence r, T (&res)[P], InputIt first, InputIt last, Callback &&callback, Args &&...args) {
    constexpr ptrdiff_t K = P - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        T tmp = res[K];
        res[K] = res[0];
        for (ptrdiff_t j = 0; j < K - 1; ++j) {
            res[j] = res[j + 1];
        }
        res[K - 1] = tmp;

        callback(it, r, res, std::forward<Args>(args)...);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            for (ptrdiff_t j = 0; j < K; ++j) {
                res[j] = res[j + 1];
            }

            forward_recur_next(r, it, res);

            callback(it, r, res, std::forward<Args>(args)...);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t P, typename InputIt>
void forward_recur(Recurrence r, T (&res)[P], InputIt first, InputIt last) {
    forward_recur(r, res, first, last, [](InputIt it, Recurrence r, const T(&res)[P]) {});
}

template <typename Recurrence, typename InputIt, typename T, ptrdiff_t KP1>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[KP1], T (&res_jac)[KP1]) {
    constexpr ptrdiff_t K = KP1 - 1;

    T coef[K];
    T coef_jac[K];
    r(it, coef, coef_jac);

    res[K] = 0;
    for (ptrdiff_t j = 0; j < K; ++j) {
        res[K] += coef[j] * res[j];
    }

    res_jac[K] = 0;
    for (ptrdiff_t j = 0; j < K; ++j) {
        res_jac[K] += coef[j] * res_jac[j] + coef_jac[j] * res[j];
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
template <typename Recurrence, typename T, ptrdiff_t P, typename InputIt, typename Callback, typename... Args>
void forward_recur(
    Recurrence r, T (&res)[P], T (&res_jac)[P], InputIt first, InputIt last, Callback &&callback, Args &&...args
) {
    constexpr ptrdiff_t K = P - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        T tmp = res[K];
        res[K] = res[0];
        for (ptrdiff_t j = 0; j < K - 1; ++j) {
            res[j] = res[j + 1];
        }
        res[K - 1] = tmp;

        T tmp_jac = res_jac[K];
        res_jac[K] = res_jac[0];
        for (ptrdiff_t j = 0; j < K - 1; ++j) {
            res_jac[j] = res_jac[j + 1];
        }
        res_jac[K - 1] = tmp_jac;

        callback(it, r, res, res_jac, std::forward<Args>(args)...);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            for (ptrdiff_t j = 0; j < K; ++j) {
                res[j] = res[j + 1];
                res_jac[j] = res_jac[j + 1];
            }

            forward_recur_next(r, it, res, res_jac);

            callback(it, r, res, res_jac, std::forward<Args>(args)...);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t P, typename InputIt>
void forward_recur(Recurrence r, T (&res)[P], T (&res_jac)[P], InputIt first, InputIt last) {
    forward_recur(r, res, res_jac, first, last, [](InputIt it, Recurrence r, const T(&res)[P], const T(&res_jac)[P]) {
    });
}

template <typename Recurrence, typename InputIt, typename T, ptrdiff_t KP1>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[KP1], T (&res_jac)[KP1], T (&res_hess)[KP1]) {
    constexpr ptrdiff_t K = KP1 - 1;

    T coef[K];
    T coef_jac[K];
    T coef_hess[K];
    r(it, coef, coef_jac, coef_hess);

    res[K] = 0;
    res_jac[K] = 0;
    res_hess[K] = 0;
    for (ptrdiff_t j = 0; j < K; ++j) {
        res[K] += coef[j] * res[j];
        res_jac[K] += coef[j] * res_jac[j] + coef_jac[j] * res[j];
        res_hess[K] += coef[j] * res_hess[j] + T(2) * coef_jac[j] * res_jac[j] + coef_hess[j] * res[j];
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
template <typename Recurrence, typename T, ptrdiff_t P, typename InputIt, typename Callback, typename... Args>
void forward_recur(
    Recurrence r, T (&res)[P], T (&res_jac)[P], T (&res_hess)[P], InputIt first, InputIt last, Callback &&callback,
    Args &&...args
) {
    constexpr ptrdiff_t K = P - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        T tmp = res[K];
        res[K] = res[0];
        for (ptrdiff_t j = 0; j < K - 1; ++j) {
            res[j] = res[j + 1];
        }
        res[K - 1] = tmp;

        T tmp_jac = res_jac[K];
        res_jac[K] = res_jac[0];
        for (ptrdiff_t j = 0; j < K - 1; ++j) {
            res_jac[j] = res_jac[j + 1];
        }
        res_jac[K - 1] = tmp_jac;

        T tmp_hess = res_hess[K];
        res_hess[K] = res_hess[0];
        for (ptrdiff_t j = 0; j < K - 1; ++j) {
            res_hess[j] = res_hess[j + 1];
        }
        res_hess[K - 1] = tmp_hess;

        callback(it, r, res, res_jac, res_hess, std::forward<Args>(args)...);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            for (ptrdiff_t j = 0; j < K; ++j) {
                res[j] = res[j + 1];
                res_jac[j] = res_jac[j + 1];
                res_hess[j] = res_hess[j + 1];
            }

            forward_recur_next(r, it, res, res_jac, res_hess);

            callback(it, r, res, res_jac, res_hess, std::forward<Args>(args)...);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t P, typename InputIt>
void forward_recur(Recurrence r, T (&res)[P], T (&res_jac)[P], T (&res_hess)[P], InputIt first, InputIt last) {
    forward_recur(
        r, res, res_jac, res_hess, first, last,
        [](InputIt it, Recurrence r, const T(&res)[P], const T(&res_jac)[P], const T(&res_hess)[P]) {}
    );
}

} // namespace special
