#pragma once

#include "config.h"

namespace special {

template <typename T>
T binom(T n, T k) {
    if (n == 1) {
        if (k == 0) {
            return 1;
        }

        if (k == 1) {
            return 1;
        }
    }

    if (n == 2) {
        if (k == 0 || k == 2) {
            return 1;
        }

        return 2;
    }

    return std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
}

template <typename Recurrence, typename InputIt, typename T, ptrdiff_t KP1, ptrdiff_t NP1>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[KP1][NP1]) {
    constexpr ptrdiff_t K = KP1 - 1;
    constexpr ptrdiff_t N = NP1 - 1;

    T coef[K][NP1];
    r(it, coef);

    res[K][0] = 0;
    for (ptrdiff_t j = 0; j < K; ++j) {
        res[K][0] += coef[j][0] * res[j][0];
    }

    for (size_t n = 1; n <= N; ++n) {
        res[K][n] = 0;
        for (size_t k = 0; k <= n; ++k) {
            for (size_t i = 0; i < K; ++i) {
                res[K][n] += binom<remove_complex_t<T>>(n, k) * coef[i][n - k] * res[i][k];
            }
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
template <
    typename Recurrence, typename T, ptrdiff_t P, ptrdiff_t N, typename InputIt, typename Callback, typename... Args>
void forward_recur(Recurrence r, T (&res)[P][N], InputIt first, InputIt last, Callback &&callback, Args &&...args) {
    constexpr ptrdiff_t K = P - 1;

    InputIt it = first;
    while (it - first != K && it != last) {
        for (ptrdiff_t k = 0; k < N; ++k) {
            T tmp = res[K][k];
            res[K][k] = res[0][k];
            for (ptrdiff_t j = 0; j < K - 1; ++j) {
                res[j][k] = res[j + 1][k];
            }
            res[K - 1][k] = tmp;
        }

        callback(it, r, res, std::forward<Args>(args)...);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            for (ptrdiff_t j = 0; j < K; ++j) {
                for (ptrdiff_t k = 0; k < N; ++k) {
                    res[j][k] = res[j + 1][k];
                }
            }

            forward_recur_next(r, it, res);

            callback(it, r, res, std::forward<Args>(args)...);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t P, size_t N, typename InputIt>
void forward_recur(Recurrence r, T (&res)[P][N], InputIt first, InputIt last) {
    forward_recur(r, res, first, last, [](InputIt it, Recurrence r, const T(&res)[P][N]) {});
}

} // namespace special
