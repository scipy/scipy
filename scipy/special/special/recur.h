#pragma once

#include "config.h"

namespace special {

template <typename T>
T binom(T n, T k) {
    return std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
}

template <typename Recurrence, typename T, size_t K, typename Callback, typename... Args>
T forward_recurrence(Recurrence r, const T (&init)[K], int n, Callback &&callback, Args &&...args) {
    T vals[K + 1] = {0};
    vals[0] = init[1];

    callback(0, r, vals, std::forward<Args>(args)...);

    if (n > 0) {
        vals[1] = vals[0];
        vals[0] = init[0];
        callback(1, r, vals, std::forward<Args>(args)...);

        for (int j = K; j <= n; ++j) {
            T coef[K];
            r(j, coef);

            for (size_t i = K; i > 0; --i) {
                vals[i] = vals[i - 1];
            }

            vals[0] = 0;
            for (size_t i = K; i >= 1; --i) {
                vals[0] += coef[i - 1] * vals[i];
            }

            callback(j, r, vals, std::forward<Args>(args)...);
        }
    }

    return vals[0];
}

template <typename Recurrence, typename T, size_t K>
T forward_recurrence(Recurrence r, const T (&init)[K], int n) {
    return forward_recurrence(r, init, n, [](int j, Recurrence r, const T(&p)[K + 1]) {});
}

template <typename Recurrence, typename InputIt, typename T, size_t P, size_t N>
void forward_recur_next(Recurrence r, InputIt it, T (&res)[P][N]) {
    constexpr int K = P - 1;

    T coef[K][N];
    r(it, coef);

    res[0][0] = 0;
    for (size_t j = K; j > 0; --j) {
        res[0][0] += coef[j - 1][0] * res[j][0];
    }

    for (size_t r2 = 1; r2 < N; ++r2) {
        res[0][r2] = 0;
        for (size_t k = 0; k <= r2; ++k) {
            res[0][r2] +=
                binom<remove_complex_t<T>>(r2, k) * (coef[0][r2 - k] * res[1][k] + coef[1][r2 - k] * res[2][k]);
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
    typename Recurrence, typename T, ptrdiff_t K, ptrdiff_t N, typename InputIt, typename Callback, typename... Args>
void forward_recur(
    Recurrence r, const T (&init)[K][N], T (&res)[K + 1][N], InputIt first, InputIt last, Callback &&callback,
    Args &&...args
) {
    InputIt it = first;
    while (it - first != K && it != last) {
        for (ptrdiff_t j = it - first; j > 0; --j) {
            for (ptrdiff_t k = 0; k < N; ++k) {
                res[j][k] = res[j - 1][k];
            }
        }

        for (ptrdiff_t k = 0; k < N; ++k) {
            res[0][k] = init[it - first][k];
        }

        callback(it, r, res, std::forward<Args>(args)...);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            for (ptrdiff_t j = K; j > 0; --j) {
                for (ptrdiff_t k = 0; k < N; ++k) {
                    res[j][k] = res[j - 1][k];
                }
            }

            forward_recur_next(r, it, res);

            callback(it, r, res, std::forward<Args>(args)...);
            ++it;
        }
    }
}

template <typename Recurrence, typename T, size_t K, size_t N, typename InputIt>
void forward_recur(Recurrence r, const T (&init)[K][N], T (&res)[K + 1][N], InputIt first, InputIt last) {
    forward_recur(r, init, res, first, last, [](InputIt it, Recurrence r, const T(&res)[K + 1][N]) {});
}

} // namespace special
