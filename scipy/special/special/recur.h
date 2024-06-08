#pragma once

#include "config.h"

namespace special {

template <typename T, size_t K>
void forward_recur_shift_left(T (&res)[K]) {
    for (size_t k = 1; k < K; ++k) {
        res[k - 1] = res[k];
    }
}

template <typename T, size_t K, size_t N>
void forward_recur_shift_left(grad<T[K], N> &res) {
    std::apply([](auto &...args) { (forward_recur_shift_left(args), ...); }, res.refs_as_tuple());
}

template <typename T, size_t K>
void forward_recur_rotate_left(T (&res)[K]) {
    T tmp = res[0];
    forward_recur_shift_left(res);
    res[K - 1] = tmp;
}

template <typename T, size_t K, size_t N>
void forward_recur_rotate_left(grad<T[K], N> &res) {
    std::apply([](auto &...args) { (forward_recur_rotate_left(args), ...); }, res.refs_as_tuple());
}

/**
 * Compute a forward recurrence that depends on N previous values.
 *
 * @param first begin iterator
 * @param last end iterator
 * @param r recurrence
 * @param res values, initialised to the leading N values
 * @param f a function to be called as f(it, res)
 */
template <typename InputIt, typename Recurrence, typename T, ssize_t N, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, T (&res)[N], Func f) {
    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        ++it;
    }

    if (last - first > N) {
        while (it != last) {
            T coef[N];
            r(it, coef);

            T res_next = 0;
            for (ssize_t n = 0; n < N; ++n) {
                res_next += coef[n] * res[n];
            }

            forward_recur_shift_left(res);
            res[N - 1] = res_next;

            f(it, res);
            ++it;
        }
    }
}

/**
 * Compute a forward recurrence that depends on N previous values.
 *
 * @param first begin iterator
 * @param last end iterator
 * @param r recurrence
 * @param res values, initialised to the leading N values
 * @param res_jac first derivatives, initialised to the leading N first derivatives
 * @param f a function to be called as f(it, res, res_jac)
 */
template <typename InputIt, typename Recurrence, typename T, ssize_t N, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, T (&res)[N], T (&res_jac)[N], Func f) {
    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate_left(res);
        forward_recur_rotate_left(res_jac);

        f(it, res, res_jac);
        ++it;
    }

    if (last - first > N) {
        while (it != last) {
            T coef[N];
            T coef_jac[N];
            r(it, coef, coef_jac);

            T res_next = 0;
            T res_next_jac = 0;
            for (ssize_t n = 0; n < N; ++n) {
                res_next += coef[n] * res[n];
                res_next_jac += coef[n] * res_jac[n] + coef_jac[n] * res[n];
            }

            forward_recur_shift_left(res);
            res[N - 1] = res_next;

            forward_recur_shift_left(res_jac);
            res_jac[N - 1] = res_next_jac;

            f(it, res, res_jac);
            ++it;
        }
    }
}

template <typename InputIt, typename Recurrence, typename T, ssize_t K, size_t N, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, grad<T[K], N> &res, Func f) {
    InputIt it = first;
    while (it - first != K && it != last) {
        forward_recur_rotate_left(res);

        f(it, res);
        ++it;
    }

    if (last - first > K) {
        while (it != last) {
            grad<T[K], N> coef;
            r(it, coef);

            grad<T, N> tmp;
            dot(coef, res, tmp);

            std::apply([](auto &...args) { (forward_recur_shift_left(args), ...); }, res.refs_as_tuple());
            std::apply([](auto &...args) { return std::tie(args[K - 1]...); }, res.refs_as_tuple()) =
                tmp.refs_as_tuple();

            f(it, res);
            ++it;
        }
    }
}

/**
 * Compute a forward recurrence that depends on N previous values.
 *
 * @param first begin iterator
 * @param last end iterator
 * @param r recurrence
 * @param res values, initialised to the leading N values
 * @param res_jac first derivatives, initialised to the leading N first derivatives
 * @param res_hess second derivatives, initialised to the leading N second derivative
 * @param f a function to be called as f(it, res, res_jac, res_hess)
 */
template <typename InputIt, typename Recurrence, typename T, ssize_t N, typename Func>
void forward_recur(InputIt first, InputIt last, Recurrence r, T (&res)[N], T (&res_jac)[N], T (&res_hess)[N], Func f) {
    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate_left(res);
        forward_recur_rotate_left(res_jac);
        forward_recur_rotate_left(res_hess);

        f(it, res, res_jac, res_hess);
        ++it;
    }

    if (last - first > N) {
        while (it != last) {
            T coef[N];
            T coef_jac[N];
            T coef_hess[N];
            r(it, coef, coef_jac, coef_hess);

            T res_next = 0;
            T res_next_jac = 0;
            T res_next_hess = 0;
            for (ssize_t n = 0; n < N; ++n) {
                res_next += coef[n] * res[n];
                res_next_jac += coef[n] * res_jac[n] + coef_jac[n] * res[n];
                res_next_hess += coef[n] * res_hess[n] + T(2) * coef_jac[n] * res_jac[n] + coef_hess[n] * res[n];
            }

            forward_recur_shift_left(res);
            res[N - 1] = res_next;

            forward_recur_shift_left(res_jac);
            res_jac[N - 1] = res_next_jac;

            forward_recur_shift_left(res_hess);
            res_hess[N - 1] = res_next_hess;

            f(it, res, res_jac, res_hess);
            ++it;
        }
    }
}

} // namespace special
