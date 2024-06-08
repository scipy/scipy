#pragma once

#include "config.h"

namespace special {

template <typename T, size_t N>
void forward_recur_shift_left(T (&res)[N]) {
    for (size_t n = 1; n < N; ++n) {
        res[n - 1] = res[n];
    }
}

template <typename T, size_t N>
void forward_recur_rotate_left(T (&res)[N]) {
    T tmp = res[0];
    forward_recur_shift_left(res);
    res[N - 1] = tmp;
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

template <typename InputIt, typename Recurrence, ssize_t N, typename T, typename Func>
void tuple_forward_recur(InputIt first, InputIt last, Recurrence r, grad_tuple<T[N], 1> &tuple_res, Func f) {
    T(&res)[N] = get<0>(tuple_res);
    T(&res_jac)[N] = get<1>(tuple_res);

    InputIt it = first;
    while (it - first != N && it != last) {
        forward_recur_rotate_left(res);
        forward_recur_rotate_left(res_jac);

        f(it, tuple_res);
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

            f(it, tuple_res);
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
