#pragma once

#include "legendre.h"
#include "numbers.h"

namespace xsf {
namespace detail {

    template <typename T>
    void sph_harm_y_next(int m, T phi, T p, complex<T> &res) {
        res = p * exp(numbers::i_v<T> * T(m) * phi);
    }

} // namespace detail

template <typename T, typename Func>
void sph_harm_y_for_each_n(int n, int m, T theta, T phi, complex<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n(n, m, theta, p, [m, phi, &res, &f](int n, const T(&p)[2]) {
        detail::sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T, typename Func>
void sph_harm_y_for_each_n_m(int n, int m, T theta, T phi, complex<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n_m(n, m, theta, p, [phi, &res, &f](int n, int m, const T(&p)[2]) {
        detail::sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T>
complex<T> sph_harm_y(int n, int m, T theta, T phi) {
    complex<T> res_n;
    sph_harm_y_for_each_n(n, m, theta, phi, res_n, [](int n, int m, const complex<T> &res_n) {});

    return res_n;
}

template <typename T, typename OutputMat>
void sph_harm_y_all(T theta, T phi, OutputMat res) {
    int n_max = res.extent(0) - 1;
    int m_max = (res.extent(1) - 1) / 2;

    complex<T> res_n_m;
    sph_harm_y_for_each_n_m(n_max, m_max, theta, phi, res_n_m, [m_max, &res](int n, int m, complex<T> &res_n_m) {
        if (m >= 0) {
            res(n, m) = res_n_m;
        } else {
            res(n, m + 2 * m_max + 1) = res_n_m;
        }
    });
}

} // namespace xsf
