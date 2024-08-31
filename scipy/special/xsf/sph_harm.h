#pragma once

#include "legendre.h"

namespace xsf {

template <typename T>
void sph_harm_y_next(int m, T phi, T p, std::complex<T> &res) {
    std::complex<T> z = exp(std::complex(T(0), T(m) * phi));
    res = p * z;
}

template <typename T, typename Func>
void sph_harm_y_for_each_n(int n, int m, T theta, T phi, std::complex<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n(n, m, theta, p, [m, phi, &res, &f](int n, const T(&p)[2]) {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T, typename Func>
void sph_harm_y_for_each_n_m(int n, int m, T theta, T phi, std::complex<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n_m(n, m, theta, p, [phi, &res, &f](int n, int m, const T(&p)[2]) {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T>
std::complex<T> sph_harm_y(int n, int m, T theta, T phi) {
    std::complex<T> res_n;
    sph_harm_y_for_each_n(n, m, theta, phi, res_n, [](int n, int m, const std::complex<T> &res_n) {});

    return res_n;
}

template <typename T, typename... OutMats>
void sph_harm_y_all(T theta, T phi, std::tuple<OutMats...> res) {
    auto &res0 = std::get<0>(res);
    int n_max = res0.extent(0) - 1;
    int m_max = (res0.extent(1) - 1) / 2;

    std::complex<T> res_n_m;
    sph_harm_y_for_each_n_m(n_max, m_max, theta, phi, res_n_m, [m_max, res](int n, int m, std::complex<T> &res_n_m) {
        if (m >= 0) {
            dual_assign_grad(res_n_m, tuples::submdspan(res, n, m));
        } else {
            dual_assign_grad(res_n_m, tuples::submdspan(res, n, m + 2 * m_max + 1));
        }
    });
}

} // namespace xsf
