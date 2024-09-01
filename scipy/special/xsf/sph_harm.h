#pragma once

#include "legendre.h"

namespace xsf {

using namespace std::complex_literals;

template <typename T>
std::complex<T> i_v;

template <>
std::complex<float> i_v<float> = 1if;

template <>
std::complex<double> i_v<double> = 1i;

template <typename T, size_t... N>
dual<std::complex<T>, N...> i_v<dual<T, N...>> = i_v<T>;

template <typename T, size_t N>
dual<dual<std::complex<T>, N>, N> i_v<dual<dual<T, N>, N>> = dual<dual<std::complex<T>, N>, N>(i_v<T>);

template <typename T>
void sph_harm_y_next(int m, T phi, T p, complex<T> &res) {
    res = complex<T>(p) * exp(i_v<T> * complex<T>(T(m) * phi));
}

template <typename T, typename Func>
void sph_harm_y_for_each_n(int n, int m, T theta, T phi, complex<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n(n, m, theta, p, [m, phi, &res, &f](int n, const T(&p)[2]) {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T, typename Func>
void sph_harm_y_for_each_n_m(int n, int m, T theta, T phi, complex<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n_m(n, m, theta, p, [phi, &res, &f](int n, int m, const T(&p)[2]) {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T>
complex<T> sph_harm_y(int n, int m, T theta, T phi) {
    complex<T> res_n;
    sph_harm_y_for_each_n(n, m, theta, phi, res_n, [](int n, int m, const complex<T> &res_n) {});

    return res_n;
}

template <typename T, typename... OutMats>
void sph_harm_y_all(T theta, T phi, std::tuple<OutMats...> res) {
    auto &res0 = std::get<0>(res);
    int n_max = res0.extent(0) - 1;
    int m_max = (res0.extent(1) - 1) / 2;

    complex<T> res_n_m;
    sph_harm_y_for_each_n_m(n_max, m_max, theta, phi, res_n_m, [m_max, res](int n, int m, complex<T> &res_n_m) {
        if (m >= 0) {
            dual_assign_grad(res_n_m, tuples::submdspan(res, n, m));
        } else {
            dual_assign_grad(res_n_m, tuples::submdspan(res, n, m + 2 * m_max + 1));
        }
    });
}

} // namespace xsf
