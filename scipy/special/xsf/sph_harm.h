#pragma once

#include "legendre.h"

namespace xsf {

template <typename T>
struct get_i;
using namespace std::complex_literals;

template <>
struct get_i<float> {

    static constexpr std::complex<float> value = 1if;
};

template <>
struct get_i<double> {

    static constexpr std::complex<double> value = 1i;
};

template <typename T, size_t N>
struct get_i<dual<T, N>> : get_i<T> {};

template <typename T>
void sph_harm_y_next(int m, T phi, T p, complex_type_t<T> &res) {
    using namespace std::complex_literals;

    complex_type_t<T> x2 = complex_type_t<T>(get_i<T>::value) * complex_type_t<T>(T(m) * phi);

    complex_type_t<T> z = exp(x2);
    res = complex_type_t<T>(p) * z;
}

template <typename T, typename Func>
void sph_harm_y_for_each_n(int n, int m, T theta, T phi, complex_type_t<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n(n, m, theta, p, [m, phi, &res, &f](int n, const T(&p)[2]) {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T, typename Func>
void sph_harm_y_for_each_n_m(int n, int m, T theta, T phi, complex_type_t<T> &res, Func f) {
    T p[2];
    sph_legendre_p_for_each_n_m(n, m, theta, p, [phi, &res, &f](int n, int m, const T(&p)[2]) {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

template <typename T>
complex_type_t<T> sph_harm_y(int n, int m, T theta, T phi) {
    complex_type_t<T> res_n;
    sph_harm_y_for_each_n(n, m, theta, phi, res_n, [](int n, int m, const complex_type_t<T> &res_n) {});

    return res_n;
}

template <typename T, typename... OutMats>
void sph_harm_y_all(T theta, T phi, std::tuple<OutMats...> res) {
    auto &res0 = std::get<0>(res);
    int n_max = res0.extent(0) - 1;
    int m_max = (res0.extent(1) - 1) / 2;

    complex_type_t<T> res_n_m;
    sph_harm_y_for_each_n_m(n_max, m_max, theta, phi, res_n_m, [m_max, res](int n, int m, complex_type_t<T> &res_n_m) {
        if (m >= 0) {
            dual_assign_grad(res_n_m, tuples::submdspan(res, n, m));
        } else {
            dual_assign_grad(res_n_m, tuples::submdspan(res, n, m + 2 * m_max + 1));
        }
    });
}

} // namespace xsf
