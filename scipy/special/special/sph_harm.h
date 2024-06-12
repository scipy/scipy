#pragma once

#include "legendre.h"

namespace special {

template <typename T>
void sph_harm_y_next(int m, T theta, std::tuple<const T (&)[2]> p, std::tuple<std::complex<T> (&)[2]> res) {
    const auto &[p0] = p;
    auto &[res0] = res;

    res0[1] = p0[1] * std::exp(std::complex(T(0), m * theta));
}

template <typename T, typename... OutputVals, typename Func>
void sph_harm_y_for_each_n(int n, int m, T theta, T phi, std::tuple<OutputVals (&)[2]...> res, Func f) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<T[2], N> p;
    sph_legendre_p_for_each_n(n, m, phi, tuples::ref(p), [m, theta, res, &f](int n, grad_tuple_t<const T(&)[2], N> p) {
        sph_harm_y_next(m, theta, p, res);

        f(n, m, res);
    });
}

template <typename T, typename... OutputVals, typename Func>
void sph_harm_y_for_each_n_m(int n, int m, T theta, T phi, std::tuple<OutputVals (&)[2]...> res, Func f) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<T[2], N> p;
    sph_legendre_p_for_each_n_m(
        n, m, phi, tuples::ref(p),
        [theta, res, &f](int n, int m, grad_tuple_t<const T(&)[2], N> p) {
            sph_harm_y_next(m, theta, p, res);

            f(n, m, res);
        }
    );
}

template <typename T, typename... OutputVals>
void sph_harm_y(int n, int m, T theta, T phi, std::tuple<OutputVals &...> res) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<std::complex<T>[2], N> res_n;
    sph_harm_y_for_each_n(
        n, m, theta, phi, tuples::ref(res_n), [](int n, int m, grad_tuple_t<const std::complex<T>(&)[2], N> res_n) {}
    );

    res = tuples::access(res_n, 1);
}

template <typename T>
std::complex<T> sph_harm_y(int n, int m, T theta, T phi) {
    std::complex<T> res;
    sph_harm_y(n, m, theta, phi, std::tie(res));

    return res;
}

template <typename T, typename... OutMats>
void sph_harm_y_all(T theta, T phi, std::tuple<OutMats...> res) {
    static constexpr size_t N = sizeof...(OutMats) - 1;

    auto &res0 = std::get<0>(res);
    int m_max = (res0.extent(0) - 1) / 2;
    int n_max = res0.extent(1) - 1;

    grad_tuple_t<std::complex<T>[2], N> res_n_m;
    sph_harm_y_for_each_n_m(
        n_max, m_max, theta, phi, tuples::ref(res_n_m),
        [m_max, res](int n, int m, grad_tuple_t<std::complex<T>(&)[2], N> res_n_m) {
            if (m >= 0) {
                tuples::call(res, m, n) = tuples::access(res_n_m, 1);
            } else {
                tuples::call(res, m + 2 * m_max + 1, n) = tuples::access(res_n_m, 1);
            }
        }
    );
}

} // namespace special
