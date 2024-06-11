#pragma once

#include "legendre.h"

namespace special {

template <typename T, typename... OutputVals, typename Func>
void sph_harm_for_each_n(int n, int m, T theta, T phi, std::tuple<OutputVals (&)[2]...> res, Func f) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<T[2], N> p;
    sph_legendre_p_for_each_n(n, m, phi, tuple_ref_each(p), [m, res, f, theta](int n, auto &p_res) {
        T p = std::get<0>(p_res)[1];

        std::get<0>(res)[1] = p * std::exp(std::complex(T(0), m * theta));

        f(n, m, res);
    });
}

template <typename T, typename... OutputVals, typename Func>
void sph_harm_for_each_n_m(int n, int m, T theta, T phi, std::tuple<OutputVals (&)[2]...> res, Func f) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<T[2], N> p;
    sph_legendre_p_for_each_n_m(
        n, m, phi, tuple_ref_each(p),
        [theta, res, &f](int n, int m, grad_tuple_t<const T(&)[2], N> p_res) {
            T p = std::get<0>(p_res)[1];

            std::get<0>(res)[1] = p * std::exp(std::complex(T(0), m * theta));

            f(n, m, res);
        }
    );
}

template <typename T, typename... OutputVals>
void sph_harm(int n, int m, T theta, T phi, std::tuple<OutputVals &...> res) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<std::complex<T>[2], N> res_n;
    sph_harm_for_each_n(
        n, m, theta, phi, tuple_ref_each(res_n), [](int n, int m, grad_tuple_t<const std::complex<T>(&)[2], N> res_n) {}
    );

    res = tuple_access_each(res_n, 1);
}

template <typename T>
std::complex<T> sph_harm(int n, int m, T theta, T phi) {
    std::complex<T> res;
    sph_harm(n, m, theta, phi, std::tie(res));

    return res;
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    int m = (y.extent(0) - 1) / 2;
    int n = y.extent(1) - 1;

    std::tuple<std::complex<T>[2]> res;
    sph_harm_for_each_n_m(
        n, m, theta, phi, tuple_ref_each(res),
        [y](int n, int m, std::tuple<const std::complex<T>(&)[2]> res) {
            if (m >= 0) {
                y(m, n) = std::get<0>(res)[1];
            } else {
                y(m + y.extent(0), n) = std::get<0>(res)[1];
            }
        }
    );
}

} // namespace special
