#pragma once

#include "legendre.h"

namespace special {

template <typename T>
void sph_harm_y_next(int m, T phi, std::tuple<const T (&)[2]> p, std::tuple<std::complex<T> &> res) {
    const auto &[p0] = p;

    std::complex<T> z = std::exp(std::complex(T(0), m * phi));
    tuples::assign(res, {p0[1] * z});
}

template <typename T>
void sph_harm_y_next(
    int m, T phi, std::tuple<const T (&)[2], const T (&)[2]> p,
    std::tuple<std::complex<T> &, std::complex<T> (&)[2]> res
) {
    const auto &[p0, p1] = p;

    std::complex<T> z = std::exp(std::complex(T(0), m * phi));
    tuples::assign(res, {p0[1] * z, {p1[1] * z, std::complex(T(0), T(m)) * p0[1] * z}});
}

template <typename T>
void sph_harm_y_next(
    int m, T phi, std::tuple<const T (&)[2], const T (&)[2], const T (&)[2]> p,
    std::tuple<std::complex<T> &, std::complex<T> (&)[2], std::complex<T> (&)[2][2]> res
) {
    const auto &[p0, p1, p2] = p;

    std::complex<T> z = std::exp(std::complex(T(0), m * phi));
    tuples::assign(
        res, {p0[1] * z,
              {p1[1] * z, std::complex(T(0), T(m)) * p0[1] * z},
              {{p2[1] * z, std::complex(T(0), T(m)) * p1[1] * z},
               {std::complex(T(0), T(m)) * p1[1] * z, -T(m * m) * p0[1] * z}}}
    );
}

template <typename T, typename... OutputVals, typename Func>
void sph_harm_y_for_each_n(int n, int m, T theta, T phi, std::tuple<OutputVals &...> res, Func f) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<T[2], N> p;
    sph_legendre_p_for_each_n(n, m, theta, tuples::ref(p), [m, phi, res, &f](int n, grad_tuple_t<const T(&)[2], N> p) {
        sph_harm_y_next(m, phi, p, res);

        f(n, m, res);
    });
}

template <typename T, typename... OutputVals, typename Func>
void sph_harm_y_for_each_n_m(int n, int m, T theta, T phi, std::tuple<OutputVals &...> res, Func f) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<T[2], N> p;
    sph_legendre_p_for_each_n_m(
        n, m, theta, tuples::ref(p),
        [phi, res, &f](int n, int m, grad_tuple_t<const T(&)[2], N> p) {
            sph_harm_y_next(m, phi, p, res);

            f(n, m, res);
        }
    );
}

template <typename T, typename... OutputVals>
void sph_harm_y(int n, int m, T theta, T phi, std::tuple<OutputVals &...> res) {
    static constexpr size_t N = sizeof...(OutputVals) - 1;

    grad_tuple_t<std::complex<T>, N, 2> res_n;
    sph_harm_y_for_each_n(
        n, m, theta, phi, tuples::ref(res_n), [](int n, int m, grad_tuple_t<std::complex<T> &, N, 2> res_n) {}
    );

    tuples::assign(res, res_n);
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
    int n_max = res0.extent(0) - 1;
    int m_max = (res0.extent(1) - 1) / 2;

    grad_tuple_t<std::complex<T>, N, 2> res_n_m;
    sph_harm_y_for_each_n_m(
        n_max, m_max, theta, phi, tuples::ref(res_n_m),
        [m_max, res](int n, int m, grad_tuple_t<std::complex<T> &, N, 2> res_n_m) {
            if (m >= 0) {
                tuples::assign(tuples::submdspan(res, n, m), res_n_m);
            } else {
                tuples::assign(tuples::submdspan(res, n, m + 2 * m_max + 1), res_n_m);
            }
        }
    );
}

} // namespace special
