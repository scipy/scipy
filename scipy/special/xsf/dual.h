#pragma once

#include <tuple>

#include <iostream>

#include "config.h"
#include "tuples.h"

namespace xsf {

inline int factorial(int n) {
    if (n == 0) {
        return 1;
    }

    return n * factorial(n - 1);
}

inline double binom2(double n, double m) { return factorial(n) / (factorial(m) * factorial(n - m)); }

template <typename T, size_t N0, size_t... N>
struct dual;

template <typename T, size_t N0, size_t... N>
struct dual_element {
    using type = dual<T, N...>;
};

template <typename T, size_t N0>
struct dual_element<T, N0> {
    using type = T;
};

template <typename T, size_t... N>
using dual_element_t = typename dual_element<T, N...>::type;

template <typename T, size_t N>
const T &dual_get(const dual<T, N> &x) {
    return x.data[0];
}

template <typename T, size_t... N>
const T &dual_get(const dual<T, N...> &x) {
    return dual_get(x.data[0]);
}

template <typename T, size_t N>
T &dual_get(dual<T, N> &x) {
    return x.data[0];
}

template <typename T, size_t... N>
T &dual_get(dual<T, N...> &x) {
    return dual_get(x.data[0]);
}

template <typename T, size_t N>
struct dual<T, N> {
  public:
    using value_type = T;
    using reference = value_type &;
    using const_reference = const value_type &;

  public:
    T data[N + 1];

  public:
    dual() = default;

    dual(T value) {
        data[0] = value;
        for (size_t i = 1; i <= N; ++i) {
            data[i] = 0;
        }
    }

    dual(std::initializer_list<T> data) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            this->data[it - data.begin()] = *it;
        }
        for (size_t i = data.size(); i <= N; ++i) {
            this->data[i] = 0;
        }
    }

    template <typename U>
    dual(const dual<U, N> &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] = other.data[i];
        }
    }

    dual &operator=(const T &other) {
        operator[](0) = other;
        for (size_t i = 1; i <= N; ++i) {
            operator[](i) = 0;
        }

        return *this;
    }

    dual &operator+=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] += other.data[i];
        }

        return *this;
    }

    dual &operator+=(const T &other) {
        data[0] += other;

        return *this;
    }

    dual &operator-=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] -= other.data[i];
        }

        return *this;
    }

    dual &operator-=(const T &other) {
        data[0] -= other;

        return *this;
    }

    dual &operator*=(const dual &other) {
        for (size_t i = N + 1; i-- > 0;) {
            data[i] *= other.data[0];
            for (size_t j = 0; j < i; ++j) {
                data[i] += T(binom2(i, i - j)) * data[j] * other.data[i - j];
            }
        }

        return *this;
    }

    dual &operator*=(const T &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] *= other;
        }

        return *this;
    }

    dual &operator/=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            for (size_t j = 1; j <= i; ++j) {
                data[i] -= T(binom2(i, i - j)) * other.data[j] * data[i - j];
            }
            data[i] /= other.data[0];
        }

        return *this;
    }

    dual &operator/=(const T &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] /= other;
        }

        return *this;
    }

    value_type &value() { return data[0]; }

    const value_type &value() const { return data[0]; }

    value_type &operator[](size_t i) { return data[i]; }

    const value_type &operator[](size_t i) const { return data[i]; }

    static constexpr size_t max_order() { return N; }

    static constexpr size_t max_order2 = N;
};

template <typename T, size_t N0, size_t N1, size_t... N>
struct dual<T, N0, N1, N...> {
  public:
    using value_type = T;
    using element_type = dual_element_t<T, N0, N1, N...>;
    using reference = value_type &;
    using const_reference = const value_type &;

  public:
    dual<T, N1, N...> data[N0 + 1];

  public:
    dual() = default;

    dual(T value) {
        data[0] = value;
        for (size_t i = 1; i <= N0; ++i) {
            data[i] = 0;
        }
    }

    dual(std::initializer_list<element_type> data) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            this->data[it - data.begin()] = *it;
        }
        for (size_t i = data.size(); i <= N0; ++i) {
            this->data[i] = 0;
        }
    }

    template <typename U>
    dual(const dual<U, N0, N1, N...> &other) {
        for (size_t i = 0; i <= N0; ++i) {
            data[i] = other.data[i];
        }
    }

    dual &operator=(const T &other) {
        operator[](0) = other;
        for (size_t i = 1; i <= N0; ++i) {
            operator[](i) = 0;
        }

        return *this;
    }

    dual &operator+=(const dual &other) {
        for (size_t i = 0; i <= N0; ++i) {
            data[i] += other.data[i];
        }

        return *this;
    }

    dual &operator+=(const T &other) {
        data[0] += other;

        return *this;
    }

    dual &operator-=(const dual &other) {
        for (size_t i = 0; i <= N0; ++i) {
            data[i] -= other.data[i];
        }

        return *this;
    }

    dual &operator-=(const T &other) {
        data[0] -= other;

        return *this;
    }

    dual &operator*=(const dual &other) {
        for (size_t i = N0 + 1; i-- > 0;) {
            data[i] *= other.data[0];
            for (size_t j = 0; j < i; ++j) {
                data[i] += T(binom2(i, i - j)) * data[j] * other.data[i - j];
            }
        }

        return *this;
    }

    dual &operator*=(const T &other) {
        for (size_t i = 0; i <= N0; ++i) {
            data[i] *= other;
        }

        return *this;
    }

    dual &operator/=(const dual &other) {
        for (size_t i = 0; i <= N0; ++i) {
            for (size_t j = 1; j <= i; ++j) {
                data[i] -= T(binom2(i, i - j)) * other.data[j] * data[i - j];
            }
            data[i] /= other.data[0];
        }

        return *this;
    }

    dual &operator/=(const T &other) {
        for (size_t i = 0; i <= N0; ++i) {
            data[i] /= other;
        }

        return *this;
    }

    value_type &value() { return dual_get(*this); }

    const value_type &value() const { return dual_get(*this); }

    element_type &operator[](size_t i) { return data[i]; }

    const element_type &operator[](size_t i) const { return data[i]; }

    static constexpr size_t max_order() { return std::max({N0, N...}); }

    static constexpr size_t max_order2 = std::max({N0, N...});
};

template <typename T>
void dual_assign_grad(const dual<T, 0> &x, std::tuple<T &> res) {
    auto &[res0] = res;

    res0 = x[0];
}

template <typename T>
void dual_assign_grad(const dual<T, 1> &x, std::tuple<T &, T &> res) {
    auto &[res0, res1] = res;

    res0 = x[0];
    res1 = x[1];
}

template <typename T>
void dual_assign_grad(const dual<T, 2> &x, std::tuple<T &, T &, T &> res) {
    auto &[res0, res1, res2] = res;

    res0 = x[0];
    res1 = x[1];
    res2 = x[2];
}

template <typename T>
void dual_assign_grad(const dual<T, 0, 0> &x, std::tuple<T &> res) {
    auto &[res0] = res;

    res0 = x[0][0];
}

template <typename T>
void dual_assign_grad(const dual<T, 1, 1> &z, std::tuple<T &, T (&)[2]> res) {
    auto &[res0, res1] = res;

    res0 = z[0][0];

    res1[0] = z[1][0];
    res1[1] = z[0][1];
}

template <typename T>
void dual_assign_grad(const dual<T, 2, 2> &z, std::tuple<T &, T (&)[2], T (&)[2][2]> res) {
    auto &[res0, res1, res2] = res;

    res0 = z[0][0];

    res1[0] = z[1][0];
    res1[1] = z[0][1];

    res2[0][0] = z[2][0];
    res2[0][1] = z[1][1];
    res2[1][0] = z[1][1];
    res2[1][1] = z[0][2];
}

template <typename T>
void dual_assign_grad(
    const dual<T, 1, 1> &z, std::tuple<T &, std::mdspan<T, std::dextents<ptrdiff_t, 1>, std::layout_stride>> res
) {
    auto &[res0, res1] = res;

    res0 = z[0][0];

    res1(0) = z[1][0];
    res1(1) = z[0][1];
}

template <typename T>
void dual_assign_grad(
    const dual<T, 2, 2> &z, std::tuple<
                                T &, std::mdspan<T, std::dextents<ptrdiff_t, 1>, std::layout_stride>,
                                std::mdspan<T, std::dextents<ptrdiff_t, 2>, std::layout_stride>>
                                res
) {
    auto &[res0, res1, res2] = res;

    res0 = z[0][0];

    res1(0) = z[1][0];
    res1(1) = z[0][1];

    res2(0, 0) = z[2][0];
    res2(0, 1) = z[1][1];
    res2(1, 0) = z[1][1];
    res2(1, 1) = z[0][2];
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> operator-(const dual<T, N0, N...> &rhs) {
    dual<T, N0, N...> out;
    for (size_t i = 0; i <= N0; ++i) {
        out[i] = -rhs[i];
    }

    return out;
}

template <typename T, size_t... N>
dual<T, N...> operator+(const dual<T, N...> &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator+(const dual<T, N...> &lhs, const T &rhs) {
    dual<T, N...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator+(const T &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator-(const dual<T, N...> &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res -= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator-(const dual<T, N...> &lhs, const T &rhs) {
    dual<T, N...> res = lhs;
    res -= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator*(const dual<T, N...> &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator*(const dual<T, N...> &lhs, const T &rhs) {
    dual<T, N...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator*(const T &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = rhs;
    res *= lhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator/(const dual<T, N...> &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator/(const dual<T, N...> &lhs, const T &rhs) {
    dual<T, N...> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator/(const T &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res /= rhs;

    return res;
}
template <typename T, size_t... N>
bool operator<(const dual<T, N...> &lhs, const dual<T, N...> &rhs) {
    return lhs.value() < rhs.value();
}

template <typename T, size_t... N, typename U>
bool operator<(const dual<T, N...> &lhs, const U &rhs) {
    return lhs.value() < rhs;
}

template <typename T, size_t... N, typename U>
bool operator==(const dual<T, N...> &lhs, const U &rhs) {
    return lhs.value() == rhs;
}

template <size_t Order, typename T>
dual<T, Order> make_dual(T value) {
    if (Order == 0) {
        return {value};
    }

    return {value, 1};
}

template <typename T, size_t N, size_t... Orders>
dual<T, Orders...> dual_taylor_series(const T (&coef)[N], const dual<T, Orders...> &x) {
    dual<T, Orders...> res = coef[0];
    if (N >= 2) {
        dual<T, Orders...> y = x - x.value();
        T denom = 1;

        res += y * coef[1];

        for (size_t i = 2; i < N; ++i) {
            y *= x - x.value();
            denom *= T(i);

            res += y * coef[i] / denom;
        }
    }

    return res;
}

template <typename T, size_t... N>
dual<T, 0, N...> sqrt(const dual<T, 0, N...> &z) {
    T coef[1] = {sqrt(z.value())};

    return dual_taylor_series(coef, z);
}

template <typename T, size_t... N>
dual<T, 1, N...> sqrt(const dual<T, 1, N...> &z) {
    T z0_sqrt = sqrt(z.value());
    T coef[2] = {z0_sqrt, T(1) / (T(2) * z0_sqrt)};

    return dual_taylor_series(coef, z);
}

template <typename T, size_t... N>
dual<T, 2, N...> sqrt(const dual<T, 2, N...> &z) {
    T z0_sqrt = sqrt(z.value());
    T coef[3] = {z0_sqrt, T(1) / (T(2) * z0_sqrt), -T(1) / (T(4) * z0_sqrt * z.value())};

    return dual_taylor_series(coef, z);
}

template <typename T, size_t... Orders>
dual<T, Orders...> sin(const dual<T, Orders...> &x) {
    static constexpr size_t MaxOrder = dual<T, Orders...>::max_order();

    T coef[MaxOrder + 1] = {sin(x.value())};
    if (MaxOrder >= 1) {
        coef[1] = cos(x.value());

        for (size_t i = 2; i <= MaxOrder; ++i) {
            coef[i] = -coef[i - 2];
        }
    }

    return dual_taylor_series(coef, x);
}

template <typename T, size_t... Orders>
dual<T, Orders...> cos(const dual<T, Orders...> &x) {
    static constexpr size_t MaxOrder = dual<T, Orders...>::max_order();

    T coef[MaxOrder + 1] = {cos(x.value())};
    if (MaxOrder >= 1) {
        coef[1] = -sin(x.value());

        for (size_t i = 2; i <= MaxOrder; ++i) {
            coef[i] = -coef[i - 2];
        }
    }

    return dual_taylor_series(coef, x);
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> real(dual<std::complex<T>, N0, N...> z) {
    dual<T, N0, N...> res;
    for (size_t i = 0; i <= N0; ++i) {
        res[i] = real(z[i]);
    }

    return res;
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> real(dual<T, N0, N...> z) {
    dual<T, N0, N...> res;
    for (size_t i = 0; i <= N0; ++i) {
        res[i] = real(z[i]);
    }

    return res;
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> imag(dual<std::complex<T>, N0, N...> z) {
    dual<T, N0, N...> res;
    for (size_t i = 0; i <= N0; ++i) {
        res[i] = imag(z[i]);
    }

    return res;
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> imag(dual<T, N0, N...> z) {
    dual<T, N0, N...> res;
    for (size_t i = 0; i <= N0; ++i) {
        res[i] = imag(z[i]);
    }

    return res;
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> exp(const dual<T, N0, N...> &z) {
    constexpr size_t M = N0 + 10;

    T coef[M] = {exp(z.value())};
    for (size_t i = 1; i < M; ++i) {
        coef[i] = coef[0];
    }

    return dual_taylor_series(coef, z);
}

template <typename T, size_t... Orders>
bool isfinite(const dual<T, Orders...> &x) {
    return isfinite(x.value());
}

template <typename T, size_t... Orders>
bool isinf(const dual<T, Orders...> &x) {
    return isinf(x.value());
}

using std::isnan;

template <typename T, size_t... Orders>
bool isnan(const dual<T, Orders...> &x) {
    return isnan(x.value());
}

using std::copysign;
using std::signbit;

template <typename T, size_t... Orders>
dual<T, Orders...> copysign(const dual<T, Orders...> &z, const dual<T, Orders...> &z2) {
    if (signbit(z.value()) == signbit(z2.value())) {
        return z;
    }

    return -z;
}

template <typename T, size_t... N>
struct complex_type<dual<T, N...>> {
    using type = dual<typename complex_type<T>::type, N...>;
};

} // namespace xsf

namespace std {

template <typename T, size_t... N>
xsf::dual<T, 0, N...> abs(xsf::dual<T, 0, N...> z) {
    return dual_taylor_series({std::abs(z.value())}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 0, N...> abs(xsf::dual<std::complex<T>, 0, N...> z) {
    return dual_taylor_series({std::abs(z.value())}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 1, N...> abs(xsf::dual<T, 1, N...> z) {
    if (z < 0) {
        return dual_taylor_series({std::abs(z.value()), T(-1)}, z);
    }

    return dual_taylor_series({std::abs(z.value()), T(1)}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 1, N...> abs(xsf::dual<std::complex<T>, 1, N...> z) {
    if (std::real(z) < 0) {
        return dual_taylor_series({std::abs(z.value()), std::real(z.value()) / std::abs(z.value())}, z);
    }

    return dual_taylor_series({std::abs(z.value()), std::real(z.value()) / std::abs(z.value())}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 2, N...> abs(xsf::dual<T, 2, N...> z) {
    if (z < 0) {
        return dual_taylor_series({std::abs(z.value()), T(-1), T(0)}, z);
    }

    return dual_taylor_series({std::abs(z.value()), T(1), T(0)}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 2, N...> abs(xsf::dual<std::complex<T>, 2, N...> z) {
    if (std::real(z.value()) < 0) {
        return dual_taylor_series({std::abs(z.value()), std::real(z.value()) / std::abs(z.value()), T(0)}, z);
    }

    return dual_taylor_series({std::abs(z.value()), std::real(z.value()) / std::abs(z.value()), T(0)}, z);
}

} // namespace std
