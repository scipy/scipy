#pragma once

#include <tuple>

#include "config.h"
#include "tuples.h"

namespace xsf {

template <typename T, size_t N, size_t... M>
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

template <typename T, size_t N, size_t... M>
struct dual {
  public:
    using value_type = T;
    using element_type = dual_element_t<T, N, M...>;
    using reference = value_type &;
    using const_reference = const value_type &;
    using dual_element = element_type;

  public:
    element_type data[N + 1];

  public:
    dual() = default;

    dual(T value) {
        data[0] = value;
        for (size_t i = 1; i <= N; ++i) {
            data[i] = 0;
        }
    }

    dual(std::initializer_list<element_type> data) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            this->data[it - data.begin()] = *it;
        }
        for (size_t i = data.size(); i <= N; ++i) {
            this->data[i] = 0;
        }
    }

    template <typename U>
    dual(const dual<U, N, M...> &other) {
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

    dual &operator+=(const dual_element &other) {
        data[0] += other;

        return *this;
    }

    dual &operator-=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] -= other.data[i];
        }

        return *this;
    }

    dual &operator-=(const dual_element &other) {
        data[0] -= other;

        return *this;
    }

    dual &operator*=(const dual &other) {
        for (size_t i = N + 1; i-- > 0;) {
            data[i] *= other.data[0];
            for (size_t j = 0; j < i; ++j) {
                data[i] += data[j] * other.data[i - j];
            }
        }

        return *this;
    }

    dual &operator*=(const dual_element &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] *= other;
        }

        return *this;
    }

    dual &operator/=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            for (size_t j = 1; j <= i; ++j) {
                data[i] -= other.data[j] * data[i - j];
            }
            data[i] /= other.data[0];
        }

        return *this;
    }

    dual &operator/=(const dual_element &other) {
        for (size_t i = 0; i <= N; ++i) {
            data[i] /= other;
        }

        return *this;
    }

    value_type &value() { return dual_get(*this); }

    const value_type &value() const { return dual_get(*this); }

    element_type &operator[](size_t i) { return data[i]; }

    const element_type &operator[](size_t i) const { return data[i]; }
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
    res2 = T(2) * x[2];
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

    res2[0][0] = T(2) * z[2][0];
    res2[0][1] = z[1][1];
    res2[1][0] = z[1][1];
    res2[1][1] = T(2) * z[0][2];
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

    res2(0, 0) = T(2) * z[2][0];
    res2(0, 1) = z[1][1];
    res2(1, 0) = z[1][1];
    res2(1, 1) = T(2) * z[0][2];
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
dual<T, N...> operator+(const dual<T, N...> &lhs, const dual_element_t<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator+(const dual_element_t<T, N...> &lhs, const dual<T, N...> &rhs) {
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
dual<T, N...> operator*(const dual<T, N...> &lhs, const dual<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator*(const dual<T, N...> &lhs, const dual_element_t<T, N...> &rhs) {
    dual<T, N...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... N>
dual<T, N...> operator*(const dual_element_t<T, N...> &lhs, const dual<T, N...> &rhs) {
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
dual<T, N...> operator/(const dual<T, N...> &lhs, const dual_element_t<T, N...> &rhs) {
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
dual<T, N...> operator/(const dual_element<T, N...> &lhs, const dual<T, N...> &rhs) {
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

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> dual_apply2(const dual_element_t<T, N0, N...> (&coef)[N0 + 1], const dual<T, N0, N...> &z) {
    dual<T, N0, N...> res = T(0);

    dual<T, N0, N...> x = z;
    x[0] = 0;

    T fac = T(1);

    res[0] = coef[0];
    for (size_t i = 1; i <= N0; ++i) {
        fac *= T(i);
        res += x * (coef[i] / fac);
        x = x * x;
    }

    return res;
}

template <typename T, size_t... N>
dual<T, 0, N...> sqrt(const dual<T, 0, N...> &z) {
    dual_element_t<T, 0, N...> coef[1] = {sqrt(z[0])};

    return dual_apply2(coef, z);
}

template <typename T, size_t... N>
dual<T, 1, N...> sqrt(const dual<T, 1, N...> &z) {
    dual_element_t<T, 1, N...> z0_sqrt = sqrt(z[0]);
    dual_element_t<T, 1, N...> coef[2] = {z0_sqrt, T(1) / (T(2) * z0_sqrt)};

    return dual_apply2(coef, z);
}

template <typename T, size_t... N>
dual<T, 2, N...> sqrt(const dual<T, 2, N...> &z) {
    dual_element_t<T, 2, N...> z0_sqrt = sqrt(z[0]);
    dual_element_t<T, 2, N...> coef[3] = {z0_sqrt, T(1) / (T(2) * z0_sqrt), -T(1) / (T(4) * z0_sqrt * z[0])};

    return dual_apply2(coef, z);
}

template <typename T, size_t... N>
dual<T, 0, N...> sin(const dual<T, 0, N...> &z) {
    return dual_apply2({sin(z[0])}, z);
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> sin(const dual<T, N0, N...> &z) {
    dual_element_t<T, N0, N...> coef[N0 + 1] = {sin(z[0]), cos(z[0])};
    for (size_t i = 2; i <= N0; ++i) {
        coef[i] = -coef[i - 2];
    }

    return dual_apply2(coef, z);
}

template <typename T, size_t... N>
dual<T, 0, N...> cos(const dual<T, 0, N...> &z) {
    return dual_apply2({cos(z[0])}, z);
}

template <typename T, size_t N0, size_t... N>
dual<T, N0, N...> cos(const dual<T, N0, N...> &z) {
    dual_element_t<T, N0, N...> coef[N0 + 1] = {cos(z[0]), -sin(z[0])};
    for (size_t i = 2; i <= N0; ++i) {
        coef[i] = -coef[i - 2];
    }

    return dual_apply2(coef, z);
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
    dual_element_t<T, N0, N...> coef[N0 + 1] = {exp(z[0])};
    for (size_t i = 1; i <= N0; ++i) {
        coef[i] = coef[0];
    }

    return dual_apply2(coef, z);
}

template <typename T, size_t... N>
bool isfinite(const dual<T, N...> &z) {
    return isfinite(z.value());
}

template <typename T, size_t... N>
bool isinf(const dual<T, N...> &z) {
    return isinf(z.value());
}

using std::isnan;

template <typename T, size_t... N>
bool isnan(const dual<T, N...> &z) {
    return isnan(z.value());
}

using std::copysign;

template <typename T, size_t... N>
dual<T, N...> copysign(const dual<T, N...> &z, const dual<T, N...> &z2) {
    if (z2 < 0) {
        return -z;
    }

    return z;
}

template <typename T, size_t... N>
struct complex_type<dual<T, N...>> {
    using type = dual<typename complex_type<T>::type, N...>;
};

} // namespace xsf

namespace std {

template <typename T, size_t... N>
xsf::dual<T, 0, N...> abs(xsf::dual<T, 0, N...> z) {
    return dual_apply2({std::abs(z[0])}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 0, N...> abs(xsf::dual<std::complex<T>, 0, N...> z) {
    return dual_apply2({std::abs(z[0])}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 1, N...> abs(xsf::dual<T, 1, N...> z) {
    if (z < 0) {
        return dual_apply2({std::abs(z[0]), T(-1)}, z);
    }

    return dual_apply2({std::abs(z[0]), T(1)}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 1, N...> abs(xsf::dual<std::complex<T>, 1, N...> z) {
    if (std::real(z) < 0) {
        return dual_apply2({std::abs(z[0]), std::real(z[0]) / std::abs(z[0])}, z);
    }

    return dual_apply2({std::abs(z[0]), std::real(z[0]) / std::abs(z[0])}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 2, N...> abs(xsf::dual<T, 2, N...> z) {
    if (z < 0) {
        return dual_apply2({std::abs(z[0]), T(-1), T(0)}, z);
    }

    return dual_apply2({std::abs(z[0]), T(1), T(0)}, z);
}

template <typename T, size_t... N>
xsf::dual<T, 2, N...> abs(xsf::dual<std::complex<T>, 2, N...> z) {
    if (std::real(z.value()) < 0) {
        return dual_apply2({std::abs(z[0]), std::real(z[0]) / std::abs(z[0]), T(0)}, z);
    }

    return dual_apply2({std::abs(z[0]), std::real(z[0]) / std::abs(z[0]), T(0)}, z);
}

} // namespace std
