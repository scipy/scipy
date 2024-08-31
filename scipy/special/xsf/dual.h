#pragma once

#include <tuple>

#include "config.h"
#include "tuples.h"

namespace xsf {

template <typename T, size_t... N>
struct dual;

template <typename T, size_t N>
struct dual<T, N> {
  public:
    using value_type = T;
    using reference = value_type &;
    using const_reference = const value_type &;

  public:
    T values[N + 1];

  public:
    dual() = default;

    dual(T value) {
        values[0] = value;
        for (size_t i = 1; i <= N; ++i) {
            values[i] = 0;
        }
    }

    dual(std::initializer_list<T> values) {
        for (auto it = values.begin(); it != values.end(); ++it) {
            this->values[it - values.begin()] = *it;
        }
        for (size_t i = values.size(); i <= N; ++i) {
            this->values[i] = 0;
        }
    }

    template <typename U>
    dual(const dual<U, N> &other) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] = other.values[i];
        }
    }

    dual &operator=(const T &other) {
        values[0] = other;
        for (size_t i = 1; i <= N; ++i) {
            values[i] = 0;
        }

        return *this;
    }

    dual &operator+=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] += other[i];
        }

        return *this;
    }

    dual &operator+=(const T &other) {
        values[0] += other;

        return *this;
    }

    dual &operator-=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] -= other[i];
        }

        return *this;
    }

    dual &operator-=(const T &other) {
        values[0] -= other;

        return *this;
    }

    dual &operator*=(const dual &other) {
        T tmp[N + 1];
        for (size_t i = 0; i <= N; ++i) {
            tmp[i] = values[i];
        }

        for (size_t i = 0; i <= N; ++i) {
            values[i] = 0;
            for (size_t j = 0; j <= i; ++j) {
                values[i] += tmp[j] * other[i - j];
            }
        }

        return *this;
    }

    dual &operator*=(const T &other) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] *= other;
        }

        return *this;
    }

    dual &operator/=(const dual &other) {
        for (size_t i = 0; i <= N; ++i) {
            for (size_t j = 1; j <= i; ++j) {
                values[i] -= other.values[j] * values[i - j];
            }
            values[i] /= other.values[0];
        }

        return *this;
    }

    dual &operator/=(const T &other) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] /= other;
        }

        return *this;
    }

    dual apply(const T (&coef)[N + 1]) const {
        dual res = T(0);

        dual x = *this;
        x.values[0] = 0;

        T fac = T(1);

        res.values[0] = coef[0];
        for (size_t i = 1; i <= N; ++i) {
            fac *= T(i);
            res += x * (coef[i] / fac);
            x = x * x;
        }

        return res;
    }

    reference front() { return values[0]; }

    const_reference front() const { return values[0]; }

    T value() const { return values[0]; }

    T &operator[](size_t i) { return values[i]; }

    const T &operator[](size_t i) const { return values[i]; }
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
void dual_assign_grad(const dual<dual<T, 0>, 0> &x, std::tuple<T &> res) {
    auto &[res0] = res;

    res0 = x[0][0];
}

template <typename T>
void dual_assign_grad(const dual<dual<T, 1>, 1> &z, std::tuple<T &, T (&)[2]> res) {
    auto &[res0, res1] = res;

    res0 = z[0][0];

    res1[0] = z[1][0];
    res1[1] = z[0][1];
}

template <typename T>
void dual_assign_grad(const dual<dual<T, 2>, 2> &z, std::tuple<T &, T (&)[2], T (&)[2][2]> res) {
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
void dual_assign_grad(const std::complex<dual<dual<T, 0>, 0>> &z, std::tuple<std::complex<T> &> res) {
    dual<dual<T, 0>, 0> z_real = real(z);
    dual<dual<T, 0>, 0> z_imag = imag(z);

    auto &[res0] = res;

    res0 = std::complex(z_real[0][0], z_imag[0][0]);
}

template <typename T>
void dual_assign_grad(
    const std::complex<dual<dual<T, 1>, 1>> &z, std::tuple<std::complex<T> &, std::complex<T> (&)[2]> res
) {
    dual<dual<T, 1>, 1> z_real = real(z);
    dual<dual<T, 1>, 1> z_imag = imag(z);

    auto &[res0, res1] = res;

    res0 = std::complex(z_real[0][0], z_imag[0][0]);

    res1[0] = std::complex(z_real[1][0], z_imag[1][0]);
    res1[1] = std::complex(z_real[0][1], z_imag[0][1]);
}

template <typename T>
void dual_assign_grad(
    const std::complex<dual<dual<T, 1>, 1>> &z,
    std::tuple<std::complex<T> &, std::mdspan<std::complex<T>, std::dextents<ptrdiff_t, 1>, std::layout_stride>> res
) {
    dual<dual<T, 1>, 1> z_real = real(z);
    dual<dual<T, 1>, 1> z_imag = imag(z);

    auto &[res0, res1] = res;

    res0 = std::complex(z_real[0][0], z_imag[0][0]);

    res1(0) = std::complex(z_real[1][0], z_imag[1][0]);
    res1(1) = std::complex(z_real[0][1], z_imag[0][1]);
}

template <typename T>
void dual_assign_grad(
    const std::complex<dual<dual<T, 2>, 2>> &z,
    std::tuple<std::complex<T> &, std::complex<T> (&)[2], std::complex<T> (&)[2][2]> res
) {
    dual<dual<T, 2>, 2> z_real = real(z);
    dual<dual<T, 2>, 2> z_imag = imag(z);

    auto &[res0, res1, res2] = res;

    res0 = std::complex(z_real[0][0], z_imag[0][0]);

    res1[0] = std::complex(z_real[1][0], z_imag[1][0]);
    res1[1] = std::complex(z_real[0][1], z_imag[0][1]);

    res2[0][0] = T(2) * std::complex(z_real[2][0], z_imag[2][0]);
    res2[0][1] = std::complex(z_real[1][1], z_imag[1][1]);
    res2[1][0] = std::complex(z_real[1][1], z_imag[1][1]);
    res2[1][1] = T(2) * std::complex(z_real[0][2], z_imag[0][2]);
}

template <typename T>
void dual_assign_grad(const dual<dual<T, 1>, 1> &z, std::tuple<T &> res) {
    auto &[res0] = res;

    res0 = z[0][0];
}

template <typename T>
void dual_assign_grad(
    const dual<dual<T, 1>, 1> &z, std::tuple<T &, std::mdspan<T, std::dextents<ptrdiff_t, 1>, std::layout_stride>> res
) {
    auto &[res0, res1] = res;

    res0 = z[0][0];

    res1(0) = z[1][0];
    res1(1) = z[0][1];
}

template <typename T>
void dual_assign_grad(
    const dual<dual<T, 2>, 2> &z, std::tuple<
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

template <typename T>
void dual_assign_grad(
    const std::complex<dual<dual<T, 2>, 2>> &z,
    std::tuple<
        std::complex<T> &, std::mdspan<std::complex<T>, std::dextents<ptrdiff_t, 1>, std::layout_stride>,
        std::mdspan<std::complex<T>, std::dextents<ptrdiff_t, 2>, std::layout_stride>>
        res
) {
    dual<dual<T, 2>, 2> z_real = real(z);
    dual<dual<T, 2>, 2> z_imag = imag(z);

    auto &[res0, res1, res2] = res;

    res0 = std::complex(z_real[0][0], z_imag[0][0]);

    res1(0) = std::complex(z_real[1][0], z_imag[1][0]);
    res1(1) = std::complex(z_real[0][1], z_imag[0][1]);

    res2(0, 0) = T(2) * std::complex(z_real[2][0], z_imag[2][0]);
    res2(0, 1) = std::complex(z_real[1][1], z_imag[1][1]);
    res2(1, 0) = std::complex(z_real[1][1], z_imag[1][1]);
    res2(1, 1) = T(2) * std::complex(z_real[0][2], z_imag[0][2]);
}

template <typename T, size_t N>
dual<T, N> operator-(const dual<T, N> &rhs) {
    dual<T, N> out;
    for (size_t i = 0; i <= N; ++i) {
        out.values[i] = -rhs.values[i];
    }

    return out;
}

template <typename T, size_t N>
dual<T, N> operator+(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator+(const dual<T, N> &lhs, const T &rhs) {
    dual<T, N> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator+(const T &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator-(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = lhs;
    res -= rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator*(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator*(const dual<T, N> &lhs, const T &rhs) {
    dual<T, N> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator*(const T &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = rhs;
    res *= lhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator/(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator/(const dual<T, N> &lhs, const T &rhs) {
    dual<T, N> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t N>
dual<T, N> operator/(const T &lhs, const dual<T, N> &rhs) {
    dual<T, N> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t N>
bool operator<(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    return lhs.front() < rhs.front();
}

template <typename T, size_t N, typename U>
bool operator<(const dual<T, N> &lhs, const U &rhs) {
    return lhs.front() < rhs;
}

template <typename T, size_t N, typename U>
bool operator==(const dual<T, N> &lhs, const U &rhs) {
    return lhs.front() == rhs;
}

template <size_t N, typename T>
dual<T, N> make_dual(T value) {
    if (N == 0) {
        return {value};
    }

    return {value, 1};
}

template <typename T>
dual<T, 0> sqrt(const dual<T, 0> &z) {
    T z0_sqrt = sqrt(z[0]);

    return z.apply({z0_sqrt});
}

template <typename T>
dual<T, 1> sqrt(const dual<T, 1> &z) {
    T z0_sqrt = sqrt(z[0]);

    return z.apply({z0_sqrt, T(1) / (T(2) * z0_sqrt)});
}

template <typename T>
dual<T, 2> sqrt(const dual<T, 2> &z) {
    T z0_sqrt = sqrt(z[0]);

    return z.apply({z0_sqrt, T(1) / (T(2) * z0_sqrt), -T(1) / (T(4) * z0_sqrt * z[0])});
}

template <typename T>
dual<T, 0> sin(const dual<T, 0> &z) {
    return z.apply({sin(z.front())});
}

template <typename T, size_t N>
dual<T, N> sin(const dual<T, N> &z) {
    T coef[N + 1] = {sin(z.front()), cos(z.front())};
    for (size_t i = 2; i <= N; ++i) {
        coef[i] = -coef[i - 2];
    }

    return z.apply(coef);
}

template <typename T>
dual<T, 0> cos(const dual<T, 0> &z) {
    return z.apply({cos(z.front())});
}

template <typename T, size_t N>
dual<T, N> cos(const dual<T, N> &z) {
    T coef[N + 1] = {cos(z.front()), -sin(z.front())};
    for (size_t i = 2; i <= N; ++i) {
        coef[i] = -coef[i - 2];
    }

    return z.apply(coef);
}

template <typename T, size_t N>
dual<T, N> real(dual<std::complex<T>, N> z) {
    dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = real(z[i]);
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> real(dual<T, N> z) {
    dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = real(z[i]);
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> imag(dual<std::complex<T>, N> z) {
    dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = imag(z[i]);
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> imag(dual<T, N> z) {
    dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = imag(z[i]);
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> exp(const dual<T, N> &z) {
    T coef[N + 1] = {exp(z.front())};
    for (size_t i = 1; i <= N; ++i) {
        coef[i] = coef[0];
    }

    return z.apply(coef);
}

template <typename T, size_t N>
bool isfinite(const dual<T, N> &z) {
    return isfinite(z.front());
}

template <typename T, size_t N>
bool isinf(const dual<T, N> &z) {
    return isinf(z.front());
}

using std::copysign;

template <typename T, size_t N>
dual<T, N> copysign(const dual<T, N> &z, const dual<T, N> &z2) {
    return z;
}

template <typename T, size_t N>
struct complex_type<dual<T, N>> {
    using type = dual<typename complex_type<T>::type, N>;
};

} // namespace xsf

namespace std {

template <typename T>
xsf::dual<T, 0> abs(xsf::dual<T, 0> z) {
    return z.apply({std::abs(z[0])});
}

template <typename T>
xsf::dual<T, 0> abs(xsf::dual<std::complex<T>, 0> z) {
    return z.apply({std::abs(z[0])});
}

template <typename T>
xsf::dual<T, 1> abs(xsf::dual<T, 1> z) {
    if (z < 0) {
        return z.apply({std::abs(z.value()), T(-1)});
    }

    return z.apply({std::abs(z.value()), T(1)});
}

template <typename T>
xsf::dual<T, 1> abs(xsf::dual<std::complex<T>, 1> z) {
    if (std::real(z[0]) < 0) {
        return z.apply({std::abs(z.value()), std::real(z[0]) / std::abs(z[0])});
    }

    return z.apply({std::abs(z.value()), std::real(z[0]) / std::abs(z[0])});
}

template <typename T>
xsf::dual<T, 2> abs(xsf::dual<T, 2> z) {
    if (z < 0) {
        return z.apply({std::abs(z.value()), T(-1), T(0)});
    }

    return z.apply({std::abs(z.value()), T(1), T(0)});
}

template <typename T>
xsf::dual<T, 2> abs(xsf::dual<std::complex<T>, 2> z) {
    if (std::real(z[0]) < 0) {
        return z.apply({std::abs(z.value()), std::real(z[0]) / std::abs(z[0]), T(0)});
    }

    return z.apply({std::abs(z.value()), std::real(z[0]) / std::abs(z[0]), T(0)});
}

} // namespace std
