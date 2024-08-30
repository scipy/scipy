#pragma once

namespace xsf {

template <typename T>
std::tuple<T> derivatives_impl(const T (&values)[1]) {
    return {values[0]};
}

template <typename T>
std::tuple<T, T> derivatives_impl(const T (&values)[2]) {
    return {values[0], values[1]};
}

template <typename T>
std::tuple<T, T, T> derivatives_impl(const T (&values)[3]) {
    return {values[0], values[1], T(2) * values[2]};
}

template <typename T, size_t N>
struct dual {
    T values[N + 1];

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

    dual &operator=(const T &rhs) {
        values[0] = rhs;
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
        dual res(0);

        dual x = *this;
        x.values[0] = 0;

        T fac = 1;

        res.values[0] = coef[0];
        for (size_t i = 1; i <= N; ++i) {
            fac *= T(i);
            res += x * (coef[i] / fac);
            x = x * x;
        }

        return res;
    }

    T value() const { return values[0]; }

    T jac() const { return values[1]; }

    T hess() const { return T(2) * values[2]; }

    decltype(auto) derivatives() const { return derivatives_impl(values); }

    T &operator[](size_t i) { return values[i]; }

    const T &operator[](size_t i) const { return values[i]; }
};

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
dual<T, N> operator<(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    return lhs[0] < rhs[0];
}

template <typename T, typename U, size_t N>
bool operator<(const dual<T, N> &lhs, const U &rhs) {
    return lhs.values[0] < rhs;
}

template <typename T, typename U, size_t N>
bool operator==(const dual<T, N> &lhs, const U &rhs) {
    return lhs.values[0] == rhs;
}

template <size_t N, typename T>
dual<T, N> make_dual(T value) {
    if (N == 0) {
        return {value};
    }

    return {value, 1};
}

using std::sqrt;

template <typename T>
dual<T, 0> sqrt(const dual<T, 0> &z) {
    T z0_sqrt = std::sqrt(z[0]);

    return z.apply({z0_sqrt});
}

template <typename T>
dual<T, 1> sqrt(const dual<T, 1> &z) {
    T z0_sqrt = std::sqrt(z[0]);

    return z.apply({z0_sqrt, T(1) / (T(2) * z0_sqrt)});
}

template <typename T>
dual<T, 2> sqrt(const dual<T, 2> &z) {
    T z0_sqrt = std::sqrt(z[0]);

    return z.apply({z0_sqrt, T(1) / (T(2) * z0_sqrt), -T(1) / (T(4) * z0_sqrt * z[0])});
}

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
    if (std::real(z) < 0) {
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
    if (std::real(z) < 0) {
        return z.apply({std::abs(z.value()), std::real(z[0]) / std::abs(z[0]), T(0)});
    }

    return z.apply({std::abs(z.value()), std::real(z[0]) / std::abs(z[0]), T(0)});
}

template <typename T, size_t N>
xsf::dual<T, N> real(xsf::dual<complex<T>, N> z) {
    xsf::dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = std::real(z[i]);
    }

    return res;
}

template <typename T, size_t N>
xsf::dual<T, N> real(xsf::dual<T, N> z) {
    xsf::dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = std::real(z[i]);
    }

    return res;
}

template <typename T, size_t N>
xsf::dual<T, N> imag(xsf::dual<complex<T>, N> z) {
    xsf::dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = std::imag(z[i]);
    }

    return res;
}

template <typename T, size_t N>
xsf::dual<T, N> imag(xsf::dual<T, N> z) {
    xsf::dual<T, N> res;
    for (size_t i = 0; i <= N; ++i) {
        res[i] = std::imag(z[i]);
    }

    return res;
}

} // namespace std
