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
        for (int i = values.size(); i <= N; ++i) {
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

    dual &operator+=(const dual &rhs) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] += rhs.values[i];
        }

        return *this;
    }

    dual &operator+=(const T &rhs) {
        values[0] += rhs;

        return *this;
    }

    dual &operator-=(const dual &rhs) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] -= rhs.values[i];
        }

        return *this;
    }

    dual &operator/=(const dual &rhs) {
        *this = (*this / rhs);

        return *this;
    }

    dual &operator/=(const T &rhs) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] /= rhs;
        }

        return *this;
    }

    dual from_coef(const T (&coef)[N + 1]) {
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

    template <typename Func>
    dual<std::result_of_t<Func(T)>, N> apply(Func f) {
        using U = std::result_of_t<Func(T)>;

        dual<U, N> res;
        for (size_t i = 0; i <= N; ++i) {
            res.values[i] = f(values[i]);
        }

        return res;
    }
};

template <typename T, size_t N>
dual<T, N> operator-(const dual<T, N> &rhs) {
    dual<T, N> out;
    for (int i = 0; i <= N; ++i) {
        out.values[i] = -rhs.values[i];
    }

    return out;
}

template <typename T, size_t N>
dual<T, N> operator+(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> out;
    for (int i = 0; i <= N; ++i) {
        out.values[i] = lhs.values[i] + rhs.values[i];
    }

    return out;
}

template <typename T, size_t N>
dual<T, N> operator-(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> out;
    for (int i = 0; i <= N; ++i) {
        out.values[i] = lhs.values[i] - rhs.values[i];
    }

    return out;
}

template <typename T, size_t N>
dual<T, N> operator*(const T &lhs, const dual<T, N> &rhs) {
    dual<T, N> res;
    for (int i = 0; i <= N; ++i) {
        res.values[i] = lhs * rhs[i];
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> operator*(const dual<T, N> &lhs, const T &rhs) {
    dual<T, N> res;
    for (int i = 0; i <= N; ++i) {
        res.values[i] = lhs.values[i] * rhs;
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> operator*(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> res;
    for (int i = 0; i <= N; ++i) {
        res.values[i] = 0;
        for (int j = 0; j <= i; ++j) {
            res.values[i] += lhs.values[j] * rhs.values[i - j];
        }
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> operator/(const dual<T, N> &lhs, const dual<T, N> &rhs) {
    dual<T, N> res;
    for (int i = 0; i <= N; ++i) {
        res.values[i] = lhs.values[i];
        for (int j = 1; j <= i; ++j) {
            res.values[i] -= rhs.values[j] * res.values[i - j];
        }
        res.values[i] /= rhs.values[0];
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> operator/(const T &lhs, const dual<T, N> &rhs) {
    dual<T, N> res;
    for (int i = 0; i <= N; ++i) {
        res[i] = lhs / rhs[i];
    }

    return res;
}

template <typename T, size_t N>
dual<T, N> operator/(const dual<T, N> &lhs, const T &rhs) {
    dual<T, N> res;
    for (int i = 0; i <= N; ++i) {
        res[i] = lhs[i] / rhs;
    }

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

} // namespace xsf

namespace std {

template <typename T>
xsf::dual<T, 0> abs(xsf::dual<T, 0> z) {
    return z.from_coef({std::abs(z[0])});
}

template <typename T>
xsf::dual<T, 0> abs(xsf::dual<std::complex<T>, 0> z) {
    return z.from_coef({std::abs(z[0])});
}

template <typename T>
xsf::dual<T, 1> abs(xsf::dual<T, 1> z) {
    if (z < 0) {
        return z.from_coef({std::abs(z.value()), T(-1)});
    }

    return z.from_coef({std::abs(z.value()), T(1)});
}

template <typename T>
xsf::dual<T, 1> abs(xsf::dual<std::complex<T>, 1> z) {
    if (std::real(z) < 0) {
        return z.from_coef({std::abs(z.value()), std::real(z[0]) / std::abs(z[0])});
    }

    return z.from_coef({std::abs(z.value()), std::real(z[0]) / std::abs(z[0])});
}

template <typename T>
xsf::dual<T, 2> abs(xsf::dual<T, 2> z) {
    if (z < 0) {
        return z.from_coef({std::abs(z.value()), T(-1), T(0)});
    }

    return z.from_coef({std::abs(z.value()), T(1), T(0)});
}

template <typename T>
xsf::dual<T, 2> abs(xsf::dual<std::complex<T>, 2> z) {
    if (std::real(z) < 0) {
        return z.from_coef({std::abs(z.value()), std::real(z[0]) / std::abs(z[0]), T(0)});
    }

    return z.from_coef({std::abs(z.value()), std::real(z[0]) / std::abs(z[0]), T(0)});
}

template <typename T>
xsf::dual<T, 0> sqrt(xsf::dual<T, 0> z) {
    T coef[1] = {std::sqrt(z.value())};

    return z.from_coef(coef);
}

template <typename T>
xsf::dual<T, 1> sqrt(xsf::dual<T, 1> z) {
    T z_sqrt = std::sqrt(z.value());

    return z.from_coef({z_sqrt, T(1) / (T(2) * z_sqrt)});
}

template <typename T>
xsf::dual<T, 2> sqrt(xsf::dual<T, 2> z) {
    T z_sqrt = std::sqrt(z.value());

    return z.from_coef({z_sqrt, T(1) / (T(2) * z_sqrt), -T(1) / (T(4) * z.value() * z_sqrt)});
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
