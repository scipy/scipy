#pragma once

#include "binom.h"
#include "config.h"
#include "numbers.h"

namespace xsf {

namespace detail {
    template <typename T>
    constexpr T small_binom_coefs[3][3] = {
        {T(1.0), T(0.0), T(0.0)}, {T(1.0), T(1.0), T(0.0)}, {T(1.0), T(2.0), T(1.0)}
    };

    /* Since we only compute derivatives up to order 2, we only need
     * Binomial coefficients with n <= 2 for use in the General
     * Leibniz rule. Get these from a lookup table. */
    template <typename T>
    T fast_binom(size_t n, size_t k) {
        if ((n <= 2) && (k <= 2)) {
            return small_binom_coefs<T>[n][k];
        }
        return T(xsf::binom(static_cast<double>(n), static_cast<double>(k)));
    }
} // namespace detail

template <typename T, size_t Order0, size_t... Orders>
class dual;

template <typename T, size_t Order>
class dual<T, Order> {
  public:
    using value_type = T;

  private:
    value_type data[Order + 1];

  public:
    dual() = default;

    dual(value_type value) {
        data[0] = value;
        for (size_t i = 1; i <= Order; ++i) {
            data[i] = 0;
        }
    }

    dual(std::initializer_list<value_type> data) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            this->data[it - data.begin()] = *it;
        }
        for (size_t i = data.size(); i <= Order; ++i) {
            this->data[i] = 0;
        }
    }

    template <typename U>
    dual(const dual<U, Order> &other) {
        for (size_t i = 0; i <= Order; ++i) {
            data[i] = other[i];
        }
    }

    dual &operator=(const value_type &other) {
        data[0] = other;
        for (size_t i = 1; i <= Order; ++i) {
            data[i] = 0;
        }

        return *this;
    }

    dual &operator+=(const dual &other) {
        for (size_t i = 0; i <= Order; ++i) {
            data[i] += other.data[i];
        }

        return *this;
    }

    dual &operator+=(const value_type &other) {
        data[0] += other;

        return *this;
    }

    dual &operator-=(const dual &other) {
        for (size_t i = 0; i <= Order; ++i) {
            data[i] -= other.data[i];
        }

        return *this;
    }

    dual &operator-=(const value_type &other) {
        data[0] -= other;

        return *this;
    }

    dual &operator*=(const dual &other) {
        for (size_t i = Order + 1; i-- > 0;) {
            data[i] *= other.data[0];
            // General Leibniz Rule
            for (size_t j = 0; j < i; ++j) {
                data[i] += detail::fast_binom<T>(i, j) * data[j] * other.data[i - j];
            }
        }

        return *this;
    }

    dual &operator*=(const value_type &other) {
        for (size_t i = 0; i <= Order; ++i) {
            data[i] *= other;
        }

        return *this;
    }

    dual &operator/=(const dual &other) {
        for (size_t i = 0; i <= Order; ++i) {
            for (size_t j = 1; j <= i; ++j) {
                data[i] -= detail::fast_binom<T>(i - 1, j) * other.data[j] * data[i - j];
            }

            data[i] /= other.data[0];
        }

        return *this;
    }

    dual &operator/=(const value_type &other) {
        for (size_t i = 0; i <= Order; ++i) {
            data[i] /= other;
        }

        return *this;
    }

    value_type &value() { return data[0]; }

    const value_type &value() const { return data[0]; }

    value_type &operator[](size_t i) { return data[i]; }

    const value_type &operator[](size_t i) const { return data[i]; }

    static constexpr size_t max_order() { return Order; }
};

template <typename T, size_t Order0, size_t Order1, size_t... Orders>
class dual<T, Order0, Order1, Orders...> {
  public:
    using value_type = T;

  private:
    dual<T, Order1, Orders...> data[Order0 + 1];

  public:
    dual() = default;

    dual(value_type value) {
        data[0] = value;
        for (size_t i = 1; i <= Order0; ++i) {
            data[i] = 0;
        }
    }

    dual(dual<value_type, Order1, Orders...> value) {
        data[0] = value;
        for (size_t i = 1; i <= Order0; ++i) {
            data[i] = 0;
        }
    }

    dual(std::initializer_list<dual<value_type, Order1, Orders...>> data) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            this->data[it - data.begin()] = *it;
        }
        for (size_t i = data.size(); i <= Order0; ++i) {
            this->data[i] = 0;
        }
    }

    template <typename U>
    dual(const dual<U, Order0, Order1, Orders...> &other) {
        for (size_t i = 0; i <= Order0; ++i) {
            data[i] = other[i];
        }
    }

    dual &operator=(const value_type &other) {
        data[0] = other;
        for (size_t i = 1; i <= Order0; ++i) {
            data[i] = 0;
        }

        return *this;
    }

    dual &operator+=(const dual &other) {
        for (size_t i = 0; i <= Order0; ++i) {
            data[i] += other.data[i];
        }

        return *this;
    }

    dual &operator+=(const value_type &other) {
        data[0] += other;

        return *this;
    }

    dual &operator-=(const dual &other) {
        for (size_t i = 0; i <= Order0; ++i) {
            data[i] -= other.data[i];
        }

        return *this;
    }

    dual &operator-=(const value_type &other) {
        data[0] -= other;

        return *this;
    }

    dual &operator*=(const dual &other) {
        for (size_t i = Order0 + 1; i-- > 0;) {
            data[i] *= other.data[0];
            // General Leibniz Rule
            for (size_t j = 0; j < i; ++j) {
                data[i] += detail::fast_binom<T>(i, j) * data[j] * other.data[i - j];
            }
        }

        return *this;
    }

    dual &operator*=(const value_type &other) {
        for (size_t i = 0; i <= Order0; ++i) {
            data[i] *= other;
        }

        return *this;
    }

    dual &operator/=(const dual &other) {
        for (size_t i = 0; i <= Order0; ++i) {
            for (size_t j = 1; j <= i; ++j) {
                data[i] -= detail::fast_binom<T>(i - 1, j) * other.data[j] * data[i - j];
            }

            data[i] /= other.data[0];
        }

        return *this;
    }

    dual &operator/=(const value_type &other) {
        for (size_t i = 0; i <= Order0; ++i) {
            data[i] /= other;
        }

        return *this;
    }

    value_type &value() {
        dual<value_type, Order1, Orders...> &data0 = data[0];

        return data0.value();
    }

    const value_type &value() const {
        const dual<value_type, Order1, Orders...> &data0 = data[0];

        return data0.value();
    }

    dual<value_type, Order1, Orders...> &operator[](size_t i) { return data[i]; }

    const dual<value_type, Order1, Orders...> &operator[](size_t i) const { return data[i]; }

    static constexpr size_t max_order() { return std::max({Order0, Order1, Orders...}); }
};

template <typename T, size_t Order0, size_t... Orders>
dual<T, Order0, Orders...> operator+(const dual<T, Order0, Orders...> &x) {
    dual<T, Order0, Orders...> res;
    for (size_t i = 0; i <= Order0; ++i) {
        res[i] = +x[i];
    }

    return res;
}

template <typename T, size_t Order0, size_t... Orders>
dual<T, Order0, Orders...> operator-(const dual<T, Order0, Orders...> &x) {
    dual<T, Order0, Orders...> res;
    for (size_t i = 0; i <= Order0; ++i) {
        res[i] = -x[i];
    }

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator+(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator+(const dual<T, Orders...> &lhs, const T &rhs) {
    dual<T, Orders...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator+(const T &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = lhs;
    res += rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator-(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = lhs;
    res -= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator-(const dual<T, Orders...> &lhs, const T &rhs) {
    dual<T, Orders...> res = lhs;
    res -= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator*(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<complex<T>, Orders...> operator*(const dual<complex<T>, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    dual<complex<T>, Orders...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<complex<T>, Orders...> operator*(const dual<T, Orders...> &lhs, const dual<complex<T>, Orders...> &rhs) {
    dual<complex<T>, Orders...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator*(const dual<T, Orders...> &lhs, const T &rhs) {
    dual<T, Orders...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<complex<T>, Orders...> operator*(const dual<T, Orders...> &lhs, const complex<T> &rhs) {
    dual<complex<T>, Orders...> res = lhs;
    res *= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator*(const T &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = rhs;
    res *= lhs;

    return res;
}

template <typename T, size_t... Orders>
dual<complex<T>, Orders...> operator*(const complex<T> &lhs, const dual<T, Orders...> &rhs) {
    dual<complex<T>, Orders...> res = rhs;
    res *= lhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator/(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator/(const dual<T, Orders...> &lhs, const T &rhs) {
    dual<T, Orders...> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> operator/(const T &lhs, const dual<T, Orders...> &rhs) {
    dual<T, Orders...> res = lhs;
    res /= rhs;

    return res;
}

template <typename T, size_t... Orders>
bool operator==(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    return lhs.value() == rhs.value();
}

template <typename T, size_t... Orders, typename U>
bool operator==(const dual<T, Orders...> &lhs, const U &rhs) {
    return lhs.value() == rhs;
}

template <typename T, size_t... Orders>
bool operator!=(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    return lhs.value() != rhs.value();
}

template <typename T, size_t... Orders>
bool operator<(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    return lhs.value() < rhs.value();
}

template <typename T, size_t... Orders>
bool operator>(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    return lhs.value() > rhs.value();
}

template <typename T, size_t... Orders>
bool operator<=(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    return lhs.value() <= rhs.value();
}

template <typename T, size_t... Orders>
bool operator>=(const dual<T, Orders...> &lhs, const dual<T, Orders...> &rhs) {
    return lhs.value() >= rhs.value();
}

template <size_t Order, typename T>
dual<T, Order> dual_var(T value, size_t dim = 0) {
    // dim must be zero

    if constexpr (Order >= 1) {
        return {value, 1};
    }

    return value;
}

template <size_t Orders0, size_t Orders1, size_t... Orders, typename T>
dual<T, Orders0, Orders1, Orders...> dual_var(T value, size_t dim = 0) {
    if (dim == 0) {
        if constexpr (Orders0 >= 1) {
            return {value, 1};
        }

        return value;
    }

    return dual_var<Orders1, Orders...>(value, dim - 1);
}

template <typename T, size_t N, size_t... Orders>
dual<T, Orders...> dual_taylor_series(const T (&coef)[N], const dual<T, Orders...> &x, T a) {
    dual<T, Orders...> res = coef[0];
    if constexpr (N >= 2) {
        dual<T, Orders...> y = x - a;
        T denom = 1; // factorial

        res += y * coef[1];

        for (size_t i = 2; i < N; ++i) {
            y *= x - a;
            denom *= T(i);

            res += y * coef[i] / denom;
        }
    }

    return res;
}

template <typename T, size_t... Orders>
dual<T, Orders...> abs(dual<T, Orders...> z) {
    if (z.value() < 0) {
        return dual_taylor_series({abs(z.value()), T(-1)}, z, z.value());
    }

    return dual_taylor_series({abs(z.value()), T(1)}, z, z.value());
}

template <typename T, size_t... Orders>
dual<T, Orders...> abs(dual<complex<T>, Orders...> z) {
    return dual_taylor_series({abs(z.value()), real(z.value()) / abs(z.value())}, z, z.value());
}

template <typename T, size_t... Orders>
dual<T, Orders...> exp(const dual<T, Orders...> &x) {
    static constexpr size_t MaxOrder = dual<T, Orders...>::max_order();

    T coef[MaxOrder + 1] = {exp(x.value())};
    for (size_t i = 1; i <= MaxOrder; ++i) {
        coef[i] = coef[0];
    }

    return dual_taylor_series(coef, x, x.value());
}

template <typename T, size_t... Orders>
dual<T, Orders...> sqrt(const dual<T, Orders...> &z) {
    static constexpr size_t MaxOrder = dual<T, Orders...>::max_order();

    T coef[MaxOrder + 1] = {sqrt(z.value())};
    if constexpr (MaxOrder >= 1) {
        coef[1] = T(1) / (T(2) * coef[0]);

        if constexpr (MaxOrder >= 2) {
            coef[2] = -T(1) / (T(4) * coef[0] * z.value());
        }
    }

    return dual_taylor_series(coef, z, z.value());
}

template <typename T, size_t... Orders>
dual<T, Orders...> sin(const dual<T, Orders...> &x) {
    static constexpr size_t MaxOrder = dual<T, Orders...>::max_order();

    T coef[MaxOrder + 1] = {sin(x.value())};
    if constexpr (MaxOrder >= 1) {
        coef[1] = cos(x.value());

        for (size_t i = 2; i <= MaxOrder; ++i) {
            coef[i] = -coef[i - 2];
        }
    }

    return dual_taylor_series(coef, x, x.value());
}

template <typename T, size_t... Orders>
dual<T, Orders...> cos(const dual<T, Orders...> &x) {
    static constexpr size_t MaxOrder = dual<T, Orders...>::max_order();

    T coef[MaxOrder + 1] = {cos(x.value())};
    if constexpr (MaxOrder >= 1) {
        coef[1] = -sin(x.value());

        for (size_t i = 2; i <= MaxOrder; ++i) {
            coef[i] = -coef[i - 2];
        }
    }

    return dual_taylor_series(coef, x, x.value());
}

template <typename T, size_t... Orders>
dual<T, Orders...> copysign(const dual<T, Orders...> &x, const dual<T, Orders...> &y) {
    if (signbit(x.value()) == signbit(y.value())) {
        return x;
    }

    return -x;
}

template <typename T, size_t... Orders>
bool isfinite(const dual<T, Orders...> &x) {
    return isfinite(x.value());
}

template <typename T, size_t... Orders>
bool isinf(const dual<T, Orders...> &x) {
    return isinf(x.value());
}

template <typename T, size_t... Orders>
bool isnan(const dual<T, Orders...> &x) {
    return isnan(x.value());
}

template <typename T, size_t Order0, size_t... Orders>
dual<T, Order0, Orders...> real(dual<std::complex<T>, Order0, Orders...> z) {
    dual<T, Order0, Orders...> res;
    for (size_t i = 0; i <= Order0; ++i) {
        res[i] = real(z[i]);
    }

    return res;
}

template <typename T, size_t Order0, size_t... Orders>
dual<T, Order0, Orders...> real(dual<T, Order0, Orders...> x) {
    dual<T, Order0, Orders...> res;
    for (size_t i = 0; i <= Order0; ++i) {
        res[i] = real(x[i]);
    }

    return res;
}

template <typename T, size_t Order0, size_t... Orders>
dual<T, Order0, Orders...> imag(dual<std::complex<T>, Order0, Orders...> z) {
    dual<T, Order0, Orders...> res;
    for (size_t i = 0; i <= Order0; ++i) {
        res[i] = imag(z[i]);
    }

    return res;
}

template <typename T, size_t Order0, size_t... Orders>
dual<T, Order0, Orders...> imag(dual<T, Order0, Orders...> x) {
    dual<T, Order0, Orders...> res;
    for (size_t i = 0; i <= Order0; ++i) {
        res[i] = imag(x[i]);
    }

    return res;
}

template <typename T, size_t... Orders>
struct complex_type<dual<T, Orders...>> {
    using type = dual<typename complex_type<T>::type, Orders...>;
};

template <typename T>
struct remove_dual {
    using type = T;
};

template <typename T, size_t... Orders>
struct remove_dual<dual<T, Orders...>> {
    using type = T;
};

template <typename T>
using remove_dual_t = typename remove_dual<T>::type;

namespace numbers {

    template <typename T, size_t... Orders>
    dual<std::complex<T>, Orders...> i_v<dual<T, Orders...>> = i_v<T>;

} // namespace numbers
} // namespace xsf
