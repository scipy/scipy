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

    dual &operator+=(const dual &rhs) {
        for (size_t i = 0; i <= N; ++i) {
            values[i] += rhs.values[i];
        }

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

    T value() const { return values[0]; }

    T jac() const { return values[1]; }

    T hess() const { return T(2) * values[2]; }

    decltype(auto) derivatives() const { return derivatives_impl(values); }
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

template <typename T, size_t N, typename U>
dual<T, N> operator*(const dual<T, N> &lhs, const U &rhs) {
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

    //    return {
    //      lhs.value() / rhs.value(), (lhs.jac() * rhs.value() - lhs.value() * rhs.jac()) / (rhs.value() * rhs.value())
    //};
}

template <size_t N, typename T>
dual<T, N> make_dual(T value) {
    return {value, 1};
}

} // namespace xsf
