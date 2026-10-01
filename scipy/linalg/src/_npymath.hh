/*
 * Helpers for dealing with c128/c64 types.
 */
#pragma once

using namespace wrapper;

namespace sp_linalg {

namespace detail {

/*
 * Numeric limits and useful constants for NumPy scalar types.
 */

template<typename T> struct numeric_limits {};

template<>
struct numeric_limits<float>{
    static constexpr float zero = 0.0f;
    static constexpr float one = 1.0f;
    static constexpr float nan = std::numeric_limits<float>::quiet_NaN();
    static constexpr float eps = std::numeric_limits<float>::epsilon();
};

template<>
struct numeric_limits<double>{
    static constexpr double zero = 0.0;
    static constexpr double one = 1.0;
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    static constexpr double eps = std::numeric_limits<double>::epsilon();
};


template<>
struct numeric_limits<c64>{
    static constexpr c64 zero = {0.0f, 0.0f};
    static constexpr c64 one = {1.0f, 0.0f};
    static constexpr c64 nan = {std::numeric_limits<float>::quiet_NaN(),
                                       std::numeric_limits<float>::quiet_NaN()};
};

template<>
struct numeric_limits<c128>{
    static constexpr c128 zero = {0.0, 0.0};
    static constexpr c128 one = {1.0, 0.0};
    static constexpr c128 nan = {std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()};
};

/*
 * Type traits for mapping NumPy types (e.g., c64) to their C++ equivalents.
 *
 * Provides:
 *   - real_type: the underlying real type (float/double)
 *   - value_type: C++ type for operations (float/double/std::complex<float>/std::complex<double>)
 *   - npy_complex_type: the corresponding NumPy complex type
 *   - typenum: NumPy type number (NPY_FLOAT, NPY_DOUBLE, etc.)
 *   - is_complex: boolean indicating whether the type is complex
 */
template<typename T> struct type_traits {};
template<> struct type_traits<float> {
    using real_type = float;
    using value_type = float;
    using npy_complex_type = c64;
    static constexpr int typenum = NPY_FLOAT;
    static constexpr bool is_complex = false;

};
template<> struct type_traits<double> {
    using real_type = double;
    using value_type = double;
    using npy_complex_type = c128;
    static constexpr int typenum = NPY_DOUBLE;
    static constexpr bool is_complex = false;
};
template<> struct type_traits<c64> {
    using real_type = float;
    using value_type = std::complex<float>;
    using npy_complex_type = c64;
    static constexpr int typenum = NPY_COMPLEX64;
    static constexpr bool is_complex = true;
};
template<> struct type_traits<c128> {
    using real_type = double;
    using value_type = std::complex<double>;
    using npy_complex_type = c128;
    static constexpr int typenum = NPY_COMPLEX128;
    static constexpr bool is_complex = true;
};

} // namespace detail

} // namespace sp_linalg
