// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <charconv>
#include <cmath>
#include <cstring>
#include <complex>
#include <limits>
#include <iomanip>
#include <type_traits>

#ifdef FMM_USE_FAST_FLOAT
#include <fast_float/fast_float.h>
#endif

#ifdef FMM_USE_DRAGONBOX
#include <dragonbox/dragonbox_to_chars.h>
#endif

#ifdef FMM_USE_RYU
#include <ryu/ryu.h>
#endif

#include "fast_matrix_market.hpp"

namespace fast_matrix_market {
    ///////////////////////////////////////////
    // Whitespace management
    ///////////////////////////////////////////

    inline const char* skip_spaces(const char* pos) {
        return pos + std::strspn(pos, " \t\r");
    }

    inline const char* skip_spaces_and_newlines(const char* pos, int64_t& line_num) {
        pos = skip_spaces(pos);
        while (*pos == '\n') {
            ++line_num;
            ++pos;
            pos = skip_spaces(pos);
        }
        return pos;
    }

    inline const char* bump_to_next_line(const char* pos, const char* end) {
        if (pos == end) {
            return pos;
        }

        // find the newline
        pos = std::strchr(pos, '\n');

        // bump to start of next line
        if (pos != end) {
            ++pos;
        }
        return pos;
    }

    ///////////////////////////////////////////
    // Integer / Floating Point Field parsers
    ///////////////////////////////////////////

#ifdef FMM_FROM_CHARS_INT_SUPPORTED
    /**
     * Parse integer using std::from_chars
     */
    template <typename IT>
    const char* read_int_from_chars(const char* pos, const char* end, IT& out) {
        std::from_chars_result result = std::from_chars(pos, end, out);
        if (result.ec != std::errc()) {
            if (result.ec == std::errc::result_out_of_range) {
                throw out_of_range("Integer out of range.");
            } else {
                throw invalid_mm("Invalid integer value.");
            }
        }
        return result.ptr;
    }
#endif

    inline const char* read_int_fallback(const char* pos, [[maybe_unused]] const char* end, long long& out) {
        errno = 0;

        char* value_end;
        out = std::strtoll(pos, &value_end, 10);
        if (errno != 0 || pos == value_end) {
            if (errno == ERANGE) {
                throw out_of_range("Integer out of range.");
            } else {
                throw invalid_mm("Invalid integer value.");
            }
        }

        return value_end;
    }

    inline const char* read_int_fallback(const char* pos, [[maybe_unused]] const char* end, unsigned long long& out) {
        errno = 0;

        char *value_end;
        out = std::strtoull(pos, &value_end, 10);
        if (errno != 0 || pos == value_end) {
            if (errno == ERANGE) {
              throw out_of_range("Integer out of range.");
            } else {
                throw invalid_mm("Invalid integer value.");
            }
        }
        return value_end;
    }

    /**
     * Parse integers using C routines.
     *
     * This is a compatibility fallback.
     */
    template <typename T>
    const char* read_int_fallback(const char* pos, const char* end, T& out) {
        long long i64;
        const char* ret = read_int_fallback(pos, end, i64);

        if (sizeof(T) < sizeof(long long)) {
            if (i64 > (long long) std::numeric_limits<T>::max() ||
                i64 < (long long) std::numeric_limits<T>::min()) {
                throw out_of_range(std::string("Integer out of range."));
            }
        }
        out = static_cast<T>(i64);
        return ret;
    }

    /**
     * Parse integer using best available method
     */
    template <typename IT>
    const char* read_int(const char* pos, const char* end, IT& out) {
#ifdef FMM_FROM_CHARS_INT_SUPPORTED
        return read_int_from_chars(pos, end, out);
#else
        return read_int_fallback(pos, end, out);
#endif
    }

#ifdef FMM_USE_FAST_FLOAT
    /**
     * Parse float or double using fast_float::from_chars
     */
    template <typename FT>
    const char* read_float_fast_float(const char* pos, const char* end, FT& out, out_of_range_behavior oorb) {
        fast_float::from_chars_result result = fast_float::from_chars(pos, end, out, fast_float::chars_format::general);

        if (result.ec != std::errc()) {
            if (result.ec == std::errc::result_out_of_range) {
                if (oorb == ThrowOutOfRange) {
                    throw out_of_range("Floating-point value out of range.");
                }
            } else {
                throw invalid_mm("Invalid floating-point value.");
            }
        }
        return result.ptr;
    }
#endif


#ifdef FMM_FROM_CHARS_DOUBLE_SUPPORTED
    /**
     * Parse float or double using std::from_chars
     */
    template <typename FT>
    const char* read_float_from_chars(const char* pos, const char* end, FT& out, out_of_range_behavior oorb) {
        std::from_chars_result result = std::from_chars(pos, end, out);
        if (result.ec != std::errc()) {
            if (result.ec == std::errc::result_out_of_range) {
                if (oorb == ThrowOutOfRange) {
                    throw out_of_range("Floating-point overflow");
                } else {
                    // std::from_chars does not return a best match on under/overflow, so fall back to strtod
                    out = static_cast<FT>(std::strtod(pos, nullptr));
                }
            } else {
                throw invalid_mm("Invalid floating-point value.");
            }
        }
        return result.ptr;
    }
#endif

    /**
     * Parse double using strtod(). This is a compatibility fallback.
     */
    inline const char* read_float_fallback(const char* pos, [[maybe_unused]] const char* end, double& out, out_of_range_behavior oorb = ThrowOutOfRange) {
        errno = 0;

        char* value_end;
        out = std::strtod(pos, &value_end);
        if (errno != 0 || pos == value_end) {
            if (errno == ERANGE) {
                if (oorb == ThrowOutOfRange) {
                    throw out_of_range("Floating-point value out of range.");
                }
            } else {
                throw invalid_mm("Invalid floating-point value.");
            }
        }
        return value_end;
    }

    /**
     * Parse float using strtof(). This is a compatibility fallback.
     */
    inline const char* read_float_fallback(const char* pos, [[maybe_unused]] const char* end, float& out, out_of_range_behavior oorb) {
        errno = 0;

        char* value_end;
        out = std::strtof(pos, &value_end);
        if (errno != 0 || pos == value_end) {
            if (errno == ERANGE) {
                if (oorb == ThrowOutOfRange) {
                    throw out_of_range("Floating-point value out of range.");
                }
            } else {
                throw invalid_mm("Invalid floating-point value.");
            }
        }
        return value_end;
    }

    template <typename FT>
    const char* read_float(const char* pos, const char* end, FT& out, out_of_range_behavior oorb) {
        constexpr bool have_fast_float =
#ifdef FMM_USE_FAST_FLOAT
        true;
#else
        false;
#endif

        if constexpr (have_fast_float && (std::is_same_v<FT, float> || std::is_same_v<FT, double>)) {
            return read_float_fast_float(pos, end, out, oorb);
        } else {
#if defined(FMM_FROM_CHARS_DOUBLE_SUPPORTED)
            return read_float_from_chars(pos, end, out, oorb);
#else
            return read_float_fallback(pos, end, out, oorb);
#endif
        }
    }

#ifdef FMM_FROM_CHARS_LONG_DOUBLE_SUPPORTED
    /**
     * Parse long double using std::from_chars
     */
    inline const char* read_float_from_chars(const char* pos, const char* end, long double& out, out_of_range_behavior oorb) {
        std::from_chars_result result = std::from_chars(pos, end, out);
        if (result.ec != std::errc()) {
            if (result.ec == std::errc::result_out_of_range) {
                if (oorb == ThrowOutOfRange) {
                    throw out_of_range("Floating-point value out of range.");
                } else {
                    // std::from_chars does not return a best match on under/overflow, so fall back to strtold
                    out = std::strtold(pos, nullptr);
                }
            } else {
                throw invalid_mm("Invalid floating-point value.");
            }
        }
        return result.ptr;
    }
#endif

    /**
     * Parse `long double` using std::strtold().
     *
     * fast_float does not support long double.
     */
    inline const char* read_float_fallback(const char* pos, [[maybe_unused]] const char* end, long double& out, out_of_range_behavior oorb) {
        errno = 0;

        char* value_end;
        out = std::strtold(pos, &value_end);
        if (errno != 0 || pos == value_end) {
            if (errno == ERANGE) {
                if (oorb == ThrowOutOfRange) {
                    throw out_of_range("Floating-point value out of range.");
                }
            } else {
                throw invalid_mm("Invalid floating-point value.");
            }
        }
        return value_end;
    }

    inline const char* read_float(const char* pos, [[maybe_unused]] const char* end, long double& out, out_of_range_behavior oorb) {
#ifdef FMM_FROM_CHARS_LONG_DOUBLE_SUPPORTED
        return read_float_from_chars(pos, end, out, oorb);
#else
        return read_float_fallback(pos, end, out, oorb);
#endif
    }

    //////////////////////////////////////
    // Read value. These evaluate to the field parsers above, depending on requested type
    //////////////////////////////////////

    /**
     * Pattern values are no-ops.
     */
    inline const char* read_value(const char* pos, [[maybe_unused]] const char* end, [[maybe_unused]] pattern_placeholder_type& out, [[maybe_unused]] const read_options& options = {}) {
        return pos;
    }

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    const char* read_value(const char* pos, const char* end, T& out, [[maybe_unused]] const read_options& options = {}) {
        return read_int(pos, end, out);
    }

    inline const char* read_value(const char* pos, const char* end, bool& out, const read_options& options = {}) {
        double parsed;
        auto ret = read_float(pos, end, parsed, options.float_out_of_range_behavior);
        out = (parsed != 0);
        return ret;
    }

    template <typename T, typename std::enable_if<std::is_floating_point_v<T>, int>::type = 0>
    const char* read_value(const char* pos, const char* end, T& out, const read_options& options = {}) {
        return read_float(pos, end, out, options.float_out_of_range_behavior);
    }

    template <typename COMPLEX, typename std::enable_if<is_complex<COMPLEX>::value, int>::type = 0>
    const char* read_value(const char* pos, const char* end, COMPLEX& out, const read_options& options = {}) {
        typename COMPLEX::value_type real, imaginary;
        pos = read_float(pos, end, real, options.float_out_of_range_behavior);
        pos = skip_spaces(pos);
        pos = read_float(pos, end, imaginary, options.float_out_of_range_behavior);

        out.real(real);
        out.imag(imaginary);

        return pos;
    }

    template <typename T, typename std::enable_if<is_complex<T>::value, int>::type = 0>
    T complex_conjugate(const T& value) {
        return T(value.real(), -value.imag());
    }

    template <typename T, typename std::enable_if<!is_complex<T>::value, int>::type = 0>
    T complex_conjugate(const T& value) {
        return value;
    }

    ////////////////////////////////////////////
    // Value to String conversions
    ////////////////////////////////////////////

#ifdef FMM_TO_CHARS_INT_SUPPORTED
    /**
     * Convert integral types to string.
     * std::to_string and std::to_chars has similar performance, however std::to_string is locale dependent and
     * therefore will cause thread serialization.
     */
    template <typename T>
    std::string int_to_string(const T& value) {
        std::string ret(20, ' ');
        std::to_chars_result result = std::to_chars(ret.data(), ret.data() + ret.size(), value);
        if (result.ec == std::errc()) {
            ret.resize(result.ptr - ret.data());
            return ret;
        } else {
            return std::to_string(value);
        }
    }
#else
    /**
     * Convert integral types to string. This is the fallback due to locale dependence (and hence thread serialization).
     */
    template <typename T>
    std::string int_to_string(const T& value) {
        return std::to_string(value);
    }
#endif

    inline std::string value_to_string([[maybe_unused]] const pattern_placeholder_type& value, [[maybe_unused]] int precision) {
        return {};
    }

    inline std::string value_to_string(const bool & value, [[maybe_unused]] int precision) {
        return value ? "1" : "0";
    }

    template <typename T, typename std::enable_if<std::is_integral_v<T>, int>::type = 0>
    std::string value_to_string(const T& value, [[maybe_unused]] int precision) {
        return int_to_string(value);
    }

    /**
     * stdlib fallback
     */
    template <typename T>
    std::string value_to_string_fallback(const T& value, int precision) {
        if (precision < 0) {
            // shortest representation
            if constexpr (std::is_floating_point_v<T>
                            && !std::is_same_v<T, float>
                            && !std::is_same_v<T, double>
                            && !std::is_same_v<T, long double>) {
                // std::to_string is faster but does not support fixed-width floating point types.
                std::ostringstream oss;
                oss << value;
                return oss.str();
            } else {
                return std::to_string(value);
            }
        } else {
            std::ostringstream oss;
            oss << std::setprecision(precision) << value;
            return oss.str();
        }
    }

// Sometimes Dragonbox and Ryu may render '1' as '1E0'
// This controls whether to truncate those suffixes.
#ifndef FMM_DROP_ENDING_E0
#define FMM_DROP_ENDING_E0 1
#endif

// Same as above, but for context where precision is specified
#ifndef FMM_DROP_ENDING_E0_PRECISION
#define FMM_DROP_ENDING_E0_PRECISION 0
#endif

#ifdef FMM_USE_DRAGONBOX

    inline std::string value_to_string_dragonbox(const float& value) {
        std::string buffer(jkj::dragonbox::max_output_string_length<jkj::dragonbox::ieee754_binary32> + 1, ' ');

        char *end_ptr = jkj::dragonbox::to_chars(value, buffer.data());
        buffer.resize(end_ptr - buffer.data());

#if FMM_DROP_ENDING_E0
        if (ends_with(buffer, "E0")) {
            buffer.resize(buffer.size() - 2);
        }
#endif
        return buffer;
    }

    inline std::string value_to_string_dragonbox(const double& value) {
        std::string buffer(jkj::dragonbox::max_output_string_length<jkj::dragonbox::ieee754_binary64> + 1, ' ');

        char *end_ptr = jkj::dragonbox::to_chars(value, buffer.data());
        buffer.resize(end_ptr - buffer.data());

#if FMM_DROP_ENDING_E0
        if (ends_with(buffer, "E0")) {
            buffer.resize(buffer.size() - 2);
        }
#endif
        return buffer;
    }
#endif

#ifdef FMM_USE_RYU
    inline std::string value_to_string_ryu(const float& value, int precision) {
        std::string ret(16, ' ');

        if (precision < 0) {
            // shortest representation
            auto len = f2s_buffered_n(value, ret.data());
            ret.resize(len);

#if FMM_DROP_ENDING_E0
            if (ends_with(ret, "E0")) {
                ret.resize(ret.size() - 2);
            }
#endif
        } else {
            // explicit precision
            if (precision > 0) {
                // d2exp_buffered_n's precision means number of places after the decimal point, but
                // we expect it to mean number of sigfigs.
                --precision;
            }
            auto len = d2exp_buffered_n(static_cast<double>(value), precision, ret.data());
            ret.resize(len);

#if FMM_DROP_ENDING_E0_PRECISION
            if (ends_with(ret, "e+00")) {
                ret.resize(ret.size() - 4);
            }
#endif
        }

        return ret;
    }

    inline std::string value_to_string_ryu(const double& value, int precision) {
        std::string ret(26, ' ');

        if (precision < 0) {
            // shortest representation
            auto len = d2s_buffered_n(value, ret.data());
            ret.resize(len);

#if FMM_DROP_ENDING_E0
            if (ends_with(ret, "E0")) {
                ret.resize(ret.size() - 2);
            }
#endif
        } else {
            // explicit precision
            if (precision > 0) {
                // d2exp_buffered_n's precision means number of places after the decimal point, but
                // we expect it to mean number of sigfigs.
                --precision;
            }
            auto len = d2exp_buffered_n(value, precision, ret.data());
            ret.resize(len);

#if FMM_DROP_ENDING_E0_PRECISION
            if (ends_with(ret, "e+00")) {
                ret.resize(ret.size() - 4);
            }
#endif
        }

        return ret;
    }
#endif

#ifdef FMM_TO_CHARS_DOUBLE_SUPPORTED
    template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
    inline std::string value_to_string_to_chars(const T& value, int precision) {
        std::string ret(100, ' ');
        std::to_chars_result result{};
        if (precision < 0) {
            // shortest representation
            result = std::to_chars(ret.data(), ret.data() + ret.size(), value);
        } else {
            // explicit precision
            result = std::to_chars(ret.data(), ret.data() + ret.size(), value, std::chars_format::general, precision);
        }
        if (result.ec == std::errc()) {
            ret.resize(result.ptr - ret.data());
            return ret;
        } else {
            return value_to_string_fallback(value, precision);
        }
    }

#endif

#ifdef FMM_TO_CHARS_LONG_DOUBLE_SUPPORTED
    inline std::string value_to_string_to_chars(const long double& value, int precision) {
        std::string ret(50, ' ');
        std::to_chars_result result{};
        if (precision < 0) {
            // shortest representation
            result = std::to_chars(ret.data(), ret.data() + ret.size(), value);
        } else {
            // explicit precision
            result = std::to_chars(ret.data(), ret.data() + ret.size(), value, std::chars_format::general, precision);
        }
        if (result.ec == std::errc()) {
            ret.resize(result.ptr - ret.data());
            return ret;
        } else {
            return value_to_string_fallback(value, precision);
        }
    }
#endif

    /**
     * long double to string.
     *
     * Preference order: to_chars, fallback.
     * Note: Ryu's generic_128 can do this on some platforms, but it is not reliable.
     * see https://github.com/ulfjack/ryu/issues/215
     */
    inline std::string value_to_string(const long double& value, int precision) {
#if defined(FMM_TO_CHARS_LONG_DOUBLE_SUPPORTED)
        return value_to_string_to_chars(value, precision);
#else
        return value_to_string_fallback(value, precision);
#endif
    }


    /**
     * floating-point to string.
     *
     * Preference order: Dragonbox (fast but no precision support), to_chars, Ryu, fallback.
     *
     * Dragonbox and Ryu only support float and double types.
     */
    template <typename T, typename std::enable_if<std::is_floating_point_v<T> && !std::is_same_v<T, long double>, int>::type = 0>
    std::string value_to_string(const T& value, int precision) {
#ifdef FMM_USE_DRAGONBOX
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            if (precision < 0) {
                // Shortest representation. Dragonbox is fastest.
                return value_to_string_dragonbox(value);
            }
        }
#endif

#ifdef FMM_TO_CHARS_DOUBLE_SUPPORTED
        return value_to_string_to_chars(value, precision);
#else
        constexpr bool have_ryu =
#ifdef FMM_USE_RYU
            true;
#else
            false;
#endif

        if constexpr (have_ryu && (std::is_same_v<T, float> || std::is_same_v<T, double>)) {
            return value_to_string_ryu(value, precision);
        } else {
            return value_to_string_fallback(value, precision);
        }
#endif
    }

    template <typename COMPLEX, typename std::enable_if<is_complex<COMPLEX>::value, int>::type = 0>
    inline std::string value_to_string(const COMPLEX& value, int precision) {
        return value_to_string(value.real(), precision) + " " + value_to_string(value.imag(), precision);
    }

    /**
     * Catchall
     */
    template <typename T, typename std::enable_if<!std::is_integral_v<T> && !std::is_floating_point_v<T> && !is_complex<T>::value, int>::type = 0>
    std::string value_to_string(const T& value, int precision) {
        return value_to_string_fallback(value, precision);
    }
}