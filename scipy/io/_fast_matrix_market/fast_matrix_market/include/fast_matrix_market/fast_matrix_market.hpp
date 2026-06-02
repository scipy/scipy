// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#ifndef FAST_MATRIX_MARKET_H
#define FAST_MATRIX_MARKET_H

#pragma once

#include <algorithm>
#include <future>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "types.hpp"

// Support std::string as a user type
#include "app/user_type_string.hpp"

namespace fast_matrix_market {

    // Version macros.
    // Keep in sync with python/pyproject.toml
#define FAST_MATRIX_MARKET_VERSION_MAJOR 1
#define FAST_MATRIX_MARKET_VERSION_MINOR 7
#define FAST_MATRIX_MARKET_VERSION_PATCH 4

    constexpr std::string_view kSpace = " ";
    constexpr std::string_view kNewline = "\n";

    /**
     *
     */
    class fmm_error : public std::exception {
    public:
        explicit fmm_error(std::string msg): msg(std::move(msg)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return msg.c_str();
        }
    protected:
        std::string msg;
    };

    /**
     * The provided stream does not represent a Matrix Market file.
     */
    class invalid_mm : public fmm_error {
    public:
        explicit invalid_mm(std::string msg): fmm_error(std::move(msg)) {}
        explicit invalid_mm(std::string msg, int64_t line_num) : fmm_error(std::move(msg)) {
            prepend_line_number(line_num);
        }

        void prepend_line_number(int64_t line_num) {
            msg = std::string("Line ") + std::to_string(line_num) + ": " + msg;
        }
    };

    /**
     * A value was encountered that does not fit the provided type.
     */
    class out_of_range : public invalid_mm {
    public:
        explicit out_of_range(std::string msg): invalid_mm(std::move(msg)) {}
    };

    /**
     * Passed in argument was not valid.
     */
    class invalid_argument : public fmm_error {
    public:
        explicit invalid_argument(std::string msg): fmm_error(std::move(msg)) {}
    };

    /**
     * Matrix Market file has complex fields but the datastructure to load into cannot handle complex values.
     */
    class complex_incompatible : public invalid_argument {
    public:
        explicit complex_incompatible(std::string msg): invalid_argument(std::move(msg)) {}
    };

    /**
     * A Matrix Market file requires a feature that has been disabled via compilation flags.
     */
    class support_not_selected : public invalid_argument {
    public:
        explicit support_not_selected(std::string msg): invalid_argument(std::move(msg)) {}
    };

    /**
     * Matrix Market file is a `vector` type, but vector support is disabled in this build.
     */
    class no_vector_support : public support_not_selected {
    public:
        explicit no_vector_support(std::string msg): support_not_selected(std::move(msg)) {}
    };

    /**
     * A value type to use for pattern matrices. Pattern Matrix Market files do not write a value column, only the
     * coordinates. Setting this as the value type signals the parser to not attempt to read a column that isn't there.
     */
    struct pattern_placeholder_type {};

    /**
     * Negation of a pattern_placeholder_type needed to support symmetry generalization.
     * Skew-symmetric symmetry negates values.
     */
    inline pattern_placeholder_type operator-(const pattern_placeholder_type& o) { return o; }

    /**
     * MSVC does not like std::negate<bool>
     */
    inline bool negate(const bool o) {
        return !o;
    }

    inline bool negate(const std::vector<bool>::reference o) {
        return !o;
    }

    template <typename T>
    T negate(const T& o) {
        return std::negate<T>()(o);
    }

    template <typename T>
    T pattern_default_value([[maybe_unused]] const T* type) {
        return 1;
    }

    /**
     * Zero generator for generalize symmetry with ExtraZeroElement.
     */
    template <typename T>
    T get_zero() {
        return {};
    }

    /**
     * Determine if a std::future is ready to return a result, i.e. finished computing.
     * @return true if the future is ready.
     */
    template<typename R>
    bool is_ready(std::future<R> const& f)
    {
        return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    /**
     * @param flags flags bitwise ORed together
     * @param flag flag bit to test for
     * @return true if the flag bit is set in flags, false otherwise
     */
    inline bool test_flag(int flags, int flag) {
        return (flags & flag) == flag;
    }

    inline bool starts_with(const std::string &str, const std::string& prefix) {
        if (prefix.size() > str.size()) {
            return false;
        }
        return std::equal(prefix.begin(), prefix.end(), str.begin());
    }

    inline bool ends_with(const std::string &str, const std::string& suffix) {
        if (suffix.size() > str.size()) {
            return false;
        }
        return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
    }

    /**
     * Trim the whitespace from both ends of a string. Returns a copy.
     */
    inline std::string trim(std::string s) {
        // ltrim
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));

        // rtrim
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());

        return s;
    }

    /**
     * Replace all instances of `from` with `to` in `str`.
     */
    inline std::string replace_all(const std::string& str, const std::string& from, const std::string& to) {
        std::string ret(str);

        if (from.empty())
            return ret;

        std::string::size_type start_pos = 0;
        while ((start_pos = ret.find(from, start_pos)) != std::string::npos) {
            ret.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }

        return ret;
    }
}

#include "field_conv.hpp"
#include "header.hpp"
#include "parse_handlers.hpp"
#include "formatters.hpp"
#include "read_body.hpp"
#include "write_body.hpp"
#include "app/array.hpp"
#include "app/doublet.hpp"
#include "app/triplet.hpp"

#endif
