// Copyright (C) 2022-2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

/*
 * Support std::string as a value type.
 *
 * Supports such things as a triplet where the value type is std::vector<std::string>, and such structure can read
 * any Matrix Market file regardless of field type.
 *
 * IMPORTANT!
 *
 * If using this as a template for your own user type, note the weird #include ordering requirements!
 *
 * ```
 * #include <fast_matrix_market/types.hpp>
 *
 * namespace fast_matrix_market {
 *   // declare your type stuff here
 * }
 *
 * // Now include the main header.
 * #include <fast_matrix_market/fast_matrix_market.hpp>
 * ```
 */

namespace fast_matrix_market {
    /**
     * (Needed for read) Parse a value.
     *
     * This method must find the end of the value, too. Likely the next newline or end of string is the end.
     *
     * @param pos starting character
     * @param end end of string. Do not dereference >end.
     * @param out out parameter of parsed value
     * @return a pointer to the next character following the value. Likely the newline.
     */
    inline const char *read_value(const char *pos, const char *end, std::string &out, [[maybe_unused]] const read_options& options = {}) {
        const char *field_start = pos;
        // find the end of the line
        while (pos != end && *pos != '\n') {
            ++pos;
        }
        out = std::string(field_start, (pos - field_start));

        return pos;
    }

    /**
     * (Optional, used for read) Declare that read_value(std::string) declared above can read complex values too.
     */
    template<> struct can_read_complex<std::string> : std::true_type {};

    /**
     * (Needed for read) Used to handle skew-symmetric files. Must be defined regardless.
     */
    inline std::string negate(const std::string& o) {
        return "-" + o;
    }

    /**
     * (Needed for read) Default value to set for patterns (if option selected). Must be defined regardless.
     */
    inline std::string pattern_default_value([[maybe_unused]] const std::string* type) {
        return "";
    }

    // If using default dense array loader, type must also work with std::plus<>.

    /**
     * (Needed for write) Used to determine what field type to set in the MatrixMarket header.
     *
     * IMPORTANT! This is the default type. If it does not match your actual data, set the appropriate type
     * in the header, then use write_options::fill_header_field_type = false.
     */
    inline field_type get_field_type([[maybe_unused]] const std::string* type) {
        return real;
    }

    /**
     * (Needed for write) Write a value.
     */
    inline std::string value_to_string(const std::string& value, [[maybe_unused]] int precision) {
        return value;
    }
}
