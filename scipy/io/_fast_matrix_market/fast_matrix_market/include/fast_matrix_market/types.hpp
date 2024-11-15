// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <complex>
#include <map>
#include <cstdint>
#include <string>

namespace fast_matrix_market {

    enum object_type {matrix, vector};
    const std::map<object_type, const std::string> object_map = {
            {matrix, "matrix"},
            {vector, "vector"},
    };

    enum format_type {array, coordinate};
    const std::map<format_type, const std::string> format_map = {
            {array, "array"},
            {coordinate, "coordinate"},
    };

    enum field_type {real, double_, complex, integer, pattern, unsigned_integer};
    const std::map<field_type, const std::string> field_map = {
            {real, "real"},
            {double_, "double"},  // non-standard
            {complex, "complex"},
            {integer, "integer"},
            {pattern, "pattern"},
            {unsigned_integer, "unsigned-integer"}, // SciPy only
    };

    enum symmetry_type {general, symmetric, skew_symmetric, hermitian};
    const std::map<symmetry_type, const std::string> symmetry_map = {
            {general, "general"},
            {symmetric, "symmetric"},
            {skew_symmetric, "skew-symmetric"},
            {hermitian, "hermitian"},
    };

    /**
     * Matrix Market header
     */
    struct matrix_market_header {
        matrix_market_header() = default;
        explicit matrix_market_header(int64_t vector_length) : object(vector), vector_length(vector_length) {}
        matrix_market_header(int64_t nrows, int64_t ncols) : nrows(nrows), ncols(ncols) {}

        object_type object = matrix;
        format_type format = coordinate;
        field_type field = real;
        symmetry_type symmetry = general;

        // Matrix dimensions
        int64_t nrows = 0;
        int64_t ncols = 0;

        // Vector dimensions
        int64_t vector_length = 0;

        // Number of non-zeros for sparse objects
        int64_t nnz = 0;

        // Comment written in the file header
        std::string comment;

        // Number of lines the header takes up. This is populated by read_header().
        int64_t header_line_count = 1;
    };

    enum storage_order {row_major = 1, col_major = 2};
    enum out_of_range_behavior {BestMatch = 1, ThrowOutOfRange = 2};

    struct read_options {
        /**
         * Chunk size for the parsing step, in bytes.
         */
        int64_t chunk_size_bytes = 2 << 20;

        /**
         * If true then any symmetries other than general are expanded out.
         * For any symmetries other than general, only entries in the lower triangular portion need be supplied.
         * symmetric: for (row, column, value), also generate (column, row, value) except if row==column
         * skew-symmetric: for (row, column, value), also generate (column, row, -value) except if row==column
         * hermitian: for (row, column, value), also generate (column, row, complex_conjugate(value)) except if row==column
         */
        bool generalize_symmetry = true;

        /**
         * If true, perform symmetry generalization in the application binding as a post-processing step.
         * If supported by the binding this method can avoid extra diagonal elements.
         * If false or unsupported, diagonals are handled according to `generalize_coordinate_diagnonal_values`.
         */
        bool generalize_symmetry_app = true;

        /**
         * Generalize Symmetry:
         * How to handle a value on the diagonal of a symmetric coordinate matrix.
         *  - DuplicateElement: Duplicate the diagonal element
         *  - ExtraZeroElement: emit a zero along with the diagonal element. The zero will appear first.
         *
         *  The extra cannot simply be omitted because the handlers work by setting already-allocated memory. This
         *  is necessary for efficient parallelization.
         *
         *  This value is ignored if the parse handler has the kAppending flag set. In that case only a single
         *  diagonal element is emitted.
         */
        enum {ExtraZeroElement, DuplicateElement} generalize_coordinate_diagnonal_values = ExtraZeroElement;

        /**
         * Whether parallel implementation is allowed.
         */
        bool parallel_ok = true;

        /**
         * Number of threads to use. 0 means std::thread::hardware_concurrency().
         */
        int num_threads = 0;

        /**
         * How to handle floating-point values that do not fit into their declared type.
         * For example, parsing 1e9999 will
         *  - BestMatch: return Infinity
         *  - ThrowOutOfRange: throw out_of_range exception
         */
        out_of_range_behavior float_out_of_range_behavior = BestMatch;
    };

    struct write_options {
        int64_t chunk_size_values = 2 << 12;

        /**
         * Whether parallel implementation is allowed.
         */
        bool parallel_ok = true;

        /**
         * Number of threads to use. 0 means std::thread::hardware_concurrency().
         */
        int num_threads = 0;

        /**
         * Floating-point formatting precision.
         * Placeholder. Currently not used due to the various supported float rendering backends.
         */
        int precision = -1;

        /**
         * Whether to always write a comment line even if comment is empty.
         */
        bool always_comment = false;

        /**
         * Whether to determine header field type based on the supplied datastructure.
         *
         * If true then set header.field using `get_field_type()`. The only exception is
         * if field == pattern then it is left unchanged.
         *
         * Possible reasons to set to false:
         *  - Using a custom type, such as std::string, where `get_field_type()` would return the wrong type
         *  - Writing integer structures as real
         */
        bool fill_header_field_type = true;
    };

    template<class T> struct is_complex : std::false_type {};
    template<class T> struct is_complex<std::complex<T>> : std::true_type {};

    template<class T> struct can_read_complex : is_complex<T> {};
}
