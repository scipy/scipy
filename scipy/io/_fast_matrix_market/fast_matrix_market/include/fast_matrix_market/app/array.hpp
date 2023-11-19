// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "../fast_matrix_market.hpp"

namespace fast_matrix_market {

#if __cplusplus >= 202002L
    // If available, use C++20 concepts for programmer clarity and better error messages.
    // If not using C++20 this shows what fast_matrix_market expects each template type to support.

    /**
     * We can read a matrix market file into any vector type that can be resized and iterated.
     * Classic example is std::vector<>.
     */
    template <typename VEC>
    concept array_read_vector = requires(VEC v) {
        v.empty();
        v.resize(1);  // assumption is that the new elements are default constructed
        v.begin();
    };

    /**
     * We can write a matrix market file from any vector type that can report its size (for error checking) and
     * can be iterated.
     * Classic example is std::vector<>.
     */
    template <typename VEC>
    concept array_write_vector = requires(VEC v) {
        v.size();
        v.begin();
    };
#else
    // Make everything still work if concepts are not available
#define array_read_vector typename
#define array_write_vector typename
#endif

    /**
     * Read a Matrix Market file into an array.
     */
    template <array_read_vector VEC>
    void read_matrix_market_array(std::istream &instream,
                                  matrix_market_header& header,
                                  VEC& values,
                                  storage_order order = row_major,
                                  const read_options& options = {}) {
        read_header(instream, header);

        if (!values.empty()) {
            values.resize(0);
        }
        values.resize(header.nrows * header.ncols);

        auto handler = dense_adding_parse_handler(values.begin(), order, header.nrows, header.ncols);
        read_matrix_market_body(instream, header, handler, 1, options);
    }

    /**
     * Convenience method that omits the header requirement if the user only cares about the dimensions.
     */
    template <array_read_vector VEC, typename DIM>
    void read_matrix_market_array(std::istream &instream,
                                  DIM& nrows, DIM& ncols,
                                  VEC& values,
                                  storage_order order = row_major,
                                  const read_options& options = {}) {
        matrix_market_header header;
        read_matrix_market_array(instream, header, values, order, options);
        nrows = header.nrows;
        ncols = header.ncols;
    }

    /**
     * Convenience method that omits the header requirement if the user only cares about the values
     * (e.g. loading a 1D vector, where the std::vector length already includes the length).
     */
    template <array_read_vector VEC>
    void read_matrix_market_array(std::istream &instream,
                                  VEC& values,
                                  storage_order order = row_major,
                                  const read_options& options = {}) {
        matrix_market_header header;
        read_matrix_market_array(instream, header, values, order, options);
    }

    /**
     * Write an array to a Matrix Market file.
     */
    template <array_write_vector VEC>
    void write_matrix_market_array(std::ostream &os,
                                   matrix_market_header header,
                                   const VEC& values,
                                   storage_order order = row_major,
                                   const write_options& options = {}) {
        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;

        if (header.nrows * header.ncols != (int64_t)values.size()) {
            throw invalid_argument("Array length does not match matrix dimensions.");
        }

        header.nnz = values.size();

        header.object = matrix;
        if (options.fill_header_field_type) {
            header.field = get_field_type((const VT *) nullptr);
        }
        header.format = array;
        header.symmetry = general;

        write_header(os, header, options);

        line_formatter<int64_t, VT> lf(header, options);
        auto formatter = array_formatter(lf, values.begin(), order, header.nrows, header.ncols);
        write_body(os, formatter, options);
    }

#if __cplusplus < 202002L
// clean up after ourselves
#undef array_read_vector
#undef array_write_vector
#endif
}