// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "../fast_matrix_market.hpp"

namespace fast_matrix_market {
#if __cplusplus >= 202002L
    // If available, use C++20 concepts for programmer clarity.
    // This shows what fast_matrix_market expects each template type to support.

    /**
     * We can read a matrix market file into any vector type that can be resized and iterated.
     * Classic example is std::vector<>.
     */
    template <typename VEC>
    concept triplet_read_vector = requires(VEC v) {
        v.resize(1);
        v.begin();
    };

    /**
     * We can write a matrix market file from any vector type that can be iterated.
     * Classic example is std::vector<>.
     */
    template <typename VEC>
    concept triplet_write_vector = requires(VEC v) {
        v.cbegin();
        v.cend();
    };
#else
    // Make everything still work if concepts are not available
#define triplet_read_vector typename
#define triplet_write_vector typename
#endif

    /**
     * Read a Matrix Market file into a triplet (i.e. row, column, value vectors).
     */
    template <triplet_read_vector IVEC, triplet_read_vector VVEC>
    void read_matrix_market_triplet(std::istream &instream,
                                    matrix_market_header& header,
                                    IVEC& rows, IVEC& cols, VVEC& values,
                                    const read_options& options = {}) {
        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;

        read_header(instream, header);

        rows.resize(get_storage_nnz(header, options));
        cols.resize(get_storage_nnz(header, options));
        values.resize(get_storage_nnz(header, options));

        auto handler = triplet_parse_handler(rows.begin(), cols.begin(), values.begin());
        read_matrix_market_body(instream, header, handler, pattern_default_value((const VT*)nullptr), options);
    }

    /**
     * Convenience method that omits the header requirement if the user only cares about the dimensions.
     */
    template <triplet_read_vector IVEC, triplet_read_vector VVEC, typename DIM>
    void read_matrix_market_triplet(std::istream &instream,
                                    DIM& nrows, DIM& ncols,
                                    IVEC& rows, IVEC& cols, VVEC& values,
                                    const read_options& options = {}) {
        matrix_market_header header;
        read_matrix_market_triplet(instream, header, rows, cols, values, options);
        nrows = header.nrows;
        ncols = header.ncols;
    }

    /**
     * Write triplets to a Matrix Market file.
     */
    template <triplet_write_vector IVEC, triplet_write_vector VVEC>
    void write_matrix_market_triplet(std::ostream &os,
                                     matrix_market_header header,
                                     const IVEC& rows,
                                     const IVEC& cols,
                                     const VVEC& values,
                                     const write_options& options = {}) {
        using IT = typename std::iterator_traits<decltype(rows.begin())>::value_type;
        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;

        header.nnz = values.size();

        header.object = matrix;
        if (header.field != pattern) {
            header.field = get_field_type((const VT *) nullptr);
        }
        header.format = coordinate;

        write_header(os, header);

        line_formatter<IT, VT> lf(header, options);
        auto formatter = triplet_formatter(lf,
                                           rows.cbegin(), rows.cend(),
                                           cols.cbegin(), cols.cend(),
                                           values.cbegin(), header.field == pattern ? values.cbegin() : values.cend());
        write_body(os, formatter, options);
    }

#if __cplusplus < 202002L
    // clean up after ourselves
#undef triplet_read_vector
#undef triplet_write_vector
#endif
}