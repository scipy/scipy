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
    concept doublet_read_vector = requires(VEC v) {
        v.resize(1);
        v.begin();
    };

    /**
     * We can write a matrix market file from any vector type that can be iterated.
     * Classic example is std::vector<>.
     */
    template <typename VEC>
    concept doublet_write_vector = requires(VEC v) {
        v.cbegin();
        v.cend();
    };
#else
    // Make everything still work if concepts are not available
#define doublet_read_vector typename
#define doublet_write_vector typename
#endif

    /**
     * Read a Matrix Market vector file into a doublet sparse vector (i.e. index, value vectors).
     *
     * Any vector-like Matrix Market file will work:
     *  - object=vector file, either dense or sparse
     *  - object=matrix file as long as nrows=1 or ncols=1
     */
    template <doublet_read_vector IVEC, doublet_read_vector VVEC>
    void read_matrix_market_doublet(std::istream &instream,
                                    matrix_market_header& header,
                                    IVEC& indices, VVEC& values,
                                    const read_options& options = {}) {
        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;

        read_header(instream, header);

        indices.resize(header.nnz);
        values.resize(get_storage_nnz(header, options));

        auto handler = doublet_parse_handler(indices.begin(), values.begin());
        read_matrix_market_body(instream, header, handler, pattern_default_value((const VT*)nullptr), options);
    }

    /**
     * Convenience method that omits the header requirement if the user only cares about the dimensions.
     */
    template <doublet_read_vector IVEC, doublet_read_vector VVEC, typename DIM>
    void read_matrix_market_doublet(std::istream &instream,
                                    DIM& length,
                                    IVEC& indices, VVEC& values,
                                    const read_options& options = {}) {
        matrix_market_header header;
        read_matrix_market_doublet(instream, header, indices, values, options);
        length = header.vector_length;
    }

    /**
     * Write doublets to a Matrix Market file.
     */
    template <doublet_write_vector IVEC, doublet_write_vector VVEC>
    void write_matrix_market_doublet(std::ostream &os,
                                     matrix_market_header header,
                                     const IVEC& indices,
                                     const VVEC& values,
                                     const write_options& options = {}) {
        using IT = typename std::iterator_traits<decltype(indices.begin())>::value_type;
        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;

        header.nnz = indices.size();

        header.object = vector;
        if (header.nnz > 0 && (values.cbegin() == values.cend())) {
            header.field = pattern;
        } else if (header.field != pattern) {
            header.field = get_field_type((const VT *) nullptr);
        }
        header.format = coordinate;

        write_header(os, header, options);

        vector_line_formatter<IT, VT> lf(header, options);
        auto formatter = triplet_formatter(lf,
                                          indices.cbegin(), indices.cend(),
                                          indices.cbegin(), indices.cend(),
                                          values.cbegin(), header.field == pattern ? values.cbegin() : values.cend());
        write_body(os, formatter, options);
    }

#if __cplusplus < 202002L
    // clean up after ourselves
#undef doublet_read_vector
#undef doublet_write_vector
#endif
}