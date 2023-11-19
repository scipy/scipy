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
     * Generalize symmetry of triplet.
     *
     * Does not duplicate diagonal elements.
     */
    template <typename IVEC, typename VVEC>
    void generalize_symmetry_triplet(IVEC& rows, IVEC& cols, VVEC& values, const symmetry_type& symmetry) {
        if (symmetry == general) {
            return;
        }

        std::size_t num_diagonal_elements = 0;

        // count how many diagonal elements there are (these do not get duplicated)
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (rows[i] == cols[i]) {
                ++num_diagonal_elements;
            }
        }

        // resize vectors
        auto orig_size = rows.size();
        auto new_size = 2*orig_size - num_diagonal_elements;
        rows.resize(new_size);
        cols.resize(new_size);
        values.resize(new_size);

        // fill in the new values
        auto row_iter = rows.begin() + orig_size;
        auto col_iter = cols.begin() + orig_size;
        auto val_iter = values.begin() + orig_size;
        for (std::size_t i = 0; i < orig_size; ++i) {
            if (rows[i] == cols[i]) {
                continue;
            }

            *row_iter = cols[i];
            *col_iter = rows[i];
            *val_iter = get_symmetric_value<typename VVEC::value_type>(values[i], symmetry);

            ++row_iter; ++col_iter; ++val_iter;
        }
    }

    template <triplet_read_vector IVEC, triplet_read_vector VVEC, typename T>
    void read_matrix_market_body_triplet(std::istream &instream,
                                         const matrix_market_header& header,
                                         IVEC& rows, IVEC& cols, VVEC& values,
                                         T pattern_value,
                                         read_options options = {}) {
        bool app_generalize = false;
        if (options.generalize_symmetry && options.generalize_symmetry_app) {
            app_generalize = true;
            options.generalize_symmetry = false;
        }

        auto nnz = get_storage_nnz(header, options);
        rows.resize(nnz);
        cols.resize(nnz);
        values.resize(nnz);

        auto handler = triplet_parse_handler(rows.begin(), cols.begin(), values.begin());
        read_matrix_market_body(instream, header, handler, pattern_value, options);

        if (app_generalize) {
            generalize_symmetry_triplet(rows, cols, values, header.symmetry);
        }
    }

    /**
     * Read a Matrix Market file into a triplet (i.e. row, column, value vectors).
     */
    template <triplet_read_vector IVEC, triplet_read_vector VVEC>
    void read_matrix_market_triplet(std::istream &instream,
                                    matrix_market_header& header,
                                    IVEC& rows, IVEC& cols, VVEC& values,
                                    const read_options& options = {}) {
        read_header(instream, header);

        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;
        read_matrix_market_body_triplet(instream, header, rows, cols, values, pattern_default_value((const VT*)nullptr), options);
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

        header.nnz = rows.size();

        header.object = matrix;
        if (header.nnz > 0 && (values.cbegin() == values.cend())) {
            header.field = pattern;
        } else if (header.field != pattern && options.fill_header_field_type) {
            header.field = get_field_type((const VT *) nullptr);
        }
        header.format = coordinate;

        write_header(os, header, options);

        line_formatter<IT, VT> lf(header, options);
        auto formatter = triplet_formatter(lf,
                                           rows.cbegin(), rows.cend(),
                                           cols.cbegin(), cols.cend(),
                                           values.cbegin(), header.field == pattern ? values.cbegin() : values.cend());
        write_body(os, formatter, options);
    }

    /**
     * Write CSC/CSR to a Matrix Market file.
     */
    template <triplet_write_vector IVEC, triplet_write_vector VVEC>
    void write_matrix_market_csc(std::ostream &os,
                                 matrix_market_header header,
                                 const IVEC& indptr,
                                 const IVEC& indices,
                                 const VVEC& values,
                                 bool is_csr,
                                 const write_options& options = {}) {
        using IT = typename std::iterator_traits<decltype(indptr.begin())>::value_type;
        using VT = typename std::iterator_traits<decltype(values.begin())>::value_type;

        header.nnz = indices.size();

        header.object = matrix;
        if (header.nnz > 0 && (values.cbegin() == values.cend())) {
            header.field = pattern;
        } else if (header.field != pattern && options.fill_header_field_type) {
            header.field = get_field_type((const VT *) nullptr);
        }
        header.format = coordinate;

        write_header(os, header, options);

        line_formatter<IT, VT> lf(header, options);
        auto formatter = csc_formatter(lf,
                                       indptr.cbegin(), indptr.cend() - 1,
                                       indices.cbegin(), indices.cend(),
                                       values.cbegin(), header.field == pattern ? values.cbegin() : values.cend(),
                                       is_csr);
        write_body(os, formatter, options);
    }

#if __cplusplus < 202002L
    // clean up after ourselves
#undef triplet_read_vector
#undef triplet_write_vector
#endif
}