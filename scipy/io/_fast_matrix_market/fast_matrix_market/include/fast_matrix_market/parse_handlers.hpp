// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <charconv>
#include <complex>
#include <functional>
#include <iterator>
#include <type_traits>

#include "fast_matrix_market.hpp"

namespace fast_matrix_market {

    /**
     * This parse handler supports parallelism.
     */
    constexpr int kParallelOk = 1;

    /**
     * Writing to the same (row, column) position a second time will affect the previous value.
     * This means that if there is a possibility of duplicate writes from different threads then this
     * parse handler is unsafe. Example: coordinate file parsed into a dense array. Coordinate files may have dupes.
     */
    constexpr int kDense = 2;

    /**
     * This parse handler can handle a variable number of elements. If this flag is not set then the memory
     * has already been allocated and must be filled. If this flag is set, then potentially fewer elements can be
     * written without loss of correctness.
     *
     * This is useful for not needing to duplicate elements on the main diagonal when generalizing symmetry.
     */
    constexpr int kAppending = 4;

    /**
     * Tuple handler. A single vector of (row, column, value) tuples.
     */
    template<typename IT, typename VT, typename ITER>
    class tuple_parse_handler {
    public:
        using coordinate_type = IT;
        using value_type = VT;
        static constexpr int flags = kParallelOk;

        using TUPLE = typename std::iterator_traits<ITER>::value_type;

        explicit tuple_parse_handler(const ITER& iter) : begin_iter(iter), iter(iter) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type value) {
            *iter = TUPLE(row, col, value);
            ++iter;
        }

        tuple_parse_handler<IT, VT, ITER> get_chunk_handler(int64_t offset_from_begin) {
            return tuple_parse_handler(begin_iter + offset_from_begin);
        }

    protected:
        ITER begin_iter;
        ITER iter;
    };

    /**
     * Triplet handler. Separate row, column, value iterators.
     */
    template<typename IT_ITER, typename VT_ITER>
    class triplet_parse_handler {
    public:
        using coordinate_type = typename std::iterator_traits<IT_ITER>::value_type;
        using value_type = typename std::iterator_traits<VT_ITER>::value_type;
        static constexpr int flags = kParallelOk;

        explicit triplet_parse_handler(const IT_ITER& rows,
                                       const IT_ITER& cols,
                                       const VT_ITER& values) : begin_rows(rows), begin_cols(cols), begin_values(values),
                                                                rows(rows), cols(cols), values(values) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type value) {
            *rows = row;
            *cols = col;
            *values = value;

            ++rows;
            ++cols;
            ++values;
        }

        void handle(const coordinate_type row, const coordinate_type col, [[maybe_unused]] const pattern_placeholder_type& pat) {
            *rows = row;
            *cols = col;

            ++rows;
            ++cols;
            ++values;
        }

        triplet_parse_handler<IT_ITER, VT_ITER> get_chunk_handler(int64_t offset_from_begin) {
            return triplet_parse_handler(begin_rows + offset_from_begin,
                                         begin_cols + offset_from_begin,
                                         begin_values + offset_from_begin);
        }

    protected:
        IT_ITER begin_rows;
        IT_ITER begin_cols;
        VT_ITER begin_values;

        IT_ITER rows;
        IT_ITER cols;
        VT_ITER values;
    };

    /**
     * Triplet handler for pattern matrices. Row and column vectors only.
     */
    template<typename IT_ITER>
    class triplet_pattern_parse_handler {
    public:
        using coordinate_type = typename std::iterator_traits<IT_ITER>::value_type;
        using value_type = pattern_placeholder_type;
        static constexpr int flags = kParallelOk;

        explicit triplet_pattern_parse_handler(const IT_ITER& rows,
                                               const IT_ITER& cols) : begin_rows(rows), begin_cols(cols),
                                                                      rows(rows), cols(cols) {}

        void handle(const coordinate_type row, const coordinate_type col, [[maybe_unused]] const value_type ignored) {
            *rows = row;
            *cols = col;

            ++rows;
            ++cols;
        }

        triplet_pattern_parse_handler<IT_ITER> get_chunk_handler(int64_t offset_from_begin) {
            return triplet_pattern_parse_handler(begin_rows + offset_from_begin,
                                                 begin_cols + offset_from_begin);
        }
    protected:
        IT_ITER begin_rows;
        IT_ITER begin_cols;

        IT_ITER rows;
        IT_ITER cols;
    };

    /**
     * Triplet handler. Separate row, column, value iterators.
     */
    template<typename IT, typename VT, typename IT_ARR, typename VT_ARR>
    class triplet_calling_parse_handler {
    public:
        using coordinate_type = IT;
        using value_type = VT;
        static constexpr int flags = kParallelOk;

        explicit triplet_calling_parse_handler(IT_ARR& rows,
                                               IT_ARR& cols,
                                               VT_ARR& values,
                                               int64_t offset = 0) : rows(rows), cols(cols), values(values), offset(offset) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type value) {
            rows(offset) = row;
            cols(offset) = col;
            values(offset) = value;

            ++offset;
        }

        triplet_calling_parse_handler<IT, VT, IT_ARR, VT_ARR> get_chunk_handler(int64_t offset_from_begin) {
            return triplet_calling_parse_handler(rows, cols, values, offset_from_begin);
        }

    protected:
        IT_ARR& rows;
        IT_ARR& cols;
        VT_ARR& values;

        int64_t offset;
    };

    /**
     * Doublet handler, for a (index, value) sparse vector.
     */
    template<typename IT_ITER, typename VT_ITER>
    class doublet_parse_handler {
    public:
        using coordinate_type = typename std::iterator_traits<IT_ITER>::value_type;
        using value_type = typename std::iterator_traits<VT_ITER>::value_type;
        static constexpr int flags = kParallelOk;

        explicit doublet_parse_handler(const IT_ITER& index,
                                       const VT_ITER& values) : begin_index(index), begin_values(values),
                                                                index(index), values(values) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type value) {
            *index = std::max(row, col);
            *values = value;

            ++index;
            ++values;
        }

        doublet_parse_handler<IT_ITER, VT_ITER> get_chunk_handler(int64_t offset_from_begin) {
            return doublet_parse_handler(begin_index + offset_from_begin,
                                         begin_values + offset_from_begin);
        }
    protected:
        IT_ITER begin_index;
        VT_ITER begin_values;

        IT_ITER index;
        VT_ITER values;
    };

    /**
     * Works with any type where `mat(row, column) += value` works.
     */
    template<typename MAT, typename IT, typename VT>
    class dense_2d_call_adding_parse_handler {
    public:
        using coordinate_type = IT;
        using value_type = VT;
        static constexpr int flags = kParallelOk | kDense;

        explicit dense_2d_call_adding_parse_handler(MAT &mat) : mat(mat) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type value) {
            mat(row, col) += std::plus<value_type>()(mat(row, col), value);
        }

        dense_2d_call_adding_parse_handler<MAT, IT, VT> get_chunk_handler([[maybe_unused]] int64_t offset_from_begin) {
            return *this;
        }

    protected:
        MAT &mat;
    };

    /**
     * Dense array handler (row-major).
     */
    template <typename VT_ITER>
    class dense_adding_parse_handler {
    public:
        using coordinate_type = int64_t;
        using value_type = typename std::iterator_traits<VT_ITER>::value_type;
        static constexpr int flags = kParallelOk | kDense;

        explicit dense_adding_parse_handler(const VT_ITER& values, storage_order order, int64_t nrows, int64_t ncols) :
        values(values), order(order), nrows(nrows), ncols(ncols) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type value) {
            int64_t offset;
            if (order == row_major) {
                offset = row * ncols + col;
            } else {
                offset = col * nrows + row;
            }
            values[offset] = std::plus<value_type>()(values[offset], value);
        }

        dense_adding_parse_handler<VT_ITER> get_chunk_handler([[maybe_unused]] int64_t offset_from_begin) {
            return *this;
        }

    protected:
        VT_ITER values;
        storage_order order;
        int64_t nrows;
        int64_t ncols;
    };
}