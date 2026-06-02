// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <algorithm>
#include <utility>

#include "fast_matrix_market.hpp"

namespace fast_matrix_market {

    /**
     * Format individual lines (matrix version).
     */
    template <typename IT, typename VT>
    class line_formatter {
    public:
        line_formatter(const matrix_market_header &header, const write_options &options) : header(header),
                                                                                           options(options) {}

        std::string coord_matrix(const IT& row, const IT& col, const VT& val) {
            if (header.format == array) {
                return array_matrix(row, col, val);
            }

            std::string line{};
            line += int_to_string(row + 1);
            line += kSpace;
            line += int_to_string(col + 1);

            if (header.field != pattern) {
                line += kSpace;
                line += value_to_string(val, options.precision);
            }
            line += kNewline;

            return line;
        }

        std::string coord_matrix_pattern(const IT& row, const IT& col) {
            std::string line{};
            line += int_to_string(row + 1);
            line += kSpace;
            line += int_to_string(col + 1);
            line += kNewline;

            return line;
        }

        std::string array_matrix(const IT& row, const IT& col, const VT& val) {
            if (header.symmetry != general) {
                if (row < col) {
                    // omit upper triangle
                    return {};
                }
                if (header.symmetry == skew_symmetric && row == col) {
                    // omit diagonal for skew-symmetric
                    return {};
                }
            }

            std::string ret = value_to_string(val, options.precision);
            ret += kNewline;
            return ret;
        }
    protected:
        const matrix_market_header& header;
        const write_options& options;
    };

    /**
     * Format individual lines (vector version).
     */
    template <typename IT, typename VT>
    class vector_line_formatter {
    public:
        vector_line_formatter(const matrix_market_header &header, const write_options &options) : header(header),
                                                                                                  options(options) {}
        std::string coord_matrix(const IT& row, [[maybe_unused]] const IT& col, const VT& val) {
            std::string line{};
            line += int_to_string(row + 1);

            if (header.field != pattern) {
                line += kSpace;
                line += value_to_string(val, options.precision);
            }
            line += kNewline;

            return line;
        }

        std::string coord_matrix_pattern(const IT& row, [[maybe_unused]] const IT& col) {
            std::string line{};
            line += int_to_string(row + 1);
            line += kNewline;

            return line;
        }

    protected:
        const matrix_market_header& header;
        const write_options& options;
    };

    /**
     * Format row, column, value vectors.
     *
     * Value range may be empty (i.e. val_begin == val_end) to omit writing values at all. Useful for pattern matrices.
     *
     * @tparam A_ITER
     * @tparam B_ITER
     * @tparam C_ITER Must be a valid iterator, but if begin==end then the values are not written.
     * @tparam COLUMN_IS_VALUE
     */
    template<typename LF, typename A_ITER, typename B_ITER, typename C_ITER>
    class triplet_formatter {
    public:
        explicit triplet_formatter(LF lf,
                                   const A_ITER row_begin, const A_ITER row_end,
                                   const B_ITER col_begin, const B_ITER col_end,
                                   const C_ITER val_begin, const C_ITER val_end) :
                                   line_formatter(lf),
                                   row_iter(row_begin), row_end(row_end),
                                   col_iter(col_begin),
                                   val_iter(val_begin), val_end(val_end) {
            if (row_end - row_begin != col_end - col_begin ||
                    (row_end - row_begin != val_end - val_begin && val_end != val_begin)) {
                throw invalid_argument("Row, column, and value ranges must have equal length.");
            }
        }

        [[nodiscard]] bool has_next() const {
            return row_iter != row_end;
        }

        class chunk {
        public:
            explicit chunk(LF lf,
                           const A_ITER row_begin, const A_ITER row_end,
                           const B_ITER col_begin,
                           const C_ITER val_begin, const C_ITER val_end) :
                    line_formatter(lf),
                    row_iter(row_begin), row_end(row_end),
                    col_iter(col_begin),
                    val_iter(val_begin), val_end(val_end) {}

            std::string operator()() {
                std::string chunk;
                chunk.reserve((row_end - row_iter)*25);

                for (; row_iter != row_end; ++row_iter, ++col_iter) {
                    if (val_iter != val_end) {
                        chunk += line_formatter.coord_matrix(*row_iter, *col_iter, *val_iter);
                        ++val_iter;
                    } else {
                        chunk += line_formatter.coord_matrix_pattern(*row_iter, *col_iter);
                    }
                }

                return chunk;
            }

            LF line_formatter;
            A_ITER row_iter, row_end;
            B_ITER col_iter;
            C_ITER val_iter, val_end;
        };

        chunk next_chunk(const write_options& options) {
            auto chunk_size = std::min(options.chunk_size_values, (int64_t)(row_end - row_iter));
            A_ITER row_chunk_end = row_iter + chunk_size;
            B_ITER col_chunk_end = col_iter + chunk_size;
            C_ITER val_chunk_end = (val_iter != val_end) ? val_iter + chunk_size: val_end;

            chunk c(line_formatter,
                    row_iter, row_chunk_end,
                    col_iter,
                    val_iter, val_chunk_end);

            row_iter = row_chunk_end;
            col_iter = col_chunk_end;
            val_iter = val_chunk_end;

            return c;
        }

    protected:
        LF line_formatter;
        A_ITER row_iter, row_end;
        B_ITER col_iter;
        C_ITER val_iter, val_end;
    };

    /**
     * Format CSC structures.
     */
    template<typename LF, typename PTR_ITER, typename IND_ITER, typename VAL_ITER>
    class csc_formatter {
    public:
        explicit csc_formatter(LF lf,
                               const PTR_ITER ptr_begin, const PTR_ITER ptr_end,
                               const IND_ITER ind_begin, const IND_ITER ind_end,
                               const VAL_ITER val_begin, const VAL_ITER val_end,
                               bool transpose = false) :
                line_formatter(lf),
                ptr_begin(ptr_begin), ptr_iter(ptr_begin), ptr_end(ptr_end),
                ind_begin(ind_begin),
                val_begin(val_begin), val_end(val_end),
                transpose(transpose) {
            if (ind_end - ind_begin != val_end - val_begin && val_end != val_begin) {
                throw invalid_argument("Index and value ranges must have equal length.");
            }

            auto num_columns = (ptr_end - ptr_iter);
            auto nnz = (ind_end - ind_begin);
            nnz_per_column = ((double)nnz) / num_columns;
        }

        [[nodiscard]] bool has_next() const {
            return ptr_iter != ptr_end;
        }

        class chunk {
        public:
            explicit chunk(LF lf,
                           const PTR_ITER ptr_begin, const PTR_ITER ptr_iter, const PTR_ITER ptr_end,
                           const IND_ITER ind_begin,
                           const VAL_ITER val_begin, const VAL_ITER val_end,
                           bool transpose) :
                    line_formatter(lf),
                    ptr_begin(ptr_begin), ptr_iter(ptr_iter), ptr_end(ptr_end),
                    ind_begin(ind_begin),
                    val_begin(val_begin), val_end(val_end),
                    transpose(transpose) {}

            std::string operator()() {
                std::string chunk;
                chunk.reserve((ptr_end - ptr_iter)*250);

                // emit the columns [ptr_iter, ptr_end)

                // iterate over assigned columns
                for (; ptr_iter != ptr_end; ++ptr_iter) {
                    auto column_number = (int64_t)(ptr_iter - ptr_begin);

                    // iterate over rows in column
                    IND_ITER row_end = ind_begin + *(ptr_iter+1);
                    IND_ITER row_iter = ind_begin + *ptr_iter;
                    VAL_ITER val_iter = val_begin;
                    if (val_begin != val_end) {
                        val_iter = val_begin + *ptr_iter;
                    }
                    for (; row_iter != row_end; ++row_iter) {

                        int64_t lf_row = *row_iter;
                        int64_t lf_col = column_number;
                        if (transpose) {
                            std::swap(lf_row, lf_col);
                        }

                        if (val_iter != val_end) {
                            chunk += line_formatter.coord_matrix(lf_row, lf_col, *val_iter);
                            ++val_iter;
                        } else {
                            chunk += line_formatter.coord_matrix_pattern(lf_row, lf_col);
                        }
                    }
                }

                return chunk;
            }

            LF line_formatter;
            PTR_ITER ptr_begin, ptr_iter, ptr_end;
            IND_ITER ind_begin;
            VAL_ITER val_begin, val_end;
            bool transpose;
        };

        chunk next_chunk(const write_options& options) {
            auto num_columns = (int64_t)(((double)options.chunk_size_values / nnz_per_column) + 1);

            num_columns = std::min(num_columns, (int64_t)(ptr_end - ptr_iter));
            PTR_ITER ptr_chunk_end = ptr_iter + num_columns;

            chunk c(line_formatter,
                    ptr_begin, ptr_iter, ptr_chunk_end,
                    ind_begin,
                    val_begin, val_end,
                    transpose);

            ptr_iter = ptr_chunk_end;

            return c;
        }

    protected:
        LF line_formatter;
        PTR_ITER ptr_begin, ptr_iter, ptr_end;
        IND_ITER ind_begin;
        VAL_ITER val_begin, val_end;
        bool transpose;
        double nnz_per_column;
    };

    /**
     * Format dense arrays.
     */
    template<typename LF, typename VT_ITER>
    class array_formatter {
    public:
        explicit array_formatter(LF lf, const VT_ITER& values, storage_order order, int64_t nrows, int64_t ncols) :
                line_formatter(lf), values(values), order(order), nrows(nrows), ncols(ncols) {}

        [[nodiscard]] bool has_next() const {
            return cur_col != ncols;
        }

        class chunk {
        public:
            explicit chunk(LF lf, const VT_ITER& values, storage_order order, int64_t nrows, int64_t ncols, int64_t cur_col) :
                    line_formatter(lf), values(values), order(order), nrows(nrows), ncols(ncols), cur_col(cur_col) {}

            std::string operator()() {
                std::string c;
                c.reserve(ncols * 15);

                for (int64_t row = 0; row < nrows; ++row) {
                    int64_t offset;
                    if (order == row_major) {
                        offset = row * ncols + cur_col;
                    } else {
                        offset = cur_col * nrows + row;
                    }

                    c += line_formatter.array_matrix(row, cur_col, *(values + offset));
                }

                return c;
            }

            LF line_formatter;
            const VT_ITER values;
            storage_order order;
            int64_t nrows, ncols;
            int64_t cur_col;
        };

        chunk next_chunk([[maybe_unused]] const write_options& options) {
            return chunk(line_formatter, values, order, nrows, ncols, cur_col++);
        }

    protected:
        LF line_formatter;
        const VT_ITER values;
        storage_order order;
        int64_t nrows, ncols;
        int64_t cur_col = 0;
    };

    /**
     * Formats any structure that has:
     * operator(row, col) - returns the value at (row, col)
     *
     * Includes Eigen Dense Matrix/Vector and NumPy arrays.
     */
    template<typename LF, typename DenseType, typename DIM>
    class dense_2d_call_formatter {
    public:
        explicit dense_2d_call_formatter(LF lf, const DenseType& mat, DIM nrows, DIM ncols) :
        line_formatter(lf), mat(mat), nrows(nrows), ncols(ncols) {}

        [[nodiscard]] bool has_next() const {
            return col_iter < ncols;
        }

        class chunk {
        public:
            explicit chunk(LF lf, const DenseType& mat, DIM nrows, DIM col_iter, DIM col_end) :
                    line_formatter(lf), mat(mat), nrows(nrows), col_iter(col_iter), col_end(col_end) {}

            std::string operator()() {
                std::string chunk;
                chunk.reserve((col_end - col_iter) * nrows * 15);

                // iterate over assigned columns
                for (; col_iter != col_end; ++col_iter) {

                    for (DIM row = 0; row < nrows; ++row)
                    {
                        chunk += line_formatter.array_matrix(row, col_iter, mat(row, col_iter));
                    }
                }

                return chunk;
            }

            LF line_formatter;
            const DenseType& mat;
            DIM nrows;
            DIM col_iter, col_end;
        };

        chunk next_chunk(const write_options& options) {
            auto num_columns = (DIM)((double)options.chunk_size_values / nrows) + 1;
            num_columns = std::min(num_columns, ncols - col_iter);

            DIM col_end = col_iter + num_columns;
            chunk c(line_formatter, mat, nrows, col_iter, col_end);
            col_iter = col_end;

            return c;
        }

    protected:
        LF line_formatter;
        const DenseType& mat;
        DIM nrows, ncols;
        DIM col_iter = 0;
    };
}
