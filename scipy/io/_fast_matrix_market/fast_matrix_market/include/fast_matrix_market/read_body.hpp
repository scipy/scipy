// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <functional>
#include "fast_matrix_market.hpp"

#include "chunking.hpp"

namespace fast_matrix_market {

    struct line_counts {
        int64_t file_line = 0;
        int64_t element_num = 0;
    };

    /**
     * A handler wrapper for easily handling pattern matrices. This forwards a fixed value. For example, write 1.0 to
     * double matrices. Avoid using zero.
     */
    template<typename FWD_HANDLER>
    class pattern_parse_adapter {
    public:
        using coordinate_type = typename FWD_HANDLER::coordinate_type;
        using value_type = pattern_placeholder_type;
        static constexpr int flags = FWD_HANDLER::flags;

        explicit pattern_parse_adapter(const FWD_HANDLER &handler, typename FWD_HANDLER::value_type fwd_value) : handler(
                handler), fwd_value(fwd_value) {}

        void handle(const coordinate_type row, const coordinate_type col, [[maybe_unused]] const value_type ignored) {
            handler.handle(row, col, fwd_value);
        }

        pattern_parse_adapter<FWD_HANDLER> get_chunk_handler(int64_t offset_from_start) {
            return pattern_parse_adapter<FWD_HANDLER>(handler.get_chunk_handler(offset_from_start), fwd_value);
        }

    protected:
        FWD_HANDLER handler;
        typename FWD_HANDLER::value_type fwd_value;
    };

    /**
     * A handler wrapper so that real/integer files can be read into std::complex matrices by setting all
     * imaginary parts to zero.
     */
    template<typename COMPLEX_HANDLER>
    class complex_parse_adapter {
    public:
        using coordinate_type = typename COMPLEX_HANDLER::coordinate_type;
        using complex_type = typename COMPLEX_HANDLER::value_type;
        using value_type = typename complex_type::value_type;
        static constexpr int flags = COMPLEX_HANDLER::flags;

        explicit complex_parse_adapter(const COMPLEX_HANDLER &handler) : handler(handler) {}

        void handle(const coordinate_type row, const coordinate_type col, const value_type real) {
            handler.handle(row, col, complex_type(real, 0));
        }

        complex_parse_adapter<COMPLEX_HANDLER> get_chunk_handler(int64_t offset_from_start) {
            return complex_parse_adapter(handler.get_chunk_handler(offset_from_start));
        }

    protected:
        COMPLEX_HANDLER handler;
    };

    ///////////////////////////////////////////////////////////////////
    // Limit bool parallelism
    // vector<bool> is specialized to use a bitfield-like scheme. This means
    // that different elements can share the same bytes, making
    // writes to this container require locking.
    // Instead, disable parallelism for bools.

    template <typename T, typename std::enable_if<std::is_same<T, bool>::value, int>::type = 0>
    bool limit_parallelism_for_value_type(bool) {
        return false;
    }

    template <typename T, typename std::enable_if<!std::is_same<T, bool>::value, int>::type = 0>
    bool limit_parallelism_for_value_type(bool parallelism_selected) {
        return parallelism_selected;
    }

    ///////////////////////////////////////////////////////////////////
    // Chunks
    ///////////////////////////////////////////////////////////////////

    template<typename HANDLER>
    line_counts read_chunk_matrix_coordinate(const std::string &chunk, const matrix_market_header &header,
                                             line_counts line, HANDLER &handler, const read_options &options) {
        const char *pos = chunk.c_str();
        const char *end = pos + chunk.size();

        while (pos != end) {
            try {
                typename HANDLER::coordinate_type row, col;
                typename HANDLER::value_type value;

                pos = skip_spaces_and_newlines(pos, line.file_line);
                if (pos == end) {
                    // empty line
                    break;
                }
                if (line.element_num >= header.nnz) {
                    throw invalid_mm("Too many lines in file (file too long)");
                }

                pos = read_int(pos, end, row);
                pos = skip_spaces(pos);
                pos = read_int(pos, end, col);
                pos = skip_spaces(pos);
                pos = read_value(pos, end, value, options);
                pos = bump_to_next_line(pos, end);

                // validate
                if (row <= 0 || static_cast<int64_t>(row) > header.nrows) {
                    throw invalid_mm("Row index out of bounds");
                }
                if (col <= 0 || static_cast<int64_t>(col) > header.ncols) {
                    throw invalid_mm("Column index out of bounds");
                }

                // Generalize symmetry
                // This appears before the regular handler call for ExtraZeroElement handling.
                if (header.symmetry != general && options.generalize_symmetry) {
                    if (col != row) {
                        switch (header.symmetry) {
                            case symmetric:
                                handler.handle(col - 1, row - 1, value);
                                break;
                            case skew_symmetric:
                                if constexpr (!std::is_unsigned_v<typename HANDLER::value_type>) {
                                    handler.handle(col - 1, row - 1, negate(value));
                                } else {
                                    throw invalid_argument("Cannot load skew-symmetric matrix into unsigned value type.");
                                }
                                break;
                            case hermitian:
                                handler.handle(col - 1, row - 1, complex_conjugate(value));
                                break;
                            case general: break;
                        }
                    } else {
                        if (!test_flag(HANDLER::flags, kAppending)) {
                            switch (options.generalize_coordinate_diagnonal_values) {
                                case read_options::ExtraZeroElement:
                                    handler.handle(row - 1, col - 1, get_zero<typename HANDLER::value_type>());
                                    break;
                                case read_options::DuplicateElement:
                                    handler.handle(row - 1, col - 1, value);
                                    break;
                            }
                        }
                    }
                }

                // Matrix Market is one-based
                handler.handle(row - 1, col - 1, value);

                ++line.file_line;
                ++line.element_num;
            } catch (invalid_mm& inv) {
                inv.prepend_line_number(line.file_line + 1);
                throw;
            }
        }
        return line;
    }

    template<typename HANDLER>
    line_counts read_chunk_vector_coordinate(const std::string &chunk, const matrix_market_header &header,
                                             line_counts line, HANDLER &handler, const read_options &options) {
        const char *pos = chunk.c_str();
        const char *end = pos + chunk.size();

        while (pos != end) {
            try {
                typename HANDLER::coordinate_type row;
                typename HANDLER::value_type value;

                pos = skip_spaces_and_newlines(pos, line.file_line);
                if (pos == end) {
                    // empty line
                    break;
                }
                if (line.element_num >= header.nnz) {
                    throw invalid_mm("Too many lines in file (file too long)");
                }
                pos = read_int(pos, end, row);
                pos = skip_spaces(pos);
                pos = read_value(pos, end, value, options);
                pos = bump_to_next_line(pos, end);

                // validate
                if (row <= 0 || static_cast<int64_t>(row) > header.vector_length) {
                    throw invalid_mm("Vector index out of bounds");
                }

                // Matrix Market is one-based
                handler.handle(row - 1, 0, value);

                ++line.file_line;
                ++line.element_num;
            } catch (invalid_mm& inv) {
                inv.prepend_line_number(line.file_line + 1);
                throw;
            }
        }
        return line;
    }

    template<typename HANDLER>
    line_counts read_chunk_array(const std::string &chunk, const matrix_market_header &header, line_counts line,
                                 HANDLER &handler, const read_options &options,
                                 typename HANDLER::coordinate_type &row,
                                 typename HANDLER::coordinate_type &col) {
        const char *pos = chunk.c_str();
        const char *end = pos + chunk.size();

        if (header.symmetry == skew_symmetric) {
            if (row == 0 && col == 0 && header.nrows > 0) {
                // skew-symmetric matrices have zero diagonals
//                if (test_flag(HANDLER::flags, kDense)) {
//                    handler.handle(row, col, get_zero<typename HANDLER::value_type>());
//                }
                row = 1;
            }
        }

        while (pos != end) {
            try {
                typename HANDLER::value_type value;

                pos = skip_spaces_and_newlines(pos, line.file_line);
                if (pos == end) {
                    // empty line
                    break;
                }
                if (static_cast<int64_t>(col) >= header.ncols) {
                    throw invalid_mm("Too many values in array (file too long)");
                }

                pos = read_value(pos, end, value, options);
                pos = bump_to_next_line(pos, end);

                handler.handle(row, col, value);

                if (row != col && options.generalize_symmetry) {
                    switch (header.symmetry) {
                        case symmetric:
                            handler.handle(col, row, value);
                            break;
                        case skew_symmetric:
                            if constexpr (!std::is_unsigned_v<typename HANDLER::value_type>) {
                                handler.handle(col, row, negate(value));
                            } else {
                                throw invalid_argument("Cannot load skew-symmetric matrix into unsigned value type.");
                            }
                            break;
                        case hermitian:
                            handler.handle(col, row, complex_conjugate(value));
                            break;
                        case general:
                            break;
                    }
                }

                // Matrix Market is column-major, advance down the column
                ++row;
                if (static_cast<int64_t>(row) == header.nrows) {
                    ++col;
                    if (header.symmetry == general) {
                        row = 0;
                    } else {
                        row = col;
                        if (header.symmetry == skew_symmetric) {
                            // skew-symmetric matrices have zero diagonals
//                            if (test_flag(HANDLER::flags, kDense)) {
//                                handler.handle(row, col, get_zero<typename HANDLER::value_type>());
//                            }
                            if (static_cast<int64_t>(row) < header.nrows-1) {
                                ++row;
                            }
                        }
                    }
                }

                ++line.file_line;
                ++line.element_num;
            } catch (invalid_mm& inv) {
                inv.prepend_line_number(line.file_line + 1);
                throw;
            }
        }
        return line;
    }

    ////////////////////////////////////////////////
    // Read Matrix Market body
    // Get chunks from file, read chunks
    ///////////////////////////////////////////////
}

#include "read_body_threads.hpp"

namespace fast_matrix_market {

    template <typename HANDLER>
    line_counts read_coordinate_body_sequential(std::istream& instream, const matrix_market_header& header,
                                                HANDLER& handler, const read_options& options = {}) {
        line_counts lc{header.header_line_count, 0};

        // Read the file in chunks
        while (instream.good()) {
            std::string chunk = get_next_chunk(instream, options);

            // parse the chunk
            if (header.object == matrix) {
                lc = read_chunk_matrix_coordinate(chunk, header, lc, handler, options);
            } else {
                lc = read_chunk_vector_coordinate(chunk, header, lc, handler, options);
            }
        }

        return lc;
    }

    template <typename HANDLER>
    line_counts read_array_body_sequential(std::istream& instream, const matrix_market_header& header,
                                           HANDLER& handler,
                                           const read_options& options = {}) {
        line_counts lc{header.header_line_count, 0};

        typename HANDLER::coordinate_type row = 0;
        typename HANDLER::coordinate_type col = 0;

        // Read the file in chunks
        while (instream.good()) {
            std::string chunk = get_next_chunk(instream, options);

            // parse the chunk
            lc = read_chunk_array(chunk, header, lc, handler, options, row, col);
        }

        return lc;
    }

    /**
     * Read the body with no automatic adaptations.
     */
    template <typename HANDLER>
    void read_matrix_market_body_no_adapters(std::istream& instream, const matrix_market_header& header,
                                             HANDLER& handler, const read_options& options = {}) {
        // Verify generalize symmetry is compatible with this file.
        if (header.symmetry != general && options.generalize_symmetry) {
            if (header.object != matrix) {
                throw invalid_mm("Invalid Symmetry: vectors cannot have symmetry. Set generalize_symmetry=false to disregard this symmetry.");
            }
        }

        line_counts lc;
        bool threads = options.parallel_ok && options.num_threads != 1 && test_flag(HANDLER::flags, kParallelOk);

        threads = limit_parallelism_for_value_type<typename HANDLER::value_type>(threads);

        if (header.symmetry != general && header.format == array) {
            // Parallel array loader does not handle symmetry
            threads = false;
        }

        if (header.format == coordinate && test_flag(HANDLER::flags, kDense)) {
            // Potential race condition if the file contains duplicates.
            threads = false;
        }

        if (threads) {
            lc = read_body_threads(instream, header, handler, options);
        } else {
            if (header.format == coordinate) {
                lc = read_coordinate_body_sequential(instream, header, handler, options);
            } else {
                lc = read_array_body_sequential(instream, header, handler, options);
            }
        }

        // verify the file is not truncated
        if (lc.element_num < header.nnz) {
            if (!(header.symmetry != general && header.format == array)) {
                throw invalid_mm(std::string("Truncated file. Expected another ") +
                                 std::to_string(header.nnz - lc.element_num) + " lines.");
            }
        }
    }

    /**
     * Read the body by adapting real files to complex HANDLER.
     */
    template <typename HANDLER, typename std::enable_if<is_complex<typename HANDLER::value_type>::value, int>::type = 0>
    void read_matrix_market_body_no_pattern(std::istream& instream, const matrix_market_header& header,
                                            HANDLER& handler, const read_options& options = {}) {
        if (header.field == complex) {
            read_matrix_market_body_no_adapters(instream, header, handler, options);
        } else {
            // the handler is expecting std::complex values, but the file is only integer/real
            // provide adapter
            auto fwd_handler = complex_parse_adapter<HANDLER>(handler);
            read_matrix_market_body_no_adapters(instream, header, fwd_handler, options);
        }
    }

    /**
     * Read the body by adapting real files to complex HANDLER.
     */
    template <typename HANDLER, typename std::enable_if<!is_complex<typename HANDLER::value_type>::value, int>::type = 0>
    void read_matrix_market_body_no_pattern(std::istream& instream, const matrix_market_header& header,
                                            HANDLER& handler, const read_options& options = {}) {
        if (header.field != complex) {
            read_matrix_market_body_no_adapters(instream, header, handler, options);
        } else {
            // the file is complex but the values are not
            throw complex_incompatible("Matrix Market file has complex fields but passed data structure cannot handle complex values.");
        }
    }

    /**
     * Main body reader entry point.
     *
     * This will handle the following adaptations automatically:
     *  - If the file is a pattern file, the pattern_value will be substituted for each element
     *  - If the HANDLER expects std::complex values but the file is not complex then imag=0 is provided for each value.
     */
    template <typename HANDLER>
    void read_matrix_market_body(std::istream& instream, const matrix_market_header& header,
                                 HANDLER& handler,
                                 typename HANDLER::value_type pattern_value,
                                 const read_options& options = {}) {
        if (header.field == pattern) {
            auto fwd_handler = pattern_parse_adapter<HANDLER>(handler, pattern_value);
            read_matrix_market_body_no_pattern(instream, header, fwd_handler, options);
        } else {
            read_matrix_market_body_no_pattern(instream, header, handler, options);
        }
    }
}