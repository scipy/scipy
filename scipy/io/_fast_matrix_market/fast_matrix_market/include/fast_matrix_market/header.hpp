// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <fast_matrix_market/fast_matrix_market.hpp>

#include "chunking.hpp"

namespace fast_matrix_market {

    /**
     * Matrix Market header starts with this string.
     */
    const std::string kMatrixMarketBanner = "%%MatrixMarket";

    /**
     * Invalid banner, but some packages emit this instead of the double %% version.
     */
    const std::string kMatrixMarketBanner2 = "%MatrixMarket";

    template <typename ENUM>
    ENUM parse_enum(const std::string& s, std::map<ENUM, const std::string> mp) {
        // Make s lowercase for a case-insensitive match
        std::string lower(s);
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        for (const auto& [key, value] : mp) {
            if (value == lower) {
                return key;
            }
        }

        std::string acceptable;
        std::string delim;
        for (const auto& [key, value] : mp) {
            acceptable += delim + std::string(value);
            delim = ", ";
        }
        throw invalid_argument(std::string("Invalid value. Must be one of: ") + acceptable);
    }

    inline bool is_line_all_spaces(const std::string& line) {
        if (line.empty()) {
            return true;
        }

        auto end = std::cend(line);

        if (line[line.size()-1] == '\n') {
            // ignore newline
            --end;
        }

        return is_all_spaces(std::cbegin(line), end);
    }

    inline void strip_trailing_cr(std::string &line) {
        if (!line.empty() && line[line.size() - 1] == '\r') {
            line.resize(line.size() - 1);
        }
    }

    /**
     * Calculate how many nonzero elements will need to be stored.
     * For general matrices this will be the same as header.nnz, but if the MatrixMarket file has symmetry and
     * generalize symmetry is selected then this function will calculate the total needed.
     */
    inline int64_t get_storage_nnz(const matrix_market_header& header, const read_options options) {
        if (header.object == vector) {
            return header.nnz;
        }

        if (header.format == coordinate) {
            if (header.symmetry != general && options.generalize_symmetry) {
                return 2 * header.nnz;
            } else {
                return header.nnz;
            }
        } else {
            auto diag_count = header.nrows;
            auto off_diag_count = header.nrows * header.ncols - diag_count;
            auto off_diag_half = off_diag_count / 2;

            if (options.generalize_symmetry) {
                if (header.symmetry == skew_symmetric) {
                    // skew-symmetric diagonals must be zero
                    return off_diag_count;
                } else {
                    return header.nnz;
                }
            } else {
                switch (header.symmetry) {
                    case symmetric:
                        return off_diag_half + diag_count;
                    case skew_symmetric:
                        return off_diag_half;
                    case hermitian:
                        return off_diag_half + diag_count;
                    case general:
                        return header.nnz;
                }
            }
        }
        throw fmm_error("Unknown configuration for get_storage_nnz().");
    }

    /**
     * Parse a Matrix Market header comment line.
     * @param header
     * @param line
     * @return
     */
    inline bool read_comment(matrix_market_header& header, const std::string& line) {
        // empty lines are allowed anywhere in the file and are to be ignored
        if (is_line_all_spaces(line)) {
            return true;
        }

        unsigned int pos = 0;

        // skip leading whitespace
        while ((pos+1) < line.size() && std::isblank(line[pos])) {
            ++pos;
        }

        if (line[pos] != '%') {
            return false;
        }

        // skip the '%'
        ++pos;

        // Line is a comment. Save it to the header.
        header.comment += line.substr(pos) + "\n";
        return true;
    }

    /**
     * Parse an enum, but with error message fitting parsing of header.
     */
    template <typename ENUM>
    ENUM parse_header_enum(const std::string& s, std::map<ENUM, const std::string> mp, int64_t line_num) {
        // Make s lowercase for a case-insensitive match
        std::string lower(s);
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        for (const auto& [key, value] : mp) {
            if (value == lower) {
                return key;
            }
        }
        throw invalid_mm(std::string("Invalid MatrixMarket header element: ") + s, line_num);
    }


    /**
     * Reads
     * @param instream stream to read from
     * @param header structure that will be filled with read header
     * @return number of lines read
     */
    inline int64_t read_header(std::istream& instream, matrix_market_header& header) {
        int64_t lines_read = 0;
        std::string line;

        // read banner
        std::getline(instream, line);
        strip_trailing_cr(line);
        lines_read++;

        if (line.find("MatrixMarket", 0) == std::string::npos) {
            // not a matrix market file because the banner is missing
            throw invalid_mm("Not a Matrix Market file. Missing banner.", lines_read);
        }

        // parse banner
        {
            std::istringstream iss(line);
            std::string banner, f_object, f_format, f_field, f_symmetry;
            iss >> banner >> f_object >> f_format >> f_field >> f_symmetry;
            if (banner != kMatrixMarketBanner && banner != kMatrixMarketBanner2) {
                // not a matrix market file because the banner is wrong
                throw invalid_mm("Not a Matrix Market file. Missing banner.", lines_read);
            }

            header.object = parse_header_enum(f_object, object_map, lines_read);
            header.format = parse_header_enum(f_format, format_map, lines_read);
            header.field = parse_header_enum(f_field, field_map, lines_read);
            header.symmetry = parse_header_enum(f_symmetry, symmetry_map, lines_read);
        }

        // Read any comments
        do {
            std::getline(instream, line);
            strip_trailing_cr(line);
            lines_read++;

            if (!instream) {
                throw invalid_mm("Invalid MatrixMarket header: Premature EOF", lines_read);
            }
        } while (read_comment(header, line));

        // trim off final comment newline
        if (ends_with(header.comment, "\n")) {
            header.comment.resize(header.comment.size() - 1);
        }

        // parse the dimension line
        {
            std::istringstream iss(line);
            int expected_length = -1;

            const char* pos = line.c_str();
            const char* end = line.c_str() + line.size();
            pos = skip_spaces(pos);

            if (header.object == vector) {
                pos = read_int(pos, end, header.vector_length);

                if (header.vector_length < 0) {
                    throw invalid_mm("Vector length can't be negative.", lines_read);
                }

                if (header.format == coordinate) {
                    pos = skip_spaces(pos);
                    pos = read_int(pos, end, header.nnz);
                    expected_length = 2;
                } else {
                    header.nnz = header.vector_length;
                    expected_length = 1;
                }

                header.nrows = header.vector_length;
                header.ncols = 1;
            } else {
                pos = read_int(pos, end, header.nrows);
                pos = skip_spaces(pos);
                pos = read_int(pos, end, header.ncols);
                if (header.nrows < 0 || header.ncols < 0) {
                    throw invalid_mm("Matrix dimensions can't be negative.", lines_read);
                }

                if (header.format == coordinate) {
                    pos = skip_spaces(pos);
                    pos = read_int(pos, end, header.nnz);
                    if (header.nnz < 0) {
                        throw invalid_mm("Matrix NNZ can't be negative.", lines_read);
                    }
                    expected_length = 3;
                } else {
                    header.nnz = header.nrows * header.ncols;
                    expected_length = 2;
                }
                if (std::min(header.nrows, header.ncols) == 1) {
                    // row or column matrix. Either one can be loaded into a vector data structure.
                    header.vector_length = std::max(header.nrows, header.ncols);
                } else {
                    header.vector_length = -1;
                }
            }

            pos = skip_spaces(pos);
            if (pos != end) {
                // The text of this message is to be able to match SciPy's message to pass their unit test.
                throw invalid_mm("Header dimension line not of length " + std::to_string(expected_length));
            }
        }

        header.header_line_count = lines_read;

        return lines_read;
    }

    inline bool write_header(std::ostream& os, const matrix_market_header& header, const write_options options = {}) {
        // Write the banner
        os << kMatrixMarketBanner << kSpace;
        os << object_map.at(header.object) << kSpace;
        os << format_map.at(header.format) << kSpace;
        os << field_map.at(header.field) << kSpace;
        os << symmetry_map.at(header.symmetry) << kNewline;

        // Write the comment
        if (!header.comment.empty()) {
            std::string write_comment = replace_all(header.comment, "\n", "\n%");

            os << "%" << write_comment << kNewline;
        } else if (options.always_comment) {
            os << "%" << kNewline;
        }

        // Write dimension line
        if (header.object == vector) {
            os << header.vector_length;
            if (header.format == coordinate) {
                os << kSpace << header.nnz;
            }
        } else {
            os << header.nrows << kSpace << header.ncols;
            if (header.format == coordinate) {
                os << kSpace << header.nnz;
            }
        }
        os << kNewline;

        return true;
    }

}
