// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <istream>
#include <string>

namespace fast_matrix_market {
    inline void get_next_chunk(std::string& chunk, std::istream &instream, const read_options &options) {
        constexpr size_t chunk_extra = 4096; // extra chunk bytes to leave room for rest of line
        size_t chunk_length = 0;

        // ensure enough space
        chunk.resize(options.chunk_size_bytes);

        // read chunk from the stream
        auto bytes_to_read = chunk.size() > chunk_extra ? (std::streamsize) (chunk.size() - chunk_extra) : 0;
        if (bytes_to_read > 0) {
            instream.read(chunk.data(), bytes_to_read);
            auto num_read = instream.gcount();
            chunk_length = num_read;

            // test for EOF
            if (num_read == 0 || instream.eof() || chunk[chunk_length - 1] == '\n') {
                chunk.resize(chunk_length);
                return;
            }
        }

        // Read rest of line and append to the chunk.
        std::string suffix;
        std::getline(instream, suffix);
        if (instream.good()) {
            suffix += "\n";
        }

        if (chunk_length + suffix.size() > chunk.size()) {
            // rest of line didn't fit in the extra space, must copy
            chunk.resize(chunk_length);
            chunk += suffix;
        } else {
            // the suffix fits in the chunk.
            std::copy(suffix.begin(), suffix.end(), chunk.begin() + (ptrdiff_t) chunk_length);
            chunk_length += suffix.size();
            chunk.resize(chunk_length);
        }
    }

    inline std::string get_next_chunk(std::istream &instream, const read_options &options) {
        // allocate chunk
        std::string chunk(options.chunk_size_bytes, ' ');
        get_next_chunk(chunk, instream, options);
        return chunk;
    }

    template <typename ITER>
    bool is_all_spaces(ITER begin, ITER end) {
        return std::all_of(begin, end, [](char c) { return c == ' ' || c == '\t' || c == '\r'; });
    }

    /**
     * Find the number of total lines and empty lines in a multiline string.
     */
    inline std::pair<int64_t, int64_t> count_lines(const std::string& chunk) {
        int64_t num_newlines = 0;
        int64_t num_empty_lines = 0;

        auto pos = std::cbegin(chunk);
        auto end = std::cend(chunk);
        auto line_start = pos;
        for (; pos != end; ++pos) {
            if (*pos == '\n') {
                ++num_newlines;
                if (is_all_spaces(line_start, pos)) {
                    ++num_empty_lines;
                }
                line_start = pos + 1;
            }
        }

        if (line_start != end) {
            // last line does not end in newline, but it might still be empty
            if (is_all_spaces(line_start, end)) {
                ++num_empty_lines;
            }
        }

        if (num_newlines == 0) {
            // single line is still a line
            if (chunk.empty()) {
                num_empty_lines = 1;
            }
            return std::make_pair(1, num_empty_lines);
        }

        if (chunk[chunk.size()-1] != '\n') {
            ++num_newlines;
        }

        return std::make_pair(num_newlines, num_empty_lines);
    }
}
