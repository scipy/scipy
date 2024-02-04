// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <future>
#include <queue>

#include "fast_matrix_market.hpp"
#include "thirdparty/task_thread_pool.hpp"

namespace fast_matrix_market {

    struct line_count_result_s {
        std::string chunk;
        line_counts counts;

        explicit line_count_result_s(std::string && c): chunk(c) {}
    };

    using line_count_result = std::shared_ptr<line_count_result_s>;

    inline line_count_result count_chunk_lines(line_count_result lcr) {
        auto [lines, empties] = count_lines(lcr->chunk);

        lcr->counts.file_line = lines;
        lcr->counts.element_num = lines - empties;
        return lcr;
    }

    template <typename HANDLER, compile_format FORMAT = compile_all>
    line_counts read_body_threads(std::istream& instream, const matrix_market_header& header,
                                  HANDLER& handler, const read_options& options = {}) {
        /*
         * Pipeline:
         * 1. Read chunk
         * 2. Calculate chunk's line count
         * 3. Parse chunk.
         *
         * The line count is needed for
         * 1. for array files the line number determines the row/column indices of the value
         * 2. for coordinate files the line number determines the chunk's offset into the result arrays
         * 3. for error messages
         *
         * Only the main threads performs I/O reads. Everything else is done by tasks in a thread pool.
         *
         * The line count is fast, but we still spawn line count tasks. The futures for these tasks are saved in a
         * queue to be retrieved in order. This enables easy tracking of the line numbers of each chunk.
         *
         * Once a chunk's line count is complete we spawn a task to parse it. We also then read another chunk from
         * the input stream and start its line count.
         *
         * The line count step is significantly faster than the parse step. As a form of backpressure we don't read
         * additional chunks if there are too many inflight chunks.
         */
        line_counts lc{header.header_line_count, 0};

        std::queue<std::future<line_count_result>> line_count_futures;
        std::queue<std::future<line_count_result>> parse_futures;
        task_thread_pool::task_thread_pool pool(options.num_threads);

        // Reuse the line_count_result objects. Each chunk would otherwise allocate a new 1MB std::string.
        // The lifetime of these strings is relatively short, but some allocators do not immediately reuse the memory.
        // This object pool can reduce overall RSS memory usage in many cases.
        std::queue<line_count_result> lcr_reuse_pool;

        int generalizing_symmetry_factor = (header.symmetry != general && options.generalize_symmetry) ? 2 : 1;

        // Number of concurrent chunks available to work on.
        // Too few may starve workers (such as due to uneven chunk splits)
        // Too many increases costs, such as storing chunk results in memory before they're written.
        const unsigned inflight_count = pool.get_num_threads() + 1;

        // Start reading chunks and counting lines.
        for (unsigned seed_i = 0; seed_i < inflight_count && instream.good(); ++seed_i) {
            line_count_result lcr = std::make_shared<line_count_result_s>(get_next_chunk(instream, options));
            line_count_futures.push(pool.submit(count_chunk_lines, lcr));
        }

        // Read chunks in order, as they become available.
        while (!line_count_futures.empty()) {

            // Wait on any parse results. This serves as backpressure.
            while (!parse_futures.empty() && (is_ready(parse_futures.front()) || parse_futures.size() > inflight_count)) {
                // This will throw any parse errors.
                auto lcr_to_reuse = parse_futures.front().get();
                parse_futures.pop();

                // save the lcr struct to reuse the memory
                lcr_reuse_pool.push(lcr_to_reuse);
            }

            // We are ready to start another parse task.
            line_count_result lcr = line_count_futures.front().get();
            line_count_futures.pop();

            // Next chunk has finished line count. Start another to replace it.
            if (instream.good()) {
                line_count_result lcr_reuse;
                // attempt to reuse the chunk string object from a previous chunk
                if (lcr_reuse_pool.empty()) {
                    lcr_reuse = std::make_shared<line_count_result_s>("");
                } else {
                    lcr_reuse = lcr_reuse_pool.front();
                    lcr_reuse_pool.pop();
                }

                get_next_chunk(lcr_reuse->chunk, instream, options);
                line_count_futures.push(pool.submit(count_chunk_lines, lcr_reuse));
            }

            // Parse it.
            if (lc.element_num > header.nnz) {
                throw invalid_mm("File too long", lc.file_line + 1);
            }
            auto chunk_handler = handler.get_chunk_handler(lc.element_num * generalizing_symmetry_factor);
            if (header.format == array) {
                if constexpr ((FORMAT & compile_array_only) == compile_array_only) {
                    // compute the starting row/column for this array chunk
                    typename HANDLER::coordinate_type row = lc.element_num % header.nrows;
                    typename HANDLER::coordinate_type col = lc.element_num / header.nrows;

                    parse_futures.push(pool.submit([=]() mutable {
                        read_chunk_array(lcr->chunk, header, lc, chunk_handler, options, row, col);
                        return lcr;
                    }));
                } else {
                    throw support_not_selected("Matrix is array but reading array files not enabled for this method.");
                }
            } else if (header.object == matrix) {
                if constexpr ((FORMAT & compile_coordinate_only) == compile_coordinate_only) {
                    parse_futures.push(pool.submit([=]() mutable {
                        read_chunk_matrix_coordinate(lcr->chunk, header, lc, chunk_handler, options);
                        return lcr;
                    }));
                } else {
                    throw support_not_selected("Matrix is coordinate but reading coordinate files not enabled for this method.");
                }
            } else {
#ifdef FMM_NO_VECTOR
                throw no_vector_support("Vector Matrix Market files not supported.");
#else
                parse_futures.push(pool.submit([=]() mutable {
                    read_chunk_vector_coordinate(lcr->chunk, header, lc, chunk_handler, options);
                    return lcr;
                }));
#endif
            }

            // Advance counts for next chunk
            lc.file_line += lcr->counts.file_line;
            lc.element_num += lcr->counts.element_num;
        }

        // Wait on any parse results. This will throw any parse errors.
        while (!parse_futures.empty()) {
            parse_futures.front().get();
            parse_futures.pop();
        }

        return lc;
    }
}
