// Copyright (C) 2022 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <queue>

#include "fast_matrix_market.hpp"
#include "thirdparty/task_thread_pool.hpp"

namespace fast_matrix_market {
    /**
     * Write Matrix Market body.
     *
     * Chunk based so that it can be made parallel. Each chunk is written by a FORMATTER class.
     * @tparam FORMATTER implementation class that writes chunks.
     */
    template <typename FORMATTER>
    void write_body_threads(std::ostream& os,
                            FORMATTER& formatter, const write_options& options = {}) {
        /*
         * Requirements:
         * Chunks must be created sequentially by the formatter.
         * Chunks can be computed in parallel (i.e. call operator()).
         * Chunk results must be written in the same order that they were created in.
         *
         * This is effectively a pipeline with a serial producer (chunk generator), parallel workers, and serial
         * consumer (writer).
         *
         * The biggest obstacle is the final requirement to write all chunks sequentially.
         *
         * We take a simple approach. The main thread handles the serial chunk generation and I/O,
         * and a thread pool performs the parallel work.
         */
        std::queue<std::future<std::string>> futures;
        task_thread_pool::task_thread_pool pool(options.num_threads);

        // Number of concurrent chunks available to work on.
        // Too few may starve workers (such as due to uneven chunk splits)
        // Too many increases costs, such as storing chunk results in memory before they're written.
        const int inflight_count = 2 * (int)pool.get_num_threads();

        // Start computing tasks.
        for (int batch_i = 0; batch_i < inflight_count && formatter.has_next(); ++batch_i) {
            // Could push the chunk directly, but MSVC.
            futures.push(pool.submit([](auto chunk){ return chunk(); }, formatter.next_chunk(options)));
//            futures.push(pool.submit(formatter.next_chunk(options)));
        }

        // Write chunks in order as they become available.
        while (!futures.empty()) {
            std::string chunk = futures.front().get();
            futures.pop();

            // Next chunk is ready. Start another to replace it.
            if (formatter.has_next()) {
                futures.push(pool.submit([](auto chunk){ return chunk(); }, formatter.next_chunk(options)));
            }

            // Write this one out.
            os.write(chunk.c_str(), (std::streamsize) chunk.size());
        }
    }
}
