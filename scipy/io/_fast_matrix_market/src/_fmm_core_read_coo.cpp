// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_fmm_core.hpp"

/**
 * Read Matrix Market body into triplets.
 */
template <typename IT, typename VT>
void read_body_coo(read_cursor& cursor, py::array_t<IT>& row, py::array_t<IT>& col, py::array_t<VT>& data) {
    if (row.size() != cursor.header.nnz || col.size() != cursor.header.nnz || data.size() != cursor.header.nnz) {
        throw std::invalid_argument("NumPy Array sizes need to equal matrix nnz");
    }
    auto row_unchecked = row.mutable_unchecked();
    auto col_unchecked = col.mutable_unchecked();
    auto data_unchecked = data.mutable_unchecked();
    auto handler = fmm::triplet_calling_parse_handler<IT, VT, decltype(row_unchecked), decltype(data_unchecked)>(
        row_unchecked, col_unchecked, data_unchecked);

    // The mmread() will only call this method if the matrix is a coordinate. Disable the code paths for reading
    // array matrices here to reduce final library size and compilation time.
#ifdef FMM_SCIPY_PRUNE
    fmm::read_matrix_market_body<decltype(handler), fmm::compile_coordinate_only>(cursor.stream(), cursor.header, handler, 1, cursor.options);
#else
    fmm::read_matrix_market_body<decltype(handler), fmm::compile_all>(cursor.stream(), cursor.header, handler, 1, cursor.options);
#endif
    cursor.close();
}


void init_read_coo(py::module_ &m) {
    m.def("read_body_coo", &read_body_coo<int32_t, int64_t>);
    m.def("read_body_coo", &read_body_coo<int32_t, uint64_t>);
    m.def("read_body_coo", &read_body_coo<int32_t, double>);
    m.def("read_body_coo", &read_body_coo<int32_t, std::complex<double>>);

    m.def("read_body_coo", &read_body_coo<int64_t, int64_t>);
    m.def("read_body_coo", &read_body_coo<int64_t, uint64_t>);
    m.def("read_body_coo", &read_body_coo<int64_t, double>);
    m.def("read_body_coo", &read_body_coo<int64_t, std::complex<double>>);

#ifndef FMM_SCIPY_PRUNE
    m.def("read_body_coo", &read_body_coo<int32_t, float>);
    m.def("read_body_coo", &read_body_coo<int32_t, long double>);
    m.def("read_body_coo", &read_body_coo<int32_t, std::complex<float>>);
    m.def("read_body_coo", &read_body_coo<int32_t, std::complex<long double>>);

    m.def("read_body_coo", &read_body_coo<int64_t, float>);
    m.def("read_body_coo", &read_body_coo<int64_t, long double>);
    m.def("read_body_coo", &read_body_coo<int64_t, std::complex<float>>);
    m.def("read_body_coo", &read_body_coo<int64_t, std::complex<long double>>);
#endif
}