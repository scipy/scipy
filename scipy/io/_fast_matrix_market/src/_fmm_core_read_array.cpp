// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_fmm_core.hpp"

/**
 * Read Matrix Market body into a numpy array.
 *
 * @param cursor Opened by open_read().
 * @param array NumPy array. Assumed to be the correct size and zeroed out.
 */
template <typename T>
void read_body_array(read_cursor& cursor, py::array_t<T>& array) {
    cursor.options.generalize_symmetry = true;
    auto unchecked = array.mutable_unchecked();
    auto handler = fmm::dense_2d_call_adding_parse_handler<decltype(unchecked), int64_t, T>(unchecked);

    // The mmread() will only call this method if the matrix is an array. Disable the code paths for reading
    // coordinate matrices here to reduce final library size and compilation time.
#ifdef FMM_SCIPY_PRUNE
    fmm::read_matrix_market_body<decltype(handler), fmm::compile_array_only>(cursor.stream(), cursor.header, handler, 1, cursor.options);
#else
    fmm::read_matrix_market_body<decltype(handler), fmm::compile_all>(cursor.stream(), cursor.header, handler, 1, cursor.options);
#endif
    cursor.close();
}


void init_read_array(py::module_ &m) {
    m.def("read_body_array", &read_body_array<int64_t>);
    m.def("read_body_array", &read_body_array<uint64_t>);
    m.def("read_body_array", &read_body_array<double>);
    m.def("read_body_array", &read_body_array<std::complex<double>>);

#ifndef FMM_SCIPY_PRUNE
    m.def("read_body_array", &read_body_array<float>);
    m.def("read_body_array", &read_body_array<long double>);
    m.def("read_body_array", &read_body_array<std::complex<float>>);
    m.def("read_body_array", &read_body_array<std::complex<long double>>);
#endif
}