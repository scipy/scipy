// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_fmm_core.hpp"

/**
 * Write numpy array to MatrixMarket file
 */
template <typename T>
void write_body_array(write_cursor& cursor, py::array_t<T>& array) {
    if (array.ndim() != 2) {
        throw std::invalid_argument("Only 2D arrays supported.");
    }

    cursor.header.nrows = array.shape(0);
    cursor.header.ncols = array.shape(1);

    cursor.header.object = fmm::matrix;
    cursor.header.field = fmm::get_field_type((const T*)nullptr);
    cursor.header.format = fmm::array;

    fmm::write_header(cursor.stream(), cursor.header, cursor.options);

    auto unchecked = array.unchecked();
    fmm::line_formatter<int64_t, T> lf(cursor.header, cursor.options);
    auto formatter = fmm::dense_2d_call_formatter<decltype(lf), decltype(unchecked), int64_t>(
        lf, unchecked, cursor.header.nrows, cursor.header.ncols);
    fmm::write_body(cursor.stream(), formatter, cursor.options);
    cursor.close();
}


void init_write_array(py::module_ &m) {
    m.def("write_body_array", &write_body_array<int32_t>);
    m.def("write_body_array", &write_body_array<uint32_t>);
    m.def("write_body_array", &write_body_array<int64_t>);
    m.def("write_body_array", &write_body_array<uint64_t>);
    m.def("write_body_array", &write_body_array<float>);
    m.def("write_body_array", &write_body_array<double>);
    m.def("write_body_array", &write_body_array<long double>);
    m.def("write_body_array", &write_body_array<std::complex<float>>);
    m.def("write_body_array", &write_body_array<std::complex<double>>);
    m.def("write_body_array", &write_body_array<std::complex<long double>>);
}