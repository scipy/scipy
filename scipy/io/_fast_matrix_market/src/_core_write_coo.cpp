// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_core.hpp"

/**
 * Write Python triplets to MatrixMarket.
 */
template <typename IT, typename VT>
void write_coo(write_cursor& cursor, const std::tuple<int64_t, int64_t>& shape,
               py::array_t<IT>& rows, py::array_t<IT>& cols, py::array_t<VT>& data) {
    if (rows.size() != cols.size()) {
        throw std::invalid_argument("len(row) must equal len(col).");
    }
    if (rows.size() != data.size() && data.size() != 0) {
        throw std::invalid_argument("len(row) must equal len(data).");
    }

    cursor.header.nrows = std::get<0>(shape);
    cursor.header.ncols = std::get<1>(shape);
    cursor.header.nnz = rows.size();

    cursor.header.object = fmm::matrix;
    cursor.header.field = (data.size() == 0 ? (cursor.header.nnz == 0 ? fmm::real : fmm::pattern) : fmm::get_field_type((const VT*)nullptr));
    cursor.header.format = fmm::coordinate;

    fmm::write_header(cursor.stream(), cursor.header);

    auto rows_unchecked = rows.unchecked();
    auto cols_unchecked = cols.unchecked();
    auto data_unchecked = data.unchecked();

    fmm::line_formatter<IT, VT> lf(cursor.header, cursor.options);
    auto formatter = fmm::triplet_formatter(lf,
                                            py_array_iterator<decltype(rows_unchecked), IT>(rows_unchecked),
                                            py_array_iterator<decltype(rows_unchecked), IT>(rows_unchecked, rows_unchecked.size()),
                                            py_array_iterator<decltype(cols_unchecked), IT>(cols_unchecked),
                                            py_array_iterator<decltype(cols_unchecked), IT>(cols_unchecked, cols_unchecked.size()),
                                            py_array_iterator<decltype(data_unchecked), VT>(data_unchecked),
                                            py_array_iterator<decltype(data_unchecked), VT>(data_unchecked, data_unchecked.size()));
    fmm::write_body(cursor.stream(), formatter, cursor.options);
}


void init_write_coo(py::module_ &m) {
    m.def("write_coo", &write_coo<int32_t, int32_t>);
    m.def("write_coo", &write_coo<int32_t, uint32_t>);
    m.def("write_coo", &write_coo<int32_t, int64_t>);
    m.def("write_coo", &write_coo<int32_t, uint64_t>);
    m.def("write_coo", &write_coo<int32_t, float>);
    m.def("write_coo", &write_coo<int32_t, double>);
    m.def("write_coo", &write_coo<int32_t, long double>);
    m.def("write_coo", &write_coo<int32_t, std::complex<float>>);
    m.def("write_coo", &write_coo<int32_t, std::complex<double>>);
    m.def("write_coo", &write_coo<int32_t, std::complex<long double>>);

    m.def("write_coo", &write_coo<int64_t, int32_t>);
    m.def("write_coo", &write_coo<int64_t, uint32_t>);
    m.def("write_coo", &write_coo<int64_t, int64_t>);
    m.def("write_coo", &write_coo<int64_t, uint64_t>);
    m.def("write_coo", &write_coo<int64_t, float>);
    m.def("write_coo", &write_coo<int64_t, double>);
    m.def("write_coo", &write_coo<int64_t, long double>);
    m.def("write_coo", &write_coo<int64_t, std::complex<float>>);
    m.def("write_coo", &write_coo<int64_t, std::complex<double>>);
    m.def("write_coo", &write_coo<int64_t, std::complex<long double>>);
}