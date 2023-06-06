// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_core.hpp"

/**
 * Write Python CSC/CSR to MatrixMarket.
 */
template <typename IT, typename VT>
void write_csc(write_cursor& cursor, const std::tuple<int64_t, int64_t>& shape,
               py::array_t<IT>& indptr, py::array_t<IT>& indices, py::array_t<VT>& data, bool is_csr) {
    if (indices.size() != data.size() && data.size() != 0) {
        throw std::invalid_argument("len(indices) must equal len(data).");
    }

    cursor.header.nrows = std::get<0>(shape);
    cursor.header.ncols = std::get<1>(shape);
    cursor.header.nnz = indices.size();

    if ((is_csr && indptr.size() != cursor.header.nrows + 1) ||
        (!is_csr && indptr.size() != cursor.header.ncols + 1)) {
        throw std::invalid_argument("indptr length does not match matrix shape.");
    }

    cursor.header.object = fmm::matrix;
    cursor.header.field = (data.size() == 0 ? (cursor.header.nnz == 0 ? fmm::real : fmm::pattern) : fmm::get_field_type((const VT*)nullptr));
    cursor.header.format = fmm::coordinate;
    cursor.header.symmetry = fmm::general;

    fmm::write_header(cursor.stream(), cursor.header);

    auto indptr_unchecked = indptr.unchecked();
    auto indices_unchecked = indices.unchecked();
    auto data_unchecked = data.unchecked();

    fmm::line_formatter<IT, VT> lf(cursor.header, cursor.options);
    auto formatter = fmm::csc_formatter(lf,
                                        py_array_iterator<decltype(indptr_unchecked), IT>(indptr_unchecked),
                                        py_array_iterator<decltype(indptr_unchecked), IT>(indptr_unchecked, indptr_unchecked.size() - 1),
                                        py_array_iterator<decltype(indices_unchecked), IT>(indices_unchecked),
                                        py_array_iterator<decltype(indices_unchecked), IT>(indices_unchecked, indices_unchecked.size()),
                                        py_array_iterator<decltype(data_unchecked), VT>(data_unchecked),
                                        py_array_iterator<decltype(data_unchecked), VT>(data_unchecked, data_unchecked.size()),
                                        is_csr);
    fmm::write_body(cursor.stream(), formatter, cursor.options);
}


void init_write_csc(py::module_ &m) {
    m.def("write_csc", &write_csc<int32_t, int32_t>);
    m.def("write_csc", &write_csc<int32_t, uint32_t>);
    m.def("write_csc", &write_csc<int32_t, int64_t>);
    m.def("write_csc", &write_csc<int32_t, uint64_t>);
    m.def("write_csc", &write_csc<int32_t, float>);
    m.def("write_csc", &write_csc<int32_t, double>);
    m.def("write_csc", &write_csc<int32_t, long double>);
    m.def("write_csc", &write_csc<int32_t, std::complex<float>>);
    m.def("write_csc", &write_csc<int32_t, std::complex<double>>);
    m.def("write_csc", &write_csc<int32_t, std::complex<long double>>);

    m.def("write_csc", &write_csc<int64_t, int32_t>);
    m.def("write_csc", &write_csc<int64_t, uint32_t>);
    m.def("write_csc", &write_csc<int64_t, int64_t>);
    m.def("write_csc", &write_csc<int64_t, uint64_t>);
    m.def("write_csc", &write_csc<int64_t, float>);
    m.def("write_csc", &write_csc<int64_t, double>);
    m.def("write_csc", &write_csc<int64_t, long double>);
    m.def("write_csc", &write_csc<int64_t, std::complex<float>>);
    m.def("write_csc", &write_csc<int64_t, std::complex<double>>);
    m.def("write_csc", &write_csc<int64_t, std::complex<long double>>);
}