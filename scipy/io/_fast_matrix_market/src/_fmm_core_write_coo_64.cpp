// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_fmm_core.hpp"

void init_write_coo_64(py::module_ &m) {
    m.def("write_body_coo", &write_body_coo<int64_t, int32_t>);
    m.def("write_body_coo", &write_body_coo<int64_t, uint32_t>);
    m.def("write_body_coo", &write_body_coo<int64_t, int64_t>);
    m.def("write_body_coo", &write_body_coo<int64_t, uint64_t>);
    m.def("write_body_coo", &write_body_coo<int64_t, float>);
    m.def("write_body_coo", &write_body_coo<int64_t, double>);
    m.def("write_body_coo", &write_body_coo<int64_t, long double>);
    m.def("write_body_coo", &write_body_coo<int64_t, std::complex<float>>);
    m.def("write_body_coo", &write_body_coo<int64_t, std::complex<double>>);
    m.def("write_body_coo", &write_body_coo<int64_t, std::complex<long double>>);
}