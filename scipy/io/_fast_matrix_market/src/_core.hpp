// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#ifdef FMM_SCIPY_PRUNE
#define FMM_NO_VECTOR
#endif

#include <fstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "pystreambuf.h"

#include <fast_matrix_market/types.hpp>
namespace fast_matrix_market {
    // Be able to set unsigned-integer field type. This type is only used by SciPy to represent uint64 values.
    field_type get_field_type([[maybe_unused]] const uint32_t* type);

    field_type get_field_type([[maybe_unused]] const uint64_t* type);
}
#include <fast_matrix_market/fast_matrix_market.hpp>

namespace py = pybind11;
using namespace pybind11::literals;
namespace fmm = fast_matrix_market;

/**
 * A structure that represents an open MatrixMarket file or stream (for reading)
 */
struct read_cursor {
    /**
     * Open a file.
     */
    read_cursor(const std::string& filename): stream_ptr(std::make_shared<std::ifstream>(filename)) {}

    /**
     * Use a Python stream. Needs to be a shared_ptr because this stream object needs to stay alive for the lifetime
     * of this cursor object.
     */
    read_cursor(std::shared_ptr<pystream::istream>& external): stream_ptr(external) {}

    std::shared_ptr<std::istream> stream_ptr;

    fmm::matrix_market_header header{};
    fmm::read_options options{};

    std::istream& stream() {
        return *stream_ptr;
    }

    /**
     * Finish using the cursor. If a file has been opened it will be closed.
     */
    void close() {
        // If stream is a std::ifstream() then close the file.
        std::ifstream* f = dynamic_cast<std::ifstream*>(stream_ptr.get());
        if (f != nullptr) {
            f->close();
        }

        // Remove this reference to the stream.
        stream_ptr.reset();
    }
};

/**
 * A structure that represents an open MatrixMarket file or stream (for writing)
 */
struct write_cursor {
    /**
     * Open a file
     * @param filename path
     */
    write_cursor(const std::string& filename): stream_ptr(std::make_unique<std::ofstream>(filename, std::ios_base::out | std::ios_base::binary)) {}

    /**
     * Use a Python stream. Needs to be a shared_ptr because this stream object needs to stay alive for the lifetime
     * of this cursor object.
     */
    write_cursor(std::shared_ptr<pystream::ostream>& external): stream_ptr(external) {}

    std::shared_ptr<std::ostream> stream_ptr;

    fmm::matrix_market_header header{};
    fmm::write_options options{};

    std::ostream& stream() {
        return *stream_ptr;
    }

    /**
     * Finish using the cursor. Flush the backing stream, and if a file has been opened it will be closed.
     */
    void close() {
        // If stream is a std::ofstream() then close the file.
        std::ofstream* f = dynamic_cast<std::ofstream*>(stream_ptr.get());
        if (f != nullptr) {
            f->close();
        } else {
            stream_ptr->flush();
        }

        // Remove this reference to the stream.
        stream_ptr.reset();
    }
};

/**
 * An iterator adapter over py::array_t numpy arrays.
 *
 * This allows using the iterator-based fast_matrix_market methods.
 */
template<typename ARR, typename T>
class py_array_iterator
{
public:
    using value_type = T;
    using difference_type = int64_t;

    py_array_iterator(ARR& array) : array(array), index(0) {}
    py_array_iterator(ARR& array, int64_t index) : array(array), index(index) {}
    py_array_iterator(const py_array_iterator &rhs) : array(rhs.array), index(rhs.index) {}
    /* py_array_iterator& operator=(Type* rhs) {index = rhs; return *this;} */
    py_array_iterator& operator=(const py_array_iterator &rhs) {index = rhs.index; return *this;}
    py_array_iterator& operator+=(difference_type rhs) {index += rhs; return *this;}
    py_array_iterator& operator-=(difference_type rhs) {index -= rhs; return *this;}
    T operator*() const {return array(index);}
//    T* operator->() const {return index;}
    T operator[](difference_type rhs) const {return array(index + rhs);}

    py_array_iterator& operator++() {++index; return *this;}
    py_array_iterator& operator--() {--index; return *this;}
    py_array_iterator operator++(int) {py_array_iterator tmp(*this); ++index; return tmp;}
    py_array_iterator operator--(int) {py_array_iterator tmp(*this); --index; return tmp;}
    /* py_array_iterator operator+(const py_array_iterator& rhs) {return py_array_iterator(array, index+rhs.index);} */
    difference_type operator-(const py_array_iterator& rhs) const {return index-rhs.index;}
    py_array_iterator operator+(difference_type rhs) const {return py_array_iterator(array, index+rhs);}
    py_array_iterator operator-(difference_type rhs) const {return py_array_iterator(array, index-rhs);}
    friend py_array_iterator operator+(difference_type lhs, const py_array_iterator& rhs) {return py_array_iterator(rhs.array, lhs+rhs.index);}
    friend py_array_iterator operator-(difference_type lhs, const py_array_iterator& rhs) {return py_array_iterator(rhs.array, lhs-rhs.index);}

    bool operator==(const py_array_iterator& rhs) const {return index == rhs.index;}
    bool operator!=(const py_array_iterator& rhs) const {return index != rhs.index;}
    bool operator>(const py_array_iterator& rhs) const {return index > rhs.index;}
    bool operator<(const py_array_iterator& rhs) const {return index < rhs.index;}
    bool operator>=(const py_array_iterator& rhs) const {return index >= rhs.index;}
    bool operator<=(const py_array_iterator& rhs) const {return index <= rhs.index;}
private:
    ARR& array;
    int64_t index;
};


/**
 * Write Python triplets to MatrixMarket.
 */
template <typename IT, typename VT>
void write_body_coo(write_cursor& cursor, const std::tuple<int64_t, int64_t>& shape,
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

    fmm::write_header(cursor.stream(), cursor.header, cursor.options);

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
    cursor.close();
}

#ifndef FMM_SCIPY_PRUNE
/**
 * Write Python CSC/CSR to MatrixMarket.
 */
template <typename IT, typename VT>
void write_body_csc(write_cursor& cursor, const std::tuple<int64_t, int64_t>& shape,
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

    fmm::write_header(cursor.stream(), cursor.header, cursor.options);

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
    cursor.close();
}
#endif

void init_read_array(py::module_ &);
void init_write_array(py::module_ &);
void init_read_coo(py::module_ &);
void init_write_coo_32(py::module_ &);
void init_write_coo_64(py::module_ &);
void init_write_csc_32(py::module_ &);
void init_write_csc_64(py::module_ &);
