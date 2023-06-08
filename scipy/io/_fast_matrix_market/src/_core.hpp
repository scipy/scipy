// Copyright (C) 2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

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

void init_read_array(py::module_ &);
void init_write_array(py::module_ &);
void init_read_coo(py::module_ &);
void init_write_coo(py::module_ &);
void init_write_csc(py::module_ &);
