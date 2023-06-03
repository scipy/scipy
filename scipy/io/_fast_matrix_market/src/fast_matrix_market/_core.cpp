// Copyright (C) 2022-2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.

#include <fstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "pystreambuf.h"

#include <fast_matrix_market/types.hpp>
namespace fast_matrix_market {
    // Be able to set unsigned-integer field type. This type is only used by SciPy to represent uint64 values.
    field_type get_field_type([[maybe_unused]] const uint32_t* type) {
        return unsigned_integer;
    }

    field_type get_field_type([[maybe_unused]] const uint64_t* type) {
        return unsigned_integer;
    }
}
#include <fast_matrix_market/fast_matrix_market.hpp>

namespace py = pybind11;
using namespace pybind11::literals;
namespace fmm = fast_matrix_market;


////////////////////////////////////////////////
//// Header methods
////////////////////////////////////////////////

std::tuple<int64_t, int64_t> get_header_shape(const fmm::matrix_market_header& header) {
    return std::make_tuple(header.nrows, header.ncols);
}

void set_header_shape(fmm::matrix_market_header& header, const std::tuple<int64_t, int64_t>& shape) {
    header.nrows = std::get<0>(shape);
    header.ncols = std::get<1>(shape);
}

std::string get_header_object(const fmm::matrix_market_header& header) {
    return fmm::object_map.at(header.object);
}
std::string get_header_format(const fmm::matrix_market_header& header) {
    return fmm::format_map.at(header.format);
}
std::string get_header_field(const fmm::matrix_market_header& header) {
    return fmm::field_map.at(header.field);
}
std::string get_header_symmetry(const fmm::matrix_market_header& header) {
    return fmm::symmetry_map.at(header.symmetry);
}

void set_header_object(fmm::matrix_market_header& header, const std::string& value) {
    header.object = fmm::parse_enum<fmm::object_type>(value, fmm::object_map);
}
void set_header_format(fmm::matrix_market_header& header, const std::string& value) {
    header.format = fmm::parse_enum<fmm::format_type>(value, fmm::format_map);
}
void set_header_field(fmm::matrix_market_header& header, const std::string& value) {
    header.field = fmm::parse_enum<fmm::field_type>(value, fmm::field_map);
}
void set_header_symmetry(fmm::matrix_market_header& header, const std::string& value) {
    header.symmetry = fmm::parse_enum<fmm::symmetry_type>(value, fmm::symmetry_map);
}

fmm::matrix_market_header create_header(const std::tuple<int64_t, int64_t>& shape, int64_t nnz,
                                        const std::string& comment,
                                        const std::string& object, const std::string& format,
                                        const std::string& field, const std::string& symmetry) {
    fmm::matrix_market_header header{};
    set_header_shape(header, shape);
    header.nnz = nnz;
    header.comment = comment;
    set_header_object(header, object);
    set_header_format(header, format);
    set_header_field(header, field);
    set_header_symmetry(header, symmetry);
    return header;
}

py::dict header_to_dict(fmm::matrix_market_header& header) {
    py::dict dict;
    dict["shape"] = py::make_tuple(header.nrows, header.ncols);
    dict["nnz"] = header.nnz;
    dict["comment"] = header.comment;
    dict["object"] = get_header_object(header);
    dict["format"] = get_header_format(header);
    dict["field"] = get_header_field(header);
    dict["symmetry"] = get_header_symmetry(header);

    return dict;
}

std::string header_repr(const fmm::matrix_market_header& header) {
    std::ostringstream oss;
    oss << "header(";
    oss << "shape=(" << header.nrows << ", " << header.ncols << "), ";
    oss << "nnz=" << header.nnz << ", ";
    oss << "comment=\"" << header.comment << "\", ";
    oss << "object=\"" << get_header_object(header) << "\", ";
    oss << "format=\"" << get_header_format(header) << "\", ";
    oss << "field=\"" << get_header_field(header) << "\", ";
    oss << "symmetry=\"" << get_header_symmetry(header) << "\"";
    oss << ")";
    return oss.str();
}

////////////////////////////////////////////////
//// Read cursor - open files/streams for reading
////////////////////////////////////////////////

/**
 * A structure that represents an open MatrixMarket file or stream
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

void open_read_rest(read_cursor& cursor) {
    // This is done later in Python to match SciPy behavior
    cursor.options.generalize_symmetry = false;

    // read header
    fmm::read_header(cursor.stream(), cursor.header);
}

read_cursor open_read_file(const std::string& filename, int num_threads) {
    read_cursor cursor(filename);
    // Set options
    cursor.options.num_threads = num_threads;
    // Python parses 1e9999 as Inf
    cursor.options.float_out_of_range_behavior = fmm::BestMatch;

    open_read_rest(cursor);
    return cursor;
}

read_cursor open_read_stream(std::shared_ptr<pystream::istream>& external, int num_threads) {
    read_cursor cursor(external);
    // Set options
    cursor.options.num_threads = num_threads;
    // Python parses 1e9999 as Inf
    cursor.options.float_out_of_range_behavior = fmm::BestMatch;

    open_read_rest(cursor);
    return cursor;
}

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
    fmm::read_matrix_market_body(cursor.stream(), cursor.header, handler, 1, cursor.options);
}

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
    fmm::read_matrix_market_body(cursor.stream(), cursor.header, handler, 1, cursor.options);
}

////////////////////////////////////////////////
//// Write cursor - open files/streams writing reading
////////////////////////////////////////////////

struct write_cursor {
    /**
     * Open a file
     * @param filename path
     */
    write_cursor(const std::string& filename): stream_ptr(std::make_unique<std::ofstream>(filename)) {}

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

write_cursor open_write_file(const std::string& filename, const fmm::matrix_market_header& header,
                             int num_threads, int precision) {
    write_cursor cursor(filename);
    // Set options
    cursor.options.num_threads = num_threads;
    cursor.options.precision = precision;
    cursor.header = header;
    return cursor;
}

write_cursor open_write_stream(std::shared_ptr<pystream::ostream>& stream, fmm::matrix_market_header& header,
                               int num_threads, int precision) {
    write_cursor cursor(stream);
    // Set options
    cursor.options.num_threads = num_threads;
    cursor.options.precision = precision;
    cursor.header = header;
    return cursor;
}

void write_header_only(write_cursor& cursor) {
    fmm::write_header(cursor.stream(), cursor.header);
}

/**
 * Write numpy array to MatrixMarket file
 */
template <typename T>
void write_array(write_cursor& cursor, py::array_t<T>& array) {
    if (array.ndim() != 2) {
        throw std::invalid_argument("Only 2D arrays supported.");
    }

    cursor.header.nrows = array.shape(0);
    cursor.header.ncols = array.shape(1);

    cursor.header.object = fmm::matrix;
    cursor.header.field = fmm::get_field_type((const T*)nullptr);
    cursor.header.format = fmm::array;

    fmm::write_header(cursor.stream(), cursor.header);

    auto unchecked = array.unchecked();
    fmm::line_formatter<int64_t, T> lf(cursor.header, cursor.options);
    auto formatter = fmm::dense_2d_call_formatter<decltype(lf), decltype(unchecked), int64_t>(
            lf, unchecked, cursor.header.nrows, cursor.header.ncols);
    fmm::write_body(cursor.stream(), formatter, cursor.options);
}

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

////////////////////////////////////////////////
//// pybind11 module definition
//// Define the _core module here, it is used by __init__.py
////////////////////////////////////////////////


PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        fast_matrix_market
    )pbdoc";

    // translate exceptions
    py::register_local_exception_translator([](std::exception_ptr p) {
        try {
            if (p) {
                std::rethrow_exception(p);
            }
        } catch (const fmm::out_of_range &e) {
            PyErr_SetString(PyExc_OverflowError, e.what());
        } catch (const fmm::fmm_error &e) {
            // Everything else we throw maps best to ValueError
            PyErr_SetString(PyExc_ValueError, e.what());
        }
    });

    py::class_<fmm::matrix_market_header>(m, "header")
    .def(py::init<>())
    .def(py::init<int64_t, int64_t>())
    .def(py::init([](std::tuple<int64_t, int64_t> shape) { return fmm::matrix_market_header{std::get<0>(shape), std::get<1>(shape)}; }))
    .def(py::init(&create_header), py::arg("shape")=std::make_tuple(0, 0), "nnz"_a=0, "comment"_a=std::string(), "object"_a="matrix", "format"_a="coordinate", "field"_a="real", "symmetry"_a="general")
    .def_readwrite("nrows", &fmm::matrix_market_header::nrows)
    .def_readwrite("ncols", &fmm::matrix_market_header::ncols)
    .def_property("shape", &get_header_shape, &set_header_shape)
    .def_readwrite("nnz", &fmm::matrix_market_header::nnz)
    .def_readwrite("comment", &fmm::matrix_market_header::comment)
    .def_property("object", &get_header_object, &set_header_object)
    .def_property("format", &get_header_format, &set_header_format)
    .def_property("field", &get_header_field, &set_header_field)
    .def_property("symmetry", &get_header_symmetry, &set_header_symmetry)
    .def("to_dict", &header_to_dict, R"pbdoc(
        Return the values in the header as a dict.
    )pbdoc")
    .def("__repr__", [](const fmm::matrix_market_header& header) { return header_repr(header); });

    m.def("write_header_only", &write_header_only);

    // Read methods
    py::class_<read_cursor>(m, "_read_cursor")
    .def_readonly("header", &read_cursor::header);

    m.def("open_read_file", &open_read_file);
    m.def("open_read_stream", &open_read_stream);

    // Read arrays
    m.def("read_body_array", &read_body_array<int64_t>);
    m.def("read_body_array", &read_body_array<uint64_t>);
    m.def("read_body_array", &read_body_array<float>);
    m.def("read_body_array", &read_body_array<double>);
    m.def("read_body_array", &read_body_array<long double>);
    m.def("read_body_array", &read_body_array<std::complex<float>>);
    m.def("read_body_array", &read_body_array<std::complex<double>>);
    m.def("read_body_array", &read_body_array<std::complex<long double>>);

    // Read triplets
    m.def("read_body_coo", &read_body_coo<int32_t, int64_t>);
    m.def("read_body_coo", &read_body_coo<int32_t, uint64_t>);
    m.def("read_body_coo", &read_body_coo<int32_t, float>);
    m.def("read_body_coo", &read_body_coo<int32_t, double>);
    m.def("read_body_coo", &read_body_coo<int32_t, long double>);
    m.def("read_body_coo", &read_body_coo<int32_t, std::complex<double>>);
    m.def("read_body_coo", &read_body_coo<int32_t, std::complex<long double>>);

    m.def("read_body_coo", &read_body_coo<int64_t, int64_t>);
    m.def("read_body_coo", &read_body_coo<int64_t, uint64_t>);
    m.def("read_body_coo", &read_body_coo<int64_t, float>);
    m.def("read_body_coo", &read_body_coo<int64_t, double>);
    m.def("read_body_coo", &read_body_coo<int64_t, long double>);
    m.def("read_body_coo", &read_body_coo<int64_t, std::complex<float>>);
    m.def("read_body_coo", &read_body_coo<int64_t, std::complex<double>>);
    m.def("read_body_coo", &read_body_coo<int64_t, std::complex<long double>>);

    // Write methods
    py::class_<write_cursor>(m, "_write_cursor")
    .def_readwrite("header", &write_cursor::header);

    m.def("open_write_file", &open_write_file);
    m.def("open_write_stream", &open_write_stream);

    // Write arrays
    m.def("write_array", &write_array<int32_t>);
    m.def("write_array", &write_array<uint32_t>);
    m.def("write_array", &write_array<int64_t>);
    m.def("write_array", &write_array<uint64_t>);
    m.def("write_array", &write_array<float>);
    m.def("write_array", &write_array<double>);
    m.def("write_array", &write_array<long double>);
    m.def("write_array", &write_array<std::complex<float>>);
    m.def("write_array", &write_array<std::complex<double>>);
    m.def("write_array", &write_array<std::complex<long double>>);

    // Write triplets
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

    // Write CSC/CSR
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

    // Module version
#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
