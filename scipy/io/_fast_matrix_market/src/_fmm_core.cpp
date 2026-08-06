// Copyright (C) 2022-2023 Adam Lugowski. All rights reserved.
// Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
// SPDX-License-Identifier: BSD-2-Clause

#include "_fmm_core.hpp"

#include <fast_matrix_market/types.hpp>
#include <cstdint>
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

#ifndef FMM_SCIPY_PRUNE
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
#endif

////////////////////////////////////////////////
//// Read cursor - open files/streams for reading
////////////////////////////////////////////////

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


////////////////////////////////////////////////
//// Write cursor - open files/streams writing reading
////////////////////////////////////////////////

write_cursor open_write_file(const std::string& filename, const fmm::matrix_market_header& header,
                             int num_threads, int precision) {
    write_cursor cursor(filename);
    // Set options
    cursor.options.num_threads = num_threads;
    cursor.options.precision = precision;
    cursor.options.always_comment = true; // scipy.io._mmio always writes a comment line, even if comment is empty.
    cursor.header = header;
    return cursor;
}

write_cursor open_write_stream(std::shared_ptr<pystream::ostream>& stream, fmm::matrix_market_header& header,
                               int num_threads, int precision) {
    write_cursor cursor(stream);
    // Set options
    cursor.options.num_threads = num_threads;
    cursor.options.precision = precision;
    cursor.options.always_comment = true; // scipy.io._mmio always writes a comment line, even if comment is empty.
    cursor.header = header;
    return cursor;
}

#ifndef FMM_SCIPY_PRUNE
void write_header_only(write_cursor& cursor) {
    fmm::write_header(cursor.stream(), cursor.header, cursor.options);
    cursor.close();
}
#endif

////////////////////////////////////////////////
//// pybind11 module definition
//// Define the _fmm_core module here, it is used by __init__.py
////////////////////////////////////////////////


PYBIND11_MODULE(_fmm_core, m, py::mod_gil_not_used()) {
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
        } catch (const fmm::support_not_selected& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        } catch (const fmm::fmm_error &e) {
            // Everything else we throw maps best to ValueError
            PyErr_SetString(PyExc_ValueError, e.what());
        }
    });

    py::class_<fmm::matrix_market_header>(m, "header", py::module_local())
#ifndef FMM_SCIPY_PRUNE
    .def(py::init<>())
    .def(py::init<int64_t, int64_t>())
    .def(py::init([](std::tuple<int64_t, int64_t> shape) { return fmm::matrix_market_header{std::get<0>(shape), std::get<1>(shape)}; }))
#endif
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
#ifndef FMM_SCIPY_PRUNE
    .def("to_dict", &header_to_dict, R"pbdoc(
        Return the values in the header as a dict.
    )pbdoc")
    .def("__repr__", [](const fmm::matrix_market_header& header) { return header_repr(header); })
#endif
    ;

#ifndef FMM_SCIPY_PRUNE
    m.def("write_header_only", &write_header_only);
#endif
    ///////////////////////////////
    // Read methods
    py::class_<read_cursor>(m, "_read_cursor", py::module_local())
    .def_readonly("header", &read_cursor::header)
    .def("close", &read_cursor::close);

    m.def("open_read_file", &open_read_file);
    m.def("open_read_stream", &open_read_stream);

    init_read_array(m);
    init_read_coo(m);

    ///////////////////////////////
    // Write methods
    py::class_<write_cursor>(m, "_write_cursor", py::module_local())
#ifndef FMM_SCIPY_PRUNE
    .def_readwrite("header", &write_cursor::header)
#endif
    ;

    m.def("open_write_file", &open_write_file);
    m.def("open_write_stream", &open_write_stream);

    init_write_array(m);
    init_write_coo_32(m);
    init_write_coo_64(m);

#ifndef FMM_SCIPY_PRUNE
    init_write_csc_32(m);
    init_write_csc_64(m);
#endif

    // Module version
#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
