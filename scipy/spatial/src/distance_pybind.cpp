#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <cassert>

#include "function_ref.h"
#include "views.h"
#include "distance_metrics.h"

#include <sstream>
#include <string>

namespace py = pybind11;

namespace {

template <typename T>
using DistanceFunc = FunctionRef<
    void(StridedView2D<T>, StridedView2D<const T>, StridedView2D<const T>)>;

template <typename T>
using WeightedDistanceFunc = FunctionRef<
    void(StridedView2D<T>, StridedView2D<const T>,
         StridedView2D<const T>, StridedView2D<const T>)>;

// Validate weights are >= 0
template <typename T>
void validate_weights(const ArrayDescriptor& w, const T* w_data) {
    intptr_t idx[NPY_MAXDIMS] = {0};
    if (w.ndim > NPY_MAXDIMS) {
        throw std::invalid_argument("Too many dimensions");
    }

    intptr_t numiter = 1;
    for (intptr_t ax = 0; ax < w.ndim - 1; ++ax) {
        numiter *= w.shape[ax];
    }

    bool is_valid = true;
    const T* row_ptr = w_data;
    const auto inner_size = w.shape[w.ndim - 1];
    const auto stride = w.strides[w.ndim - 1];

    while (is_valid && numiter > 0) {
        for (intptr_t i = 0; i < inner_size; ++i) {
            if (row_ptr[i * stride] < 0) {
                is_valid = false;
            }
        }

        for (intptr_t ax = w.ndim - 2; ax >= 0; --ax) {
            if (idx[ax] + 1 < w.shape[ax]) {
                ++idx[ax];
                row_ptr += w.strides[ax];
                break;
            } else {
                row_ptr -= idx[ax] * w.strides[ax];
                idx[ax] = 0;
            }
        }
        --numiter;
    }

    if (!is_valid) {
        throw std::invalid_argument("Input weights should be all non-negative");
    }
}

template <typename T>
void pdist_impl(ArrayDescriptor out, T* out_data,
                ArrayDescriptor x, const T* in_data,
                DistanceFunc<T> f) {
    const intptr_t num_rows = x.shape[0], num_cols = x.shape[1];

    StridedView2D<T> out_view;
    out_view.strides = {out.strides[0], 0};
    out_view.shape = {x.shape[0] - 1, x.shape[1]};
    out_view.data = out_data;

    StridedView2D<const T> x_view;
    x_view.strides = {x.strides[0], x.strides[1]};
    x_view.shape = {out_view.shape[0], num_cols};
    x_view.data = in_data + x.strides[0];

    StridedView2D<const T> y_view;
    y_view.strides = {0, x.strides[1]};
    y_view.shape = {out_view.shape[0], num_cols};
    y_view.data = in_data;

    for (intptr_t i = 0; i < num_rows - 1; ++i) {
        f(out_view, x_view, y_view);

        out_view.data += out_view.shape[0] * out_view.strides[0];
        out_view.shape[0] -= 1;
        x_view.shape[0] = y_view.shape[0] = out_view.shape[0];
        x_view.data += x.strides[0];
        y_view.data += x.strides[0];
    }
}

template <typename T>
void pdist_weighted_impl(ArrayDescriptor out, T* out_data,
                         ArrayDescriptor x, const T* x_data,
                         ArrayDescriptor w, const T* w_data,
                         WeightedDistanceFunc<T> f) {
    if (x.ndim != 2) {
        throw std::invalid_argument("x must be 2-dimensional");
    }

    StridedView2D<T> out_view;
    out_view.strides = {out.strides[0], 0};
    out_view.shape = {x.shape[0] - 1, x.shape[1]};
    out_view.data = out_data;

    StridedView2D<const T> w_view;
    w_view.strides = {0, w.strides[0]};
    w_view.shape = out_view.shape;
    w_view.data = w_data;

    StridedView2D<const T> x_view;
    x_view.strides = {x.strides[0], x.strides[1]};
    x_view.shape = out_view.shape;
    x_view.data = x_data + x.strides[0];

    StridedView2D<const T> y_view;
    y_view.strides = {0, x.strides[1]};
    y_view.shape = out_view.shape;
    y_view.data = x_data;

    const intptr_t num_rows = x.shape[0];
    for (intptr_t i = 0; i < num_rows - 1; ++i) {
        f(out_view, x_view, y_view, w_view);

        out_view.data += out_view.shape[0] * out_view.strides[0];
        out_view.shape[0] -= 1;
        x_view.shape[0] = y_view.shape[0] = w_view.shape[0] = out_view.shape[0];
        x_view.data += x.strides[0];
        y_view.data += x.strides[0];
    }
}

template <typename T>
void cdist_impl(ArrayDescriptor out, T* out_data,
                ArrayDescriptor x, const T* x_data,
                ArrayDescriptor y, const T* y_data,
                DistanceFunc<T> f) {

    const auto num_rowsX = x.shape[0];
    const auto num_rowsY = y.shape[0];
    const auto num_cols = x.shape[1];

    StridedView2D<T> out_view;
    out_view.strides = {out.strides[1], 0};
    out_view.shape = {num_rowsY, num_cols};
    out_view.data = out_data;

    StridedView2D<const T> x_view;
    x_view.strides = {0, x.strides[1]};
    x_view.shape = {num_rowsY, num_cols};
    x_view.data = x_data;

    StridedView2D<const T> y_view;
    y_view.strides = {y.strides[0], y.strides[1]};
    y_view.shape = {out_view.shape[0], num_cols};
    y_view.data = y_data;

    for (intptr_t i = 0; i < num_rowsX; ++i) {
        f(out_view, x_view, y_view);

        out_view.data += out.strides[0];
        x_view.data += x.strides[0];
    }
}

template <typename T>
void cdist_weighted_impl(ArrayDescriptor out, T* out_data,
                         ArrayDescriptor x, const T* x_data,
                         ArrayDescriptor y, const T* y_data,
                         ArrayDescriptor w, const T* w_data,
                         WeightedDistanceFunc<T> f) {

    const auto num_rowsX = x.shape[0];
    const auto num_rowsY = y.shape[0];
    const auto num_cols = x.shape[1];

    StridedView2D<T> out_view;
    out_view.strides = {out.strides[1], 0};
    out_view.shape = {num_rowsY, num_cols};
    out_view.data = out_data;

    StridedView2D<const T> x_view;
    x_view.strides = {0, x.strides[1]};
    x_view.shape = {num_rowsY, num_cols};
    x_view.data = x_data;

    StridedView2D<const T> y_view;
    y_view.strides = {y.strides[0], y.strides[1]};
    y_view.shape = {num_rowsY, num_cols};
    y_view.data = y_data;

    StridedView2D<const T> w_view;
    w_view.strides = {0, w.strides[0]};
    w_view.shape = {num_rowsY, num_cols};
    w_view.data = w_data;

    for (intptr_t i = 0; i < num_rowsX; ++i) {
        f(out_view, x_view, y_view, w_view);

        out_view.data += out.strides[0];
        x_view.data += x.strides[0];
    }
}

// Extract shape and stride information from NumPy array. Converts byte-strides
// to element strides, and avoids an extra pointer indirection on access.
ArrayDescriptor get_descriptor(const py::array& arr) {
    const auto ndim = arr.ndim();
    ArrayDescriptor desc(ndim);

    const auto arr_shape = arr.shape();
    desc.shape.assign(arr_shape, arr_shape + ndim);

    desc.element_size = arr.itemsize();
    const auto arr_strides = arr.strides();
    desc.strides.assign(arr_strides, arr_strides + ndim);
    for (intptr_t i = 0; i < ndim; ++i) {
        if (desc.strides[i] % desc.element_size != 0) {
            throw std::runtime_error("Arrays must be aligned");
        }
        desc.strides[i] /= desc.element_size;
    }
    return desc;
}

// Cast python object to NumPy array of data type T.
// flags can be any NumPy array constructor flags.
template <typename T>
py::array_t<T> npy_asarray(const py::handle& obj, int flags = 0) {
    auto descr = reinterpret_cast<PyArray_Descr*>(
        py::dtype::of<T>().release().ptr());
    auto* arr = PyArray_FromAny(obj.ptr(), descr, 0, 0, flags, nullptr);
    if (arr == nullptr) {
        throw py::error_already_set();
    }
    return py::reinterpret_steal<py::array_t<T>>(arr);
}

// Cast python object to NumPy array with unspecified dtype.
// flags can be any NumPy array constructor flags.
py::array npy_asarray(const py::handle& obj, int flags = 0) {
    auto* arr = PyArray_FromAny(obj.ptr(), nullptr, 0, 0, flags, nullptr);
    if (arr == nullptr) {
        throw py::error_already_set();
    }
    return py::reinterpret_steal<py::array>(arr);
}

template <typename scalar_t>
py::array pdist_unweighted(const py::array& out_obj, const py::array& x_obj,
                           DistanceFunc<scalar_t> f) {
    auto x = npy_asarray<scalar_t>(x_obj,
                                   NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out = py::cast<py::array_t<scalar_t>>(out_obj);
    auto out_desc = get_descriptor(out);
    auto out_data = out.mutable_data();
    auto x_desc = get_descriptor(x);
    auto x_data = x.data();
    {
        py::gil_scoped_release guard;
        pdist_impl(out_desc, out_data, x_desc, x_data, f);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array pdist_weighted(
        const py::array& out_obj, const py::array& x_obj,
        const py::array& w_obj, WeightedDistanceFunc<scalar_t> f) {
    auto x = npy_asarray<scalar_t>(x_obj,
                                   NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto w = npy_asarray<scalar_t>(w_obj,
                                   NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out = py::cast<py::array_t<scalar_t>>(out_obj);
    auto out_desc = get_descriptor(out);
    auto out_data = out.mutable_data();
    auto x_desc = get_descriptor(x);
    auto x_data = x.data();
    auto w_desc = get_descriptor(w);
    auto w_data = w.data();
    {
        py::gil_scoped_release guard;
        validate_weights(w_desc, w_data);
        pdist_weighted_impl(
            out_desc, out_data, x_desc, x_data, w_desc, w_data, f);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array cdist_unweighted(const py::array& out_obj, const py::array& x_obj,
                           const py::array& y_obj, DistanceFunc<scalar_t> f) {
    auto x = npy_asarray<scalar_t>(x_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto y = npy_asarray<scalar_t>(y_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out = py::cast<py::array_t<scalar_t>>(out_obj);

    auto out_desc = get_descriptor(out);
    auto out_data = out.mutable_data();
    auto x_desc = get_descriptor(x);
    auto x_data = x.data();
    auto y_desc = get_descriptor(y);
    auto y_data = y.data();
    {
        py::gil_scoped_release guard;
        cdist_impl(out_desc, out_data, x_desc, x_data, y_desc, y_data, f);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array cdist_weighted(
        const py::array& out_obj, const py::array& x_obj,
        const py::array& y_obj, const py::array& w_obj,
        WeightedDistanceFunc<scalar_t> f) {
    auto x = npy_asarray<scalar_t>(x_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto y = npy_asarray<scalar_t>(y_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto w = npy_asarray<scalar_t>(w_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out = py::cast<py::array_t<scalar_t>>(out_obj);

    auto out_desc = get_descriptor(out);
    auto out_data = out.mutable_data();
    auto x_desc = get_descriptor(x);
    auto x_data = x.data();
    auto y_desc = get_descriptor(y);
    auto y_data = y.data();
    auto w_desc = get_descriptor(w);
    auto w_data = w.data();
    {
        py::gil_scoped_release guard;
        validate_weights(w_desc, w_data);
        cdist_weighted_impl(
            out_desc, out_data, x_desc, x_data, y_desc, y_data, w_desc, w_data, f);
    }
    return std::move(out);
}

py::dtype npy_promote_types(const py::dtype& type1, const py::dtype& type2) {
    PyArray_Descr* descr = PyArray_PromoteTypes(
        reinterpret_cast<PyArray_Descr*>(type1.ptr()),
        reinterpret_cast<PyArray_Descr*>(type2.ptr()));
    if (descr == nullptr) {
        throw py::error_already_set();
    }
    return py::reinterpret_steal<py::dtype>(reinterpret_cast<PyObject*>(descr));
}

template <typename Container>
py::array prepare_out_argument(const py::object& obj, const py::dtype& dtype,
                               const Container& out_shape) {
    if (obj.is_none()) {
        return py::array(dtype, out_shape);
    }

    if (!py::isinstance<py::array>(obj)) {
        throw py::type_error("out argument must be an ndarray");
    }

    py::array out = py::cast<py::array>(obj);
    const auto ndim = out.ndim();
    const auto shape = out.shape();
    auto pao = reinterpret_cast<PyArrayObject*>(out.ptr());

    if (ndim != static_cast<intptr_t>(out_shape.size()) ||
        !std::equal(shape, shape + ndim, out_shape.begin())) {
        throw std::invalid_argument("Output array has incorrect shape.");
    }
    if (!PyArray_ISCONTIGUOUS(pao)) {
        throw std::invalid_argument("Output array must be C-contiguous");
    }
    if (out.dtype().not_equal(dtype)) {
        const py::handle& handle = dtype;
        throw std::invalid_argument("wrong out dtype, expected " +
                                    std::string(py::str(handle)));
    }
    if (!PyArray_ISBEHAVED(pao)) {
        throw std::invalid_argument(
            "out array must be aligned, writable and native byte order");
    }
    return out;
}

py::array prepare_single_weight(const py::object& obj, intptr_t len) {
    py::array weight = npy_asarray(obj);
    if (weight.ndim() != 1) {
        throw std::invalid_argument("Weights must be a vector (ndim = 1)");
    } else if (weight.shape(0) != len) {
        std::stringstream msg;
        msg << "Weights must have same size as input vector. ";
        msg << weight.shape(0) << " vs. " << len << ".";
        throw std::invalid_argument(msg.str());
    }
    return weight;
}

py::dtype common_type(py::dtype type) {
    return type;
}

template <typename... Args>
py::dtype common_type(const py::dtype& type1, const py::dtype& type2,
                      const Args&... tail) {
    return common_type(npy_promote_types(type1, type2), tail...);
}

int dtype_num(const py::dtype& dtype) {
    return reinterpret_cast<const PyArray_Descr*>(
        dtype.ptr())->type_num;
}

py::dtype promote_type_real(const py::dtype& dtype) {
    switch (dtype.kind()) {
    case 'b':
    case 'i':
    case 'u': {
        // Promote integral and boolean types to double
        return py::dtype::template of<double>();
    }
    case 'f': {
        if (dtype_num(dtype) == NPY_LONGDOUBLE) {
            return dtype;
        } else {
            // TODO: Allow float32 output
            return py::dtype::template of<double>();
        }
    }

    default: {
        return dtype;
    }
    }
}

// From a NumPy dtype, run "expression" with scalar_t aliasing the C++ type
#define DISPATCH_DTYPE(dtype, expression)                               \
    do {                                                                \
        const py::dtype& type_obj = dtype;                              \
        switch (dtype_num(type_obj)) {                                  \
        case NPY_HALF:                                                  \
        case NPY_FLOAT: /* TODO: Enable scalar_t=float dispatch */      \
        case NPY_DOUBLE: {                                              \
            using scalar_t = double;                                    \
            expression();                                               \
            break;                                                      \
        }                                                               \
        case NPY_LONGDOUBLE: {                                          \
            using scalar_t = long double;                               \
            expression();                                               \
            break;                                                      \
        }                                                               \
        default: {                                                      \
            const py::handle& handle = type_obj;                        \
            throw std::invalid_argument(                                \
                "Unsupported dtype " + std::string(py::str(handle)));   \
        }                                                               \
        }                                                               \
    } while (0)

template <typename Func>
py::array pdist(const py::object& out_obj, const py::object& x_obj,
                const py::object& w_obj, Func&& f) {
    auto x = npy_asarray(x_obj);
    if (x.ndim() != 2) {
        throw std::invalid_argument("x must be 2-dimensional");
    }

    const intptr_t m = x.shape(1);
    const intptr_t n = x.shape(0);
    std::array<intptr_t, 1> out_shape{{(n * (n - 1)) / 2}};
    if (w_obj.is_none()) {
        auto dtype = promote_type_real(x.dtype());
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            pdist_unweighted<scalar_t>(out, x, f);
        });
        return out;
    }

    auto w = prepare_single_weight(w_obj, m);
    auto dtype = promote_type_real(common_type(x.dtype(), w.dtype()));
    auto out = prepare_out_argument(out_obj, dtype, out_shape);
    DISPATCH_DTYPE(dtype, [&]{
        pdist_weighted<scalar_t>(out, x, w, f);
    });
    return out;
}

template <typename Func>
py::array cdist(const py::object& out_obj, const py::object& x_obj,
                const py::object& y_obj, const py::object& w_obj, Func&& f) {
    auto x = npy_asarray(x_obj);
    auto y = npy_asarray(y_obj);
    if (x.ndim() != 2) {
        throw std::invalid_argument("XA must be a 2-dimensional array.");
    }
    if (y.ndim() != 2) {
        throw std::invalid_argument("XB must be a 2-dimensional array.");
    }
    const intptr_t m = x.shape(1);
    if (m != y.shape(1)) {
        throw std::invalid_argument(
            "XA and XB must have the same number of columns "
            "(i.e. feature dimension).");
    }

    std::array<intptr_t, 2> out_shape{{x.shape(0), y.shape(0)}};
    if (w_obj.is_none()) {
        auto dtype = promote_type_real(common_type(x.dtype(), y.dtype()));
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            cdist_unweighted<scalar_t>(out, x, y, f);
        });
        return out;
    }

    auto w = prepare_single_weight(w_obj, m);
    auto dtype = promote_type_real(
        common_type(x.dtype(), y.dtype(), w.dtype()));
    auto out = prepare_out_argument(out_obj, dtype, out_shape);
    DISPATCH_DTYPE(dtype, [&]{
        cdist_weighted<scalar_t>(out, x, y, w, f);
    });
    return out;
}

PYBIND11_MODULE(_distance_pybind, m) {
    if (_import_array() != 0) {
        throw py::error_already_set();
    }
    using namespace pybind11::literals;
    m.def("pdist_canberra",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, CanberraDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_chebyshev",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, ChebyshevDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_cityblock",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, CityBlockDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_euclidean",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, EuclideanDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_minkowski",
          [](py::object x, py::object w, py::object out, double p) {
              if (p == 1.0) {
                  return pdist(out, x, w, CityBlockDistance{});
              } else if (p == 2.0) {
                  return pdist(out, x, w, EuclideanDistance{});
              } else if (std::isinf(p)) {
                  return pdist(out, x, w, ChebyshevDistance{});
              } else {
                  return pdist(out, x, w, MinkowskiDistance{p});
              }
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none(), "p"_a=2.0);
    m.def("pdist_sqeuclidean",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, SquareEuclideanDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_braycurtis",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, BraycurtisDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_canberra",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, CanberraDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_chebyshev",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, ChebyshevDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_cityblock",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, CityBlockDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_euclidean",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, EuclideanDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_minkowski",
          [](py::object x, py::object y, py::object w, py::object out,
             double p) {
              if (p == 1.0) {
                  return cdist(out, x, y, w, CityBlockDistance{});
              } else if (p == 2.0) {
                  return cdist(out, x, y, w, EuclideanDistance{});
              } else if (std::isinf(p)) {
                  return cdist(out, x, y, w, ChebyshevDistance{});
              } else {
                  return cdist(out, x, y, w, MinkowskiDistance{p});
              }
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none(), "p"_a=2.0);
    m.def("cdist_sqeuclidean",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, SquareEuclideanDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_braycurtis",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, BraycurtisDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
}

}  // namespace (anonymous)
