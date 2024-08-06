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

enum PreprocessingType {
    None, Normalise, CentraliseAndNormalise, ComputeVariance, ScaleInputForJensenshannon_
};

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
ALWAYS_INLINE void preprocess_inputs(
    ArrayDescriptor x, const T* in_data, PreprocessingType type) {
    if( type == PreprocessingType::None ) {
        return ;
    }

    T* in_data_ = const_cast<T*>(in_data);

    StridedView2D<T> x_view;
    x_view.strides = {x.strides[0], x.strides[1]};
    x_view.shape = {x.shape[0], x.shape[1]};
    x_view.data = in_data_;

    switch( type ) {
        case PreprocessingType::Normalise: {
            NormaliseInput{}(x_view);
            break;
        }
        case PreprocessingType::CentraliseAndNormalise: {
            CentraliseInput{}(x_view);
            NormaliseInput{}(x_view);
            break;
        }
        case PreprocessingType::ScaleInputForJensenshannon_: {
            ScaleInputForJensenshannon{}(x_view);
            break;
        }
        default: {
        }
    }
}

template <typename T>
ALWAYS_INLINE void preprocess_inputs(
    ArrayDescriptor x, const T* in_data,
    StridedView2D<const T> w_view, PreprocessingType type) {
    if( type == PreprocessingType::None ) {
        return ;
    }

    T* in_data_ = const_cast<T*>(in_data);

    StridedView2D<T> x_view;
    x_view.strides = {x.strides[0], x.strides[1]};
    x_view.shape = {x.shape[0], x.shape[1]};
    x_view.data = in_data_;

    switch( type ) {
        case PreprocessingType::Normalise: {
            NormaliseInput{}(x_view, w_view);
            break;
        }
        case PreprocessingType::CentraliseAndNormalise: {
            CentraliseInput{}(x_view, w_view);
            NormaliseInput{}(x_view, w_view);
            break;
        }
        default: {
        }
    }
}

template <typename T>
void pdist_impl(ArrayDescriptor out, T* out_data,
                ArrayDescriptor x, const T* in_data,
                DistanceFunc<T> f,
                PreprocessingType type) {
    const intptr_t num_rows = x.shape[0], num_cols = x.shape[1];

    preprocess_inputs(x, in_data, type);

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
                         WeightedDistanceFunc<T> f,
                         PreprocessingType type) {
    if (x.ndim != 2) {
        throw std::invalid_argument("x must be 2-dimensional");
    }

    const intptr_t num_rows = x.shape[0];

    StridedView2D<T> out_view;
    out_view.strides = {out.strides[0], 0};
    out_view.shape = {x.shape[0] - 1, x.shape[1]};
    out_view.data = out_data;

    StridedView2D<const T> w_view;
    w_view.strides = {0, w.strides[0]};
    w_view.shape = out_view.shape;
    w_view.data = w_data;

    preprocess_inputs(x, x_data, w_view, type);

    StridedView2D<const T> x_view;
    x_view.strides = {x.strides[0], x.strides[1]};
    x_view.shape = out_view.shape;
    x_view.data = x_data + x.strides[0];

    StridedView2D<const T> y_view;
    y_view.strides = {0, x.strides[1]};
    y_view.shape = out_view.shape;
    y_view.data = x_data;

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
ALWAYS_INLINE void _compute_variance_across_rows(
    T* var, ArrayDescriptor x, const T* x_data,
    ArrayDescriptor y, const T* y_data) {
    const intptr_t num_rowsX = x.shape[0], num_rowsY = y.shape[0];
    const intptr_t num_cols = x.shape[1];

    for( intptr_t j = 0; j < num_cols; j++ ) {
        T sum[4] = {0, 0, 0, 0};
        intptr_t i;
        for( i = 0; i + 3 < num_rowsX; i += 4 ) {
            const intptr_t index0 = i * x.strides[0] + j * x.strides[1];
            const intptr_t index1 = (i + 1) * x.strides[0] + j * x.strides[1];
            const intptr_t index2 = (i + 2) * x.strides[0] + j * x.strides[1];
            const intptr_t index3 = (i + 3) * x.strides[0] + j * x.strides[1];
            sum[0] += x_data[index0];
            sum[1] += x_data[index1];
            sum[2] += x_data[index2];
            sum[3] += x_data[index3];
        }
        for( ; i < num_rowsX; i++ ) {
            const intptr_t index = i * x.strides[0] + j * x.strides[1];
            sum[0] += x_data[index];
        }

        for( i = 0; i + 3 < num_rowsY; i += 4 ) {
            const intptr_t index0 = i * y.strides[0] + j * y.strides[1];
            const intptr_t index1 = (i + 1) * y.strides[0] + j * y.strides[1];
            const intptr_t index2 = (i + 2) * y.strides[0] + j * y.strides[1];
            const intptr_t index3 = (i + 3) * y.strides[0] + j * y.strides[1];
            sum[0] += y_data[index0];
            sum[1] += y_data[index1];
            sum[2] += y_data[index2];
            sum[3] += y_data[index3];
        }
        for( ; i < num_rowsY; i++ ) {
            const intptr_t index = i * y.strides[0] + j * y.strides[1];
            sum[0] += y_data[index];
        }
        var[j] = (sum[0] + sum[1] + sum[2] + sum[3])/(num_rowsX + num_rowsY);
    }

    for( intptr_t j = 0; j < num_cols; j++ ) {
        T sum[2] = {0, 0};
        T mean = var[j];
        intptr_t i;
        for( i = 0; i + 1 < num_rowsX; i += 2 ) {
            const intptr_t index0 = i * x.strides[0] + j * x.strides[1];
            const intptr_t index1 = (i + 1) * x.strides[0] + j * x.strides[1];
            sum[0] += (x_data[index0] - mean) * (x_data[index0] - mean);
            sum[1] += (x_data[index1] - mean) * (x_data[index1] - mean);
        }
        for( ; i < num_rowsX; i++ ) {
            const intptr_t index = i * x.strides[0] + j * x.strides[1];
            sum[0] += (x_data[index] - mean) * (x_data[index] - mean);
        }
        for( i = 0; i + 1 < num_rowsY; i += 2 ) {
            const intptr_t index0 = i * y.strides[0] + j * y.strides[1];
            const intptr_t index1 = (i + 1) * y.strides[0] + j * y.strides[1];
            sum[0] += (y_data[index0] - mean) * (y_data[index0] - mean);
            sum[1] += (y_data[index1] - mean) * (y_data[index1] - mean);
        }
        for( ; i < num_rowsY; i++ ) {
            const intptr_t index = i * y.strides[0] + j * y.strides[1];
            sum[0] += (y_data[index] - mean) * (y_data[index] - mean);
        }
        var[j] = (num_rowsX + num_rowsY - 1)/(sum[0] + sum[1]);
    }
}

template <typename T>
ALWAYS_INLINE void _compute_variance_across_rows(
    T* var, ArrayDescriptor x, const T* x_data) {
    const intptr_t num_rowsX = x.shape[0];
    const intptr_t num_cols = x.shape[1];

    for( intptr_t j = 0; j < num_cols; j++ ) {
        T sum[4] = {0, 0, 0, 0};
        intptr_t i;
        for( i = 0; i + 3 < num_rowsX; i += 4 ) {
            const intptr_t index0 = i * x.strides[0] + j * x.strides[1];
            const intptr_t index1 = (i + 1) * x.strides[0] + j * x.strides[1];
            const intptr_t index2 = (i + 2) * x.strides[0] + j * x.strides[1];
            const intptr_t index3 = (i + 3) * x.strides[0] + j * x.strides[1];
            sum[0] += x_data[index0];
            sum[1] += x_data[index1];
            sum[2] += x_data[index2];
            sum[3] += x_data[index3];
        }
        for( ; i < num_rowsX; i++ ) {
            const intptr_t index = i * x.strides[0] + j * x.strides[1];
            sum[0] += x_data[index];
        }

        var[j] = (sum[0] + sum[1] + sum[2] + sum[3])/num_rowsX;
    }

    for( intptr_t j = 0; j < num_cols; j++ ) {
        T sum[2] = {0, 0};
        T mean = var[j];
        intptr_t i;
        for( i = 0; i + 1 < num_rowsX; i += 2 ) {
            const intptr_t index0 = i * x.strides[0] + j * x.strides[1];
            const intptr_t index1 = (i + 1) * x.strides[0] + j * x.strides[1];
            sum[0] += (x_data[index0] - mean) * (x_data[index0] - mean);
            sum[1] += (x_data[index1] - mean) * (x_data[index1] - mean);
        }
        for( ; i < num_rowsX; i++ ) {
            const intptr_t index = i * x.strides[0] + j * x.strides[1];
            sum[0] += (x_data[index] - mean) * (x_data[index] - mean);
        }
        var[j] = (num_rowsX - 1)/(sum[0] + sum[1]);
    }
}

template <typename T>
void cdist_impl(ArrayDescriptor out, T* out_data,
                ArrayDescriptor x, const T* x_data,
                ArrayDescriptor y, const T* y_data,
                DistanceFunc<T> f, PreprocessingType type) {

    const auto num_rowsX = x.shape[0];
    const auto num_rowsY = y.shape[0];
    const auto num_cols = x.shape[1];

    preprocess_inputs(x, x_data, type);
    preprocess_inputs(y, y_data, type);

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
                         WeightedDistanceFunc<T> f,
                         PreprocessingType type) {

    const auto num_rowsX = x.shape[0];
    const auto num_rowsY = y.shape[0];
    const auto num_cols = x.shape[1];

    StridedView2D<const T> w_view;
    w_view.strides = {0, w.strides[0]};
    w_view.shape = {num_rowsY, num_cols};
    w_view.data = w_data;

    preprocess_inputs(x, x_data, w_view, type);
    preprocess_inputs(y, y_data, w_view, type);

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
        if (arr_shape[i] <= 1) {
            // Under NumPy's relaxed stride checking, dimensions with
            // 1 or fewer elements are ignored.
            desc.strides[i] = 0;
            continue;
        }

        if (desc.strides[i] % desc.element_size != 0) {
            std::stringstream msg;
            msg << "Arrays must be aligned to element size, but found stride of ";
            msg << desc.strides[i] << " bytes for elements of size " << desc.element_size;
            throw std::runtime_error(msg.str());
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
                           DistanceFunc<scalar_t> f,
                           PreprocessingType type) {
    auto x = npy_asarray<scalar_t>(x_obj,
                                   NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out = py::cast<py::array_t<scalar_t>>(out_obj);
    auto out_desc = get_descriptor(out);
    auto out_data = out.mutable_data();
    auto x_desc = get_descriptor(x);
    auto x_data = x.data();
    {
        py::gil_scoped_release guard;
        pdist_impl(out_desc, out_data, x_desc, x_data, f, type);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array pdist_weighted(
        const py::array& out_obj, const py::array& x_obj,
        const py::array& w_obj, WeightedDistanceFunc<scalar_t> f,
        PreprocessingType type) {
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
            out_desc, out_data, x_desc, x_data, w_desc, w_data, f, type);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array pdist_weighted(
        const py::array& out_obj, const py::array& x_obj,
        WeightedDistanceFunc<scalar_t> f) {
    auto x = npy_asarray<scalar_t>(x_obj,
                                   NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out = py::cast<py::array_t<scalar_t>>(out_obj);
    auto out_desc = get_descriptor(out);
    auto out_data = out.mutable_data();
    auto x_desc = get_descriptor(x);
    auto x_data = x.data();

    scalar_t* w_data = new scalar_t[x_desc.shape[1]];
    _compute_variance_across_rows(w_data, x_desc, x_data);
    ArrayDescriptor w_desc(1);
    w_desc.element_size = x_desc.element_size;
    w_desc.shape[0] = x_desc.shape[1];
    w_desc.strides[0] = 1;

    {
        py::gil_scoped_release guard;
        pdist_weighted_impl(
            out_desc, out_data, x_desc, x_data, w_desc,
            w_data, f, PreprocessingType::None);
    }

    delete [] w_data;
    return std::move(out);
}

template <typename scalar_t>
py::array cdist_unweighted(const py::array& out_obj, const py::array& x_obj,
                           const py::array& y_obj, DistanceFunc<scalar_t> f,
                           PreprocessingType type) {
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
        cdist_impl(out_desc, out_data, x_desc, x_data, y_desc, y_data, f, type);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array cdist_weighted(
        const py::array& out_obj, const py::array& x_obj,
        const py::array& y_obj, const py::array& w_obj,
        WeightedDistanceFunc<scalar_t> f,
        PreprocessingType type) {
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
            out_desc, out_data, x_desc, x_data, y_desc, y_data,
            w_desc, w_data, f, type);
    }
    return std::move(out);
}

template <typename scalar_t>
py::array cdist_weighted(
        const py::array& out_obj, const py::array& x_obj,
        const py::array& y_obj, WeightedDistanceFunc<scalar_t> f) {
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

    scalar_t* w_data = new scalar_t[x_desc.shape[1]];
    _compute_variance_across_rows(w_data, x_desc, x_data, y_desc, y_data);
    ArrayDescriptor w_desc(1);
    w_desc.element_size = x_desc.element_size;
    w_desc.shape[0] = x_desc.shape[1];
    w_desc.strides[0] = 1;

    {
        py::gil_scoped_release guard;
        cdist_weighted_impl(
            out_desc, out_data, x_desc, x_data, y_desc, y_data,
            w_desc, w_data, f, PreprocessingType::None);
    }

    delete [] w_data;
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
                const py::object& w_obj, Func&& f,
                PreprocessingType type=PreprocessingType::None) {
    auto x = npy_asarray(x_obj);
    if (x.ndim() != 2) {
        throw std::invalid_argument("x must be 2-dimensional");
    }

    const intptr_t m = x.shape(1);
    const intptr_t n = x.shape(0);
    std::array<intptr_t, 1> out_shape{{(n * (n - 1)) / 2}};

    if( type == PreprocessingType::ComputeVariance ) {
        auto dtype = promote_type_real(common_type(x.dtype()));
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            pdist_weighted<scalar_t>(out, x, f);
        });
        return out;
    }

    if (w_obj.is_none()) {
        auto dtype = promote_type_real(x.dtype());
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            pdist_unweighted<scalar_t>(out, x, f, type);
        });
        return out;
    }

    auto w = prepare_single_weight(w_obj, m);
    auto dtype = promote_type_real(common_type(x.dtype(), w.dtype()));
    auto out = prepare_out_argument(out_obj, dtype, out_shape);
    DISPATCH_DTYPE(dtype, [&]{
        pdist_weighted<scalar_t>(out, x, w, f, type);
    });
    return out;
}

template <typename Func>
py::array cdist(const py::object& out_obj, const py::object& x_obj,
                const py::object& y_obj, const py::object& w_obj, Func&& f,
                PreprocessingType type=PreprocessingType::None) {
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

    if( type == PreprocessingType::ComputeVariance ) {
        auto dtype = promote_type_real(common_type(x.dtype(), y.dtype()));
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            cdist_weighted<scalar_t>(out, x, y, f);
        });
        return out;
    }

    if (w_obj.is_none()) {
        auto dtype = promote_type_real(common_type(x.dtype(), y.dtype()));
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            cdist_unweighted<scalar_t>(out, x, y, f, type);
        });
        return out;
    }

    auto w = prepare_single_weight(w_obj, m);
    auto dtype = promote_type_real(
        common_type(x.dtype(), y.dtype(), w.dtype()));
    auto out = prepare_out_argument(out_obj, dtype, out_shape);
    DISPATCH_DTYPE(dtype, [&]{
        cdist_weighted<scalar_t>(out, x, y, w, f, type);
    });
    return out;
}

template <typename scalar_t>
ALWAYS_INLINE void cdist_mahalanobis_impl(ArrayDescriptor x, const scalar_t* x_data,
    ArrayDescriptor y, const scalar_t* y_data, ArrayDescriptor vi, const scalar_t* vi_data,
    ArrayDescriptor out, scalar_t* out_data, intptr_t vi_col_offset) {
    const intptr_t num_rowsX = x.shape[0], num_rowsY = y.shape[0];
    const intptr_t num_cols = x.shape[1];
    scalar_t* uv_diff1 = new scalar_t[2*num_cols];
    scalar_t* uv_diff2 = uv_diff1 + num_cols;

    for( intptr_t i = 0; i < num_rowsX; i++ ) {
        const scalar_t* u = &x_data[i * x.strides[0]];
        for( intptr_t j = 0; j < num_rowsY; j++ ) {
            const scalar_t* v = &y_data[j * y.strides[0]];

            for( intptr_t k = 0; k < num_cols; k++ ) {
                uv_diff1[k] = u[k*x.strides[1]] - v[k*y.strides[1]];
            }

            for( intptr_t k = 0; k < num_cols; k++ ) {
                uv_diff2[k] = 0;
                for( intptr_t l = 0; l < num_cols; l++ ) {
                    uv_diff2[k] += uv_diff1[l] * vi_data[k*vi.strides[0] + (l + vi_col_offset)*vi.strides[1]];
                }
            }

            scalar_t s = 0;
            for( intptr_t k = 0; k < num_cols; k++ ) {
                s += uv_diff1[k] * uv_diff2[k];
            }

            out_data[i*out.strides[0] + j*out.strides[1]] = std::sqrt(s);
        }
    }

    delete [] uv_diff1;
}

template <typename scalar_t>
ALWAYS_INLINE void _compute_inverse_matrix(ArrayDescriptor vi, scalar_t* vi_data) {

    const intptr_t num_cols = vi.shape[0];

    for( intptr_t i = 0; i < num_cols; i++ ) {
        for( intptr_t j = num_cols; j < 2*num_cols; j++ ) {
            vi_data[i*vi.strides[0] + j*vi.strides[1]] = 0;
        }
    }

    for( intptr_t i = 0; i < num_cols; i++ ) {
        for( intptr_t j = 0; j < 2*num_cols; j++ ) {
            if( (j - i) == num_cols ) {
                vi_data[i*vi.strides[0] + j*vi.strides[1]] = 1;
            }
        }
    }

    for( intptr_t i = num_cols - 1; i >= 1; i-- ) {
        if( vi_data[(i - 1)*vi.strides[0]] < vi_data[i*vi.strides[0]] ) {
            for( intptr_t j = 0; j < 2*num_cols; j++ ) {
                scalar_t tmp = vi_data[i*vi.strides[0] + j*vi.strides[1]];
                vi_data[i*vi.strides[0] + j*vi.strides[1]] = vi_data[(i - 1)*vi.strides[0] + j*vi.strides[1]];
                vi_data[(i - 1)*vi.strides[0] + j*vi.strides[1]] = tmp;
            }
        }
    }

    for( intptr_t i = 0; i < num_cols; i++ ) {
        for( intptr_t j = 0; j < num_cols; j++ ) {
            if( j != i ) {
                const intptr_t ji = j*vi.strides[0] + i*vi.strides[1];
                const intptr_t ii = i*vi.strides[0] + i*vi.strides[1];
                scalar_t factor = vi_data[ji] / vi_data[ii];
                for( intptr_t k = 0; k < 2*num_cols; k++ ) {
                    const intptr_t jk = j*vi.strides[0] + k*vi.strides[1];
                    const intptr_t ik = i*vi.strides[0] + k*vi.strides[1];
                    vi_data[jk] -= vi_data[ik] * factor;
                }
            }
        }
    }

    for( intptr_t i = 0; i < num_cols; i++ ) {
        const intptr_t ii = i*vi.strides[0] + i*vi.strides[1];
        scalar_t factor = vi_data[ii];
        for( intptr_t j = 0; j < 2*num_cols; j++ ) {
            const intptr_t ij = i*vi.strides[0] + j*vi.strides[1];
            vi_data[ij] /= factor;
        }
    }
}

template <typename scalar_t>
ALWAYS_INLINE scalar_t _compute_mean(
    ArrayDescriptor x[], const scalar_t* x_data[], intptr_t nx, intptr_t j) {
    scalar_t sum[4] = {0, 0, 0, 0};
    intptr_t total_values = 0;

    for( intptr_t ix = 0; ix < nx; ix++ ) {
        const intptr_t num_rows = x[ix].shape[0];
        total_values += num_rows;
        intptr_t i;
        for( i = 0; i + 3 < num_rows; i += 4 ) {
            auto index0 = i*x[ix].strides[0] + j*x[ix].strides[1];
            auto index1 = (i + 1)*x[ix].strides[0] + j*x[ix].strides[1];
            auto index2 = (i + 2)*x[ix].strides[0] + j*x[ix].strides[1];
            auto index3 = (i + 3)*x[ix].strides[0] + j*x[ix].strides[1];
            sum[0] += x_data[ix][index0];
            sum[1] += x_data[ix][index1];
            sum[2] += x_data[ix][index2];
            sum[3] += x_data[ix][index3];
        }
        for( ; i < num_rows; i++ ) {
            auto index = i*x[ix].strides[0] + j*x[ix].strides[1];
            sum[0] += x_data[ix][index];
        }
    }
    return (sum[0] + sum[1] + sum[2] + sum[3])/total_values;
}

template <typename scalar_t>
ALWAYS_INLINE scalar_t _compute_covariance(
    ArrayDescriptor x[], const scalar_t* x_data[],
    intptr_t nx, intptr_t j, intptr_t k) {
    scalar_t j_mean = _compute_mean(x, x_data, nx, j);
    scalar_t k_mean = _compute_mean(x, x_data, nx, k);
    scalar_t cov_jk[4] = {0, 0, 0, 0};
    intptr_t total_values = 0;

    for( intptr_t ix = 0; ix < nx; ix++ ) {
        const intptr_t num_rows = x[ix].shape[0];
        total_values += num_rows;
        intptr_t i;
        for( i = 0; i + 3 < num_rows; i += 4 ) {
            auto jindex0 = i*x[ix].strides[0] + j*x[ix].strides[1];
            auto kindex0 = i*x[ix].strides[0] + k*x[ix].strides[1];
            auto jindex1 = (i + 1)*x[ix].strides[0] + j*x[ix].strides[1];
            auto kindex1 = (i + 1)*x[ix].strides[0] + k*x[ix].strides[1];
            auto jindex2 = (i + 2)*x[ix].strides[0] + j*x[ix].strides[1];
            auto kindex2 = (i + 2)*x[ix].strides[0] + k*x[ix].strides[1];
            auto jindex3 = (i + 3)*x[ix].strides[0] + j*x[ix].strides[1];
            auto kindex3 = (i + 3)*x[ix].strides[0] + k*x[ix].strides[1];
            cov_jk[0] += (x_data[ix][jindex0] - j_mean)*(x_data[ix][kindex0] - k_mean);
            cov_jk[1] += (x_data[ix][jindex1] - j_mean)*(x_data[ix][kindex1] - k_mean);
            cov_jk[2] += (x_data[ix][jindex2] - j_mean)*(x_data[ix][kindex2] - k_mean);
            cov_jk[3] += (x_data[ix][jindex3] - j_mean)*(x_data[ix][kindex3] - k_mean);
        }
        for( ; i < num_rows; i++ ) {
            auto jindex = i*x[ix].strides[0] + j*x[ix].strides[1];
            auto kindex = i*x[ix].strides[0] + k*x[ix].strides[1];
            cov_jk[0] += (x_data[ix][jindex] - j_mean)*(x_data[ix][kindex] - k_mean);
        }
    }
    return (cov_jk[0] + cov_jk[1] + cov_jk[2] + cov_jk[3])/(total_values - 1);
}

template <typename scalar_t>
ALWAYS_INLINE scalar_t _compute_covariance(
    ArrayDescriptor x[], const scalar_t* x_data[],
    intptr_t nx, intptr_t j) {
    scalar_t j_mean = _compute_mean(x, x_data, nx, j);
    scalar_t cov_jj[4] = {0, 0, 0, 0};
    intptr_t total_values = 0;

    for( intptr_t ix = 0; ix < nx; ix++ ) {
        const intptr_t num_rows = x[ix].shape[0];
        total_values += num_rows;
        intptr_t i;
        for( i = 0; i + 3 < num_rows; i += 4 ) {
            auto jindex0 = i*x[ix].strides[0] + j*x[ix].strides[1];
            auto jindex1 = (i + 1)*x[ix].strides[0] + j*x[ix].strides[1];
            auto jindex2 = (i + 2)*x[ix].strides[0] + j*x[ix].strides[1];
            auto jindex3 = (i + 3)*x[ix].strides[0] + j*x[ix].strides[1];
            cov_jj[0] += (x_data[ix][jindex0] - j_mean)*(x_data[ix][jindex0] - j_mean);
            cov_jj[1] += (x_data[ix][jindex1] - j_mean)*(x_data[ix][jindex1] - j_mean);
            cov_jj[2] += (x_data[ix][jindex2] - j_mean)*(x_data[ix][jindex2] - j_mean);
            cov_jj[3] += (x_data[ix][jindex3] - j_mean)*(x_data[ix][jindex3] - j_mean);
        }
        for( ; i < num_rows; i++ ) {
            auto jindex = i*x[ix].strides[0] + j*x[ix].strides[1];
            cov_jj[0] += (x_data[ix][jindex] - j_mean)*(x_data[ix][jindex] - j_mean);
        }
    }
    return (cov_jj[0] + cov_jj[1] + cov_jj[2] + cov_jj[3])/(total_values - 1);
}

template <typename scalar_t>
ALWAYS_INLINE void _compute_covariance_matrix(ArrayDescriptor x[], const scalar_t* x_data[],
    intptr_t nx, ArrayDescriptor vi, scalar_t* vi_data) {
    const intptr_t num_cols = x[0].shape[1];
    for( intptr_t j = 0; j < num_cols; j++ ) {
        for( intptr_t k = j + 1; k < num_cols; k++ ) {
            scalar_t vi_ik_kj = _compute_covariance(x, x_data, nx, j, k);
            vi_data[j*vi.strides[0] + k*vi.strides[1]] = vi_ik_kj;
            vi_data[k*vi.strides[0] + j*vi.strides[1]] = vi_ik_kj;
        }
    }

    for( intptr_t j = 0; j < num_cols; j++ ) {
        scalar_t vi_ik_kj = _compute_covariance(x, x_data, nx, j);
        vi_data[j*vi.strides[0] + j*vi.strides[1]] = vi_ik_kj;
    }
}

template <typename scalar_t>
py::array cdist_mahalanobis_without_vi(
    const py::array& out_obj, const py::array& x_obj, const py::array& y_obj) {
    auto x_array = npy_asarray<scalar_t>(x_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto y_array = npy_asarray<scalar_t>(y_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out_array = py::cast<py::array_t<scalar_t>>(out_obj);
    auto x = get_descriptor(x_array);
    auto y = get_descriptor(y_array);
    auto out = get_descriptor(out_array);
    const scalar_t* x_data = static_cast<const scalar_t*>(x_array.data());
    const scalar_t* y_data = static_cast<const scalar_t*>(y_array.data());
    scalar_t* out_data = static_cast<scalar_t*>(out_array.mutable_data());

    const intptr_t num_cols = x.shape[1];

    ArrayDescriptor vi(2);
    vi.element_size = x.element_size;
    vi.shape = {num_cols, 2*num_cols};
    vi.strides = {2*num_cols, 1};
    scalar_t* vi_data = new scalar_t[2*num_cols*num_cols];

    ArrayDescriptor X[2] = {x, y};
    const scalar_t* X_data[2] = {x_data, y_data};

    _compute_covariance_matrix(X, X_data, 2, vi, vi_data);
    _compute_inverse_matrix(vi, vi_data);

    cdist_mahalanobis_impl(x, x_data, y, y_data, vi, vi_data, out, out_data, num_cols);

    delete [] vi_data;

    return out_array;
}

template <typename scalar_t>
py::array cdist_mahalanobis_with_vi(const py::array& out_obj, const py::array& x_obj,
                                    const py::array& y_obj, const py::array& vi_obj) {
    auto x_array = npy_asarray<scalar_t>(x_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto y_array = npy_asarray<scalar_t>(y_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto vi_array = npy_asarray<scalar_t>(vi_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out_array = py::cast<py::array_t<scalar_t>>(out_obj);
    auto x = get_descriptor(x_array);
    auto y = get_descriptor(y_array);
    auto vi = get_descriptor(vi_array);
    auto out = get_descriptor(out_array);
    const scalar_t* x_data = static_cast<const scalar_t*>(x_array.data());
    const scalar_t* y_data = static_cast<const scalar_t*>(y_array.data());
    const scalar_t* vi_data = static_cast<const scalar_t*>(vi_array.data());
    scalar_t* out_data = static_cast<scalar_t*>(out_array.mutable_data());

    cdist_mahalanobis_impl(x, x_data, y, y_data, vi, vi_data, out, out_data, 0);

    return out_array;
}

py::array cdist_mahalanobis(py::object& out_obj, const py::object& x_obj,
                const py::object& y_obj, const py::object& vi_obj) {
    auto x = npy_asarray(x_obj);
    auto y = npy_asarray(y_obj);
    auto vi = npy_asarray(vi_obj);
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

    if (vi_obj.is_none()) {
        auto m = x.shape(0) + y.shape(0);
        auto n = x.shape(1);
        if( m <= n ) {
            throw std::invalid_argument(
                std::string("The number of observations (") +
                std::to_string(m) + std::string(") is too ") +
                "small; the covariance matrix is " +
                "singular. For observations with " + std::to_string(n) +
                " dimensions, at least " + std::to_string(n + 1) +
                std::string(" observations are required."));
        }
        auto dtype = promote_type_real(common_type(x.dtype(), y.dtype()));
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            cdist_mahalanobis_without_vi<scalar_t>(out, x, y);
        });
        return out;
    }

    if( m != vi.shape(0) || m != vi.shape(1) ) {
        throw std::invalid_argument(
            "Inverse covariance matrix must have the same "
            "number of columns as XA and XB "
            "(i.e. feature dimension).");
    }

    auto dtype = promote_type_real(
        common_type(x.dtype(), y.dtype(), vi.dtype()));
    auto out = prepare_out_argument(out_obj, dtype, out_shape);
    DISPATCH_DTYPE(dtype, [&]{
        cdist_mahalanobis_with_vi<scalar_t>(out, x, y, vi);
    });
    return out;
}

template <typename scalar_t>
ALWAYS_INLINE void pdist_mahalanobis_impl(
    ArrayDescriptor x, const scalar_t* x_data,
    ArrayDescriptor vi, const scalar_t* vi_data,
    ArrayDescriptor out, scalar_t* out_data, intptr_t vi_col_offset) {
    const intptr_t num_rowsX = x.shape[0];
    const intptr_t num_cols = x.shape[1];
    scalar_t* uv_diff1 = new scalar_t[2*num_cols];
    scalar_t* uv_diff2 = uv_diff1 + num_cols;
    intptr_t o = 0;

    for( intptr_t i = 0; i < num_rowsX; i++ ) {
        const scalar_t* u = &x_data[i * x.strides[0]];
        for( intptr_t j = i + 1; j < num_rowsX; j++, o++ ) {
            const scalar_t* v = &x_data[j * x.strides[0]];

            for( intptr_t k = 0; k < num_cols; k++ ) {
                uv_diff1[k] = u[k*x.strides[1]] - v[k*x.strides[1]];
            }

            for( intptr_t k = 0; k < num_cols; k++ ) {
                uv_diff2[k] = 0;
                for( intptr_t l = 0; l < num_cols; l++ ) {
                    uv_diff2[k] += uv_diff1[l] * vi_data[k*vi.strides[0] + (l + vi_col_offset)*vi.strides[1]];
                }
            }

            scalar_t s = 0;
            for( intptr_t k = 0; k < num_cols; k++ ) {
                s += uv_diff1[k] * uv_diff2[k];
            }

            out_data[o] = std::sqrt(s);
        }
    }

    delete [] uv_diff1;
}

template <typename scalar_t>
py::array pdist_mahalanobis_without_vi(
    const py::array& out_obj, const py::array& x_obj) {
    auto x_array = npy_asarray<scalar_t>(x_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out_array = py::cast<py::array_t<scalar_t>>(out_obj);
    auto x = get_descriptor(x_array);
    auto out = get_descriptor(out_array);
    const scalar_t* x_data = static_cast<const scalar_t*>(x_array.data());
    scalar_t* out_data = static_cast<scalar_t*>(out_array.mutable_data());

    const intptr_t num_cols = x.shape[1];

    ArrayDescriptor vi(2);
    vi.element_size = x.element_size;
    vi.shape = {num_cols, 2*num_cols};
    vi.strides = {2*num_cols, 1};
    scalar_t* vi_data = new scalar_t[2*num_cols*num_cols];

    ArrayDescriptor X[1] = {x};
    const scalar_t* X_data[1] = {x_data};

    _compute_covariance_matrix(X, X_data, 1, vi, vi_data);
    _compute_inverse_matrix(vi, vi_data);

    pdist_mahalanobis_impl(x, x_data, vi, vi_data, out, out_data, num_cols);

    delete [] vi_data;

    return out_array;
}

template <typename scalar_t>
py::array pdist_mahalanobis_with_vi(
    const py::array& out_obj, const py::array& x_obj, const py::array& vi_obj) {
    auto x_array = npy_asarray<scalar_t>(x_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto vi_array = npy_asarray<scalar_t>(vi_obj,
                                 NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED);
    auto out_array = py::cast<py::array_t<scalar_t>>(out_obj);
    auto x = get_descriptor(x_array);
    auto vi = get_descriptor(vi_array);
    auto out = get_descriptor(out_array);
    const scalar_t* x_data = static_cast<const scalar_t*>(x_array.data());
    const scalar_t* vi_data = static_cast<const scalar_t*>(vi_array.data());
    scalar_t* out_data = static_cast<scalar_t*>(out_array.mutable_data());

    pdist_mahalanobis_impl(x, x_data, vi, vi_data, out, out_data, 0);

    return out_array;
}

py::array pdist_mahalanobis(const py::object& out_obj,
    const py::object& x_obj, const py::object& vi_obj) {
    auto x = npy_asarray(x_obj);
    auto vi = npy_asarray(vi_obj);
    if (x.ndim() != 2) {
        throw std::invalid_argument("x must be 2-dimensional");
    }

    const intptr_t n = x.shape(1);
    const intptr_t m = x.shape(0);
    std::array<intptr_t, 1> out_shape{{(m * (m - 1)) / 2}};

    if (vi_obj.is_none()) {
        if( m <= n ) {
            throw std::invalid_argument(
                std::string("The number of observations (") +
                std::to_string(m) + std::string(") is too ") +
                "small; the covariance matrix is " +
                "singular. For observations with " + std::to_string(n) +
                " dimensions, at least " + std::to_string(n + 1) +
                std::string(" observations are required."));
        }
        auto dtype = promote_type_real(x.dtype());
        auto out = prepare_out_argument(out_obj, dtype, out_shape);
        DISPATCH_DTYPE(dtype, [&]{
            pdist_mahalanobis_without_vi<scalar_t>(out, x);
        });
        return out;
    }

    if( n != vi.shape(0) || n != vi.shape(1) ) {
        throw std::invalid_argument(
            "Inverse covariance matrix must have the same "
            "number of columns as X (i.e. feature dimension).");
    }

    auto dtype = promote_type_real(common_type(x.dtype(), vi.dtype()));
    auto out = prepare_out_argument(out_obj, dtype, out_shape);
    DISPATCH_DTYPE(dtype, [&]{
        pdist_mahalanobis_with_vi<scalar_t>(out, x, vi);
    });
    return out;
}

PYBIND11_MODULE(_distance_pybind, m, py::mod_gil_not_used()) {
    if (_import_array() != 0) {
        throw py::error_already_set();
    }
    using namespace pybind11::literals;
    m.def("pdist_canberra",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, CanberraDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_hamming",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, HammingDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_dice",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, DiceDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_jaccard",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, JaccardDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_kulczynski1",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, Kulczynski1Distance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_rogerstanimoto",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, RogerstanimotoDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_russellrao",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, RussellRaoDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_sokalmichener",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, SokalmichenerDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_sokalsneath",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, SokalsneathDistance{});
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_yule",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, YuleDistance{});
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
    m.def("pdist_seuclidean",
          [](py::object x, py::object w, py::object out) {
              if( w.is_none() ) {
                return pdist(out, x, w, EuclideanDistance{},
                             PreprocessingType::ComputeVariance);
              } else {
                return pdist(out, x, w, EuclideanDistance{});
              }
          },
          "x"_a, "w"_a, "out"_a=py::none());
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
    m.def("pdist_mahalanobis",
          [](py::object x, py::object vi, py::object out) {
              return pdist_mahalanobis(out, x, vi);
          },
          "x"_a, "vi"_a=py::none(), "out"_a=py::none());
    m.def("pdist_cosine",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, CosineDistance{}, PreprocessingType::Normalise);
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_correlation",
          [](py::object x, py::object w, py::object out) {
              return pdist(out, x, w, CosineDistance{},
                           PreprocessingType::CentraliseAndNormalise);
          },
          "x"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("pdist_jensenshannon",
          [](py::object x, py::object out) {
              return pdist(out, x, py::none(), JensenshannonDistance{},
                           PreprocessingType::ScaleInputForJensenshannon_);
          },
          "x"_a, "out"_a=py::none());
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
    m.def("cdist_dice",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, DiceDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_jaccard",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, JaccardDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_kulczynski1",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, Kulczynski1Distance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_hamming",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, HammingDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_rogerstanimoto",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, RogerstanimotoDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_russellrao",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, RussellRaoDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_sokalmichener",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, SokalmichenerDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_sokalsneath",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, SokalsneathDistance{});
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_yule",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, YuleDistance{});
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
    m.def("cdist_seuclidean",
          [](py::object x, py::object y, py::object w, py::object out) {
              if( w.is_none() ) {
                  return cdist(out, x, y, w, EuclideanDistance{},
                               PreprocessingType::ComputeVariance);
              } else {
                  return cdist(out, x, y, w, EuclideanDistance{});
              }
          },
          "x"_a, "y"_a, "w"_a, "out"_a=py::none());
    m.def("cdist_mahalanobis",
          [](py::object x, py::object y, py::object vi, py::object out) {
              return cdist_mahalanobis(out, x, y, vi);
          },
          "x"_a, "y"_a, "vi"_a=py::none(), "out"_a=py::none());
    m.def("cdist_cosine",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, CosineDistance{},
                           PreprocessingType::Normalise);
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_correlation",
          [](py::object x, py::object y, py::object w, py::object out) {
              return cdist(out, x, y, w, CosineDistance{},
                           PreprocessingType::CentraliseAndNormalise);
          },
          "x"_a, "y"_a, "w"_a=py::none(), "out"_a=py::none());
    m.def("cdist_jensenshannon",
          [](py::object x, py::object y, py::object out) {
              return cdist(out, x, y, py::none(), JensenshannonDistance{},
                           PreprocessingType::ScaleInputForJensenshannon_);
          },
          "x"_a, "y"_a, "out"_a=py::none());
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
