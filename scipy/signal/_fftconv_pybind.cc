#include <string_view>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "_fftconv.h"

namespace py = pybind11;

static fftconv::ConvMode parse_mode(std::string_view m) {
    if (m == "full")  return fftconv::ConvMode::Full;
    if (m == "same")  return fftconv::ConvMode::Same;
    if (m == "valid") return fftconv::ConvMode::Valid;
    throw std::invalid_argument("acceptable mode flags are 'valid', 'same', or 'full'");
}

// Parse the axes argument to a single non-negative conv axis, or return -1 to
// decline. None -> axis 0 only when 1-D (None on N-D means all axes = multi-axis).
static long parse_single_axis(py::object axes, py::ssize_t ndim) {
    if (axes.is_none()) return (ndim == 1) ? 0 : -1;
    auto norm = [ndim](long a) -> long { if (a < 0) a += ndim; return (a >= 0 && a < ndim) ? a : -1; };
    if (py::isinstance<py::int_>(axes)) return norm(axes.cast<long>());
    try {
        auto v = axes.cast<std::vector<long>>();
        if (v.size() != 1) return -1;
        return norm(v[0]);
    } catch (...) { return -1; }
}

// Non-contiguous inputs of the accepted dtype are handled by an internal
// contiguous copy (the `c_style | forcecast` conversion below), not declined.
// dtype/ndim/axes are already validated in `oaconvolve` before this point, so
// forcecast here only normalizes memory layout; it never silently changes dtype.
template<typename T>
static py::object run_nd(py::array_t<T, py::array::c_style | py::array::forcecast> a,
                         py::array_t<T, py::array::c_style | py::array::forcecast> b,
                         long axis, fftconv::ConvMode mode) {
    const py::ssize_t ndim = a.ndim();
    size_t outer = 1, inner = 1;
    for (py::ssize_t d = 0; d < axis; ++d)        outer *= static_cast<size_t>(a.shape(d));
    for (py::ssize_t d = axis + 1; d < ndim; ++d) inner *= static_cast<size_t>(a.shape(d));
    const size_t n = static_cast<size_t>(a.shape(axis));
    const size_t m = static_cast<size_t>(b.shape(axis));
    size_t out_n;
    switch (mode) {
        case fftconv::ConvMode::Full:  out_n = n + m - 1; break;
        case fftconv::ConvMode::Same:  out_n = n; break;
        default: { size_t big = std::max(n, m), sm = std::min(n, m); out_n = big - sm + 1; }
    }
    std::vector<py::ssize_t> oshape(a.shape(), a.shape() + ndim);
    oshape[axis] = static_cast<py::ssize_t>(out_n);
    py::array_t<T> out(oshape);
    bool ok;
    {
        py::gil_scoped_release rel;
        if (outer == 1 && inner == 1)
            ok = fftconv::oaconvolve_1d<T>(a.data(), n, b.data(), m,
                                           out.mutable_data(), mode, 1);
        else
            ok = fftconv::oaconvolve_axis<T>(a.data(), outer, n, inner, b.data(), m,
                                             out.mutable_data(), mode, 1);
    }
    if (!ok) return py::none();  // decline -> Python fallback
    return out;
}

static py::object oaconvolve(py::object in1, py::object in2,
                             std::string_view mode_str, py::object axes) {
    auto a = py::array::ensure(in1);
    auto b = py::array::ensure(in2);
    if (!a || !b) return py::none();
    if (a.ndim() != b.ndim()) return py::none();
    const long axis = parse_single_axis(axes, a.ndim());
    if (axis < 0) return py::none();  // not a single-axis request -> decline
    // Non-conv dims must match between the inputs (broadcasting -> decline).
    for (py::ssize_t d = 0; d < a.ndim(); ++d)
        if (d != axis && a.shape(d) != b.shape(d)) return py::none();
    const auto mode = parse_mode(mode_str);
    // Decline the block-size fallback here, before run_nd's forcecast
    // conversion (which contiguous-copies non-contiguous inputs) and output
    // allocation can do speculative work that a decline would throw away.
    // The same check inside oaconvolve_1d/oaconvolve_axis is then unreachable
    // from Python but kept so the C++ core stands on its own.
    if (fftconv::calc_oa_lens(static_cast<size_t>(a.shape(axis)),
                              static_cast<size_t>(b.shape(axis))).fallback)
        return py::none();

    if (py::isinstance<py::array_t<double>>(a) && py::isinstance<py::array_t<double>>(b))
        return run_nd<double>(in1, in2, axis, mode);
    if (py::isinstance<py::array_t<float>>(a) && py::isinstance<py::array_t<float>>(b))
        return run_nd<float>(in1, in2, axis, mode);
    if (py::isinstance<py::array_t<std::complex<double>>>(a)
            && py::isinstance<py::array_t<std::complex<double>>>(b))
        return run_nd<std::complex<double>>(in1, in2, axis, mode);
    if (py::isinstance<py::array_t<std::complex<float>>>(a)
            && py::isinstance<py::array_t<std::complex<float>>>(b))
        return run_nd<std::complex<float>>(in1, in2, axis, mode);
    return py::none();  // mixed / other dtype -> Python fallback
}

// Mirrors _signaltools._calc_oa_lens: (block_size, overlap, in1_step,
// in2_step), with overlap=None signaling the single-FFT fallback.
static py::tuple calc_oa_lens(size_t s1, size_t s2) {
    fftconv::OaLens r = fftconv::calc_oa_lens(s1, s2);
    py::object overlap = r.fallback ? py::none() : py::cast(r.overlap);
    return py::make_tuple(r.block_size, overlap, r.in1_step, r.in2_step);
}

PYBIND11_MODULE(_fftconv, m) {
    m.doc() = "C++ overlap-add convolution (duccfft backend)";
    m.def("oaconvolve", &oaconvolve,
          py::arg("in1"), py::arg("in2"), py::arg("mode"), py::arg("axes"));
    // Internal helpers exposed for testing only.
    m.def("_lambertw_wm1", &fftconv::lambertw_wm1, py::arg("x"));
    m.def("_calc_oa_lens", &calc_oa_lens, py::arg("s1"), py::arg("s2"));
}
