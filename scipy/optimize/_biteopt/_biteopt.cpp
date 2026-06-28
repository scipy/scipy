#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <exception>
#include <limits>
#include <cmath>

#include <numpy/random/bitgen.h>

#include "biteopt.h"

namespace py = pybind11;

namespace {

// Adapter matching biteopt's ``biteopt_rng`` typedef
// (``uint32_t (*)(void* rng_data)``). biteopt's CBiteRnd calls this twice per
// 64-bit draw, combining two 32-bit outputs, so returning a single uint32 per
// call is exactly what is expected. ``rng_data`` is a NumPy ``bitgen_t*``.
uint32_t numpy_rng_adapter(void* rng_data) {
    auto* bitgen = static_cast<bitgen_t*>(rng_data);
    return bitgen->next_uint32(bitgen->state);
}

// Context threaded through biteopt as the objective's ``void* data``.
struct CallbackContext {
    py::object func;        // the user-supplied Python objective
    std::exception_ptr eptr;  // first exception raised by the objective, if any
};

// Adapter matching ``biteopt_func``. It is the single bridge between biteopt's
// C++ optimization loop and the Python objective. Crucially, it never lets an
// exception unwind into biteopt (which is not exception-safe): any exception is
// captured and re-raised after biteopt has returned cleanly. While unwinding is
// pending it returns NaN, which biteopt sanitizes to a large penalty cost
// internally (fixCostNaN), so the affected points are simply never selected.
double trampoline(int N, const double* x, void* data) {
    auto* ctx = static_cast<CallbackContext*>(data);

    // Short-circuit once an error has been seen: do not re-enter Python.
    // biteopt has no clean mid-run abort, so it keeps iterating to the end of
    // its budget; every remaining call lands here and returns NaN cheaply
    // (no Python), and the captured exception is re-raised once biteopt_minimize
    // returns. For very large maxfun this delays the raise, but the wasted
    // iterations do no real work.
    if (ctx->eptr) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    try {
        // Acquire the GIL for the duration of the Python call
        py::gil_scoped_acquire gil;

        // Hand the objective a fresh copy of the parameter vector so it cannot
        // corrupt biteopt's internal buffers. The cost is passed straight
        // through; biteopt sanitizes any NaN itself.
        py::array_t<double> arr(static_cast<size_t>(N), x);
        return ctx->func(arr).cast<double>();
    } catch (...) {
        // Capture the live Python exception (or any C++ exception, e.g. a
        // failed cast) and stop feeding biteopt real values.
        ctx->eptr = std::current_exception();
        return std::numeric_limits<double>::quiet_NaN();
    }
}

py::object minimize(
    py::object func,
    py::array_t<double, py::array::c_style | py::array::forcecast> lb,
    py::array_t<double, py::array::c_style | py::array::forcecast> ub,
    int iter,
    int depth,
    int attc,
    py::object bit_generator,
    py::object f_min
) {
    const int N = static_cast<int>(lb.shape(0));

    // Obtain NumPy's bitgen_t* from the BitGenerator capsule.
    auto* bitgen = static_cast<bitgen_t*>(
        PyCapsule_GetPointer(bit_generator.ptr(), "BitGenerator")
    );
    if (bitgen == nullptr) {
        throw py::error_already_set();
    }

    // Optional early-stopping threshold: biteopt stops once the best objective
    // value reaches f_min. A null pointer disables the criterion.
    double f_min_value = f_min.cast<double>();
    double* f_minp = nullptr;
    if (std::isfinite(f_min_value)) {
        f_minp = &f_min_value;
    }

    CallbackContext ctx;
    ctx.func = func;

    // Let biteopt write the best parameters straight into the output array's
    // buffer, avoiding an intermediate vector and a copy.
    py::array_t<double> x_out(static_cast<size_t>(N));
    double best_f = 0.0;

    // biteopt_minimize sets up CBiteOptDeep internally, seeds its PRNG from our
    // NumPy adapter (rf/rdata), runs `attc` attempts of `iter` iterations each,
    // and returns the total number of objective evaluations performed. The
    // convergence criterion (stopc) is left disabled: on smooth problems
    // biteopt keeps making small improving steps that reset its plateau
    // counter, so it rarely fires and is not a reliable success signal.
    // The GIL is intentionally *not* released here. biteopt's own per-iteration
    // C++ work is only ~100 ns, so overlapping it across threads yields no
    // measurable throughput on a GIL build (and adds per-callback acquire/release
    // overhead). Real parallelism comes from a free-threaded interpreter or from
    // objectives that release the GIL themselves; neither needs a release here.
    const int evals = biteopt_minimize(
        N, &trampoline, &ctx, lb.data(), ub.data(),
        x_out.mutable_data(), &best_f, iter, depth, attc,
        /*stopc=*/0, &numpy_rng_adapter, bitgen, f_minp
    );

    // biteopt has returned and its objects have destructed normally; now it is
    // safe to re-raise. pybind11 restores the original Python exception type
    // and traceback for an ``error_already_set``.
    if (ctx.eptr) {
        std::rethrow_exception(ctx.eptr);
    }

    py::dict result;
    result["x"] = x_out;
    result["fun"] = best_f;
    result["nfev"] = evals;
    return result;
}

}  // namespace

PYBIND11_MODULE(_biteopt, m, py::mod_gil_not_used()) {
    m.def(
        "minimize", &minimize,
        py::arg("func"), py::arg("lb"), py::arg("ub"),
        py::arg("iter"), py::arg("depth"), py::arg("attc"),
        py::arg("bit_generator"), py::arg("f_min")
    );
    // GIL safety: the trampoline acquires the GIL (py::gil_scoped_acquire)
    // before every Python call, so this module is safe to load without the
    // GIL. The minimize() entry point is always called from Python and
    // therefore has the GIL held by pybind11's normal dispatch mechanism.
}
