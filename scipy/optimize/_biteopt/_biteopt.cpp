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
//
// Callback design: the optional ``callback`` is invoked once per objective
// *evaluation* (from the trampoline), passing the point that was just
// evaluated. This deliberately does NOT implement SciPy's newer
// ``callback(intermediate_result=...)`` convention: that interface fires once
// per *iteration* with the best-so-far incumbent, whereas biteopt exposes no
// per-iteration hook and the point handed to the trampoline is a trial
// candidate that is usually not the incumbent. Feeding those trial points into
// an ``intermediate_result`` would report a non-monotonic ``fun`` and an ``x``
// that is not "best so far", so we keep the simpler ``callback(x)``
// contract (matching ``scipy.optimize.direct``).
struct CallbackContext {
    py::object func;        // the user-supplied Python objective
    py::object callback;    // optional user callback called per evaluation
    bool callback_stopped = false;  // callback raised StopIteration
    // Mutable stop threshold read by biteopt each iteration (passed as f_minp).
    // NaN disables it (all comparisons with NaN are false); a finite value
    // enables early stop at that cost; the trampoline flips it to +inf to
    // request a graceful stop after a StopIteration from the callback.
    double f_stop_value = std::numeric_limits<double>::quiet_NaN();
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
    if (ctx->eptr || ctx->callback_stopped) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    try {
        // minimize() holds the GIL for the whole run, so it is already held
        // here -- no acquire needed. (Re-acquire it here if a gil_scoped_release
        // is ever added around biteopt_minimize.)
        //
        // Hand the objective a fresh copy of the parameter vector so it cannot
        // corrupt biteopt's internal buffers. The cost is passed straight
        // through; biteopt sanitizes any NaN itself.
        py::array_t<double> arr(static_cast<size_t>(N), x);
        const double fx = ctx->func(arr).cast<double>();

        // Invoke the optional user callback with the just-evaluated point.
        // StopIteration is the sanctioned "please stop" signal across SciPy
        // optimizers, so we treat it specially: instead of letting it propagate
        // as an error (which the outer catch would capture and re-raise), we
        // record the request and arm the stop threshold so biteopt winds down
        // gracefully. This mirrors SciPy's ``_call_callback_maybe_halt``, which
        // converts a callback StopIteration into a clean early termination
        // (surfaced as ``success=False`` on the Python side).
        // Any OTHER exception from the callback is a genuine error: it falls
        // through to the outer catch and is re-raised verbatim, exactly like an
        // exception raised by the objective itself.
        if (!ctx->callback.is_none()) {
            try {
                ctx->callback(arr);
            } catch (py::error_already_set& e) {
                if (e.matches(PyExc_StopIteration)) {
                    // Graceful stop: remember it and set the early-stop
                    // threshold to +inf so biteopt's next
                    // ``getBestCost() <= f_minp`` check trips and the run exits
                    // without unwinding an exception through biteopt's
                    // (non-exception-safe) core. Later trampoline entries
                    // short-circuit at the top of this function.
                    ctx->callback_stopped = true;
                    ctx->f_stop_value = std::numeric_limits<double>::infinity();
                } else {
                    throw;
                }
            }
        }

        return fx;
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
    py::object f_min,
    py::object callback
) {
    const int N = static_cast<int>(lb.shape(0));

    // Obtain NumPy's bitgen_t* from the BitGenerator capsule.
    auto* bitgen = static_cast<bitgen_t*>(
        PyCapsule_GetPointer(bit_generator.ptr(), "BitGenerator")
    );
    if (bitgen == nullptr) {
        throw py::error_already_set();
    }

    CallbackContext ctx;
    ctx.func = func;
    ctx.callback = callback;
    // A finite f_min enables biteopt's early-stop threshold; otherwise it stays
    // NaN (disabled). If callback(x) raises StopIteration the trampoline flips
    // ctx.f_stop_value to +inf, so biteopt stops at its next threshold check
    // without unwinding exceptions through biteopt internals.
    double f_min_value = f_min.cast<double>();
    if (std::isfinite(f_min_value)) {
        ctx.f_stop_value = f_min_value;
    }

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
        /*stopc=*/0, &numpy_rng_adapter, bitgen, &ctx.f_stop_value
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
    result["callback_stopped"] = ctx.callback_stopped;
    return result;
}

}  // namespace

PYBIND11_MODULE(_biteopt, m, py::mod_gil_not_used()) {
    m.def(
        "minimize", &minimize,
        py::arg("func"), py::arg("lb"), py::arg("ub"),
        py::arg("iter"), py::arg("depth"), py::arg("attc"),
        py::arg("bit_generator"), py::arg("f_min"), py::arg("callback")
    );
    // Free-threading: minimize() uses only per-call state and holds the bit
    // generator's lock around the run (see _biteopt_py.py), so it carries no
    // GIL-dependent shared state and is safe under py::mod_gil_not_used(). It is
    // always entered from Python with an attached thread state that is held for
    // the whole run (the GIL is never released), so the objective callback can
    // call into Python directly.
}
