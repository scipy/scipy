/*
 * Compiled routines for the NumPy backend of `scipy.interpolate.RBFInterpolator`.
 *
 * Every routine that evaluates an RBF takes a `kernel` argument that is either
 * the name of a built-in RBF (str) or a PyCapsule named "double (double)" that
 * holds a user-supplied kernel (the capsule of a `LowLevelCallable`). The
 * built-in kernels are template arguments of the loops so that the compiler can
 * inline them; only user kernels are called through a function pointer.
 *
 * The loops run without the GIL. A user kernel is called without the GIL, too;
 * a ctypes callback that wraps a Python function takes the GIL itself.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <vector>

#include "numpy/arrayobject.h"


extern "C" {
typedef double rbf_kernel_func(double);
}


namespace {

const char KERNEL_SIGNATURE[] = "double (double)";

constexpr double SQRT3 = 1.7320508075688772935274463415059;
constexpr double SQRT5 = 2.2360679774997896964091736687313;


/*
 * Built-in kernels, as functions of the scaled distance r = epsilon*|x - y|.
 * Keep the order of `BUILTIN_KERNELS` below in sync with `KernelId`.
 */
struct Linear {
    static double eval(double r) noexcept { return -r; }
};
struct ThinPlateSpline {
    static double eval(double r) noexcept { return r == 0.0 ? 0.0 : r * r * std::log(r); }
};
struct Cubic {
    static double eval(double r) noexcept { return r * r * r; }
};
struct Quintic {
    static double eval(double r) noexcept { return -(r * r) * (r * r) * r; }
};
struct Multiquadric {
    static double eval(double r) noexcept { return -std::sqrt(r * r + 1.0); }
};
struct InverseMultiquadric {
    static double eval(double r) noexcept { return 1.0 / std::sqrt(r * r + 1.0); }
};
struct InverseQuadratic {
    static double eval(double r) noexcept { return 1.0 / (r * r + 1.0); }
};
struct Gaussian {
    static double eval(double r) noexcept { return std::exp(-(r * r)); }
};
struct Matern1_2 {
    static double eval(double r) noexcept { return std::exp(-r); }
};
struct Matern3_2 {
    static double eval(double r) noexcept {
        const double term = SQRT3 * r;
        return (1.0 + term) * std::exp(-term);
    }
};
struct Matern5_2 {
    static double eval(double r) noexcept {
        const double term = SQRT5 * r;
        return (1.0 + term + 5.0 * (r * r) / 3.0) * std::exp(-term);
    }
};


// Plain function with C calling convention for a built-in kernel, for `_kernel_capsule`.
template <class K>
double kernel_as_function(double r) noexcept { return K::eval(r); }


enum KernelId {
    LINEAR,
    THIN_PLATE_SPLINE,
    CUBIC,
    QUINTIC,
    MULTIQUADRIC,
    INVERSE_MULTIQUADRIC,
    INVERSE_QUADRATIC,
    GAUSSIAN,
    MATERN1_2,
    MATERN3_2,
    MATERN5_2,
    NUM_BUILTIN_KERNELS,
    USER_KERNEL = NUM_BUILTIN_KERNELS
};


struct BuiltinKernel {
    const char *name;
    rbf_kernel_func *func;
};


const BuiltinKernel BUILTIN_KERNELS[NUM_BUILTIN_KERNELS] = {
    {"linear", &kernel_as_function<Linear>},
    {"thin_plate_spline", &kernel_as_function<ThinPlateSpline>},
    {"cubic", &kernel_as_function<Cubic>},
    {"quintic", &kernel_as_function<Quintic>},
    {"multiquadric", &kernel_as_function<Multiquadric>},
    {"inverse_multiquadric", &kernel_as_function<InverseMultiquadric>},
    {"inverse_quadratic", &kernel_as_function<InverseQuadratic>},
    {"gaussian", &kernel_as_function<Gaussian>},
    {"matern1_2", &kernel_as_function<Matern1_2>},
    {"matern3_2", &kernel_as_function<Matern3_2>},
    {"matern5_2", &kernel_as_function<Matern5_2>},
};


// The `kernel` argument after parsing; `func` is only used for USER_KERNEL.
struct Kernel {
    int id;
    rbf_kernel_func *func;
};


template <class K>
struct Inlined {
    double operator()(double r) const noexcept { return K::eval(r); }
};


struct ThroughPointer {
    rbf_kernel_func *func;
    double operator()(double r) const { return func(r); }
};


// Call `f` with a callable object for the kernel `k`.
template <class F>
void with_kernel(const Kernel &k, F &&f)
{
    switch (k.id) {
        case LINEAR:               f(Inlined<Linear>{});              break;
        case THIN_PLATE_SPLINE:    f(Inlined<ThinPlateSpline>{});     break;
        case CUBIC:                f(Inlined<Cubic>{});               break;
        case QUINTIC:              f(Inlined<Quintic>{});             break;
        case MULTIQUADRIC:         f(Inlined<Multiquadric>{});        break;
        case INVERSE_MULTIQUADRIC: f(Inlined<InverseMultiquadric>{}); break;
        case INVERSE_QUADRATIC:    f(Inlined<InverseQuadratic>{});    break;
        case GAUSSIAN:             f(Inlined<Gaussian>{});            break;
        case MATERN1_2:            f(Inlined<Matern1_2>{});           break;
        case MATERN3_2:            f(Inlined<Matern3_2>{});           break;
        case MATERN5_2:            f(Inlined<Matern5_2>{});           break;
        default:                   f(ThroughPointer{k.func});         break;
    }
}


// "O&" converter for the `kernel` argument.
int kernel_converter(PyObject *obj, void *out)
{
    Kernel *k = static_cast<Kernel *>(out);
    if (PyUnicode_Check(obj)) {
        const char *name = PyUnicode_AsUTF8(obj);
        if (name == NULL) { return 0; }
        for (int i = 0; i < NUM_BUILTIN_KERNELS; i++) {
            if (std::strcmp(name, BUILTIN_KERNELS[i].name) == 0) {
                k->id = i;
                k->func = BUILTIN_KERNELS[i].func;
                return 1;
            }
        }
        PyErr_Format(PyExc_ValueError, "Unknown RBF kernel '%s'.", name);
        return 0;
    }
    if (PyCapsule_CheckExact(obj)) {
        // Sets ValueError if the capsule has another signature.
        void *ptr = PyCapsule_GetPointer(obj, KERNEL_SIGNATURE);
        if (ptr == NULL) { return 0; }
        k->id = USER_KERNEL;
        k->func = reinterpret_cast<rbf_kernel_func *>(ptr);
        return 1;
    }
    PyErr_SetString(PyExc_TypeError,
                    "`kernel` must be the name of a built-in RBF or a PyCapsule "
                    "with signature \"double (double)\".");
    return 0;
}


struct ArrayDeleter {
    void operator()(PyArrayObject *p) const noexcept { Py_XDECREF(p); }
};


/*
 * Owning reference to a NumPy array; `release()` hands the reference over, e.g. to
 * the "N" format of Py_BuildValue.
 */
using OwnedArray = std::unique_ptr<PyArrayObject, ArrayDeleter>;


/*
 * `obj` as an aligned, C-contiguous, native-endian array of `typenum` with `ndim`
 * dimensions. Copies only when needed; casts only when the cast is safe.
 */
OwnedArray as_c_array(PyObject *obj, int typenum, int ndim, const char *name)
{
    OwnedArray arr(reinterpret_cast<PyArrayObject *>(
        PyArray_FROMANY(obj, typenum, 0, 0, NPY_ARRAY_IN_ARRAY)));
    if (arr && PyArray_NDIM(arr.get()) != ndim) {
        PyErr_Format(PyExc_ValueError, "`%s` must be a %d-D array, got %d-D.",
                     name, ndim, PyArray_NDIM(arr.get()));
        arr.reset();
    }
    return arr;
}


inline const double *f64_data(const OwnedArray &a)
{
    return static_cast<const double *>(PyArray_DATA(a.get()));
}


inline double *f64_data_mut(const OwnedArray &a)
{
    return static_cast<double *>(PyArray_DATA(a.get()));
}


inline npy_intp dim(const OwnedArray &a, int i) { return PyArray_DIM(a.get(), i); }


// Releases the GIL for its lifetime; restores it during unwinding as well.
class GILRelease {
  public:
    GILRelease() noexcept : state_(PyEval_SaveThread()) {}
    ~GILRelease() { PyEval_RestoreThread(state_); }
    GILRelease(const GILRelease &) = delete;
    GILRelease &operator=(const GILRelease &) = delete;
  private:
    PyThreadState *state_;
};


// Euclidean distance between the `n`-vectors at `a` and `b`.
inline double distance(const double *a, const double *b, npy_intp n) noexcept
{
    double s = 0.0;
    for (npy_intp k = 0; k < n; k++) {
        const double t = a[k] - b[k];
        s += t * t;
    }
    return std::sqrt(s);
}


// The monomial with exponents `pw` evaluated at the `n`-vector at `x`.
inline double monomial(const double *x, const int64_t *pw, npy_intp n) noexcept
{
    double v = 1.0;
    for (npy_intp k = 0; k < n; k++) {
        v *= std::pow(x[k], static_cast<double>(pw[k]));
    }
    return v;
}


/*
 * out[i*rs + j*cs] = kern(|x_i - x_j|) for the `p` rows of the C-contiguous (p, n)
 * array `x`. The strides pick C (rs=p, cs=1) or Fortran (rs=1, cs=ld) order.
 */
template <class K>
void kernel_matrix(K kern, const double *x, npy_intp p, npy_intp n,
                   double *out, npy_intp rs, npy_intp cs)
{
    for (npy_intp i = 0; i < p; i++) {
        const double *xi = x + i * n;
        for (npy_intp j = 0; j <= i; j++) {
            const double v = kern(distance(xi, x + j * n, n));
            out[i * rs + j * cs] = v;
            out[j * rs + i * cs] = v;
        }
    }
}


// The (q, n) array `x` scaled by `epsilon`, C order.
std::vector<double> scaled(const double *x, npy_intp q, npy_intp n, double epsilon)
{
    std::vector<double> out(static_cast<size_t>(q * n));
    for (npy_intp i = 0; i < q * n; i++) { out[i] = x[i] * epsilon; }
    return out;
}


// The (q, n) array `x` shifted and scaled column-wise, C order.
std::vector<double> normalized(const double *x, npy_intp q, npy_intp n,
                               const double *shift, const double *scale)
{
    std::vector<double> out(static_cast<size_t>(q * n));
    for (npy_intp i = 0; i < q; i++) {
        for (npy_intp k = 0; k < n; k++) {
            out[i * n + k] = (x[i * n + k] - shift[k]) / scale[k];
        }
    }
    return out;
}


// Translate the C++ exception being handled into a Python exception; call from `catch`.
PyObject *set_cpp_error()
{
    try {
        throw;
    } catch (const std::bad_alloc &) {
        return PyErr_NoMemory();
    } catch (...) {
        PyErr_SetString(PyExc_SystemError, "Unexpected C++ exception in _rbfinterp_cpp.");
        return NULL;
    }
}


PyDoc_STRVAR(doc_kernel_matrix,
"_kernel_matrix(x, kernel)\n"
"\n"
"Return the RBFs with centers at `x`, evaluated at `x`.\n"
"\n"
"Parameters\n"
"----------\n"
"x : (P, N) float ndarray\n"
"    Points.\n"
"kernel : str or PyCapsule\n"
"    Name of a built-in RBF, or a capsule with signature \"double (double)\".\n"
"\n"
"Returns\n"
"-------\n"
"(P, P) float ndarray\n");


PyObject *py_kernel_matrix(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *x_obj;
    Kernel kernel = {LINEAR, NULL};
    if (!PyArg_ParseTuple(args, "OO&", &x_obj, kernel_converter, &kernel)) {
        return NULL;
    }
    OwnedArray x = as_c_array(x_obj, NPY_DOUBLE, 2, "x");
    if (!x) { return NULL; }
    const npy_intp p = dim(x, 0), n = dim(x, 1);

    npy_intp dims[2] = {p, p};
    OwnedArray out(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(2, dims, NPY_DOUBLE, 0)));
    if (!out) { return NULL; }

    const double *xd = f64_data(x);
    double *od = f64_data_mut(out);
    try {
        GILRelease nogil;
        with_kernel(kernel, [&](auto kern) { kernel_matrix(kern, xd, p, n, od, p, 1); });
    } catch (...) {
        return set_cpp_error();
    }
    return reinterpret_cast<PyObject *>(out.release());
}


PyDoc_STRVAR(doc_polynomial_matrix,
"_polynomial_matrix(x, powers)\n"
"\n"
"Return the monomials with exponents from `powers`, evaluated at `x`.\n"
"\n"
"Parameters\n"
"----------\n"
"x : (Q, N) float ndarray\n"
"    Points.\n"
"powers : (R, N) int ndarray\n"
"    The exponents for each monomial.\n"
"\n"
"Returns\n"
"-------\n"
"(Q, R) float ndarray\n");


PyObject *py_polynomial_matrix(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *x_obj, *powers_obj;
    if (!PyArg_ParseTuple(args, "OO", &x_obj, &powers_obj)) { return NULL; }
    OwnedArray x = as_c_array(x_obj, NPY_DOUBLE, 2, "x");
    if (!x) { return NULL; }
    OwnedArray powers = as_c_array(powers_obj, NPY_INT64, 2, "powers");
    if (!powers) { return NULL; }
    const npy_intp q = dim(x, 0), n = dim(x, 1), r = dim(powers, 0);
    if (r > 0 && dim(powers, 1) != n) {
        PyErr_SetString(PyExc_ValueError, "`x` and `powers` must have the same number of columns.");
        return NULL;
    }

    npy_intp dims[2] = {q, r};
    OwnedArray out(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(2, dims, NPY_DOUBLE, 0)));
    if (!out) { return NULL; }

    const double *xd = f64_data(x);
    const int64_t *pd = static_cast<const int64_t *>(PyArray_DATA(powers.get()));
    double *od = f64_data_mut(out);
    {
        GILRelease nogil;
        for (npy_intp i = 0; i < q; i++) {
            for (npy_intp k = 0; k < r; k++) {
                od[i * r + k] = monomial(xd + i * n, pd + k * n, n);
            }
        }
    }
    return reinterpret_cast<PyObject *>(out.release());
}


PyDoc_STRVAR(doc_build_system,
"_build_system(y, d, smoothing, kernel, epsilon, powers)\n"
"\n"
"Build the system used to solve for the RBF interpolant coefficients.\n"
"\n"
"Parameters\n"
"----------\n"
"y : (P, N) float ndarray\n"
"    Data point coordinates.\n"
"d : (P, S) float ndarray\n"
"    Data values at `y`.\n"
"smoothing : (P,) float ndarray\n"
"    Smoothing parameter for each data point.\n"
"kernel : str or PyCapsule\n"
"    Name of a built-in RBF, or a capsule with signature \"double (double)\".\n"
"epsilon : float\n"
"    Shape parameter.\n"
"powers : (R, N) int ndarray\n"
"    The exponents for each monomial in the polynomial.\n"
"\n"
"Returns\n"
"-------\n"
"lhs : (P + R, P + R) float ndarray, Fortran order\n"
"    Left-hand side matrix.\n"
"rhs : (P + R, S) float ndarray, Fortran order\n"
"    Right-hand side matrix.\n"
"shift : (N,) float ndarray\n"
"    Domain shift used to create the polynomial matrix.\n"
"scale : (N,) float ndarray\n"
"    Domain scaling used to create the polynomial matrix.\n");


PyObject *py_build_system(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *y_obj, *d_obj, *smoothing_obj, *powers_obj;
    Kernel kernel = {LINEAR, NULL};
    double epsilon;
    if (!PyArg_ParseTuple(args, "OOOO&dO", &y_obj, &d_obj, &smoothing_obj,
                          kernel_converter, &kernel, &epsilon, &powers_obj)) {
        return NULL;
    }
    OwnedArray y = as_c_array(y_obj, NPY_DOUBLE, 2, "y");
    if (!y) { return NULL; }
    OwnedArray d = as_c_array(d_obj, NPY_DOUBLE, 2, "d");
    if (!d) { return NULL; }
    OwnedArray smoothing = as_c_array(smoothing_obj, NPY_DOUBLE, 1, "smoothing");
    if (!smoothing) { return NULL; }
    OwnedArray powers = as_c_array(powers_obj, NPY_INT64, 2, "powers");
    if (!powers) { return NULL; }

    const npy_intp p = dim(y, 0), n = dim(y, 1), s = dim(d, 1), r = dim(powers, 0);
    if (dim(d, 0) != p || dim(smoothing, 0) != p) {
        PyErr_SetString(PyExc_ValueError,
                        "`y`, `d` and `smoothing` must have the same number of rows.");
        return NULL;
    }
    if (r > 0 && dim(powers, 1) != n) {
        PyErr_SetString(PyExc_ValueError, "`y` and `powers` must have the same number of columns.");
        return NULL;
    }

    const npy_intp m = p + r;
    npy_intp lhs_dims[2] = {m, m}, rhs_dims[2] = {m, s}, vec_dims[1] = {n};
    OwnedArray lhs(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(2, lhs_dims, NPY_DOUBLE, 1)));
    OwnedArray rhs(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(2, rhs_dims, NPY_DOUBLE, 1)));
    OwnedArray shift(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(1, vec_dims, NPY_DOUBLE, 0)));
    OwnedArray scale(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(1, vec_dims, NPY_DOUBLE, 0)));
    if (!lhs || !rhs || !shift || !scale) { return NULL; }

    const double *yd = f64_data(y), *dd = f64_data(d), *sm = f64_data(smoothing);
    const int64_t *pd = static_cast<const int64_t *>(PyArray_DATA(powers.get()));
    double *L = f64_data_mut(lhs), *R = f64_data_mut(rhs);
    double *sh = f64_data_mut(shift), *sc = f64_data_mut(scale);

    // Shift and scale the polynomial domain to [-1, 1]. Like np.min/np.max, a nan
    // in a column propagates to its shift and scale.
    for (npy_intp k = 0; k < n; k++) {
        double lo = p > 0 ? yd[k] : 0.0, hi = lo;
        for (npy_intp i = 1; i < p; i++) {
            const double v = yd[i * n + k];
            if (v < lo || std::isnan(v)) { lo = v; }
            if (v > hi || std::isnan(v)) { hi = v; }
        }
        sh[k] = (hi + lo) / 2.0;
        sc[k] = (hi - lo) / 2.0;
        // A single point, or points that share a coordinate, give zero scale.
        if (sc[k] == 0.0) { sc[k] = 1.0; }
    }

    try {
        const std::vector<double> yeps = scaled(yd, p, n, epsilon);
        const std::vector<double> yhat = normalized(yd, p, n, sh, sc);
        GILRelease nogil;
        // lhs[:p, :p]: the kernel matrix, Fortran order.
        with_kernel(kernel, [&](auto kern) { kernel_matrix(kern, yeps.data(), p, n, L, 1, m); });
        // lhs[:p, p:] and lhs[p:, :p]: the polynomial matrix and its transpose.
        for (npy_intp k = 0; k < r; k++) {
            for (npy_intp i = 0; i < p; i++) {
                const double v = monomial(yhat.data() + i * n, pd + k * n, n);
                L[i + (p + k) * m] = v;
                L[(p + k) + i * m] = v;
            }
        }
        // lhs[p:, p:] = 0
        for (npy_intp b = p; b < m; b++) {
            for (npy_intp a = p; a < m; a++) { L[a + b * m] = 0.0; }
        }
        for (npy_intp i = 0; i < p; i++) { L[i + i * m] += sm[i]; }

        // rhs[:p] = d, rhs[p:] = 0
        for (npy_intp j = 0; j < s; j++) {
            for (npy_intp i = 0; i < p; i++) { R[i + j * m] = dd[i * s + j]; }
            for (npy_intp i = p; i < m; i++) { R[i + j * m] = 0.0; }
        }
    } catch (...) {
        return set_cpp_error();
    }

    return Py_BuildValue("NNNN", lhs.release(), rhs.release(), shift.release(), scale.release());
}


PyDoc_STRVAR(doc_build_evaluation_coefficients,
"_build_evaluation_coefficients(x, y, kernel, epsilon, powers, shift, scale)\n"
"\n"
"Construct the coefficients needed to evaluate the RBF.\n"
"\n"
"Parameters\n"
"----------\n"
"x : (Q, N) float ndarray\n"
"    Evaluation point coordinates.\n"
"y : (P, N) float ndarray\n"
"    Data point coordinates.\n"
"kernel : str or PyCapsule\n"
"    Name of a built-in RBF, or a capsule with signature \"double (double)\".\n"
"epsilon : float\n"
"    Shape parameter.\n"
"powers : (R, N) int ndarray\n"
"    The exponents for each monomial in the polynomial.\n"
"shift : (N,) float ndarray\n"
"    Shifts the polynomial domain for numerical stability.\n"
"scale : (N,) float ndarray\n"
"    Scales the polynomial domain for numerical stability.\n"
"\n"
"Returns\n"
"-------\n"
"(Q, P + R) float ndarray\n");


PyObject *py_build_evaluation_coefficients(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *x_obj, *y_obj, *powers_obj, *shift_obj, *scale_obj;
    Kernel kernel = {LINEAR, NULL};
    double epsilon;
    if (!PyArg_ParseTuple(args, "OOO&dOOO", &x_obj, &y_obj, kernel_converter, &kernel,
                          &epsilon, &powers_obj, &shift_obj, &scale_obj)) {
        return NULL;
    }
    OwnedArray x = as_c_array(x_obj, NPY_DOUBLE, 2, "x");
    if (!x) { return NULL; }
    OwnedArray y = as_c_array(y_obj, NPY_DOUBLE, 2, "y");
    if (!y) { return NULL; }
    OwnedArray powers = as_c_array(powers_obj, NPY_INT64, 2, "powers");
    if (!powers) { return NULL; }
    OwnedArray shift = as_c_array(shift_obj, NPY_DOUBLE, 1, "shift");
    if (!shift) { return NULL; }
    OwnedArray scale = as_c_array(scale_obj, NPY_DOUBLE, 1, "scale");
    if (!scale) { return NULL; }

    const npy_intp q = dim(x, 0), n = dim(x, 1), p = dim(y, 0), r = dim(powers, 0);
    if (dim(y, 1) != n || dim(shift, 0) != n || dim(scale, 0) != n
            || (r > 0 && dim(powers, 1) != n)) {
        PyErr_SetString(PyExc_ValueError,
                        "`x`, `y`, `powers`, `shift` and `scale` must have the same "
                        "number of columns.");
        return NULL;
    }

    const npy_intp m = p + r;
    npy_intp dims[2] = {q, m};
    OwnedArray out(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(2, dims, NPY_DOUBLE, 0)));
    if (!out) { return NULL; }

    const double *xd = f64_data(x), *yd = f64_data(y);
    const int64_t *pd = static_cast<const int64_t *>(PyArray_DATA(powers.get()));
    double *od = f64_data_mut(out);

    try {
        const std::vector<double> xeps = scaled(xd, q, n, epsilon);
        const std::vector<double> yeps = scaled(yd, p, n, epsilon);
        const std::vector<double> xhat = normalized(xd, q, n, f64_data(shift), f64_data(scale));
        GILRelease nogil;
        with_kernel(kernel, [&](auto kern) {
            for (npy_intp i = 0; i < q; i++) {
                const double *xi = xeps.data() + i * n;
                double *row = od + i * m;
                for (npy_intp j = 0; j < p; j++) {
                    row[j] = kern(distance(xi, yeps.data() + j * n, n));
                }
            }
        });
        for (npy_intp i = 0; i < q; i++) {
            double *row = od + i * m + p;
            for (npy_intp k = 0; k < r; k++) {
                row[k] = monomial(xhat.data() + i * n, pd + k * n, n);
            }
        }
    } catch (...) {
        return set_cpp_error();
    }
    return reinterpret_cast<PyObject *>(out.release());
}


PyDoc_STRVAR(doc_kernel_capsule,
"_kernel_capsule(name)\n"
"\n"
"Return a built-in RBF as a PyCapsule with signature \"double (double)\".\n"
"\n"
"The capsule can be wrapped in a `scipy.LowLevelCallable`; it exists to test\n"
"the code path for user-supplied kernels with compiled code.\n");


PyObject *py_kernel_capsule(PyObject *Py_UNUSED(self), PyObject *args)
{
    Kernel kernel = {LINEAR, NULL};
    if (!PyArg_ParseTuple(args, "O&", kernel_converter, &kernel)) { return NULL; }
    return PyCapsule_New(reinterpret_cast<void *>(kernel.func), KERNEL_SIGNATURE, NULL);
}


PyMethodDef rbfinterp_methods[] = {
    {"_kernel_matrix", py_kernel_matrix, METH_VARARGS, doc_kernel_matrix},
    {"_polynomial_matrix", py_polynomial_matrix, METH_VARARGS, doc_polynomial_matrix},
    {"_build_system", py_build_system, METH_VARARGS, doc_build_system},
    {"_build_evaluation_coefficients", py_build_evaluation_coefficients, METH_VARARGS,
     doc_build_evaluation_coefficients},
    {"_kernel_capsule", py_kernel_capsule, METH_VARARGS, doc_kernel_capsule},
    {NULL, NULL, 0, NULL}
};


int rbfinterp_module_exec(PyObject *Py_UNUSED(module))
{
    if (_import_array() < 0) { return -1; }
    return 0;
}


PyModuleDef_Slot rbfinterp_module_slots[] = {
    {Py_mod_exec, (void *)rbfinterp_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};


// Designated initializers require C++20; MSVC /std:c++17 rejects them.
PyModuleDef rbfinterp_module = {
    PyModuleDef_HEAD_INIT,                /* m_base */
    "_rbfinterp_cpp",                     /* m_name */
    NULL,                                 /* m_doc */
    0,                                    /* m_size */
    rbfinterp_methods,                    /* m_methods */
    rbfinterp_module_slots,               /* m_slots */
    NULL,                                 /* m_traverse */
    NULL,                                 /* m_clear */
    NULL                                  /* m_free */
};

}  // namespace


PyMODINIT_FUNC
PyInit__rbfinterp_cpp(void)
{
    return PyModuleDef_Init(&rbfinterp_module);
}
