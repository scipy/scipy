#include <cassert>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

#include "sf_error.h"
#include "special/mdspan.h"

// Initializes Python and NumPy.
inline bool SpecFun_Initialize() {
    Py_Initialize();

    import_array();
    if (PyErr_Occurred() != nullptr) {
        return false; // import array failed
    }

    import_umath();
    if (PyErr_Occurred() != nullptr) {
        return false; // import umath failed
    }

    return true;
}

// Deduces the number of arguments of a callable F.
template <typename Func>
struct arity_of;

template <typename Res, typename... Args>
struct arity_of<Res(Args...)> {
    static constexpr size_t value = sizeof...(Args);
};

template <typename Func>
constexpr size_t arity_of_v = arity_of<Func>::value;

template <typename F, size_t I>
struct argument_rank;

template <typename Res, typename... Args, size_t I>
struct argument_rank<Res(Args...), I> {
    static constexpr size_t value = 0;
};

template <typename Res, typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy, typename... Args>
struct argument_rank<Res(std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>, Args...), 0> {
    static constexpr size_t value = Extents::rank();
};

template <typename Res, typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy, typename... Args,
          size_t I>
struct argument_rank<Res(std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>, Args...), I> {
    static constexpr size_t value = Extents::rank() + argument_rank<Res(Args...), I - 1>::value;
};

template <typename T, size_t N>
constexpr size_t argument_rank_v = argument_rank<T, N>::value;

// Maps a C++ type to a NumPy type identifier.
template <typename T>
struct npy_type;

template <>
struct npy_type<bool> {
    static constexpr int value = NPY_BOOL;
};

template <>
struct npy_type<long int> {
    static constexpr int value = NPY_LONG;
};

template <>
struct npy_type<float> {
    static constexpr int value = NPY_FLOAT32;
};

template <>
struct npy_type<double> {
    static constexpr int value = NPY_FLOAT64;
};

template <>
struct npy_type<int> {
    static constexpr int value = NPY_INT;
};

template <typename T>
struct npy_type<T *> {
    static constexpr int value = npy_type<T>::value;
};

template <>
struct npy_type<std::complex<float>> {
    static constexpr int value = NPY_COMPLEX64;
};

template <>
struct npy_type<std::complex<double>> {
    static constexpr int value = NPY_COMPLEX128;
};

template <>
struct npy_type<npy_cdouble> {
    static constexpr int value = NPY_COMPLEX128;
};

template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
struct npy_type<std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>> {
    static constexpr int value = npy_type<T>::value;
};

// Sets the value dst to be the value of type T at src
template <typename T>
void from_pointer(char *src, T &dst, const npy_intp *dimensions, const npy_intp *steps) {
    dst = *reinterpret_cast<T *>(src);
}

template <typename T, typename Extents, typename AccessorPolicy>
void from_pointer(char *src, std::mdspan<T, Extents, std::layout_stride, AccessorPolicy> &dst,
                  const npy_intp *dimensions, const npy_intp *steps) {
    std::array<ptrdiff_t, Extents::rank()> strides;
    for (npy_uintp i = 0; i < strides.size(); ++i) {
        strides[i] = steps[i] / sizeof(T);
    }

    std::array<ptrdiff_t, Extents::rank()> exts;
    for (npy_uintp i = 0; i < exts.size(); ++i) {
        exts[i] = dimensions[i];
    }

    dst = {reinterpret_cast<T *>(src), {exts, strides}};
}

// Sets the pointer dst to be the pointer of type T at src (helps for out arguments)
template <typename T>
void from_pointer(char *src, T *&dst, const npy_intp *dimensions, const npy_intp *steps) {
    dst = reinterpret_cast<T *>(src);
}

template <typename T>
T from_pointer(char *src, const npy_intp *dimensions, const npy_intp *steps) {
    T dst;
    from_pointer(src, dst, dimensions, steps);

    return dst;
}

struct SpecFun_UFuncData {
    void *func;
    const char *name;
};

template <typename Func, typename Indices = std::make_index_sequence<arity_of_v<Func>>>
struct ufunc_traits;

template <typename Res, typename... Args, size_t... I>
struct ufunc_traits<Res(Args...), std::index_sequence<I...>> {
    static constexpr char types[sizeof...(Args) + 1] = {npy_type<Args>::value..., npy_type<Res>::value};

    static void loop_func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        Res (*func)(Args...) = reinterpret_cast<Res (*)(Args...)>(static_cast<SpecFun_UFuncData *>(data)->func);
        const char *func_name = static_cast<SpecFun_UFuncData *>(data)->name;

        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            *reinterpret_cast<Res *>(args[sizeof...(Args)]) =
                func(from_pointer<Args>(args[I], dimensions + 1, steps + sizeof...(Args) + 1)...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
            args[sizeof...(Args)] += steps[sizeof...(Args)]; // output
        }

        sf_error_check_fpe(func_name);
    }
};

template <typename... Args, size_t... I>
struct ufunc_traits<void(Args...), std::index_sequence<I...>> {
    static constexpr char types[sizeof...(Args)] = {npy_type<Args>::value...};

    static void loop_func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        void (*func)(Args...) = reinterpret_cast<void (*)(Args...)>(static_cast<SpecFun_UFuncData *>(data)->func);
        const char *func_name = static_cast<SpecFun_UFuncData *>(data)->name;

        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            func(from_pointer<Args>(args[I], dimensions + 1,
                                    steps + sizeof...(Args) + argument_rank_v<void(Args...), I>)...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
        }

        sf_error_check_fpe(func_name);
    }
};

class SpecFun_Func {
    bool m_has_return;
    int m_nin_and_nout;
    void *m_func;
    PyUFuncGenericFunction m_loop_func;
    std::unique_ptr<char[]> m_types;

  public:
    template <typename Res, typename... Args>
    SpecFun_Func(Res (*func)(Args... args))
        : m_has_return(!std::is_void_v<Res>), m_nin_and_nout(sizeof...(Args) + m_has_return),
          m_func(reinterpret_cast<void *>(func)), m_loop_func(ufunc_traits<Res(Args...)>::loop_func),
          m_types(new char[m_nin_and_nout]) {
        std::copy(ufunc_traits<Res(Args...)>::types, ufunc_traits<Res(Args...)>::types + m_nin_and_nout, m_types.get());
    }

    int nin_and_nout() const { return m_nin_and_nout; }

    bool has_return() const { return m_has_return; }

    void *func() const { return m_func; }

    PyUFuncGenericFunction loop_func() const { return m_loop_func; }

    char *types() const { return m_types.get(); }
};

class SpecFun_UFunc {
    int m_ntypes;
    int m_nin_and_nout;
    std::unique_ptr<PyUFuncGenericFunction[]> m_func;
    std::unique_ptr<void *[]> m_data;
    std::unique_ptr<char[]> m_types;
    std::unique_ptr<SpecFun_UFuncData[]> m_data_alloc;

  public:
    SpecFun_UFunc(std::initializer_list<SpecFun_Func> func, const char *name)
        : m_ntypes(func.size()), m_nin_and_nout(func.begin()->nin_and_nout()),
          m_func(new PyUFuncGenericFunction[m_ntypes]), m_data(new void *[m_ntypes]),
          m_types(new char[m_ntypes * m_nin_and_nout]), m_data_alloc(new SpecFun_UFuncData[m_ntypes]) {
        for (auto it = func.begin(); it != func.end(); ++it) {
            size_t i = it - func.begin();

            m_func[i] = it->loop_func();
            m_data[i] = m_data_alloc.get() + i;
            std::copy(it->types(), it->types() + m_nin_and_nout, m_types.get() + i * m_nin_and_nout);

            m_data_alloc[i] = {it->func(), name};
        }
    }

    int nin_and_nout() const { return m_nin_and_nout; }

    int ntypes() const { return m_ntypes; }

    PyUFuncGenericFunction *func() const { return m_func.get(); }

    void **data() const { return m_data.get(); }

    char *types() const { return m_types.get(); }
};

PyObject *SpecFun_NewUFunc(std::initializer_list<SpecFun_Func> func, int nout, const char *name, const char *doc) {
    static std::vector<SpecFun_UFunc> ufuncs;

    for (auto it = func.begin(); it != func.end(); ++it) {
        if (it->nin_and_nout() != func.begin()->nin_and_nout()) {
            PyErr_SetString(PyExc_RuntimeError, "all functions must have the same number of arguments");
            return nullptr;
        }
        if (it->has_return() != func.begin()->has_return()) {
            PyErr_SetString(PyExc_RuntimeError, "all functions must be void if any function is");
            return nullptr;
        }
    }

    const SpecFun_UFunc &ufunc = ufuncs.emplace_back(func, name);
    return PyUFunc_FromFuncAndData(ufunc.func(), ufunc.data(), ufunc.types(), ufunc.ntypes(),
                                   ufunc.nin_and_nout() - nout, nout, PyUFunc_None, name, doc, 0);
}

PyObject *SpecFun_NewUFunc(std::initializer_list<SpecFun_Func> func, const char *name, const char *doc) {
    return SpecFun_NewUFunc(func, func.begin()->has_return(), name, doc);
}

PyObject *SpecFun_NewGUFunc(std::initializer_list<SpecFun_Func> func, int nout, const char *name, const char *doc,
                            const char *signature) {
    static std::vector<SpecFun_UFunc> ufuncs;

    for (auto it = func.begin(); it != func.end(); ++it) {
        if (it->nin_and_nout() != func.begin()->nin_and_nout()) {
            PyErr_SetString(PyExc_RuntimeError, "all functions must have the same number of arguments");
            return nullptr;
        }
        if (it->has_return() != func.begin()->has_return()) {
            PyErr_SetString(PyExc_RuntimeError, "all functions must be void if any function is");
            return nullptr;
        }
    }

    const SpecFun_UFunc &ufunc = ufuncs.emplace_back(func, name);
    return PyUFunc_FromFuncAndDataAndSignature(ufunc.func(), ufunc.data(), ufunc.types(), ufunc.ntypes(),
                                               ufunc.nin_and_nout() - nout, nout, PyUFunc_None, name, doc, 0,
                                               signature);
}

PyObject *SpecFun_NewGUFunc(std::initializer_list<SpecFun_Func> func, const char *name, const char *doc,
                            const char *signature) {
    return SpecFun_NewGUFunc(func, func.begin()->has_return(), name, doc, signature);
}
