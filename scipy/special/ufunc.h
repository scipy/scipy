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

// Maps a C++ type to a NumPy type
template <typename T>
struct npy_type;

template <>
struct npy_type<bool> {
    using type = npy_bool;
};

template <>
struct npy_type<char> {
    using type = npy_byte;
};

template <>
struct npy_type<short> {
    using type = npy_short;
};

template <>
struct npy_type<int> {
    using type = npy_int;
};

template <>
struct npy_type<long> {
    using type = npy_long;
};

template <>
struct npy_type<long long> {
    using type = npy_longlong;
};

template <>
struct npy_type<unsigned char> {
    using type = npy_ubyte;
};

template <>
struct npy_type<unsigned short> {
    using type = npy_ushort;
};

template <>
struct npy_type<unsigned int> {
    using type = npy_uint;
};

template <>
struct npy_type<unsigned long> {
    using type = npy_ulong;
};

template <>
struct npy_type<unsigned long long> {
    using type = npy_ulonglong;
};

template <>
struct npy_type<float> {
    using type = npy_float;
};

template <>
struct npy_type<double> {
    using type = npy_double;
};

template <>
struct npy_type<long double> {
    using type = npy_longdouble;
};

template <>
struct npy_type<std::complex<float>> {
    using type = npy_cfloat;
};

template <>
struct npy_type<std::complex<double>> {
    using type = npy_cdouble;
};

template <>
struct npy_type<std::complex<long double>> {
    using type = npy_clongdouble;
};

template <typename T>
using npy_type_t = typename npy_type<T>::type;

// Maps a C++ type to a NumPy type number
template <typename T>
struct npy_typenum {
    static constexpr int value = npy_typenum<npy_type_t<T>>::value;
};

// We need to specialise for bool as npy_bool is defined as npy_ubyte
template <>
struct npy_typenum<bool> {
    static constexpr int value = NPY_BOOL;
};

template <>
struct npy_typenum<npy_byte> {
    static constexpr int value = NPY_BYTE;
};

template <>
struct npy_typenum<npy_short> {
    static constexpr int value = NPY_SHORT;
};

template <>
struct npy_typenum<npy_int> {
    static constexpr int value = NPY_INT;
};

template <>
struct npy_typenum<npy_long> {
    static constexpr int value = NPY_LONG;
};

template <>
struct npy_typenum<npy_longlong> {
    static constexpr int value = NPY_LONGLONG;
};

template <>
struct npy_typenum<npy_ubyte> {
    static constexpr int value = NPY_UBYTE;
};

template <>
struct npy_typenum<npy_ushort> {
    static constexpr int value = NPY_USHORT;
};

template <>
struct npy_typenum<npy_uint> {
    static constexpr int value = NPY_UINT;
};

template <>
struct npy_typenum<npy_ulong> {
    static constexpr int value = NPY_ULONG;
};

template <>
struct npy_typenum<npy_ulonglong> {
    static constexpr int value = NPY_ULONGLONG;
};

template <>
struct npy_typenum<npy_float> {
    static constexpr int value = NPY_FLOAT;
};

template <>
struct npy_typenum<npy_double> {
    static constexpr int value = NPY_DOUBLE;
};

// When NPY_SIZEOF_LONGDOUBLE == NPY_SIZEOF_DOUBLE, npy_longdouble is defined as npy_double
// See https://github.com/numpy/numpy/blob/main/numpy/_core/include/numpy/npy_common.h
#if (NPY_SIZEOF_LONGDOUBLE != NPY_SIZEOF_DOUBLE)
template <>
struct npy_typenum<npy_longdouble> {
    static constexpr int value = NPY_LONGDOUBLE;
};
#endif

template <>
struct npy_typenum<npy_cfloat> {
    static constexpr int value = NPY_CFLOAT;
};

template <>
struct npy_typenum<npy_cdouble> {
    static constexpr int value = NPY_CDOUBLE;
};

template <>
struct npy_typenum<npy_clongdouble> {
    static constexpr int value = NPY_CLONGDOUBLE;
};

template <typename T>
struct npy_typenum<T *> {
    static constexpr int value = npy_typenum<T>::value;
};

template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
struct npy_typenum<std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>> {
    static constexpr int value = npy_typenum<T>::value;
};

template <typename T>
inline constexpr int npy_typenum_v = npy_typenum<T>::value;

// Sets the value dst to be the value of type T at src
template <typename T>
void from_npy_bytes(char *src, T &dst, const npy_intp *dimensions, const npy_intp *steps) {
    dst = *reinterpret_cast<npy_type_t<T> *>(src);
}

template <typename T>
void from_npy_bytes(char *src, std::complex<T> &dst, const npy_intp *dimensions, const npy_intp *steps) {
    dst.real(*reinterpret_cast<npy_type_t<T> *>(src));
    dst.imag(*reinterpret_cast<npy_type_t<T> *>(src + sizeof(T)));
}

// Sets the pointer dst to be the pointer of type T at src (helps for out arguments)
template <typename T>
void from_npy_bytes(char *src, T *&dst, const npy_intp *dimensions, const npy_intp *steps) {
    static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

    dst = reinterpret_cast<T *>(src);
}

template <typename T, typename Extents, typename AccessorPolicy>
void from_npy_bytes(char *src, std::mdspan<T, Extents, std::layout_stride, AccessorPolicy> &dst,
                    const npy_intp *dimensions, const npy_intp *steps) {
    static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

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

template <typename T>
T from_npy_bytes(char *src, const npy_intp *dimensions, const npy_intp *steps) {
    T dst;
    from_npy_bytes(src, dst, dimensions, steps);

    return dst;
}

template <typename Func>
struct ufunc_data {
    const char *name;
    Func func;
};

template <typename Func, typename Indices = std::make_index_sequence<arity_of_v<Func>>>
struct ufunc_traits;

template <typename Res, typename... Args, size_t... I>
struct ufunc_traits<Res(Args...), std::index_sequence<I...>> {
    static constexpr char types[sizeof...(Args) + 1] = {npy_typenum_v<Args>..., npy_typenum_v<Res>};

    static void loop_func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        Res (*func)(Args...) = static_cast<ufunc_data<Res (*)(Args...)> *>(data)->func;
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            *reinterpret_cast<Res *>(args[sizeof...(Args)]) =
                func(from_npy_bytes<Args>(args[I], dimensions + 1, steps + sizeof...(Args) + 1)...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
            args[sizeof...(Args)] += steps[sizeof...(Args)]; // output
        }

        const char *name = static_cast<ufunc_data<Res (*)(Args...)> *>(data)->name;
        sf_error_check_fpe(name);
    }
};

template <typename... Args, size_t... I>
struct ufunc_traits<void(Args...), std::index_sequence<I...>> {
    static constexpr char types[sizeof...(Args)] = {npy_typenum_v<Args>...};

    static void loop_func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        void (*func)(Args...) = static_cast<ufunc_data<void (*)(Args...)> *>(data)->func;
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            func(from_npy_bytes<Args>(args[I], dimensions + 1,
                                      steps + sizeof...(Args) + argument_rank_v<void(Args...), I>)...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
        }

        const char *name = static_cast<ufunc_data<void (*)(Args...)> *>(data)->name;
        sf_error_check_fpe(name);
    }
};

class SpecFun_Func {
    bool m_has_return;
    int m_nin_and_nout;
    PyUFuncGenericFunction m_loop_func;
    std::unique_ptr<char[]> m_types;
    size_t m_data_size;
    std::unique_ptr<std::byte[]> m_data;

  public:
    template <typename Res, typename... Args>
    SpecFun_Func(Res (*func)(Args... args))
        : m_has_return(!std::is_void_v<Res>), m_nin_and_nout(sizeof...(Args) + m_has_return),
          m_loop_func(ufunc_traits<Res(Args...)>::loop_func), m_types(new char[m_nin_and_nout]),
          m_data_size(sizeof(ufunc_data<Res (*)(Args...)>)), m_data(new std::byte[m_data_size]) {
        std::copy(ufunc_traits<Res(Args...)>::types, ufunc_traits<Res(Args...)>::types + m_nin_and_nout, m_types.get());
        reinterpret_cast<ufunc_data<Res (*)(Args...)> *>(m_data.get())->func = func;
    }

    int nin_and_nout() const { return m_nin_and_nout; }

    bool has_return() const { return m_has_return; }

    PyUFuncGenericFunction loop_func() const { return m_loop_func; }

    void *data() const { return reinterpret_cast<void *>(m_data.get()); }

    size_t data_size() const { return m_data_size; }

    char *types() const { return m_types.get(); }
};

class SpecFun_UFunc {
    int m_ntypes;
    int m_nin_and_nout;
    std::unique_ptr<PyUFuncGenericFunction[]> m_func;
    std::unique_ptr<std::byte *[]> m_data;
    std::unique_ptr<char[]> m_types;

  public:
    SpecFun_UFunc(std::initializer_list<SpecFun_Func> func, const char *name)
        : m_ntypes(func.size()), m_nin_and_nout(func.begin()->nin_and_nout()),
          m_func(new PyUFuncGenericFunction[m_ntypes]), m_data(new std::byte *[m_ntypes]),
          m_types(new char[m_ntypes * m_nin_and_nout]) {
        for (auto it = func.begin(); it != func.end(); ++it) {
            size_t i = it - func.begin();

            m_func[i] = it->loop_func();

            m_data[i] = new std::byte[it->data_size()];
            std::memcpy(m_data[i], it->data(), it->data_size());
            std::memcpy(m_data[i], &name, sizeof(const char *)); // overwrite the first member with the name pointer

            std::copy(it->types(), it->types() + m_nin_and_nout, m_types.get() + i * m_nin_and_nout);
        }
    }

    SpecFun_UFunc(SpecFun_UFunc &&) = default;

    ~SpecFun_UFunc() {
        if (m_data) {
            for (int i = 0; i < m_ntypes; ++i) {
                delete[] m_data[i];
            }
        }
    }

    int nin_and_nout() const { return m_nin_and_nout; }

    int ntypes() const { return m_ntypes; }

    PyUFuncGenericFunction *func() const { return m_func.get(); }

    void **data() const { return reinterpret_cast<void **>(m_data.get()); }

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
