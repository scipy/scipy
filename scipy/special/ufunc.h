#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cassert>
#include <cstring>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

#include "sf_error.h"
#include "special/third_party/kokkos/mdspan.hpp"


// This is std::accumulate, but that is not constexpr until C++20
template <typename InputIt, typename T>
constexpr T initializer_accumulate(InputIt first, InputIt last, T init) {
    for (InputIt it = first; it != last; ++it) {
        init = std::move(init) + *it;
    }

    return init;
}

// Deduces the number of arguments of a callable F.
template <typename Func>
struct arity_of;

template <typename Res, typename... Args>
struct arity_of<Res (*)(Args...)> {
    static constexpr size_t value = sizeof...(Args);
};

template <typename Func>
constexpr size_t arity_of_v = arity_of<Func>::value;

template <typename Func>
struct has_return;

template <typename Res, typename... Args>
struct has_return<Res (*)(Args...)> {
    static constexpr bool value = true;
};

template <typename... Args>
struct has_return<void (*)(Args...)> {
    static constexpr bool value = false;
};

template <typename Func>
constexpr size_t has_return_v = has_return<Func>::value;

template <typename T>
struct rank_of {
    static constexpr size_t value = 0;
};

template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
struct rank_of<std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>> {
    static constexpr size_t value = Extents::rank();
};

template <typename T>
inline constexpr size_t rank_of_v = rank_of<T>::value;

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
    static constexpr NPY_TYPES value = npy_typenum<npy_type_t<T>>::value;
};

// We need to specialise for bool as npy_bool is defined as npy_ubyte
template <>
struct npy_typenum<bool> {
    static constexpr NPY_TYPES value = NPY_BOOL;
};

template <>
struct npy_typenum<npy_byte> {
    static constexpr NPY_TYPES value = NPY_BYTE;
};

template <>
struct npy_typenum<npy_short> {
    static constexpr NPY_TYPES value = NPY_SHORT;
};

template <>
struct npy_typenum<npy_int> {
    static constexpr NPY_TYPES value = NPY_INT;
};

template <>
struct npy_typenum<npy_long> {
    static constexpr NPY_TYPES value = NPY_LONG;
};

template <>
struct npy_typenum<npy_longlong> {
    static constexpr NPY_TYPES value = NPY_LONGLONG;
};

template <>
struct npy_typenum<npy_ubyte> {
    static constexpr NPY_TYPES value = NPY_UBYTE;
};

template <>
struct npy_typenum<npy_ushort> {
    static constexpr NPY_TYPES value = NPY_USHORT;
};

template <>
struct npy_typenum<npy_uint> {
    static constexpr NPY_TYPES value = NPY_UINT;
};

template <>
struct npy_typenum<npy_ulong> {
    static constexpr NPY_TYPES value = NPY_ULONG;
};

template <>
struct npy_typenum<npy_ulonglong> {
    static constexpr NPY_TYPES value = NPY_ULONGLONG;
};

template <>
struct npy_typenum<npy_float> {
    static constexpr NPY_TYPES value = NPY_FLOAT;
};

template <>
struct npy_typenum<npy_double> {
    static constexpr NPY_TYPES value = NPY_DOUBLE;
};

// When NPY_SIZEOF_LONGDOUBLE == NPY_SIZEOF_DOUBLE, npy_longdouble is defined as npy_double
// See https://github.com/numpy/numpy/blob/main/numpy/_core/include/numpy/npy_common.h
#if (NPY_SIZEOF_LONGDOUBLE != NPY_SIZEOF_DOUBLE)
template <>
struct npy_typenum<npy_longdouble> {
    static constexpr NPY_TYPES value = NPY_LONGDOUBLE;
};
#endif

template <>
struct npy_typenum<npy_cfloat> {
    static constexpr NPY_TYPES value = NPY_CFLOAT;
};

template <>
struct npy_typenum<npy_cdouble> {
    static constexpr NPY_TYPES value = NPY_CDOUBLE;
};

template <>
struct npy_typenum<npy_clongdouble> {
    static constexpr NPY_TYPES value = NPY_CLONGDOUBLE;
};

template <typename T>
struct npy_typenum<T *> {
    static constexpr NPY_TYPES value = npy_typenum<T>::value;
};

template <typename T>
struct npy_typenum<T &> {
    static constexpr NPY_TYPES value = npy_typenum<T>::value;
};

template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
struct npy_typenum<std::mdspan<T, Extents, LayoutPolicy, AccessorPolicy>> {
    static constexpr NPY_TYPES value = npy_typenum<T>::value;
};

template <typename T>
inline constexpr NPY_TYPES npy_typenum_v = npy_typenum<T>::value;

template <typename T>
struct npy_traits {
    static T get(char *src, const npy_intp *dimensions, const npy_intp *steps) { return *reinterpret_cast<T *>(src); }

    static void set(char *dst, const T &src) { *reinterpret_cast<npy_type_t<T> *>(dst) = src; }
};

template <typename T>
struct npy_traits<std::complex<T>> {
    static std::complex<T> get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
        return *reinterpret_cast<std::complex<T> *>(src);
    }

    static void set(char *dst, const std::complex<T> &src) {
        *reinterpret_cast<npy_type_t<T> *>(dst) = std::real(src);
        *reinterpret_cast<npy_type_t<T> *>(dst + sizeof(T)) = std::imag(src);
    }
};

template <typename T>
struct npy_traits<T *> {
    static T *get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
        static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

        return reinterpret_cast<T *>(src);
    }
};

template <typename T>
struct npy_traits<T &> {
    static T &get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
        static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

        return *reinterpret_cast<T *>(src);
    }
};

template <typename T, typename Extents, typename AccessorPolicy>
struct npy_traits<std::mdspan<T, Extents, std::layout_stride, AccessorPolicy>> {
    static std::mdspan<T, Extents, std::layout_stride, AccessorPolicy>
    get(char *src, const npy_intp *dimensions, const npy_intp *steps) {
        static_assert(sizeof(T) == sizeof(npy_type_t<T>), "NumPy type has different size than argument type");

        std::array<ptrdiff_t, Extents::rank()> strides;
        for (npy_uintp i = 0; i < strides.size(); ++i) {
            strides[i] = steps[i] / sizeof(T);
        }

        std::array<ptrdiff_t, Extents::rank()> exts;
        for (npy_uintp i = 0; i < exts.size(); ++i) {
            exts[i] = dimensions[i];
        }

        return {reinterpret_cast<T *>(src), {exts, strides}};
    }
};

struct base_ufunc_data {
    const char *name;
};

template <typename Func>
struct ufunc_data : base_ufunc_data {
    Func func;
};

template <typename Func, typename Indices = std::make_index_sequence<arity_of_v<Func>>>
struct ufunc_traits;

template <typename Res, typename... Args, size_t... I>
struct ufunc_traits<Res (*)(Args...), std::index_sequence<I...>> {
    static constexpr char types[sizeof...(Args) + 1] = {npy_typenum_v<Args>..., npy_typenum_v<Res>};

    static constexpr size_t ranks[sizeof...(Args) + 1] = {rank_of_v<Args>..., rank_of_v<Res>};

    static constexpr size_t steps_offsets[sizeof...(Args) + 1] = {
        initializer_accumulate(ranks, ranks + I, sizeof...(Args) + 1)...,
        initializer_accumulate(ranks, ranks + sizeof...(Args) + 1, sizeof...(Args) + 1)
    };

    static void loop(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        Res (*func)(Args...) = static_cast<ufunc_data<Res (*)(Args...)> *>(data)->func;
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            Res res = func(npy_traits<Args>::get(args[I], dimensions + 1, steps + steps_offsets[I])...);
            npy_traits<Res>::set(args[sizeof...(Args)], res); // assign to the output pointer

            for (npy_uintp j = 0; j <= sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
        }

        const char *name = static_cast<ufunc_data<Res (*)(Args...)> *>(data)->name;
        sf_error_check_fpe(name);
    }
};

template <typename... Args, size_t... I>
struct ufunc_traits<void (*)(Args...), std::index_sequence<I...>> {
    static constexpr char types[sizeof...(Args)] = {npy_typenum_v<Args>...};

    static constexpr size_t ranks[sizeof...(Args)] = {rank_of_v<Args>...};

    static constexpr size_t steps_offsets[sizeof...(Args)] = {
        initializer_accumulate(ranks, ranks + I, sizeof...(Args))...
    };

    static void loop(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        void (*func)(Args...) = static_cast<ufunc_data<void (*)(Args...)> *>(data)->func;
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            func(npy_traits<Args>::get(args[I], dimensions + 1, steps + steps_offsets[I])...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
        }

        const char *name = static_cast<ufunc_data<void (*)(Args...)> *>(data)->name;
        sf_error_check_fpe(name);
    }
};

class SpecFun_UFunc {
  public:
    using data_handle_type = void *;
    using data_deleter_type = void (*)(void *);

  private:
    // This is an internal class designed only to help construction from an initializer list of functions
    struct SpecFun_Func {
        bool has_return;
        int nin_and_nout;
        PyUFuncGenericFunction func;
        data_handle_type data;
        data_deleter_type data_deleter;
        const char *types;

        template <typename Func>
        SpecFun_Func(Func func)
            : has_return(has_return_v<Func>), nin_and_nout(arity_of_v<Func> + has_return),
              func(ufunc_traits<Func>::loop), data(new ufunc_data<Func>{{nullptr}, func}),
              data_deleter([](void *ptr) { delete static_cast<ufunc_data<Func> *>(ptr); }),
              types(ufunc_traits<Func>::types) {}
    };

    int m_ntypes;
    bool m_has_return;
    int m_nin_and_nout;
    std::unique_ptr<PyUFuncGenericFunction[]> m_func;
    std::unique_ptr<data_handle_type[]> m_data;
    std::unique_ptr<data_deleter_type[]> m_data_deleters;
    std::unique_ptr<char[]> m_types;

  public:
    SpecFun_UFunc(std::initializer_list<SpecFun_Func> func)
        : m_ntypes(func.size()), m_has_return(func.begin()->has_return), m_nin_and_nout(func.begin()->nin_and_nout),
          m_func(new PyUFuncGenericFunction[m_ntypes]), m_data(new data_handle_type[m_ntypes]),
          m_data_deleters(new data_deleter_type[m_ntypes]), m_types(new char[m_ntypes * m_nin_and_nout]) {
        for (auto it = func.begin(); it != func.end(); ++it) {
            if (it->nin_and_nout != m_nin_and_nout) {
                PyErr_SetString(PyExc_RuntimeError, "all functions must have the same number of arguments");
            }
            if (it->has_return != m_has_return) {
                PyErr_SetString(PyExc_RuntimeError, "all functions must be void if any function is");
            }

            size_t i = it - func.begin();
            m_func[i] = it->func;
            m_data[i] = it->data;
            m_data_deleters[i] = it->data_deleter;
            std::memcpy(m_types.get() + i * m_nin_and_nout, it->types, m_nin_and_nout);
        }
    }

    SpecFun_UFunc(SpecFun_UFunc &&other) = default;

    ~SpecFun_UFunc() {
        if (m_data) {
            for (int i = 0; i < m_ntypes; ++i) {
                data_deleter_type data_deleter = m_data_deleters[i];
                data_deleter(m_data[i]);
            }
        }
    }

    int ntypes() const { return m_ntypes; }

    bool has_return() const { return m_has_return; }

    int nin_and_nout() const { return m_nin_and_nout; }

    PyUFuncGenericFunction *func() const { return m_func.get(); }

    data_handle_type *data() const { return m_data.get(); }

    char *types() const { return m_types.get(); }

    void set_name(const char *name) {
        for (int i = 0; i < m_ntypes; ++i) {
            static_cast<base_ufunc_data *>(m_data[i])->name = name;
        }
    }
};

PyObject *SpecFun_NewUFunc(SpecFun_UFunc func, int nout, const char *name, const char *doc) {
    static std::vector<SpecFun_UFunc> ufuncs;

    if (PyErr_Occurred()) {
        return nullptr;
    }

    SpecFun_UFunc &ufunc = ufuncs.emplace_back(std::move(func));
    ufunc.set_name(name);

    return PyUFunc_FromFuncAndData(
        ufunc.func(), ufunc.data(), ufunc.types(), ufunc.ntypes(), ufunc.nin_and_nout() - nout, nout, PyUFunc_None,
        name, doc, 0
    );
}

PyObject *SpecFun_NewUFunc(SpecFun_UFunc func, const char *name, const char *doc) {
    int nout = func.has_return();

    return SpecFun_NewUFunc(std::move(func), nout, name, doc);
}

PyObject *SpecFun_NewGUFunc(SpecFun_UFunc func, int nout, const char *name, const char *doc, const char *signature) {
    static std::vector<SpecFun_UFunc> ufuncs;

    if (PyErr_Occurred()) {
        return nullptr;
    }

    SpecFun_UFunc &ufunc = ufuncs.emplace_back(std::move(func));
    ufunc.set_name(name);

    return PyUFunc_FromFuncAndDataAndSignature(
        ufunc.func(), ufunc.data(), ufunc.types(), ufunc.ntypes(), ufunc.nin_and_nout() - nout, nout, PyUFunc_None,
        name, doc, 0, signature
    );
}

PyObject *SpecFun_NewGUFunc(SpecFun_UFunc func, const char *name, const char *doc, const char *signature) {
    int nout = func.has_return();

    return SpecFun_NewGUFunc(std::move(func), nout, name, doc, signature);
}
