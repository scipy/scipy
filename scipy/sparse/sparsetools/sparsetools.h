#ifndef SPARSETOOLS_H
#define SPARSETOOLS_H

#include <Python.h>
#include "numpy/ndarrayobject.h"

#include <array>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "bool_ops.h"
#include "complex_ops.h"
#include "util.h"

using SparsetoolsKernelFn = PY_LONG_LONG (*)(void **args);

enum class SparsetoolsDispatchKind : unsigned char {
    kNone = 0,       // Non-templated routine (1 kernel)
    kIndex = 1,      // Templated on <I>
    kIndexData = 2,  // Templated on <I, T>
};

struct SparsetoolsRoutine {
    const char *name;
    char ret_spec;
    SparsetoolsDispatchKind dispatch_kind;
    const char *arg_spec;
    const SparsetoolsKernelFn *kernels;
};

struct SparsetoolsRoutineSpan {
    const SparsetoolsRoutine *data;
    std::size_t size;
};

SparsetoolsRoutineSpan get_bsr_routines();
SparsetoolsRoutineSpan get_csr_routines();
SparsetoolsRoutineSpan get_csc_routines();
SparsetoolsRoutineSpan get_other_routines();

namespace sparsetools_detail {

// Sentinel types used only in unevaluated `decltype(&fn<ProbeI, ProbeT>)`
// expressions to deduce parameter/return specs at compile time without
// instantiating function bodies.
struct ProbeI {};
struct ProbeT {};

template <typename... Ts>
struct TypeList {};

// Template trick to get properly comma-separated type lists.
template <typename Sentinel, typename... Ts>
using DropLeadingType = TypeList<Ts...>;

#define SPTOOLS_APPEND_TYPE(T) , T
using IndexTypes =
    DropLeadingType<void SPTOOLS_FOR_EACH_INDEX_TYPE(SPTOOLS_APPEND_TYPE)>;
using DataTypes =
    DropLeadingType<void SPTOOLS_FOR_EACH_DATA_TYPE(SPTOOLS_APPEND_TYPE)>;
#undef SPTOOLS_APPEND_TYPE

template <typename List>
struct TypeListSize;

template <typename... Ts>
struct TypeListSize<TypeList<Ts...>> {
    static constexpr std::size_t value = sizeof...(Ts);
};

template <typename T, typename List>
struct TypeListIndex;

template <typename T, typename... Ts>
struct TypeListIndex<T, TypeList<Ts...>> {
    static constexpr int find() {
        constexpr bool matches[] = {std::is_same_v<T, Ts>...};
        int found = -1;
        int count = 0;
        for (std::size_t i = 0; i < sizeof...(Ts); ++i) {
            if (matches[i]) {
                found = static_cast<int>(i);
                ++count;
            }
        }
        return count == 1 ? found : -1;
    }
    static constexpr int value = find();
    static_assert(value >= 0, "Type must appear exactly once in TypeList");
};

inline constexpr std::size_t kNumIndexTypes = TypeListSize<IndexTypes>::value;
inline constexpr std::size_t kNumDataTypes = TypeListSize<DataTypes>::value;

template <typename I>
inline constexpr int kIndexTypeIndex = TypeListIndex<I, IndexTypes>::value;

template <typename T>
inline constexpr int kDataTypeIndex = TypeListIndex<T, DataTypes>::value;

template <typename R>
constexpr char deduce_ret_spec() {
    static_assert(std::is_void_v<R> || std::is_same_v<R, ProbeI> ||
                      std::is_integral_v<R>,
                  "Unsupported sparsetools return type");
    return std::is_void_v<R> ? 'v' : 'i';
}

struct SpecFrag {
    char chars[2];
    std::size_t len;
};

template <typename Arg>
constexpr SpecFrag deduce_arg_frag() {
    if constexpr (std::is_same_v<Arg, ProbeI>) {
        return {{'i', '\0'}, 1};
    } else if constexpr (std::is_same_v<Arg, npy_int64>) {
        return {{'l', '\0'}, 1};
    } else if constexpr (std::is_same_v<Arg, int>) {
        return {{'i', '\0'}, 1};
    } else if constexpr (std::is_same_v<Arg, const ProbeI *>) {
        return {{'I', '\0'}, 1};
    } else if constexpr (std::is_same_v<Arg, ProbeI *>) {
        return {{'*', 'I'}, 2};
    } else if constexpr (std::is_same_v<Arg, const ProbeT *>) {
        return {{'T', '\0'}, 1};
    } else if constexpr (std::is_same_v<Arg, ProbeT *>) {
        return {{'*', 'T'}, 2};
    } else if constexpr (std::is_same_v<Arg, const npy_bool_wrapper *>) {
        return {{'B', '\0'}, 1};
    } else if constexpr (std::is_same_v<Arg, npy_bool_wrapper *>) {
        return {{'*', 'B'}, 2};
    } else if constexpr (std::is_same_v<Arg, std::vector<ProbeI> *>) {
        return {{'*', 'V'}, 2};
    } else if constexpr (std::is_same_v<Arg, std::vector<ProbeT> *>) {
        return {{'*', 'W'}, 2};
    } else {
        static_assert(sizeof(Arg) == 0,
                      "Unsupported sparsetools parameter type");
        return {{'\0', '\0'}, 0};
    }
}

template <typename FnPtr>
struct SignatureTraits;

template <typename R, typename... Args>
struct SignatureTraits<R (*)(Args...)> {
    static constexpr char ret_spec = deduce_ret_spec<R>();
    static constexpr std::size_t spec_len =
        (std::size_t{0} + ... + deduce_arg_frag<Args>().len);

    static constexpr std::array<char, spec_len + 1> make_arg_spec() {
        std::array<char, spec_len + 1> out{};
        const SpecFrag frags[sizeof...(Args) == 0 ? 1 : sizeof...(Args)] = {
            deduce_arg_frag<Args>()...};
        std::size_t pos = 0;
        for (std::size_t i = 0; i < sizeof...(Args); ++i) {
            out[pos++] = frags[i].chars[0];
            if (frags[i].len == 2) {
                out[pos++] = frags[i].chars[1];
            }
        }
        out[pos] = '\0';
        return out;
    }

    static constexpr std::array<char, spec_len + 1> arg_spec = make_arg_spec();
};

template <typename Arg>
inline Arg extract_arg(void *p) {
    if constexpr (std::is_pointer_v<Arg>) {
        return static_cast<Arg>(p);
    } else {
        static_assert(std::is_integral_v<Arg>,
                      "Expected integral scalar parameter");
        return static_cast<Arg>(*static_cast<const npy_int64 *>(p));
    }
}

template <typename R, typename... Args, std::size_t... Idx>
inline PY_LONG_LONG call_with_void_args(
    R (*fn)(Args...), void **a, std::index_sequence<Idx...>) {
    if constexpr (std::is_void_v<R>) {
        fn(extract_arg<Args>(a[Idx])...);
        return 0;
    } else {
        return static_cast<PY_LONG_LONG>(fn(extract_arg<Args>(a[Idx])...));
    }
}

template <typename R, typename... Args>
inline PY_LONG_LONG dispatch_fn_ptr(R (*fn)(Args...), void **a) {
    return call_with_void_args(fn, a, std::index_sequence_for<Args...>{});
}

template <auto Fn>
PY_LONG_LONG invoke_kernel(void **a) {
    return dispatch_fn_ptr(Fn, a);
}

template <template <typename> class Binder, typename... Is>
constexpr std::array<SparsetoolsKernelFn, sizeof...(Is)> make_i_kernels_impl(
    TypeList<Is...>) {
    return {&invoke_kernel<Binder<Is>::fn_ptr>...};
}

template <template <typename> class Binder>
constexpr auto make_i_kernels() {
    return make_i_kernels_impl<Binder>(IndexTypes{});
}

template <template <typename, typename> class Binder, typename I, typename... Ts>
constexpr std::array<SparsetoolsKernelFn, sizeof...(Ts)> make_it_row(
    TypeList<Ts...>) {
    return {&invoke_kernel<Binder<I, Ts>::fn_ptr>...};
}

template <template <typename, typename> class Binder, typename... Is>
constexpr std::array<SparsetoolsKernelFn, kNumIndexTypes * kNumDataTypes>
make_it_kernels_impl(TypeList<Is...>) {
    constexpr std::array<std::array<SparsetoolsKernelFn, kNumDataTypes>,
                         sizeof...(Is)>
        rows = {make_it_row<Binder, Is>(DataTypes{})...};
    std::array<SparsetoolsKernelFn, kNumIndexTypes * kNumDataTypes> flat{};
    for (std::size_t i = 0; i < sizeof...(Is); ++i) {
        for (std::size_t t = 0; t < kNumDataTypes; ++t) {
            flat[i * kNumDataTypes + t] = rows[i][t];
        }
    }
    return flat;
}

template <template <typename, typename> class Binder>
constexpr auto make_it_kernels() {
    return make_it_kernels_impl<Binder>(IndexTypes{});
}

}  // namespace sparsetools_detail

#define SPTOOLS_DECL_TABLE_0(fn)                                        \
  constexpr std::array<SparsetoolsKernelFn, 1> kKernels_##fn = {        \
      &::sparsetools_detail::invoke_kernel<&fn>,                        \
  };

#define SPTOOLS_DECL_TABLE_I(fn)                                        \
  template <typename I>                                                 \
  struct Binder_##fn {                                                  \
      static constexpr auto fn_ptr = &fn<I>;                            \
  };                                                                    \
  constexpr auto kKernels_##fn =                                        \
      ::sparsetools_detail::make_i_kernels<Binder_##fn>();

#define SPTOOLS_DECL_TABLE_IT(fn)                                       \
  template <typename I, typename T>                                     \
  struct Binder_##fn {                                                  \
      static constexpr auto fn_ptr = &fn<I, T>;                         \
  };                                                                    \
  constexpr auto kKernels_##fn =                                        \
      ::sparsetools_detail::make_it_kernels<Binder_##fn>();

#define SPTOOLS_ENTRY_0(fn)                                             \
  {#fn,                                                                 \
   ::sparsetools_detail::SignatureTraits<decltype(&fn)>::ret_spec,      \
   SparsetoolsDispatchKind::kNone,                                      \
   ::sparsetools_detail::SignatureTraits<decltype(&fn)>::arg_spec.data(), \
   kKernels_##fn.data()},

#define SPTOOLS_ENTRY_I(fn)                                             \
  {#fn,                                                                 \
   ::sparsetools_detail::SignatureTraits<                               \
       decltype(&fn<::sparsetools_detail::ProbeI>)>::ret_spec,          \
   SparsetoolsDispatchKind::kIndex,                                     \
   ::sparsetools_detail::SignatureTraits<                               \
       decltype(&fn<::sparsetools_detail::ProbeI>)>::arg_spec.data(),   \
   kKernels_##fn.data()},

#define SPTOOLS_ENTRY_IT(fn)                                            \
  {#fn,                                                                 \
   ::sparsetools_detail::SignatureTraits<                               \
       decltype(&fn<::sparsetools_detail::ProbeI,                       \
                    ::sparsetools_detail::ProbeT>)>::ret_spec,          \
   SparsetoolsDispatchKind::kIndexData,                                 \
   ::sparsetools_detail::SignatureTraits<                               \
       decltype(&fn<::sparsetools_detail::ProbeI,                       \
                    ::sparsetools_detail::ProbeT>)>::arg_spec.data(),   \
   kKernels_##fn.data()},

#define SPTOOLS_DEFINE_UNIT(unit_name, ROUTINE_LIST)                    \
  namespace {                                                           \
  ROUTINE_LIST(SPTOOLS_DECL_TABLE_0, SPTOOLS_DECL_TABLE_I,              \
               SPTOOLS_DECL_TABLE_IT)                                   \
  constexpr SparsetoolsRoutine kRoutines_##unit_name[] = {              \
      ROUTINE_LIST(SPTOOLS_ENTRY_0, SPTOOLS_ENTRY_I, SPTOOLS_ENTRY_IT)  \
  };                                                                    \
  }  /* namespace */                                                    \
  SparsetoolsRoutineSpan get_##unit_name##_routines() {                 \
      return {kRoutines_##unit_name,                                    \
              std::size(kRoutines_##unit_name)};                        \
  }

#endif
