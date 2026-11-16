/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasApiBinding.h
 * @brief Provides useful classes for handling external dependencies
 */
#ifndef HIGHS_EXTRAS_API_BINDING_H_
#define HIGHS_EXTRAS_API_BINDING_H_

#include <stddef.h>

#include <tuple>


// provide metadata info for each feature
struct HighsExtrasFeatureInfo {
  HighsExtrasFeatureInfo(const char* provider_ = nullptr,
                         const char* version_ = nullptr,
                         const char* license_ = nullptr, bool enabled_ = false)
      : provider(provider_),
        version(version_),
        license(license_),
        enabled(enabled_) {}

  const char* provider;
  const char* version;
  const char* license;
  const bool enabled;
};

namespace HighsExtras {

template <class... Features>
struct require {};

// compile-time count of nested features in a list,
// e.g., feature_count<require<amd, blas>, rcm>::value == 3
template <class T>
struct feature_count : std::integral_constant<size_t, 1> {};

template <class... Fs>
struct sum_feature_counts;

template <>
struct sum_feature_counts<> : std::integral_constant<size_t, 0> {};

template <class F, class... Rest>
struct sum_feature_counts<F, Rest...>
    : std::integral_constant<size_t, feature_count<F>::value +
                                         sum_feature_counts<Rest...>::value> {};

template <class... Features>
struct feature_count<require<Features...>> : sum_feature_counts<Features...> {};

template <class Methods>
struct feature_api;

// convenience wrapper to access the HighsExtrasApi storage
template <class Family>
struct wrapper_storage {
  template <class Methods>
  static feature_api<Methods>& getApi();

  static const HighsExtrasFeatureInfo* getInfo() { return nullptr; };
};

template <class FeatureFamily, int FeatureIndex>
struct feature_base {
  using family = FeatureFamily;

  static const HighsExtrasFeatureInfo* getInfo() {
    const HighsExtrasFeatureInfo* ptr = wrapper_storage<family>::getInfo();
    return ptr ? ptr + FeatureIndex : nullptr;
  }
};

// storage for the function pointer, given a decltype(&...)
template <class Desc>
struct method_storage {
  Desc value;
  method_storage() : value(nullptr) {}
};

// builds a struct of function pointers, given tuple<decltype(&...), ...>
template <class... Desc>
struct feature_api<std::tuple<Desc...>> : method_storage<Desc>... {
  using methods_type = std::tuple<Desc...>;

  // access method by index at compile-time, e.g., api->method<0>(...)
  template <std::size_t Index>
  typename std::tuple_element<Index, methods_type>::type& method() {
    using desc_type = typename std::tuple_element<Index, methods_type>::type;
    return static_cast<method_storage<desc_type>&>(*this).value;
  }
};

// access function pointer by index, e.g., api::fn<0>()(...)
template <class Family, class Methods>
struct feature_wrapper {
  template <std::size_t Index>
  static typename std::tuple_element<Index, Methods>::type& fn() {
    return wrapper_storage<Family>::template getApi<Methods>()
        .template method<Index>();
  }
};

template <class Methods, class Tuple, std::size_t Index,
          std::size_t Count = std::tuple_size<Methods>::value>
struct bind_from_tuple_impl {
  static void apply(feature_api<Methods>& api, const Tuple& funcs) {
    api.template method<Index>() = std::get<Index>(funcs);
    bind_from_tuple_impl<Methods, Tuple, Index + 1, Count>::apply(api, funcs);
  }
};

template <class Methods, class Tuple, std::size_t Count>
struct bind_from_tuple_impl<Methods, Tuple, Count, Count> {
  static void apply(feature_api<Methods>&, const Tuple&) {}
};

}  // namespace HighsExtras

#endif  // HIGHS_EXTRAS_API_BINDING_H_