/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasApi.cpp
 * @brief C-style API implementation for the HiGHS library dependencies
 */

#include "HighsExtrasApi.h"

#include "HighsExtrasExternalDeps.h"

using namespace HighsExtras;

namespace {

template <class Methods, class... Fn>
void bind_api(HighsExtrasApi* api, const std::tuple<Fn...>& funcs) {
  static_assert(std::tuple_size<Methods>::value == sizeof...(Fn),
                "bind_from_tuple requires one function per API entry");
  bind_from_tuple_impl<Methods, std::tuple<Fn...>, 0>::apply(
      api->template as<Methods>(), funcs);
}

// Get feature name by index at runtime via template recursion
template <size_t N>
struct feature_name {
  static const char* get(size_t index) {
    return (index == N - 1) ? extras_feature<N - 1>::name()
                            : feature_name<N - 1>::get(index);
  }
};

template <>
struct feature_name<0> {
  static const char* get(size_t) { return nullptr; }
};
}  // namespace

extern "C" {

HIGHS_EXTRAS_API const char* HighsExtras_getVersion(void) {
  return HIGHS_EXTRAS_VERSION;
}

HIGHS_EXTRAS_API size_t HighsExtras_getFeatureCount(void) {
  return feature_count<extrasAll>::value;
}

HIGHS_EXTRAS_API const char* HighsExtras_getFeatureName(size_t index) {
  if (index >= feature_count<extrasAll>::value) return nullptr;
  return feature_name<feature_count<extrasAll>::value>::get(index);
}

HIGHS_EXTRAS_API const HighsExtrasFeatureInfo* HighsExtras_getFeatureInfo() {
  return extras_family::getInfo();
}

HIGHS_EXTRAS_API bool HighsExtras_getApi(HighsExtrasApi* api) {
  if (!api) return false;

  *api = HighsExtrasApi{};

  // set function pointers for each feature API, if available
  // i.e., use direct descriptor to set the storage value
#ifdef HIPO

  bind_api<amd_methods>(api,
                        std::make_tuple(&Highs_amd_defaults, &Highs_amd_order));

  bind_api<blas_methods>(
      api,
      std::make_tuple(&cblas_daxpy, &cblas_dcopy, &cblas_dscal, &cblas_dswap,
                      &cblas_dgemv, &cblas_dtpsv, &cblas_dtrsv, &cblas_dger,
                      &cblas_dgemm, &cblas_dsyrk, &cblas_dtrsm,
                      &highs_openblas_set_num_threads));

  bind_api<metis_methods>(api, std::make_tuple(&Highs_METIS_SetDefaultOptions,
                                               &Highs_METIS_NodeND));

  bind_api<rcm_methods>(api, std::make_tuple(&Highs_genrcm));

#endif

  return true;
}

}  // extern "C"
