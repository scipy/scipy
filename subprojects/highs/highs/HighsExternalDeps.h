/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExternalDeps.h
 * @brief Defines the set of external features available
 */
#ifndef HIGHS_EXTERNAL_DEPS_H_
#define HIGHS_EXTERNAL_DEPS_H_

#include "HighsExtrasApi.h"

namespace HighsExtras {
struct highs_family {
  static const HighsExtrasFeatureInfo* getInfo();
};

template <>
inline const HighsExtrasFeatureInfo* wrapper_storage<highs_family>::getInfo() {
  return highs_family::getInfo();
}

template <int Index>
struct highs_feature : feature_base<highs_family, Index> {
  static const char* name() {
    return std::get<Index>(std::make_tuple("pdqsort", "zstr", "zlib", "cuda"));
  }
};

using pdqsort = highs_feature<0>;
using zstr = highs_feature<1>;
using zlib = highs_feature<2>;
using cuda = highs_feature<3>;

// define feature sets
using all = require<extrasAll, pdqsort, zstr, zlib, cuda>;

}  // namespace HighsExtras

#endif  // HIGHS_EXTERNAL_DEPS_H_