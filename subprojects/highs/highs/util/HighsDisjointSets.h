/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_DISJOINT_SETS_H_
#define HIGHS_UTIL_DISJOINT_SETS_H_

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include "util/HighsType.h"

template <bool kMinimalRepresentative = false>
class HighsDisjointSets {
  std::vector<HighsInt> sizes;
  std::vector<HighsInt> sets;
  std::vector<HighsInt> linkCompressionStack;

 public:
  HighsDisjointSets() = default;
  HighsDisjointSets(HighsInt numItems) { reset(numItems); }

  void reset(HighsInt numItems) {
    sizes.assign(numItems, 1);
    sets.resize(numItems);
    std::iota(sets.begin(), sets.end(), 0);
  }

  HighsInt getSet(HighsInt item) {
    assert(item >= 0 && item < (HighsInt)sets.size());
    HighsInt repr = sets[item];
    assert(repr >= 0 && repr < (HighsInt)sets.size());

    if (repr != sets[repr]) {
      do {
        linkCompressionStack.push_back(item);
        item = repr;
        repr = sets[repr];
      } while (repr != sets[repr]);

      do {
        HighsInt i = linkCompressionStack.back();
        linkCompressionStack.pop_back();
        sets[i] = repr;
      } while (!linkCompressionStack.empty());

      sets[item] = repr;
    }

    return repr;
  }

  HighsInt getSetSize(HighsInt set) const {
    assert(sets[set] == set);
    return sizes[set];
  }

  void merge(HighsInt item1, HighsInt item2) {
    assert(item1 >= 0 && item1 < (HighsInt)sets.size());
    assert(item2 >= 0 && item2 < (HighsInt)sets.size());

    HighsInt repr1 = getSet(item1);
    assert(sets[repr1] == repr1);
    assert(repr1 >= 0 && repr1 < (HighsInt)sets.size());

    HighsInt repr2 = getSet(item2);
    assert(sets[repr2] == repr2);
    assert(repr2 >= 0 && repr2 < (HighsInt)sets.size());
    assert(sizes.size() == sets.size());

    if (repr1 == repr2) return;

    if (kMinimalRepresentative) {
      if (repr2 > repr1) {
        sets[repr2] = repr1;
        sizes[repr1] += sizes[repr2];
      } else {
        sets[repr1] = repr2;
        sizes[repr2] += sizes[repr1];
      }
    } else {
      if (sizes[repr1] > sizes[repr2]) {
        sets[repr2] = repr1;
        sizes[repr1] += sizes[repr2];
      } else {
        sets[repr1] = repr2;
        sizes[repr2] += sizes[repr1];
      }
    }
  }
};

#endif
