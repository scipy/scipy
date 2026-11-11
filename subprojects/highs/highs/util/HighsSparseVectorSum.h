/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SPARSE_VECTOR_SUM_H_
#define HIGHS_SPARSE_VECTOR_SUM_H_

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

#include "util/HighsCDouble.h"
#include "util/HighsType.h"

class HighsSparseVectorSum {
 public:
  std::vector<HighsCDouble> values;
  std::vector<HighsInt> nonzeroinds;
  HighsSparseVectorSum() = default;

  HighsSparseVectorSum(HighsInt dimension) { setDimension(dimension); }

  void setDimension(HighsInt dimension) {
    values.resize(dimension);
    nonzeroinds.reserve(dimension);
  }

  void add(HighsInt index, double value) {
    assert(index >= 0 && index < (HighsInt)values.size());
    if (values[index] != 0.0) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroinds.push_back(index);
    }

    if (values[index] == 0.0)
      values[index] = (std::numeric_limits<double>::min)();
  }

  void add(HighsInt index, HighsCDouble value) {
    if (values[index] != 0.0) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroinds.push_back(index);
    }

    if (values[index] == 0.0)
      values[index] = (std::numeric_limits<double>::min)();
  }

  const std::vector<HighsInt>& getNonzeros() const { return nonzeroinds; }

  template <typename T = double>
  T getValue(HighsInt index) const {
    return static_cast<T>(values[index]);
  }

  void clear() {
    if (10 * nonzeroinds.size() < 3 * values.size())
      for (HighsInt i : nonzeroinds) values[i] = 0.0;
    else
      values.assign(values.size(), false);

    nonzeroinds.clear();
  }

  template <typename Pred>
  HighsInt partition(Pred&& pred) {
    return std::partition(nonzeroinds.begin(), nonzeroinds.end(), pred) -
           nonzeroinds.begin();
  }

  template <typename IsZero>
  void cleanup(IsZero&& isZero) {
    HighsInt numNz = static_cast<HighsInt>(nonzeroinds.size());

    for (HighsInt i = numNz - 1; i >= 0; --i) {
      HighsInt pos = nonzeroinds[i];
      double val = double(values[pos]);

      if (isZero(pos, val)) {
        values[pos] = 0.0;
        --numNz;
        std::swap(nonzeroinds[numNz], nonzeroinds[i]);
      }
    }

    nonzeroinds.resize(numNz);
  }
};

#endif
