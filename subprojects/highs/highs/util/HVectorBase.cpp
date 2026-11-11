/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HVectorBase.cpp
 * @brief
 */

#include "util/HVectorBase.h"

#include <cassert>
#include <cmath>

#include "lp_data/HConst.h"
#include "stdio.h"  //Just for temporary printf
#include "util/HighsCDouble.h"

template <typename Real>
void HVectorBase<Real>::setup(HighsInt size_) {
  /*
   * Initialise an HVector instance
   */
  size = size_;
  count = 0;
  index.resize(size);
  array.assign(size, Real{0});
  cwork.assign(size + 6400, false);  // MAX invert
  iwork.assign(size * 4, 0);

  packCount = 0;
  packIndex.resize(size);
  packValue.resize(size);

  // Initialise three values that are initialised in clear(), but
  // weren't originally initialised in setup(). Probably doesn't
  // matter, since clear() is usually called before a vector is
  // (re-)used.
  packFlag = false;
  synthetic_tick = 0;
  next = 0;
}

template <typename Real>
void HVectorBase<Real>::clear() {
  /*
   * Clear an HVector instance
   */
  // Standard HVector to clear
  HighsInt dense_clear = count < 0 || count > size * 0.3;
  if (dense_clear) {
    // Treat the array as full if there are no indices or too many
    // indices
    array.assign(size, Real{0});
  } else {
    // Zero according to the indices of (possible) nonzeros
    for (HighsInt i = 0; i < count; i++) {
      array[index[i]] = 0;
    }
  }
  this->clearScalars();
}

template <typename Real>
void HVectorBase<Real>::clearScalars() {
  /*
   * Clear scalars in an HVector instance
   */
  // Reset the flag to indicate when to pack
  this->packFlag = false;
  // Zero the number of stored indices
  this->count = 0;
  // Zero the synthetic clock for operations with this vector
  this->synthetic_tick = 0;
  // Initialise the next value
  this->next = 0;
}

template <typename Real>
void HVectorBase<Real>::tight() {
  /*
   * Zero values in Vector.array that do not exceed kHighsTiny in
   * magnitude, maintaining index if it is well defined
   */
  HighsInt totalCount = 0;
  using std::abs;
  if (count < 0) {
    for (auto& val : array)
      if (abs(val) < kHighsTiny) val = 0;
  } else {
    for (HighsInt i = 0; i < count; i++) {
      const HighsInt my_index = index[i];
      const Real& value = array[my_index];
      if (abs(value) >= kHighsTiny) {
        index[totalCount++] = my_index;
      } else {
        array[my_index] = Real{0};
      }
    }
    count = totalCount;
  }
}

template <typename Real>
void HVectorBase<Real>::pack() {
  /*
   * Packing (if packFlag set): Pack values/indices in Vector.array
   * into packValue/Index
   */
  if (!packFlag) return;
  packFlag = false;
  packCount = 0;
  for (HighsInt i = 0; i < count; i++) {
    const HighsInt ipack = index[i];
    packIndex[packCount] = ipack;
    packValue[packCount] = array[ipack];
    packCount++;
  }
}

template <typename Real>
void HVectorBase<Real>::reIndex() {
  /*
   * Possibly determine the indices from scratch by passing through
   * the array
   */
  // Don't do it if there are relatively few nonzeros
  if (count >= 0 && count <= size * 0.1) return;
  count = 0;
  for (HighsInt i = 0; i < size; i++)
    if ((double)array[i]) index[count++] = i;
}

template <typename Real>
template <typename FromReal>
void HVectorBase<Real>::copy(const HVectorBase<FromReal>* from) {
  /*
   * Copy from another HVector structure to this instance
   * The real type of "from" does not need to be the same, but must be
   * convertible to this HVector's real type.
   */
  clear();
  synthetic_tick = from->synthetic_tick;
  const HighsInt fromCount = count = from->count;
  const HighsInt* fromIndex = &from->index[0];
  const FromReal* fromArray = &from->array[0];
  for (HighsInt i = 0; i < fromCount; i++) {
    const HighsInt iFrom = fromIndex[i];
    const FromReal xFrom = fromArray[iFrom];
    index[i] = iFrom;
    array[iFrom] = Real(xFrom);
  }
}
template <typename Real>
Real HVectorBase<Real>::norm2() const {
  /*
   * Compute the squared 2-norm of the vector
   */
  const HighsInt workCount = count;
  const HighsInt* workIndex = &index[0];
  const Real* workArray = &array[0];

  Real result = Real{0};
  for (HighsInt i = 0; i < workCount; i++) {
    Real value = workArray[workIndex[i]];
    result += value * value;
  }
  return result;
}

template <typename Real>
template <typename RealPivX, typename RealPiv>
void HVectorBase<Real>::saxpy(const RealPivX pivotX,
                              const HVectorBase<RealPiv>* pivot) {
  /*
   * Add a multiple pivotX of *pivot into this vector, maintaining
   * indices of nonzeros but not tracking cancellation.
   * The real types may all be different but must mix in operations and be
   * convertible to this HVector's real type.
   */
  HighsInt workCount = count;
  HighsInt* workIndex = &index[0];
  Real* workArray = &array[0];

  const HighsInt pivotCount = pivot->count;
  const HighsInt* pivotIndex = &pivot->index[0];
  const RealPiv* pivotArray = &pivot->array[0];

  using std::abs;
  for (HighsInt k = 0; k < pivotCount; k++) {
    const HighsInt iRow = pivotIndex[k];
    const Real x0 = workArray[iRow];
    const Real x1 = Real(x0 + pivotX * pivotArray[iRow]);
    if (x0 == Real{0}) workIndex[workCount++] = iRow;
    workArray[iRow] = (abs(x1) < kHighsTiny) ? kHighsZero : x1;
  }
  count = workCount;
}

template <typename Real>
bool HVectorBase<Real>::isEqual(const HVectorBase<Real>& v0) {
  if (this->size != v0.size) return false;
  if (this->count != v0.count) return false;
  if (this->index != v0.index) return false;
  if (this->array != v0.array) return false;
  //  if (this->index.size() != v0.index.size()) return false;
  //  for (HighsInt el = 0; el < (HighsInt)this->index.size(); el++)
  //    if (this->index[el] != v0.index[el]) return false;
  if (this->synthetic_tick != v0.synthetic_tick) return false;
  return true;
}

// explicitly instantiate HVectorBase<T> T=double
template class HVectorBase<double>;

// explicitly instantiate template member function "copy" for
// HVectorBase<double> with the type of the copied HVectorBase being either
// double or HighsCDouble
template void HVectorBase<double>::copy(const HVectorBase<double>*);
template void HVectorBase<double>::copy(const HVectorBase<HighsCDouble>*);

// explicitly instantiate template member function "saxpy" for
// HVectorBase<double> with all four combinations of types for pivot and pivotX:
// (double double), (double HighsCDouble), (HighsCDouble double), (HighsCDouble
// HighsCDouble)
template void HVectorBase<double>::saxpy(const double,
                                         const HVectorBase<double>*);
template void HVectorBase<double>::saxpy(const double,
                                         const HVectorBase<HighsCDouble>*);
template void HVectorBase<double>::saxpy(const HighsCDouble,
                                         const HVectorBase<double>*);
template void HVectorBase<double>::saxpy(const HighsCDouble,
                                         const HVectorBase<HighsCDouble>*);

// explicitly instantiate HVectorBase<T> T=HighsCDouble
template class HVectorBase<HighsCDouble>;

// explicitly instantiate template member function "copy" for
// HVectorBase<HighsCDouble> with the type of the copied HVectorBase being
// either double or HighsCDouble
template void HVectorBase<HighsCDouble>::copy(const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::copy(const HVectorBase<HighsCDouble>*);

// explicitly instantiate template member function "saxpy" for
// HVectorBase<HighsCDouble> with all four combinations of types for pivot and
// pivotX: (double double), (double HighsCDouble), (HighsCDouble double),
// (HighsCDouble HighsCDouble)
template void HVectorBase<HighsCDouble>::saxpy(const double,
                                               const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::saxpy(
    const double, const HVectorBase<HighsCDouble>*);
template void HVectorBase<HighsCDouble>::saxpy(const HighsCDouble,
                                               const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::saxpy(
    const HighsCDouble, const HVectorBase<HighsCDouble>*);

#if 0
// Todo: Additionally we could add the two specializations to allow pivotX as
// int so that integer constants can be used, but I think we can avoid the
// additional code bloat and I changed the single place where this is used to
// call with 1.0 instead of 1
template void HVectorBase<HighsCDouble>::saxpy(const int,
                                               const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::saxpy(
    const int, const HVectorBase<HighsCDouble>*);
template void HVectorBase<double>::saxpy(const int, const HVectorBase<double>*);
template void HVectorBase<double>::saxpy(const int,
                                         const HVectorBase<HighsCDouble>*);
#endif
