/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HVector.h
 * @brief Vector structure for HiGHS
 */
#ifndef UTIL_HVECTOR_BASE_H_
#define UTIL_HVECTOR_BASE_H_

#include <vector>

#include "util/HighsType.h"

// using std::map;
using std::vector;

/**
 * @brief Class for the vector structure for HiGHS
 */
template <typename Real>
class HVectorBase {
 public:
  /**
   * @brief Initialise a vector
   */
  void setup(HighsInt size_  //!< Dimension of the vector to be initialised
  );

  /**
   * @brief Clear the vector - or just its scalars
   *
   */
  void clear();
  void clearScalars();

  HighsInt size;           //!< Dimension of the vector
  HighsInt count;          //!< Number of nonzeros
  vector<HighsInt> index;  //!< Packed indices of nonzeros
  vector<Real> array;      //!< Full-length array of values

  double synthetic_tick;  //!< Synthetic clock for operations with this vector

  // For update
  vector<HighsBool> cwork;  //!< Working buffer for UPDATE
  vector<HighsInt> iwork;   //!< integer working buffer for UPDATE
  HVectorBase<Real>* next;  //!< Allows vectors to be linked for PAMI

  /*
   * Zero values in Vector.array that exceed kHighsTiny in magnitude
   */
  void tight();

  /**
   * @brief Packing (if packFlag set): Pack values/indices in Vector.array into
   * packValue/Index
   *
   */
  void pack();

  /**
   * @brief Possibly determine the indices from scratch by passing
   * through the array
   *
   */
  void reIndex();

  bool packFlag = false;       //!< Flag to indicate whether to pack or not
  HighsInt packCount;          //!< Number of nonzeros packed
  vector<HighsInt> packIndex;  //!< Packed indices
  vector<Real> packValue;      //!< Packed values

  /**
   * @brief Copy from another HVector structure to this instance
   */
  template <typename FromReal>
  void copy(const HVectorBase<FromReal>*
                from  //!< Source of HVector structure to be copied
  );

  /**
   * @brief Compute the squared 2-norm of the vector
   */
  Real norm2() const;

  /**
   * @brief Add a multiple pivotX of *pivot into this vector,
   * maintaining indices of nonzeros but not tracking cancellation
   */
  template <typename RealPivX, typename RealPiv>
  void saxpy(const RealPivX pivotX,  //!< The multiple of *pivot to be added
             const HVectorBase<RealPiv>*
                 pivot  //!< The vector whose multiple is to be added
  );

  bool isEqual(const HVectorBase<Real>& v0);
};

#endif /* UTIL_HVECTOR_BASE_H_ */
