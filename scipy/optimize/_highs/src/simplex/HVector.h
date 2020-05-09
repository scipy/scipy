/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HVector.h
 * @brief Vector structure for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HVECTOR_H_
#define SIMPLEX_HVECTOR_H_

#include <map>
#include <vector>

using std::map;
using std::vector;

/**
 * @brief Class for the vector structure for HiGHS
 */
class HVector {
 public:
  /**
   * @brief Initialise a vector
   */
  void setup(int size_  //!< Dimension of the vector to be initialised
  );

  /**
   * @brief Clear the vector
   *
   */
  void clear();

  int size;              //!< Dimension of the vector
  int count;             //!< Number of nonzeros
  vector<int> index;     //!< Packed indices of nonzeros
  vector<double> array;  //!< Full-length array of values

  double syntheticTick;  //!< Synthetic clock for operations with this vector

  // For update
  vector<char> cwork;  //!< char working buffer for UPDATE
  vector<int> iwork;   //!< integer working buffer for UPDATE
  HVector* next;       //!< Allows vectors to be linked for PAMI

  /**
   * @brief Packing: Zero values in Vector.array which exceed HIGHS_CONST_TINY
   * in magnitude
   *
   */
  void tight();

  /**
   * @brief Packing (if packFlag set): Pack values/indices in Vector.array into
   * packValue/Index
   *
   */
  void pack();
  bool packFlag;             //!< Flag to indicate whether to pack or not
  int packCount;             //!< Number of nonzeros packed
  vector<int> packIndex;     //!< Packed indices
  vector<double> packValue;  //!< Packed values

  /**
   * @brief Copy from another HVector structure to this instance
   */
  void copy(const HVector* from  //!< Source of HVector structure to be copied
  );

  /**
   * @brief Compute the squared 2-norm of the vector
   */
  double norm2();

  /**
   * @brief Add a multiple pivotX of *pivot into this vector,
   * maintaining indices of nonzeros but not tracking cancellation
   */
  void saxpy(const double pivotX,  //!< The multiple of *pivot to be added
             const HVector* pivot  //!< The vector whose multiple is to be added
  );
};

typedef HVector* HVector_ptr;

#endif /* SIMPLEX_HVECTOR_H_ */
