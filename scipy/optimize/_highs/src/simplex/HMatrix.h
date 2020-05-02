/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HMatrix.h
 * @brief Column-wise and partitioned row-wise representaiton of the constraint
 * matrix for the simplex solver
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HMATRIX_H_
#define SIMPLEX_HMATRIX_H_

#include <vector>

#include "HConfig.h"

class HVector;

/**
 * @brief Column-wise and partitioned row-wise representation of the
 * constraint matrix for the simplex solver
 */
class HMatrix {
 public:
  /**
   * @brief For a given basis, sets up the column-wise and partitioned
   * row-wise representation of the constraint matrix for the simplex
   * solver
   */
  void setup(
      int numCol,              //!< Number of columns in the source matrix
      int numRow,              //!< Number of rows in the source matrix
      const int* Astart,       //!< Pointer to the starts of the source matrix
      const int* Aindex,       //!< Pointer to the indices of the source matrix
      const double* Avalue,    //!< Pointer to the values of the source matrix
      const int* nonbasicFlag  //!< Pointer to the flags indicating which
                               //!< columns are basic and nonbasic
  );
  /**
   * @brief For a logical basis, sets up the column-wise and
   * partitioned row-wise representation of the constraint matrix for
   * the simplex solver
   */
  void setup_lgBs(
      int numCol,           //!< Number of columns in the source matrix
      int numRow,           //!< Number of rows in the source matrix
      const int* Astart,    //!< Pointer to the starts of the source matrix
      const int* Aindex,    //!< Pointer to the indices of the source matrix
      const double* Avalue  //!< Pointer to the values of the source matrix
  );
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T A \f$
   * column-wise, maintaining indices of nonzeros in the result
   */
  void priceByColumn(
      HVector& row_ap,               //!< Vector \f$ \mathbf{y}\f$
      const HVector& row_ep) const;  //!< Vector \f$ \mathbf{x}\f$
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise, maintaining indices of nonzeros in result
   */
  void priceByRowSparseResult(
      HVector& row_ap,               //!< Vector \f$ \mathbf{y}\f$
      const HVector& row_ep) const;  //!< Vector \f$ \mathbf{x}\f$
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise from a given row, maintaining indices of nonzeros in
   * result, with possible switches to standard row-wise PRICE either
   * immediately based on historical density of the PRICE results, or
   * during PRICE if there is too much fill-in
   */
  void priceByRowSparseResultWithSwitch(
      HVector& row_ap,            //!< Vector \f$ \mathbf{y}\f$
      const HVector& row_ep,      //!< Vector \f$ \mathbf{x}\f$
      double historical_density,  //!< Historical density of PRICE results to be
                                  //!< used
      int from_i,                 //!< Index of row to work from
      double switch_density) const;  //!< Density for switch to not maintaining
                                     //!< indices of nonzeros
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise from a given row, not maintaining indices of nonzeros in
   * result
   */
  void priceByRowDenseResult(
      HVector& row_ap,        //!< Vector \f$ \mathbf{y}\f$
      const HVector& row_ep,  //!< Vector \f$ \mathbf{x}\f$
      int from_i) const;      //!< Index of row to work from
  /**
   * @brief Remove indices of zeros from vector \f$ \mathbf{y}\f$ created by
   * cancellation in PRICE
   */
  void priceByRowSparseResultRemoveCancellation(
      HVector& row_ap) const;  //!< Vector \f$ \mathbf{y} \f$
  /**
   * @brief Update the partitioned row-wise representation according
   * to columns coming in and out of the set of indices of basic
   * variables
   */
  void update(
      int columnIn,  //!< Column entering the set of indices of basic variables
      int columnOut  //!< Column leaving the set of indices of basic variables
  );
  /**
   * @brief Compute the dot product between a vector and particular
   * column of the constraint matrix: \f$ \mathbf{x}^T\mathbf{a}_i \f$
   */
  double compute_dot(HVector& vector,  //!< Vector \f$ \mathbf{x} \f$
                     int iCol) const;  //!< Index  \f$ i\f$ of column
  /**
   * @brief Add into a vector, a multiple of a particular column of
   * the constraint matrix \f$ \mathbf{x} := \mathbf{x} + \mu \mathbf{a}_i \f$
   */
  void collect_aj(HVector& vector,           //!< Vector \f$ \mathbf{x} \f$
                  int iCol,                  //!< Index  \f$ i\f$ of column
                  double multiplier) const;  //!< Multiplier \f$ \mu \f$

  /**
   * @brief Get the pointer to the starts of the column-wise matrix
   */
  const int* getAstart() const { return &Astart[0]; }

  /**
   * @brief Get the pointer to the indices of the column-wise matrix
   */
  const int* getAindex() const { return &Aindex[0]; }

  /**
   * @brief Get the pointer to the values of the column-wise matrix
   */
  const double* getAvalue() const { return &Avalue[0]; }

#ifdef HiGHSDEV
  bool setup_ok(const int* nonbasicFlag);
  bool price_er_ck(HVector& row_ap, const HVector& row_ep) const;
  bool price_er_ck_core(HVector& row_ap, const HVector& row_ep) const;
#endif

  /**
   * @brief Density of result at which it is not worth maintaing
   * indices of nonzeros
   */
  const double hyperPRICE = 0.10;

 private:
  int numCol;
  int numRow;
  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;

  std::vector<int> ARstart;
  std::vector<int> AR_Nend;
  std::vector<int> ARindex;
  std::vector<double> ARvalue;
};

#endif /* SIMPLEX_HMATRIX_H_ */
