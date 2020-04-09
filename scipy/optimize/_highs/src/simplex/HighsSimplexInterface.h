/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexInterface.h
 * @brief Return or report data from simplex solves and interface that data with
 * changes to models
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSIMPLEXINTERFACE_H_
#define SIMPLEX_HIGHSSIMPLEXINTERFACE_H_

#include <vector>

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HSimplex.h"
#include "simplex/HVector.h"

// class HFactor;

/**
 * @brief Return or report data from simplex solves and interface that data with
 * changes to models
 */
class HighsSimplexInterface {
 public:
  HighsSimplexInterface(HighsModelObject& highs_model_object)
      : highs_model_object(highs_model_object) {}

  HighsModelObject& highs_model_object;

  /**
   * @brief Add a contiguous set of columns to the model data---making them
   * nonbasic
   */
  HighsStatus addCols(int XnumCol, const double* XcolCost,
                      const double* XcolLower, const double* XcolUpper,
                      int XnumNZ, const int* XAstart, const int* XAindex,
                      const double* XAvalue);

  HighsStatus deleteCols(int from_col, int to_col);

  HighsStatus deleteCols(int num_set_entries, const int* col_set);

  HighsStatus deleteCols(int* col_mask);

  HighsStatus deleteColsGeneral(bool interval, int from_col, int to_col,
                                bool set, int num_set_entries,
                                const int* col_set, bool mask, int* col_mask);

  HighsStatus getCols(const int from_col, const int to_col, int& num_col,
                      double* col_cost, double* col_lower, double* col_upper,
                      int& num_nz, int* col_matrix_start, int* col_matrix_index,
                      double* col_matrix_value);

  HighsStatus getCols(const int num_set_entries, const int* col_set,
                      int& num_col, double* col_cost, double* col_lower,
                      double* col_upper, int& num_nz, int* col_matrix_start,
                      int* col_matrix_index, double* col_matrix_value);

  HighsStatus getCols(const int* col_mask, int& num_col, double* col_cost,
                      double* col_lower, double* col_upper, int& num_nz,
                      int* col_matrix_start, int* col_matrix_index,
                      double* col_matrix_value);

  HighsStatus getColsGeneral(const bool interval, const int from_col,
                             const int to_col, const bool set,
                             const int num_set_entries, const int* col_set,
                             const bool mask, const int* col_mask, int& num_col,
                             double* col_cost, double* col_lower,
                             double* col_upper, int& num_nz,
                             int* col_matrix_start, int* col_matrix_index,
                             double* col_matrix_value);

  HighsStatus getRows(const int from_row, const int to_row, int& num_row,
                      double* row_lower, double* row_upper, int& num_nz,
                      int* row_matrix_start, int* row_matrix_index,
                      double* row_matrix_value);

  HighsStatus getRows(const int num_set_entries, const int* row_set,
                      int& num_row, double* row_lower, double* row_upper,
                      int& num_nz, int* row_matrix_start, int* row_matrix_index,
                      double* row_matrix_value);

  HighsStatus getRows(const int* row_mask, int& num_row, double* row_lower,
                      double* row_upper, int& num_nz, int* row_matrix_start,
                      int* row_matrix_index, double* row_matrix_value);

  HighsStatus getRowsGeneral(const bool interval, const int from_row,
                             const int to_row, const bool set,
                             const int num_set_entries, const int* row_set,
                             const bool mask, const int* row_mask, int& num_row,
                             double* row_lower, double* row_upper, int& num_nz,
                             int* row_matrix_start, int* row_matrix_index,
                             double* row_matrix_value);

  HighsStatus getCoefficient(const int Xrow, const int Xcol, double& value);

  /**
   * @brief Add a contiguous set of rows to the model data---making them basic
   */
  HighsStatus addRows(int XnumNewRow, const double* XrowLower,
                      const double* XrowUpper, int XnumNewNZ,
                      const int* XARstart, const int* XARindex,
                      const double* XARvalue);

  HighsStatus deleteRows(int from_row, int to_row);

  HighsStatus deleteRows(int num_set_entries, const int* row_set);

  HighsStatus deleteRows(int* row_mask);

  HighsStatus deleteRowsGeneral(bool interval, int from_row, int to_row,
                                bool set, int num_set_entries,
                                const int* row_set, bool mask, int* row_mask);

  HighsStatus changeCoefficient(const int Xrow, const int Xcol,
                                const double XnewValue);

  // Shift the objective
  void shiftObjectiveValue(const double Xshift);

  // Utilities to get/change costs and bounds
  // Change the objective sense
  HighsStatus changeObjectiveSense(const ObjSense Xsense);

  // Change the costs for an interval of columns
  HighsStatus changeCosts(int from_col, int to_col, const double* usr_col_cost);

  // Change the costs from an ordered set of indices
  HighsStatus changeCosts(int num_set_entries, const int* col_set,
                          const double* usr_col_cost);

  // Change the costs with a mask
  HighsStatus changeCosts(const int* col_mask, const double* usr_col_cost);

  HighsStatus changeCostsGeneral(bool interval, int from_col, int to_col,
                                 bool set, int num_set_entries,
                                 const int* col_set, bool mask,
                                 const int* col_mask,
                                 const double* usr_col_cost);

  // Change the bounds for an interval of columns
  HighsStatus changeColBounds(int from_col, int to_col,
                              const double* usr_col_lower,
                              const double* usr_col_upper);

  // Change the bounds from an ordered set of indices
  HighsStatus changeColBounds(int num_set_entries, const int* col_set,
                              const double* usr_col_lower,
                              const double* usr_col_upper);

  // Change the bounds with a mask
  HighsStatus changeColBounds(const int* col_mask, const double* usr_col_lower,
                              const double* usr_col_upper);

  HighsStatus changeColBoundsGeneral(bool interval, int from_col, int to_col,
                                     bool set, int num_set_entries,
                                     const int* col_set, bool mask,
                                     const int* col_mask,
                                     const double* usr_col_lower,
                                     const double* usr_col_upper);

  // Change the bounds for an interval of rows
  HighsStatus changeRowBounds(int from_row, int to_row,
                              const double* usr_row_lower,
                              const double* usr_row_upper);

  // Change the bounds from an ordered set of indices
  HighsStatus changeRowBounds(int num_set_entries, const int* row_set,
                              const double* usr_row_lower,
                              const double* usr_row_upper);

  // Change the bounds with a mask
  HighsStatus changeRowBounds(const int* row_mask, const double* usr_row_lower,
                              const double* usr_row_upper);

  HighsStatus changeRowBoundsGeneral(bool interval, int from_row, int to_row,
                                     bool set, int num_set_entries,
                                     const int* row_set, bool mask,
                                     const int* row_mask,
                                     const double* usr_row_lower,
                                     const double* usr_row_upper);

  HighsStatus basisSolve(const vector<double>& rhs, double* solution,
                         int* solution_num_nz, int* solution_nz_indices,
                         bool transpose = false);

#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void change_update_method(int updateMethod);
#endif

  /**
   * @brief Convert a SCIP baseStat for columns and rows to HiGHS basis
   * Postive  return value k implies invalid basis status for column k-1
   * Negative return value k implies invalid basis status for row   -k-1
   */
  int convertBaseStatToHighsBasis(const int* cstat,  //!> Column baseStat
                                  const int* rstat   //!> Row baseStat
  );

  /**
   * @brief Convert a HiGHS basis to SCIP baseStat for columns and rows
   * Postive  return value k implies invalid basis status for column k-1
   * Negative return value k implies invalid basis status for row   -k-1
   */
  int convertHighsBasisToBaseStat(int* cstat,  //!> Column baseStat
                                  int* rstat   //!> Row baseStat
  );

  /**
   * @brief Convert a simplex basis to a HiGHS basis
   */
  void convertSimplexToHighsBasis();

  /**
   * @brief Convert a HiGHS basis to a simplex basis
   */
  void convertHighsToSimplexBasis();
  /**
   * @brief Convert a simplex solution to a HiGHS solution
   */
  void convertSimplexToHighsSolution();

  /**
   * @brief Get the indices of the basic variables for SCIP
   */
  int get_basic_indices(int* bind  //!> Indices of basic variables
  );
};

#endif /* SIMPLEX_HIGHSSIMPLEXINTERFACE_H_ */
