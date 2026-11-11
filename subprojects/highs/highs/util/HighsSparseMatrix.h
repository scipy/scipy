/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSparseMatrix.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_SPARSE_MATRIX_H_
#define LP_DATA_HIGHS_SPARSE_MATRIX_H_

#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"  //For  HighsScale
#include "lp_data/HighsStatus.h"
#include "simplex/SimplexStruct.h"  //For SimplexScale until scaling is HighsScale
#include "util/HVector.h"
#include "util/HVectorBase.h"
#include "util/HighsSparseVectorSum.h"
#include "util/HighsUtils.h"

const double kHyperPriceDensity = 0.1;
const HighsInt kDebugReportOff = -2;
const HighsInt kDebugReportAll = -1;

class HighsSparseMatrix {
 public:
  HighsSparseMatrix() { clear(); }
  std::string name_;
  // Explicitly default the special member functions so that the move
  // operations (relied on by getRowwise, getColwise, ensureRowwise,
  // ensureColwise) can't silently vanish if a copy operation or the destructor
  // is ever customized in the future.
  HighsSparseMatrix(const HighsSparseMatrix&) = default;
  HighsSparseMatrix& operator=(const HighsSparseMatrix&) = default;
  HighsSparseMatrix(HighsSparseMatrix&&) = default;
  HighsSparseMatrix& operator=(HighsSparseMatrix&&) = default;
  ~HighsSparseMatrix() = default;

  MatrixFormat format_;
  HighsInt num_col_;
  HighsInt num_row_;
  std::vector<HighsInt> start_;
  std::vector<HighsInt> p_end_;
  std::vector<HighsInt> index_;
  std::vector<double> value_;

  bool operator==(const HighsSparseMatrix& matrix) const;
  bool equivalent(const HighsSparseMatrix& matrix) const;
  void clear();
  void exactResize();
  bool formatOk() const { return (this->isColwise() || this->isRowwise()); };
  bool isRowwise() const;
  bool isColwise() const;
  HighsInt numNz() const;
  void range(double& min_value, double& max_value) const;
  void setFormat(const MatrixFormat desired_format);
  void ensureColwise();
  void ensureRowwise();
  void setEmptyColwise(const HighsInt num_col);
  HighsSparseMatrix getRowwise() const;
  HighsSparseMatrix getColwise() const;
  HighsSparseMatrix oppositeFormat() const;
  void addVec(const HighsInt num_nz, const HighsInt* index, const double* value,
              const double multiple = 1);
  void addCols(const HighsSparseMatrix& new_cols,
               const int8_t* in_partition = NULL);
  void addRows(const HighsSparseMatrix& new_rows,
               const int8_t* in_partition = NULL);
  void getRow(const HighsInt iRow, HighsInt& num_nz, HighsInt* index,
              double* value) const;
  void getCol(const HighsInt iCol, HighsInt& num_nz, HighsInt* index,
              double* value) const;
  void deleteCols(const HighsIndexCollection& index_collection);
  void deleteRows(const HighsIndexCollection& index_collection);
  // Used only in Highs::writeLocalModel but (ToDo) why can't this
  // used assess?
  HighsStatus assessStart(const HighsLogOptions& log_options);
  HighsStatus assessIndexBounds(const HighsLogOptions& log_options);

  HighsInt numVec() const {
    return this->isColwise() ? this->num_col_ : this->num_row_;
  }
  HighsInt vecDim() const {
    return this->isColwise() ? this->num_row_ : this->num_col_;
  }
  HighsStatus assess(const HighsOptions& options,
                     const bool structure_only = false);
  HighsStatus assess(const HighsLogOptions& log_options,
                     const double small_matrix_value,
                     const double large_matrix_value, const bool sum_duplicates,
                     const std::string* col_names = nullptr,
                     const std::string* row_names = nullptr);
  HighsStatus assessStructure(const HighsOptions& options) {
    return this->assess(options, true);
  }
  HighsStatus assessExtremeValues(const std::string& name,
                                  const HighsOptions& options) const;
  bool hasExtremeValues(const double small_matrix_value,
                        const double large_matrix_value, double& min_abs_value,
                        double& max_abs_value) const;
  void considerColScaling(const HighsInt max_scale_factor_exponent,
                          double* col_scale);
  void considerRowScaling(const HighsInt max_scale_factor_exponent,
                          double* row_scale);
  void scaleCol(const HighsInt col, const double colScale);
  void scaleRow(const HighsInt row, const double rowScale);
  void applyScale(const HighsScale& scale);
  void applyRowScale(const HighsScale& scale);
  void applyColScale(const HighsScale& scale);
  void unapplyScale(const HighsScale& scale);
  void createSlice(const HighsSparseMatrix& matrix, const HighsInt from_col,
                   const HighsInt to_col);
  void createColwise(const HighsSparseMatrix& matrix);
  void createRowwise(const HighsSparseMatrix& matrix);
  void alphaProductPlusY(const double alpha, const std::vector<double>& x,
                         std::vector<double>& y,
                         const bool transpose = false) const;
  void product(vector<double>& result, const vector<double>& x) const;
  void productTranspose(vector<double>& result, const vector<double>& x) const;
  void productQuad(vector<double>& result, const vector<double>& x,
                   const HighsInt debug_report = kDebugReportOff) const;
  void productTransposeQuad(
      vector<double>& result_value, const vector<double>& x,
      const HighsInt debug_report = kDebugReportOff) const;
  void productTransposeQuad(
      vector<double>& result_value, vector<HighsInt>& result_index,
      const HVector& x, const HighsInt debug_report = kDebugReportOff) const;
  // Methods for PRICE, including the creation and updating of the
  // partitioned row-wise matrix
  void createRowwisePartitioned(const HighsSparseMatrix& matrix,
                                const int8_t* in_partition = NULL);
  bool debugPartitionOk(const int8_t* in_partition) const;
  void priceByColumn(const bool quad_precision, HVector& result,
                     const HVector& column,
                     const HighsInt debug_report = kDebugReportOff) const;
  void priceByRow(const bool quad_precision, HVector& result,
                  const HVector& column,
                  const HighsInt debug_report = kDebugReportOff) const;
  void priceByRowWithSwitch(
      const bool quad_precision, HVector& result, const HVector& column,
      const double expected_density, const HighsInt from_index,
      const double switch_density,
      const HighsInt debug_report = kDebugReportOff) const;
  void update(const HighsInt var_in, const HighsInt var_out,
              const HighsSparseMatrix& matrix);
  double computeDot(const HVector& column, const HighsInt use_col) const {
    return computeDot(column.array, use_col);
  }

  double computeDot(const std::vector<double>& array,
                    const HighsInt use_col) const;
  void collectAj(HVector& column, const HighsInt use_col,
                 const double multiplier) const;
  bool legalDimensions(const HighsLogOptions& log_options) const;
  HighsStatus trimDimensions(const HighsLogOptions& log_options);

 private:
  void priceByRowDenseResult(
      std::vector<double>& result, const HVector& column,
      const HighsInt from_index,
      const HighsInt debug_report = kDebugReportOff) const;
  void priceByRowDenseResult(
      std::vector<HighsCDouble>& result, const HVector& column,
      const HighsInt from_index,
      const HighsInt debug_report = kDebugReportOff) const;
  void debugReportRowPrice(const HighsInt iRow, const double multiplier,
                           const HighsInt to_iEl,
                           const vector<double>& result) const;
  void debugReportRowPrice(const HighsInt iRow, const double multiplier,
                           const HighsInt to_iEl,
                           const vector<HighsCDouble>& result) const;
  void debugReportRowPrice(const HighsInt iRow, const double multiplier,
                           const HighsInt to_iEl,
                           HighsSparseVectorSum& sum) const;
};

#endif
