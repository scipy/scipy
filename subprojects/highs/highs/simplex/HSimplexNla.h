/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSimplexNla.h
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#ifndef HSIMPLEXNLA_H_
#define HSIMPLEXNLA_H_

#include "lp_data/HighsDebug.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/SimplexStruct.h"
#include "util/HFactor.h"

const HighsInt kReportItemLimit = 25;
const double kDensityForIndexing = 0.4;

struct ProductFormUpdate {
  bool valid_ = false;
  HighsInt num_row_;
  HighsInt update_count_;
  vector<HighsInt> pivot_index_;
  vector<double> pivot_value_;
  vector<HighsInt> start_;
  vector<HighsInt> index_;
  vector<double> value_;
  void clear();
  void setup(const HighsInt num_row, const double expected_density);
  HighsInt update(HVector* aq, HighsInt* iRow);
  void btran(HVector& rhs) const;
  void ftran(HVector& rhs) const;
};

struct SimplexIterate {
  bool valid_ = false;
  SimplexBasis basis_;
  InvertibleRepresentation invert_;
  std::vector<double> dual_edge_weight_;
  void clear();
};

class HSimplexNla {
 private:
  void setup(const HighsLp* lp, HighsInt* basic_index,
             const HighsOptions* options, HighsTimer* timer,
             HighsSimplexAnalysis* analysis,
             const HighsSparseMatrix* factor_a_matrix,
             const double factor_pivot_threshold);

  void setLpAndScalePointers(const HighsLp* lp);
  void setBasicIndexPointers(HighsInt* basic_index);
  void setPointers(const HighsLp* for_lp,
                   const HighsSparseMatrix* factor_a_matrix = NULL,
                   HighsInt* basic_index = NULL,
                   const HighsOptions* options = NULL, HighsTimer* timer = NULL,
                   HighsSimplexAnalysis* analysis = NULL);
  void clear();
  HighsInt invert();
  void btran(HVector& rhs, const double expected_density,
             HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void ftran(HVector& rhs, const double expected_density,
             HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranInScaledSpace(
      HVector& rhs, const double expected_density,
      HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void ftranInScaledSpace(
      HVector& rhs, const double expected_density,
      HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void update(HVector* aq, HVector* ep, HighsInt* iRow, HighsInt* hint);

  void putInvert();
  void getInvert();

  void transformForUpdate(HVector* column, HVector* row_ep,
                          const HighsInt variable_in, const HighsInt row_out);
  double variableScaleFactor(const HighsInt iVar) const;
  double basicColScaleFactor(const HighsInt iCol) const;
  double pivotInScaledSpace(const HVector* aq, const HighsInt variable_in,
                            const HighsInt row_out) const;
  void setPivotThreshold(const double new_pivot_threshold);

  void passLpPointer(const HighsLp* lp);
  void passScalePointer(const HighsScale* scale);
  void applyBasisMatrixColScale(HVector& rhs) const;
  void applyBasisMatrixRowScale(HVector& rhs) const;
  void unapplyBasisMatrixRowScale(HVector& rhs) const;
  double rowEp2NormInScaledSpace(const HighsInt iRow,
                                 const HVector& row_ep) const;
  void addCols(const HighsLp* updated_lp);
  void addRows(const HighsLp* updated_lp, HighsInt* basic_index,
               const HighsSparseMatrix* scaled_ar_matrix);
  bool sparseLoopStyle(const HighsInt count, const HighsInt dim,
                       HighsInt& to_entry) const;
  void reportVector(const std::string& message, const HighsInt num_index,
                    const vector<double>& vector_value,
                    const vector<HighsInt>& vector_index,
                    const bool force) const;
  void reportArray(const std::string& message, const HVector* vector,
                   const bool force = false) const;
  void reportArray(const std::string& message, const HighsInt offset,
                   const HVector* vector, const bool force = false) const;
  void reportArraySparse(const std::string& message, const HVector* vector,
                         const bool force = false) const;
  void reportArraySparse(const std::string& message, const HighsInt offset,
                         const HVector* vector, const bool force = false) const;
  void reportPackValue(const std::string& message, const HVector* vector,
                       const bool force = false) const;
  // Debug methods
  HighsDebugStatus debugCheckData(const std::string& message = "") const;
  HighsDebugStatus debugCheckInvert(const std::string& message,
                                    const HighsInt alt_debug_level = -1) const;
  double debugInvertResidualError(const bool transposed,
                                  const HVector& solution,
                                  HVector& residual) const;
  HighsDebugStatus debugReportInvertSolutionError(const bool transposed,
                                                  const HVector& true_solution,
                                                  const HVector& solution,
                                                  HVector& residual,
                                                  const bool force) const;
  HighsDebugStatus debugReportInvertSolutionError(
      const std::string& source, const bool transposed,
      const double solve_error_norm, const double residual_error_norm,
      const bool force) const;

  // References:
  //
  // Pointers:

  // Class data members
  const HighsLp* lp_;
  const HighsScale* scale_;
  HighsInt* basic_index_;
  const HighsOptions* options_;
  HighsTimer* timer_;
  HighsSimplexAnalysis* analysis_;

  HFactor factor_;

  bool report_;
  double build_synthetic_tick_;

  ProductFormUpdate update_;

  // Simplex iterate data
  SimplexIterate simplex_iterate_;

  friend class HEkk;
  friend class HEkkPrimal;
  friend class HEkkDual;
};

#endif /* HSIMPLEXNLA_H_ */
