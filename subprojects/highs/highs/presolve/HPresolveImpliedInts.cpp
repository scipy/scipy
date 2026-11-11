/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::StatusResult HPresolve::isImpliedIntegral(HighsInt col) {
  // check if the integer constraint on a variable is implied by the model
  assert(model->integrality_[col] == HighsVarType::kInteger);

  bool runDualDetection = true;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    // if not all other columns are integer, skip row and also do not try the
    // dual detection in the second loop as it must hold for all rows
    if (rowsize[nz.index()] < 2 ||
        rowsizeInteger[nz.index()] < rowsize[nz.index()]) {
      runDualDetection = false;
      continue;
    }

    double rowLower = isImpliedEquationAtUpper(nz.index())
                          ? model->row_upper_[nz.index()]
                          : model->row_lower_[nz.index()];

    double rowUpper = isImpliedEquationAtLower(nz.index())
                          ? model->row_lower_[nz.index()]
                          : model->row_upper_[nz.index()];

    if (rowUpper == rowLower) {
      // if there is an equation the dual detection does not need to be tried
      runDualDetection = false;
      double scale = 1.0 / nz.value();

      if (!rowCoefficientsIntegral(nz.index(), scale)) continue;

      if (fractionality(rowLower * scale) > primal_feastol)
        return StatusResult(Result::kPrimalInfeasible);

      return StatusResult(true);
    }
  }

  if (!runDualDetection) return StatusResult(false);

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    // if row coefficients are not integral, variable is not (implied) integral
    if (!rowCoefficientsIntegral(nz.index(), scale)) return StatusResult(false);
    if (model->row_upper_[nz.index()] != kHighsInf) {
      // right-hand side: scale, round down and unscale again
      double rUpper =
          std::abs(nz.value()) *
          std::floor(model->row_upper_[nz.index()] * std::abs(scale) +
                     primal_feastol);
      // check if modification is large enough
      if (std::abs(model->row_upper_[nz.index()] - rUpper) >
          options->small_matrix_value)
        changeRowUpper(nz.index(), rUpper);
    }
    if (model->row_lower_[nz.index()] != -kHighsInf) {
      // left-hand side: scale, round up and unscale again
      double rLower =
          std::abs(nz.value()) *
          std::ceil(model->row_lower_[nz.index()] * std::abs(scale) -
                    primal_feastol);
      // check if modification is large enough
      if (std::abs(model->row_lower_[nz.index()] - rLower) >
          options->small_matrix_value)
        changeRowLower(nz.index(), rLower);
    }
  }

  return StatusResult(true);
}

HPresolve::StatusResult HPresolve::isImpliedInteger(HighsInt col) const {
  // check if a continuous variable is implied integer
  assert(model->integrality_[col] == HighsVarType::kContinuous);

  bool runDualDetection = true;

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    // if not all other columns are integer, skip row and also do not try the
    // dual detection in the second loop as it must hold for all rows
    if (rowsize[nz.index()] < 2 ||
        rowsizeInteger[nz.index()] + rowsizeImplInt[nz.index()] <
            rowsize[nz.index()] - 1) {
      runDualDetection = false;
      continue;
    }

    double rowLower = isImpliedEquationAtUpper(nz.index())
                          ? model->row_upper_[nz.index()]
                          : model->row_lower_[nz.index()];

    double rowUpper = isImpliedEquationAtLower(nz.index())
                          ? model->row_lower_[nz.index()]
                          : model->row_upper_[nz.index()];

    if (rowUpper == rowLower) {
      // if there is an equation the dual detection does not need to be tried
      runDualDetection = false;
      double scale = 1.0 / nz.value();

      if (fractionality(rowLower * scale) > primal_feastol) continue;

      if (!rowCoefficientsIntegral(nz.index(), scale)) continue;

      return StatusResult(true);
    }
  }

  if (!runDualDetection) return StatusResult(false);

  if ((model->col_lower_[col] != -kHighsInf &&
       fractionality(model->col_lower_[col]) > options->small_matrix_value) ||
      (model->col_upper_[col] != kHighsInf &&
       fractionality(model->col_upper_[col]) > options->small_matrix_value))
    return StatusResult(false);

  for (const HighsSliceNonzero& nz : getColumnVector(col)) {
    double scale = 1.0 / nz.value();
    if (model->row_upper_[nz.index()] != kHighsInf &&
        fractionality(model->row_upper_[nz.index()] * scale) > primal_feastol)
      return StatusResult(false);

    if (model->row_lower_[nz.index()] != -kHighsInf &&
        fractionality(model->row_lower_[nz.index()] * scale) > primal_feastol)
      return StatusResult(false);

    if (!rowCoefficientsIntegral(nz.index(), scale)) return StatusResult(false);
  }

  return StatusResult(true);
}

HPresolve::StatusResult HPresolve::convertImpliedInteger(HighsInt col,
                                                         HighsInt row,
                                                         bool skipInputChecks) {
  // return if column was deleted
  if (colDeleted[col]) return StatusResult(false);

  // return if column is not continuous or cannot be converted to an implied
  // integer
  if (!skipInputChecks) {
    if (model->integrality_[col] != HighsVarType::kContinuous)
      return StatusResult(false);
    StatusResult impliedInteger = isImpliedInteger(col);
    if (!impliedInteger) return impliedInteger;
  }

  // convert to implied integer
  model->integrality_[col] = HighsVarType::kImplicitInteger;

  if (row != -1) {
    // use row index supplied by caller (e.g. singleton)
    ++rowsizeImplInt[row];
  } else {
    // iterate over rows
    for (const HighsSliceNonzero& nonzero : getColumnVector(col))
      ++rowsizeImplInt[nonzero.index()];
  }

  // Potentially strengthen existing bound using implied bounds
  // If not done then there may be a stronger fractional
  // implied bound for a non-continuous column, which causes
  // errors in rules downstream that assume integrality.
  // changeColBounds will perform rounding
  double newLower = model->col_lower_[col];
  double newUpper = model->col_upper_[col];
  if (implColLower[col] > newLower + primal_feastol)
    newLower = implColLower[col];
  if (implColUpper[col] < newUpper - primal_feastol)
    newUpper = implColUpper[col];
  // round and update bounds
  return StatusResult(changeColBounds(col, newLower, newUpper));
}

HPresolve::Result HPresolve::detectImpliedIntegers() {
  for (HighsInt col = 0; col != model->num_col_; ++col)
    HPRESOLVE_CHECKED_CALL(static_cast<Result>(convertImpliedInteger(col)));
  return Result::kOk;
}

}  // namespace presolve
