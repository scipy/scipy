/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsLpAggregator.h"

#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"

HighsLpAggregator::HighsLpAggregator(const HighsLpRelaxation& lprelaxation)
    : lprelaxation(lprelaxation) {
  vectorsum.setDimension(lprelaxation.getLp().num_row_ +
                         lprelaxation.getLp().num_col_);
}

void HighsLpAggregator::addRow(HighsInt row, double weight) {
  HighsInt len;
  const double* vals;
  const HighsInt* inds;
  lprelaxation.getRow(row, len, inds, vals);

  for (HighsInt i = 0; i != len; ++i) vectorsum.add(inds[i], weight * vals[i]);

  vectorsum.add(lprelaxation.getLp().num_col_ + row, -weight);
}

void HighsLpAggregator::getCurrentAggregation(std::vector<HighsInt>& inds,
                                              std::vector<double>& vals,
                                              bool negate) {
  const double droptol =
      lprelaxation.getMipSolver().options_mip_->small_matrix_value;
  const HighsDomain& domain = lprelaxation.getMipSolver().mipdata_->getDomain();
  const HighsInt numCol = lprelaxation.numCols();
  vectorsum.cleanup([&](const HighsInt col, const double val) {
    // only drop values for columns and not for slack variables of rows as the
    // former are the only ones which might be subject to numerical error. The
    // values for row slack variables are exact copies of the weights used
    // (assuming that separators only add each row once). Also all the
    // separators will only allow rows that do have meaningful coefficient
    // contributions with the used weight
    return col < numCol && termIsZero(val, domain.col_lower_[col],
                                      domain.col_upper_[col], droptol);
  });

  inds = vectorsum.getNonzeros();
  HighsInt len = inds.size();
  vals.resize(len);

  if (negate)
    for (HighsInt i = 0; i != len; ++i) vals[i] = -vectorsum.getValue(inds[i]);
  else
    for (HighsInt i = 0; i != len; ++i) vals[i] = vectorsum.getValue(inds[i]);
}

void HighsLpAggregator::clear() { vectorsum.clear(); }
