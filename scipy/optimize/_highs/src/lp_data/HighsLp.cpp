/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLp.h"

bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis) {
  bool consistent = true;
  consistent = isBasisRightSize(lp, basis) && consistent;
  int num_basic_variables = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    if (basis.col_status[iCol] == HighsBasisStatus::BASIC)
      num_basic_variables++;
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    if (basis.row_status[iRow] == HighsBasisStatus::BASIC)
      num_basic_variables++;
  }
  bool right_num_basic_variables = num_basic_variables == lp.numRow_;
  consistent = right_num_basic_variables && consistent;
  return consistent;
}

bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution) {
  bool right_size = true;
  right_size = (int)solution.col_value.size() == lp.numCol_ && right_size;
  right_size = (int)solution.col_dual.size() == lp.numCol_ && right_size;
  right_size = (int)solution.row_value.size() == lp.numRow_ && right_size;
  right_size = (int)solution.row_dual.size() == lp.numRow_ && right_size;
  return right_size;
}

bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis) {
  bool right_size = true;
  right_size = (int)basis.col_status.size() == lp.numCol_ && right_size;
  right_size = (int)basis.row_status.size() == lp.numRow_ && right_size;
  return right_size;
}

void clearSolutionUtil(HighsSolution& solution) {
  solution.col_dual.clear();
  solution.col_value.clear();
  solution.row_dual.clear();
  solution.row_value.clear();
}

void clearBasisUtil(HighsBasis& basis) {
  basis.row_status.clear();
  basis.col_status.clear();
  basis.valid_ = false;
}

void clearLp(HighsLp& lp) {
  lp.Astart_.clear();
  lp.Aindex_.clear();
  lp.Avalue_.clear();

  lp.col_names_.clear();
  lp.row_names_.clear();

  lp.sense_ = ObjSense::MINIMIZE;

  lp.colCost_.clear();
  lp.colLower_.clear();
  lp.colUpper_.clear();

  lp.integrality_.clear();
}

bool HighsLp::equalButForNames(const HighsLp& lp) {
  bool equal = true;
  equal = this->numCol_ == lp.numCol_ && equal;
  equal = this->numRow_ == lp.numRow_ && equal;
  equal = this->sense_ == lp.sense_ && equal;
  equal = this->offset_ == lp.offset_ && equal;
  equal = this->model_name_ == lp.model_name_ && equal;
  equal = this->colCost_ == lp.colCost_ && equal;
  equal = this->colUpper_ == lp.colUpper_ && equal;
  equal = this->colLower_ == lp.colLower_ && equal;
  equal = this->rowUpper_ == lp.rowUpper_ && equal;
  equal = this->rowLower_ == lp.rowLower_ && equal;
  equal = this->Astart_ == lp.Astart_ && equal;
  equal = this->Aindex_ == lp.Aindex_ && equal;
  equal = this->Avalue_ == lp.Avalue_ && equal;
  return equal;
}

bool HighsLp::operator==(const HighsLp& lp) {
  bool equal = equalButForNames(lp);
  equal = this->row_names_ == lp.row_names_ && equal;
  equal = this->col_names_ == lp.col_names_ && equal;
  return equal;
}
