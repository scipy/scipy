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

bool HighsLp::operator==(const HighsLp& lp) {
  bool equal = equalButForNames(lp);
  equal = this->row_names_ == lp.row_names_ && equal;
  equal = this->col_names_ == lp.col_names_ && equal;
  return equal;
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

void HighsLp::clear() {
  this->numCol_ = 0;
  this->numRow_ = 0;

  this->Astart_.clear();
  this->Aindex_.clear();
  this->Avalue_.clear();
  this->colCost_.clear();
  this->colLower_.clear();
  this->colUpper_.clear();
  this->rowLower_.clear();
  this->rowUpper_.clear();

  this->sense_ = ObjSense::MINIMIZE;
  this->offset_ = 0;

  this->model_name_ = "";
  this->lp_name_ = "";

  this->col_names_.clear();
  this->row_names_.clear();

  this->integrality_.clear();
}
