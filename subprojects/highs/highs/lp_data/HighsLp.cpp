/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 */
#include "lp_data/HighsLp.h"

#include <cassert>

#include "lp_data/HighsLpUtils.h"
#include "util/HighsMatrixUtils.h"

bool HighsLp::isMip() const {
  size_t integrality_size = this->integrality_.size();
  if (integrality_size) {
    assert(static_cast<HighsInt>(integrality_size) == this->num_col_);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
      if (this->integrality_[iCol] != HighsVarType::kContinuous) return true;
  }
  return false;
}

bool HighsLp::hasInfiniteCost(const double infinite_cost) const {
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    if (this->col_cost_[iCol] >= infinite_cost) return true;
    if (this->col_cost_[iCol] <= -infinite_cost) return true;
  }
  return false;
}

bool HighsLp::hasSemiVariables() const {
  size_t integrality_size = this->integrality_.size();
  if (integrality_size) {
    assert(static_cast<HighsInt>(integrality_size) == this->num_col_);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
      if (this->integrality_[iCol] == HighsVarType::kSemiContinuous ||
          this->integrality_[iCol] == HighsVarType::kSemiInteger)
        return true;
  }
  return false;
}

bool HighsLp::hasMods() const {
  return this->mods_.save_non_semi_variable_index.size() > 0 ||
         this->mods_.save_inconsistent_semi_variable_index.size() > 0 ||
         this->mods_.save_relaxed_semi_variable_lower_bound_index.size() > 0 ||
         this->mods_.save_tightened_semi_variable_upper_bound_index.size() >
             0 ||
         this->mods_.save_inf_cost_variable_index.size() > 0;
}

bool HighsLp::needsMods(const double infinite_cost) const {
  assert(this->has_infinite_cost_ == this->hasInfiniteCost(infinite_cost));
  return this->has_infinite_cost_ || this->hasSemiVariables();
}

bool HighsLp::operator==(const HighsLp& lp) const {
  bool equal = equalButForNames(lp);
  equal = equalNames(lp) && equal;
  return equal;
}

bool HighsLp::equalButForNames(const HighsLp& lp) const {
  bool equal = equalButForScalingAndNames(lp);
  equal = equalScaling(lp) && equal;
  return equal;
}

bool HighsLp::equalButForScalingAndNames(const HighsLp& lp) const {
  const bool equal_vectors = equalVectors(lp);
  const bool equal_matrix = this->a_matrix_ == lp.a_matrix_;
#ifndef NDEBUG
  if (!equal_matrix) printf("HighsLp::equalButForNames: Unequal matrix\n");
#endif
  return equal_vectors && equal_matrix;
}

bool HighsLp::equalVectors(const HighsLp& lp) const {
  bool equal_vectors = true;
  equal_vectors = this->num_col_ == lp.num_col_ && equal_vectors;
  equal_vectors = this->num_row_ == lp.num_row_ && equal_vectors;
  equal_vectors = this->sense_ == lp.sense_ && equal_vectors;
  equal_vectors = this->offset_ == lp.offset_ && equal_vectors;
  equal_vectors = this->model_name_ == lp.model_name_ && equal_vectors;
  equal_vectors = this->col_cost_ == lp.col_cost_ && equal_vectors;
  equal_vectors = this->col_upper_ == lp.col_upper_ && equal_vectors;
  equal_vectors = this->col_lower_ == lp.col_lower_ && equal_vectors;
  equal_vectors = this->row_upper_ == lp.row_upper_ && equal_vectors;
  equal_vectors = this->row_lower_ == lp.row_lower_ && equal_vectors;
#ifndef NDEBUG
  if (!equal_vectors) printf("HighsLp::equalButForNames: Unequal vectors\n");
#endif
  return equal_vectors;
}

bool HighsLp::equalNames(const HighsLp& lp) const {
  bool equal = true;
  equal = this->objective_name_ == lp.objective_name_ && equal;
  equal = this->row_names_ == lp.row_names_ && equal;
  equal = this->col_names_ == lp.col_names_ && equal;
  return equal;
}

bool HighsLp::equalScaling(const HighsLp& lp) const {
  bool equal = true;
  equal = this->scale_.strategy == lp.scale_.strategy && equal;
  equal = this->scale_.has_scaling == lp.scale_.has_scaling && equal;
  equal = this->scale_.num_col == lp.scale_.num_col && equal;
  equal = this->scale_.num_row == lp.scale_.num_row && equal;
  equal = this->scale_.cost == lp.scale_.cost && equal;
  equal = this->scale_.col == lp.scale_.col && equal;
  equal = this->scale_.row == lp.scale_.row && equal;
#ifndef NDEBUG
  if (!equal) printf("HighsLp::equalScaling: Unequal scaling\n");
#endif
  return equal;
}

double HighsLp::objectiveValue(const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->num_col_);
  return double(objectiveCDoubleValue(solution));
}

HighsCDouble HighsLp::objectiveCDoubleValue(
    const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->num_col_);
  HighsCDouble objective_function_value = this->offset_;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    objective_function_value +=
        static_cast<HighsCDouble>(this->col_cost_[iCol]) * solution[iCol];
  return objective_function_value;
}

void HighsLp::setMatrixDimensions() {
  this->a_matrix_.num_col_ = this->num_col_;
  this->a_matrix_.num_row_ = this->num_row_;
}

void HighsLp::resetScale() {
  // Should allow user-supplied scale to be retained
  //  const bool dimensions_ok =
  //    this->scale_.num_col_ == this->num_col_ &&
  //    this->scale_.num_row_ == this->num_row_;
  this->clearScale();
}

void HighsLp::setFormat(const MatrixFormat format) {
  this->a_matrix_.setFormat(format);
}

void HighsLp::exactResize() {
  this->col_cost_.resize(this->num_col_);
  this->col_lower_.resize(this->num_col_);
  this->col_upper_.resize(this->num_col_);
  this->row_lower_.resize(this->num_row_);
  this->row_upper_.resize(this->num_row_);
  this->a_matrix_.exactResize();

  if ((int)this->col_names_.size()) this->col_names_.resize(this->num_col_);
  if ((int)this->row_names_.size()) this->row_names_.resize(this->num_row_);
  if ((int)this->integrality_.size()) this->integrality_.resize(this->num_col_);
}

bool HighsLp::okNames() const {
  if (this->col_names_.size() != static_cast<size_t>(this->num_col_))
    return false;
  if (this->row_names_.size() != static_cast<size_t>(this->num_row_))
    return false;
  // Look for names that are empty or have spaces
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    const std::string& name = this->col_names_[iCol];
    if (HighsInt(name.length()) == 0) return false;
    size_t space_pos = name.find(' ');
    if (space_pos != std::string::npos) return false;
  }
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    const std::string& name = this->row_names_[iRow];
    if (HighsInt(name.length()) == 0) return false;
    size_t space_pos = name.find(' ');
    if (space_pos != std::string::npos) return false;
  }
  return true;
}

void HighsLp::clear() {
  this->num_col_ = 0;
  this->num_row_ = 0;

  this->col_cost_.clear();
  this->col_lower_.clear();
  this->col_upper_.clear();
  this->row_lower_.clear();
  this->row_upper_.clear();

  this->a_matrix_.clear();

  this->sense_ = ObjSense::kMinimize;
  this->offset_ = 0;

  this->model_name_ = "";
  this->origin_name_ = "";
  this->objective_name_ = "";

  this->col_name_prefix_ = "";
  this->row_name_prefix_ = "";
  this->col_name_suffix_ = 0;
  this->row_name_suffix_ = 0;
  this->col_names_.clear();
  this->row_names_.clear();

  this->integrality_.clear();

  this->col_hash_.clear();
  this->row_hash_.clear();

  this->clearScale();
  this->is_scaled_ = false;
  this->is_moved_ = false;
  this->cost_row_location_ = -1;
  this->has_infinite_cost_ = false;
  this->fme_obj_col_ = -1;
  this->mods_.clear();
}

void HighsLp::clearScale() {
  this->scale_.strategy = kSimplexScaleStrategyOff;
  this->scale_.has_scaling = false;
  this->scale_.num_col = 0;
  this->scale_.num_row = 0;
  this->scale_.cost = 0;
  this->scale_.col.clear();
  this->scale_.row.clear();
}

void HighsLp::clearScaling() {
  this->unapplyScale();
  this->clearScale();
}

void HighsLp::applyScale() {
  // Ensure that any scaling is applied
  const HighsScale& scale = this->scale_;
  if (this->is_scaled_) {
    // Already scaled - so check that there is scaling and return
    assert(scale.has_scaling);
    return;
  }
  // No scaling currently applied
  this->is_scaled_ = false;
  if (scale.has_scaling) {
    // Apply the scaling to the bounds, costs and matrix, and record
    // that it has been applied
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      this->col_lower_[iCol] /= scale.col[iCol];
      this->col_upper_[iCol] /= scale.col[iCol];
      this->col_cost_[iCol] *= scale.col[iCol];
    }
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      this->row_lower_[iRow] *= scale.row[iRow];
      this->row_upper_[iRow] *= scale.row[iRow];
    }
    this->a_matrix_.applyScale(scale);
    this->is_scaled_ = true;
  }
}

void HighsLp::unapplyScale() {
  // Ensure that any scaling is not applied
  const HighsScale& scale = this->scale_;
  if (!this->is_scaled_) {
    // Not scaled so return
    return;
  }
  // Already scaled - so check that there is scaling
  assert(scale.has_scaling);
  // Unapply the scaling to the bounds, costs and matrix, and record
  // that it has been unapplied
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    this->col_lower_[iCol] *= scale.col[iCol];
    this->col_upper_[iCol] *= scale.col[iCol];
    this->col_cost_[iCol] /= scale.col[iCol];
  }
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    this->row_lower_[iRow] /= scale.row[iRow];
    this->row_upper_[iRow] /= scale.row[iRow];
  }
  this->a_matrix_.unapplyScale(scale);
  this->is_scaled_ = false;
}

void HighsLp::moveBackLpAndUnapplyScaling(HighsLp& lp) {
  assert(this->is_moved_ == true);
  *this = std::move(lp);
  this->unapplyScale();
  assert(this->is_moved_ == false);
}

void HighsLp::addColNames(const std::string& name, const HighsInt num_new_col) {
  // Don't add names if there are no columns being added
  if (this->num_col_ == 0) return;
  HighsInt col_names_size = this->col_names_.size();
  if (col_names_size <= 0) return;
  assert(col_names_size == this->num_col_ + num_new_col);
  // Handle the addition of user-defined names later
  assert(name == "");
  // Blank names for the new columns were added in
  // appendColsToLpVectors
}

void HighsLp::addRowNames(const std::string& name, const HighsInt num_new_row) {
  // Don't add names if there are no rows being added
  if (this->num_row_ == 0) return;
  HighsInt row_names_size = this->row_names_.size();
  if (row_names_size <= 0) return;
  assert(row_names_size == this->num_row_ + num_new_row);
  // Handle the addition of user-defined names later
  assert(name == "");
  // Blank names for the new rows were added in appendRowsToLpVectors
}

void HighsLp::deleteColsFromVectors(
    HighsInt& new_num_col, const HighsIndexCollection& index_collection) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  // Initialise new_num_col in case none is removed due to from_k > to_k
  new_num_col = this->num_col_;
  if (from_k > to_k) return;

  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  HighsInt col_dim = this->num_col_;
  new_num_col = 0;
  bool have_names = (this->col_names_.size() != 0);
  bool have_integrality = (this->integrality_.size() != 0);
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    // Account for the initial columns being kept
    if (k == from_k) new_num_col = delete_from_col;
    if (delete_to_col >= col_dim - 1) break;
    assert(delete_to_col < col_dim);
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      this->col_cost_[new_num_col] = this->col_cost_[col];
      this->col_lower_[new_num_col] = this->col_lower_[col];
      this->col_upper_[new_num_col] = this->col_upper_[col];
      if (have_names) this->col_names_[new_num_col] = this->col_names_[col];
      if (have_integrality)
        this->integrality_[new_num_col] = this->integrality_[col];
      new_num_col++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  this->col_cost_.resize(new_num_col);
  this->col_lower_.resize(new_num_col);
  this->col_upper_.resize(new_num_col);
  if (have_integrality) this->integrality_.resize(new_num_col);
  if (have_names) this->col_names_.resize(new_num_col);
}

void HighsLp::deleteRowsFromVectors(
    HighsInt& new_num_row, const HighsIndexCollection& index_collection) {
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  // Initialise new_num_row in case none is removed due to from_k > to_k
  new_num_row = this->num_row_;
  if (from_k > to_k) return;

  HighsInt delete_from_row;
  HighsInt delete_to_row;
  HighsInt keep_from_row;
  HighsInt keep_to_row = -1;
  HighsInt current_set_entry = 0;

  HighsInt row_dim = this->num_row_;
  new_num_row = 0;
  bool have_names = (HighsInt)this->row_names_.size() > 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_row, delete_to_row,
                     keep_from_row, keep_to_row, current_set_entry);
    if (k == from_k) {
      // Account for the initial rows being kept
      new_num_row = delete_from_row;
    }
    if (delete_to_row >= row_dim - 1) break;
    assert(delete_to_row < row_dim);
    for (HighsInt row = keep_from_row; row <= keep_to_row; row++) {
      this->row_lower_[new_num_row] = this->row_lower_[row];
      this->row_upper_[new_num_row] = this->row_upper_[row];
      if (have_names) this->row_names_[new_num_row] = this->row_names_[row];
      new_num_row++;
    }
    if (keep_to_row >= row_dim - 1) break;
  }
  this->row_lower_.resize(new_num_row);
  this->row_upper_.resize(new_num_row);
  if (have_names) this->row_names_.resize(new_num_row);
}

void HighsLp::deleteCols(const HighsIndexCollection& index_collection) {
  HighsInt new_num_col;
  this->deleteColsFromVectors(new_num_col, index_collection);
  this->a_matrix_.deleteCols(index_collection);
  this->num_col_ = new_num_col;
}

void HighsLp::deleteRows(const HighsIndexCollection& index_collection) {
  HighsInt new_num_row;
  this->deleteRowsFromVectors(new_num_row, index_collection);
  this->a_matrix_.deleteRows(index_collection);
  this->num_row_ = new_num_row;
}

bool HighsLp::legalDimensions(const std::string& message,
                              const HighsLogOptions& log_options) const {
  bool ok = true;
  const HighsInt num_col = this->num_col_;
  const HighsInt num_row = this->num_row_;
  if (!(num_col >= 0))
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on num_col = %d >= 0\n",
                 message.c_str(), (int)num_col);
  ok = num_col >= 0 && ok;
  if (!(num_row >= 0))
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on num_row = %d >= 0\n",
                 message.c_str(), (int)num_row);
  ok = num_row >= 0 && ok;
  if (!ok) return ok;

  HighsInt col_cost_size = this->col_cost_.size();
  HighsInt col_lower_size = this->col_lower_.size();
  HighsInt col_upper_size = this->col_upper_.size();
  HighsInt col_names_size = this->col_names_.size();
  HighsInt integrality_size = this->integrality_.size();
  bool legal_col_cost_size = col_cost_size >= num_col;
  bool legal_col_lower_size = col_lower_size >= num_col;
  bool legal_col_upper_size = col_upper_size >= num_col;
  bool legal_col_names_size = col_names_size == 0 || col_names_size >= num_col;
  bool legal_integrality_size =
      integrality_size == 0 || integrality_size >= num_col;
  if (!legal_col_cost_size)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on col_cost_.size() = %d < "
        "%d = num_col\n",
        message.c_str(), (int)col_cost_size, (int)num_col);
  ok = legal_col_cost_size && ok;

  if (!legal_col_lower_size)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on col_lower_.size() = %d "
                 "< %d = num_col\n",
                 message.c_str(), (int)col_lower_size, (int)num_col);
  ok = legal_col_lower_size && ok;

  if (!legal_col_upper_size)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on col_upper_.size() = %d "
                 "< %d = num_col\n",
                 message.c_str(), (int)col_upper_size, (int)num_col);
  ok = legal_col_upper_size && ok;

  if (!legal_integrality_size)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on integrality_.size() = %d "
        "< %d = num_col\n",
        message.c_str(), (int)integrality_size, (int)num_col);
  ok = legal_integrality_size && ok;

  if (!legal_col_names_size)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on col_names_.size() = %d "
                 "< %d = num_col\n",
                 message.c_str(), (int)col_names_size, (int)num_col);
  ok = legal_col_names_size && ok;

  bool legal_format = this->a_matrix_.format_ == MatrixFormat::kColwise ||
                      this->a_matrix_.format_ == MatrixFormat::kRowwise;
  if (!legal_format)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on a_matrix_.format\n",
                 message.c_str());
  ok = legal_format && ok;
  HighsInt num_vec;
  if (this->a_matrix_.isColwise()) {
    num_vec = num_col;
  } else {
    num_vec = num_row;
  }
  bool legal_matrix_dimensions = this->a_matrix_.legalDimensions(log_options);
  if (!legal_matrix_dimensions)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on a_matrix dimensions\n",
                 message.c_str());
  ok = legal_matrix_dimensions && ok;

  HighsInt row_lower_size = this->row_lower_.size();
  HighsInt row_upper_size = this->row_upper_.size();
  HighsInt row_names_size = this->row_names_.size();
  bool legal_row_lower_size = row_lower_size >= num_row;
  bool legal_row_upper_size = row_upper_size >= num_row;
  bool legal_row_names_size = row_names_size == 0 || row_names_size >= num_row;
  if (!legal_row_lower_size)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on row_lower_.size() = %d "
                 "< %d = num_row\n",
                 message.c_str(), (int)row_lower_size, (int)num_row);
  ok = legal_row_lower_size && ok;

  if (!legal_row_upper_size)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on row_upper_.size() = %d "
                 "< %d = num_row\n",
                 message.c_str(), (int)row_upper_size, (int)num_row);
  ok = legal_row_upper_size && ok;

  if (!legal_row_names_size)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on row_names_.size() = %d "
                 "< %d = num_row\n",
                 message.c_str(), (int)row_names_size, (int)num_row);
  ok = legal_row_names_size && ok;

  bool legal_a_matrix_num_col = this->a_matrix_.num_col_ == num_col;
  bool legal_a_matrix_num_row = this->a_matrix_.num_row_ == num_row;
  if (!legal_a_matrix_num_col)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on a_matrix.num_col_ = %d "
                 "!= %d = num_col\n",
                 message.c_str(), (int)this->a_matrix_.num_col_, (int)num_col);
  ok = legal_a_matrix_num_col && ok;

  if (!legal_a_matrix_num_row)
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails on a_matrix.num_row_ = %d "
                 "!= %d = num_row\n",
                 message.c_str(), (int)this->a_matrix_.num_row_, (int)num_row);
  ok = legal_a_matrix_num_row && ok;

  HighsInt scale_strategy = (HighsInt)this->scale_.strategy;
  bool legal_scale_strategy = scale_strategy >= 0;
  if (!legal_scale_strategy)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on scale_.scale_strategy\n",
        message.c_str());
  ok = legal_scale_strategy && ok;
  HighsInt scale_row_size = (HighsInt)this->scale_.row.size();
  HighsInt scale_col_size = (HighsInt)this->scale_.col.size();
  bool legal_scale_num_col = false;
  bool legal_scale_num_row = false;
  bool legal_scale_row_size = false;
  bool legal_scale_col_size = false;
  if (this->scale_.has_scaling) {
    legal_scale_num_col = this->scale_.num_col == num_col;
    legal_scale_num_row = this->scale_.num_row == num_row;
    legal_scale_row_size = scale_row_size >= num_row;
    legal_scale_col_size = scale_col_size >= num_col;
  } else {
    legal_scale_num_col = this->scale_.num_col == 0;
    legal_scale_num_row = this->scale_.num_row == 0;
    legal_scale_row_size = scale_row_size == 0;
    legal_scale_col_size = scale_col_size == 0;
  }
  if (!legal_scale_num_col)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on scale_.num_col = %d != %d\n",
        message.c_str(), (int)this->scale_.num_col,
        (int)(this->scale_.has_scaling ? num_col : 0));
  ok = legal_scale_num_col && ok;
  if (!legal_scale_num_row)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on scale_.num_row = %d != %d\n",
        message.c_str(), (int)this->scale_.num_row,
        (int)(this->scale_.has_scaling ? num_row : 0));
  ok = legal_scale_num_row && ok;
  if (!legal_scale_col_size)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on scale_.col.size() = %d %s %d\n",
        message.c_str(), (int)scale_col_size,
        this->scale_.has_scaling ? ">=" : "==",
        (int)(this->scale_.has_scaling ? num_col : 0));
  ok = legal_scale_col_size && ok;
  if (!legal_scale_row_size)
    highsLogUser(
        log_options, HighsLogType::kError,
        "LP dimension validation (%s) fails on scale_.row.size() = %d %s %d\n",
        message.c_str(), (int)scale_row_size,
        this->scale_.has_scaling ? ">=" : "==",
        (int)(this->scale_.has_scaling ? num_row : 0));
  ok = legal_scale_row_size && ok;
  if (!ok) {
    highsLogUser(log_options, HighsLogType::kError,
                 "LP dimension validation (%s) fails\n", message.c_str());
  }

  return ok;
}

HighsStatus HighsLp::trimDimensions(const HighsLogOptions& log_options) {
  const HighsInt num_col = this->num_col_;
  const HighsInt num_row = this->num_row_;

  HighsInt trimmed = 0;
  HighsStatus status = HighsStatus::kOk;

  std::string lp_name = "";
  if (this->model_name_ != "") {
    lp_name = this->model_name_;
    if (this->origin_name_ != "") {
      lp_name += "+" + this->origin_name_;
    }
    lp_name += " ";
  } else if (this->origin_name_ != "") {
    lp_name = this->origin_name_ + " ";
  }

  auto logTrimming = [&](const std::string& name) {
    if (trimmed == 0) return;
    assert(trimmed > 0);
    highsLogUser(log_options, HighsLogType::kWarning,
                 "Trimmed %d excess entr%s from %sHighsLp data member %s\n",
                 int(trimmed), highsIntToPlural(trimmed, true).c_str(),
                 lp_name.c_str(), name.c_str());
    status = HighsStatus::kWarning;
  };

  trimmed = this->col_cost_.size() - num_col;
  if (trimmed) {
    this->col_cost_.resize(num_col);
    logTrimming("col_cost_");
  }

  trimmed = this->col_lower_.size() - num_col;
  if (trimmed) {
    this->col_lower_.resize(num_col);
    logTrimming("col_lower_");
  }

  trimmed = this->col_upper_.size() - num_col;
  if (trimmed) {
    this->col_upper_.resize(num_col);
    logTrimming("col_upper_");
  }

  trimmed = this->row_lower_.size() - num_row;
  if (trimmed) {
    this->row_lower_.resize(num_row);
    logTrimming("row_lower_");
  }

  trimmed = this->row_upper_.size() - num_row;
  if (trimmed) {
    this->row_upper_.resize(num_row);
    logTrimming("row_upper_");
  }

  if (this->col_names_.size()) {
    trimmed = this->col_names_.size() - num_col;
    if (trimmed) {
      this->col_names_.resize(num_col);
      logTrimming("col_names_");
    }
  }

  if (this->row_names_.size()) {
    trimmed = this->row_names_.size() - num_row;
    if (trimmed) {
      this->row_names_.resize(num_row);
      logTrimming("row_names_");
    }
  }

  if (this->integrality_.size()) {
    trimmed = this->integrality_.size() - num_col;
    if (trimmed) {
      this->integrality_.resize(num_col);
      logTrimming("integrality_");
    }
  }

  HighsStatus matrix_trim_status = this->a_matrix_.trimDimensions(log_options);

  if (matrix_trim_status == HighsStatus::kWarning) return matrix_trim_status;
  return status;
}

HighsStatus HighsLp::assess(const HighsOptions& options) {
  // Just in case it's not yet been done
  this->setMatrixDimensions();
  return assessLp(*this, options);
}

void HighsLp::unapplyMods() {
  // Restore any non-semi types
  const HighsInt num_non_semi = this->mods_.save_non_semi_variable_index.size();
  for (HighsInt k = 0; k < num_non_semi; k++) {
    HighsInt iCol = this->mods_.save_non_semi_variable_index[k];
    assert(this->integrality_[iCol] == HighsVarType::kContinuous ||
           this->integrality_[iCol] == HighsVarType::kInteger);
    if (this->integrality_[iCol] == HighsVarType::kContinuous) {
      this->integrality_[iCol] = HighsVarType::kSemiContinuous;
    } else {
      this->integrality_[iCol] = HighsVarType::kSemiInteger;
    }
  }
  // Restore any inconsistent semi variables
  const HighsInt num_inconsistent_semi =
      this->mods_.save_inconsistent_semi_variable_index.size();
  if (!num_inconsistent_semi) {
    assert(
        !this->mods_.save_inconsistent_semi_variable_lower_bound_value.size());
    assert(
        !this->mods_.save_inconsistent_semi_variable_upper_bound_value.size());
    assert(!this->mods_.save_inconsistent_semi_variable_type.size());
  }
  for (HighsInt k = 0; k < num_inconsistent_semi; k++) {
    HighsInt iCol = this->mods_.save_inconsistent_semi_variable_index[k];
    this->col_lower_[iCol] =
        this->mods_.save_inconsistent_semi_variable_lower_bound_value[k];
    this->col_upper_[iCol] =
        this->mods_.save_inconsistent_semi_variable_upper_bound_value[k];
    this->integrality_[iCol] =
        this->mods_.save_inconsistent_semi_variable_type[k];
  }
  // Restore any relaxed lower bounds
  std::vector<HighsInt>& relaxed_semi_variable_lower_index =
      this->mods_.save_relaxed_semi_variable_lower_bound_index;
  std::vector<double>& relaxed_semi_variable_lower_value =
      this->mods_.save_relaxed_semi_variable_lower_bound_value;
  const HighsInt num_lower_bound = relaxed_semi_variable_lower_index.size();
  if (!num_lower_bound) {
    assert(!relaxed_semi_variable_lower_value.size());
  }
  for (HighsInt k = 0; k < num_lower_bound; k++) {
    HighsInt iCol = relaxed_semi_variable_lower_index[k];
    assert(this->integrality_[iCol] == HighsVarType::kSemiContinuous ||
           this->integrality_[iCol] == HighsVarType::kSemiInteger);
    this->col_lower_[iCol] = relaxed_semi_variable_lower_value[k];
  }
  // Restore any tightened upper bounds
  std::vector<HighsInt>& tightened_semi_variable_upper_bound_index =
      this->mods_.save_tightened_semi_variable_upper_bound_index;
  std::vector<double>& tightened_semi_variable_upper_bound_value =
      this->mods_.save_tightened_semi_variable_upper_bound_value;
  const HighsInt num_upper_bound =
      tightened_semi_variable_upper_bound_index.size();
  if (!num_upper_bound) {
    assert(!tightened_semi_variable_upper_bound_value.size());
  }
  for (HighsInt k = 0; k < num_upper_bound; k++) {
    HighsInt iCol = tightened_semi_variable_upper_bound_index[k];
    assert(this->integrality_[iCol] == HighsVarType::kSemiContinuous ||
           this->integrality_[iCol] == HighsVarType::kSemiInteger);
    this->col_upper_[iCol] = tightened_semi_variable_upper_bound_value[k];
  }

  this->mods_.clear();
}

void HighsLpMods::clear() {
  this->save_non_semi_variable_index.clear();
  this->save_inconsistent_semi_variable_index.clear();
  this->save_inconsistent_semi_variable_lower_bound_value.clear();
  this->save_inconsistent_semi_variable_upper_bound_value.clear();
  this->save_inconsistent_semi_variable_type.clear();
  this->save_relaxed_semi_variable_lower_bound_index.clear();
  this->save_relaxed_semi_variable_lower_bound_value.clear();
  this->save_tightened_semi_variable_upper_bound_index.clear();
  this->save_tightened_semi_variable_upper_bound_value.clear();
  this->save_inf_cost_variable_index.clear();
  this->save_inf_cost_variable_cost.clear();
  this->save_inf_cost_variable_lower.clear();
  this->save_inf_cost_variable_upper.clear();
}

bool HighsLpMods::isClear() {
  if (this->save_non_semi_variable_index.size()) return false;
  if (this->save_inconsistent_semi_variable_index.size()) return false;
  if (this->save_inconsistent_semi_variable_lower_bound_value.size())
    return false;
  if (this->save_inconsistent_semi_variable_upper_bound_value.size())
    return false;
  if (this->save_inconsistent_semi_variable_type.size()) return false;
  if (this->save_relaxed_semi_variable_lower_bound_index.size()) return false;
  if (this->save_relaxed_semi_variable_lower_bound_value.size()) return false;
  if (this->save_tightened_semi_variable_upper_bound_index.size()) return false;
  if (this->save_tightened_semi_variable_upper_bound_value.size()) return false;
  return true;
}

void HighsNameHash::form(const std::vector<std::string>& names) {
  this->clear();
  for (size_t index = 0; index < names.size(); index++) {
    auto emplace_result =
        this->name2index.emplace(names[index], static_cast<int>(index));
    const bool duplicate = !emplace_result.second;
    if (duplicate) {
      // Find the original
      auto& search = emplace_result.first;
      if (search->second >= 0) {
        // search->second is is the index associated with the first
        // occurrence of the name (not the position in name2index), so
        // initialise the duplicate count
        search->second = -1;
      } else {
        // search->second < 0, then -search->second is the duplicate
        // count, so increase it
        search->second--;
      }
    }
  }
  assert(this->ok(names));
}

bool HighsNameHash::hasDuplicate(const std::vector<std::string>& names) {
  this->clear();
  bool has_duplicate = false;
  for (size_t index = 0; index < names.size(); index++) {
    auto emplace_result =
        this->name2index.emplace(names[index], static_cast<int>(index));
    has_duplicate = !emplace_result.second;
    if (has_duplicate) break;
  }
  this->clear();
  return has_duplicate;
}

void HighsNameHash::addName(const HighsInt index, const std::string& name) {
  if (!this->size()) return;
  addNameFindsDuplicate(index, name);
}

bool HighsNameHash::addNameFindsDuplicate(const HighsInt index,
                                          const std::string& name) {
  auto emplace_result = this->name2index.emplace(name, static_cast<int>(index));
  const bool duplicate = !emplace_result.second;
  if (duplicate) {
    // Find the original and mark it as duplicate
    auto& search = emplace_result.first;
    if (search->second >= 0) {
      search->second = -1;
    } else {
      search->second--;
    }
    return true;
  }
  return false;
}

bool HighsNameHash::updateFindsDuplicate(const HighsInt index,
                                         const std::string& old_name,
                                         const std::string& new_name) {
  // Update duplication data for old_name, erasing it if it disappears
  auto find_result = this->name2index.find(old_name);
  if (find_result != this->name2index.end()) {
    // Old name is present, so find out if it's a duplicate
    if (find_result->second < 0) {
      // At least one duplicate of old_name
      find_result->second++;
    } else {
      // Old name was unique, so erase it
      this->name2index.erase(old_name);
    }
  }
  return addNameFindsDuplicate(index, new_name);
}

bool HighsNameHash::ok(const std::vector<std::string>& names) {
  HighsInt num_name = names.size();
  HighsNameHash helper;
  for (HighsInt index = 0; index < num_name; index++)
    helper.name2index.emplace(names[index], 0);
  for (auto& entry : helper.name2index) {
    std::string unique_name = entry.first;
    auto find_result = this->name2index.find(unique_name);
    if (find_result == this->name2index.end()) return false;
    HighsInt count = 0;
    for (HighsInt index = 0; index < num_name; index++)
      if (names[index] == unique_name) count++;
    if ((count == 1 && find_result->second < 0) ||
        (count > 1 && count != 1 - find_result->second))
      return false;
  }
  return true;
}
