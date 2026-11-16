/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsHessian.cpp
 * @brief
 */
#include "model/HighsHessian.h"

#include <cassert>
#include <cstdio>

#include "model/HighsHessianUtils.h"
#include "util/HighsMatrixUtils.h"
#include "util/HighsRandom.h"

void HighsHessian::clear() {
  this->name_ = "Hessian";
  this->dim_ = 0;
  this->start_.clear();
  this->index_.clear();
  this->value_.clear();
  this->format_ = HessianFormat::kTriangular;
  this->start_.assign(1, 0);
  this->oracle_.clear();
}

void HighsHessian::exactResize() {
  if (this->dim_) {
    this->start_.resize(this->dim_ + 1);
    HighsInt num_nz = this->start_[this->dim_];
    this->index_.resize(num_nz);
    this->value_.resize(num_nz);
  } else if (!this->isOracle()) {
    this->clear();
  }
}

void HighsHessian::deleteCols(const HighsIndexCollection& index_collection) {
  if (this->dim_ == 0 || this->isOracle()) return;
  // Can't handle non-triangular matrices yet
  assert(this->format_ == HessianFormat::kTriangular);
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;
  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  // Initial pass creates a look-up to for the new index of columns
  // being retained, and -1 for columns being deleted
  std::vector<HighsInt> new_index;
  new_index.assign(this->dim_, -1);
  HighsInt new_dim = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      for (HighsInt iCol = 0; iCol < delete_from_col; iCol++)
        new_index[iCol] = new_dim++;
    }
    for (HighsInt iCol = keep_from_col; iCol <= keep_to_col; iCol++)
      new_index[iCol] = new_dim++;
    // When using a mask, to_k = this->dim_, but consecutive
    // keep/delete entries are accumulated, so may not need all passes
    if (keep_to_col >= this->dim_ - 1) break;
  }
  assert(new_dim < this->dim_);
  // Now perform the pass that deletes rows/columns
  keep_to_col = -1;
  current_set_entry = 0;
  // Have to accumulate new number of entries from the outset, as
  // entries may be lost from columns being kept before any are
  // deleted. Also keep a count of the number of nonzeros, in case the new
  HighsInt check_new_dim = new_dim;
  new_dim = 0;
  HighsInt new_num_nz = 0;
  HighsInt new_num_entries = 0;
  std::vector<HighsInt> save_start = this->start_;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      for (HighsInt iCol = 0; iCol < delete_from_col; iCol++) {
        assert(new_index[iCol] >= 0);
        for (HighsInt iEl = save_start[iCol]; iEl < save_start[iCol + 1];
             iEl++) {
          HighsInt iRow = new_index[this->index_[iEl]];
          if (iRow < 0) continue;
          this->index_[new_num_entries] = iRow;
          this->value_[new_num_entries] = this->value_[iEl];
          if (this->value_[new_num_entries]) new_num_nz++;
          new_num_entries++;
        }
        new_dim++;
        this->start_[new_dim] = new_num_entries;
      }
      assert(new_dim == delete_from_col);
    }
    for (HighsInt iCol = keep_from_col; iCol <= keep_to_col; iCol++) {
      assert(new_index[iCol] >= 0);
      for (HighsInt iEl = save_start[iCol]; iEl < save_start[iCol + 1]; iEl++) {
        HighsInt iRow = new_index[this->index_[iEl]];
        if (iRow < 0) continue;
        this->index_[new_num_entries] = iRow;
        this->value_[new_num_entries] = this->value_[iEl];
        if (this->value_[new_num_entries]) new_num_nz++;
        new_num_entries++;
      }
      new_dim++;
      this->start_[new_dim] = new_num_entries;
    }
    if (keep_to_col >= this->dim_ - 1) break;
  }
  assert(new_dim == check_new_dim);
  this->dim_ = new_dim;
  if (!new_num_nz) {
    this->clear();
  } else {
    this->exactResize();
  }
}

bool HighsHessian::scaleOk(const HighsInt hessian_scale,
                           const double small_matrix_value,
                           const double large_matrix_value) const {
  if (!this->dim_) return true;
  assert(!this->isOracle());
  if (this->isOracle()) return false;
  double hessian_scale_value = std::pow(2, hessian_scale);
  for (HighsInt iEl = 0; iEl < this->start_[this->dim_]; iEl++) {
    double abs_new_value = std::abs(this->value_[iEl] * hessian_scale_value);
    if (abs_new_value >= large_matrix_value) return false;
    if (abs_new_value <= small_matrix_value) return false;
  }
  return true;
}

HighsInt HighsHessian::dim() const {
  return this->isOracle() ? this->oracle_.dim_ : this->dim_;
}

HighsInt HighsHessian::numNz() const {
  if (this->isOracle()) return 0;
  assert(this->formatOk());
  assert((HighsInt)this->start_.size() >= this->dim_ + 1);
  return this->start_[this->dim_];
}

void HighsHessian::print(const std::string& message) const {
  if (this->isOracle()) {
    printf("Hessian is represented via an oracle\n");
    return;
  }
  HighsInt num_nz = this->numNz();
  printf("%s Hessian of dimension %" HIGHSINT_FORMAT " and %" HIGHSINT_FORMAT
         " entries: %s\n",
         this->format_ == HessianFormat::kTriangular ? "Triangular" : "Square",
         dim_, num_nz, message.c_str());
  printf("Start; Index; Value of sizes %d; %d; %d\n", (int)this->start_.size(),
         (int)this->index_.size(), (int)this->value_.size());
  if (dim_ <= 0) return;
  printf(" Row|");
  for (int iCol = 0; iCol < dim_; iCol++) printf(" %4d", iCol);
  printf("\n");
  printf("-----");
  for (int iCol = 0; iCol < dim_; iCol++) printf("-----");
  printf("\n");
  std::vector<std::string> col;
  col.assign(dim_, "");
  for (HighsInt iCol = 0; iCol < dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1]; iEl++)
      col[this->index_[iEl]] = highsFormatToString("%4g", this->value_[iEl]);
    printf("%4d|", (int)iCol);
    for (int iRow = 0; iRow < dim_; iRow++) printf(" %4s", col[iRow].c_str());
    printf("\n");
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1]; iEl++)
      col[this->index_[iEl]] = "";
  }
}
bool HighsHessian::operator==(const HighsHessian& hessian) const {
  bool equal = true;
  equal = this->dim_ == hessian.dim_ && equal;
  equal = this->start_ == hessian.start_ && equal;
  equal = this->index_ == hessian.index_ && equal;
  equal = this->value_ == hessian.value_ && equal;
  return equal;
}

void HighsHessian::formFromOracle() {
  assert(this->isOracle());
  assert(this->oracle_.isValid());
  this->start_.clear();
  this->index_.clear();
  this->value_.clear();
  HighsInt col_nnz;
  HighsInt dim = this->oracle_.dim_;
  std::vector<HighsInt> index(dim, 0);
  std::vector<double> value(dim, 0);
  HighsInt hessian_nnz = 0;
  this->dim_ = dim;
  this->format_ = HessianFormat::kTriangular;
  this->start_.push_back(hessian_nnz);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    this->oracle_.getPackedColumn(iCol, col_nnz, index.data(), value.data());
    // Set the iagonal value to zero in case it's not set to a nonzero
    // value by passing through the packed values
    HighsInt diag_el = hessian_nnz++;
    this->index_.push_back(iCol);
    this->value_.push_back(0.0);
    for (HighsInt iX = 0; iX < col_nnz; iX++) {
      HighsInt iRow = index[iX];
      double entry = value[iX];
      value[iX] = 0;
      if (iRow == iCol) {
        this->value_[diag_el] = entry;
      } else if (iRow > iCol) {
        this->index_.push_back(iRow);
        this->value_.push_back(entry);
        hessian_nnz++;
      }
    }
    this->start_.push_back(hessian_nnz);
  }
}

void HighsHessian::product(const std::vector<double>& solution,
                           std::vector<double>& product) const {
  HighsInt dim = this->dim();
  if (dim <= 0) return;
  product.assign(dim, 0);
  if (this->isOracle()) {
    this->oracle_.product(solution, product);
    return;
  }
  const bool triangular = this->format_ == HessianFormat::kTriangular;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      const HighsInt iRow = this->index_[iEl];
      product[iRow] += this->value_[iEl] * solution[iCol];
      if (triangular && iRow != iCol)
        product[iCol] += this->value_[iEl] * solution[iRow];
    }
  }
}

void HighsHessian::alphaProductPlusY(const double alpha,
                                     const std::vector<double>& x,
                                     std::vector<double>& y) const {
  assert(!this->isOracle());
  /*
    HighsInt dim = this->oracle_.dim_;
    assert(static_cast<size_t>(dim) == y.size());
    assert(static_cast<size_t>(dim) == x.size());
    std::vector<double> q_x(dim, 0);
    this->oracle_.product(x, q_x);
    for (HighsInt iCol = 0; iCol < dim; iCol++) y[iCol] += alpha * q_x[iCol];
    assert(111 == 888);
    return;
  }
  */
  if (this->dim_ <= 0) return;

  const bool triangular = this->format_ == HessianFormat::kTriangular;
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      const HighsInt iRow = this->index_[iEl];
      y[iRow] += alpha * this->value_[iEl] * x[iCol];
      if (triangular && iRow != iCol)
        y[iCol] += alpha * this->value_[iEl] * x[iRow];
    }
  }
}

double HighsHessian::objectiveValue(const std::vector<double>& solution) const {
  HighsInt dim = this->dim();
  std::vector<double> q_solution(dim, 0);
  this->product(solution, q_solution);
  double objective_function_value = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    objective_function_value += solution[iCol] * q_solution[iCol];
  objective_function_value *= 0.5;
  return objective_function_value;
}

HighsCDouble HighsHessian::objectiveCDoubleValue(
    const std::vector<double>& solution) const {
  if (this->isOracle()) {
    assert(111 == 888);
    return static_cast<HighsCDouble>(this->objectiveValue(solution));
  }
  HighsCDouble objective_function_value = HighsCDouble(0);
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    // Assumes Hessian format is triangular
    assert(this->index_[iEl] == iCol);
    objective_function_value +=
        0.5 * solution[iCol] * this->value_[iEl] * solution[iCol];
    for (HighsInt iEl = this->start_[iCol] + 1; iEl < this->start_[iCol + 1];
         iEl++)
      objective_function_value +=
          solution[iCol] * this->value_[iEl] * solution[this->index_[iEl]];
  }
  return objective_function_value;
}

bool HighsHessian::isDiagonal() const {
  if (this->isOracle()) {
    assert(111 == 888);
    return false;
  }
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      const HighsInt iRow = this->index_[iEl];
      if (iRow != iCol) return false;
    }
  }
  return true;
}

double HighsHessian::diag(HighsInt i) const { return this->diagonal(i); }
double HighsHessian::diagonal(HighsInt i) const {
  if (this->isOracle()) return this->oracle_.diagonal(i);
  // Assumes that the diagonal entry is always first, possibly with explicit
  // zero value
  assert(i < dim_);
  assert(index_[start_[i]] == i);
  return value_[start_[i]];
}

HighsHessian HighsHessian::toSquare() const {
  if (this->isOracle()) {
    assert(111 == 888);
    return *this;
  }
  if (this->format_ == HessianFormat::kSquare) return *this;
  assert(this->format_ == HessianFormat::kTriangular);
  std::vector<HighsInt> iwork(this->dim_, 0);
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    HighsInt iRow = this->index_[iEl];
    assert(iRow == iCol);
    iEl++;
    for (; iEl < this->start_[iCol + 1]; iEl++) {
      HighsInt iRow = this->index_[iEl];
      assert(iRow > iCol);
      iwork[iRow]++;
    }
  }
  HighsInt triangular_hessian_off_diagonal = this->numNz() - this->dim_;
  HighsHessian square_hessian;
  square_hessian.format_ = HessianFormat::kSquare;
  square_hessian.dim_ = this->dim_;
  square_hessian.start_[0] = 0;
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt square_hessian_col_nnz =
        iwork[iCol] + this->start_[iCol + 1] - this->start_[iCol];
    square_hessian.start_.push_back(square_hessian.start_[iCol] +
                                    square_hessian_col_nnz);
    iwork[iCol] = square_hessian.start_[iCol] + 1;
  }
  HighsInt square_hessian_nnz = square_hessian.numNz();
  assert(square_hessian_nnz ==
         this->dim_ + 2 * triangular_hessian_off_diagonal);
  square_hessian.index_.resize(square_hessian_nnz);
  square_hessian.value_.resize(square_hessian_nnz);
  for (HighsInt iCol = 0; iCol < this->dim_; iCol++) {
    HighsInt iEl = this->start_[iCol];
    HighsInt iRow = this->index_[iEl];
    HighsInt square_iEl = square_hessian.start_[iCol];
    square_hessian.index_[square_iEl] = iCol;
    square_hessian.value_[square_iEl] = this->value_[iEl];
    iEl++;
    for (; iEl < this->start_[iCol + 1]; iEl++) {
      HighsInt iRow = this->index_[iEl];
      double value = this->value_[iEl];
      HighsInt square_iEl = iwork[iCol];
      square_hessian.index_[square_iEl] = iRow;
      square_hessian.value_[square_iEl] = value;
      iwork[iCol]++;
      square_iEl = iwork[iRow];
      square_hessian.index_[square_iEl] = iCol;
      square_hessian.value_[square_iEl] = value;
      iwork[iRow]++;
    }
  }
  return square_hessian;
}

bool HighsHessian::legalDimensions(const HighsLogOptions& log_options) const {
  assert(!this->isOracle());
  // Assess the Hessian dimensions and vector sizes
  std::vector<HighsInt> hessian_p_end;
  const bool partitioned = false;
  return legalMatrixDimensions(log_options, "HighsHessian", this->name_,
                               this->dim_, partitioned, this->start_,
                               hessian_p_end, this->index_, this->value_);
}

HighsStatus HighsHessian::trimDimensions(const HighsLogOptions& log_options) {
  assert(!this->isOracle());
  std::vector<HighsInt> hessian_p_end;
  const bool partitioned = false;
  return trimMatrixDimensions(log_options, "HighsHessian", this->name_,
                              this->dim_, partitioned, this->start_,
                              hessian_p_end, this->index_, this->value_);
}

HighsStatus HighsHessian::assess(const HighsOptions& options,
                                 const std::string* col_names) {
  return assessHessian(*this, options, col_names);
}

HighsStatus HighsHessian::checkOracle(const HighsLogOptions& log_options,
                                      const bool exit_on_first_error) const {
  if (!this->isOracle()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Hessian oracle is not defined\n");
    return HighsStatus::kError;
  }
  const HessianOracle& oracle = this->oracle_;
  HighsInt dim = oracle.dim_;

  HighsInt num_nz;
  std::vector<HighsInt> index(dim, 0);
  std::vector<double> value(dim, 0);
  std::vector<double> oracle_value(dim, 0);

  num_nz = 0;
  // A Hessian oracle product call is necessary
  if (!oracle.hasProductCall()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Hessian oracle has no product call\n");
    return HighsStatus::kError;
  }

  bool warning_found = false;
  bool error_found = false;
  // Set up a square Hessian corresponding to the oracle (naturally)
  // assuming that oracle.entry is implemented correctly
  HighsHessian hessian;
  hessian.dim_ = dim;
  hessian.format_ = HessianFormat::kSquare;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    for (HighsInt iRow = 0; iRow < dim; iRow++) {
      double value = oracle.entry(iCol, iRow);
      if (value) {
        hessian.index_.push_back(iRow);
        hessian.value_.push_back(value);
      }
    }
    hessian.start_.push_back(static_cast<HighsInt>(hessian.index_.size()));
  }

  // Check for symmetry
  double asymmetry;
  HighsInt num_illegal_asymmetry = 0;
  double min_illegal_asymmetry = kHighsInf;
  double max_illegal_asymmetry = 0;
  HighsInt num_ok_asymmetry = 0;
  double min_ok_asymmetry = kHighsInf;
  double max_ok_asymmetry = 0;

  auto asymmetryCheck = [&]() {
    if (asymmetry > kSquareHessianAsymmetryTolerance) {
      num_illegal_asymmetry++;
      min_illegal_asymmetry = std::min(asymmetry, min_illegal_asymmetry);
      max_illegal_asymmetry = std::max(asymmetry, max_illegal_asymmetry);
    } else if (asymmetry) {
      num_ok_asymmetry++;
      min_ok_asymmetry = std::min(asymmetry, min_ok_asymmetry);
      max_ok_asymmetry = std::max(asymmetry, max_ok_asymmetry);
    }
  };

  auto columnValuesEqual = [&](const HighsInt iCol, const HighsInt iRow,
                               const double oracle_v, const double true_v) {
    if (oracle_v == true_v) return true;
    highsLogUser(log_options, HighsLogType::kError,
                 "Hessian oracle yields column %d value %d of %g, not %g\n",
                 int(iCol), int(iRow), oracle_v, true_v);
    return false;
  };

  auto productValuesClose = [&](const HighsInt iCol, const HighsInt iRow,
                                const double oracle_v, const double true_v,
                                const double tolerance = 0.0) {
    double delta = std::fabs(oracle_v - true_v) / (1.0 + std::fabs(true_v));
    if (delta <= tolerance) return true;
    highsLogUser(log_options, HighsLogType::kError,
                 "Hessian oracle product %d yields value %d of %g, not %g, a "
                 "relative difference of %g\n",
                 int(iCol), int(iRow), oracle_v, true_v, delta);
    return false;
  };

  std::vector<double> column;
  column.assign(dim, 0);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    // Scatter the entries below the diagonal
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      if (iRow > iCol) column[iRow] = hessian.value_[iEl];
    }
    // Inspect the entries above the diagonal in the row corresponding
    // to this column
    for (HighsInt iX = iCol + 1; iX < dim; iX++) {
      for (HighsInt iEl = hessian.start_[iX]; iEl < hessian.start_[iX + 1];
           iEl++) {
        HighsInt iRow = hessian.index_[iEl];
        if (iRow == iCol) {
          // Found an entry above the diagonal in the row
          // corresponding to this column
          asymmetry = std::fabs(column[iX] - hessian.value_[iEl]);
          asymmetryCheck();
          column[iX] = 0;
          break;
        }
      }
    }
    // Check for missing values in the row corresponding to this column
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      if (iRow <= iCol) continue;
      asymmetry = std::fabs(column[iRow]);
      asymmetryCheck();
    }
  }
  if (num_ok_asymmetry) {
    highsLogUser(log_options, HighsLogType::kWarning,
                 "Hessian oracle contains %d non-symmetr%s in [%.2g, %.2g] "
                 "within tolerance of %.1g\n",
                 int(num_ok_asymmetry), num_ok_asymmetry == 1 ? "y" : "ies",
                 min_ok_asymmetry, max_ok_asymmetry,
                 kSquareHessianAsymmetryTolerance);
    warning_found = true;
  }
  if (num_illegal_asymmetry) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Hessian oracle contains %d non-symmetr%s in [%.2g, %.2g] "
                 "exceeding tolerance of %.1g\n",
                 int(num_illegal_asymmetry),
                 num_illegal_asymmetry == 1 ? "y" : "ies",
                 min_illegal_asymmetry, max_illegal_asymmetry,
                 kSquareHessianAsymmetryTolerance);
    if (exit_on_first_error) return HighsStatus::kError;
    error_found = true;
  }

  std::vector<double> check(dim, 0);
  HighsRandom random;
  // Test a full matrix-vector product
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    //    index[iCol] = iCol;
    value[iCol] = random.fraction();
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
         iEl++)
      check[hessian.index_[iEl]] += hessian.value_[iEl] * value[iCol];
  }
  // Firstly in one sensible call
  oracle.product(dim, nullptr, value.data(), oracle_value.data());
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    if (!productValuesClose(-1, iCol, oracle_value[iCol], check[iCol])) {
      if (exit_on_first_error) return HighsStatus::kError;
      error_found = true;
    }
    oracle_value[iCol] = 0;
  }
  // Secondly in another sensible call
  oracle.product(value.data(), oracle_value.data());
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    if (!productValuesClose(-1, iCol, oracle_value[iCol], check[iCol])) {
      if (exit_on_first_error) return HighsStatus::kError;
      error_found = true;
    }
    check[iCol] = 0;
    value[iCol] = 0;
    oracle_value[iCol] = 0;
  }
  // Test column extraction, then use the column to check product
  HighsInt oracle_num_el;
  std::vector<HighsInt> oracle_index(dim, 0);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    bool column_error_found = false;
    // Get the scattered column from the local Hessian and the oracle
    num_nz = 0;
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
         iEl++) {
      HighsInt iRow = hessian.index_[iEl];
      index[num_nz] = iRow;
      value[iRow] = hessian.value_[iEl];
      num_nz++;
    }
    oracle.getScatteredColumn(iCol, oracle_num_el, oracle_index.data(),
                              oracle_value.data());
    // NB Oracle may return explicit zeros, so cannot expect
    // oracle_num_el to equal num_nz
    //
    // Check oracle values are the same as column values (ie oracle
    // hasn't inserted any values)
    for (HighsInt iX = 0; iX < oracle_num_el; iX++) {
      HighsInt iRow = oracle_index[iX];
      if (!columnValuesEqual(iCol, iRow, oracle_value[iRow], value[iRow])) {
        if (exit_on_first_error) return HighsStatus::kError;
        error_found = true;
        column_error_found = true;
      }
    }
    // Check column values are the same as oracle values (ie oracle
    // hasn't missed any values), zeroing oracle_value
    for (HighsInt iX = 0; iX < num_nz; iX++) {
      HighsInt iRow = index[iX];
      if (!columnValuesEqual(iCol, iRow, oracle_value[iRow], value[iRow])) {
        if (exit_on_first_error) return HighsStatus::kError;
        error_found = true;
        column_error_found = true;
      }
      oracle_value[iRow] = 0;
    }
    // Check remaining oracle values are zero
    for (HighsInt iRow = 0; iRow < dim; iRow++) {
      if (!columnValuesEqual(iCol, iRow, oracle_value[iRow], 0.0)) {
        if (exit_on_first_error) return HighsStatus::kError;
        error_found = true;
        column_error_found = true;
      }
    }
    // Reinstate oracle nonzeros, scaling with a random number for
    // product experiement and zero value
    for (HighsInt iX = 0; iX < num_nz; iX++) {
      HighsInt iRow = index[iX];
      oracle_value[iRow] = value[iRow] * random.fraction();
      value[iRow] = 0;
    }
    if (!column_error_found) {
      // Compute the product in check directly,
      for (HighsInt iX = 0; iX < num_nz; iX++) {
        HighsInt iCol1 = index[iX];
        for (HighsInt iEl = hessian.start_[iCol1];
             iEl < hessian.start_[iCol1 + 1]; iEl++)
          check[hessian.index_[iEl]] +=
              hessian.value_[iEl] * oracle_value[iCol1];
      }
      // Compute the product in value using the oracle
      oracle.product(oracle_num_el, oracle_index.data(), oracle_value.data(),
                     value.data());
      // Check results are close
      for (HighsInt iRow = 0; iRow < dim; iRow++) {
        if (!productValuesClose(iCol, iRow, value[iRow], check[iRow], 1e-12)) {
          if (exit_on_first_error) return HighsStatus::kError;
          error_found = true;
        }
      }
      // Zero value and check
      for (HighsInt iRow = 0; iRow < dim; iRow++) {
        value[iRow] = 0;
        check[iRow] = 0;
      }
    }
    // Zero oracle entries
    for (HighsInt iX = 0; iX < num_nz; iX++) oracle_value[index[iX]] = 0;
  }

  if (error_found) return HighsStatus::kError;
  if (warning_found) return HighsStatus::kWarning;
  return HighsStatus::kOk;
}

void HessianOracle::clear() {
  this->dim_ = 0;
  this->multiplier_ = 1;
  this->shift_ = 0;
  this->call_ = nullptr;
  this->data_ = nullptr;
}

bool HessianOracle::hasProductCall() const {
  assert(this->call_);
  if (!this->call_) return false;
  // A Hessian oracle is valid if it has product call
  HighsInt num_nz = 0;
  HighsInt index = 0;
  double value = 0.0;
  return this->call_(kHessianOracleCallTypeProduct, &num_nz, &index, &value,
                     nullptr, nullptr, &value, this->data_) == 0;
}

double HessianOracle::diagonal(const HighsInt i) const { return entry(i, i); }

double HessianOracle::entry(const HighsInt i, const HighsInt j) const {
  assert(this->call_);
  double entry;
  HighsInt hessian_x_index = j;
  HighsInt status =
      this->call_(kHessianOracleCallTypeEntry, nullptr, &i, nullptr, nullptr,
                  &hessian_x_index, &entry, this->data_);
  if (status) {
    HighsInt num_entries;
    std::vector<HighsInt> index(this->dim_, 0);
    std::vector<double> value(this->dim_, 0);
    this->getScatteredColumn(i, num_entries, index.data(), value.data());
    entry = value[j];
  }
  entry *= this->multiplier_;
  if (i == j) entry += this->shift_;
  return entry;
}

void HessianOracle::getPackedColumn(const HighsInt col,
                                    HighsInt& col_num_entries,
                                    HighsInt* col_index,
                                    double* col_value) const {
  assert(col >= 0 && col < this->dim_);
  HighsInt status =
      this->call_(kHessianOracleCallTypeColumn, nullptr, &col, nullptr,
                  &col_num_entries, col_index, col_value, this->data_);
  if (status) {
    std::vector<double> x_value(this->dim_, 0);
    x_value[col] = 1.0;
    this->product(1, &col, x_value.data(), col_value);
    col_num_entries = 0;
    for (HighsInt iRow = 0; iRow < this->dim_; iRow++) {
      if (!col_value[iRow]) continue;
      col_value[col_num_entries] = col_value[iRow];
      col_index[col_num_entries] = iRow;
      col_num_entries++;
    }
  }
  if (this->multiplier_ != 1.0) {
    for (HighsInt iX = 0; iX < col_num_entries; iX++)
      col_value[iX] *= this->multiplier_;
  }
  if (this->shift_) {
    bool found = false;
    for (HighsInt iX = 0; iX < col_num_entries; iX++) {
      if (col_index[iX] == col) {
        col_value[iX] += this->shift_;  // Exploits value = 1
        found = true;
        break;
      }
    }
    if (!found) {
      col_index[col_num_entries] = col;
      col_value[col_num_entries] = this->shift_;
      col_num_entries++;
    }
  }
}

void HessianOracle::getScatteredColumn(const HighsInt col,
                                       HighsInt& col_num_entries,
                                       HighsInt* col_index,
                                       double* col_value) const {
  // Get the packed column values using col_value to avoid allocating
  // a full vector
  this->getPackedColumn(col, col_num_entries, col_index, col_value);
  // Copy col_value into col_packed and zero col_value
  std::vector<double> col_packed(col_num_entries);
  for (HighsInt iEl = 0; iEl < col_num_entries; iEl++) {
    col_packed[iEl] = col_value[iEl];
    col_value[iEl] = 0;
  }
  // Now scatter into col_value
  for (HighsInt iEl = 0; iEl < col_num_entries; iEl++) {
    HighsInt iRow = col_index[iEl];
    col_value[iRow] = col_packed[iEl];
  }
}

void HessianOracle::product(const std::vector<double>& x_value,
                            std::vector<double>& hessian_x_value) const {
  HighsInt dim = this->dim_;
  assert(static_cast<size_t>(dim) == x_value.size());
  assert(static_cast<size_t>(dim) == hessian_x_value.size());
  this->product(x_value.data(), hessian_x_value.data());
}

// For full x
void HessianOracle::product(const double* x_value,
                            double* hessian_x_value) const {
  assert(this->call_);
  this->call_(kHessianOracleCallTypeProduct, nullptr, nullptr, x_value, nullptr,
              nullptr, hessian_x_value, this->data_);
  this->scaleAndShift(nullptr, nullptr, x_value, hessian_x_value);
}

// For scattered, sparse, x. Also reasonable to call this with x_index
// as null pointer for full x
void HessianOracle::product(const HighsInt x_num_entries,
                            const HighsInt* x_index, const double* x_value,
                            double* hessian_x_value) const {
  assert(this->call_);
  if (x_num_entries == 0) return;
  if (x_index) {
    // Scattered, sparse, x
    assert(x_num_entries > 0);
    this->call_(kHessianOracleCallTypeProduct, &x_num_entries, x_index, x_value,
                nullptr, nullptr, hessian_x_value, this->data_);
  } else {
    // Passed x_index as null pointer, so assume full x
    this->product(x_value, hessian_x_value);
  }
  this->scaleAndShift(&x_num_entries, x_index, x_value, hessian_x_value);
}

void HessianOracle::scaleAndShift(const HighsInt* x_num_entries,
                                  const HighsInt* x_index,
                                  const double* x_value,
                                  double* hessian_x_value) const {
  // Compute multiplier_*Qx + scale_*x
  if (this->multiplier_ != 1.0) {
    for (HighsInt iRow = 0; iRow < this->dim_; iRow++)
      hessian_x_value[iRow] *= this->multiplier_;
  }
  if (this->shift_ != 0.0) {
    if (x_index != nullptr) {
      assert(x_num_entries != nullptr);
      assert(*x_num_entries > 0);
      for (HighsInt iX = 0; iX < *x_num_entries; iX++) {
        HighsInt iRow = x_index[iX];
        hessian_x_value[iRow] += this->shift_ * x_value[iRow];
      }
    } else {
      for (HighsInt iRow = 0; iRow < this->dim_; iRow++)
        hessian_x_value[iRow] += this->shift_ * x_value[iRow];
    }
  }
}
