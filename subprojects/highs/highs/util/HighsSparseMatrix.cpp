/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSparseMatrix.cpp
 * @brief
 */
#include "util/HighsSparseMatrix.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "util/HighsCDouble.h"
#include "util/HighsMatrixUtils.h"
#include "util/HighsSort.h"
#include "util/HighsSparseVectorSum.h"

using std::fabs;
using std::max;
using std::min;
using std::swap;
using std::vector;

bool HighsSparseMatrix::operator==(const HighsSparseMatrix& matrix) const {
  bool equal = true;
  equal = this->format_ == matrix.format_ && equal;
  equal = this->num_col_ == matrix.num_col_ && equal;
  equal = this->num_row_ == matrix.num_row_ && equal;
  equal = this->start_ == matrix.start_ && equal;
  equal = this->index_ == matrix.index_ && equal;
  equal = this->value_ == matrix.value_ && equal;
  return equal;
}

bool HighsSparseMatrix::equivalent(const HighsSparseMatrix& matrix) const {
  // Use this to test when the matrices may be in opposite format
  if (*this == matrix) return true;
  // Creat a local copy in of matrix in the same format as this and
  // test for it being equal to this
  HighsSparseMatrix lc_matrix;
  if (this->isColwise()) {
    lc_matrix = matrix.getColwise();
    assert(lc_matrix.isColwise());
    // lc_matrix.ensureColwise(); is now redundant
  } else {
    lc_matrix = matrix.getRowwise();
    assert(lc_matrix.isRowwise());
    // lc_matrix.ensureRowwise(); is now redundant
  }
  return *this == lc_matrix;
}

void HighsSparseMatrix::clear() {
  this->name_ = "";
  this->num_col_ = 0;
  this->num_row_ = 0;
  this->start_.clear();
  this->p_end_.clear();
  this->index_.clear();
  this->value_.clear();
  this->format_ = MatrixFormat::kColwise;
  this->start_.assign(1, 0);
}

void HighsSparseMatrix::exactResize() {
  if (this->isColwise()) {
    this->start_.resize(this->num_col_ + 1);
  } else {
    this->start_.resize(this->num_row_ + 1);
  }
  const HighsInt num_nz = this->isColwise() ? this->start_[this->num_col_]
                                            : this->start_[this->num_row_];
  if (this->format_ == MatrixFormat::kRowwisePartitioned) {
    this->p_end_.resize(this->num_row_);
  } else {
    assert((int)this->p_end_.size() == 0);
    this->p_end_.clear();
  }
  this->index_.resize(num_nz);
  this->value_.resize(num_nz);
}

bool HighsSparseMatrix::isRowwise() const {
  return this->format_ == MatrixFormat::kRowwise ||
         this->format_ == MatrixFormat::kRowwisePartitioned;
}

bool HighsSparseMatrix::isColwise() const {
  return this->format_ == MatrixFormat::kColwise;
}

HighsInt HighsSparseMatrix::numNz() const {
  assert(this->formatOk());
  if (this->isColwise()) {
    assert((HighsInt)this->start_.size() >= this->num_col_ + 1);
    return this->start_[this->num_col_];
  } else {
    assert((HighsInt)this->start_.size() >= this->num_row_ + 1);
    return this->start_[this->num_row_];
  }
}

void HighsSparseMatrix::range(double& min_value, double& max_value) const {
  assert(this->formatOk());
  for (HighsInt iEl = 0; iEl < this->numNz(); iEl++) {
    double value = fabs(this->value_[iEl]);
    min_value = min(min_value, value);
    max_value = max(max_value, value);
  }
}

void HighsSparseMatrix::setFormat(const MatrixFormat desired_format) {
  assert(this->formatOk());
  if (desired_format == MatrixFormat::kColwise) {
    this->ensureColwise();
  } else {
    this->ensureRowwise();
  }
  assert(this->format_ == desired_format);
}

void HighsSparseMatrix::ensureColwise() {
  // If a change of format is needed, the result of oppositeFormat() is built
  // in a temporary and move-assigned into *this, so no extra copy is made.
  assert(formatOk());

  if (isColwise()) return;

  HighsInt num_nz = numNz();
  (void)num_nz;  // avoid unused variable warning

  *this = oppositeFormat();

  assert(isColwise());
  assert(numNz() == num_nz);
  assert(static_cast<HighsInt>(start_.size()) == num_col_ + 1);
  assert(static_cast<HighsInt>(index_.size()) == numNz());
  assert(static_cast<HighsInt>(value_.size()) == numNz());
}
void HighsSparseMatrix::ensureRowwise() {
  // If a change of format is needed, the result of oppositeFormat() is built
  // in a temporary and move-assigned into *this, so no extra copy is made.
  assert(formatOk());

  if (isRowwise()) return;

  HighsInt num_nz = numNz();
  (void)num_nz;  // avoid unused variable warning

  *this = oppositeFormat();

  assert(isRowwise());
  assert(numNz() == num_nz);
  assert(static_cast<HighsInt>(start_.size()) == num_row_ + 1);
  assert(static_cast<HighsInt>(index_.size()) == numNz());
  assert(static_cast<HighsInt>(value_.size()) == numNz());
}

void HighsSparseMatrix::setEmptyColwise(const HighsInt num_col) {
  this->clear();
  assert(this->isColwise());
  this->start_.assign(num_col + 1, 0);
  this->num_col_ = num_col;
  this->num_row_ = 0;
}

void HighsSparseMatrix::addVec(const HighsInt num_nz, const HighsInt* index,
                               const double* value, const double multiple) {
  HighsInt num_vec = 0;
  if (this->isColwise()) {
    num_vec = this->num_col_;
  } else {
    num_vec = this->num_row_;
  }
  assert((int)this->start_.size() == num_vec + 1);
  assert((int)this->index_.size() == this->numNz());
  assert((int)this->value_.size() == this->numNz());
  for (HighsInt iEl = 0; iEl < num_nz; iEl++) {
    this->index_.push_back(index[iEl]);
    this->value_.push_back(multiple * value[iEl]);
  }
  this->start_.push_back(this->start_[num_vec] + num_nz);
  if (this->isColwise()) {
    this->num_col_++;
  } else {
    this->num_row_++;
  }
}

void HighsSparseMatrix::addCols(const HighsSparseMatrix& new_cols,
                                const int8_t* in_partition) {
  assert(new_cols.isColwise());
  const HighsInt num_new_col = new_cols.num_col_;
  const HighsInt num_new_nz = new_cols.numNz();
  const vector<HighsInt>& new_matrix_start = new_cols.start_;
  const vector<HighsInt>& new_matrix_index = new_cols.index_;
  const vector<double>& new_matrix_value = new_cols.value_;

  assert(this->formatOk());
  // Adding columns to a row-wise partitioned matrix needs the
  // partition information
  const bool partitioned = this->format_ == MatrixFormat::kRowwisePartitioned;
  // Cannot handle the row-wise partitioned case
  assert(!partitioned);
  if (partitioned) {
    //    if (in_partition == NULL) { printf("in_partition == NULL\n"); }
    assert(in_partition != NULL);
  }
  assert(num_new_col >= 0);
  assert(num_new_nz >= 0);
  if (num_new_col == 0) {
    // No columns are being added, so check that no nonzeros are being
    // added
    assert(num_new_nz == 0);
    return;
  }
  // Adding a positive number of columns to a matrix
  if (num_new_nz) {
    // Nonzeros are being added, so ensure that non-null data are
    // being passed
    assert(!new_matrix_start.empty());
    assert(!new_matrix_index.empty());
    assert(!new_matrix_value.empty());
  }
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  HighsInt num_nz = this->numNz();
  // Check that nonzeros aren't being appended to a matrix with no rows
  assert(num_new_nz <= 0 || num_row > 0);

  // If matrix is currently a standard row-wise matrix and there are
  // more new nonzeros than current nonzeros so flip column-wise
  if (this->format_ == MatrixFormat::kRowwise && num_new_nz > num_nz)
    this->ensureColwise();

  // Determine the new number of columns and nonzeros in the matrix
  HighsInt new_num_col = num_col + num_new_col;
  HighsInt new_num_nz = num_nz + num_new_nz;

  if (this->isColwise()) {
    // Matrix is column-wise
    this->start_.resize(new_num_col + 1);
    // Append the starts of the new columns
    if (num_new_nz) {
      // Nontrivial number of nonzeros being added, so use new_matrix_start
      for (HighsInt iNewCol = 0; iNewCol < num_new_col; iNewCol++)
        this->start_[num_col + iNewCol] = num_nz + new_matrix_start[iNewCol];
    } else {
      // No nonzeros being added, so new_matrix_start may be null, but entries
      // of zero are implied.
      for (HighsInt iNewCol = 0; iNewCol < num_new_col; iNewCol++)
        this->start_[num_col + iNewCol] = num_nz;
    }
    this->start_[num_col + num_new_col] = new_num_nz;
    // Update the number of columns
    this->num_col_ += num_new_col;
    // If no nonzeros are being added then there's nothing else to do
    if (num_new_nz <= 0) return;
    // Adding a non-trivial matrix: resize the column-wise matrix arrays
    // accordingly
    this->index_.resize(new_num_nz);
    this->value_.resize(new_num_nz);
    // Copy in the new indices and values
    for (HighsInt iEl = 0; iEl < num_new_nz; iEl++) {
      this->index_[num_nz + iEl] = new_matrix_index[iEl];
      this->value_[num_nz + iEl] = new_matrix_value[iEl];
    }
  } else {
    // Matrix is row-wise
    if (num_new_nz) {
      // Adding a positive number of nonzeros
      this->index_.resize(new_num_nz);
      this->value_.resize(new_num_nz);
      // Determine the row lengths of the new columns being added
      std::vector<HighsInt> new_row_length;
      new_row_length.assign(num_row, 0);
      for (HighsInt iEl = 0; iEl < num_new_nz; iEl++)
        new_row_length[new_matrix_index[iEl]]++;
      // Now shift the indices and values to make space
      HighsInt entry_offset = num_new_nz;
      HighsInt to_original_el = this->start_[num_row];
      this->start_[num_row] = new_num_nz;
      for (HighsInt iRow = num_row - 1; iRow >= 0; iRow--) {
        entry_offset -= new_row_length[iRow];
        HighsInt from_original_el = this->start_[iRow];
        // Can now use this new_row_length to store the start for the
        // new entries
        new_row_length[iRow] = to_original_el + entry_offset;
        for (HighsInt iEl = to_original_el - 1; iEl >= from_original_el;
             iEl--) {
          this->index_[iEl + entry_offset] = this->index_[iEl];
          this->value_[iEl + entry_offset] = this->value_[iEl];
        }
        to_original_el = from_original_el;
        this->start_[iRow] = entry_offset + from_original_el;
      }
      // Now insert the indices and values for the new columns
      for (HighsInt iCol = 0; iCol < num_new_col; iCol++) {
        for (HighsInt iEl = new_matrix_start[iCol];
             iEl < new_matrix_start[iCol + 1]; iEl++) {
          HighsInt iRow = new_matrix_index[iEl];
          this->index_[new_row_length[iRow]] = num_col + iCol;
          this->value_[new_row_length[iRow]] = new_matrix_value[iEl];
          new_row_length[iRow]++;
        }
      }
    }
    // Have to increase the number of columns, even if no nonzeros are being
    // added
    this->num_col_ += num_new_col;
  }
}

void HighsSparseMatrix::addRows(const HighsSparseMatrix& new_rows,
                                const int8_t* in_partition) {
  assert(new_rows.isRowwise());
  const HighsInt num_new_row = new_rows.num_row_;
  const HighsInt num_new_nz = new_rows.numNz();
  const vector<HighsInt>& new_matrix_start = new_rows.start_;
  const vector<HighsInt>& new_matrix_index = new_rows.index_;
  const vector<double>& new_matrix_value = new_rows.value_;

  assert(this->formatOk());
  // Adding rows to a row-wise partitioned matrix needs the
  // partition information
  const bool partitioned = this->format_ == MatrixFormat::kRowwisePartitioned;
  if (partitioned) {
    assert(1 == 0);
    assert(in_partition != NULL);
  }
  assert(num_new_row >= 0);
  assert(num_new_nz >= 0);
  if (num_new_row == 0) {
    // No rows are being added, so check that no nonzeros are being
    // added
    assert(num_new_nz == 0);
    return;
  }
  // Adding a positive number of rows to a matrix
  if (num_new_nz) {
    // Nonzeros are being added, so ensure that non-null data are
    // being passed
    assert(!new_matrix_start.empty());
    assert(!new_matrix_index.empty());
    assert(!new_matrix_value.empty());
  }
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  HighsInt num_nz = this->numNz();
  // Check that nonzeros aren't being appended to a matrix with no columns
  assert(num_new_nz <= 0 || num_col > 0);

  if (this->isColwise()) {
    // Matrix is currently a standard col-wise matrix, so flip
    // row-wise if there are more new nonzeros than current nonzeros
    if (num_new_nz > num_nz) this->ensureRowwise();
  }
  // Determine the new number of rows and nonzeros in the matrix
  HighsInt new_num_nz = num_nz + num_new_nz;
  HighsInt new_num_row = num_row + num_new_row;

  if (this->isRowwise()) {
    // Matrix is row-wise
    this->start_.resize(new_num_row + 1);
    // Append the starts of the new rows
    if (num_new_nz) {
      // Nontrivial number of nonzeros being added, so use new_matrix_start
      for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++)
        this->start_[num_row + iNewRow] = num_nz + new_matrix_start[iNewRow];
    } else {
      // No nonzeros being added, so new_matrix_start may be NULL, but entries
      // of zero are implied.
      for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++)
        this->start_[num_row + iNewRow] = num_nz;
    }
    this->start_[new_num_row] = new_num_nz;
    if (num_new_nz > 0) {
      // Adding a non-trivial matrix: resize the matrix arrays accordingly
      this->index_.resize(new_num_nz);
      this->value_.resize(new_num_nz);
      // Copy in the new indices and values
      if (partitioned) {
        // Insert the entries in the partition
        assert(1 == 0);
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
          HighsInt iRow = num_row + iNewRow;
          for (HighsInt iNewEl = new_matrix_start[iNewRow];
               iNewEl < new_matrix_start[iNewRow + 1]; iNewEl++) {
            HighsInt iCol = new_matrix_index[iNewEl];
            if (in_partition[iCol]) {
              HighsInt iEl = this->start_[iRow];
              this->index_[iEl] = new_matrix_index[iNewEl];
              this->value_[iEl] = new_matrix_value[iNewEl];
              this->start_[iRow]++;
            }
          }
        }
        // Use the incremented starts to initialise p_end, save these
        // values and reset the starts
        vector<HighsInt> save_p_end;
        save_p_end.resize(num_new_row);
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
          HighsInt iRow = num_row + iNewRow;
          this->start_[iRow] = num_nz + new_matrix_start[iNewRow];
          this->p_end_[iRow] = this->start_[iRow];
          save_p_end[iNewRow] = this->p_end_[iRow];
        }
        // Insert the entries not in the partition
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
          HighsInt iRow = num_row + iNewRow;
          for (HighsInt iNewEl = new_matrix_start[iNewRow];
               iNewEl < new_matrix_start[iNewRow + 1]; iNewEl++) {
            HighsInt iCol = new_matrix_index[iNewEl];
            if (!in_partition[iCol]) {
              HighsInt iEl = this->p_end_[iRow];
              this->index_[iEl] = new_matrix_index[iNewEl];
              this->value_[iEl] = new_matrix_value[iNewEl];
              this->p_end_[iRow]++;
            }
          }
        }
        // Reset p_end using the saved values
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++)
          this->p_end_[num_row + iNewRow] = save_p_end[iNewRow];
      } else {
        for (HighsInt iNewEl = 0; iNewEl < num_new_nz; iNewEl++) {
          this->index_[num_nz + iNewEl] = new_matrix_index[iNewEl];
          this->value_[num_nz + iNewEl] = new_matrix_value[iNewEl];
        }
      }
    }
  } else {
    // Storing the matrix column-wise, so have to insert the new rows
    assert(this->isColwise());
    if (num_new_nz) {
      vector<HighsInt> length;
      length.assign(num_col, 0);
      for (HighsInt iEl = 0; iEl < num_new_nz; iEl++)
        length[new_matrix_index[iEl]]++;
      // Determine the new number of nonzeros and resize the column-wise matrix
      // arrays
      this->index_.resize(new_num_nz);
      this->value_.resize(new_num_nz);
      // Append the new rows
      // Shift the existing columns to make space for the new entries
      HighsInt new_iEl = new_num_nz;
      for (HighsInt iCol = num_col - 1; iCol >= 0; iCol--) {
        HighsInt start_col_plus_1 = new_iEl;
        new_iEl -= length[iCol];
        for (HighsInt iEl = this->start_[iCol + 1] - 1;
             iEl >= this->start_[iCol]; iEl--) {
          new_iEl--;
          this->index_[new_iEl] = this->index_[iEl];
          this->value_[new_iEl] = this->value_[iEl];
        }
        this->start_[iCol + 1] = start_col_plus_1;
      }
      assert(new_iEl == 0);
      // Insert the new entries
      for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
        HighsInt first_el = new_matrix_start[iNewRow];
        HighsInt last_el =
            (iNewRow < num_new_row - 1 ? new_matrix_start[iNewRow + 1]
                                       : num_new_nz);
        for (HighsInt iEl = first_el; iEl < last_el; iEl++) {
          HighsInt iCol = new_matrix_index[iEl];
          new_iEl = this->start_[iCol + 1] - length[iCol];
          length[iCol]--;
          this->index_[new_iEl] = num_row + iNewRow;
          this->value_[new_iEl] = new_matrix_value[iEl];
        }
      }
    }
  }
  // Update the number of rows
  this->num_row_ += num_new_row;
}

void HighsSparseMatrix::getCol(const HighsInt iCol, HighsInt& num_nz,
                               HighsInt* index, double* value) const {
  assert(iCol >= 0 && iCol < this->num_col_);
  num_nz = 0;
  if (this->isColwise()) {
    for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
         iEl++) {
      index[num_nz] = this->index_[iEl];
      value[num_nz] = this->value_[iEl];
      num_nz++;
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        if (this->index_[iEl] == iCol) {
          index[num_nz] = iRow;
          value[num_nz] = this->value_[iEl];
          num_nz++;
          break;
        }
      }
    }
  }
}

void HighsSparseMatrix::getRow(const HighsInt iRow, HighsInt& num_nz,
                               HighsInt* index, double* value) const {
  assert(iRow >= 0 && iRow < this->num_row_);
  num_nz = 0;
  if (this->isRowwise()) {
    for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
         iEl++) {
      index[num_nz] = this->index_[iEl];
      value[num_nz] = this->value_[iEl];
      num_nz++;
    }
  } else {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        if (this->index_[iEl] == iRow) {
          index[num_nz] = iCol;
          value[num_nz] = this->value_[iEl];
          num_nz++;
          break;
        }
      }
    }
  }
}

void HighsSparseMatrix::deleteCols(
    const HighsIndexCollection& index_collection) {
  assert(this->formatOk());
  // Can't handle rowwise matrices yet
  assert(!this->isRowwise());
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

  HighsInt col_dim = this->num_col_;
  HighsInt new_num_col = 0;
  HighsInt new_num_nz = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, delete_from_col, delete_to_col,
                     keep_from_col, keep_to_col, current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      new_num_col = delete_from_col;
      new_num_nz = this->start_[delete_from_col];
    }
    // Ensure that the starts of the deleted columns are zeroed to
    // avoid redundant start information for columns whose indices
    // aren't used after the deletion takes place. In particular, if
    // all columns are deleted then something must be done to ensure
    // that the matrix isn't magically recreated by increasing the
    // number of columns from zero when there are no rows in the
    // matrix.
    for (HighsInt col = delete_from_col; col <= delete_to_col; col++)
      this->start_[col] = 0;
    // Shift the starts - both in place and value - to account for the
    // columns and nonzeros removed
    const HighsInt keep_from_el = this->start_[keep_from_col];
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      this->start_[new_num_col] = new_num_nz + this->start_[col] - keep_from_el;
      new_num_col++;
    }
    for (HighsInt el = keep_from_el; el < this->start_[keep_to_col + 1]; el++) {
      this->index_[new_num_nz] = this->index_[el];
      this->value_[new_num_nz] = this->value_[el];
      new_num_nz++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  // Ensure that the start of the spurious last column is zeroed so
  // that it doesn't give a positive number of matrix entries if the
  // number of columns in the matrix is increased when there are no
  // rows in the matrix.
  this->start_[this->num_col_] = 0;
  this->start_[new_num_col] = new_num_nz;
  this->start_.resize(new_num_col + 1);
  this->index_.resize(new_num_nz);
  this->value_.resize(new_num_nz);
  // Update the number of columns
  this->num_col_ = new_num_col;
}

void HighsSparseMatrix::deleteRows(
    const HighsIndexCollection& index_collection) {
  assert(this->formatOk());
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k > to_k) return;

  HighsInt delete_from_row;
  HighsInt delete_to_row;
  HighsInt keep_from_row;
  HighsInt row_dim = this->num_row_;
  HighsInt keep_to_row = -1;
  HighsInt current_set_entry = 0;

  // Set up a row mask to indicate the new row index of kept rows and
  // -1 for deleted rows so that the kept entries in the column-wise
  // matrix can be identified and have their correct row index.
  vector<HighsInt> new_index;
  new_index.resize(this->num_row_);
  HighsInt new_num_row = 0;
  bool mask = index_collection.is_mask_;
  const vector<HighsInt>& row_mask = index_collection.mask_;
  if (!mask) {
    keep_to_row = -1;
    current_set_entry = 0;
    for (HighsInt k = from_k; k <= to_k; k++) {
      updateOutInIndex(index_collection, delete_from_row, delete_to_row,
                       keep_from_row, keep_to_row, current_set_entry);
      if (k == from_k) {
        // Account for any initial rows being kept
        for (HighsInt row = 0; row < delete_from_row; row++) {
          new_index[row] = new_num_row;
          new_num_row++;
        }
      }
      for (HighsInt row = delete_from_row; row <= delete_to_row; row++) {
        new_index[row] = -1;
      }
      for (HighsInt row = keep_from_row; row <= keep_to_row; row++) {
        new_index[row] = new_num_row;
        new_num_row++;
      }
      if (keep_to_row >= row_dim - 1) break;
    }
  } else {
    for (HighsInt row = 0; row < this->num_row_; row++) {
      if (row_mask[row]) {
        new_index[row] = -1;
      } else {
        new_index[row] = new_num_row;
        new_num_row++;
      }
    }
  }
  HighsInt new_num_nz = 0;
  for (HighsInt col = 0; col < this->num_col_; col++) {
    HighsInt from_el = this->start_[col];
    this->start_[col] = new_num_nz;
    for (HighsInt el = from_el; el < this->start_[col + 1]; el++) {
      HighsInt row = this->index_[el];
      HighsInt new_row = new_index[row];
      if (new_row >= 0) {
        this->index_[new_num_nz] = new_row;
        this->value_[new_num_nz] = this->value_[el];
        new_num_nz++;
      }
    }
  }
  this->start_[this->num_col_] = new_num_nz;
  this->start_.resize(this->num_col_ + 1);
  this->index_.resize(new_num_nz);
  this->value_.resize(new_num_nz);
  // Update the number of rows
  this->num_row_ = new_num_row;
}

HighsStatus HighsSparseMatrix::assessStart(const HighsLogOptions& log_options) {
  // Identify main dimensions
  const HighsInt vec_dim = this->vecDim();
  const HighsInt num_vec = this->numVec();
  if (this->start_[0]) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Matrix start[0] = %d, not 0\n", int(this->start_[0]));
    return HighsStatus::kError;
  }
  HighsInt num_nz = this->numNz();
  for (HighsInt iVec = 1; iVec < num_vec; iVec++) {
    if (this->start_[iVec] < this->start_[iVec - 1]) {
      highsLogUser(log_options, HighsLogType::kError,
                   "Matrix start[%d] = %d > %d = start[%d]\n", int(iVec),
                   int(this->start_[iVec]), int(this->start_[iVec - 1]),
                   int(iVec - 1));
      return HighsStatus::kError;
    }
    if (this->start_[iVec] > num_nz) {
      highsLogUser(log_options, HighsLogType::kError,
                   "Matrix start[%d] = %d > %d = number of nonzeros\n",
                   int(iVec), int(this->start_[iVec]), int(num_nz));
      return HighsStatus::kError;
    }
  }
  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::assessIndexBounds(
    const HighsLogOptions& log_options) {
  // Identify main dimensions
  const HighsInt vec_dim = this->vecDim();
  for (HighsInt iEl = 0; iEl < this->numNz(); iEl++) {
    if (this->index_[iEl] < 0 || this->index_[iEl] >= vec_dim) {
      highsLogUser(log_options, HighsLogType::kError,
                   "Matrix index[%d] = %d is not in legal range of [0, %d)\n",
                   int(iEl), int(this->index_[iEl]), vec_dim);
      return HighsStatus::kError;
    }
  }
  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::assess(const HighsOptions& options,
                                      const bool structure_only) {
  const bool sum_duplicates = false;
  const double small_matrix_value =
      structure_only ? 0 : options.small_matrix_value;
  const double large_matrix_value =
      structure_only ? kHighsInf : options.large_matrix_value;
  return this->assess(options.log_options, small_matrix_value,
                      large_matrix_value, sum_duplicates);
}

HighsStatus HighsSparseMatrix::assess(const HighsLogOptions& log_options,
                                      const double small_matrix_value,
                                      const double large_matrix_value,
                                      const bool sum_duplicates,
                                      const std::string* col_names,
                                      const std::string* row_names) {
  assert(this->formatOk());
  return assessMatrix(log_options, "HighsSparseMatrix", this->name_,
                      this->vecDim(), this->numVec(),
                      this->format_ == MatrixFormat::kRowwisePartitioned,
                      this->start_, this->p_end_, this->index_, this->value_,
                      small_matrix_value, large_matrix_value, sum_duplicates,
                      col_names, row_names);
}

HighsStatus HighsSparseMatrix::assessExtremeValues(
    const std::string& name, const HighsOptions& options) const {
  // Returns HighsStatus::kError if the matrix contains zeros or large
  // values, and HighsStatus::kWarning if the matrix contains small
  // nonzero values, otherwise returns HighsStatus::kOk
  double min_abs_value;
  double max_abs_value;
  if (!this->hasExtremeValues(options.small_matrix_value,
                              options.large_matrix_value, min_abs_value,
                              max_abs_value))
    return HighsStatus::kOk;
  bool error_found = false;
  bool warning_found = false;
  if (min_abs_value <= options.small_matrix_value) {
    warning_found = true;
    error_found = min_abs_value == 0;
    highsLogUser(options.log_options,
                 error_found ? HighsLogType::kError : HighsLogType::kWarning,
                 "Matrix %shas smallest |value| of %.2g at/below "
                 "tolerance of %g\n",
                 name.c_str(), min_abs_value, options.small_matrix_value);
  }
  if (max_abs_value >= options.large_matrix_value) {
    error_found = true;
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Matrix %shas largest |value| of %.2g at/above "
                 "tolerance of %g\n",
                 name.c_str(), max_abs_value, options.large_matrix_value);
  }
  const HighsInt num_values = this->value_.size();
  analyseVectorValues(&options.log_options, "Extreme values in matrix",
                      num_values, this->value_, false, "");
  assert(warning_found || error_found);
  return error_found ? HighsStatus::kError : HighsStatus::kWarning;
}

bool HighsSparseMatrix::hasExtremeValues(const double small_matrix_value,
                                         const double large_matrix_value,
                                         double& min_abs_value,
                                         double& max_abs_value) const {
  min_abs_value = kHighsInf;
  max_abs_value = 0;
  const HighsInt num_values = this->value_.size();
  for (HighsInt iX = 0; iX < num_values; iX++) {
    min_abs_value = std::min(std::abs(this->value_[iX]), min_abs_value);
    max_abs_value = std::max(std::abs(this->value_[iX]), max_abs_value);
  }
  if (min_abs_value > small_matrix_value && max_abs_value < large_matrix_value)
    return false;
  return true;
}

void HighsSparseMatrix::considerColScaling(
    const HighsInt max_scale_factor_exponent, double* col_scale) {
  const double log2 = log(2.0);
  const double max_allow_scale = pow(2.0, max_scale_factor_exponent);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_col_scale = min_allow_scale;
  const double max_allow_col_scale = max_allow_scale;

  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      double col_max_value = 0;
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        col_max_value = max(fabs(this->value_[iEl]), col_max_value);
      if (col_max_value) {
        double col_scale_value = 1 / col_max_value;
        // Convert the col scale factor to the nearest power of two, and
        // ensure that it is not excessively large or small
        col_scale_value = pow(2.0, floor(log(col_scale_value) / log2 + 0.5));
        col_scale_value =
            min(max(min_allow_col_scale, col_scale_value), max_allow_col_scale);
        col_scale[iCol] = col_scale_value;
        // Scale the column
        for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
             iEl++)
          this->value_[iEl] *= col_scale[iCol];
      } else {
        // Empty column
        col_scale[iCol] = 1;
      }
    }
  } else {
    assert(1 == 0);
  }
}

void HighsSparseMatrix::considerRowScaling(
    const HighsInt max_scale_factor_exponent, double* row_scale) {
  const double log2 = log(2.0);
  const double max_allow_scale = pow(2.0, max_scale_factor_exponent);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_row_scale = min_allow_scale;
  const double max_allow_row_scale = max_allow_scale;

  if (this->isRowwise()) {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      double row_max_value = 0;
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        row_max_value = max(fabs(this->value_[iEl]), row_max_value);
      if (row_max_value) {
        double row_scale_value = 1 / row_max_value;
        // Convert the row scale factor to the nearest power of two, and
        // ensure that it is not excessively large or small
        row_scale_value = pow(2.0, floor(log(row_scale_value) / log2 + 0.5));
        row_scale_value =
            min(max(min_allow_row_scale, row_scale_value), max_allow_row_scale);
        row_scale[iRow] = row_scale_value;
        // Scale the row
        for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
             iEl++)
          this->value_[iEl] *= row_scale[iRow];
      } else {
        // Empty row
        row_scale[iRow] = 1;
      }
    }
  } else {
    assert(1 == 0);
  }
}

void HighsSparseMatrix::scaleCol(const HighsInt col, const double colScale) {
  assert(this->formatOk());
  assert(col >= 0);
  assert(col < this->num_col_);
  assert(colScale);

  if (this->isColwise()) {
    for (HighsInt iEl = this->start_[col]; iEl < this->start_[col + 1]; iEl++)
      this->value_[iEl] *= colScale;
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        if (this->index_[iEl] == col) this->value_[iEl] *= colScale;
      }
    }
  }
}

void HighsSparseMatrix::scaleRow(const HighsInt row, const double rowScale) {
  assert(this->formatOk());
  assert(row >= 0);
  assert(row < this->num_row_);
  assert(rowScale);

  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        if (this->index_[iEl] == row) this->value_[iEl] *= rowScale;
      }
    }
  } else {
    for (HighsInt iEl = this->start_[row]; iEl < this->start_[row + 1]; iEl++)
      this->value_[iEl] *= rowScale;
  }
}

void HighsSparseMatrix::applyScale(const HighsScale& scale) {
  assert(this->formatOk());
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        HighsInt iRow = this->index_[iEl];
        this->value_[iEl] *= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        HighsInt iCol = this->index_[iEl];
        this->value_[iEl] *= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  }
}

void HighsSparseMatrix::applyColScale(const HighsScale& scale) {
  assert(this->formatOk());
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        this->value_[iEl] *= scale.col[iCol];
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        this->value_[iEl] *= scale.col[this->index_[iEl]];
    }
  }
}

void HighsSparseMatrix::applyRowScale(const HighsScale& scale) {
  assert(this->formatOk());
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        this->value_[iEl] *= scale.row[this->index_[iEl]];
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        this->value_[iEl] *= scale.row[iRow];
    }
  }
}

void HighsSparseMatrix::unapplyScale(const HighsScale& scale) {
  assert(this->formatOk());
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        HighsInt iRow = this->index_[iEl];
        this->value_[iEl] /= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        HighsInt iCol = this->index_[iEl];
        this->value_[iEl] /= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  }
}

void HighsSparseMatrix::createSlice(const HighsSparseMatrix& matrix,
                                    const HighsInt from_col,
                                    const HighsInt to_col) {
  assert(matrix.formatOk());
  assert(matrix.isColwise());
  assert(this->formatOk());
  HighsInt num_row = matrix.num_row_;
  const vector<HighsInt>& a_start = matrix.start_;
  const vector<HighsInt>& a_index = matrix.index_;
  const vector<double>& a_value = matrix.value_;
  vector<HighsInt>& slice_start = this->start_;
  vector<HighsInt>& slice_index = this->index_;
  vector<double>& slice_value = this->value_;
  HighsInt slice_num_col = to_col + 1 - from_col;
  HighsInt slice_num_nz = a_start[to_col + 1] - a_start[from_col];
  slice_start.resize(slice_num_col + 1);
  slice_index.resize(slice_num_nz);
  slice_value.resize(slice_num_nz);
  HighsInt from_col_start = a_start[from_col];
  for (HighsInt iCol = from_col; iCol < to_col + 1; iCol++)
    slice_start[iCol - from_col] = a_start[iCol] - from_col_start;
  slice_start[slice_num_col] = slice_num_nz;
  for (HighsInt iEl = a_start[from_col]; iEl < a_start[to_col + 1]; iEl++) {
    slice_index[iEl - from_col_start] = a_index[iEl];
    slice_value[iEl - from_col_start] = a_value[iEl];
  }
  this->num_col_ = slice_num_col;
  this->num_row_ = num_row;
  this->format_ = MatrixFormat::kColwise;
}

void HighsSparseMatrix::createRowwise(const HighsSparseMatrix& matrix) {
  assert(matrix.formatOk());
  assert(matrix.isColwise());
  assert(this->formatOk());

  HighsInt num_col = matrix.num_col_;
  HighsInt num_row = matrix.num_row_;
  HighsInt num_nz = matrix.numNz();
  const vector<HighsInt>& a_start = matrix.start_;
  const vector<HighsInt>& a_index = matrix.index_;
  const vector<double>& a_value = matrix.value_;
  vector<HighsInt>& ar_start = this->start_;
  vector<HighsInt>& ar_index = this->index_;
  vector<double>& ar_value = this->value_;

  // Use ar_end to compute lengths, which are then transformed into
  // the ends of the inserted entries
  std::vector<HighsInt> ar_end;
  ar_start.resize(num_row + 1);
  ar_end.assign(num_row, 0);
  // Count the nonzeros in each row
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
      HighsInt iRow = a_index[iEl];
      ar_end[iRow]++;
    }
  }
  // Compute the starts and turn the lengths into ends
  ar_start[0] = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    ar_start[iRow + 1] = ar_start[iRow] + ar_end[iRow];
    ar_end[iRow] = ar_start[iRow];
  }
  ar_index.resize(num_nz);
  ar_value.resize(num_nz);
  // Insert the entries
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
      HighsInt iRow = a_index[iEl];
      HighsInt iToEl = ar_end[iRow]++;
      ar_index[iToEl] = iCol;
      ar_value[iToEl] = a_value[iEl];
    }
  }
  this->format_ = MatrixFormat::kRowwise;
  this->num_col_ = num_col;
  this->num_row_ = num_row;
}

void HighsSparseMatrix::createColwise(const HighsSparseMatrix& matrix) {
  assert(matrix.formatOk());
  assert(matrix.isRowwise());
  assert(this->formatOk());

  HighsInt num_col = matrix.num_col_;
  HighsInt num_row = matrix.num_row_;
  HighsInt num_nz = matrix.numNz();
  const vector<HighsInt>& ar_start = matrix.start_;
  const vector<HighsInt>& ar_index = matrix.index_;
  const vector<double>& ar_value = matrix.value_;
  vector<HighsInt>& a_start = this->start_;
  vector<HighsInt>& a_index = this->index_;
  vector<double>& a_value = this->value_;

  // Use a_end to compute lengths, which are then transformed into
  // the ends of the inserted entries
  std::vector<HighsInt> a_end;
  a_start.resize(num_col + 1);
  a_end.assign(num_col, 0);
  // Count the nonzeros in each col
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    for (HighsInt iEl = ar_start[iRow]; iEl < ar_start[iRow + 1]; iEl++) {
      HighsInt iCol = ar_index[iEl];
      a_end[iCol]++;
    }
  }
  // Compute the starts and turn the lengths into ends
  a_start[0] = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    a_start[iCol + 1] = a_start[iCol] + a_end[iCol];
    a_end[iCol] = a_start[iCol];
  }
  a_index.resize(num_nz);
  a_value.resize(num_nz);
  // Insert the entries
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    for (HighsInt iEl = ar_start[iRow]; iEl < ar_start[iRow + 1]; iEl++) {
      HighsInt iCol = ar_index[iEl];
      HighsInt iToEl = a_end[iCol]++;
      a_index[iToEl] = iRow;
      a_value[iToEl] = ar_value[iEl];
    }
  }
  this->format_ = MatrixFormat::kColwise;
  this->num_col_ = num_col;
  this->num_row_ = num_row;
}

void HighsSparseMatrix::alphaProductPlusY(const double alpha,
                                          const std::vector<double>& x,
                                          std::vector<double>& y,
                                          const bool transpose) const {
  assert(x.size() >=
         static_cast<size_t>(transpose ? this->num_row_ : this->num_col_));
  assert(y.size() >=
         static_cast<size_t>(transpose ? this->num_col_ : this->num_row_));
  if (this->isColwise()) {
    if (transpose) {
      for (int iCol = 0; iCol < this->num_col_; iCol++)
        for (int iEl = this->start_[iCol]; iEl < this->start_[iCol + 1]; iEl++)
          y[iCol] += alpha * this->value_[iEl] * x[this->index_[iEl]];
    } else {
      for (int iCol = 0; iCol < this->num_col_; iCol++)
        for (int iEl = this->start_[iCol]; iEl < this->start_[iCol + 1]; iEl++)
          y[this->index_[iEl]] += alpha * this->value_[iEl] * x[iCol];
    }
  } else {
    if (transpose) {
      for (int iRow = 0; iRow < this->num_row_; iRow++)
        for (int iEl = this->start_[iRow]; iEl < this->start_[iRow + 1]; iEl++)
          y[this->index_[iEl]] += alpha * this->value_[iEl] * x[iRow];
    } else {
      for (int iRow = 0; iRow < this->num_row_; iRow++)
        for (int iEl = this->start_[iRow]; iEl < this->start_[iRow + 1]; iEl++)
          y[iRow] += alpha * this->value_[iEl] * x[this->index_[iEl]];
    }
  }
}

void HighsSparseMatrix::product(vector<double>& result,
                                const vector<double>& row) const {
  assert(this->formatOk());
  assert((int)row.size() >= this->num_col_);
  result.assign(this->num_row_, 0.0);
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        result[this->index_[iEl]] += row[iCol] * this->value_[iEl];
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        result[iRow] += row[this->index_[iEl]] * this->value_[iEl];
    }
  }
}

void HighsSparseMatrix::productTranspose(vector<double>& result,
                                         const vector<double>& col) const {
  assert(this->formatOk());
  assert((int)col.size() >= this->num_row_);
  result.assign(this->num_col_, 0.0);
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        result[iCol] += col[this->index_[iEl]] * this->value_[iEl];
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        result[this->index_[iEl]] += col[iRow] * this->value_[iEl];
    }
  }
}

void HighsSparseMatrix::productQuad(vector<double>& result,
                                    const vector<double>& row,
                                    const HighsInt debug_report) const {
  assert(this->formatOk());
  assert((int)row.size() >= this->num_col_);
  result.assign(this->num_row_, 0.0);
  if (this->isColwise()) {
    std::vector<HighsCDouble> value(this->num_row_, 0);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        value[this->index_[iEl]] += row[iCol] * this->value_[iEl];
    }
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++)
      result[iRow] = double(value[iRow]);
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      HighsCDouble value = 0.0;
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        value += row[this->index_[iEl]] * this->value_[iEl];
      result[iRow] = double(value);
    }
  }
}

void HighsSparseMatrix::productTransposeQuad(
    vector<double>& result, const vector<double>& row,
    const HighsInt debug_report) const {
  assert(this->formatOk());
  assert((int)row.size() >= this->num_row_);
  result.assign(this->num_col_, 0.0);
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      HighsCDouble value = 0.0;
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        value += row[this->index_[iEl]] * this->value_[iEl];
      result[iCol] = double(value);
    }
  } else {
    std::vector<HighsCDouble> value(this->num_col_, 0);
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        value[this->index_[iEl]] += row[iRow] * this->value_[iEl];
    }
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
      result[iCol] = double(value[iCol]);
  }
}

void HighsSparseMatrix::productTransposeQuad(
    vector<double>& result_value, vector<HighsInt>& result_index,
    const HVector& column, const HighsInt debug_report) const {
  if (debug_report >= kDebugReportAll)
    printf("\nHighsSparseMatrix::productTranspose:\n");
  if (this->isColwise()) {
    result_value.reserve(this->num_col_);
    result_index.reserve(this->num_col_);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      HighsCDouble value = 0.0;
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        value += column.array[this->index_[iEl]] * this->value_[iEl];

      if (abs(value) - kHighsTiny > 0.0) {
        result_value.push_back(double(value));
        result_index.push_back(iCol);
      }
    }
  } else {
    HighsSparseVectorSum sum(num_col_);
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      double multiplier = column.array[iRow];
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++)
        sum.add(this->index_[iEl], multiplier * this->value_[iEl]);
    }
    if (debug_report >= kDebugReportAll) {
      HighsSparseVectorSum report_sum(num_col_);
      for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
        double multiplier = column.array[iRow];
        if (debug_report == kDebugReportAll || debug_report == iRow)
          debugReportRowPrice(iRow, multiplier, this->start_[iRow + 1],
                              report_sum);
      }
    }

    sum.cleanup([](HighsInt, double x) { return std::abs(x) <= kHighsTiny; });
    result_index = std::move(sum.nonzeroinds);
    HighsInt result_num_nz = static_cast<HighsInt>(result_index.size());
    result_value.reserve(result_num_nz);
    for (HighsInt i = 0; i < result_num_nz; ++i)
      result_value.push_back(sum.getValue(result_index[i]));
  }
}

void HighsSparseMatrix::createRowwisePartitioned(
    const HighsSparseMatrix& matrix, const int8_t* in_partition) {
  assert(matrix.formatOk());
  assert(matrix.isColwise());
  assert(this->formatOk());
  const bool all_in_partition = in_partition == NULL;

  HighsInt num_col = matrix.num_col_;
  HighsInt num_row = matrix.num_row_;
  HighsInt num_nz = matrix.numNz();
  const vector<HighsInt>& a_start = matrix.start_;
  const vector<HighsInt>& a_index = matrix.index_;
  const vector<double>& a_value = matrix.value_;
  vector<HighsInt>& ar_start = this->start_;
  vector<HighsInt>& ar_p_end = this->p_end_;
  vector<HighsInt>& ar_index = this->index_;
  vector<double>& ar_value = this->value_;

  // Use ar_p_end and ar_end to compute lengths, which are then transformed into
  // the p_ends and ends of the inserted entries
  std::vector<HighsInt> ar_end;
  ar_start.resize(num_row + 1);
  ar_p_end.assign(num_row, 0);
  ar_end.assign(num_row, 0);
  // Count the nonzeros of nonbasic and basic columns in each row
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (all_in_partition || in_partition[iCol]) {
      for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
        HighsInt iRow = a_index[iEl];
        ar_p_end[iRow]++;
      }
    } else {
      for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
        HighsInt iRow = a_index[iEl];
        ar_end[iRow]++;
      }
    }
  }
  // Compute the starts and turn the lengths into ends
  ar_start[0] = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    ar_start[iRow + 1] = ar_start[iRow] + ar_p_end[iRow] + ar_end[iRow];
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    ar_end[iRow] = ar_start[iRow] + ar_p_end[iRow];
    ar_p_end[iRow] = ar_start[iRow];
  }
  // Insert the entries
  ar_index.resize(num_nz);
  ar_value.resize(num_nz);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (all_in_partition || in_partition[iCol]) {
      for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
        HighsInt iRow = a_index[iEl];
        HighsInt iToEl = ar_p_end[iRow]++;
        ar_index[iToEl] = iCol;
        ar_value[iToEl] = a_value[iEl];
      }
    } else {
      for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
        HighsInt iRow = a_index[iEl];
        HighsInt iToEl = ar_end[iRow]++;
        ar_index[iToEl] = iCol;
        ar_value[iToEl] = a_value[iEl];
      }
    }
  }
  this->format_ = MatrixFormat::kRowwisePartitioned;
  this->num_col_ = num_col;
  this->num_row_ = num_row;
}

bool HighsSparseMatrix::debugPartitionOk(const int8_t* in_partition) const {
  assert(this->format_ == MatrixFormat::kRowwisePartitioned);
  bool ok = true;
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    for (HighsInt iEl = this->start_[iRow]; iEl < this->p_end_[iRow]; iEl++) {
      if (!in_partition[this->index_[iEl]]) {
        ok = false;
        break;
      }
    }
    if (!ok) break;
    for (HighsInt iEl = this->p_end_[iRow]; iEl < this->start_[iRow + 1];
         iEl++) {
      if (in_partition[this->index_[iEl]]) {
        ok = false;
        break;
      }
    }
    if (!ok) break;
  }
  return ok;
}

void HighsSparseMatrix::priceByColumn(const bool quad_precision,
                                      HVector& result, const HVector& column,
                                      const HighsInt debug_report) const {
  assert(this->isColwise());
  if (debug_report >= kDebugReportAll)
    printf("\nHighsSparseMatrix::priceByColumn:\n");
  result.count = 0;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    double value = 0;
    if (quad_precision) {
      HighsCDouble quad_value = 0.0;
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        quad_value += column.array[this->index_[iEl]] * this->value_[iEl];
      value = (double)quad_value;
    } else {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++)
        value += column.array[this->index_[iEl]] * this->value_[iEl];
    }
    if (fabs(value) > kHighsTiny) {
      result.array[iCol] = value;
      result.index[result.count++] = iCol;
    }
  }
}

void HighsSparseMatrix::priceByRow(const bool quad_precision, HVector& result,
                                   const HVector& column,
                                   const HighsInt debug_report) const {
  assert(this->isRowwise());
  if (debug_report >= kDebugReportAll)
    printf("\nHighsSparseMatrix::priceByRow:\n");
  // Vanilla hyper-sparse row-wise PRICE. Set up parameters so that
  // priceByRowWithSwitch runs as vanilla hyper-sparse PRICE
  // Expected density always forces hyper-sparse PRICE
  const double expected_density = -kHighsInf;
  // Always start from first index of column
  HighsInt from_index = 0;
  // Never switch to standard row-wise PRICE
  const double switch_density = kHighsInf;
  this->priceByRowWithSwitch(quad_precision, result, column, expected_density,
                             from_index, switch_density);
}

void HighsSparseMatrix::priceByRowWithSwitch(
    const bool quad_precision, HVector& result, const HVector& column,
    const double expected_density, const HighsInt from_index,
    const double switch_density, const HighsInt debug_report) const {
  assert(this->isRowwise());
  HighsSparseVectorSum sum;
  // todo @Julian: Setting up the sparse vector sum is equivalent to calling
  // HVector::setup() I think there should instead be overloads where the result
  // vector is of type HVectorQuad for the future and not the boolean parameter
  // quad_precision. Then the buffer can be maintained similar to row_ap.
  if (quad_precision) sum.setDimension(num_col_);
  if (debug_report >= kDebugReportAll)
    printf("\nHighsSparseMatrix::priceByRowWithSwitch\n");
  // (Continue) hyper-sparse row-wise PRICE with possible switches to
  // standard row-wise PRICE either immediately based on historical
  // density or during hyper-sparse PRICE if there is too much fill-in
  HighsInt next_index = from_index;
  // Possibly don't perform hyper-sparse PRICE based on historical density
  //
  // Ensure that result was set up for this number of columns, and
  // that result.index is still of correct size
  assert(HighsInt(result.size) == this->num_col_);
  assert(HighsInt(result.index.size()) == this->num_col_);
  if (expected_density <= kHyperPriceDensity) {
    double inv_num_col = 1.0 / this->num_col_;
    for (HighsInt ix = next_index; ix < column.count; ix++) {
      HighsInt iRow = column.index[ix];
      // Determine whether p_end_ or the next start_ ends the loop
      HighsInt to_iEl;
      if (this->format_ == MatrixFormat::kRowwisePartitioned) {
        to_iEl = this->p_end_[iRow];
      } else {
        to_iEl = this->start_[iRow + 1];
      }
      // Possibly switch to standard row-wise price
      HighsInt row_num_nz = to_iEl - this->start_[iRow];
      double local_density = (1.0 * result.count) * inv_num_col;
      bool switch_to_dense = result.count + row_num_nz >= this->num_col_ ||
                             local_density > switch_density;
      if (switch_to_dense) break;
      double multiplier = column.array[iRow];
      if (debug_report == kDebugReportAll || debug_report == iRow)
        debugReportRowPrice(iRow, multiplier, to_iEl, result.array);
      if (multiplier) {
        if (quad_precision) {
          for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
            sum.add(this->index_[iEl], multiplier * this->value_[iEl]);
          }
        } else {
          for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
            HighsInt iCol = this->index_[iEl];
            double value0 = result.array[iCol];
            double value1 = value0 + multiplier * this->value_[iEl];
            if (value0 == 0) result.index[result.count++] = iCol;
            result.array[iCol] =
                (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
          }
        }
      }
      next_index = ix + 1;
    }
  }
  if (quad_precision)
    sum.cleanup([](HighsInt, double x) { return std::abs(x) <= kHighsTiny; });
  if (next_index < column.count) {
    // PRICE is not complete: finish without maintaining nonzeros of result
    if (quad_precision) {
      std::vector<HighsCDouble> result_array = std::move(sum.values);
      this->priceByRowDenseResult(result_array, column, next_index);
      // Determine indices of nonzeros in result
      result.count = 0;
      for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
        double value1 = (double)result_array[iCol];
        if (fabs(value1) < kHighsTiny) {
          result.array[iCol] = 0;
        } else {
          result.array[iCol] = value1;
          result.index[result.count++] = iCol;
        }
      }
    } else {
      this->priceByRowDenseResult(result.array, column, next_index);
      // Determine indices of nonzeros in result
      result.count = 0;
      for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
        double value1 = result.array[iCol];
        if (fabs(value1) < kHighsTiny) {
          result.array[iCol] = 0;
        } else {
          result.index[result.count++] = iCol;
        }
      }
    }
  } else {
    if (quad_precision) {
      // HVector result should have result.index of size this->num_col_
      // by virtue of result.setup. However, it will generally lose
      // this property by virtue of the following move
      result.index = std::move(sum.nonzeroinds);
      HighsInt result_num_nz = static_cast<HighsInt>(result.index.size());
      // Restore the size of result.index
      result.index.resize(this->num_col_);
      result.count = result_num_nz;
      for (HighsInt i = 0; i < result_num_nz; ++i) {
        HighsInt iRow = result.index[i];
        result.array[iRow] = sum.getValue(iRow);
      }
    } else {
      // PRICE is complete maintaining nonzeros of result
      // Remove small values
      result.tight();
    }
  }
}

void HighsSparseMatrix::update(const HighsInt var_in, const HighsInt var_out,
                               const HighsSparseMatrix& matrix) {
  assert(matrix.format_ == MatrixFormat::kColwise);
  assert(this->format_ == MatrixFormat::kRowwisePartitioned);
  if (var_in < this->num_col_) {
    for (HighsInt iEl = matrix.start_[var_in]; iEl < matrix.start_[var_in + 1];
         iEl++) {
      HighsInt iRow = matrix.index_[iEl];
      HighsInt iFind = this->start_[iRow];
      HighsInt iSwap = --this->p_end_[iRow];
      while (this->index_[iFind] != var_in) iFind++;
      // todo @ Julian : this assert can fail
      assert(iFind >= 0 && iFind < int(this->index_.size()));
      assert(iSwap >= 0 && iSwap < int(this->value_.size()));
      swap(this->index_[iFind], this->index_[iSwap]);
      swap(this->value_[iFind], this->value_[iSwap]);
    }
  }

  if (var_out < this->num_col_) {
    for (HighsInt iEl = matrix.start_[var_out];
         iEl < matrix.start_[var_out + 1]; iEl++) {
      HighsInt iRow = matrix.index_[iEl];
      HighsInt iFind = this->p_end_[iRow];
      HighsInt iSwap = this->p_end_[iRow]++;
      while (this->index_[iFind] != var_out) iFind++;
      swap(this->index_[iFind], this->index_[iSwap]);
      swap(this->value_[iFind], this->value_[iSwap]);
    }
  }
}

double HighsSparseMatrix::computeDot(const std::vector<double>& array,
                                     const HighsInt use_col) const {
  assert(this->isColwise());
  double result = 0;
  if (use_col < this->num_col_) {
    for (HighsInt iEl = this->start_[use_col]; iEl < this->start_[use_col + 1];
         iEl++)
      result += array[this->index_[iEl]] * this->value_[iEl];
  } else {
    result = array[use_col - this->num_col_];
  }
  return result;
}

void HighsSparseMatrix::collectAj(HVector& column, const HighsInt use_col,
                                  const double multiplier) const {
  assert(this->isColwise());
  if (use_col < this->num_col_) {
    for (HighsInt iEl = this->start_[use_col]; iEl < this->start_[use_col + 1];
         iEl++) {
      HighsInt iRow = this->index_[iEl];
      double value0 = column.array[iRow];
      double value1 = value0 + multiplier * this->value_[iEl];
      if (value0 == 0) column.index[column.count++] = iRow;
      column.array[iRow] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
    }
  } else {
    HighsInt iRow = use_col - this->num_col_;
    double value0 = column.array[iRow];
    double value1 = value0 + multiplier;
    if (value0 == 0) column.index[column.count++] = iRow;
    column.array[iRow] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
  }
}

HighsSparseMatrix HighsSparseMatrix::oppositeFormat() const {
  // Given a matrix A, doing
  //    HighsSparseMatrix B = A.oppositeFormat();
  // is more efficient than doing
  //    HighsSparseMatrix B = A;
  //    B.ensureRowwise();
  // because the latter requires a copy (B = A) and an internal copy in
  // ensureRowWise(). Similarly for ensureColwise.

  HighsSparseMatrix T;

  const bool colwise = isColwise();
  const HighsInt num_major = colwise ? num_col_ : num_row_;
  const HighsInt num_minor = colwise ? num_row_ : num_col_;

  T.format_ = colwise ? MatrixFormat::kRowwise : MatrixFormat::kColwise;
  T.num_row_ = num_row_;
  T.num_col_ = num_col_;
  T.p_end_ = p_end_;

  highsSparseTranspose(num_minor, num_major, start_, index_, value_, T.start_,
                       T.index_, T.value_);

  return T;
}

HighsSparseMatrix HighsSparseMatrix::getRowwise() const {
  // To get a row-wise copy B of A, calling
  //    HighsSparseMatrix B = A;
  //    B.ensureRowwise();
  // makes one copy if A is already rowwise and two otherwise (one is hidden
  // inside ensureRowwise).
  // Calling
  //    HighsSparseMatrix B = A.getRowwise();
  // makes one copy if A is already rowwise; otherwise the rowwise matrix is
  // built by oppositeFormat() directly into B, with no additional copy (the
  // return value is elided, or at worst moved).
  if (isRowwise()) return *this;
  return oppositeFormat();
}
HighsSparseMatrix HighsSparseMatrix::getColwise() const {
  if (isColwise()) return *this;
  return oppositeFormat();
}

void HighsSparseMatrix::priceByRowDenseResult(
    std::vector<double>& result, const HVector& column,
    const HighsInt from_index, const HighsInt debug_report) const {
  // Assumes that result is zeroed beforehand - in case continuing
  // priceByRow after switch from sparse
  assert(this->isRowwise());
  for (HighsInt ix = from_index; ix < column.count; ix++) {
    HighsInt iRow = column.index[ix];
    double multiplier = column.array[iRow];
    // Determine whether p_end_ or the next start_ should be used to end the
    // loop
    HighsInt to_iEl;
    if (this->format_ == MatrixFormat::kRowwisePartitioned) {
      to_iEl = this->p_end_[iRow];
    } else {
      to_iEl = this->start_[iRow + 1];
    }
    if (debug_report == kDebugReportAll || debug_report == iRow)
      debugReportRowPrice(iRow, multiplier, to_iEl, result);
    for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
      HighsInt iCol = this->index_[iEl];
      double value0 = result[iCol];
      double value1 = value0 + multiplier * this->value_[iEl];
      result[iCol] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
    }
  }
}

void HighsSparseMatrix::priceByRowDenseResult(
    std::vector<HighsCDouble>& result, const HVector& column,
    const HighsInt from_index, const HighsInt debug_report) const {
  // Assumes that result is zeroed beforehand - in case continuing
  // priceByRow after switch from sparse
  assert(this->isRowwise());
  for (HighsInt ix = from_index; ix < column.count; ix++) {
    HighsInt iRow = column.index[ix];
    double multiplier = column.array[iRow];
    // Determine whether p_end_ or the next start_ should be used to end the
    // loop
    HighsInt to_iEl;
    if (this->format_ == MatrixFormat::kRowwisePartitioned) {
      to_iEl = this->p_end_[iRow];
    } else {
      to_iEl = this->start_[iRow + 1];
    }
    for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
      HighsInt iCol = this->index_[iEl];
      HighsCDouble value0 = result[iCol];
      HighsCDouble value1 = value0 + multiplier * this->value_[iEl];
      result[iCol] = (fabs((double)value1) < kHighsTiny) ? kHighsZero : value1;
    }
  }
}

void HighsSparseMatrix::debugReportRowPrice(
    const HighsInt iRow, const double multiplier, const HighsInt to_iEl,
    const vector<double>& result) const {
  if (this->start_[iRow] >= to_iEl) return;
  printf("Row %d: value = %11.4g", (int)iRow, multiplier);
  HighsInt num_print = 0;
  for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
    HighsInt iCol = this->index_[iEl];
    double value0 = result[iCol];
    double value1 = value0 + multiplier * this->value_[iEl];
    double value2 = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
    if (num_print % 5 == 0) printf("\n");
    printf("[%4d %11.4g] ", (int)(iCol), value2);
    num_print++;
  }
  printf("\n");
}

void HighsSparseMatrix::debugReportRowPrice(
    const HighsInt iRow, const double multiplier, const HighsInt to_iEl,
    const vector<HighsCDouble>& result) const {
  if (this->start_[iRow] >= to_iEl) return;
  printf("Row %d: value = %11.4g", (int)iRow, multiplier);
  HighsInt num_print = 0;
  for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
    HighsInt iCol = this->index_[iEl];
    double value0 = (double)result[iCol];
    double value1 = value0 + multiplier * this->value_[iEl];
    double value2 = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
    if (num_print % 5 == 0) printf("\n");
    printf("[%4d %11.4g] ", (int)(iCol), value2);
    num_print++;
  }
  printf("\n");
}

void HighsSparseMatrix::debugReportRowPrice(const HighsInt iRow,
                                            const double multiplier,
                                            const HighsInt to_iEl,
                                            HighsSparseVectorSum& sum) const {
  if (this->start_[iRow] >= to_iEl) return;
  if (!multiplier) return;
  printf("Row %d: value = %11.4g", (int)iRow, multiplier);
  HighsInt num_print = 0;
  for (HighsInt iEl = this->start_[iRow]; iEl < to_iEl; iEl++) {
    HighsInt iCol = this->index_[iEl];
    sum.add(iCol, multiplier * this->value_[iEl]);
    if (num_print % 5 == 0) printf("\n");
    printf("[%4d %11.4g] ", (int)(iCol), sum.getValue(iCol));
    num_print++;
  }
  printf("\n");
}

bool HighsSparseMatrix::legalDimensions(
    const HighsLogOptions& log_options) const {
  return legalMatrixDimensions(
      log_options, "HighsSparseMatrix", this->name_, this->numVec(),
      this->format_ == MatrixFormat::kRowwisePartitioned, this->start_,
      this->p_end_, this->index_, this->value_);
}

HighsStatus HighsSparseMatrix::trimDimensions(
    const HighsLogOptions& log_options) {
  return trimMatrixDimensions(
      log_options, "HighsSparseMatrix", this->name_, this->numVec(),
      this->format_ == MatrixFormat::kRowwisePartitioned, this->start_,
      this->p_end_, this->index_, this->value_);
}
