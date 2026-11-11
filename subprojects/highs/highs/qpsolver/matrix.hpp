/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_MATRIX_HPP__
#define __SRC_LIB_MATRIX_HPP__

#include <cassert>
#include <vector>

#include "lp_data/HStruct.h"
#include "qpvector.hpp"

struct MatrixBase {
  HighsInt num_row;
  HighsInt num_col;
  std::vector<HighsInt> start;
  std::vector<HighsInt> index;
  std::vector<double> value;
  HessianOracle oracle_;

  void clear() {
    num_row = 0;
    num_col = 0;
    start.clear();
    index.clear();
    value.clear();
    oracle_.clear();
  }

  bool isOracle() const { return oracle_.call_ != nullptr; }

  bool testOracle() const { return this->num_col > 0 && this->isOracle(); }

  bool doubleRelEqual(const double v_check, const double v_true) const {
    double rel_v = std::fabs(v_check - v_true) / (1.0 + std::fabs(v_true));
    return rel_v < 1e-5;
  }

  bool doubleVectorRelEqual(const double* v_check, const double* v_true) const {
    for (HighsInt iRow = 0; iRow < this->num_col; iRow++)
      if (!doubleRelEqual(v_check[iRow], v_true[iRow])) return false;
    return true;
  }

  bool intVectorEqual(const HighsInt* v_check, const HighsInt* v_true) const {
    for (HighsInt iRow = 0; iRow < this->num_col; iRow++)
      if (v_check[iRow] != v_true[iRow]) return false;
    return true;
  }

  bool packedVectorRelEqual(const HighsInt v_check_num_entries,
                            const HighsInt* v_check_index,
                            const double* v_check_value,
                            const HighsInt v_true_num_entries,
                            const HighsInt* v_true_index,
                            const double* v_true_value) const {
    if (v_check_num_entries != v_true_num_entries) return false;
    std::vector<double> v_check(this->num_col, 0);
    std::vector<double> v_true(this->num_col, 0);
    for (HighsInt iX = 0; iX < v_true_num_entries; iX++) {
      v_check[v_check_index[iX]] = v_check_value[iX];
      v_true[v_true_index[iX]] = v_true_value[iX];
    }
    return scatteredVectorRelEqual(v_check_num_entries, v_check_index,
                                   v_check.data(), v_true_num_entries,
                                   v_true_index, v_true.data());
  }

  bool scatteredVectorRelEqual(const HighsInt v_check_num_entries,
                               const HighsInt* v_check_index,
                               const double* v_check_value,
                               const HighsInt v_true_num_entries,
                               const HighsInt* v_true_index,
                               const double* v_true_value) const {
    if (v_check_num_entries != v_true_num_entries) return false;
    std::vector<HighsInt> v_check(this->num_col, -1);
    std::vector<HighsInt> v_true(this->num_col, -1);
    for (HighsInt iX = 0; iX < v_true_num_entries; iX++) {
      HighsInt check_row = v_check_index[iX];
      if (v_check[check_row] != -1) return false;
      v_true[v_true_index[iX]] = v_true_index[iX];
      v_check[check_row] = v_check_index[iX];
    }
    return intVectorEqual(v_check.data(), v_true.data()) &&
           doubleVectorRelEqual(v_check_value, v_true_value);
  }

  void callLog(const std::string& message, const HighsInt id = 0) const {
    return;
    printf("%s: %d\n", message.c_str(), int(id));
  }

  double diagonal(const HighsInt iCol) const {
    double diagonal_value = 0;
    if (this->num_col > 0) {
      for (HighsInt iEl = this->start[iCol]; iEl < this->start[iCol + 1]; iEl++)
        if (this->index[iEl] == iCol) {
          diagonal_value = this->value[iEl];
          break;
        }
    }
    if (this->isOracle()) {
      callLog("MatrixBase::diagonal", iCol);
      double oracle_diagonal_value = this->oracle_.diagonal(iCol);
      if (testOracle())
        assert(doubleRelEqual(oracle_diagonal_value, diagonal_value));
      return oracle_diagonal_value;
    }
    return diagonal_value;
  }

  void getColumn(const HighsInt iCol, HighsInt& num_entries, HighsInt* index,
                 double* value) {
    if (this->num_col > 0) {
      num_entries = 0;
      for (HighsInt iEl = this->start[iCol]; iEl < this->start[iCol + 1];
           iEl++) {
        index[num_entries] = this->index[iEl];
        value[num_entries] = this->value[iEl];
        num_entries++;
      }
    }
    if (this->isOracle()) {
      callLog("MatrixBase::getColumn", iCol);
      if (testOracle()) {
        HighsInt oracle_num_entries;
        std::vector<HighsInt> oracle_index(this->oracle_.dim_);
        std::vector<double> oracle_value(this->oracle_.dim_);
        this->oracle_.getPackedColumn(iCol, oracle_num_entries,
                                      oracle_index.data(), oracle_value.data());
        assert(packedVectorRelEqual(oracle_num_entries, oracle_index.data(),
                                    oracle_value.data(), num_entries, index,
                                    value));
      } else {
        this->oracle_.getPackedColumn(iCol, num_entries, index, value);
      }
      return;
    }
  }

  QpVector& mat_vec(const QpVector& other, QpVector& target) const {
    return mat_vec_seq(other, target);
  }

  QpVector& mat_vec_seq(const QpVector& other, QpVector& target) const {
    target.reset();
    if (this->num_col > 0) {
      for (HighsInt i = 0; i < other.num_nz; i++) {
        HighsInt col = other.index[i];
        for (HighsInt idx = start[col]; idx < start[col + 1]; idx++) {
          HighsInt row = index[idx];
          target.value[row] += value[idx] * other.value[col];
        }
      }
    }
    if (isOracle()) {
      callLog("MatrixBase::mat_vec_seq");
      // product packs values in other, which yields nullptr
      // if other.num_nz = 0, and then assert in the oracle, so just
      // return target - which is correct since it's zero
      if (other.num_nz == 0) return target;
      // Need local copy of dim_ since mat_vec_seq is const
      HighsInt dim = oracle_.dim_;
      if (testOracle()) {
        std::vector<double> oracle_value(dim);
        oracle_.product(other.num_nz, other.index.data(), other.value.data(),
                        oracle_value.data());
        assert(doubleVectorRelEqual(oracle_value.data(), target.value.data()));
      } else {
        oracle_.product(other.num_nz, other.index.data(), other.value.data(),
                        target.value.data());
      }
    }
    target.resparsify();
    return target;
  }

  QpVector mat_vec(const QpVector& other) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->isOracle());
    HighsInt dim = this->isOracle() ? this->oracle_.dim_ : this->num_row;
    QpVector result(dim);
    mat_vec(other, result);
    return result;
  }

  QpVector vec_mat(HighsInt* idx, double* val, HighsInt nnz) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->isOracle());
    HighsInt dim = this->isOracle() ? this->oracle_.dim_ : this->num_col;
    QpVector result(dim);
    for (HighsInt i = 0; i < num_col; i++) {
      double dot = 0.0;
      // HighsInt idx_other = 0;
      // HighsInt idx_this = start[i];
      // while (idx_this < start[i+1] && idx_other < nnz) {
      //    if (idx[idx_other] == index[idx_this]) {
      //       dot += val[idx_other] * value[idx_this];
      //    } else if (idx[idx_other] < index[idx_this]) {
      //       idx_other++;
      //    } else {
      //       idx_this++;
      //    }
      // }

      for (HighsInt j = start[i]; j < start[i + 1]; j++) {
        // does the vector have an entry for index index[j]?
        double other_value = 0.0;
        for (HighsInt k = 0; k < nnz; k++) {
          if (idx[k] == index[j]) {
            other_value = val[k];
            break;
          }
        }

        dot += other_value * value[j];
      }

      if (dot != 0.0) {
        result.value[i] = dot;
        result.index[result.num_nz] = i;
        result.num_nz++;
      }
    }
    return result;
  }

  QpVector& vec_mat(const QpVector& other, QpVector& target) const {
    return vec_mat_1(other, target);
  }

  QpVector& vec_mat_1(const QpVector& other, QpVector& target) const {
    target.reset();
    if (other.num_nz == 0) return target;
    if (num_col > 0) {
      for (HighsInt col = 0; col < num_col; col++) {
        double dot = 0.0;
        for (HighsInt j = start[col]; j < start[col + 1]; j++) {
          dot += other.value[index[j]] * value[j];
        }
        target.value[col] = dot;
      }
    }
    if (this->isOracle()) {
      callLog("MatrixBase::vec_mat_1");
      HighsInt target_dim = other.dim;
      if (testOracle()) {
        HighsInt dim = oracle_.dim_;
        std::vector<double> oracle_value(dim);
        this->oracle_.product(other.num_nz, other.index.data(),
                              other.value.data(), oracle_value.data());
        assert(doubleVectorRelEqual(oracle_value.data(), target.value.data()));
      } else {
        this->oracle_.product(other.num_nz, other.index.data(),
                              other.value.data(), target.value.data());
      }
    }
    target.resparsify();
    return target;
  }

  QpVector vec_mat(const QpVector& other) const {
    // Check to see whether this is needed for a Hessian oracle
    HighsInt dim = this->isOracle() ? this->oracle_.dim_ : this->num_col;
    QpVector result(dim);

    return vec_mat(other, result);
  }

  // computes this * mat, where "this" is a transposed matrix
  MatrixBase tran_mat_(const MatrixBase& other) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->isOracle());
    MatrixBase res;
    res.num_row = num_col;
    res.num_col = other.num_col;

    res.start.push_back(0);
    QpVector buffer_col(other.num_row);
    QpVector buffer_col_res(num_col);
    for (HighsInt r = 0; r < other.num_col; r++) {
      other.extractcol(r, buffer_col);

      vec_mat(buffer_col, buffer_col_res);
      for (HighsInt i = 0; i < buffer_col_res.num_nz; i++) {
        res.index.push_back(buffer_col_res.index[i]);
        res.value.push_back(buffer_col_res.value[buffer_col_res.index[i]]);
      }
      res.start.push_back(res.start[r] + buffer_col_res.num_nz);
    }

    return res;
  }

  QpVector& extractcol(HighsInt col, QpVector& target) const {
    HighsInt row_dim = this->isOracle() ? this->oracle_.dim_ : this->num_row;
    HighsInt col_dim = this->isOracle() ? this->oracle_.dim_ : this->num_col;
    assert(target.dim == row_dim);
    target.reset();

    if (col >= col_dim) {
      assert(!this->isOracle());
      target.index[0] = col - num_col;
      target.value[col - num_col] = 1.0;
      target.num_nz = 1;
    } else {
      if (this->num_col > 0) {
        for (HighsInt i = 0; i < start[col + 1] - start[col]; i++) {
          target.index[i] = index[start[col] + i];
          target.value[target.index[i]] = value[start[col] + i];
        }
        target.num_nz = start[col + 1] - start[col];
      }
      if (this->isOracle()) {
        callLog("MatrixBase::extractcol");
        if (testOracle()) {
          HighsInt oracle_num_entries;
          std::vector<HighsInt> oracle_index(this->oracle_.dim_);
          std::vector<double> oracle_value(this->oracle_.dim_);
          this->oracle_.getScatteredColumn(col, oracle_num_entries,
                                           oracle_index.data(),
                                           oracle_value.data());
          assert(scatteredVectorRelEqual(
              oracle_num_entries, oracle_index.data(), oracle_value.data(),
              target.num_nz, target.index.data(), target.value.data()));
        }
        this->oracle_.getScatteredColumn(
            col, target.num_nz, target.index.data(), target.value.data());
      }
    }

    return target;
  }

  QpVector extractcol(HighsInt col) const {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->isOracle());
    HighsInt dim = this->isOracle() ? this->oracle_.dim_ : this->num_row;
    QpVector res(dim);

    return extractcol(col, res);
  }

  HighsInt numNz() const {
    if (this->isOracle()) {
      HighsInt nnz = 0;
      for (HighsInt iCol = 0; iCol < this->oracle_.dim_; iCol++) {
        for (HighsInt iRow = 0; iRow < iCol; iRow++)
          if (this->oracle_.entry(iRow, iCol)) nnz += 2;
        if (this->oracle_.entry(iCol, iCol)) nnz++;
      }
      return nnz;
    } else {
      assert(this->start.size() >= static_cast<size_t>(this->num_col));
      return this->start[this->num_col];
    }
  }

  bool isDiagonal() const {
    if (this->isOracle()) {
      HighsInt nnz = 0;
      for (HighsInt iCol = 0; iCol < this->oracle_.dim_; iCol++) {
        for (HighsInt iRow = 0; iRow < iCol; iRow++)
          if (this->oracle_.entry(iRow, iCol)) return false;
      }
      return true;
    }
    // Only relevant if matrix is square
    assert(this->num_col == this->num_row);
    return this->numNz() == this->num_col;
  }
};

struct Matrix {
 private:
  MatrixBase tran;
  bool has_transpose = false;

  void transpose() {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    if (!has_transpose) {
      std::vector<std::vector<HighsInt>> row_indices(mat.num_row);
      std::vector<std::vector<double>> row_values(mat.num_row);

      for (HighsInt col = 0; col < mat.num_col; col++) {
        for (HighsInt entry = mat.start[col]; entry < mat.start[col + 1];
             entry++) {
          HighsInt row = mat.index[entry];
          double val = mat.value[entry];
          row_indices[row].push_back(col);
          row_values[row].push_back(val);
        }
      }
      tran.start.clear();
      tran.index.clear();
      tran.value.clear();
      tran.start.reserve(mat.num_row + 1);
      tran.index.reserve(mat.index.size());
      tran.value.reserve(mat.value.size());

      tran.start.push_back(0);
      for (HighsInt row = 0; row < mat.num_row; row++) {
        tran.index.insert(tran.index.end(), row_indices[row].begin(),
                          row_indices[row].end());
        tran.value.insert(tran.value.end(), row_values[row].begin(),
                          row_values[row].end());

        tran.start.push_back(tran.start[row] + row_indices[row].size());
      }

      tran.num_col = mat.num_row;
      tran.num_row = mat.num_col;
    }
  }

 public:
  MatrixBase mat;

  Matrix(HighsInt nr, HighsInt nc) {
    mat.num_row = nr;
    mat.num_col = nc;
  };

  Matrix(const MatrixBase& m, bool needstran) {
    mat = m;
    // if (needstran) {
    //    transpose();
    // }
  }

  void append(const QpVector& vec) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    if (mat.num_col == 0 && mat.start.size() == 0) {
      mat.start.push_back(0);
    }
    for (HighsInt i = 0; i < vec.num_nz; i++) {
      mat.index.push_back(vec.index[i]);
      mat.value.push_back(vec.value[vec.index[i]]);
    }
    mat.start.push_back(mat.start[mat.num_col] + vec.num_nz);
    mat.num_col++;
    has_transpose = false;
  }

  void append(HighsInt* idx, double* val, HighsInt nnz) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    if (mat.num_col == 0 && mat.start.size() == 0) {
      mat.start.push_back(0);
    }
    for (HighsInt i = 0; i < nnz; i++) {
      mat.index.push_back(idx[i]);
      mat.value.push_back(val[i]);
    }
    mat.start.push_back(mat.start[mat.num_col] + nnz);
    mat.num_col++;
    has_transpose = false;
  }

  void append(HighsInt num_nz, HighsInt* index, double* value) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    if (mat.num_col == 0 && mat.start.size() == 0) {
      mat.start.push_back(0);
    }
    for (HighsInt i = 0; i < num_nz; i++) {
      mat.index.push_back(index[i]);
      mat.value.push_back(value[i]);
    }
    mat.start.push_back(mat.start[mat.num_col] + num_nz);
    mat.num_col++;
    has_transpose = false;
  }

  void dropcol(HighsInt col) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    assert(col < mat.num_col);
    has_transpose = false;

    mat.index.erase(mat.index.begin() + mat.start[col],
                    mat.index.begin() + mat.start[col + 1]);
    mat.value.erase(mat.value.begin() + mat.start[col],
                    mat.value.begin() + mat.start[col + 1]);

    HighsInt num_elements_in_col = mat.start[col + 1] - mat.start[col];
    for (; col < mat.num_col; col++) {
      mat.start[col] = mat.start[col + 1] - num_elements_in_col;
    }
    mat.start.pop_back();
    mat.num_col--;
  }

  MatrixBase& t() {
    if (!has_transpose) {
      transpose();
      has_transpose = true;
    }
    return tran;
  }

  Matrix mat_mat(Matrix& other) {
    Matrix res(mat.num_row, 0);

    QpVector buffer(other.mat.num_row);
    QpVector buffer2(mat.num_col);
    for (HighsInt col = 0; col < other.mat.num_col; col++) {
      res.append(vec_mat(other.mat.extractcol(col, buffer), buffer2));
    }

    return res;
  }

  Matrix tran_mat(Matrix& other) {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    Matrix res(mat.num_col, 0);

    QpVector buffer(other.mat.num_row);
    QpVector buffer2(mat.num_row);
    for (HighsInt col = 0; col < other.mat.num_col; col++) {
      res.append(mat_vec(other.mat.extractcol(col, buffer), buffer2));
    }
    return res;
  }

  QpVector& mat_vec(const QpVector& other, QpVector& target) {
    return mat.mat_vec(other, target);
  }

  QpVector mat_vec(const QpVector& other) { return mat.mat_vec(other); }

  QpVector vec_mat(const QpVector& other) const { return mat.vec_mat(other); }

  QpVector& vec_mat(const QpVector& other, QpVector& target) const {
    return mat.vec_mat(other, target);
  }

  QpVector vec_mat(HighsInt* index, double* value, HighsInt num_nz) {
    return mat.vec_mat(index, value, num_nz);
  }

  void report(const std::string& name = "") const {
    // Check to see whether this is needed for a Hessian oracle
    assert(!this->mat.isOracle());
    if (name != "") {
      printf("%s:", name.c_str());
    }
    printf("[%" HIGHSINT_FORMAT " x %" HIGHSINT_FORMAT "]\n", mat.num_row,
           mat.num_col);
    printf("start: ");
    for (HighsInt i : mat.start) {
      printf("%" HIGHSINT_FORMAT " ", i);
    }
    printf("\n");

    printf("index: ");
    for (HighsInt i : mat.index) {
      printf("%" HIGHSINT_FORMAT " ", i);
    }
    printf("\n");

    printf("value: ");
    for (double d : mat.value) {
      printf("%lf ", d);
    }
    printf("\n");
  }
};

#endif
