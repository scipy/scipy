#ifndef HIPO_MODEL_H
#define HIPO_MODEL_H

#include <limits>
#include <string>
#include <vector>

#include "PreProcess.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsLp.h"
#include "model/HighsHessian.h"
#include "util/HighsSparseMatrix.h"

namespace hipo {

// Optimization problem:
// min  c^T * x + 1/2 x^T * Q * x
// s.t. A * x = b
//      x - xl = lower
//      x + xu = upper
//      xl, xu >= 0
//
// A is of size num_con x num_var, stored in CSC format using ptr, rows, vals.
//
// See Schork, Gondzio "Implementation of an interior point method with basis
// preconditioning", Math. Prog. Comput. 12, 2020
//

class Model {
 private:
  // data of original problem
  Int n_orig_{};
  Int m_orig_{};
  const HighsLp* lp_orig_ = nullptr;
  double offset_;
  ObjSense sense_ = ObjSense::kMinimize;
  std::vector<double> rhs_orig_;
  std::vector<char> constraints_orig_;

  // data of reformulated problem
  Int n_{};
  Int m_{};
  std::vector<double> c_{};
  std::vector<double> b_{};
  std::vector<double> lower_{};
  std::vector<double> upper_{};
  HighsSparseMatrix A_{};
  HighsHessian Q_{};
  std::vector<char> constraints_{};

  Int64 AS_nz_, NE_nz_lb_, NE_nz_ub_;

  std::vector<double> colscale_, rowscale_;

  double norm_unscaled_rhs_, norm_scaled_rhs_, norm_unscaled_obj_,
      norm_scaled_obj_;

  bool ready_ = false;

  Preprocessor preprocessor_;

  // norms of rows and cols of A
  std::vector<double> one_norm_cols_, one_norm_rows_, inf_norm_cols_,
      inf_norm_rows_;

  std::vector<Int> free_variables_;

  std::vector<double> best_lower_, best_upper_;

  void preprocess();
  Int checkData() const;
  void computeNorms();

 public:
  Int init(const HighsLp& lp, const HighsHessian& Q);

  void print(const Logger& logger) const;

  void printDense() const;

  void postprocess(std::vector<double>& x, std::vector<double>& xl,
                   std::vector<double>& xu, std::vector<double>& slack,
                   std::vector<double>& y, std::vector<double>& zl,
                   std::vector<double>& zu, const Iterate& it) const;

  double normScaledRhs() const { return norm_scaled_rhs_; }
  double normScaledObj() const { return norm_scaled_obj_; }
  double normUnscaledObj() const { return norm_unscaled_obj_; }
  double normUnscaledRhs() const { return norm_unscaled_rhs_; }

  void nzBounds();

  void adjustFreeVars(std::vector<double>& x, std::vector<double>& xl,
                      std::vector<double>& xu, const Logger& logger);

  void saveBounds();
  void resetBounds();

  bool hasLb(Int j) const { return std::isfinite(lower_[j]); }
  bool hasUb(Int j) const { return std::isfinite(upper_[j]); }

  Int m() const { return m_; }
  Int n() const { return n_; }
  Int n_orig() const { return n_orig_; }
  Int m_orig() const { return m_orig_; }
  const HighsLp* lpOrig() const { return lp_orig_; }
  const std::vector<double>& rhsOrig() const { return rhs_orig_; }
  const std::vector<char>& constraintsOrig() const { return constraints_orig_; }
  bool qp() const { return !Q_.empty(); }
  bool nonSeparableQp() const { return qp() && !Q_.isDiagonal(); }
  double sense() const { return (double)sense_; }
  const HighsSparseMatrix& A() const { return A_; }
  const HighsHessian& Q() const { return Q_; }
  const std::vector<double>& b() const { return b_; }
  const std::vector<double>& c() const { return c_; }
  double lb(Int i) const { return lower_[i]; }
  double ub(Int i) const { return upper_[i]; }
  char constraint(Int i) const { return constraints_[i]; }
  double colScale(Int i) const { return scaled() ? colscale_[i] : 1.0; }
  double rowScale(Int i) const { return scaled() ? rowscale_[i] : 1.0; }
  bool ready() const { return ready_; }
  bool scaled() const { return !colscale_.empty(); }
  double offset() const { return offset_; }
  double oneNormRows(Int i) const { return one_norm_rows_[i]; }
  double oneNormCols(Int i) const { return one_norm_cols_[i]; }
  double infNormRows(Int i) const { return inf_norm_rows_[i]; }
  double infNormCols(Int i) const { return inf_norm_cols_[i]; }
  Int64 nzAS() const { return AS_nz_; }
  Int64 nzNElb() const { return NE_nz_lb_; }
  Int64 nzNEub() const { return NE_nz_ub_; }
  const std::vector<Int>& freeVars() const { return free_variables_; }

  Int loadIntoIpx(ipx::LpSolver& lps) const;

  // classes for preprocessing
  friend struct PreprocessEmptyRows;
  friend struct PreprocessFixedVars;
  friend struct PreprocessScaling;
  friend struct PreprocessFormulation;
  friend struct PreprocessFreeVars;
};

}  // namespace hipo

#endif
