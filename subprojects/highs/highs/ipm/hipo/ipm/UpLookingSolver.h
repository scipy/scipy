#ifndef HIPO_UP_LOOKING_SOLVER_H
#define HIPO_UP_LOOKING_SOLVER_H

#include "Info.h"
#include "IpmData.h"
#include "KktMatrix.h"
#include "LinearSolver.h"
#include "Model.h"

// Up-looking factorisation
// Based on Tim Davis "Direct Methods for Sparse Linear Systems", cholmod, qdldl

namespace hipo {

class UpLookingSolver : public LinearSolver {
  KktMatrix& kkt_;

  // reference to matrix that is actually used (AS or NE)
  std::vector<Int>& ptr_;
  std::vector<Int>& rows_;
  std::vector<double>& val_;

  const Int n_;

  std::vector<Int64> ptrL_;
  std::vector<Int> rowsL_;
  std::vector<double> valL_;

  std::vector<double> regularisation_;
  std::vector<Int> signs_;

  const double reg_threshold_ = 1e-13;
  const double reg_value_ = 2e-7;

  std::vector<Int> parent_;
  double flops_;

  std::vector<double> as_buffer_;

  Info& info_;
  IpmData& data_;
  const Regularisation& regul_;
  const Model& model_;

  void etreeAndCounts(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                      std::vector<Int64>& colcount);
  void factor(const std::vector<Int>& ptr, const std::vector<Int>& rows,
              const std::vector<double>& val);
  void solve(std::vector<double>& x);

 public:
  UpLookingSolver(KktMatrix& kkt, Info& info, IpmData& data,
                  const Regularisation& regul, const Model& model);

  Int factorAS(const std::vector<double>& scaling) override;
  Int factorNE(const std::vector<double>& scaling) override;
  Int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  Int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  Int setup() override;
  void clear() override { valid_ = false; }
  double flops() const override { return flops_; }
  double spops() const override { return 0; }
  double nz() const override { return valL_.size(); }
  void getReg(std::vector<double>& reg) override;

  LinearSolverType type() const override { return kUpLookingType; }
};

}  // namespace hipo

#endif