#ifndef HIPO_FACTORHIGHS_SOLVER_H
#define HIPO_FACTORHIGHS_SOLVER_H

#include <algorithm>
#include <atomic>

#include "Info.h"
#include "IpmData.h"
#include "KktMatrix.h"
#include "LinearSolver.h"
#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/factorhighs/FactorHighs.h"

namespace hipo {

class FactorHighsSolver : public LinearSolver {
  FHsolver FH_;
  KktMatrix& kkt_;
  const Regularisation& regul_;
  Info& info_;
  IpmData& data_;
  const Logger& logger_;
  const Model& model_;
  Options& options_;

  std::vector<double> as_buffer_;

  std::string ordering_AS_ = "none";
  std::string ordering_NE_ = "none";

  Int chooseNla();
  Int setNla();
  void setParallelBeforeSymbolic();
  void setParallelAfterSymbolic();
  void printParallel() const;
  Int chooseOrdering(const std::vector<Int>& rows, const std::vector<Int>& ptr,
                     const std::vector<Int>& signs, Symbolic& S,
                     std::string& ordering, const std::string& nla);
  Int analyseAS(Symbolic& S);
  Int analyseNE(Symbolic& S);

 public:
  FactorHighsSolver(KktMatrix& kkt, Options& options, const Model& model,
                    const Regularisation& regul, Info& info, IpmData& record,
                    const Logger& logger);

  Int factorAS(const std::vector<double>& scaling) override;
  Int factorNE(const std::vector<double>& scaling) override;
  Int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  Int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  Int setup() override;
  void clear() override;
  double flops() const override;
  double spops() const override;
  double nz() const override;
  void getReg(std::vector<double>& reg) override;
  LinearSolverType type() const override { return kFactorHighsType; }
};

}  // namespace hipo

#endif
