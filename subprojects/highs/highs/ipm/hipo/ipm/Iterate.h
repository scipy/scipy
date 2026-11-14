#ifndef HIPO_ITERATE_H
#define HIPO_ITERATE_H

#include <vector>

#include "Info.h"
#include "IpmData.h"
#include "LinearSolver.h"
#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Holds the Newton direction Delta(x,y,xl,xu,zl,zu)
struct NewtonDir {
  std::vector<double> x{};
  std::vector<double> y{};
  std::vector<double> xl{};
  std::vector<double> xu{};
  std::vector<double> zl{};
  std::vector<double> zu{};

  NewtonDir(Int m, Int n);
  void clear();
  void add(const NewtonDir& d);
};

struct Residuals {
  std::vector<double> r1, r2, r3, r4, r5, r6;
};

struct Iterate {
  Model& model;
  Info& info;
  IpmData data;
  std::vector<double> x, xl, xu, y, zl, zu;
  Residuals res;
  NewtonDir delta;
  double pobj, dobj, pinf, dinf, pdgap;
  double mu;
  std::vector<double> scaling;
  double best_mu;
  const Regularisation& regul;
  std::vector<double> total_reg;
  double* Rp;
  double* Rd;
  Residuals ires;  // residuals for iterative refinement

  // statistics for stagnation
  Int bad_iter_{};
  double largest_dx_x_{}, largest_dy_y_{};
  double best_pinf_ = kHighsInf, best_dinf_ = kHighsInf;

  double best_violation = kHighsInf;
  Int best_iter = -1;
  std::vector<double> best_x, best_xl, best_xu, best_y, best_zl, best_zu;

  // ===================================================================================
  // Functions to construct, clear and check for nan or inf
  // ===================================================================================
  Iterate(Model& model_input, Info& info_input, Regularisation& r);

  // clear existing data
  void clearIter();
  void clearRes();
  void clearDir();
  void clearIres();

  // check if any component is nan or infinite
  bool isNan() const;
  bool isInf() const;
  bool isResNan() const;
  bool isResInf() const;
  bool isDirNan(const NewtonDir& d) const;
  bool isDirInf(const NewtonDir& d) const;

  // ===================================================================================
  // Compute:
  //  mu = \sum xl(i) * zl(i) + \sum xu(j) * zu(j)
  // for variables: i with a finite lower bound
  //                j with a finite upper bound
  // ===================================================================================
  void computeMu();

  // ===================================================================================
  // Compute diagonal scaling Theta^{-1}
  //  Theta^{-1}_{ii} = zl(i) / xl(i) + zu(i) / xu(i)
  // Theta^{-1} only considers the terms above if the corresponding upper/lower
  // bound is finite.
  // ===================================================================================
  void computeScaling();

  // ===================================================================================
  // Compute convergence indicators
  // - primal and dual objectives, and primal-dual relative gap
  // - primal and dual infeasibilities (either scaled or unscaled)
  // - complementairy products
  // ===================================================================================
  void indicators();

  // objectives
  void primalObj();
  void dualObj();
  void pdGap();

  // infeasibilities
  void primalInfeas();
  void dualInfeas();
  void primalInfeasUnscaled();
  void dualInfeasUnscaled();
  double infeasAfterDropping() const;

  // complementarity products
  void products();

  // ===================================================================================
  // Compute:
  //  res1 = rhs - A * x
  //  res2 = lower - x + xl
  //  res3 = upper - x - xu
  //  res4 = c - A^T * y - zl + zu + Q * x
  // Components of residuals 2,3 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void residual1234();

  // ===================================================================================
  // Compute:
  //  res5 = sigma * mu * e - Xl * Zl * e
  //  res6 = sigma * mu * e - Xu * Zu * e
  // Components of residuals 5,6 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void residual56(double sigma);

  // ===================================================================================
  // Compute:
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  // ===================================================================================
  std::vector<double> residual7(const Residuals& r) const;

  // ===================================================================================
  // Compute:
  //  res8 = res1 + A * Theta * res7
  // ===================================================================================
  std::vector<double> residual8(const Residuals& r,
                                const std::vector<double>& res7) const;

  // ===================================================================================
  // Construct a complementary point (x,y,z), such that for each j, either xj is
  // at one of the bounds (lower or upper), or zj is zero.
  // ===================================================================================
  void dropToComplementarity(std::vector<double>& x_cmp,
                             std::vector<double>& y_cmp,
                             std::vector<double>& z_cmp) const;

  void finalResiduals() const;

  // ===================================================================================
  // Compute residual of 6x6 linear system for iterative refinement.
  // ===================================================================================
  void residuals6x6(const NewtonDir& d);

  void makeStep(double alpha_primal, double alpha_dual);

  void saveBest(double feas_tol, double opt_tol, Int iter);
  bool resetBest(Int iter);

  bool stagnation(std::stringstream& log_stream);

  void getReg(LinearSolver& LS, const std::string& nla);

  void assertConsistency(Int n, Int m) const;
};

}  // namespace hipo

#endif