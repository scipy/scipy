#ifndef HIPO_SOLVER_H
#define HIPO_SOLVER_H

#include <string>

#include "Control.h"
#include "FactorHighsSolver.h"
#include "Info.h"
#include "Iterate.h"
#include "LinearSolver.h"
#include "Model.h"
#include "Options.h"
#include "Parameters.h"
#include "Status.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "ipm/hipo/factorhighs/FactorHighs.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsCallback.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "util/HighsSparseMatrix.h"
#include "util/HighsTimer.h"

namespace hipo {

using isFailure = bool;
using shouldTerminate = bool;
using isSuccess = bool;

class Solver {
  Model model_;
  std::unique_ptr<LinearSolver> LS_;
  std::unique_ptr<KktMatrix> kkt_;
  std::unique_ptr<Iterate> it_;
  Int m_{}, n_{};
  Int iter_{};
  double alpha_primal_{}, alpha_dual_{};
  double sigma_{};
  Regularisation regul_{};
  Info info_;
  Control control_;
  Options options_{};
  Options options_orig_{};
  HighsOptions Hoptions_{};
  ipx::LpSolver ipx_lps_;
  Logger logger_;
  double start_time_;

  // Status of ipm iterations up to pd feas solution found (potentially using
  // ipx to refine, in which case the ipx status_ipm is converted to hipo
  // status).
  Status status_phase1 = kStatusNotSet;

  // Status of ipm iterations after pd feas solution found, or status of
  // crossover if it is run, in which case ipx status_crossover is converted to
  // hipo status.
  Status status_phase2 = kStatusNotSet;

 public:
  // ===================================================================================
  // Load an LP or QP:
  //
  //  min   obj^T * x + 1/2 * x^T * Q * x
  //  s.t.  Ax {<=,=,>=} rhs
  //        lower <= x <= upper
  //
  // Transform constraints in equalities by adding slacks to inequalities:
  //  <= : add slack    0 <= s_i <= +inf
  //  >= : add slack -inf <= s_i <=    0
  // ===================================================================================
  Int load(const HighsLp& lp, const HighsHessian& Q);

  // ===================================================================================
  // Specify options, callback and timer.
  // ===================================================================================
  void setOptions(const HighsOptions& highs_options);
  void setCallback(HighsCallback& callback);
  void setTimer(const HighsTimer& timer);

  // ===================================================================================
  // Solve the LP or QP
  // ===================================================================================
  void solve();

  // ===================================================================================
  // Extract information
  // ===================================================================================
  void getInteriorSolution(std::vector<double>& x, std::vector<double>& xl,
                           std::vector<double>& xu, std::vector<double>& slack,
                           std::vector<double>& y, std::vector<double>& zl,
                           std::vector<double>& zu) const;
  Int getBasicSolution(std::vector<double>& x, std::vector<double>& slack,
                       std::vector<double>& y, std::vector<double>& z,
                       Int* cbasis, Int* vbasis) const;
  const Info& getInfo() const;
  void getOriginalDims(Int& num_row, Int& num_col) const;

  bool solved() const;
  bool stopped() const;
  bool failed() const;

  void setIpxTimerOffset(const double offset) {
    this->ipx_lps_.setTimerOffset(offset);
  }

 private:
  // Functions to run the various stages of the ipm
  void doSolve();
  void runIpm();
  isFailure initialise();
  shouldTerminate prepareIter();
  isFailure predictor();
  isFailure correctors();

  // ===================================================================================
  // Interface with IPX
  // ===================================================================================
  isFailure prepareIpx();
  isFailure prepareIpxStartingPoint();
  void runIpx();
  void refineWithIpx();
  void crossoverWithIpx();

  // ===================================================================================
  // Determine the maximum number of correctors to use, based on the relative
  // cost of factorisation and solve. Based on the heuristic in "Multiple
  // Centrality Corrections in a Primal-Dual Method for Linear Programming".
  // ===================================================================================
  void chooseNumberOfCorrectors();

  // ===================================================================================
  // Solve:
  //
  // ___Augmented system___
  //
  //      [ -Theta^{-1}-Q  A^T ] [ Deltax ] = [ res7 ]
  //      [ A              0   ] [ Deltay ] = [ res1 ]
  //
  // with:
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  //  Theta^{-1} = diag( scaling )
  //
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  //
  // OR
  //
  // ___Normal equations___
  //
  //      A * (Theta^{-1}+Q)^{-1} * A^T * Deltay = res8
  //      Delta x = (Theta^{-1}+Q)^{-1} * (A^T* Deltay - res7)
  //
  // with:
  //  res8 = res1 + A * (Theta^{-1}+Q)^{-1} * res7
  //
  // NB: normal equations available only if Q is zero or diagonal.
  // ===================================================================================
  isFailure solveNewtonSystem(NewtonDir& delta);
  void solve2x2(NewtonDir& delta, const Residuals& rhs);
  void solve6x6(NewtonDir& delta, const Residuals& rhs);

  // ===================================================================================
  // Reconstruct the solution of the full Newton system:
  //
  //  Deltaxl = Deltax - res2
  //  Deltaxu = res3 - Deltax
  //  Deltazl = Xl^{-1} * (res5 - zl * Deltaxl)
  //  Deltazu = Xu^{-1} * (res6 - zu * Deltaxu)
  // ===================================================================================
  void recoverDirection(NewtonDir& delta, const Residuals& rhs);

  // ===================================================================================
  // Functions for iterative refinement on the large 6x6 system
  // ===================================================================================
  void refine(NewtonDir& delta);
  double computeOmega(const NewtonDir& delta) const;

  // ===================================================================================
  // Steps to boundary are computed so that
  //
  //  x  + alpha_primal * Deltax
  //  xl + alpha_primal * Deltaxl > 0     (if lower bound finite)
  //  xu + alpha_primal * Deltaxu > 0     (if upper bound finite)
  //
  //  y  + alpha_dual * Deltay
  //  zl + alpha_dual * Deltazl > 0       (if lower bound finite)
  //  zu + alpha_dual * Deltazu > 0       (if upper bound finite)
  //
  // If cor is valid, the direction used is delta + weight * cor
  // If block is valid, the blocking index is returned.
  // ===================================================================================

  // step to boundary for single direction
  double stepToBoundary(const std::vector<double>& x,
                        const std::vector<double>& dx,
                        const std::vector<double>* cor, double weight, bool lo,
                        Int* block = nullptr) const;

  // primal and dual steps to boundary
  void stepsToBoundary(double& alpha_primal, double& alpha_dual,
                       const NewtonDir& delta, const NewtonDir* cor = nullptr,
                       double weight = 1.0) const;

  // ===================================================================================
  // Find stepsizes using Mehrotra heuristic.
  // Given the steps to the boundary for xl, xu, zl, zu, and the blocking
  // indices, find stepsizes so that the primal (resp dual) blocking variable
  // produces a complementarity product not too far from the mu that would be
  // obtained using the steps to the boundary.
  // ===================================================================================
  void stepSizes();

  // ===================================================================================
  // Make the step in the Newton direction with appropriate stepsizes.
  // ===================================================================================
  void makeStep();

  // ===================================================================================
  // Compute the Mehrotra starting point.
  // This requires to solve two linear systems with matrix A*A^T.
  // ===================================================================================
  isFailure startingPoint();

  // ===================================================================================
  // Compute the sigma to use for affine scaling direction or correctors, based
  // on the smallest stepsize of the previous iteration.
  // If stepsize was large, use small sigma to reduce mu.
  // If stepsize was small, use large sigma to re-centre.
  //
  //  alpha | sigma  |   sigma    |
  //        | affine | correctors |
  //  1.0   |--------|------------|
  //        |        |    0.01    |
  //  0.5   |        |------------|
  //        |        |    0.10    |
  //  0.2   |        |------------|
  //        |  0.01  |    0.25    |
  //  0.1   |        |------------|
  //        |        |    0.50    |
  //  0.05  |        |------------|
  //        |        |    0.90    |
  //  0.0   |--------|------------|
  //
  // ===================================================================================
  void sigmaAffine();
  void sigmaCorrectors();

  // ===================================================================================
  // Compute the residuals for the computation of multiple centrality
  // correctors.
  // ===================================================================================
  void residualsMcc();

  // ===================================================================================
  // Iteratively compute correctors, until they improve the stepsizes.
  // Based on Gondzio, "Multiple centrality corrections in a primal-dual method
  // for linear programming" and Colombo, Gondzio, "Further Development of
  // Multiple Centrality Correctors for Interior Point Methods".
  // ===================================================================================
  isFailure centralityCorrectors();

  // ===================================================================================
  // Given the current direction delta and the latest corrector, compute the
  // best primal and dual weights, that maximize the primal and dual stepsize.
  // ===================================================================================
  void bestWeight(const NewtonDir& delta, const NewtonDir& corrector,
                  double& wp, double& wd, double& alpha_p,
                  double& alpha_d) const;

  // ===================================================================================
  // If the current iterate is nan or inf, abort the iterations.
  // ===================================================================================
  shouldTerminate checkIterate();

  // ===================================================================================
  // Stop if detection is detected, or if the problem is primal or dual
  // infeasible.
  // ===================================================================================
  shouldTerminate checkStagnation();
  shouldTerminate checkBadIter();

  // ===================================================================================
  // Check the termination criterion:
  //  - primal infeasibility < tolerance
  //  - dual infeasiblity    < tolerance
  //  - relative dual gap    < tolerance
  // ===================================================================================
  shouldTerminate checkTermination();
  isSuccess checkTerminationKkt();

  // ===================================================================================
  // Check for user interrupt or time limit
  // ===================================================================================
  shouldTerminate checkInterrupt();

  // ===================================================================================
  // Check and set status
  // ===================================================================================
  bool statusNeedsRefinement() const;
  bool statusAllowsCrossover() const;
  bool refinementIsOn() const;
  bool crossoverIsOn() const;
  bool errorOrInterrupt() const;

  void setStatus1(Status status);
  void setStatus2(Status status);
  void setStatus(Status status);
  Status getStatus1() const;
  Status getStatus2() const;
  Status getStatus() const;
  void finaliseStatus();

  // ===================================================================================
  // Print to screen
  // ===================================================================================
  void printInfo() const;
  void printHeader() const;
  void printOutput(bool reset = false) const;
  void printSummary() const;

  void resetOptions();
  void reset();
  void resetToBestIter(Int iter, bool print = true);

  isFailure initialiseLinearSolver();
  isSuccess switchToMultifrontal();

  void chooseAllowedParallelism(const HighsOptions& highs_options);
};

}  // namespace hipo

#endif
