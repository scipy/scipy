#include "Solver.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "UpLookingSolver.h"
#include "ipm/IpxWrapper.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "lp_data/HighsSolution.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Int Solver::load(const HighsLp& lp, const HighsHessian& Q) {
  if (model_.init(lp, Q)) {
    logger_.printInfo("Error with model\n");
    return kErrorModel;
  }
  return kOk;
}

void Solver::chooseAllowedParallelism(const HighsOptions& Hoptions) {
  for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i)
    options_.parallel[i] = ParallelType::kChoose;

  // set default based on option `parallel`
  options_.setParallel(ParallelTechnique::kAnalyse, ParallelType::kOn);
  options_.setParallel(ParallelTechnique::kOrderNE, ParallelType::kOn);
  options_.setParallel(ParallelTechnique::kOrderAS, ParallelType::kOn);
  if (Hoptions.parallel == kHighsOffString) {
    options_.setParallel(ParallelTechnique::kTree, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kNode, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kForwardSolve, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kDiagonalSolve, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kBackwardSolve, ParallelType::kOff);
  }

  // override with option `hipo_parallel_type`
  if (Hoptions.parallel == kHighsOnString) {
    if (Hoptions.hipo_parallel_type == kHipoTreeString) {
      options_.setParallel(ParallelTechnique::kTree, ParallelType::kOn);
      options_.setParallel(ParallelTechnique::kNode, ParallelType::kOff);
    } else if (Hoptions.hipo_parallel_type == kHipoNodeString) {
      options_.setParallel(ParallelTechnique::kTree, ParallelType::kOff);
      options_.setParallel(ParallelTechnique::kNode, ParallelType::kOn);
    } else if (Hoptions.hipo_parallel_type == kHipoBothString) {
      options_.setParallel(ParallelTechnique::kTree, ParallelType::kOn);
      options_.setParallel(ParallelTechnique::kNode, ParallelType::kOn);
    }
  }

  // override if threads is 1
  if (highs::parallel::num_threads() == 1) {
    for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i)
      options_.parallel[i] = ParallelType::kOff;
  }

  // override with option `hipo_parallel_force` or `hipo_parallel_forbid`
  for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i) {
    bool force = testParallelBit(Hoptions.hipo_parallel_force, i);
    bool forbid = testParallelBit(Hoptions.hipo_parallel_forbid, i);
    if (force && forbid) continue;
    if (force) options_.parallel[i] = ParallelType::kOn;
    if (forbid) options_.parallel[i] = ParallelType::kOff;
  }
}

void Solver::setOptions(const HighsOptions& highs_options) {
  options_.display = highs_options.output_flag;
  options_.log_options = &highs_options.log_options;

  // Debug option is already considered through log_options.log_dev_level in
  // hipo::Logger::debug

  options_.timeless_log = highs_options.timeless_log;
  options_.feasibility_tol =
      std::min(highs_options.primal_feasibility_tolerance,
               highs_options.dual_feasibility_tolerance);
  options_.optimality_tol = highs_options.ipm_optimality_tolerance;
  options_.crossover_tol = highs_options.start_crossover_tolerance;

  if (highs_options.kkt_tolerance != kDefaultKktTolerance) {
    options_.feasibility_tol = highs_options.kkt_tolerance;
    options_.optimality_tol = 1e-1 * highs_options.kkt_tolerance;
    options_.crossover_tol = 1e-1 * highs_options.kkt_tolerance;
  }

  // hipo uses same timer as highs, so it is fine to pass the same time limit
  options_.time_limit = highs_options.time_limit;

  options_.max_iter = highs_options.ipm_iteration_limit;
  options_.crossover = highs_options.run_crossover;
  options_.nla = highs_options.hipo_system;
  options_.ordering = highs_options.hipo_ordering;
  options_.factor = highs_options.hipo_factor;
  options_.block_size = highs_options.hipo_block_size;
  options_.random_seed = highs_options.random_seed + 42;
  chooseAllowedParallelism(highs_options);

  options_orig_ = options_;
  Hoptions_ = highs_options;
  resetOptions();
}

void Solver::resetOptions() {
  options_ = options_orig_;
  if (options_.display) logger_.setOptions(options_.log_options);
  control_.setOptions(options_);
}
void Solver::setCallback(HighsCallback& callback) {
  control_.setCallback(callback);
}
void Solver::setTimer(const HighsTimer& timer) { control_.setTimer(timer); }

void Solver::reset() {
  m_ = model_.m();
  n_ = model_.n();

  info_ = Info{};

  info_.ipx_used = false;
  info_.m_solver = m_;
  info_.n_solver = n_;
  info_.m_original = model_.m_orig();
  info_.n_original = model_.n_orig();

  resetOptions();

  iter_ = 0;
  alpha_primal_ = 0.0;
  alpha_dual_ = 0.0;
  sigma_ = 0.0;
  regul_ = Regularisation{};

  LS_.reset();
  it_.reset();
}

void Solver::solve() {
  doSolve();
  info_.ipm_iter = iter_;
  finaliseStatus();
  printSummary();
}

void Solver::doSolve() {
  if (!model_.ready()) {
    info_.error = kErrorModel;
    return;
  }

  reset();

  // iterate object needs to be initialised before potentially interrupting
  it_.reset(new Iterate(model_, info_, regul_));
  if (!it_) {
    info_.error = kErrorFailedAllocation;
    return;
  }

  if (checkInterrupt()) return;

  printInfo();

  runIpm();
  if (errorOrInterrupt()) return;

  runIpx();
  if (errorOrInterrupt()) return;

  it_->finalResiduals();
}

void Solver::runIpm() {
  if (initialise()) return;

  while (iter_ < options_.max_iter) {
    if (prepareIter()) break;
    if (predictor()) break;
    if (correctors()) break;
    makeStep();
    printOutput();
  }
}

isFailure Solver::initialise() {
  Clock clock;

  start_time_ = control_.elapsed();

  kkt_.reset(new KktMatrix(model_, regul_, info_, logger_, options_));
  if (!kkt_) {
    info_.error = kErrorFailedAllocation;
    return true;
  }

  if (initialiseLinearSolver()) return true;

  it_->data.append();
  it_->data.back().factorisation_used = LS_->type();

  if (checkInterrupt()) return true;

  chooseNumberOfCorrectors();

  if (startingPoint()) return true;

  it_->residual1234();
  it_->computeMu();
  it_->indicators();

  printOutput();

  if (checkInterrupt()) return true;

  info_.times[kInitialiseTime] += clock.stop();
  return false;
}

shouldTerminate Solver::prepareIter() {
  Clock clock;

  it_->saveBest(options_.feasibility_tol, options_.optimality_tol, iter_);
  if (checkIterate()) return true;
  if (checkBadIter()) return true;
  if (checkTermination()) return true;
  if (checkInterrupt()) return true;

  ++iter_;

  model_.adjustFreeVars(it_->x, it_->xl, it_->xu, logger_);

  it_->clearDir();

  LS_->clear();

  it_->data.append();
  it_->data.back().factorisation_used = LS_->type();

  it_->computeScaling();

  info_.times[kPrepareTime] += clock.stop();
  return false;
}

isFailure Solver::predictor() {
  Clock clock;

  if (checkInterrupt()) return true;

  sigmaAffine();
  it_->residual56(sigma_);

  if (solveNewtonSystem(it_->delta)) return true;

  info_.times[kPredictorTime] += clock.stop();
  return false;
}

isFailure Solver::correctors() {
  Clock clock;

  if (checkInterrupt()) return true;

  sigmaCorrectors();
  if (centralityCorrectors()) return true;

  info_.times[kCorrectorsTime] += clock.stop();
  return false;
}

isFailure Solver::prepareIpx() {
  assert(!model_.qp());

  ipx::Parameters ipx_param;
  ipx_param.display = options_.display;
  ipx_param.highs_logging = true;
  ipx_param.log_options = options_.log_options;
  ipx_param.dualize = 0;

  if (options_.crossover == kHighsOnString)
    ipx_param.run_crossover = 1;
  else if (options_.crossover == kHighsOffString)
    ipx_param.run_crossover = 0;
  else
    ipx_param.run_crossover = -1;

  ipx_param.ipm_feasibility_tol = options_.feasibility_tol;
  ipx_param.ipm_optimality_tol = options_.optimality_tol;
  ipx_param.start_crossover_tol = options_.crossover_tol;
  ipx_param.time_limit = options_.time_limit;
  ipx_param.timeless_log = options_.timeless_log;
  ipx_param.ipm_maxiter = options_.max_iter - iter_;
  ipx_lps_.SetParameters(ipx_param);

  ipx_lps_.SetCallback(control_.callback());

  Int load_status = model_.loadIntoIpx(ipx_lps_);

  if (load_status) {
    logger_.printInfo("Error loading model into IPX\n");
    return true;
  }
  return false;
}

isFailure Solver::prepareIpxStartingPoint() {
  std::vector<double> x, xl, xu, slack, y, zl, zu;
  getInteriorSolution(x, xl, xu, slack, y, zl, zu);

  Int start_point_status = ipx_lps_.LoadIPMStartingPoint(
      x.data(), xl.data(), xu.data(), slack.data(), y.data(), zl.data(),
      zu.data());

  if (start_point_status) {
    logger_.printInfo("Error loading starting point into IPX\n");
    return true;
  }

  return false;
}

void Solver::runIpx() {
  if (checkInterrupt()) return;

  if (statusNeedsRefinement() && refinementIsOn()) {
    logger_.print("\nRestarting with IPX\n");
    refineWithIpx();
  } else if (statusAllowsCrossover() && crossoverIsOn()) {
    logger_.print("\nRunning crossover with IPX\n");
    crossoverWithIpx();
  }
}

void Solver::refineWithIpx() {
  if (prepareIpx()) return;
  if (prepareIpxStartingPoint()) return;
  ipx_lps_.Solve();
  info_.ipx_info = ipx_lps_.GetInfo();

  // Hipo status_1 is no progress. If ipx fails, leave current hipo status,
  // because no progress is better than failed status.

  const bool ipx_status_usable =
      info_.ipx_info.status_ipm != IPX_STATUS_not_run &&
      info_.ipx_info.status_ipm != IPX_STATUS_failed;

  if (ipx_status_usable) {
    info_.ipx_used = true;
    setStatus1(IpxToHipoStatus(info_.ipx_info.status_ipm));

    const bool crossover_status_usable =
        info_.ipx_info.status_crossover != IPX_STATUS_not_run &&
        info_.ipx_info.status_crossover != IPX_STATUS_failed;

    if (crossover_status_usable)
      setStatus2(IpxToHipoStatus(info_.ipx_info.status_crossover));
  }
}

void Solver::crossoverWithIpx() {
  // at the moment this is almost identical to refineWithIpx, but in the
  // future it will use ipx_lps_.CrossoverFromStartingPoint
  if (prepareIpx()) return;
  if (prepareIpxStartingPoint()) return;
  ipx_lps_.Solve();
  info_.ipx_info = ipx_lps_.GetInfo();

  // Hipo status_1 is optimal. If ipx fails, leave current hipo status,
  // because optimal without crossover is better than failed status.

  const bool ipx_status_usable =
      info_.ipx_info.status_ipm != IPX_STATUS_not_run &&
      info_.ipx_info.status_ipm != IPX_STATUS_failed;

  if (ipx_status_usable) {
    info_.ipx_used = true;
    setStatus2(IpxToHipoStatus(info_.ipx_info.status_ipm));

    const bool crossover_status_usable =
        info_.ipx_info.status_crossover != IPX_STATUS_not_run &&
        info_.ipx_info.status_crossover != IPX_STATUS_failed;

    if (crossover_status_usable)
      setStatus2(IpxToHipoStatus(info_.ipx_info.status_crossover));
  }
}

isFailure Solver::solveNewtonSystem(NewtonDir& delta) {
  solve6x6(delta, it_->res);
  bool terminate = info_.error;

  if (!terminate) {
    refine(delta);
    if (it_->isDirNan(delta)) {
      logger_.printInfo("Direction is nan\n");
      terminate = true;
      info_.error = kErrorNan;
    } else if (it_->isDirInf(delta)) {
      logger_.printInfo("Direction is inf\n");
      terminate = true;
      info_.error = kErrorNan;
    }
  }

  if (terminate) {
    if (switchToMultifrontal()) {
      // try FactorHighs factorisation before giving up
      it_->data.back().factorisation_used = LS_->type();
      logger_.printInfo("Switching to FactorHighs because of bad direction\n");
      resetToBestIter(iter_ - 1);

      // Solve again. This does not create an infinite loop, as
      // switchToMultifrontal() is false in the second solve.
      delta.clear();
      it_->data.back().omega = 0.0;
      terminate = solveNewtonSystem(delta);
    }
  }

  return terminate;
}

void Solver::solve2x2(NewtonDir& delta, const Residuals& rhs) {
  Clock clock;

  std::vector<double>& theta_inv = it_->scaling;

  std::vector<double> res7 = it_->residual7(rhs);

  // NORMAL EQUATIONS
  if (options_.nla == kHipoNormalEqString) {
    std::vector<double> res8 = it_->residual8(rhs, res7);

    if (!LS_->valid_) {
      info_.error = LS_->factorNE(theta_inv);
      if (info_.error) {
        logger_.printInfo("Error while factorising normal equations\n");
        return;
      }
      it_->getReg(*LS_, options_.nla);
    }

    info_.error = LS_->solveNE(res8, delta.y);
    if (info_.error) {
      logger_.printInfo("Error while solving normal equations\n");
      return;
    }

    // Compute delta.x
    // Deltax = A^T * Deltay - res7;
    delta.x = res7;
    model_.A().alphaProductPlusY(-1.0, delta.y, delta.x, true);
    vectorScale(delta.x, -1.0);

    // Deltax = (Theta^-1+Rp+Q)^-1 * Deltax
    for (Int i = 0; i < n_; ++i) {
      double denom = theta_inv[i] + regul_.primal;
      if (model_.qp()) denom += model_.sense() * model_.Q().diag(i);
      delta.x[i] /= denom;
    }

  }

  // AUGMENTED SYSTEM
  else {
    if (!LS_->valid_) {
      info_.error = LS_->factorAS(theta_inv);
      if (info_.error) {
        logger_.printInfo("Error while factorising augmented system\n");
        return;
      }
      it_->getReg(*LS_, options_.nla);
    }

    info_.error = LS_->solveAS(res7, rhs.r1, delta.x, delta.y);
    if (info_.error) {
      logger_.printInfo("Error while solving augmented system\n");
      return;
    }
  }

  info_.times[kSolve2x2Time] += clock.stop();
}

void Solver::solve6x6(NewtonDir& delta, const Residuals& rhs) {
  solve2x2(delta, rhs);
  if (info_.error) return;
  recoverDirection(delta, rhs);
}

void Solver::recoverDirection(NewtonDir& delta, const Residuals& rhs) {
  Clock clock;

  const std::vector<double>& xl = it_->xl;
  const std::vector<double>& xu = it_->xu;
  const std::vector<double>& zl = it_->zl;
  const std::vector<double>& zu = it_->zu;
  const std::vector<double>& res2 = rhs.r2;
  const std::vector<double>& res3 = rhs.r3;
  const std::vector<double>& res4 = rhs.r4;
  const std::vector<double>& res5 = rhs.r5;
  const std::vector<double>& res6 = rhs.r6;

  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      delta.xl[i] = delta.x[i] - res2[i];
      delta.zl[i] = (res5[i] - zl[i] * delta.xl[i]) / xl[i];
    } else {
      delta.xl[i] = 0.0;
      delta.zl[i] = 0.0;
    }
  }
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      delta.xu[i] = res3[i] - delta.x[i];
      delta.zu[i] = (res6[i] - zu[i] * delta.xu[i]) / xu[i];
    } else {
      delta.xu[i] = 0.0;
      delta.zu[i] = 0.0;
    }
  }

  info_.times[kRecoverTime] += clock.stop();
}

double Solver::stepToBoundary(const std::vector<double>& x,
                              const std::vector<double>& dx,
                              const std::vector<double>* cor, double weight,
                              bool lo, Int* block) const {
  // Compute the largest alpha s.t. x + alpha * dx >= 0.
  // If cor is valid, consider x + alpha * (dx + w * cor) instead.
  // Use lo=1 for xl and zl, lo=0 for xu and zu.
  // Return the blocking index in block.

  const double damp = 1.0 - 1e3 * std::numeric_limits<double>::epsilon();

  double alpha = 1.0;
  Int bl = -1;

  for (Int i = 0; i < n_; ++i) {
    if ((lo && model_.hasLb(i)) || (!lo && model_.hasUb(i))) {
      double c = (cor ? (*cor)[i] * weight : 0.0);
      if (x[i] + alpha * (dx[i] + c) < 0.0) {
        alpha = -(x[i] * damp) / (dx[i] + c);
        bl = i;
      }
    }
  }
  if (block) *block = bl;
  return alpha;
}

void Solver::stepsToBoundary(double& alpha_primal, double& alpha_dual,
                             const NewtonDir& delta, const NewtonDir* cor,
                             double weight) const {
  // compute primal and dual steps to boundary, given direction, corrector and
  // weight.

  Int block;
  double axl = stepToBoundary(it_->xl, delta.xl, cor ? &(cor->xl) : nullptr,
                              weight, true);
  double axu = stepToBoundary(it_->xu, delta.xu, cor ? &(cor->xu) : nullptr,
                              weight, false);
  double azl = stepToBoundary(it_->zl, delta.zl, cor ? &(cor->zl) : nullptr,
                              weight, true);
  double azu = stepToBoundary(it_->zu, delta.zu, cor ? &(cor->zu) : nullptr,
                              weight, false);

  alpha_primal = std::min(axl, axu);
  alpha_primal = std::min(alpha_primal, 1.0);
  alpha_dual = std::min(azl, azu);
  alpha_dual = std::min(alpha_dual, 1.0);
}

void Solver::stepSizes() {
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& dxl = it_->delta.xl;
  std::vector<double>& dxu = it_->delta.xu;
  std::vector<double>& dzl = it_->delta.zl;
  std::vector<double>& dzu = it_->delta.zu;

  // parameters for Mehrotra heuristic
  const double gamma_f = 0.9;
  const double gamma_a = 1.0 / (1.0 - gamma_f);

  // compute stepsizes and blocking components
  Int block_xl, block_xu, block_zl, block_zu;
  double alpha_xl = stepToBoundary(xl, dxl, nullptr, 0, true, &block_xl);
  double alpha_xu = stepToBoundary(xu, dxu, nullptr, 0, false, &block_xu);
  double alpha_zl = stepToBoundary(zl, dzl, nullptr, 0, true, &block_zl);
  double alpha_zu = stepToBoundary(zu, dzu, nullptr, 0, false, &block_zu);

  double max_p = std::min(alpha_xl, alpha_xu);
  double max_d = std::min(alpha_zl, alpha_zu);

  // compute mu with current stepsizes
  double mu_full = 0.0;
  Int num_finite = 0;
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu_full += (xl[i] + max_p * dxl[i]) * (zl[i] + max_d * dzl[i]);
      ++num_finite;
    }
    if (model_.hasUb(i)) {
      mu_full += (xu[i] + max_p * dxu[i]) * (zu[i] + max_d * dzu[i]);
      ++num_finite;
    }
  }
  mu_full /= num_finite;
  mu_full /= gamma_a;

  // compute new stepsizes based on Mehrotra heuristic

  // primal
  double alpha_p = 1.0;
  Int block_p = -1;
  if (max_p < 1.0) {
    if (alpha_xl <= alpha_xu) {
      block_p = block_xl;
      double temp = mu_full / (zl[block_p] + max_d * dzl[block_p]);
      alpha_p = (temp - xl[block_p]) / dxl[block_p];
    } else {
      block_p = block_xu;
      double temp = mu_full / (zu[block_p] + max_d * dzu[block_p]);
      alpha_p = (temp - xu[block_p]) / dxu[block_p];
    }
    alpha_p = std::max(alpha_p, gamma_f * max_p);
    alpha_p = std::min(alpha_p, 1.0);
    assert(block_p >= 0);
  }

  // dual
  double alpha_d = 1.0;
  Int block_d = -1;
  if (max_d < 1.0) {
    if (alpha_zl <= alpha_zu) {
      block_d = block_zl;
      double temp = mu_full / (xl[block_d] + max_p * dxl[block_d]);
      alpha_d = (temp - zl[block_d]) / dzl[block_d];
    } else {
      block_d = block_zu;
      double temp = mu_full / (xu[block_d] + max_p * dxu[block_d]);
      alpha_d = (temp - zu[block_d]) / dzu[block_d];
    }
    alpha_d = std::max(alpha_d, gamma_f * max_d);
    alpha_d = std::min(alpha_d, 1.0);
    assert(block_d >= 0);
  }

  alpha_primal_ = std::min(alpha_p, 1.0 - 1e-4);
  alpha_dual_ = std::min(alpha_d, 1.0 - 1e-4);

  assert(alpha_primal_ > 0 && alpha_primal_ < 1 && alpha_dual_ > 0 &&
         alpha_dual_ < 1);
}

void Solver::makeStep() {
  Clock clock;
  stepSizes();
  it_->makeStep(alpha_primal_, alpha_dual_);
  it_->residual1234();
  it_->computeMu();
  it_->indicators();
  info_.times[kStepTime] += clock.stop();
}

isFailure Solver::startingPoint() {
  std::vector<double>& x = it_->x;
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& y = it_->y;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;

  // *********************************************************************
  // x starting point
  // *********************************************************************
  // compute feasible x
  for (Int i = 0; i < n_; ++i) {
    x[i] = 0.0;
    x[i] = std::max(x[i], model_.lb(i));
    x[i] = std::min(x[i], model_.ub(i));
  }

  const std::vector<double> empty_scaling;
  std::vector<double> temp_m(m_);

  // store b in y
  y = model_.b();
  if (norm2(y) < 1e-6) vectorAdd(y, 1e-3);

  if (options_.nla == kHipoNormalEqString) {
    // use y to store b-A*x
    model_.A().alphaProductPlusY(-1.0, x, y);

    // solve A*A^T * dx = b-A*x with factorisation and store the result in
    // temp_m

    // factorise A*A^T
    info_.error = LS_->factorNE(empty_scaling);
    if (info_.error) {
      logger_.printInfo("Error while factorising normal equations\n");
      return true;
    }

    info_.error = LS_->solveNE(y, temp_m);
    if (info_.error) {
      logger_.printInfo("Error while solving normal equations\n");
      return true;
    }

  } else if (options_.nla == kHipoAugmentedString) {
    // obtain solution of A*A^T * dx = b-A*x by solving
    // [ -I  A^T] [...] = [ -x]
    // [  A   0 ] [ dx] = [ b ]

    info_.error = LS_->factorAS(empty_scaling);
    if (info_.error) {
      logger_.printInfo("Error while factorising augmented system\n");
      return true;
    }

    std::vector<double> rhs_x(n_);
    for (Int i = 0; i < n_; ++i) rhs_x[i] = -x[i];
    std::vector<double> lhs_x(n_);
    info_.error = LS_->solveAS(rhs_x, y, lhs_x, temp_m);
    if (info_.error) {
      logger_.printInfo("Error while solving augmented system\n");
      return true;
    }
  }

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  xl.assign(n_, 0.0);
  model_.A().alphaProductPlusY(1.0, temp_m, xl, true);

  // x += dx;
  vectorAdd(x, 1.0, xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
  {
    double violation{};
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xl[i] = x[i] - model_.lb(i);
        violation = std::min(violation, xl[i]);
      } else {
        xl[i] = 0.0;
      }
      if (model_.hasUb(i)) {
        xu[i] = model_.ub(i) - x[i];
        violation = std::min(violation, xu[i]);
      } else {
        xu[i] = 0.0;
      }
    }

    // shift to be positive
    violation = 1.0 + std::max(0.0, -1.5 * violation);
    vectorAdd(xl, violation);
    vectorAdd(xu, violation);
  }
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************

  std::vector<double> cPlusQx = model_.c();
  if (norm2(cPlusQx) < 1e-6) vectorAdd(cPlusQx, 1e-3);
  if (model_.qp()) model_.Q().alphaProductPlusY(model_.sense(), x, cPlusQx);

  double max_norm_x = std::max(infNorm(x), std::max(infNorm(xl), infNorm(xu)));
  if (infNorm(cPlusQx) > max_norm_x * 1e3) cPlusQx = model_.c();

  if (options_.nla == kHipoNormalEqString) {
    // compute A*(c+Qx)
    std::fill(temp_m.begin(), temp_m.end(), 0.0);
    model_.A().alphaProductPlusY(1.0, cPlusQx, temp_m);

    info_.error = LS_->solveNE(temp_m, y);
    if (info_.error) {
      logger_.printInfo("Error while solvingF normal equations\n");
      return true;
    }

  } else if (options_.nla == kHipoAugmentedString) {
    // obtain solution of A*A^T * y = A*(c+Qx) by solving
    // [ -I  A^T] [...] = [ c+Qx ]
    // [  A   0 ] [ y ] = [   0  ]

    std::vector<double> rhs_y(m_, 0.0);
    std::vector<double> lhs_x(n_);

    info_.error = LS_->solveAS(cPlusQx, rhs_y, lhs_x, y);
    if (info_.error) {
      logger_.printInfo("Error while solving augmented system\n");
      return true;
    }
  }
  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c + Q * x - A^T * y and store in zl
  zl = cPlusQx;
  model_.A().alphaProductPlusY(-1.0, y, zl, true);

  // split result between zl and zu
  {
    double violation = 0.0;
    for (Int i = 0; i < n_; ++i) {
      double val = zl[i];
      zl[i] = 0.0;
      zu[i] = 0.0;

      if (model_.hasLb(i) && model_.hasUb(i)) {
        zl[i] = 0.5 * val;
        zu[i] = -0.5 * val;
      } else if (model_.hasLb(i)) {
        zl[i] = val;
      } else if (model_.hasUb(i)) {
        zu[i] = -val;
      }

      violation = std::min(violation, zl[i]);
      violation = std::min(violation, zu[i]);
    }

    // shift to be positive

    violation = 1.0 + std::max(0.0, -1.5 * violation);
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        zl[i] += violation;
      }
      if (model_.hasUb(i)) {
        zu[i] += violation;
      }
    }
  }
  // *********************************************************************

  // *********************************************************************
  // improve centrality
  // *********************************************************************
  {
    double xsum{1.0};
    double zsum{1.0};
    double mu{1.0};

    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xsum += xl[i];
        zsum += zl[i];
        mu += xl[i] * zl[i];
      }
      if (model_.hasUb(i)) {
        xsum += xu[i];
        zsum += zu[i];
        mu += xu[i] * zu[i];
      }
    }

    double dx = 0.5 * mu / zsum;
    double dz = 0.5 * mu / xsum;

    vectorAdd(xl, dx);
    vectorAdd(xu, dx);
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        zl[i] += dz;
      }
      if (model_.hasUb(i)) {
        zu[i] += dz;
      }
    }
  }
  // *********************************************************************

  return false;
}

void Solver::sigmaAffine() {
  sigma_ = kSigmaAffine;

  it_->data.back().sigma_aff = sigma_;
}

void Solver::sigmaCorrectors() {
  if ((alpha_primal_ > 0.5 && alpha_dual_ > 0.5) || iter_ == 1) {
    sigma_ = 0.01;
  } else if (alpha_primal_ > 0.2 && alpha_dual_ > 0.2) {
    sigma_ = 0.1;
  } else if (alpha_primal_ > 0.1 && alpha_dual_ > 0.1) {
    sigma_ = 0.25;
  } else if (alpha_primal_ > 0.05 && alpha_dual_ > 0.05) {
    sigma_ = 0.5;
  } else {
    sigma_ = 0.9;
  }

  it_->data.back().sigma = sigma_;
}

void Solver::residualsMcc() {
  Clock clock;

  // compute right-hand side for multiple centrality correctors
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& res5 = it_->res.r5;
  std::vector<double>& res6 = it_->res.r6;
  double& mu = it_->mu;

  it_->clearRes();

  // stepsizes of current direction
  double alpha_p, alpha_d;
  stepsToBoundary(alpha_p, alpha_d, it_->delta);

  // compute increased stepsizes
  alpha_p = std::max(1.0, alpha_p + kMccIncreaseAlpha);
  alpha_d = std::max(1.0, alpha_d + kMccIncreaseAlpha);

  // compute trial point
  std::vector<double> xlt = xl;
  std::vector<double> xut = xu;
  std::vector<double> zlt = zl;
  std::vector<double> zut = zu;
  vectorAdd(xlt, 1.0, it_->delta.xl, alpha_p);
  vectorAdd(xut, 1.0, it_->delta.xu, alpha_p);
  vectorAdd(zlt, 1.0, it_->delta.zl, alpha_d);
  vectorAdd(zut, 1.0, it_->delta.zu, alpha_d);

  // compute right-hand side for mcc
  for (Int i = 0; i < n_; ++i) {
    // res5
    if (model_.hasLb(i)) {
      double prod = xlt[i] * zlt[i];
      if (prod < sigma_ * mu * kGammaCorrector) {
        // prod is small, we add something positive to res5

        double temp = sigma_ * mu * kGammaCorrector - prod;
        res5[i] += temp;

      } else if (prod > sigma_ * mu / kGammaCorrector) {
        // prod is large, we may subtract something large from res5.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu / kGammaCorrector);
        res5[i] += temp;
      }
    } else {
      res5[i] = 0.0;
    }

    // res6
    if (model_.hasUb(i)) {
      double prod = xut[i] * zut[i];
      if (prod < sigma_ * mu * kGammaCorrector) {
        // prod is small, we add something positive to res6

        double temp = sigma_ * mu * kGammaCorrector - prod;
        res6[i] += temp;

      } else if (prod > sigma_ * mu / kGammaCorrector) {
        // prod is large, we may subtract something large from res6.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu / kGammaCorrector);
        res6[i] += temp;
      }
    } else {
      res6[i] = 0.0;
    }
  }

  info_.times[kResidualsTime] += clock.stop();
}

isFailure Solver::centralityCorrectors() {
  // compute stepsizes of current direction
  double alpha_p_old, alpha_d_old;
  stepsToBoundary(alpha_p_old, alpha_d_old, it_->delta);

  Int cor;
  for (cor = 0; cor < info_.correctors; ++cor) {
    residualsMcc();

    NewtonDir corrector(m_, n_);
    if (solveNewtonSystem(corrector)) return true;

    double alpha_p, alpha_d;
    double wp = alpha_p_old * alpha_d_old;
    double wd = wp;
    bestWeight(it_->delta, corrector, wp, wd, alpha_p, alpha_d);

    if (alpha_p < alpha_p_old + kMccIncreaseAlpha * kMccIncreaseMin &&
        alpha_d < alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // reject corrector
      break;
    }

    if (alpha_p >= alpha_p_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // accept primal corrector
      vectorAdd(it_->delta.x, 1.0, corrector.x, wp);
      vectorAdd(it_->delta.xl, 1.0, corrector.xl, wp);
      vectorAdd(it_->delta.xu, 1.0, corrector.xu, wp);
      alpha_p_old = alpha_p;
    }
    if (alpha_d >= alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // accept dual corrector
      vectorAdd(it_->delta.y, 1.0, corrector.y, wd);
      vectorAdd(it_->delta.zl, 1.0, corrector.zl, wd);
      vectorAdd(it_->delta.zu, 1.0, corrector.zu, wd);
      alpha_d_old = alpha_d;
    }

    if (alpha_p > 0.95 && alpha_d > 0.95) {
      // stepsizes are large enough, stop
      ++cor;
      break;
    }

    // else, keep computing correctors
  }

  it_->data.back().correctors = cor;

  return false;
}

void Solver::bestWeight(const NewtonDir& delta, const NewtonDir& corrector,
                        double& wp, double& wd, double& alpha_p,
                        double& alpha_d) const {
  // Find the best primal and dual weights for the corrector in the interval
  // [alpha_p_old * alpha_d_old, 1].
  // Upon return, wp and wd are the optimal weights, alpha_p and alpha_d are
  // the corresponding stepsizes.

  // keep track of best stepsizes
  alpha_p = 0.0;
  alpha_d = 0.0;

  // initial weight
  double w = wp;

  // divide interval into 9 points
  const double step = (1.0 - w) / 8;

  // for each weight, compute stepsizes and save best ones
  for (; w <= 1.0; w += step) {
    double ap, ad;
    stepsToBoundary(ap, ad, delta, &corrector, w);
    if (ap > alpha_p) {
      alpha_p = ap;
      wp = w;
    }
    if (ad > alpha_d) {
      alpha_d = ad;
      wd = w;
    }

    if (step == 0.0) break;
  }
}

shouldTerminate Solver::checkIterate() {
  bool terminate = false;

  if (it_->isNan()) {
    logger_.printInfo("\nIterate is nan\n");
    info_.error = kErrorNan;
    terminate = true;
  } else if (it_->isInf()) {
    logger_.printInfo("\nIterate is inf\n");
    info_.error = kErrorNan;
    bool terminate = true;
  } else {
    for (Int i = 0; i < n_; ++i) {
      if ((model_.hasLb(i) && it_->xl[i] < 0) ||
          (model_.hasLb(i) && it_->zl[i] < 0) ||
          (model_.hasUb(i) && it_->xu[i] < 0) ||
          (model_.hasUb(i) && it_->zu[i] < 0)) {
        logger_.printInfo("\nIterate has negative component\n");
        info_.error = kErrorNegativeComponent;
        terminate = true;
      }
    }
  }

  return terminate;
}

shouldTerminate Solver::checkStagnation() {
  std::stringstream log_stream;
  bool stagnation = it_->stagnation(log_stream);
  logger_.printInfo(log_stream.str().c_str());
  return stagnation;
}

shouldTerminate Solver::checkBadIter() {
  bool terminate = false;
  bool stagnation = iter_ > 0 ? checkStagnation() : false;

  // check for infeasibility
  bool mu_is_large =
      it_->best_mu > 0.0 ? it_->mu > it_->best_mu * kDivergeTol : false;
  bool pobj_is_very_large =
      it_->pobj <
      -std::max(std::abs(it_->dobj) * kDivergeTol * kDivergeTol, 1.0);
  bool dobj_is_very_large =
      it_->dobj >
      std::max(std::abs(it_->pobj) * kDivergeTol * kDivergeTol, 1.0);
  bool clearly_infeasible =
      iter_ > 5 && (pobj_is_very_large || dobj_is_very_large);

  if (stagnation || mu_is_large || clearly_infeasible) {
    bool pobj_is_larger =
        it_->pobj < -std::max(std::abs(it_->dobj) * kDivergeTol, 1.0);
    bool dobj_is_larger =
        it_->dobj > std::max(std::abs(it_->pobj) * kDivergeTol, 1.0);

    if (pobj_is_larger) {
      // problem is likely to be primal unbounded, i.e. dual infeasible
      logger_.print("=== Dual infeasible\n");
      setStatus(kStatusDualInfeasible);
      terminate = true;
    } else if (dobj_is_larger) {
      // problem is likely to be dual unbounded, i.e. primal infeasible
      logger_.print("=== Primal infeasible\n");
      setStatus(kStatusPrimalInfeasible);
      terminate = true;
    } else if (stagnation) {
      if (switchToMultifrontal()) {
        logger_.printInfo("Switching to FactorHighs because of no progress\n");
        resetToBestIter(iter_);
      } else {
        if (getStatus1() == kStatusOptimal) {
          assert(crossoverIsOn());
          setStatus2(kStatusNoProgress);
          resetToBestIter(iter_);
        } else {
          if (checkTerminationKkt()) {
            logger_.printw(
                "HiPO stagnated but HiGHS considers the solution "
                "acceptable\n");
            logger_.print("=== Primal-dual feasible point found\n");
            setStatus1(kStatusOptimal);
          } else {
            setStatus1(kStatusNoProgress);
            resetToBestIter(iter_);
          }
        }

        terminate = true;
      }
    }
  }

  return terminate;
}

shouldTerminate Solver::checkTermination() {
  bool feasible = it_->pinf < options_.feasibility_tol &&
                  it_->dinf < options_.feasibility_tol;
  bool optimal = it_->pdgap < options_.optimality_tol;

  bool terminate = false;

  if (feasible && optimal) {
    if (crossoverIsOn()) {
      if (getStatus1() != kStatusOptimal) {
        logger_.print("=== Primal-dual feasible point found\n");
        setStatus1(kStatusOptimal);
      }
      bool ready_for_crossover =
          it_->infeasAfterDropping() < options_.crossover_tol;
      if (ready_for_crossover) {
        logger_.print("=== Ready for crossover\n");
        terminate = true;
      }

    } else {
      terminate = model_.qp() ? true : checkTerminationKkt();
      if (terminate) {
        logger_.print("=== Primal-dual feasible point found\n");
        setStatus1(kStatusOptimal);
      }
    }
  }

  return terminate;
}

isSuccess Solver::checkTerminationKkt() {
  if (model_.qp()) {
    // Not yet implemented for QP
    logger_.printInfo("kktCheck skipped for QP\n");
    return false;

  } else {
    // LP check
    logger_.printInfo("Solution may be optimal, perform kktCheck\n");
    HighsModelStatus model_status = HighsModelStatus::kOptimal;
    HighsInfo highs_info;
    HighsSolution highs_solution;

    // Allow kktCheck to print only in debug mode (this is a copy of highs
    // options, not the original)
    Hoptions_.output_flag = logger_.debug(2);

    if (!model_.lpOrig()) return false;

    getHipoNonVertexSolution(Hoptions_, *model_.lpOrig(), model_.n_orig(),
                             model_.m_orig(), model_.rhsOrig(),
                             model_.constraintsOrig(), *this, model_status,
                             highs_solution);

    lpNoBasisKktCheck(model_status, highs_info, *model_.lpOrig(),
                      highs_solution, Hoptions_, "During HiPO solve");

    if (model_status == HighsModelStatus::kOptimal) {
      logger_.printInfo("Check successful\n");
      return true;
    } else
      logger_.printInfo("Check failed\n");
  }

  return false;
}

shouldTerminate Solver::checkInterrupt() {
  bool terminate = false;
  Int status = control_.interruptCheck(iter_);
  if (status) {
    setStatus((Status)status);
    terminate = true;
  }
  return terminate;
}

void Solver::resetToBestIter(Int iter, bool print) {
  bool did_reset = it_->resetBest(iter);
  if (did_reset && print) printOutput(true);
}

isFailure Solver::initialiseLinearSolver() {
  LS_.reset(new FactorHighsSolver(*kkt_, options_, model_, regul_, info_,
                                  it_->data, logger_));
  if (!LS_) {
    info_.error = kErrorFailedAllocation;
    return true;
  }
  info_.error = LS_->setup();
  if (info_.error) return true;

  LS_->clear();

  // Use uplooking factorisation instead of multifrontal if
  // - very few operations needed or
  // - very few nz per column in L or
  // - multifrontal is dominated by assembly operations
  const bool should_switch_to_uplooking =
      LS_->flops() < kUplookFlopsThresh ||
      LS_->nz() < kUplookNzPerColLower * kkt_->n() ||
      LS_->flops() < kUplookSpopsRatioLower * LS_->spops() ||
      (LS_->nz() < kUplookNzPerColUpper * kkt_->n() &&
       LS_->flops() < kUplookSpopsRatioUpper * LS_->spops());

  bool do_switch;
  if (options_.factor == kHipoFactorMultifrontal)
    do_switch = false;
  else if (options_.factor == kHipoFactorUplooking)
    do_switch = true;
  else
    do_switch = should_switch_to_uplooking;

  if (do_switch) {
    LS_.reset(new UpLookingSolver(*kkt_, info_, it_->data, regul_, model_));
    if (!LS_) {
      info_.error = kErrorFailedAllocation;
      return true;
    }
    LS_->setup();
    LS_->clear();
  }

  return false;
}

isSuccess Solver::switchToMultifrontal() {
  bool switch_successful = false;

  if (LS_->type() == kUpLookingType && options_.factor == kHighsChooseString) {
    LS_.reset(new FactorHighsSolver(*kkt_, options_, model_, regul_, info_,
                                    it_->data, logger_));
    if (!LS_) {
      info_.error = kErrorFailedAllocation;
    } else {
      LS_->clear();
      it_->bad_iter_ = 0;
      it_->best_pinf_ = kHighsInf;
      it_->best_dinf_ = kHighsInf;
      it_->largest_dx_x_ = 0;
      it_->largest_dy_y_ = 0;
      info_.error = kOk;
      switch_successful = true;
    }
  }

  return switch_successful;
}

void Solver::printHeader() const {
  if (iter_ % 20 == 0) {
    logger_.print(
        " iter       primal obj         dual obj       pinf       dinf       "
        "gap");
    if (!options_.timeless_log) logger_.print("    time");
    if (logger_.debug(1)) {
      logger_.print(
          "     alpha p/d   sigma af/co   cor  solv  fact    static reg p/d  "
          "  "
          " minT     maxT  (xj * zj / mu)_range_&_num   max_res");
    }
    logger_.print("\n");
  }
}

void Solver::printOutput(bool reset) const {
  printHeader();

  if (!reset) {
    logger_.print("%5d ", iter_);
  } else {
    logger_.print(" Resetting to iteration %d\n", it_->best_iter);
    logger_.print("%5d ", it_->best_iter);
  }

  logger_.print("%16.8e %16.8e %10.2e %10.2e %9.2e", it_->pobj, it_->dobj,
                it_->pinf, it_->dinf, it_->pdgap);
  if (!options_.timeless_log) logger_.print(" %7.1f", control_.elapsed());
  if (logger_.debug(1) && !reset) {
    const IpmIterData& data = it_->data.back();
    logger_.print(
        " %6.2f %6.2f %6.2f %6.2f %5d %5d %5c %8.1e %8.1e %8.1e %8.1e %8.1e "
        "%8.1e %4d %4d %9.1e",
        alpha_primal_, alpha_dual_, data.sigma_aff, data.sigma, data.correctors,
        data.num_solves,
        LinearSolverTypeToString((LinearSolverType)data.factorisation_used)
            .front(),
        regul_.primal, regul_.dual, data.min_theta, data.max_theta,
        data.min_prod, data.max_prod, data.num_small_prod, data.num_large_prod,
        data.omega);
  }
  logger_.print("\n");
}

void Solver::printInfo() const {
  std::stringstream log_stream;
  log_stream << "\nRunning HiPO\n";
  log_stream << textline("Threads:") << highs::parallel::num_threads() << '\n';
  logger_.print(log_stream.str().c_str());

  model_.print(logger_);
}

void Solver::printSummary() const {
  std::stringstream log_stream;

  log_stream << "\nSummary\n";
  if (!options_.timeless_log)
    log_stream << textline("HiPO runtime:")
               << fix(control_.elapsed() - start_time_, 0, 2) << "\n";

  log_stream << textline("Status:") << statusString(info_.status) << "\n";
  log_stream << textline("Status Crossover:")
             << statusString(IpxToHipoStatus(info_.ipx_info.status_crossover))
             << "\n";
  if (info_.error)
    log_stream << textline("Error:") << errorString((Error)info_.error) << "\n";

  log_stream << textline("HiPO iterations:") << integer(iter_) << "\n";
  if (info_.ipx_used)
    log_stream << textline("IPX iterations:") << integer(info_.ipx_info.iter)
               << "\n";

  if (logger_.debug(1)) {
    log_stream << textline("Correctors:") << integer(info_.correctors) << '\n';
    log_stream << textline("Factorisations:") << integer(info_.factor_number)
               << '\n';
    log_stream << textline("Solves:") << integer(info_.solve_number) << '\n';

    if (!options_.timeless_log) {
      log_stream << textline("Avg time per factorisation:")
                 << sci(info_.times[kFactoriseTime] / info_.factor_number, 0, 2)
                 << '\n';
      log_stream << textline("Avg time per solve:")
                 << sci(info_.times[kSolveTime] / info_.solve_number, 0, 2)
                 << '\n';
    }
  }
  logger_.print(log_stream.str().c_str());

  if (logger_.debug(1)) {
    logger_.print("\nTime profile\n");

    logger_.print(" Initialise           %.2f\n", info_.times[kInitialiseTime]);
    logger_.print(" Prepare iter         %.2f\n", info_.times[kPrepareTime]);
    logger_.print(" Predictor            %.2f\n", info_.times[kPredictorTime]);
    logger_.print(" Correctors           %.2f\n", info_.times[kCorrectorsTime]);
    logger_.print(" Step                 %.2f\n", info_.times[kStepTime]);
    logger_.print("\n");

    logger_.print(" 2x2 system           %.2f\n", info_.times[kSolve2x2Time]);
    logger_.print(" Recover direction    %.2f\n", info_.times[kRecoverTime]);
    logger_.print(" Residuals            %.2f\n", info_.times[kResidualsTime]);
    logger_.print("\n");

    logger_.print(" Structure NE         %.2f\n",
                  info_.times[kMatrixStructureTime_NE]);
    logger_.print(" Structure AS         %.2f\n",
                  info_.times[kMatrixStructureTime_AS]);
    logger_.print(" Matrix values        %.2f\n",
                  info_.times[kMatrixValuesTime]);

    logger_.print(" Analyse NE           %.2f\n", info_.times[kAnalyseTime_NE]);
    logger_.print(" Analyse AS           %.2f\n", info_.times[kAnalyseTime_AS]);

    logger_.print(" Factorise            %.2f\n", info_.times[kFactoriseTime]);

    logger_.print(" Solve                %.2f\n", info_.times[kSolveTime]);
    logger_.print(" Residual             %.2f\n",
                  info_.times[kRefinementResTime]);
    logger_.print(" Omega                %.2f\n",
                  info_.times[kRefinementOmegaTime]);
    logger_.print(" Insert/split         %.2f\n",
                  info_.times[kInsertSplitTime]);
    logger_.print("\n");
  }

  logger_.print("\n");
}

void Solver::getOriginalDims(Int& num_row, Int& num_col) const {
  num_row = model_.m_orig();
  num_col = model_.n_orig();
}
const Info& Solver::getInfo() const { return info_; }
void Solver::getInteriorSolution(
    std::vector<double>& x, std::vector<double>& xl, std::vector<double>& xu,
    std::vector<double>& slack, std::vector<double>& y, std::vector<double>& zl,
    std::vector<double>& zu) const {
  // prepare and return solution with internal format

  if (!info_.ipx_used) {
    model_.postprocess(x, xl, xu, slack, y, zl, zu, *it_);
  } else {
    ipx_lps_.GetInteriorSolution(x.data(), xl.data(), xu.data(), slack.data(),
                                 y.data(), zl.data(), zu.data());
  }
}

Int Solver::getBasicSolution(std::vector<double>& x, std::vector<double>& slack,
                             std::vector<double>& y, std::vector<double>& z,
                             Int* cbasis, Int* vbasis) const {
  return ipx_lps_.GetBasicSolution(x.data(), slack.data(), y.data(), z.data(),
                                   cbasis, vbasis);
}

void Solver::chooseNumberOfCorrectors() {
  if (kMaxCorrectors > 0) {
    // Compute estimate of effort to factorise and solve

    // Effort to factorise depends on the number of flops
    double fact_effort = LS_->flops() + 100 * LS_->spops();

    // Effort to solve depends on the number of nonzeros of L multiplied by 2,
    // because there are two sweeps through L (forward and backward).
    double solv_effort = 2.0 * LS_->nz();

    // The factorise phase uses BLAS-3 and can be parallelised, the solve
    // phase uses BLAS-2 and cannot be parallelised. To account for this, the
    // factorisation effort is multiplied by a coefficient < 1, estimated
    // empirically.
    double alpha = 1.0 / 112.0;

    double ratio = alpha * fact_effort / solv_effort;
    if (std::isnan(ratio)) ratio = 0;

    // At each ipm iteration, there are up to (1+k) directions computed, where
    // k is the number of correctors. Each direction requires up (1+f) solves,
    // where f is the number of refinement steps. So, up to (1+k)(1+f) solves
    // are performed per iteration. However, not all refinement steps are used
    // all the time, so use f/2.
    // Therefore, we want (1+k)(1+f/2) < ratio.

    double thresh = ratio / (1.0 + kRefineMaxIter / 2.0) - 1;

    info_.correctors = std::floor(thresh);
    info_.correctors = std::max(info_.correctors, (Int)1);
    info_.correctors = std::min(info_.correctors, kMaxCorrectors);

  } else {
    info_.correctors = -kMaxCorrectors;
  }
}

bool Solver::solved() const { return info_.status >= kStatusTypeSolved; }
bool Solver::stopped() const {
  return info_.status > kStatusTypeStopped && info_.status < kStatusTypeFailed;
}
bool Solver::failed() const {
  return info_.status > kStatusTypeFailed && info_.status < kStatusTypeSolved;
}
bool Solver::statusAllowsCrossover() const {
  return getStatus1() == kStatusOptimal;
}
bool Solver::statusNeedsRefinement() const {
  return getStatus1() == kStatusNoProgress;
}
bool Solver::refinementIsOn() const {
  return options_.refine_with_ipx && !model_.qp();
}
bool Solver::crossoverIsOn() const {
  return options_.crossover == kHighsOnString && !model_.qp();
}

bool Solver::errorOrInterrupt() const {
  return getStatus() == kStatusTimeLimit ||
         getStatus() == kStatusUserInterrupt || info_.error;
}

// status1 should only be set if status2 is empty.
// status2 should only be set if status1 is optimal or imprecise.
// final status should consider status2 if non empty, or status 1 otherwise.

void Solver::setStatus1(Status status) {
  assert(getStatus2() == kStatusNotSet);
  status_phase1 = status;
}
void Solver::setStatus2(Status status) {
  assert(getStatus1() == kStatusOptimal || getStatus1() == kStatusImprecise);
  status_phase2 = status;
}
void Solver::setStatus(Status status) {
  if (getStatus1() == kStatusNotSet)
    setStatus1(status);
  else
    setStatus2(status);
}

Status Solver::getStatus1() const { return status_phase1; }
Status Solver::getStatus2() const { return status_phase2; }
Status Solver::getStatus() const {
  if (getStatus2() != kStatusNotSet)
    return getStatus2();
  else
    return getStatus1();
}

void Solver::finaliseStatus() {
  info_.status = getStatus();
  if (info_.error) info_.status = kStatusError;
  if (info_.status == kStatusNotSet) info_.status = kStatusMaxIter;
}

}  // namespace hipo
