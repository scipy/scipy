/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_SETTINGS_HPP__
#define __SRC_LIB_SETTINGS_HPP__

#include "eventhandler.hpp"
#include "qpconst.hpp"
#include "statistics.hpp"

enum class RatiotestStrategy { TwoPass, Textbook };

enum class PricingStrategy { SteepestEdge, DantzigWolfe, Devex };

enum class Phase1Strategy { HIGHS, QUASS, BOUNDED };

struct Settings {
  RatiotestStrategy ratiotest = RatiotestStrategy::TwoPass;
  double ratiotest_t = 1e-9;
  double ratiotest_d = 1e-8;
  double primal_feasibility_tolerance = 1e-7;

  PricingStrategy pricing = PricingStrategy::Devex;

  double pnorm_zero_threshold =
      1e-11;  // if ||p|| < this threshold, p is determined to not be an
              // improving search direction
  double improvement_zero_threshold =
      1e-4;  // if p^t gradient < this threshold, p is determined to not be an
             // improving search direction
  double d_zero_threshold = 1e-12;  // minimal value for pivot, will declare
                                    // degeneracy if no larger pivot is found
  double lambda_zero_threshold =
      1e-9;  // used for pricing / optimality checking
  double pQp_zero_threshold =
      1e-7;  // if p'Qp < this, p is determined to not have curvature, a
             // simplex-like iteration is performed.

  bool hessian_regularization_on =
      false;  // if true, a small multiple of the identity matrix will be added
              // to the Hessian.
              //
              // This is always false, so regularization in
              // regularize(Runtime& rt) never happens, but
              // regularization in solveqp is _always_ performed
              //
              // This explains the "perturbed" solutions of problems
              // in TestQpSolver.cpp
  double hessian_regularization_value =
      1e-7;  // multiple of identity matrix added to Hessian in case of
             // regularization

  bool allow_hot_start = false;
  Phase1Strategy phase1strategy = Phase1Strategy::HIGHS;
  bool phase1movefreevarsbasic = false;
  bool phase1boundfreevars = false;

  HighsInt reportingfequency = 1;
  Eventhandler<HighsInt&> iteration_log_header;
  Eventhandler<Statistics&> iteration_log;
  Eventhandler<QpModelStatus&> qp_model_status_log;
  Eventhandler<HighsInt&> nullspace_limit_log;
  Eventhandler<std::pair<HighsInt, double>&> degeneracy_fail_log;

  HighsInt nullspace_limit = 4000;

  HighsInt reinvertfrequency = 1000;
  HighsInt gradientrecomputefrequency = 100;
  HighsInt reducedgradientrecomputefrequency = kHighsIInf;
  HighsInt reducedhessianrecomputefrequency = kHighsIInf;

  HighsInt iteration_limit = kHighsIInf;
  double time_limit = kHighsInf;

  bool rowscaling = true;
  bool varscaling = true;

  bool perturbation = false;
};

#endif
