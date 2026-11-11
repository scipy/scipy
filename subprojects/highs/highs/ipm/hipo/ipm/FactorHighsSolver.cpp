#include "FactorHighsSolver.h"

#include <limits>

#include "HighsExternalApi.h"
#include "Status.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "parallel/HighsParallel.h"
#include "util/stringutil.h"

namespace hipo {

FactorHighsSolver::FactorHighsSolver(KktMatrix& kkt, Options& options,
                                     const Model& model,
                                     const Regularisation& regul, Info& info,
                                     IpmData& record, const Logger& logger)
    : FH_{},
      kkt_{kkt},
      regul_{regul},
      info_{info},
      data_{record},
      logger_{logger},
      model_{model},
      options_{options} {
  FH_.setBlockSize(options.block_size);
  FH_.setLogger(&logger);
}

void FactorHighsSolver::clear() {
  valid_ = false;
  FH_.newIter();
}

// =========================================================================
// Analyse phase
// =========================================================================

Int FactorHighsSolver::analyseAS(Symbolic& S) {
  // Perform analyse phase of augmented system and return symbolic factorisation
  // in object S and the status.

  if (!kkt_.AS_ready) return kErrorAnalyse;

  const Int m = model_.A().num_row_;
  const Int n = model_.A().num_col_;

  std::vector<Int> pivot_signs(n + m, -1);
  for (Int i = 0; i < m; ++i) pivot_signs[n + i] = 1;

  logger_.printInfo("Performing AS analyse phase\n");

  Clock clock;
  Int status = chooseOrdering(kkt_.rowsAS, kkt_.ptrAS, pivot_signs, S,
                              ordering_AS_, "AS");
  info_.times[kAnalyseTime_AS] += clock.stop();

  return status;
}

Int FactorHighsSolver::analyseNE(Symbolic& S) {
  // Perform analyse phase of normal equations and return symbolic factorisation
  // in object S and the status. Structure of the matrix must be already
  // computed.

  if (!kkt_.NE_ready) return kErrorAnalyse;

  std::vector<Int> pivot_signs(model_.A().num_row_, 1);

  logger_.printInfo("Performing NE analyse phase\n");

  Clock clock;
  Int status = chooseOrdering(kkt_.rowsNE, kkt_.ptrNE, pivot_signs, S,
                              ordering_NE_, "NE");
  info_.times[kAnalyseTime_NE] += clock.stop();

  return status;
}

// =========================================================================
// Factorise phase
// =========================================================================

Int FactorHighsSolver::factorAS(const std::vector<double>& scaling) {
  assert(!this->valid_);
  if (!kkt_.isAS()) return kErrorFactorise;

  kkt_.buildASvalues(scaling);

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  Clock clock;
  if (FH_.factorise(kkt_.S, kkt_.n(), kkt_.nz(), kkt_.rowsAS.data(),
                    kkt_.ptrAS.data(), kkt_.valAS.data()))
    return kErrorFactorise;
  info_.times[kFactoriseTime] += clock.stop();
  info_.factor_number++;

  this->valid_ = true;
  return kOk;
}

Int FactorHighsSolver::factorNE(const std::vector<double>& scaling) {
  assert(!this->valid_);
  if (!kkt_.isNE()) return kErrorFactorise;

  kkt_.buildNEvalues(scaling);

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  Clock clock;
  if (FH_.factorise(kkt_.S, kkt_.n(), kkt_.nz(), kkt_.rowsNE.data(),
                    kkt_.ptrNE.data(), kkt_.valNE.data()))
    return kErrorFactorise;
  info_.times[kFactoriseTime] += clock.stop();
  info_.factor_number++;

  this->valid_ = true;
  return kOk;
}

// =========================================================================
// Solve phase
// =========================================================================

Int FactorHighsSolver::solveAS(const std::vector<double>& rhs_x,
                               const std::vector<double>& rhs_y,
                               std::vector<double>& lhs_x,
                               std::vector<double>& lhs_y) {
  assert(this->valid_);
  assert(kkt_.isAS());

  Int n = rhs_x.size();

  Clock clock;
  as_buffer_.resize(rhs_x.size() + rhs_y.size());
  std::copy(rhs_x.begin(), rhs_x.end(), as_buffer_.begin());
  std::copy(rhs_y.begin(), rhs_y.end(), as_buffer_.begin() + n);
  info_.times[kInsertSplitTime] += clock.stop();

  clock.start();
  if (FH_.solve(as_buffer_.data())) return kErrorSolve;
  info_.times[kSolveTime] += clock.stop();
  info_.solve_number++;

  data_.back().num_solves++;

  clock.start();
  lhs_x.assign(as_buffer_.begin(), as_buffer_.begin() + n);
  lhs_y.assign(as_buffer_.begin() + n, as_buffer_.end());
  info_.times[kInsertSplitTime] += clock.stop();

  return kOk;
}

Int FactorHighsSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  assert(this->valid_);
  assert(kkt_.isNE());

  lhs = rhs;

  Clock clock;
  if (FH_.solve(lhs.data())) return kErrorSolve;
  info_.times[kSolveTime] += clock.stop();
  info_.solve_number++;

  data_.back().num_solves++;

  return kOk;
}

// =========================================================================
// Automatic selections
// =========================================================================

Int FactorHighsSolver::setup() {
  if (kkt_.S.empty()) {
    Clock clock;
    setParallelBeforeSymbolic();
    if (Int status = setNla()) return status;
    setParallelAfterSymbolic();
    printParallel();
    info_.times[kAnalyseTime] += clock.stop();

    if (!options_.timeless_log) {
      std::stringstream log_stream;
      log_stream << textline("Analyse time:")
                 << fix(info_.times[kAnalyseTime], 0, 2) << '\n';
      logger_.print(log_stream.str().c_str());
    }

    kkt_.S.print(logger_, logger_.debug(1));

    if (kkt_.S.storage() > kParallelLargeStorageGB * 1024 * 1024 * 1024) {
      logger_.printw("Large amount of memory required\n");
    }

    logger_.print("\n");
  }

  return kOk;
}

Int FactorHighsSolver::chooseNla() {
  Symbolic symb_NE{};
  Symbolic symb_AS{};
  bool failure_NE = false;
  bool failure_AS = false;
  bool overflow_NE = false;
  bool overflow_AS = false;

  const bool AS_possible = true;
  const bool NE_possible = !model_.nonSeparableQp();

  const bool expect_AS_much_cheaper =
      model_.nzNElb() > model_.nzAS() * kSkipSystemNzBoundsRatio;
  const bool expect_NE_much_cheaper =
      model_.nzAS() > model_.nzNEub() * kSkipSystemNzBoundsRatio;

  const bool skip_AS = NE_possible && expect_NE_much_cheaper;
  const bool skip_NE = AS_possible && expect_AS_much_cheaper;

  bool allow_skip_AS = true;
  bool allow_skip_NE = true;

  auto run_structure_NE = [&]() {
    if (!NE_possible || (skip_NE && allow_skip_NE)) {
      failure_NE = true;
      logger_.printInfo("NE skipped\n");
    } else {
      Int status = kkt_.buildNEstructure();
      if (status) {
        failure_NE = true;
        if (status == kErrorOverflow) {
          logger_.printInfo("Integer overflow forming NE matrix\n");
          overflow_NE = true;
        }
        return;
      }
    }
  };

  auto run_analyse_NE = [&]() {
    if (failure_NE) return;
    Int NE_status = analyseNE(symb_NE);
    if (NE_status) failure_NE = true;
  };

  auto run_analyse_AS = [&]() {
    if (!AS_possible || (skip_AS && allow_skip_AS)) {
      failure_AS = true;
      logger_.printInfo("AS skipped\n");
    } else {
      Int AS_status = kkt_.buildASstructure();
      if (!AS_status) AS_status = analyseAS(symb_AS);
      if (AS_status) failure_AS = true;
      if (AS_status == kErrorOverflow) {
        logger_.printInfo("Integer overflow forming AS matrix\n");
        overflow_AS = true;
      }

      // If NE has more nonzeros than the factor of AS, then it's likely that AS
      // will be preferred, so stop computation of NE.
      Int64 NE_nz_limit = symb_AS.nz() * kSystemSymbNzMult;
      if (failure_AS || NE_nz_limit > kHighsIInf) NE_nz_limit = kHighsIInf;
      kkt_.NE_nz_limit.store(NE_nz_limit, std::memory_order_relaxed);
    }
  };

  // In parallel, run AS analyse and build NE structure. NE analyse runs only
  // after AS analyse is finished, so that it can be skipped based on the number
  // of nz of NE matrix and AS factor.
  if (options_.getParallel(ParallelTechnique::kAnalyse)) {
    highs::parallel::TaskGroup tg;
    tg.spawn([&]() { run_analyse_AS(); });
    tg.spawn([&]() { run_structure_NE(); });
    tg.taskWait();
  } else {
    run_analyse_AS();
    run_structure_NE();
  }

  // if NE was skipped but AS failed, use NE
  if (skip_NE && failure_AS) {
    allow_skip_NE = false;
    run_structure_NE();
  }

  run_analyse_NE();

  // if AS was skipped but NE failed, use AS
  if (skip_AS && failure_NE) {
    allow_skip_AS = false;
    run_analyse_AS();
  }

  Int status = kOk;

  std::stringstream log_stream;

  if (failure_NE && !failure_AS) {
    options_.nla = kHipoAugmentedString;
    log_stream << textline("Newton system:") << "AS preferred (NE failed)\n";
  } else if (failure_AS && !failure_NE) {
    options_.nla = kHipoNormalEqString;
    log_stream << textline("Newton system:") << "NE preferred (AS failed)\n";
  } else if (failure_AS && failure_NE) {
    if (overflow_AS && overflow_NE)
      status = kErrorOverflow;
    else
      status = kErrorAnalyse;

    logger_.printe("Both NE and AS failed analyse phase\n");
  } else {
    // Total number of operations, given by dense flops and sparse indexing
    // operations, weighted with an empirical factor
    double ops_NE = symb_NE.flops() + symb_NE.spops() * kSystemSpopsWeight;
    double ops_AS = symb_AS.flops() + symb_AS.spops() * kSystemSpopsWeight;

    double sn_size_NE = (double)symb_NE.size() / symb_NE.sn();
    double sn_size_AS = (double)symb_AS.size() / symb_AS.sn();

    double ratio_ops = ops_NE / ops_AS;
    double ratio_sn = sn_size_AS / sn_size_NE;

    bool NE_much_more_expensive = ratio_ops > kSystemRatioOpsThresh;
    bool AS_not_too_expensive = ratio_ops > 1.0 / kSystemRatioOpsThresh;
    bool sn_AS_larger_than_NE = ratio_sn > kSystemRatioSnThresh;

    if (NE_much_more_expensive ||
        (sn_AS_larger_than_NE && AS_not_too_expensive)) {
      options_.nla = kHipoAugmentedString;
      log_stream << textline("Newton system:") << "AS preferred\n";
    } else {
      options_.nla = kHipoNormalEqString;
      log_stream << textline("Newton system:") << "NE preferred\n";
    }
  }

  logger_.print(log_stream.str().c_str());

  if (status == kOk) {
    if (options_.nla == kHipoAugmentedString) {
      kkt_.S = std::move(symb_AS);
      kkt_.freeNEmemory();
    } else {
      kkt_.S = std::move(symb_NE);
      kkt_.freeASmemory();
    }
  }

  return status;
}

Int FactorHighsSolver::chooseOrdering(const std::vector<Int>& rows,
                                      const std::vector<Int>& ptr,
                                      const std::vector<Int>& signs,
                                      Symbolic& S, std::string& ordering,
                                      const std::string& nla) {
  // Run analyse phase.
  // - If ordering is "amd", "metis", "rcm" run only the ordering requested.
  // - If ordering is "choose", run "amd", "metis", and choose the best.

  std::vector<std::string> orderings_to_try;
  if (options_.ordering != kHighsChooseString)
    orderings_to_try.push_back(options_.ordering);
  else {
    orderings_to_try.push_back(kHipoAmdString);
    orderings_to_try.push_back(kHipoMetisString);

    // rcm is much worse in general, so no point in trying for now
  }

  const Int k = orderings_to_try.size();

  std::vector<HighsBool> failure(k, false);

  if (nla == "NE") {
    if (ptr.back() >= kkt_.NE_nz_limit.load(std::memory_order_relaxed)) {
      logger_.printInfo("NE interrupted\n");
      return kErrorAnalyse;
    }
  }

  // compute full-format matrix without diagonal entries
  std::vector<Int> full_ptr, full_rows;
  fullFromLower(ptr.size() - 1, rows.size(), ptr.data(), rows.data(), full_ptr,
                full_rows);
  Int n = full_ptr.size() - 1;
  std::vector<Int> perm(n), iperm(n);

  std::vector<std::vector<Int>> permutations(k, std::vector<Int>(n));

  std::vector<Symbolic> symbolics(k, S);

  auto run_ordering_and_analyse = [&](Int i) {
    logger_.printInfo("Running %s for %s\n", orderings_to_try[i].c_str(),
                      nla.c_str());

    if (n > 0) {
      if (orderings_to_try[i] == kHipoMetisString) {
        Int status = FH_.reorderMetis(n, rows.size(), full_rows.data(),
                                      full_ptr.data(), permutations[i].data(),
                                      true, options_.random_seed);
        if (status) failure[i] = true;

      } else if (orderings_to_try[i] == kHipoAmdString) {
        Int status =
            FH_.reorderAmd(n, rows.size(), full_rows.data(), full_ptr.data(),
                           permutations[i].data(), true);
        if (status) failure[i] = true;

      } else if (orderings_to_try[i] == kHipoRcmString) {
        Int status =
            FH_.reorderRcm(n, rows.size(), full_rows.data(), full_ptr.data(),
                           permutations[i].data(), true);
        if (status) failure[i] = true;

      } else {
        assert(1 == 0);
      }
    }

    logger_.printInfo("Finished %s for %s\n", orderings_to_try[i].c_str(),
                      nla.c_str());

    if (failure[i]) {
      logger_.printInfo("Error with %s for %s\n", orderings_to_try[i].c_str(),
                        nla.c_str());
      return;
    }

    failure[i] =
        FH_.analyse(symbolics[i], ptr.size() - 1, rows.size(), rows.data(),
                    ptr.data(), signs.data(), permutations[i].data());

    if (failure[i] && logger_.debug(2)) {
      logger_.print("Failed symbolic:");
      symbolics[i].print(logger_, true);
    }
  };

  const bool parallel_ordering =
      nla == "NE" ? options_.getParallel(ParallelTechnique::kOrderNE)
                  : options_.getParallel(ParallelTechnique::kOrderAS);

  if (parallel_ordering) {
    highs::parallel::for_each(
        0, k, [&](Int start, Int end) { run_ordering_and_analyse(start); }, 1);
  } else {
    for (Int i = 0; i < k; ++i) run_ordering_and_analyse(i);
  }

  Int num_success = 0;
  for (bool b : failure) {
    if (!b) ++num_success;
  }

  if (num_success > 0) {
    for (Int i = 0; i < k; ++i) {
      if (!failure[i])
        logger_.printInfo(
            "%20s for %s: %.2e %.2f\n", orderings_to_try[i].c_str(),
            nla.c_str(), symbolics[i].flops(),
            static_cast<double>(symbolics[i].size()) / symbolics[i].sn());
    }

    // find the ordering with best flops
    double best_flops = kHighsInf;
    for (Int i = 0; i < k; ++i) {
      if (!failure[i] && symbolics[i].flops() < best_flops) {
        best_flops = symbolics[i].flops();
      }
    }

    // find orderings with flops within kOrderingFlopsThresh of the best
    std::vector<Int> consider;
    for (Int i = 0; i < k; ++i) {
      if (!failure[i] &&
          symbolics[i].flops() <= kOrderingFlopsThresh * best_flops) {
        consider.push_back(i);
      }
    }

    // among them, find the one with largest supernodes
    double best_sn_size = 0.0;
    Int chosen = -1;
    for (Int i : consider) {
      double sn_avg_size =
          static_cast<double>(symbolics[i].size()) / symbolics[i].sn();
      if (sn_avg_size > best_sn_size) {
        best_sn_size = sn_avg_size;
        chosen = i;
      }
    }

    // fix selection if one or more require too much memory
    const double bytes_thresh = kParallelLargeStorageGB * 1024 * 1024 * 1024;
    double best_memory = kHighsInf;
    Int ind_best_memory = -1;
    for (Int i = 0; i < k; ++i) {
      if (symbolics[i].storage() < best_memory) {
        best_memory = symbolics[i].storage();
        ind_best_memory = i;
      }
    }

    if (n == 0) chosen = 0;

    assert(chosen >= 0 && chosen < k);

    S = std::move(symbolics[chosen]);
    ordering = orderings_to_try[chosen];
  }

  return num_success > 0 ? kOk : kErrorAnalyse;
}

Int FactorHighsSolver::setNla() {
  std::stringstream log_stream;

  if (options_.nla == kHipoNormalEqString && model_.nonSeparableQp()) {
    logger_.printw("Normal equations not available for non-separable QP\n");
    options_.nla = kHighsChooseString;
  }

  if (options_.nla == kHipoAugmentedString) {
    Int status = kkt_.buildASstructure();
    if (!status) status = analyseAS(kkt_.S);
    if (status == kErrorOverflow) {
      logger_.printe("AS requested, integer overflow\n");
      return kErrorOverflow;
    } else if (status) {
      logger_.printe("AS requested, failed analyse phase\n");
      return kErrorAnalyse;
    }
    log_stream << textline("Newton system:") << "AS requested\n";

  } else if (options_.nla == kHipoNormalEqString) {
    Int status = kkt_.buildNEstructure();
    if (!status) status = analyseNE(kkt_.S);
    if (status == kErrorOverflow) {
      logger_.printe("NE requested, integer overflow\n");
      return kErrorOverflow;
    } else if (status) {
      logger_.printe("NE requested, failed analyse phase\n");
      return kErrorAnalyse;
    }
    log_stream << textline("Newton system:") << "NE requested\n";

  } else if (options_.nla == kHighsChooseString) {
    if (Int status = chooseNla()) return status;

  } else
    assert(1 == 0);

  if (logger_.debug(1))
    log_stream << textline("Ordering:")
               << (options_.nla == kHipoAugmentedString ? ordering_AS_
                                                        : ordering_NE_)
               << '\n';
  logger_.print(log_stream.str().c_str());

  return kOk;
}

void FactorHighsSolver::setParallelBeforeSymbolic() {
  const bool parallel_analyse_default = true;
  options_.chooseParallel(ParallelTechnique::kAnalyse,
                          parallel_analyse_default);

  const bool parallel_order_NE_default = true;
  options_.chooseParallel(ParallelTechnique::kOrderNE,
                          parallel_order_NE_default);

  const bool parallel_order_AS_default = true;
  options_.chooseParallel(ParallelTechnique::kOrderNE,
                          parallel_order_AS_default);

  const double A_nz_per_col = (double)model_.A().numNz() / model_.A().num_col_;
  const double A_nz_per_row = (double)model_.A().numNz() / model_.A().num_row_;
  const bool A_is_dense = A_nz_per_col > kParallelNEnzPerColThresh ||
                          A_nz_per_row > kParallelNEnzPerRowThresh;
  const bool A_is_large = model_.A().num_row_ > kParallelNEsizeThresh ||
                          model_.A().num_col_ > kParallelNEsizeThresh;

  const bool parallel_NE_struct_default = A_is_dense && A_is_large;
  options_.chooseParallel(ParallelTechnique::kNEStruct,
                          parallel_NE_struct_default);

  const bool parallel_NE_values_default = A_is_large;
  options_.chooseParallel(ParallelTechnique::kNEValues,
                          parallel_NE_values_default);
}

void FactorHighsSolver::setParallelAfterSymbolic() {
  bool parallel_tree = false;
  bool parallel_node = false;

  if (is_substring_case(HighsExtras::blas::getInfo()->provider, "apple")) {
    // Blas on Apple do not work well with parallel_node, but parallel_tree
    // seems to always be beneficial.
    parallel_node = false;
    parallel_tree = true;
  } else {
    // Otherwise, parallel_node is active because it is triggered only if the
    // frontal matrix is large enough anyway.
    parallel_node = true;

    // parallel_tree instead is chosen with a heuristic

    double tree_speedup = kkt_.S.flops() / kkt_.S.critops();
    double sn_size = (double)kkt_.S.size() / kkt_.S.sn();

    bool enough_sn = kkt_.S.sn() > kParallelMinNumberSn;
    bool enough_flops = kkt_.S.flops() > kParallelLargeFlopsThresh;
    bool speedup_is_large = tree_speedup > kParallelLargeSpeedupThresh;
    bool sn_are_large = sn_size > kParallelLargeSnThresh;
    bool sn_are_not_small = sn_size > kParallelSmallSnThresh;

    // parallel_tree is active if the supernodes are large, or if there is a
    // large expected speedup and the supernodes are not too small, provided
    // that the number of flops and supernodes is not too small.
    if (enough_sn && enough_flops &&
        (sn_are_large || (speedup_is_large && sn_are_not_small))) {
      parallel_tree = true;
    }
  }

  // If serial memory is too large, switch off tree parallelism to avoid
  // running out of memory
  double num_GB = kkt_.S.storage() / 1024 / 1024 / 1024;
  if (num_GB > kParallelLargeStorageGB) {
    parallel_tree = false;
  }

  // switch off tree parallelism if depth of recursion is too large
  if (kkt_.S.depth() > kParallelMaxTreeDepth) parallel_tree = false;

  options_.chooseParallel(ParallelTechnique::kTree, parallel_tree);
  options_.chooseParallel(ParallelTechnique::kNode, parallel_node);

  // choose parallelism for solve phase
  bool parallel_forward = false;
  bool parallel_backward = false;
  bool parallel_diag = false;
  if (kkt_.S.size() > kParallelSolveMinSize) {
    parallel_diag = true;

    if (kkt_.S.solveTreeSpeedup() > kParallelForwardMinSpeedup)
      parallel_forward = true;

    if (kkt_.S.solveTreeSpeedup() > kParallelBackwardMinSpeedup)
      parallel_backward = true;
  }

  options_.chooseParallel(ParallelTechnique::kForwardSolve, parallel_forward);
  options_.chooseParallel(ParallelTechnique::kBackwardSolve, parallel_backward);
  options_.chooseParallel(ParallelTechnique::kDiagonalSolve, parallel_diag);

  for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i)
    assert(options_.parallel[i] == ParallelType::kOn ||
           options_.parallel[i] == ParallelType::kOff);

  FH_.setParallel(options_.getParallel(ParallelTechnique::kTree),
                  options_.getParallel(ParallelTechnique::kNode));

  FH_.setParallelSolve(options_.getParallel(ParallelTechnique::kForwardSolve),
                       options_.getParallel(ParallelTechnique::kBackwardSolve),
                       options_.getParallel(ParallelTechnique::kDiagonalSolve));
}

void FactorHighsSolver::printParallel() const {
  std::stringstream log_stream;
  log_stream
      << textline("Parallelism:")
      << (options_.getParallel(ParallelTechnique::kAnalyse) ? "A" : "_")
      << (options_.getParallel(ParallelTechnique::kOrderNE) ? "O" : "_")
      << (options_.getParallel(ParallelTechnique::kOrderAS) ? "O" : "_") << "|"
      << (options_.getParallel(ParallelTechnique::kNEStruct) ? "S" : "_")
      << (options_.getParallel(ParallelTechnique::kNEValues) ? "V" : "_") << "|"
      << (options_.getParallel(ParallelTechnique::kTree) ? "T" : "_")
      << (options_.getParallel(ParallelTechnique::kNode) ? "N" : "_") << "|"
      << (options_.getParallel(ParallelTechnique::kForwardSolve) ? "F" : "_")
      << (options_.getParallel(ParallelTechnique::kDiagonalSolve) ? "D" : "_")
      << (options_.getParallel(ParallelTechnique::kBackwardSolve) ? "B" : "_")
      << '\n';
  logger_.printInfo(log_stream.str().c_str());
}

// =========================================================================
// Other stuff
// =========================================================================

double FactorHighsSolver::flops() const { return kkt_.S.flops(); }
double FactorHighsSolver::spops() const { return kkt_.S.spops(); }
double FactorHighsSolver::nz() const { return (double)kkt_.S.nz(); }
void FactorHighsSolver::getReg(std::vector<double>& reg) {
  FH_.getRegularisation(reg.data());
}

}  // namespace hipo
