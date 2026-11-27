/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSolverSelect.h
 * @brief Choose an LP algorithm (simplex / IPX / HiPO) from cheap
 *        structural features of the (presolved) model.
 *
 * Regenerating the models in HighsSolverSelectModels.h
 *
 * The models are fitted offline by the scripts in select/ of the HiGHS
 * test repo:
 *
 * 1. presolve_set.sh writes the presolved model of each instance.
 * 2. run_solvers.sh <version> LP|MIP times dual simplex, IPX and HiPO on
 *    the presolved LPs and on the LP relaxations of presolved MIPs.
 * 3. extract_features/ applies LP presolve, as solveLp sees the LP, and
 *    writes computeLpFeatures() of each model to a CSV.
 * 4. collect_results.py --version <version> writes select.csv.
 * 5. visualize.py --cpp-out <HiGHS>/highs/util/HighsSolverSelectModels.h
 *    fits the gradient boosted trees and writes the header.
 *
 * The models index highsLpFeatureVector() by position: changing the order
 * of the features, or what computeLpFeatures() computes, needs a refit.
 * The test solver-select-feature-order fixes the order.
 */
#ifndef UTIL_HIGHS_SOLVER_SELECT_H_
#define UTIL_HIGHS_SOLVER_SELECT_H_

#include <string>
#include <utility>
#include <vector>

#include "util/HighsInt.h"

class HighsLp;

enum class HighsSolverSelect { kDualSimplex, kIpx, kHipo };

// Effect of the options solver_select_strategy and
// solver_select_require_basis
//
// solver_select_strategy is a SolverSelectStrategy. It only acts where the
// user has left the choice of LP solver open:
//
// * LP (solveLp): when solver is "choose", the strategy is positive and
//   run_centring is false. Any other value of solver is used as given. The
//   choice is made for the LP after presolve. If there is a useful basis,
//   or the LP has no constraints or no nonzeros, then dual simplex is used,
//   whatever the strategy
//
// * MIP (HighsLpRelaxation::run): solver plays no part. An LP with a valid
//   basis is solved by simplex. Otherwise mip_lp_solver "simplex", "ipm",
//   "ipx" or "hipo" is used as given, whatever the strategy. Only
//   mip_lp_solver "choose" defers to a positive strategy, and not when
//   simplex_strategy is primal. LPs solved using mip_ipm_solver are
//   unaffected
//
// If HiPO is selected but isn't available in the build, IPX is used
//
// A selected IPM solver has fallbacks that an explicit one doesn't
// (solveLpWithSelection): IPX if HiPO runs out of memory, and simplex if
// IPM fails, returns an imprecise result, can't tell unbounded from
// infeasible (unless allow_unbounded_or_infeasible), or yields no basis
// when one is needed. Simplex isn't run again if solveLp has already used
// it to clean up. In the MIP solver, IPM is limited to 200 iterations for
// LPs with a simplex iteration limit
//
// Presolve precedes selection, so if the strategy may choose IPM then
// presolve allows for dual postsolve, as it does for solver "ipm"
//
// solver_select_require_basis makes a basis part of the result: a claim
// of optimality, infeasibility or unboundedness by a selected IPM solver
// without a basis, or IPM reaching its iteration limit, is followed by a
// simplex solve, and simplex is selected if run_crossover is "off". The
// MIP solver sets it,
// since it needs a basis. When false, a basis is only needed for an
// optimal solution when run_crossover is "on"
//
// Only strategies 0, 1 and 2 are documented for users. The rest are for
// development, but can still be set
enum class SolverSelectStrategy {
  kSimplexOnly = 0,
  kColsRowsNnz = 1,
  kRankedFeatures = 2,
  kHeuristic = 3,
  // Always choose one solver, to compare them
  kDualSimplex = 4,
  kIpx = 5,
  kHipo = 6,
  // Dual simplex below the fitted models' min_num_nz, otherwise always
  // one IPM solver
  kIpxIfLarge = 7,
  kHipoIfLarge = 8,
};

// Whether LPs go through solver selection with this strategy
inline bool useSolverSelection(const HighsInt strategy) {
  return strategy > static_cast<HighsInt>(SolverSelectStrategy::kSimplexOnly);
}

// Whether solver selection with this strategy may choose IPM. Presolve
// must then allow for dual postsolve
inline bool solverSelectionMayChooseIpm(const HighsInt strategy) {
  return useSolverSelection(strategy) &&
         strategy != static_cast<HighsInt>(SolverSelectStrategy::kDualSimplex);
}

// All features are doubles, since that's how the fitted models read them
struct HighsLpFeatures {
  // --- Raw size ---
  double num_row = 0;
  double num_col = 0;
  double num_nz = 0;
  double num_integer_col = 0;  // Set only by offline feature extraction

  // --- How constrained is the problem? (barrier likes an interior) ---
  double relative_num_equalities = 0.0;           // equality rows / #rows
  double relative_num_inequalities = 0.0;         // one-sided rows / #rows
  double relative_num_ranged_rows = 0.0;          // l < u, both finite
  double relative_num_free_rows = 0.0;            // -inf..inf rows / #rows
  double relative_num_cols_without_upper = 0.0;   // u_j = +inf / #cols
  double relative_num_cols_without_lower = 0.0;   // l_j = -inf / #cols
  double relative_num_free_cols = 0.0;            // -inf..inf cols / #cols
  double relative_num_boxed_cols = 0.0;           // both bounds finite
  double relative_num_singly_bounded_cols = 0.0;  // exactly one finite
  double relative_num_fixed_cols = 0.0;           // l_j == u_j

  // --- How degenerate is the problem? (simplex can get stuck) ---
  double relative_num_almost_identical_nonzeros = 0.0;
  double relative_num_identical_obj_values = 0.0;  // repeated c_j
  double relative_num_identical_rhs_values = 0.0;  // repeated rhs

  // --- How stable is the problem? (conditioning) ---
  double relative_max_matrix_coefficient = 0.0;  // max|a| / min|a|
  double relative_max_obj_coefficient = 0.0;     // max|c| / min|c|
  double relative_max_rhs_value = 0.0;           // max|rhs| / min|rhs|
  double max_abs_matrix_coefficient = 0.0;
  double min_abs_matrix_coefficient = 0.0;
  double max_abs_obj_coefficient = 0.0;
  double min_abs_obj_coefficient = 0.0;
  double max_abs_rhs_value = 0.0;
  double min_abs_rhs_value = 0.0;

  // --- How dense is the problem? (simplex thrives on sparsity) ---
  double matrix_nonzero_density = 0.0;   // #nz / (#rows * #cols)
  double avg_nonzeros_per_col = 0.0;     // #nz / #cols
  double avg_nonzeros_per_row = 0.0;     // #nz / #rows
  double relative_num_dense_rows = 0.0;  // "dense" rows / #rows
  double max_nonzeros_in_row = 0.0;
  double max_nonzeros_in_col = 0.0;
};

HighsLpFeatures computeLpFeatures(const HighsLp& lp);

// Flatten `features` into (name, value) pairs, in the stable order that
// the fitted models index.  `highsLpFeatureNames()` returns just the names
// (same order) so a header row can be emitted without a features instance.
std::vector<std::pair<std::string, double>> highsLpFeatureVector(
    const HighsLpFeatures& features);
std::vector<std::string> highsLpFeatureNames();

// The solver that the strategy chooses for the LP
HighsSolverSelect selectSolver(const HighsLp& lp,
                               const SolverSelectStrategy strategy);

#endif /* UTIL_HIGHS_SOLVER_SELECT_H_ */
