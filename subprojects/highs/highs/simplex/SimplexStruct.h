/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/SimplexStruct.h
 * @brief Structs for HiGHS simplex solvers
 */
#ifndef SIMPLEX_SIMPLEXSTRUCT_H_
#define SIMPLEX_SIMPLEXSTRUCT_H_

#include <cstdint>
#include <vector>

// #include "lp_data/HighsLp.h"
#include "lp_data/HConst.h"
#include "simplex/SimplexConst.h"

struct SimplexBasis {
  // The basis for the simplex method consists of basicIndex,
  // nonbasicFlag and nonbasicMove. If HighsSimplexStatus has_basis
  // is true then it is assumed that basicIndex_ and nonbasicFlag_ are
  // self-consistent and correspond to the dimensions of an associated
  // HighsLp, but the basis matrix B is not necessarily nonsingular.
  std::vector<HighsInt> basicIndex_;
  std::vector<int8_t> nonbasicFlag_;
  std::vector<int8_t> nonbasicMove_;
  //  std::vector<double> debug_dual;
  uint64_t hash;
  HighsInt debug_id = -1;
  HighsInt debug_update_count = -1;
  std::string debug_origin_name = "";
  void clear();
  void setup(const HighsInt num_col, const HighsInt num_row);
};

struct HighsSimplexStatus {
  // Status of LP solved by the simplex method and its data
  bool initialised_for_new_lp = false;
  bool is_dualized = false;
  bool is_permuted = false;
  bool initialised_for_solve = false;
  bool has_basis = false;      // The simplex LP has a valid simplex basis
  bool has_ar_matrix = false;  // HEkk has the row-wise matrix
  bool has_nla = false;        // SimplexNla is set up
  bool has_dual_steepest_edge_weights = false;  // The DSE weights are known
  bool has_invert =
      false;  // The representation of B^{-1} corresponds to the current basis
  bool has_fresh_invert = false;  // The representation of B^{-1} corresponds to
                                  // the current basis and is fresh
  bool has_fresh_rebuild = false;  // The data are fresh from rebuild
  bool has_dual_objective_value =
      false;  // The dual objective function value is known
  bool has_primal_objective_value =
      false;  // The dual objective function value is known
};

struct HighsSimplexInfo {
  // Simplex information regarding primal solution, dual solution and
  // objective for this Highs Model Object. This is information which
  // should be retained from one run to the next in order to provide
  // hot starts.
  //
  // Part of working model which are assigned and populated as much as
  // possible when a model is being defined

  // workCost: Originally just costs from the model but, in solve(), may
  // be perturbed or set to alternative values in Phase I??
  //
  // workDual: Values of the dual variables corresponding to
  // workCost. Latter not known until solve() is called since B^{-1}
  // is required to compute them.
  //
  // workShift: Values added to workCost in order that workDual
  // remains feasible, thereby remaining dual feasible in phase 2
  //
  std::vector<double> workCost_;
  std::vector<double> workDual_;
  std::vector<double> workShift_;

  // workLower/workUpper: Originally just lower (upper) bounds from
  // the model but, in solve(), may be perturbed or set to
  // alternative values in Phase I??
  //
  // workRange: Distance between lower and upper bounds
  //
  // workValue: Values of the nonbasic variables corresponding to
  // workLower/workUpper and the basis. Always known.
  //
  std::vector<double> workLower_;
  std::vector<double> workUpper_;
  std::vector<double> workRange_;
  std::vector<double> workValue_;
  std::vector<double> workLowerShift_;
  std::vector<double> workUpperShift_;
  //
  // baseLower/baseUpper/baseValue: Lower and upper bounds on the
  // basic variables and their values. Latter not known until solve()
  // is called since B^{-1} is required to compute them.
  //
  std::vector<double> baseLower_;
  std::vector<double> baseUpper_;
  std::vector<double> baseValue_;
  //
  // Vectors of random reals for column cost perturbation, a random
  // permutation of all indices for CHUZR and a random permutation of
  // column indices for permuting the columns
  std::vector<double> numTotRandomValue_;
  std::vector<HighsInt> numTotPermutation_;
  std::vector<HighsInt> numColPermutation_;

  std::vector<HighsInt> devex_index_;

  // Records of the row chosen by dual simplex or column chosen by
  // primal simplex, plus the pivot values - since last revinversion
  std::vector<HighsInt> index_chosen_;
  std::vector<double> pivot_;

  // Data for backtracking in the event of a singular basis
  HighsInt phase1_backtracking_test_done = false;
  HighsInt phase2_backtracking_test_done = false;
  bool backtracking_ = false;
  bool valid_backtracking_basis_ = false;
  SimplexBasis backtracking_basis_;
  HighsInt backtracking_basis_costs_shifted_;
  HighsInt backtracking_basis_costs_perturbed_;
  HighsInt backtracking_basis_bounds_shifted_;
  HighsInt backtracking_basis_bounds_perturbed_;
  std::vector<double> backtracking_basis_workShift_;
  std::vector<double> backtracking_basis_workLowerShift_;
  std::vector<double> backtracking_basis_workUpperShift_;
  std::vector<double> backtracking_basis_edge_weight_;

  // Options from HighsOptions for the simplex solver
  HighsInt simplex_strategy;
  HighsInt dual_edge_weight_strategy;
  HighsInt primal_edge_weight_strategy;
  HighsInt price_strategy;

  double dual_simplex_cost_perturbation_multiplier;
  double primal_simplex_phase1_cost_perturbation_multiplier = 1;
  double primal_simplex_bound_perturbation_multiplier;
  double factor_pivot_threshold;
  HighsInt update_limit;

  // Simplex control parameters from HSA
  HighsInt control_iteration_count0;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  double col_steepest_edge_density;
  double col_basic_feasibility_change_density;
  double row_basic_feasibility_change_density;
  double col_BFRT_density;
  double primal_col_density;
  double dual_col_density;
  // For control of switch from DSE to Devex in dual simplex
  bool allow_dual_steepest_edge_to_devex_switch;
  double dual_steepest_edge_weight_log_error_threshold;
  double costly_DSE_frequency;
  HighsInt num_costly_DSE_iteration;
  double costly_DSE_measure;

  double average_log_low_DSE_weight_error;
  double average_log_high_DSE_weight_error;
  // Needed globally??

  // Internal options - can't be changed externally
  bool run_quiet = false;
  bool store_squared_primal_infeasibility = false;

  //  bool analyse_lp_solution = true;
  // Options for reporting timing
  bool report_simplex_inner_clock = false;
  bool report_simplex_outer_clock = false;
  bool report_simplex_phases_clock = false;
  bool report_HFactor_clock = false;
  // Option for analysing the LP simplex iterations, INVERT time and rebuild
  // time
  bool analyse_lp = false;
  bool analyse_iterations = false;
  bool analyse_invert_form = false;
  //  bool analyse_invert_condition = false;
  //  bool analyse_invert_time = false;
  //  bool analyse_rebuild_time = false;

  // Simplex runtime information
  bool allow_cost_shifting = true;
  bool allow_cost_perturbation = true;
  bool allow_bound_perturbation = true;
  bool costs_shifted = false;
  bool costs_perturbed = false;
  bool bounds_shifted = false;
  bool bounds_perturbed = false;

  HighsInt num_primal_infeasibilities = -1;
  double max_primal_infeasibility;
  double sum_primal_infeasibilities;
  HighsInt num_dual_infeasibilities = -1;
  double max_dual_infeasibility;
  double sum_dual_infeasibilities;

  // Records of cumulative iteration counts - updated at the end of a phase
  HighsInt dual_phase1_iteration_count = 0;
  HighsInt dual_phase2_iteration_count = 0;
  HighsInt primal_phase1_iteration_count = 0;
  HighsInt primal_phase2_iteration_count = 0;
  HighsInt primal_bound_swap = 0;

  // Starting values for use in reportSimplexPhaseIterations
  HighsInt iteration_count0 = 0;
  HighsInt dual_phase1_iteration_count0 = 0;
  HighsInt dual_phase2_iteration_count0 = 0;
  HighsInt primal_phase1_iteration_count0 = 0;
  HighsInt primal_phase2_iteration_count0 = 0;
  HighsInt primal_bound_swap0 = 0;

  HighsInt min_concurrency = 1;
  HighsInt num_concurrency = 1;
  HighsInt max_concurrency = kSimplexConcurrencyLimit;

  // Info on PAMI iterations
  HighsInt multi_iteration = 0;

  // Number of UPDATE operations performed - should be zeroed when INVERT is
  // performed
  HighsInt update_count;
  // Value of dual objective - only set when computed from scratch in dual
  // rebuild()
  double dual_objective_value;
  // Value of primal objective - only set when computed from scratch in primal
  // rebuild()
  double primal_objective_value;

  // Value of dual objective that is updated in dual simplex solver
  double updated_dual_objective_value;
  // Value of primal objective that is updated in primal simplex solver
  double updated_primal_objective_value;
  // Number of logical variables in the basis
  HighsInt num_basic_logicals;
};

struct HighsSimplexBadBasisChangeRecord {
  bool taboo;
  HighsInt row_out;
  HighsInt variable_out;
  HighsInt variable_in;
  BadBasisChangeReason reason;
  double save_value = 0.0;
};

struct HighsRayRecord {
  HighsInt index;
  HighsInt sign;
  std::vector<double> value;
  HighsRayRecord getRayRecord() const;
  void setRayRecord(const HighsRayRecord& from_record);
  void clear();
};
#endif /* SIMPLEX_SIMPLEXSTRUCT_H_ */
