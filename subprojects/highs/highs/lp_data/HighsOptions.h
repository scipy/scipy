/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_OPTIONS_H_
#define LP_DATA_HIGHS_OPTIONS_H_

#include <cstring>  // For strlen
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"
#include "simplex/SimplexConst.h"
#include "util/HFactorConst.h"

using std::string;

enum class OptionStatus { kOk = 0, kUnknownOption, kIllegalValue };

const bool kAdvancedInDocumentation = false;

class OptionRecord {
 public:
  HighsOptionType type;
  std::string name;
  std::string description;
  bool advanced;

  OptionRecord(HighsOptionType Xtype, const std::string& Xname,
               const std::string& Xdescription, bool Xadvanced) {
    this->type = Xtype;
    this->name = Xname;
    this->description = Xdescription;
    this->advanced = Xadvanced;
  }

  virtual ~OptionRecord() {}
};

class OptionRecordBool : public OptionRecord {
 public:
  bool* value;
  bool default_value;
  OptionRecordBool(const std::string& Xname, const std::string& Xdescription,
                   bool Xadvanced, bool* Xvalue_pointer, bool Xdefault_value)
      : OptionRecord(HighsOptionType::kBool, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  void assignvalue(bool Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordBool() {}
};

class OptionRecordInt : public OptionRecord {
 public:
  HighsInt* value;
  HighsInt lower_bound;
  HighsInt default_value;
  HighsInt upper_bound;
  OptionRecordInt(const std::string& Xname, const std::string& Xdescription,
                  bool Xadvanced, HighsInt* Xvalue_pointer,
                  HighsInt Xlower_bound, HighsInt Xdefault_value,
                  HighsInt Xupper_bound)
      : OptionRecord(HighsOptionType::kInt, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    lower_bound = Xlower_bound;
    default_value = Xdefault_value;
    upper_bound = Xupper_bound;
    *value = default_value;
  }

  void assignvalue(HighsInt Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordInt() {}
};

class OptionRecordDouble : public OptionRecord {
 public:
  double* value;
  double lower_bound;
  double upper_bound;
  double default_value;
  OptionRecordDouble(const std::string& Xname, const std::string& Xdescription,
                     bool Xadvanced, double* Xvalue_pointer,
                     double Xlower_bound, double Xdefault_value,
                     double Xupper_bound)
      : OptionRecord(HighsOptionType::kDouble, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    lower_bound = Xlower_bound;
    default_value = Xdefault_value;
    upper_bound = Xupper_bound;
    *value = default_value;
  }

  void assignvalue(double Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordDouble() {}
};

class OptionRecordString : public OptionRecord {
 public:
  std::string* value;
  std::string default_value;
  OptionRecordString(const std::string& Xname, const std::string& Xdescription,
                     bool Xadvanced, std::string* Xvalue_pointer,
                     const std::string& Xdefault_value)
      : OptionRecord(HighsOptionType::kString, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  void assignvalue(const std::string& Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordString() {}
};

void highsOpenLogFile(HighsLogOptions& log_options,
                      std::vector<OptionRecord*>& option_records,
                      const std::string& log_file);

bool optionOffChooseOnOk(const HighsLogOptions& report_log_options,
                         const string& name, const string& value);
bool optionOffOnOk(const HighsLogOptions& report_log_options,
                   const string& name, const string& value);
bool optionSolverOk(const HighsLogOptions& report_log_options,
                    const string& value);
bool optionMipLpSolverOk(const HighsLogOptions& report_log_options,
                         const string& value);
bool optionMipIpmSolverOk(const HighsLogOptions& report_log_options,
                          const string& value);
bool optionHipoParallelTypeOk(const HighsLogOptions& report_log_options,
                              const string& value);
bool optionHipoSystemOk(const HighsLogOptions& report_log_options,
                        const string& value);
bool optionHipoOrderingOk(const HighsLogOptions& report_log_options,
                          const string& value);
bool optionHipoFactorOk(const HighsLogOptions& report_log_options,
                        const string& value);

bool boolFromString(std::string value, bool& bool_value);

OptionStatus getOptionIndex(const HighsLogOptions& report_log_options,
                            const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            HighsInt& index);

OptionStatus checkOptions(const HighsLogOptions& report_log_options,
                          const std::vector<OptionRecord*>& option_records);
OptionStatus checkOption(const HighsLogOptions& report_log_options,
                         const OptionRecordInt& option);
OptionStatus checkOption(const HighsLogOptions& report_log_options,
                         const OptionRecordDouble& option);

OptionStatus checkOptionValue(const HighsLogOptions& report_log_options,
                              OptionRecordInt& option_records,
                              const HighsInt value);
OptionStatus checkOptionValue(const HighsLogOptions& report_log_options,
                              OptionRecordDouble& option_records,
                              const double value);
OptionStatus checkOptionValue(const HighsLogOptions& report_log_options,
                              OptionRecordString& option_records,
                              const std::string& value);

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 std::vector<OptionRecord*>& option_records,
                                 const bool value);

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 std::vector<OptionRecord*>& option_records,
                                 const HighsInt value);
#ifdef HIGHSINT64
inline OptionStatus setLocalOptionValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    std::vector<OptionRecord*>& option_records, const int value) {
  return setLocalOptionValue(report_log_options, name, option_records,
                             HighsInt{value});
}
#endif
OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 std::vector<OptionRecord*>& option_records,
                                 const double value);
OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 HighsLogOptions& log_options,
                                 std::vector<OptionRecord*>& option_records,
                                 const std::string& value);
OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 HighsLogOptions& log_options,
                                 std::vector<OptionRecord*>& option_records,
                                 const char* value);

OptionStatus setLocalOptionValue(OptionRecordBool& option, const bool value);
OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 OptionRecordInt& option, const HighsInt value);
OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 OptionRecordDouble& option,
                                 const double value);
OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 OptionRecordString& option,
                                 const std::string& value);

void possibleLowerCaseOptionValue(const std::string& name, std::string& value);

OptionStatus passLocalOptions(const HighsLogOptions& report_log_options,
                              const HighsOptions& from_options,
                              HighsOptions& to_options);

OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<OptionRecord*>& option_records, bool* current_value,
    bool* default_value = nullptr);
OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<OptionRecord*>& option_records, HighsInt* current_value,
    HighsInt* min_value = nullptr, HighsInt* max_value = nullptr,
    HighsInt* default_value = nullptr);
OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<OptionRecord*>& option_records, double* current_value,
    double* min_value = nullptr, double* max_value = nullptr,
    double* default_value = nullptr);
OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<OptionRecord*>& option_records,
    std::string* current_value, std::string* default_value = nullptr);

OptionStatus getLocalOptionType(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<OptionRecord*>& option_records,
    HighsOptionType* type = nullptr);

void resetLocalOptions(std::vector<OptionRecord*>& option_records);

HighsStatus writeOptionsToFile(
    FILE* file, const HighsLogOptions& report_log_options,
    const std::vector<OptionRecord*>& option_records,
    const bool report_only_deviations = false,
    const HighsFileType file_type = HighsFileType::kFull);
void reportOptions(FILE* file, const HighsLogOptions& report_log_options,
                   const std::vector<OptionRecord*>& option_records,
                   const bool report_only_deviations = false,
                   const HighsFileType file_type = HighsFileType::kFull);
void reportOption(FILE* file, const HighsLogOptions& report_log_options,
                  const OptionRecordBool& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type);
void reportOption(FILE* file, const HighsLogOptions& report_log_options,
                  const OptionRecordInt& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type);
void reportOption(FILE* file, const HighsLogOptions& report_log_options,
                  const OptionRecordDouble& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type);
void reportOption(FILE* file, const HighsLogOptions& report_log_options,
                  const OptionRecordString& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type);

const string kHighsRunLogFile = "Highs.log";

const string kSimplexString = "simplex";
const string kIpmString = "ipm";
const string kHipoString = "hipo";
const string kIpxString = "ipx";
const string kPdlpString = "pdlp";
const string kQpAsmString = "qpasm";
const string kHiPdlpString = "hipdlp";

const HighsInt kKeepNRowsDeleteRows = -1;
const HighsInt kKeepNRowsDeleteEntries = 0;
const HighsInt kKeepNRowsKeepRows = 1;

// Strings for command line options
const string kModelFileString = "model_file";
const string kReadBasisFileString = "read_basis_file";
const string kWriteBasisFileString = "write_basis_file";
const string kPresolveString = "presolve";
const string kSolverString = "solver";
const string kParallelString = "parallel";
const string kThreadsString = "threads";
const string kRunCrossoverString = "run_crossover";
const string kTimeLimitString = "time_limit";
const string kOptionsFileString = "options_file";
const string kRandomSeedString = "random_seed";
const string kWriteSolutionFileString = "solution_file";
const string kRangingString = "ranging";
const string kVersionString = "version";
const string kWriteModelFileString = "write_model_file";
const string kWritePresolvedModelFileString = "write_presolved_model_file";
const string kWriteIisModelFileString = "write_iis_model_file";
const string kReadSolutionFileString = "read_solution_file";

// String for HiGHS log file option
const string kLogFileString = "log_file";

// String for other file names
#ifdef HIGHS_DEBUGSOL
const string kMipDebugSolutionFileString = "mip_debug_solution_file";
#endif
const string kMipImprovingSolutionFileString = "mip_improving_solution_file";

// Strings for HiPO system option
const string kHipoSystemString = "hipo_system";
const string kHipoAugmentedString = "augmented";
const string kHipoNormalEqString = "normaleq";

// Strings for MIP LP/IPM options
const string kMipLpSolverString = "mip_lp_solver";
const string kMipIpmSolverString = "mip_ipm_solver";
// Strings for HiPO parallel method
const string kHipoParallelString = "hipo_parallel_type";
const string kHipoTreeString = "tree";
const string kHipoNodeString = "node";
const string kHipoBothString = "both";

// Strings for HiPO matrix reordering
const string kHipoOrderingString = "hipo_ordering";
const string kHipoMetisString = "metis";
const string kHipoAmdString = "amd";
const string kHipoRcmString = "rcm";

const string kHipoFactorString = "hipo_factor";
const string kHipoFactorMultifrontal = "multifrontal";
const string kHipoFactorUplooking = "uplooking";

struct HighsOptionsStruct {
  // Run-time options read from the command line
  std::string presolve;
  std::string solver;
  std::string parallel;
  HighsInt threads;
  std::string run_crossover;
  double time_limit;
  std::string read_solution_file;
  std::string read_basis_file;
  std::string write_model_file;
  std::string solution_file;
  std::string write_basis_file;
  HighsInt random_seed;
  std::string ranging;

  // Options read from the file
  double infinite_cost;
  double infinite_bound;
  double small_matrix_value;
  double large_matrix_value;
  double kkt_tolerance;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  double primal_residual_tolerance;
  double dual_residual_tolerance;
  double optimality_tolerance;
  double objective_bound;
  double objective_target;
  HighsInt user_objective_scale;
  HighsInt user_bound_scale;
  HighsInt highs_debug_level;
  HighsInt highs_analysis_level;
  HighsInt simplex_strategy;
  HighsInt simplex_scale_strategy;
  HighsInt simplex_crash_strategy;
  HighsInt simplex_dual_edge_weight_strategy;
  HighsInt simplex_primal_edge_weight_strategy;
  HighsInt simplex_iteration_limit;
  HighsInt simplex_update_limit;
  HighsInt simplex_min_concurrency;
  HighsInt simplex_max_concurrency;

  std::string log_file;
  // Three bools are deprecated: remove in V2.0
  bool write_model_to_file;
  bool write_presolved_model_to_file;
  bool write_solution_to_file;

  HighsInt write_solution_style;
  HighsInt glpsol_cost_row_location;
  std::string write_presolved_model_file;
  std::string write_iis_model_file;

  // Control of HiGHS log
  bool output_flag;
  bool log_to_console;
  bool timeless_log;

  // Options for IPM solver
  double ipm_optimality_tolerance;
  HighsInt ipm_iteration_limit;
  std::string hipo_system;
  std::string hipo_parallel_type;
  std::string hipo_ordering;
  std::string hipo_factor;
  HighsInt hipo_block_size;
  HighsInt hipo_parallel_force;
  HighsInt hipo_parallel_forbid;
  bool hipo_serial_blas;

  // Options for PDLP solver
  HighsInt pdlp_features_off;
  HighsInt pdlp_iteration_limit;
  HighsInt pdlp_scaling_mode;
  HighsInt pdlp_ruiz_iterations;
  HighsInt pdlp_restart_strategy;
  HighsInt pdlp_cupdlpc_restart_method;
  HighsInt pdlp_step_size_strategy;
  double pdlp_optimality_tolerance;

  // Options for QP solver
  bool test_qp_oracle;
  bool qp_allow_hot_start;
  HighsInt qp_iteration_limit;
  HighsInt qp_nullspace_limit;
  double qp_regularization_value;

  // Options for IIS calculation
  HighsInt iis_strategy;
  double iis_time_limit;

  // Option for multi-objective optimization
  bool blend_multi_objectives;

  // Advanced options
  HighsInt log_dev_level;
  bool log_githash;
  bool solve_relaxation;
  bool allow_unbounded_or_infeasible;
  bool use_implied_bounds_from_presolve;
  bool lp_presolve_requires_basis_postsolve;
  bool mps_parser_type_free;
  bool use_warm_start;
  bool write_matrix_image;
  bool write_hessian_image;
  HighsInt keep_n_rows;
  HighsInt cost_scale_factor;
  HighsInt allowed_matrix_scale_factor;
  HighsInt allowed_cost_scale_factor;
  HighsInt ipx_dualize_strategy;
  HighsInt simplex_dualize_strategy;
  HighsInt simplex_permute_strategy;
  HighsInt max_dual_simplex_cleanup_level;
  HighsInt max_dual_simplex_phase1_cleanup_level;
  HighsInt simplex_price_strategy;
  HighsInt simplex_unscaled_solution_strategy;
  HighsInt presolve_mode;
  HighsInt presolve_reduction_limit;
  HighsInt restart_presolve_reduction_limit;
  HighsInt presolve_substitution_maxfillin;
  HighsInt presolve_rule_off;
  HighsInt presolve_rule_test;
  bool presolve_rule_logging;
  bool presolve_remove_slacks;
  bool presolve_fm_obj_reformulation;
  bool no_unnecessary_rebuild_refactor;
  double rebuild_refactor_solution_error_tolerance;
  double dual_steepest_edge_weight_error_tolerance;
  double dual_steepest_edge_weight_log_error_threshold;
  double dual_simplex_cost_perturbation_multiplier;
  double primal_simplex_bound_perturbation_multiplier;
  double dual_simplex_pivot_growth_tolerance;
  double presolve_pivot_threshold;
  double factor_pivot_threshold;
  double factor_pivot_tolerance;
  double start_crossover_tolerance;
  bool less_infeasible_DSE_check;
  bool less_infeasible_DSE_choose_row;
  bool use_original_HFactor_logic;
  //  bool allow_pdlp_cleanup;
  bool run_centring;
  HighsInt max_centring_steps;
  double centring_ratio_tolerance;

  // Options for iCrash
  bool icrash;
  bool icrash_dualize;
  std::string icrash_strategy;
  double icrash_starting_weight;
  HighsInt icrash_iterations;
  HighsInt icrash_approx_iter;
  bool icrash_exact;
  bool icrash_breakpoints;

  // Options for MIP solver
  bool mip_detect_symmetry;
  bool mip_allow_restart;
  HighsInt mip_max_nodes;
  HighsInt mip_max_stall_nodes;
  HighsInt mip_max_start_nodes;
  HighsInt mip_max_leaves;
  HighsInt mip_max_improving_sols;
  HighsInt mip_lp_age_limit;
  HighsInt mip_pool_age_limit;
  HighsInt mip_pool_soft_limit;
  HighsInt mip_pscost_minreliable;
  HighsInt mip_min_cliquetable_entries_for_parallelism;
  HighsInt mip_report_level;
  double mip_feasibility_tolerance;
  double mip_rel_gap;
  double mip_abs_gap;
  double mip_heuristic_effort;
  bool mip_heuristic_run_feasibility_jump;
  bool mip_heuristic_run_rins;
  bool mip_heuristic_run_rens;
  bool mip_heuristic_run_root_reduced_cost;
  bool mip_heuristic_run_zi_round;
  bool mip_heuristic_run_shifting;
  double mip_min_logging_interval;
  std::string mip_lp_solver;
  std::string mip_ipm_solver;
#ifdef HIGHS_DEBUGSOL
  std::string mip_debug_solution_file;
#endif
  bool mip_improving_solution_save;
  bool mip_improving_solution_report_sparse;
  std::string mip_improving_solution_file;
  bool mip_root_presolve_only;
  HighsInt mip_lifting_for_probing;
  bool mip_search_simulate_concurrency;
  bool mip_allow_cut_separation_at_nodes;

  // Logging callback identifiers
  HighsLogOptions log_options;
  virtual ~HighsOptionsStruct() {}

  HighsOptionsStruct()
      : presolve(""),
        solver(""),
        parallel(""),
        threads(0),
        run_crossover(""),
        time_limit(0.0),
        read_solution_file(""),
        read_basis_file(""),
        write_model_file(""),
        solution_file(""),
        write_basis_file(""),
        random_seed(0),
        ranging(""),
        infinite_cost(0.0),
        infinite_bound(0.0),
        small_matrix_value(0.0),
        large_matrix_value(0.0),
        kkt_tolerance(0.0),
        primal_feasibility_tolerance(0.0),
        dual_feasibility_tolerance(0.0),
        primal_residual_tolerance(0.0),
        dual_residual_tolerance(0.0),
        optimality_tolerance(0.0),
        objective_bound(0.0),
        objective_target(0.0),
        user_objective_scale(0),
        user_bound_scale(0),
        highs_debug_level(0),
        highs_analysis_level(0),
        simplex_strategy(0),
        simplex_scale_strategy(0),
        simplex_crash_strategy(0),
        simplex_dual_edge_weight_strategy(0),
        simplex_primal_edge_weight_strategy(0),
        simplex_iteration_limit(0),
        simplex_update_limit(0),
        simplex_min_concurrency(0),
        simplex_max_concurrency(0),
        log_file(""),
        write_model_to_file(false),
        write_presolved_model_to_file(false),
        write_solution_to_file(false),
        write_solution_style(0),
        glpsol_cost_row_location(0),
        write_presolved_model_file(""),
        write_iis_model_file(""),
        output_flag(false),
        log_to_console(false),
        timeless_log(false),
        ipm_optimality_tolerance(0.0),
        ipm_iteration_limit(0),
        hipo_system(""),
        hipo_parallel_type(""),
        hipo_ordering(""),
        hipo_factor(""),
        hipo_block_size(0),
        hipo_parallel_force(0),
        hipo_parallel_forbid(0),
        hipo_serial_blas(true),
        pdlp_features_off(0),
        pdlp_iteration_limit(0),
        pdlp_scaling_mode(0),
        pdlp_ruiz_iterations(0),
        pdlp_restart_strategy(0),
        pdlp_cupdlpc_restart_method(0),
        pdlp_step_size_strategy(0),
        pdlp_optimality_tolerance(0.0),
        test_qp_oracle(false),
        qp_allow_hot_start(false),
        qp_iteration_limit(0),
        qp_nullspace_limit(0),
        qp_regularization_value(0),
        iis_strategy(0),
        iis_time_limit(0.0),
        blend_multi_objectives(false),
        log_dev_level(0),
        log_githash(false),
        solve_relaxation(false),
        allow_unbounded_or_infeasible(false),
        use_implied_bounds_from_presolve(false),
        lp_presolve_requires_basis_postsolve(false),
        mps_parser_type_free(false),
        use_warm_start(true),
        write_matrix_image(false),
        write_hessian_image(false),
        keep_n_rows(0),
        cost_scale_factor(0),
        allowed_matrix_scale_factor(0),
        allowed_cost_scale_factor(0),
        ipx_dualize_strategy(0),
        simplex_dualize_strategy(0),
        simplex_permute_strategy(0),
        max_dual_simplex_cleanup_level(0),
        max_dual_simplex_phase1_cleanup_level(0),
        simplex_price_strategy(0),
        simplex_unscaled_solution_strategy(0),
        presolve_mode(0),
        presolve_reduction_limit(0),
        restart_presolve_reduction_limit(0),
        presolve_substitution_maxfillin(0),
        presolve_rule_off(0),
        presolve_rule_test(0),
        presolve_rule_logging(false),
        presolve_remove_slacks(false),
        presolve_fm_obj_reformulation(false),
        no_unnecessary_rebuild_refactor(false),
        rebuild_refactor_solution_error_tolerance(0.0),
        dual_steepest_edge_weight_error_tolerance(0.0),
        dual_steepest_edge_weight_log_error_threshold(0.0),
        dual_simplex_cost_perturbation_multiplier(0.0),
        primal_simplex_bound_perturbation_multiplier(0.0),
        dual_simplex_pivot_growth_tolerance(0.0),
        presolve_pivot_threshold(0.0),
        factor_pivot_threshold(0.0),
        factor_pivot_tolerance(0.0),
        start_crossover_tolerance(0.0),
        less_infeasible_DSE_check(false),
        less_infeasible_DSE_choose_row(false),
        use_original_HFactor_logic(false),
        //        allow_pdlp_cleanup(false),
        run_centring(false),
        max_centring_steps(0),
        centring_ratio_tolerance(0.0),
        icrash(false),
        icrash_dualize(false),
        icrash_strategy(""),
        icrash_starting_weight(0.0),
        icrash_iterations(0),
        icrash_approx_iter(0),
        icrash_exact(false),
        icrash_breakpoints(false),
        mip_detect_symmetry(false),
        mip_allow_restart(false),
        mip_max_nodes(0),
        mip_max_stall_nodes(0),
        mip_max_start_nodes(0),
        mip_max_leaves(0),
        mip_max_improving_sols(0),
        mip_lp_age_limit(0),
        mip_pool_age_limit(0),
        mip_pool_soft_limit(0),
        mip_pscost_minreliable(0),
        mip_min_cliquetable_entries_for_parallelism(0),
        mip_report_level(0),
        mip_feasibility_tolerance(0.0),
        mip_rel_gap(0.0),
        mip_abs_gap(0.0),
        mip_heuristic_effort(0.0),
        mip_heuristic_run_feasibility_jump(false),
        mip_heuristic_run_rins(false),
        mip_heuristic_run_rens(false),
        mip_heuristic_run_root_reduced_cost(false),
        mip_heuristic_run_zi_round(false),
        mip_heuristic_run_shifting(false),
        mip_min_logging_interval(0.0),
        mip_lp_solver(""),
        mip_ipm_solver(""),
#ifdef HIGHS_DEBUGSOL
        mip_debug_solution_file(""),
#endif
        mip_improving_solution_save(false),
        mip_improving_solution_report_sparse(false),
        mip_improving_solution_file(""),
        mip_root_presolve_only(false),
        mip_lifting_for_probing(-1),
        mip_search_simulate_concurrency(false),
        // clang-format off
        mip_allow_cut_separation_at_nodes(true) {};
  // clang-format on
};

// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily. The options below
// are just what has been used to parse options from argv.
// todo: when creating the new options don't forget underscores for class
// variables but no underscores for struct
class HighsOptions : public HighsOptionsStruct {
 public:
  HighsOptions() {
    initRecords();
    setLogOptions();
  }

  HighsOptions(const HighsOptions& options) {
    initRecords();
    HighsOptionsStruct::operator=(options);
    setLogOptions();
  }

  const HighsOptions& operator=(const HighsOptions& other) {
    if (&other != this) {
      if ((HighsInt)records.size() == 0) initRecords();
      HighsOptionsStruct::operator=(other);
      this->log_options.log_stream = other.log_options.log_stream;
      setLogOptions();
    }
    return *this;
  }

  // there is deliberately no move constructor or move assignment operator.
  // each option record holds a pointer to the member of the object that
  // created it, so records cannot be taken over from another object: a move
  // constructor that steals them leaves the new object's records pointing at
  // the members of the moved-from object (and dangling once it is destroyed).
  // without move operations, rvalues bind to the copy constructor and copy
  // assignment operator above, which give each object its own records. do
  // not declare the move operations as deleted, since this would make moving
  // (e.g. via std::move or returning by value) fail to compile.

  virtual ~HighsOptions() {
    if (records.size() > 0) deleteRecords();
  }

 private:
  void initRecords() {
    OptionRecordBool* record_bool;
    OptionRecordInt* record_int;
    OptionRecordDouble* record_double;
    OptionRecordString* record_string;
    bool advanced = false;
    const bool now_advanced = true;
    // Options read from the command line
    record_string = new OptionRecordString(
        kPresolveString, "Presolve: \"off\", \"choose\" or \"on\"", advanced,
        &presolve, kHighsChooseString);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kSolverString,
        "LP/QP solver: \"choose\", \"simplex\", "
        "\"ipm\", \"ipx\", \"hipo\", \"pdlp\", \"qpasm\" or \"hipdlp\", ",
        advanced, &solver, kHighsChooseString);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kParallelString, "Parallel: \"off\", \"choose\" or \"on\"", advanced,
        &parallel, kHighsChooseString);
    records.push_back(record_string);

    record_int = new OptionRecordInt(
        kThreadsString,
        "Maximum number of threads used by HiGHS (0: automatic)", advanced,
        &threads, 0, 0, kHighsIInf);
    records.push_back(record_int);

    record_string = new OptionRecordString(
        kRunCrossoverString, "Run IPM crossover: \"off\", \"choose\" or \"on\"",
        advanced, &run_crossover, kHighsOnString);
    records.push_back(record_string);

    record_double =
        new OptionRecordDouble(kTimeLimitString, "Time limit (seconds)",
                               advanced, &time_limit, 0, kHighsInf, kHighsInf);
    records.push_back(record_double);

    record_string =
        new OptionRecordString(kRangingString,
                               "Compute cost, bound, RHS and basic solution "
                               "ranging: \"off\" or \"on\"",
                               advanced, &ranging, kHighsOffString);
    records.push_back(record_string);
    //
    // Options read from the file
    record_double = new OptionRecordDouble(
        "infinite_cost",
        "Limit on |cost coefficient|: values greater than or equal to "
        "this will be treated as infinite",
        advanced, &infinite_cost, 1e15, 1e20, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "infinite_bound",
        "Limit on |variable/constraint bound|: values greater than or equal to "
        "this will be treated as infinite",
        advanced, &infinite_bound, 1e15, 1e20, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "small_matrix_value",
        "Lower limit on |matrix entries|: values less than or equal to this "
        "will be "
        "treated as zero",
        advanced, &small_matrix_value, 1e-12, 1e-9, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "large_matrix_value",
        "Upper limit on |matrix entries|: values greater than or equal to "
        "this will be treated as infinite",
        advanced, &large_matrix_value, 1e0, 1e15, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "kkt_tolerance",
        "If changed from its default value, this tolerance is used for all "
        "feasibility and optimality (KKT) measures",
        advanced, &kkt_tolerance, kMinimumKktTolerance, kDefaultKktTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "primal_feasibility_tolerance", "Primal feasibility tolerance",
        advanced, &primal_feasibility_tolerance, kMinimumKktTolerance,
        kDefaultKktTolerance, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_feasibility_tolerance", "Dual feasibility tolerance", advanced,
        &dual_feasibility_tolerance, kMinimumKktTolerance, kDefaultKktTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "primal_residual_tolerance", "Primal residual tolerance", advanced,
        &primal_residual_tolerance, kMinimumKktTolerance, kDefaultKktTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_residual_tolerance", "Dual residual tolerance", advanced,
        &dual_residual_tolerance, kMinimumKktTolerance, kDefaultKktTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "optimality_tolerance", "Optimality tolerance", advanced,
        &optimality_tolerance, kMinimumKktTolerance, kDefaultKktTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "objective_bound",
        "Objective bound for termination of the dual simplex solver", advanced,
        &objective_bound, -kHighsInf, kHighsInf, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "objective_target",
        "Objective target for termination of the MIP solver", advanced,
        //"primal simplex and "
        &objective_target, -kHighsInf, -kHighsInf, kHighsInf);
    records.push_back(record_double);

    record_int =
        new OptionRecordInt(kRandomSeedString, "Random seed used in HiGHS",
                            advanced, &random_seed, 0, 0, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "user_objective_scale",
        "Exponent of power-of-two objective scaling for model", advanced,
        &user_objective_scale, -kHighsIInf, 0, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "user_bound_scale", "Exponent of power-of-two bound scaling for model",
        advanced, &user_bound_scale, -kHighsIInf, 0, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt("highs_debug_level",
                                     "Debugging level in HiGHS", now_advanced,
                                     &highs_debug_level, kHighsDebugLevelMin,
                                     kHighsDebugLevelMin, kHighsDebugLevelMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "highs_analysis_level", "Analysis level in HiGHS", now_advanced,
        &highs_analysis_level, kHighsAnalysisLevelMin,
        kHighsAnalysisLevelMin,  // kHighsAnalysisLevelMipTime,  //
        kHighsAnalysisLevelMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_strategy",
        "Strategy for simplex solver 0 => Choose; 1 => Dual (serial); 2 => "
        "Dual (SIP); 3 => Dual (PAMI); 4 => Primal",
        advanced, &simplex_strategy, kSimplexStrategyMin, kSimplexStrategyDual,
        kSimplexStrategyMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_scale_strategy",
        "Simplex scaling strategy: off / choose / equilibration (default) / "
        "forced "
        "equilibration / max value (0/1/2/3/4)",
        advanced, &simplex_scale_strategy, kSimplexScaleStrategyMin,
        kSimplexScaleStrategyEquilibration, kSimplexScaleStrategyMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_crash_strategy",
        "Strategy for simplex crash: off / LTSSF / Bixby (0/1/2)", now_advanced,
        &simplex_crash_strategy, kSimplexCrashStrategyMin,
        kSimplexCrashStrategyOff, kSimplexCrashStrategyMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_dual_edge_weight_strategy",
        "Strategy for simplex dual edge weights: Choose / "
        "Dantzig / Devex / Steepest "
        "Edge (-1/0/1/2)",
        advanced, &simplex_dual_edge_weight_strategy,
        kSimplexEdgeWeightStrategyMin, kSimplexEdgeWeightStrategyChoose,
        kSimplexEdgeWeightStrategyMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_primal_edge_weight_strategy",
        "Strategy for simplex primal edge weights: Choose "
        "/ Dantzig / Devex / Steepest "
        "Edge (-1/0/1/2)",
        advanced, &simplex_primal_edge_weight_strategy,
        kSimplexEdgeWeightStrategyMin, kSimplexEdgeWeightStrategyChoose,
        kSimplexEdgeWeightStrategyMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_iteration_limit",
        "Iteration limit for simplex solver when solving LPs, but not "
        "subproblems in the MIP solver",
        advanced, &simplex_iteration_limit, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_update_limit",
        "Limit on the number of simplex UPDATE operations", advanced,
        &simplex_update_limit, 0, 5000, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_min_concurrency",
        "Minimum level of concurrency in parallel simplex", now_advanced,
        &simplex_min_concurrency, 1, 1, kSimplexConcurrencyLimit);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("simplex_max_concurrency",
                            "Maximum level of concurrency in parallel simplex",
                            advanced, &simplex_max_concurrency, 1,
                            kSimplexConcurrencyLimit, kSimplexConcurrencyLimit);
    records.push_back(record_int);

    record_bool =
        new OptionRecordBool("output_flag", "Enables or disables solver output",
                             advanced, &output_flag, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("log_to_console",
                                       "Enables or disables console logging",
                                       advanced, &log_to_console, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "timeless_log", "Suppression of time-based data in logging", true,
        &timeless_log, false);
    records.push_back(record_bool);

    record_string = new OptionRecordString(kLogFileString, "Log file", advanced,
                                           &log_file, "");
    records.push_back(record_string);

    record_bool =
        new OptionRecordBool("write_model_to_file", "Write the model to a file",
                             advanced, &write_model_to_file, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "write_presolved_model_to_file", "Write the presolved model to a file",
        advanced, &write_presolved_model_to_file, false);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("write_solution_to_file",
                             "Write the primal and dual solution to a file",
                             advanced, &write_solution_to_file, false);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "write_solution_style",
        "Style of solution file (raw = computer-readable, "
        "pretty = human-readable): "
        "-1 => HiGHS old raw (deprecated); 0 => HiGHS raw; "
        "1 => HiGHS pretty; 2 => Glpsol raw; 3 => Glpsol pretty; "
        "4 => HiGHS sparse raw; 5 => HiGHS partial raw",
        advanced, &write_solution_style, kSolutionStyleMin, kSolutionStyleRaw,
        kSolutionStyleMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "glpsol_cost_row_location",
        "Location of cost row for Glpsol file: "
        "-2 => Last; -1 => None; 0 => None if empty, otherwise data file "
        "location; 1 <= n <= num_row => Location n; n > "
        "num_row => Last",
        advanced, &glpsol_cost_row_location, kGlpsolCostRowLocationMin, 0,
        kHighsIInf);
    records.push_back(record_int);

    record_bool = new OptionRecordBool("icrash", "Run iCrash", now_advanced,
                                       &icrash, false);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("icrash_dualize", "Dualize strategy for iCrash",
                             now_advanced, &icrash_dualize, false);
    records.push_back(record_bool);

    record_string =
        new OptionRecordString("icrash_strategy", "Strategy for iCrash",
                               now_advanced, &icrash_strategy, "ICA");
    records.push_back(record_string);

    record_double = new OptionRecordDouble(
        "icrash_starting_weight", "iCrash starting weight", now_advanced,
        &icrash_starting_weight, 1e-10, 1e-3, 1e50);
    records.push_back(record_double);

    record_int =
        new OptionRecordInt("icrash_iterations", "iCrash iterations",
                            now_advanced, &icrash_iterations, 0, 30, 200);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "icrash_approx_iter", "iCrash approximate minimization iterations",
        now_advanced, &icrash_approx_iter, 0, 50, 100);
    records.push_back(record_int);

    record_bool = new OptionRecordBool("icrash_exact",
                                       "Exact subproblem solution for iCrash",
                                       now_advanced, &icrash_exact, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "icrash_breakpoints", "Exact subproblem solution for iCrash",
        now_advanced, &icrash_breakpoints, false);
    records.push_back(record_bool);

    record_string = new OptionRecordString(
        kReadSolutionFileString, "Read solution file", advanced,
        &read_solution_file, kHighsFilenameDefault);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kReadBasisFileString, "Read basis file", advanced, &read_basis_file,
        kHighsFilenameDefault);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kWriteModelFileString, "Write model file", advanced, &write_model_file,
        kHighsFilenameDefault);
    records.push_back(record_string);

    record_string =
        new OptionRecordString(kWriteSolutionFileString, "Write solution file",
                               advanced, &solution_file, kHighsFilenameDefault);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kWriteBasisFileString, "Write basis file", advanced, &write_basis_file,
        kHighsFilenameDefault);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kWritePresolvedModelFileString, "Write presolved model file", advanced,
        &write_presolved_model_file, kHighsFilenameDefault);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kWriteIisModelFileString, "Write IIS model file", advanced,
        &write_iis_model_file, kHighsFilenameDefault);
    records.push_back(record_string);

    record_int = new OptionRecordInt(
        "presolve_mode",
        "Presolve mode 0 => default; 1 => initial sweep; 2 => light; 3 => full",
        advanced, &presolve_mode, kPresolveModeMin, kPresolveModeDefault,
        kPresolveModeMax);
    records.push_back(record_int);

    record_bool = new OptionRecordBool(
        "mip_detect_symmetry", "Whether MIP symmetry should be detected",
        advanced, &mip_detect_symmetry, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("mip_allow_restart",
                                       "Whether MIP restart is permitted",
                                       advanced, &mip_allow_restart, true);
    records.push_back(record_bool);

    record_int = new OptionRecordInt("mip_max_nodes",
                                     "MIP solver max number of nodes", advanced,
                                     &mip_max_nodes, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_max_stall_nodes",
        "MIP solver max number of nodes where estimate is above cutoff bound",
        advanced, &mip_max_stall_nodes, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_max_start_nodes",
        "MIP solver max number of nodes when completing a partial MIP start",
        advanced, &mip_max_start_nodes, 0, 500, kHighsIInf);
    records.push_back(record_int);

#ifdef HIGHS_DEBUGSOL
    record_string = new OptionRecordString(
        kMipDebugSolutionFileString,
        "Solution file for debug solution of the MIP solver", advanced,
        &mip_debug_solution_file, kHighsFilenameDefault);
    records.push_back(record_string);
#endif

    record_bool =
        new OptionRecordBool("mip_improving_solution_save",
                             "Whether improving MIP solutions should be saved",
                             advanced, &mip_improving_solution_save, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "mip_improving_solution_report_sparse",
        "Whether improving MIP solutions should be reported in sparse format",
        advanced, &mip_improving_solution_report_sparse, false);
    records.push_back(record_bool);

    record_string = new OptionRecordString(
        kMipImprovingSolutionFileString,
        "File for reporting improving MIP solutions: not reported for an empty "
        "string \\\"\\\"",
        advanced, &mip_improving_solution_file, kHighsFilenameDefault);
    records.push_back(record_string);

    record_bool = new OptionRecordBool(
        "mip_root_presolve_only",
        "Whether MIP presolve is only applied at the root node", advanced,
        &mip_root_presolve_only, false);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "mip_lifting_for_probing", "Level of lifting for probing that is used",
        advanced, &mip_lifting_for_probing, -1, -1, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_max_leaves", "MIP solver max number of leaf nodes", advanced,
        &mip_max_leaves, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_max_improving_sols",
        "Limit on the number of improving solutions found to stop the MIP "
        "solver prematurely",
        advanced, &mip_max_improving_sols, 1, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_lp_age_limit",
        "Maximal age of dynamic LP rows before "
        "they are removed from the LP relaxation in the MIP solver",
        advanced, &mip_lp_age_limit, 0, 10,
        (std::numeric_limits<int16_t>::max)());
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_pool_age_limit",
        "Maximal age of rows in the MIP solver cutpool before they are deleted",
        advanced, &mip_pool_age_limit, 0, 30, 1000);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_pool_soft_limit",
        "Soft limit on the number of rows in the "
        "MIP solver cutpool for dynamic age adjustment",
        advanced, &mip_pool_soft_limit, 1, 10000, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_pscost_minreliable",
        "Minimal number of observations before "
        "MIP solver pseudo costs are considered reliable",
        advanced, &mip_pscost_minreliable, 0, 8, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "mip_min_cliquetable_entries_for_parallelism",
        "Minimal number of entries in the MIP solver cliquetable before "
        "neighbourhood "
        "queries of the conflict graph use parallel processing",
        advanced, &mip_min_cliquetable_entries_for_parallelism, 0, 100000,
        kHighsIInf);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("mip_report_level", "MIP solver reporting level",
                            now_advanced, &mip_report_level, 0, 1, 2);
    records.push_back(record_int);

    record_double = new OptionRecordDouble(
        "mip_feasibility_tolerance", "MIP integrality tolerance", advanced,
        &mip_feasibility_tolerance, kMinimumMipTolerance, kDefaultMipTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "mip_heuristic_effort", "Effort spent for MIP heuristics", advanced,
        &mip_heuristic_effort, 0.0, 0.05, 1.0);
    records.push_back(record_double);

    record_bool =
        new OptionRecordBool("mip_heuristic_run_feasibility_jump",
                             "Use the feasibility jump heuristic", advanced,
                             &mip_heuristic_run_feasibility_jump, true);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("mip_heuristic_run_rins", "Use the RINS heuristic",
                             advanced, &mip_heuristic_run_rins, true);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("mip_heuristic_run_rens", "Use the RENS heuristic",
                             advanced, &mip_heuristic_run_rens, true);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("mip_heuristic_run_root_reduced_cost",
                             "Use the rootReducedCost heuristic", advanced,
                             &mip_heuristic_run_root_reduced_cost, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("mip_heuristic_run_zi_round",
                                       "Use the ZI Round heuristic", advanced,
                                       &mip_heuristic_run_zi_round, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("mip_heuristic_run_shifting",
                                       "Use the Shifting heuristic", advanced,
                                       &mip_heuristic_run_shifting, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "mip_allow_cut_separation_at_nodes",
        "Whether cut separation at nodes other than the root node is permitted",
        advanced, &mip_allow_cut_separation_at_nodes, true);
    records.push_back(record_bool);

    record_double = new OptionRecordDouble(
        "mip_rel_gap",
        "Tolerance on relative gap, |ub-lb|/|ub|, to determine whether "
        "optimality has been reached for a MIP instance",
        advanced, &mip_rel_gap, 0.0, 1e-4, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "mip_abs_gap",
        "Tolerance on absolute gap of MIP, |ub-lb|, to determine whether "
        "optimality has been reached for a MIP instance",
        advanced, &mip_abs_gap, 0.0, 1e-6, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "mip_min_logging_interval", "MIP minimum logging interval", advanced,
        &mip_min_logging_interval, 0, 5, kHighsInf);
    records.push_back(record_double);

    record_string =
        new OptionRecordString(kMipLpSolverString,
                               "MIP LP solver: \"choose\", \"simplex\", "
                               "\"ipm\", \"ipx\" or \"hipo\"",
                               advanced, &mip_lp_solver, kHighsChooseString);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kMipIpmSolverString, "MIP IPM solver: \"choose\", \"ipx\" or \"hipo\"",
        advanced, &mip_ipm_solver, kHighsChooseString);
    records.push_back(record_string);

    record_double = new OptionRecordDouble(
        "ipm_optimality_tolerance", "IPM optimality tolerance", advanced,
        &ipm_optimality_tolerance, kMinimumIpmTolerance, kDefaultIpmTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_bool = new OptionRecordBool(
        "mip_search_simulate_concurrency",
        "Simulate MIP search concurrency on a single thread", advanced,
        &mip_search_simulate_concurrency, false);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "ipm_iteration_limit", "Iteration limit for IPM solver", advanced,
        &ipm_iteration_limit, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_string = new OptionRecordString(
        kHipoSystemString,
        "HiPO Newton system: \"choose\", \"augmented\" or \"normaleq\"",
        advanced, &hipo_system, kHighsChooseString);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kHipoParallelString,
        "HiPO parallelism: \"tree\", \"node\", \"both\" or \"choose\"",
        advanced, &hipo_parallel_type, kHighsChooseString);
    records.push_back(record_string);

    record_int = new OptionRecordInt(
        "hipo_parallel_force", "Bit mask to force parallel techniques in HiPO",
        advanced, &hipo_parallel_force, 0, 0,
        static_cast<int>(hipo::ParallelTechnique::kMaxSum));
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("hipo_parallel_forbid",
                            "Bit mask to forbid parallel techniques in HiPO",
                            advanced, &hipo_parallel_forbid, 0, 0,
                            static_cast<int>(hipo::ParallelTechnique::kMaxSum));
    records.push_back(record_int);

    record_string = new OptionRecordString(
        kHipoOrderingString,
        "HiPO matrix reordering: \"choose\", \"metis\", \"amd\" or \"rcm\"",
        advanced, &hipo_ordering, kHighsChooseString);
    records.push_back(record_string);

    record_string = new OptionRecordString(
        kHipoFactorString,
        "HiPO matrix factorisation: \"choose\", \"multifrontal\" "
        "or \"uplooking\"",
        advanced, &hipo_factor, kHighsChooseString);
    records.push_back(record_string);

    record_int = new OptionRecordInt(
        "hipo_block_size", "Block size for dense linear algebra within HiPO",
        advanced, &hipo_block_size, 0, 128, kHighsIInf);
    records.push_back(record_int);

    record_bool = new OptionRecordBool(
        "hipo_serial_blas",
        "Force OpenBLAS to run in serial when used by HiPO. This setting is "
        "process-wide, so affects any other use of OpenBLAS",
        advanced, &hipo_serial_blas, true);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "pdlp_iteration_limit", "Iteration limit for PDLP solver", advanced,
        &pdlp_iteration_limit, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("pdlp_scaling_mode",
                            "Scaling mode for PDLP solver (default = "
                            "5): 1 => Ruiz; 2 => L2; 4 => PC",
                            advanced, &pdlp_scaling_mode, kPdlpScalingMin,
                            kPdlpScalingRuiz + kPdlpScalingPC, kPdlpScalingMax);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("pdlp_ruiz_iterations",
                            "Number of Ruiz scaling iteraitons for PDLP solver",
                            advanced, &pdlp_ruiz_iterations, 0, 10, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "pdlp_restart_strategy",
        "Restart strategy for PDLP solver: 0 => off; "
        "1 => fixed; 2 => adaptive; 3 => Halpern",
        advanced, &pdlp_restart_strategy, kPdlpRestartStrategyMin,
        kPdlpRestartStrategyAdaptive, kPdlpRestartStrategyMax);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("pdlp_cupdlpc_restart_method",
                            "Restart mode for cuPDLP-C solver: 0 => none; "
                            "1 => GPU (default); 2 => CPU ",
                            advanced, &pdlp_cupdlpc_restart_method, 0, 1, 2);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "pdlp_step_size_strategy",
        "Step size strategy for PDLP solver: 0 => fixed; "
        "1 => adaptive; 2 => Malitsky-Pock; 3 => PID",
        advanced, &pdlp_step_size_strategy, kPdlpStepSizeStrategyMin,
        kPdlpStepSizeStrategyAdaptive, kPdlpStepSizeStrategyMax);
    records.push_back(record_int);

    record_double = new OptionRecordDouble(
        "pdlp_optimality_tolerance", "PDLP optimality tolerance", advanced,
        &pdlp_optimality_tolerance, kMinimumKktTolerance, kDefaultKktTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_bool = new OptionRecordBool("test_qp_oracle",
                                       "Use an oracle for the QP Hessian",
                                       advanced, &test_qp_oracle, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "qp_allow_hot_start", "Allow the active set QP solver to hot start",
        advanced, &qp_allow_hot_start, false);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "qp_iteration_limit", "Iteration limit for the active set QP solver",
        advanced, &qp_iteration_limit, 0, kHighsIInf, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "qp_nullspace_limit", "Nullspace limit for the active set QP solver",
        advanced, &qp_nullspace_limit, 0, 4000, kHighsIInf);
    records.push_back(record_int);

    record_double = new OptionRecordDouble(
        "qp_regularization_value",
        "Regularization value added to the Hessian in the active set QP solver",
        advanced, &qp_regularization_value, 0, kHessianRegularizationValue,
        kHighsInf);
    records.push_back(record_double);

    record_int = new OptionRecordInt(
        "iis_strategy",
        "Strategy for IIS calculation: "
        "0 => Light test; 1 => Try dual ray; "
        "2 => Try elastic LP; 4 => Find true IIS; "
        "8 => Prioritise columns; 16 => Find relaxation IIS for MIP",
        advanced, &iis_strategy, kIisStrategyMin, kIisStrategyDefault,
        kIisStrategyMax);
    records.push_back(record_int);

    record_double = new OptionRecordDouble(
        "iis_time_limit", "Time limit for computing IIS (seconds)", advanced,
        &iis_time_limit, 0, kHighsInf, kHighsInf);
    records.push_back(record_double);

    record_bool = new OptionRecordBool(
        "blend_multi_objectives",
        "Blend multiple objectives or apply lexicographically", advanced,
        &blend_multi_objectives, true);
    records.push_back(record_bool);

    // Fix the number of user settable options
    num_user_settable_options_ = static_cast<HighsInt>(records.size());

    // Advanced options
    advanced = true;

    record_int =
        new OptionRecordInt("log_dev_level",
                            "Output development messages: 0 => none; 1 => "
                            "info; 2 => detailed; 3 => verbose",
                            advanced, &log_dev_level, kHighsLogDevLevelMin,
                            kHighsLogDevLevelNone, kHighsLogDevLevelMax);
    records.push_back(record_int);

    record_bool = new OptionRecordBool("log_githash", "Log the githash",
                                       advanced, &log_githash, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "solve_relaxation", "Solve the relaxation of discrete model components",
        advanced, &solve_relaxation, false);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("allow_unbounded_or_infeasible",
                             "Allow ModelStatus::kUnboundedOrInfeasible",
                             advanced, &allow_unbounded_or_infeasible, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "use_implied_bounds_from_presolve",
        "Use relaxed implied bounds from presolve: redundant option!", advanced,
        &use_implied_bounds_from_presolve, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "lp_presolve_requires_basis_postsolve",
        "Prevents LP presolve steps for which postsolve cannot maintain a "
        "basis: redundant option!",
        advanced, &lp_presolve_requires_basis_postsolve, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("mps_parser_type_free",
                                       "Use the free format MPS file reader",
                                       advanced, &mps_parser_type_free, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("use_warm_start",
                                       "Use any warm start that is available",
                                       advanced, &use_warm_start, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "write_matrix_image",
        "Write an image of the constraint matrix to a file", advanced,
        &write_matrix_image, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "write_hessian_image", "Write an image of the Hessian to a file",
        advanced, &write_hessian_image, false);
    records.push_back(record_bool);

    record_int =
        new OptionRecordInt("keep_n_rows",
                            "For multiple N-rows in MPS files: delete rows / "
                            "delete entries / keep rows (-1/0/1)",
                            advanced, &keep_n_rows, kKeepNRowsDeleteRows,
                            kKeepNRowsDeleteRows, kKeepNRowsKeepRows);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("cost_scale_factor", "Scaling factor for costs",
                            advanced, &cost_scale_factor, -20, 0, 20);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "allowed_matrix_scale_factor",
        "Largest power-of-two factor permitted when "
        "scaling the constraint matrix",
        advanced, &allowed_matrix_scale_factor, 0,
        kDefaultAllowedMatrixPow2Scale, kMaxAllowedMatrixPow2Scale);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "allowed_cost_scale_factor",
        "Largest power-of-two factor permitted when scaling the costs",
        advanced, &allowed_cost_scale_factor, 0, 0, 20);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "ipx_dualize_strategy", "Strategy for dualizing before IPX", advanced,
        &ipx_dualize_strategy, kIpxDualizeStrategyMin, kIpxDualizeStrategyLukas,
        kIpxDualizeStrategyMax);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_dualize_strategy", "Strategy for dualizing before simplex",
        advanced, &simplex_dualize_strategy, kHighsOptionOff, kHighsOptionOff,
        kHighsOptionOn);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_permute_strategy", "Strategy for permuting before simplex",
        advanced, &simplex_permute_strategy, kHighsOptionOff, kHighsOptionOff,
        kHighsOptionOn);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "max_dual_simplex_cleanup_level", "Max level of dual simplex cleanup",
        advanced, &max_dual_simplex_cleanup_level, 0, 1, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "max_dual_simplex_phase1_cleanup_level",
        "Max level of dual simplex phase 1 cleanup", advanced,
        &max_dual_simplex_phase1_cleanup_level, 0, 2, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_price_strategy", "Strategy for PRICE in simplex", advanced,
        &simplex_price_strategy, kSimplexPriceStrategyMin,
        kSimplexPriceStrategyRowSwitchColSwitch, kSimplexPriceStrategyMax);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("simplex_unscaled_solution_strategy",
                            "Strategy for solving unscaled LP in simplex",
                            advanced, &simplex_unscaled_solution_strategy,
                            kSimplexUnscaledSolutionStrategyMin,
                            kSimplexUnscaledSolutionStrategyRefine,
                            kSimplexUnscaledSolutionStrategyMax);
    records.push_back(record_int);

    record_bool = new OptionRecordBool(
        "no_unnecessary_rebuild_refactor",
        "No unnecessary refactorization on simplex rebuild", advanced,
        &no_unnecessary_rebuild_refactor, true);
    records.push_back(record_bool);

    record_double = new OptionRecordDouble(
        "rebuild_refactor_solution_error_tolerance",
        "Tolerance on solution error when considering refactorization on "
        "simplex rebuild",
        advanced, &rebuild_refactor_solution_error_tolerance, -kHighsInf, 1e-8,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_steepest_edge_weight_error_tolerance",
        "Tolerance on dual steepest edge weight errors", advanced,
        &dual_steepest_edge_weight_error_tolerance, 0.0, kHighsInf, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_steepest_edge_weight_log_error_threshold",
        "Threshold on dual steepest edge weight errors for Devex switch",
        advanced, &dual_steepest_edge_weight_log_error_threshold, 1.0, 1e1,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_simplex_cost_perturbation_multiplier",
        "Dual simplex cost perturbation multiplier: 0 => no perturbation",
        advanced, &dual_simplex_cost_perturbation_multiplier, 0.0, 1.0,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "primal_simplex_bound_perturbation_multiplier",
        "Primal simplex bound perturbation multiplier: 0 => no perturbation",
        advanced, &primal_simplex_bound_perturbation_multiplier, 0.0, 1.0,
        kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_simplex_pivot_growth_tolerance",
        "Dual simplex pivot growth tolerance", advanced,
        &dual_simplex_pivot_growth_tolerance, 1e-12, 1e-9, kHighsInf);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "presolve_pivot_threshold",
        "Matrix factorization pivot threshold for substitutions in presolve",
        advanced, &presolve_pivot_threshold, kMinPivotThreshold, 0.01,
        kMaxPivotThreshold);
    records.push_back(record_double);

    record_int = new OptionRecordInt(
        "presolve_reduction_limit",
        "Limit on number of presolve reductions -1 => no limit", advanced,
        &presolve_reduction_limit, -1, -1, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "restart_presolve_reduction_limit",
        "Limit on number of further presolve reductions on restart in MIP "
        "solver -1 => no limit, otherwise, must be positive",
        advanced, &restart_presolve_reduction_limit, -1, -1, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "presolve_rule_off", "Bit mask of presolve rules that are not allowed",
        advanced, &presolve_rule_off, 0, 0, kHighsIInf);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "presolve_rule_test", "Presolve rule to test - DEV only!", advanced,
        &presolve_rule_test, 0, 0, kPresolveRuleMax);
    records.push_back(record_int);

    record_bool = new OptionRecordBool(
        "presolve_rule_logging", "Log effectiveness of presolve rules for LP",
        advanced, &presolve_rule_logging, false);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("presolve_remove_slacks",
                             "Remove slacks after presolve: redundant option!",
                             advanced, &presolve_remove_slacks, false);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("presolve_fm_obj_reformulation",
                             "Fourier-Motzkin objective reformulation",
                             advanced, &presolve_fm_obj_reformulation, false);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "presolve_substitution_maxfillin",
        "Maximal fillin allowed for substitutions in presolve", advanced,
        &presolve_substitution_maxfillin, 0, 10, kHighsIInf);
    records.push_back(record_int);

    record_double = new OptionRecordDouble(
        "factor_pivot_threshold", "Matrix factorization pivot threshold",
        advanced, &factor_pivot_threshold, kMinPivotThreshold,
        kDefaultPivotThreshold, kMaxPivotThreshold);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "factor_pivot_tolerance", "Matrix factorization pivot tolerance",
        advanced, &factor_pivot_tolerance, kMinPivotTolerance,
        kDefaultPivotTolerance, kMaxPivotTolerance);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "start_crossover_tolerance",
        "Tolerance to be satisfied before IPM crossover will start", advanced,
        &start_crossover_tolerance, kMinimumIpmTolerance, kDefaultIpmTolerance,
        kHighsInf);
    records.push_back(record_double);

    record_bool = new OptionRecordBool(
        "use_original_HFactor_logic",
        "Use original HFactor logic for sparse vs hyper-sparse TRANs", advanced,
        &use_original_HFactor_logic, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool(
        "less_infeasible_DSE_check", "Check whether LP is candidate for LiDSE",
        advanced, &less_infeasible_DSE_check, true);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("less_infeasible_DSE_choose_row",
                             "Use LiDSE if LP has right properties", advanced,
                             &less_infeasible_DSE_choose_row, true);
    records.push_back(record_bool);

    record_bool =
        new OptionRecordBool("run_centring", "Perform centring steps or not",
                             advanced, &run_centring, false);
    records.push_back(record_bool);

    /*
    record_bool = new OptionRecordBool("allow_pdlp_cleanup",
                                       "Allow PDLP to be used to clean up "
                                       "model with unknown status and no basis",
                                       advanced, &allow_pdlp_cleanup, true);
    records.push_back(record_bool);
    */

    record_int =
        new OptionRecordInt("max_centring_steps",
                            "Maximum number of steps to use (default = 5) "
                            "when computing the analytic centre",
                            advanced, &max_centring_steps, 0, 5, kHighsIInf);
    records.push_back(record_int);

    record_double = new OptionRecordDouble(
        "centring_ratio_tolerance",
        "Centring stops when the ratio max(x_j*s_j) / min(x_j*s_j) is below "
        "this tolerance (default = 100)",
        advanced, &centring_ratio_tolerance, 0, 100, kHighsInf);
    records.push_back(record_double);

    // Set up the log_options aliases
    log_options.clear();
    log_options.log_stream =
        log_file.empty() ? nullptr : fopen(log_file.c_str(), "w");
    log_options.output_flag = &output_flag;
    log_options.log_to_console = &log_to_console;
    log_options.log_dev_level = &log_dev_level;
  }

  void deleteRecords() {
    for (auto record : records) delete record;
  }

 public:
  std::vector<OptionRecord*> records;
  HighsInt num_user_settable_options_;
  void setLogOptions();
  double getMipPrimalFeasibilityTolerance() const {
    return this->mip_feasibility_tolerance;
  }
  double getMipDualFeasibilityTolerance() const {
    return this->mip_feasibility_tolerance * 0.1;
  }
};

void warnSolverInvalid(const HighsOptions& options,
                       const std::string& problem_type);

bool solverValidForLp(const std::string& solver);
bool solverValidForMip(const std::string& solver);
bool solverValidForQp(const std::string& solver);

bool useIpm(const std::string& solver);
bool usePdlp(const std::string& solver);
bool mayRequireDualPostsolve(const HighsOptions* options);
bool mayRequireBasisPostsolve(const HighsOptions* options);
#endif
