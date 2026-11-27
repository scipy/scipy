/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HStruct.h
 * @brief Structs for HiGHS
 */
#ifndef LP_DATA_HSTRUCT_H_
#define LP_DATA_HSTRUCT_H_

#include <functional>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsTimer.h"

struct HighsSolution {
  bool value_valid = false;
  bool dual_valid = false;
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
  bool hasUndefined() const;
  void invalidate();
  void clear();
  void print(const std::string& prefix = "",
             const std::string& message = "") const;
};

struct HighsObjectiveSolution {
  double objective;
  std::vector<double> col_value;
  void clear();
};

struct RefactorInfo {
  bool use = false;
  std::vector<HighsInt> pivot_row;
  std::vector<HighsInt> pivot_var;
  std::vector<int8_t> pivot_type;
  double build_synthetic_tick;
  void clear();
};

// Unused, but retained since there is a const reference to this in a
// deprecated method
struct HotStart {
  bool valid = false;
  RefactorInfo refactor_info;
  std::vector<int8_t> nonbasicMove;
};

struct HighsBasis {
  // Logical flags for a HiGHS basis:
  //
  // valid: has been factored by HiGHS
  //
  // alien: a basis that's been set externally, so cannot be assumed
  // to even have the right number of basic and nonbasic variables
  //
  // useful: a basis that may be useful
  //
  // Need useful since, by default, a basis is alien but not useful
  bool valid = false;
  bool alien = true;
  bool useful = false;
  bool was_alien = true;
  HighsInt debug_id = -1;
  HighsInt debug_update_count = -1;
  std::string debug_origin_name = "None";
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
  void print(const std::string& prefix = "",
             const std::string& message = "") const;
  void printScalars(const std::string& prefix = "",
                    const std::string& message = "") const;
  void invalidate();
  void clear();
};

struct HighsScale {
  HighsInt strategy;
  bool has_scaling;
  HighsInt num_col;
  HighsInt num_row;
  double cost;
  std::vector<double> col;
  std::vector<double> row;
};

struct HighsLpMods {
  // Semi-variables with zero lower bound that are treated as non-semi
  std::vector<HighsInt> save_non_semi_variable_index;

  // Semi-variables with inconsistent bounds that are fixed at zero
  std::vector<HighsInt> save_inconsistent_semi_variable_index;
  std::vector<double> save_inconsistent_semi_variable_lower_bound_value;
  std::vector<double> save_inconsistent_semi_variable_upper_bound_value;
  std::vector<HighsVarType> save_inconsistent_semi_variable_type;

  // Semi-variables whose lower bound is ignored when solving the
  // relaxation
  std::vector<HighsInt> save_relaxed_semi_variable_lower_bound_index;
  std::vector<double> save_relaxed_semi_variable_lower_bound_value;

  // Semi-variables whose upper bound is too large to be used as a
  // big-M when converting them to an integer variables plus an
  // integer/continuous variables as appropriate
  std::vector<HighsInt> save_tightened_semi_variable_upper_bound_index;
  std::vector<double> save_tightened_semi_variable_upper_bound_value;

  // Variables with infinite costs that are fixed during solve
  std::vector<HighsInt> save_inf_cost_variable_index;
  std::vector<double> save_inf_cost_variable_cost;
  std::vector<double> save_inf_cost_variable_lower;
  std::vector<double> save_inf_cost_variable_upper;

  void clear();
  bool isClear();
};

struct HighsNameHash {
  std::unordered_map<std::string, int> name2index;
  size_t size() { return this->name2index.size(); }
  void form(const std::vector<std::string>& names);
  bool hasDuplicate(const std::vector<std::string>& names);
  void addName(const HighsInt index, const std::string& name);
  bool addNameFindsDuplicate(const HighsInt index, const std::string& name);
  bool updateFindsDuplicate(const HighsInt index, const std::string& old_name,
                            const std::string& new_name);
  bool ok(const std::vector<std::string>& names);
  void clear() { this->name2index.clear(); }
};

struct HighsPresolveRuleLog {
  HighsInt call;
  HighsInt col_removed;
  HighsInt row_removed;
};

struct HighsPresolveLog {
  std::vector<HighsPresolveRuleLog> rule;
  void clear();
};

struct HighsIllConditioningRecord {
  HighsInt index;
  double multiplier;
};

struct HighsIllConditioning {
  std::vector<HighsIllConditioningRecord> record;
  void clear();
};

struct HighsLinearObjective {
  double weight = 0;
  double offset = 0;
  std::vector<double> coefficients;
  double abs_tolerance = -1;
  double rel_tolerance = -1;
  HighsInt priority = 0;
  void clear();
};

struct HighsProfilingRecord {
  std::vector<HighsInt> num_call;
  std::vector<double> run_time;
  std::vector<double> start_time;
};

struct HighsProfiling {
  HighsTimer* timer = nullptr;
  bool multi_threaded = true;
  std::string model_name_ = "";
  bool sub_solver_ = false;
  bool mip_ = false;
  HighsInt num_profiling_clock_ = -1;
  std::vector<std::string> name;
  // These vectors are over threads
  std::vector<HighsBool> submip;
  std::vector<HighsProfilingRecord> record;
  std::vector<HighsProfilingRecord> submip_record;
  bool initialized = false;

  void initialize(HighsTimer& timer_, const bool subsolver_profiling,
                  const bool mip_profiling = false);
  void clear();
  HighsInt numThread();
  HighsInt myThread();
  void setSubMip(const bool submip);
  bool isSubMip();
  HighsProfilingRecord* getHighsProfilingRecord(
      const HighsInt record_type = kChooseRecord);
  void start(const HighsInt profiling_clock, const bool restart = false);
  void stop(const HighsInt profiling_clock);
  double read(const HighsInt profiling_clock,
              const HighsInt record_type = kChooseRecord);
  bool running(const HighsInt profiling_clock,
               const HighsInt record_type = kChooseRecord);
  HighsInt numCall(const HighsInt profiling_clock,
                   const HighsInt record_type = kChooseRecord);
  void solveCall(const std::string& model, const bool submip);
  //  HighsInt getSepaClockIndex(const std::string& name);
};

// Starts a profiling clock and stops it when the scope is left, on every exit
using HighsProfilingScope = HighsClockScope<HighsProfiling>;

struct HighsSimplexStats {
  bool valid;
  HighsInt iteration_count;
  HighsInt num_invert;
  HighsInt last_invert_num_el;
  HighsInt last_factored_basis_num_el;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  void report(FILE* file, const std::string& message = "") const;
  void initialise(const HighsInt iteration_count_ = 0);
};

struct HighsUserScaleData {
  HighsInt user_objective_scale = 0;
  HighsInt user_bound_scale = 0;
  double infinite_cost;
  double infinite_bound;
  double small_matrix_value;
  double large_matrix_value;
  HighsInt num_infinite_costs;
  HighsInt num_infinite_hessian_values;
  HighsInt num_infinite_col_bounds;
  HighsInt num_infinite_row_bounds;
  HighsInt num_small_matrix_values;
  HighsInt num_large_matrix_values;
  HighsInt suggested_user_objective_scale;
  HighsInt suggested_user_bound_scale;
  bool applied = false;
  void initialise(const HighsInt& user_objective_scale_,
                  const HighsInt& user_bound_scale_,
                  const double& infinite_cost_, const double& infinite_bound_,
                  const double& small_matrix_value_,
                  const double& large_matrix_value_);
  bool scaleError(std::string& message) const;
  bool scaleWarning(std::string& message) const;
};

using HighsHessianFunctionType = std::function<HighsInt(
    const HighsInt call_type, const HighsInt* x_num_entries,
    const HighsInt* x_index, const double* x_value,
    HighsInt* hessian_x_num_entries, HighsInt* hessian_x_index,
    double* hessian_x_value, void*)>;

struct HessianOracle {
  // Oracle to obtain values of Q' = multiplier_*Q + shift_*I and form
  // products with Q', by calling a user-supplied function (call_) to
  // form products with Q and then apply any non-trivial multiplier_
  // or shift_
  HighsInt dim_ = 0;
  double multiplier_ = 1.0;  // Minimize
  double shift_ = 0.0;       // No regularization
  HighsHessianFunctionType call_ = nullptr;
  void* data_ = nullptr;
  void clear();
  bool isValid() const { return hasProductCall(); }
  bool hasProductCall() const;
  void formFromOracle();
  double diagonal(const HighsInt i) const;
  double entry(const HighsInt i, const HighsInt j) const;
  void getPackedColumn(const HighsInt col, HighsInt& col_num_entries,
                       HighsInt* col_index, double* col_value) const;
  void getScatteredColumn(const HighsInt col, HighsInt& col_num_entries,
                          HighsInt* col_index, double* col_value) const;
  void product(const std::vector<double>& x_value,
               std::vector<double>& hessian_x_value) const;
  void product(const double* x_value, double* hessian_x_value) const;
  void product(const HighsInt x_num_entries, const HighsInt* x_index,
               const double* x_value, double* hessian_x_value) const;
  void scaleAndShift(const HighsInt* x_num_entries, const HighsInt* x_index,
                     const double* x_value, double* hessian_x_value) const;
};

struct StatusString {
  std::string full_;
  std::string s2_;
};

#endif /* LP_DATA_HSTRUCT_H_ */
