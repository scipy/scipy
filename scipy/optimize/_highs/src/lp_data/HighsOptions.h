/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_OPTIONS_H_
#define LP_DATA_HIGHS_OPTIONS_H_

#include <cstring>  // For strlen
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"
#include "simplex/SimplexConst.h"

using std::string;

enum class OptionStatus { OK = 0, NO_FILE, UNKNOWN_OPTION, ILLEGAL_VALUE };

class OptionRecord {
 public:
  HighsOptionType type;
  std::string name;
  std::string description;
  bool advanced;

  OptionRecord(HighsOptionType Xtype, std::string Xname,
               std::string Xdescription, bool Xadvanced) {
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
  OptionRecordBool(std::string Xname, std::string Xdescription, bool Xadvanced,
                   bool* Xvalue_pointer, bool Xdefault_value)
      : OptionRecord(HighsOptionType::BOOL, Xname, Xdescription, Xadvanced) {
    advanced = Xadvanced;
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  void assignvalue(bool Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordBool() {}
};

class OptionRecordInt : public OptionRecord {
 public:
  int* value;
  int lower_bound;
  int default_value;
  int upper_bound;
  OptionRecordInt(std::string Xname, std::string Xdescription, bool Xadvanced,
                  int* Xvalue_pointer, int Xlower_bound, int Xdefault_value,
                  int Xupper_bound)
      : OptionRecord(HighsOptionType::INT, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    lower_bound = Xlower_bound;
    default_value = Xdefault_value;
    upper_bound = Xupper_bound;
    *value = default_value;
  }

  void assignvalue(int Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordInt() {}
};

class OptionRecordDouble : public OptionRecord {
 public:
  double* value;
  double lower_bound;
  double upper_bound;
  double default_value;
  OptionRecordDouble(std::string Xname, std::string Xdescription,
                     bool Xadvanced, double* Xvalue_pointer,
                     double Xlower_bound, double Xdefault_value,
                     double Xupper_bound)
      : OptionRecord(HighsOptionType::DOUBLE, Xname, Xdescription, Xadvanced) {
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
  OptionRecordString(std::string Xname, std::string Xdescription,
                     bool Xadvanced, std::string* Xvalue_pointer,
                     std::string Xdefault_value)
      : OptionRecord(HighsOptionType::STRING, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  void assignvalue(std::string Xvalue) { *value = Xvalue; }

  virtual ~OptionRecordString() {}
};

inline const char* bool2string(bool b);

bool commandLineOffChooseOnOk(FILE* logfile, const string& value);
bool commandLineSolverOk(FILE* logfile, const string& value);

bool boolFromString(const std::string value, bool& bool_value);

OptionStatus getOptionIndex(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            int& index);

OptionStatus checkOptions(FILE* logfile,
                          const std::vector<OptionRecord*>& option_records);
OptionStatus checkOption(FILE* logfile, const OptionRecordInt& option);
OptionStatus checkOption(FILE* logfile, const OptionRecordDouble& option);

OptionStatus checkOptionValue(FILE* logfile,
                              std::vector<OptionRecord*>& option_records,
                              const int value);
OptionStatus checkOptionValue(FILE* logfile,
                              std::vector<OptionRecord*>& option_records,
                              const double value);
OptionStatus checkOptionValue(FILE* logfile,
                              std::vector<OptionRecord*>& option_records,
                              const std::string value);

OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const bool value);
OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const int value);
OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const double value);
OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const std::string value);
OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const char* value);

OptionStatus setOptionValue(OptionRecordBool& option, const bool value);
OptionStatus setOptionValue(FILE* logfile, OptionRecordInt& option,
                            const int value);
OptionStatus setOptionValue(FILE* logfile, OptionRecordDouble& option,
                            const double value);
OptionStatus setOptionValue(FILE* logfile, OptionRecordString& option,
                            std::string const value);

OptionStatus passOptions(FILE* logfile, const HighsOptions& from_options,
                         HighsOptions& to_options);

OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            bool& value);
OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            int& value);
OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            double& value);
OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            std::string& value);

void resetOptions(std::vector<OptionRecord*>& option_records);

HighsStatus writeOptionsToFile(FILE* file,
                               const std::vector<OptionRecord*>& option_records,
                               const bool report_only_non_default_values = true,
                               const bool html = false);
void reportOptions(FILE* file, const std::vector<OptionRecord*>& option_records,
                   const bool report_only_non_default_values = true,
                   const bool html = false);
void reportOption(FILE* file, const OptionRecordBool& option,
                  const bool report_only_non_default_values, const bool html);
void reportOption(FILE* file, const OptionRecordInt& option,
                  const bool report_only_non_default_values, const bool html);
void reportOption(FILE* file, const OptionRecordDouble& option,
                  const bool report_only_non_default_values, const bool html);
void reportOption(FILE* file, const OptionRecordString& option,
                  const bool report_only_non_default_values, const bool html);

const string simplex_string = "simplex";
const string ipm_string = "ipm";
const int KEEP_N_ROWS_DELETE_ROWS = -1;
const int KEEP_N_ROWS_DELETE_ENTRIES = 0;
const int KEEP_N_ROWS_KEEP_ROWS = 1;

// Strings for command line options
const string model_file_string = "model_file";
const string presolve_string = "presolve";
const string solver_string = "solver";
const string parallel_string = "parallel";
const string time_limit_string = "time_limit";
const string options_file_string = "options_file";

// enum objSense { OBJSENSE_MINIMIZE = 1, OBJSENSE_MAXIMIZE = -1 };

struct HighsOptionsStruct {
  // Options read from the command line
  std::string model_file;
  std::string presolve;
  std::string solver;
  std::string parallel;
  double time_limit;
  std::string options_file;

  // Options read from the file
  double infinite_cost;
  double infinite_bound;
  double small_matrix_value;
  double large_matrix_value;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  double ipm_optimality_tolerance;
  double dual_objective_value_upper_bound;
  int highs_debug_level;
  int simplex_strategy;
  int simplex_scale_strategy;
  int simplex_crash_strategy;
  int simplex_dual_edge_weight_strategy;
  int simplex_primal_edge_weight_strategy;
  int simplex_iteration_limit;
  int simplex_update_limit;
  int ipm_iteration_limit;
  int highs_min_threads;
  int highs_max_threads;
  int message_level;
  std::string solution_file;
  bool write_solution_to_file;
  bool write_solution_pretty;

  // Advanced options
  bool run_crossover;
  bool mps_parser_type_free;
  int keep_n_rows;
  int allowed_simplex_matrix_scale_factor;
  int allowed_simplex_cost_scale_factor;
  int simplex_dualise_strategy;
  int simplex_permute_strategy;
  int dual_simplex_cleanup_strategy;
  int simplex_price_strategy;
  int dual_chuzc_sort_strategy;
  bool simplex_initial_condition_check;
  double simplex_initial_condition_tolerance;
  double dual_steepest_edge_weight_log_error_threshhold;
  double dual_simplex_cost_perturbation_multiplier;
  bool less_infeasible_DSE_check;
  bool less_infeasible_DSE_choose_row;
  bool use_original_HFactor_logic;

  // Options for MIP solver
  int mip_max_nodes;
  int mip_report_level;

  // Switch for MIP solver
  bool mip;

  // Options for HighsPrintMessage and HighsLogMessage
  FILE* logfile = stdout;
  FILE* output = stdout;

  void (*printmsgcb)(int level, const char* msg, void* msgcb_data) = NULL;
  void (*logmsgcb)(HighsMessageType type, const char* msg,
                   void* msgcb_data) = NULL;
  void* msgcb_data = NULL;

  virtual ~HighsOptionsStruct() {}
};

// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily. The options below
// are just what has been used to parse options from argv.
// todo: when creating the new options don't forget underscores for class
// variables but no underscores for struct
class HighsOptions : public HighsOptionsStruct {
 public:
  HighsOptions() { initRecords(); }

  HighsOptions(const HighsOptions& options) {
    initRecords();
    HighsOptionsStruct::operator=(options);
  }

  HighsOptions(HighsOptions&& options) {
    records = std::move(options.records);
    HighsOptionsStruct::operator=(std::move(options));
  }

  const HighsOptions& operator=(const HighsOptions& other) {
    if (&other != this) {
      if ((int)records.size() == 0) initRecords();
      HighsOptionsStruct::operator=(other);
    }
    return *this;
  }

  const HighsOptions& operator=(HighsOptions&& other) {
    if (&other != this) {
      if ((int)records.size() == 0) initRecords();
      HighsOptionsStruct::operator=(other);
    }
    return *this;
  }

  virtual ~HighsOptions() {
    if (records.size() > 0) deleteRecords();
  }

 private:
  void initRecords() {
    OptionRecordBool* record_bool;
    OptionRecordInt* record_int;
    OptionRecordDouble* record_double;
    OptionRecordString* record_string;
    bool advanced;
    advanced = false;
    // Options read from the command line
    record_string =
        new OptionRecordString(model_file_string, "Model file", advanced,
                               &model_file, FILENAME_DEFAULT);
    records.push_back(record_string);
    record_string = new OptionRecordString(
        presolve_string, "Presolve option: \"off\", \"choose\" or \"on\"",
        advanced, &presolve, choose_string);
    records.push_back(record_string);
    record_string = new OptionRecordString(
        solver_string, "Solver option: \"simplex\", \"choose\" or \"ipm\"",
        advanced, &solver, choose_string);
    records.push_back(record_string);
    record_string = new OptionRecordString(
        parallel_string, "Parallel option: \"off\", \"choose\" or \"on\"",
        advanced, &parallel, choose_string);
    records.push_back(record_string);
    record_double = new OptionRecordDouble(time_limit_string, "Time limit",
                                           advanced, &time_limit, 0,
                                           HIGHS_CONST_INF, HIGHS_CONST_INF);
    records.push_back(record_double);
    record_string =
        new OptionRecordString(options_file_string, "Options file", advanced,
                               &options_file, FILENAME_DEFAULT);
    records.push_back(record_string);
    // Options read from the file
    record_double =
        new OptionRecordDouble("infinite_cost",
                               "Limit on cost coefficient: values larger than "
                               "this will be treated as infinite",
                               advanced, &infinite_cost, 1e15, 1e20, 1e25);
    records.push_back(record_double);

    record_double =
        new OptionRecordDouble("infinite_bound",
                               "Limit on |constraint bound|: values larger "
                               "than this will be treated as infinite",
                               advanced, &infinite_bound, 1e15, 1e20, 1e25);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "small_matrix_value",
        "Lower limit on |matrix entries|: values smaller than this will be "
        "treated as zero",
        advanced, &small_matrix_value, 1e-12, 1e-9, HIGHS_CONST_INF);
    records.push_back(record_double);

    record_double =
        new OptionRecordDouble("large_matrix_value",
                               "Upper limit on |matrix entries|: values larger "
                               "than this will be treated as infinite",
                               advanced, &large_matrix_value, 1e0, 1e15, 1e20);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "primal_feasibility_tolerance", "Primal feasibility tolerance",
        advanced, &primal_feasibility_tolerance, 1e-10, 1e-7, HIGHS_CONST_INF);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_feasibility_tolerance", "Dual feasibility tolerance", advanced,
        &dual_feasibility_tolerance, 1e-10, 1e-7, HIGHS_CONST_INF);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "ipm_optimality_tolerance", "IPM optimality tolerance", advanced,
        &ipm_optimality_tolerance, 1e-12, 1e-8, HIGHS_CONST_INF);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_objective_value_upper_bound",
        "Upper bound on objective value for dual simplex: algorithm terminates "
        "if reached",
        advanced, &dual_objective_value_upper_bound, -HIGHS_CONST_INF,
        HIGHS_CONST_INF, HIGHS_CONST_INF);
    records.push_back(record_double);

    record_int =
        new OptionRecordInt("highs_debug_level", "Debugging level in HiGHS",
                            advanced, &highs_debug_level, HIGHS_DEBUG_LEVEL_MIN,
                            HIGHS_DEBUG_LEVEL_MIN, HIGHS_DEBUG_LEVEL_MAX);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("simplex_strategy", "Strategy for simplex solver",
                            advanced, &simplex_strategy, SIMPLEX_STRATEGY_MIN,
                            SIMPLEX_STRATEGY_DUAL, SIMPLEX_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_scale_strategy",
        "Strategy for scaling before simplex solver: off / on (0/1)", advanced,
        &simplex_scale_strategy, SIMPLEX_SCALE_STRATEGY_MIN,
        SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED, SIMPLEX_SCALE_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_crash_strategy",
        "Strategy for simplex crash: off / LTSSF / Bixby (0/1/2)", advanced,
        &simplex_crash_strategy, SIMPLEX_CRASH_STRATEGY_MIN,
        SIMPLEX_CRASH_STRATEGY_OFF, SIMPLEX_CRASH_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_dual_edge_weight_strategy",
        "Strategy for simplex dual edge weights: Dantzig / Devex / Steepest "
        "Edge (0/1/2)",
        advanced, &simplex_dual_edge_weight_strategy,
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN,
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH,
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_primal_edge_weight_strategy",
        "Strategy for simplex primal edge weights: Dantzig / Devex (0/1)",
        advanced, &simplex_primal_edge_weight_strategy,
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN,
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt("simplex_iteration_limit",
                                     "Iteration limit for simplex solver",
                                     advanced, &simplex_iteration_limit, 0,
                                     HIGHS_CONST_I_INF, HIGHS_CONST_I_INF);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_update_limit",
        "Limit on the number of simplex UPDATE operations", advanced,
        &simplex_update_limit, 0, 5000, HIGHS_CONST_I_INF);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "ipm_iteration_limit", "Iteration limit for IPM solver", advanced,
        &ipm_iteration_limit, 0, HIGHS_CONST_I_INF, HIGHS_CONST_I_INF);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "highs_min_threads", "Minimum number of threads in parallel execution",
        advanced, &highs_min_threads, 1, 1, HIGHS_THREAD_LIMIT);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "highs_max_threads", "Maximum number of threads in parallel execution",
        advanced, &highs_max_threads, 1, HIGHS_THREAD_LIMIT,
        HIGHS_THREAD_LIMIT);
    records.push_back(record_int);

    record_int = new OptionRecordInt("message_level",
                                     "HiGHS message level: bit-mask 1 => "
                                     "VERBOSE; 2 => DETAILED 4 => MINIMAL",
                                     advanced, &message_level, ML_MIN,
                                     ML_MINIMAL, ML_MAX);
    records.push_back(record_int);

    record_string =
        new OptionRecordString("solution_file", "Solution file", advanced,
                               &solution_file, FILENAME_DEFAULT);
    records.push_back(record_string);

    record_bool =
        new OptionRecordBool("write_solution_to_file",
                             "Write the primal and dual solution to a file",
                             advanced, &write_solution_to_file, false);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("write_solution_pretty",
                                       "Write the primal and dual solution in "
                                       "a pretty (human-readable) format",
                                       advanced, &write_solution_pretty, false);
    records.push_back(record_bool);

    record_int = new OptionRecordInt(
        "mip_max_nodes", "MIP solver max number of nodes", advanced,
        &mip_max_nodes, 0, HIGHS_CONST_I_INF, HIGHS_CONST_I_INF);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("mip_report_level", "MIP solver reporting level",
                            advanced, &mip_report_level, 0, 1, 2);
    records.push_back(record_int);

    // Advanced options
    advanced = true;

    record_bool = new OptionRecordBool("run_crossover",
                                       "Run the crossover routine for IPX",
                                       advanced, &run_crossover, true);
    records.push_back(record_bool);

    record_bool = new OptionRecordBool("mps_parser_type_free",
                                       "Use the free format MPS file reader",
                                       advanced, &mps_parser_type_free, true);
    records.push_back(record_bool);
    record_int =
        new OptionRecordInt("keep_n_rows",
                            "For multiple N-rows in MPS files: delete rows / "
                            "delete entries / keep rows (-1/0/1)",
                            advanced, &keep_n_rows, KEEP_N_ROWS_DELETE_ROWS,
                            KEEP_N_ROWS_DELETE_ROWS, KEEP_N_ROWS_KEEP_ROWS);
    records.push_back(record_int);
    record_int = new OptionRecordInt(
        "allowed_simplex_matrix_scale_factor",
        "Largest power-of-two factor permitted when scaling the constraint "
        "matrix for the simplex solver",
        advanced, &allowed_simplex_matrix_scale_factor, 0, 10, 20);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "allowed_simplex_cost_scale_factor",
        "Largest power-of-two factor permitted when scaling the costs for the "
        "simplex solver",
        advanced, &allowed_simplex_cost_scale_factor, 0, 0, 20);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_dualise_strategy", "Strategy for dualising before simplex",
        advanced, &simplex_dualise_strategy, OPTION_OFF, OPTION_OFF, OPTION_ON);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_permute_strategy", "Strategy for permuting before simplex",
        advanced, &simplex_permute_strategy, OPTION_OFF, OPTION_OFF, OPTION_ON);
    records.push_back(record_int);

    record_int =
        new OptionRecordInt("dual_simplex_cleanup_strategy",
                            "Strategy for cleanup in dual simplex solver: none "
                            "/ HPrimal / HQPrimal (0/1/2)",
                            advanced, &dual_simplex_cleanup_strategy,
                            DUAL_SIMPLEX_CLEANUP_STRATEGY_MIN,
                            DUAL_SIMPLEX_CLEANUP_STRATEGY_HPRIMAL,
                            DUAL_SIMPLEX_CLEANUP_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "simplex_price_strategy", "Strategy for PRICE in simplex", advanced,
        &simplex_price_strategy, SIMPLEX_PRICE_STRATEGY_MIN,
        SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH,
        SIMPLEX_PRICE_STRATEGY_MAX);
    records.push_back(record_int);

    record_int = new OptionRecordInt(
        "dual_chuzc_sort_strategy", "Strategy for CHUZC sort in dual simplex",
        advanced, &dual_chuzc_sort_strategy, SIMPLEX_DUAL_CHUZC_STRATEGY_MIN,
        SIMPLEX_DUAL_CHUZC_STRATEGY_CHOOSE, SIMPLEX_DUAL_CHUZC_STRATEGY_MAX);
    records.push_back(record_int);

    record_bool =
        new OptionRecordBool("simplex_initial_condition_check",
                             "Perform initial basis condition check in simplex",
                             advanced, &simplex_initial_condition_check, true);
    records.push_back(record_bool);

    record_double = new OptionRecordDouble(
        "simplex_initial_condition_tolerance",
        "Tolerance on initial basis condition in simplex", advanced,
        &simplex_initial_condition_tolerance, 1.0, 1e14, HIGHS_CONST_INF);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_steepest_edge_weight_log_error_threshhold",
        "Threshhold on dual steepest edge weight errors for Devex switch",
        advanced, &dual_steepest_edge_weight_log_error_threshhold, 1.0, 1e1,
        HIGHS_CONST_INF);
    records.push_back(record_double);

    record_double = new OptionRecordDouble(
        "dual_simplex_cost_perturbation_multiplier",
        "Dual simplex cost perturbation multiplier: 0 => no perturbation",
        advanced, &dual_simplex_cost_perturbation_multiplier, 0.0, 1.0,
        HIGHS_CONST_INF);
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
        new OptionRecordBool("mip", "Use mip solver.", advanced, &mip, false);
    records.push_back(record_bool);
  }

  void deleteRecords() {
    for (unsigned int i = 0; i < records.size(); i++) delete records[i];
  }

 public:
  std::vector<OptionRecord*> records;
};

#endif
