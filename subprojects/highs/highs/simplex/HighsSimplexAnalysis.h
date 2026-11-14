/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexAnalysis.h
 * @brief Analyse simplex iterations, both for run-time control and data
 * gathering
 */
#ifndef SIMPLEX_HIGHSSIMPLEXANALYSIS_H_
#define SIMPLEX_HIGHSSIMPLEXANALYSIS_H_

#include <cassert>
#include <memory>
#include <sstream>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "simplex/SimplexConst.h"
#include "util/HFactor.h"

enum TRAN_STAGE {
  TRAN_STAGE_FTRAN_LOWER = 0,
  TRAN_STAGE_FTRAN_UPPER_FT,
  TRAN_STAGE_FTRAN_UPPER,
  TRAN_STAGE_BTRAN_UPPER,
  TRAN_STAGE_BTRAN_UPPER_FT,
  TRAN_STAGE_BTRAN_LOWER,
  NUM_TRAN_STAGE_TYPE,
};

struct TranStageAnalysis {
  std::string name_;
  HighsScatterData rhs_density_;
  HighsInt num_decision_;
  HighsInt num_wrong_original_sparse_decision_;
  HighsInt num_wrong_original_hyper_decision_;
  HighsInt num_wrong_new_sparse_decision_;
  HighsInt num_wrong_new_hyper_decision_;
};

const HighsInt kAnIterTraceMaxNumRec = 20;
const HighsLogType kIterationReportLogType = HighsLogType::kVerbose;

/**
 * @brief Analyse simplex iterations, both for run-time control and data
 * gathering
 */
class HighsSimplexAnalysis {
 public:
  HighsSimplexAnalysis()
      : timer_(nullptr),
        pointer_serial_factor_clocks(nullptr),
        numRow(0),
        numCol(0),
        numTot(0),
        model_name_(""),
        lp_name_(""),
        analyse_lp_data(false),
        analyse_simplex_summary_data(false),
        analyse_simplex_runtime_data(false),
        analyse_simplex_time(false),
        analyse_factor_data(false),
        analyse_factor_time(false),
        analyse_simplex_data(false),
        simplex_strategy(0),
        edge_weight_mode(EdgeWeightMode::kSteepestEdge),
        solve_phase(0),
        simplex_iteration_count(0),
        devex_iteration_count(0),
        pivotal_row_index(0),
        leaving_variable(0),
        entering_variable(0),
        rebuild_reason(0),
        rebuild_reason_string(""),
        reduced_rhs_value(0.0),
        reduced_cost_value(0.0),
        edge_weight(0.0),
        edge_weight_error(0.0),
        primal_delta(0.0),
        primal_step(0.0),
        dual_step(0.0),
        pivot_value_from_column(0.0),
        pivot_value_from_row(0.0),
        factor_pivot_threshold(0.0),
        numerical_trouble(0.0),
        objective_value(0.0),
        num_primal_infeasibility(0),
        num_dual_infeasibility(0),
        sum_primal_infeasibility(0.0),
        sum_dual_infeasibility(0.0),
        num_dual_phase_1_lp_dual_infeasibility(0),
        max_dual_phase_1_lp_dual_infeasibility(0.0),
        sum_dual_phase_1_lp_dual_infeasibility(0.0),
        num_devex_framework(0),
        col_aq_density(0.0),
        row_ep_density(0.0),
        row_ap_density(0.0),
        row_DSE_density(0.0),
        col_steepest_edge_density(0.0),
        col_basic_feasibility_change_density(0.0),
        row_basic_feasibility_change_density(0.0),
        col_BFRT_density(0.0),
        primal_col_density(0.0),
        dual_col_density(0.0),
        num_costly_DSE_iteration(0),
        costly_DSE_measure(0.0),
        multi_iteration_count(0),
        multi_chosen(0),
        multi_finished(0),
        min_concurrency(0),
        num_concurrency(0),
        max_concurrency(0),
        num_col_price(0),
        num_row_price(0),
        num_row_price_with_switch(0),
        num_primal_cycling_detections(0),
        num_dual_cycling_detections(0),
        num_quad_chuzc(0),
        num_heap_chuzc(0),
        sum_quad_chuzc_size(0.0),
        sum_heap_chuzc_size(0.0),
        max_quad_chuzc_size(0),
        max_heap_chuzc_size(0),
        num_improve_choose_column_row_call(0),
        num_remove_pivot_from_pack(0),
        num_correct_dual_primal_flip(0),
        min_correct_dual_primal_flip_dual_infeasibility(kHighsInf),
        max_correct_dual_primal_flip(0.0),
        num_correct_dual_cost_shift(0),
        max_correct_dual_cost_shift_dual_infeasibility(0.0),
        max_correct_dual_cost_shift(0.0),
        net_num_single_cost_shift(0),
        num_single_cost_shift(0),
        max_single_cost_shift(0.0),
        sum_single_cost_shift(0.0),
        num_dual_steepest_edge_weight_check(0),
        num_dual_steepest_edge_weight_reject(0),
        num_wrong_low_dual_steepest_edge_weight(0),
        num_wrong_high_dual_steepest_edge_weight(0),
        average_frequency_low_dual_steepest_edge_weight(0.0),
        average_frequency_high_dual_steepest_edge_weight(0.0),
        average_log_low_dual_steepest_edge_weight_error(0.0),
        average_log_high_dual_steepest_edge_weight_error(0.0),
        max_average_frequency_low_dual_steepest_edge_weight(0.0),
        max_average_frequency_high_dual_steepest_edge_weight(0.0),
        max_sum_average_frequency_extreme_dual_steepest_edge_weight(0.0),
        max_average_log_low_dual_steepest_edge_weight_error(0.0),
        max_average_log_high_dual_steepest_edge_weight_error(0.0),
        max_sum_average_log_extreme_dual_steepest_edge_weight_error(0.0),
        num_invert_report_since_last_header(-1),
        num_iteration_report_since_last_header(-1),
        highs_run_time(0.0),
        last_user_log_time(-kHighsInf),
        delta_user_log_time(1e0),
        timeless_log(false),
        average_concurrency(0.0),
        average_fraction_of_possible_minor_iterations_performed(0.0),
        sum_multi_chosen(0),
        sum_multi_finished(0),
        num_invert(0),
        num_kernel(0),
        num_major_kernel(0),
        max_kernel_dim(0.0),
        sum_kernel_dim(0.0),
        running_average_kernel_dim(0.0),
        sum_invert_fill_factor(0.0),
        sum_kernel_fill_factor(0.0),
        sum_major_kernel_fill_factor(0.0),
        running_average_invert_fill_factor(1.0),
        running_average_kernel_fill_factor(1.0),
        running_average_major_kernel_fill_factor(1.0),
        AnIterIt0(0),
        AnIterPrevIt(0),
        AnIterOp{},
        AnIterTraceNumRec(0),
        AnIterTraceIterDl(0),
        AnIterTrace{},
        AnIterNumInvert{},
        AnIterNumEdWtIt{} {}

  // Pointer to timer
  HighsTimer* timer_;

  void setup(const std::string& lp_name, const HighsLp& lp,
             const HighsOptions& options,
             const HighsInt simplex_iteration_count);
  void setupSimplexTime(const HighsOptions& options);
  void setupFactorTime(const HighsOptions& options);
  void messaging(const HighsLogOptions& log_options_);
  void iterationReport();
  void invertReport();
  void invertReport(const bool header);
  void userInvertReport(const bool force);
  void userInvertReport(const bool header, const bool force);
  bool predictEndDensity(const HighsInt tran_stage_id,
                         const double start_density, double& end_density) const;
  void afterTranStage(const HighsInt tran_stage_id, const double start_density,
                      const double end_density, const double historical_density,
                      const double predicted_end_density,
                      const bool use_solve_sparse_original_HFactor_logic,
                      const bool use_solve_sparse_new_HFactor_logic);

  void simplexTimerStart(const HighsInt simplex_clock,
                         const HighsInt thread_id = 0);
  void simplexTimerStop(const HighsInt simplex_clock,
                        const HighsInt thread_id = 0);
  bool simplexTimerRunning(const HighsInt simplex_clock,
                           const HighsInt thread_id = 0) const;
  HighsInt simplexTimerNumCall(const HighsInt simplex_clock,
                               const HighsInt thread_id = 0) const;
  double simplexTimerRead(const HighsInt simplex_clock,
                          const HighsInt thread_id = 0) const;

  HighsTimerClock* getThreadFactorTimerClockPointer();

  const std::vector<HighsTimerClock>& getThreadSimplexTimerClocks() {
    return thread_simplex_clocks;
  }
  HighsTimerClock* getThreadSimplexTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_simplex_clocks.size());
    return &thread_simplex_clocks[i];
  }

  const std::vector<HighsTimerClock>& getThreadFactorTimerClocks() {
    return thread_factor_clocks;
  }
  HighsTimerClock* getThreadFactorTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_factor_clocks.size());
    return &thread_factor_clocks[i];
  }

  void iterationRecord();
  void iterationRecordMajor();
  void operationRecordBefore(const HighsInt operation_type,
                             const HVector& vector,
                             const double historical_density);
  void operationRecordBefore(const HighsInt operation_type,
                             const HighsInt current_count,
                             const double historical_density);
  void operationRecordAfter(const HighsInt operation_type,
                            const HVector& vector);
  void operationRecordAfter(const HighsInt operation_type,
                            const HighsInt result_count);
  void summaryReport();
  void summaryReportFactor() const;
  void reportSimplexTimer() const;
  void reportFactorTimer();
  void updateInvertFormData(const HFactor& factor);
  void reportInvertFormData() const;
  HighsInt numInvert() const { return num_invert; }

  // Control methods to be moved to HEkkControl
  void dualSteepestEdgeWeightError(const double computed_edge_weight,
                                   const double updated_edge_weight);
  //  bool switchToDevex();

  std::vector<HighsTimerClock> thread_simplex_clocks;
  std::vector<HighsTimerClock> thread_factor_clocks;
  HighsTimerClock* pointer_serial_factor_clocks;

  // Local copies of LP data
  HighsInt numRow;
  HighsInt numCol;
  HighsInt numTot;
  std::string model_name_;
  std::string lp_name_;

  // Local copies of IO data
  HighsLogOptions log_options;

  // Interpreted shortcuts from bit settings in highs_analysis_level
  bool analyse_lp_data;
  bool analyse_simplex_summary_data;
  bool analyse_simplex_runtime_data;
  bool analyse_simplex_time;
  bool analyse_factor_data;
  bool analyse_factor_time;
  bool analyse_simplex_data;

  // Control parameters moving to info
  //  bool allow_dual_steepest_edge_to_devex_switch;
  //  double dual_steepest_edge_weight_log_error_threshold;

  // Local copies of simplex data for reporting
  HighsInt simplex_strategy;
  EdgeWeightMode edge_weight_mode;
  HighsInt solve_phase;
  HighsInt simplex_iteration_count;
  HighsInt devex_iteration_count;
  HighsInt pivotal_row_index;
  HighsInt leaving_variable;
  HighsInt entering_variable;
  HighsInt rebuild_reason;
  std::string rebuild_reason_string;
  double reduced_rhs_value;
  double reduced_cost_value;
  double edge_weight;
  double edge_weight_error;
  double primal_delta;
  double primal_step;
  double dual_step;
  double pivot_value_from_column;
  double pivot_value_from_row;
  double factor_pivot_threshold;
  double numerical_trouble;
  double objective_value;
  HighsInt num_primal_infeasibility;
  HighsInt num_dual_infeasibility;
  double sum_primal_infeasibility;
  double sum_dual_infeasibility;
  // This triple is an original infeasibility record, so it includes max,
  // but it's only used for reporting
  HighsInt num_dual_phase_1_lp_dual_infeasibility;
  double max_dual_phase_1_lp_dual_infeasibility;
  double sum_dual_phase_1_lp_dual_infeasibility;
  HighsInt num_devex_framework;
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
  HighsInt num_costly_DSE_iteration;
  double costly_DSE_measure;

  // Local copies of parallel simplex data for reporting
  HighsInt multi_iteration_count;
  HighsInt multi_chosen;
  HighsInt multi_finished;
  HighsInt min_concurrency;
  HighsInt num_concurrency;
  HighsInt max_concurrency;

  // Unused
  //  HighsInt multi_num = 0; // Useless
  //  double basis_condition = 0; // Maybe useful

  // Records of how pivotal row PRICE was done
  HighsInt num_col_price;
  HighsInt num_row_price;
  HighsInt num_row_price_with_switch;

  HighsValueDistribution before_ftran_upper_sparse_density;
  HighsValueDistribution ftran_upper_sparse_density;
  HighsValueDistribution before_ftran_upper_hyper_density;
  HighsValueDistribution ftran_upper_hyper_density;
  HighsValueDistribution cost_perturbation1_distribution;
  HighsValueDistribution cost_perturbation2_distribution;
  HighsValueDistribution cleanup_dual_change_distribution;
  HighsValueDistribution cleanup_primal_step_distribution;
  HighsValueDistribution cleanup_dual_step_distribution;
  HighsValueDistribution cleanup_primal_change_distribution;

  HighsInt num_primal_cycling_detections;
  HighsInt num_dual_cycling_detections;

  HighsInt num_quad_chuzc;
  HighsInt num_heap_chuzc;
  double sum_quad_chuzc_size;
  double sum_heap_chuzc_size;
  HighsInt max_quad_chuzc_size;
  HighsInt max_heap_chuzc_size;

  HighsInt num_improve_choose_column_row_call;
  HighsInt num_remove_pivot_from_pack;

  HighsInt num_correct_dual_primal_flip;
  double min_correct_dual_primal_flip_dual_infeasibility;
  double max_correct_dual_primal_flip;
  HighsInt num_correct_dual_cost_shift;
  double max_correct_dual_cost_shift_dual_infeasibility;
  double max_correct_dual_cost_shift;
  HighsInt net_num_single_cost_shift;
  HighsInt num_single_cost_shift;
  double max_single_cost_shift;
  double sum_single_cost_shift;

  // Tolerances for analysis of TRAN stages - could be needed for
  // control if this is ever used again!
  vector<double> original_start_density_tolerance;
  vector<double> new_start_density_tolerance;
  vector<double> historical_density_tolerance;
  vector<double> predicted_density_tolerance;
  vector<TranStageAnalysis> tran_stage;

  std::unique_ptr<std::stringstream> analysis_log;

 private:
  void iterationReport(const bool header);
  void reportAlgorithmPhase(const bool header);
  void reportIterationObjective(const bool header);
  void reportInfeasibility(const bool header);
  void reportThreads(const bool header);
  void reportMulti(const bool header);
  void reportOneDensity(const double density);
  void printOneDensity(const double density) const;
  void reportDensity(const bool header);
  void reportInvert(const bool header);
  //  void reportCondition(const bool header);
  void reportIterationData(const bool header);
  void reportRunTime(const bool header, const double run_time);
  void reportFreeListSize(const bool header);
  HighsInt intLog10(const double v) const;
  bool dualAlgorithm() const;

  //  double AnIterCostlyDseFq;  //!< Frequency of iterations when DSE is costly
  //  double AnIterCostlyDseMeasure;

  HighsInt num_dual_steepest_edge_weight_check;
  HighsInt num_dual_steepest_edge_weight_reject;
  HighsInt num_wrong_low_dual_steepest_edge_weight;
  HighsInt num_wrong_high_dual_steepest_edge_weight;
  double average_frequency_low_dual_steepest_edge_weight;
  double average_frequency_high_dual_steepest_edge_weight;
  double average_log_low_dual_steepest_edge_weight_error;
  double average_log_high_dual_steepest_edge_weight_error;
  double max_average_frequency_low_dual_steepest_edge_weight;
  double max_average_frequency_high_dual_steepest_edge_weight;
  double max_sum_average_frequency_extreme_dual_steepest_edge_weight;
  double max_average_log_low_dual_steepest_edge_weight_error;
  double max_average_log_high_dual_steepest_edge_weight_error;
  double max_sum_average_log_extreme_dual_steepest_edge_weight_error;

  HighsInt num_invert_report_since_last_header;
  HighsInt num_iteration_report_since_last_header;
  double highs_run_time;
  double last_user_log_time;
  double delta_user_log_time;
  bool timeless_log;

  double average_concurrency;
  double average_fraction_of_possible_minor_iterations_performed;
  HighsInt sum_multi_chosen;
  HighsInt sum_multi_finished;

  // Analysis of INVERT form
  HighsInt num_invert;
  HighsInt num_kernel;
  HighsInt num_major_kernel;
  double max_kernel_dim;
  double sum_kernel_dim;
  double running_average_kernel_dim;
  double sum_invert_fill_factor;
  double sum_kernel_fill_factor;
  double sum_major_kernel_fill_factor;
  double running_average_invert_fill_factor;
  double running_average_kernel_fill_factor;
  double running_average_major_kernel_fill_factor;

  HighsInt AnIterIt0;
  HighsInt AnIterPrevIt;

  // Major operation analysis struct
  struct AnIterOpRec {
    double AnIterOpHyperCANCEL;
    double AnIterOpHyperTRAN;
    HighsInt AnIterOpRsDim;
    HighsInt AnIterOpNumCa;
    HighsInt AnIterOpNumHyperOp;
    HighsInt AnIterOpNumHyperRs;
    double AnIterOpSumLog10RsDensity;
    HighsInt AnIterOpRsMxNNZ;
    std::string AnIterOpName;
    HighsValueDistribution AnIterOp_density;
  };
  AnIterOpRec AnIterOp[kNumSimplexNlaOperation];

  struct AnIterTraceRec {
    double AnIterTraceTime;
    double AnIterTraceMulti;
    double AnIterTraceDensity[kNumSimplexNlaOperation];
    double AnIterTraceCostlyDse;
    HighsInt AnIterTraceIter;
    HighsInt AnIterTrace_simplex_strategy;
    HighsInt AnIterTrace_edge_weight_mode;
  };

  HighsInt AnIterTraceNumRec;
  HighsInt AnIterTraceIterDl;
  AnIterTraceRec AnIterTrace[1 + kAnIterTraceMaxNumRec + 1];

  HighsInt AnIterNumInvert[kRebuildReasonCount];
  HighsInt AnIterNumEdWtIt[(HighsInt)EdgeWeightMode::kCount];

  HighsValueDistribution primal_step_distribution;
  HighsValueDistribution dual_step_distribution;
  HighsValueDistribution simplex_pivot_distribution;
  HighsValueDistribution numerical_trouble_distribution;
  HighsValueDistribution factor_pivot_threshold_distribution;
  HighsValueDistribution edge_weight_error_distribution;
};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
