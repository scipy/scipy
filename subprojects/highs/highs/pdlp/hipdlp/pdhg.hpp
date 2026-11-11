/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/pdhg.hpp
 * @brief
 */
#ifndef PDHG_HPP
#define PDHG_HPP

#include "HConfig.h"
#if defined(CUPDLP_GPU) || defined(HIPDLP_GPU)
#include "gpu_backend.hpp"
#endif

#define PDLP_PROFILE (0)

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "linalg.hpp"
#include "logger.hpp"
#if PDLP_PROFILE
#include "pdlp/HiPdlpTimer.h"
#endif
#include "restart.hpp"
#include "scaling.hpp"
#include "solver_results.hpp"

const bool kUseCupdlpx = true;
const bool kTempSetting = true;

// GPU error-checking macros are defined in gpu_backend.hpp as
// GPU_CHECK / GPU_BLAS_CHECK / GPU_SPARSE_CHECK.

// Forward declarations
struct StepSizeConfig;

class PDLPSolver {
 public:
  // --- setup & Main Interface ---
  HighsStatus setup(const HighsOptions& options, HighsTimer& timer);
  void passLp(const HighsLp* lp) { original_lp_ = lp; }
  void preprocessLp();
  void scaleProblem();

  // Main entry point
  void solve(std::vector<double>& x, std::vector<double>& y);

  void unscaleSolution(std::vector<double>& x, std::vector<double>& y);
  PostSolveRetcode postprocess(HighsSolution& solution);
  void logSummary();

  // --- Getters ---
  TerminationStatus getTerminationCode() const { return results_.term_code; }
  HighsInt getIterationCount() const { return final_iter_count_; }
  HighsInt getNumCol() const { return lp_.num_col_; }
  HighsInt getNumRow() const { return lp_.num_row_; }

#if PDLP_DEBUG_LOG
  FILE* debug_pdlp_log_file_ = nullptr;
  DebugPdlpData debug_pdlp_data_;
  void closeDebugLog();
#endif
#if PDLP_PROFILE
  void reportHipdlpTimer();
#endif

 private:
  // --- Core Algorithm Logic ---
  void initialize();
  void solveReturn(const TerminationStatus term_code);

  // Returns true if solver should terminate, false if it should continue
  // Updates 'status' if returning true.
  bool runConvergenceCheck(size_t iter, std::vector<double>& output_x,
                           std::vector<double>& output_y,
                           TerminationStatus& status);
  bool checkRestartCriteria(double current_fpe, double initial_fpe,
                            double last_trial_fpe, int halpern_iteration,
                            int final_iter_count);
  void performPdhgStep();
  void performHalpernStep();
  void accumulateAverages(size_t iter);
  void prepareNextIteration();  // Swaps pointers/vectors

  // --- Convergence & Math Helpers ---
  double powerMethod();
#if defined(CUPDLP_GPU) || defined(HIPDLP_GPU)
  double powerMethodGpu();
#endif
  void initializeStepSizes();

  // Convergence checks
  bool checkConvergence(const HighsInt iter, const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& ax_vector,
                        const std::vector<double>& aty_vector, double epsilon,
                        SolverResults& results, const char* type,
                        std::vector<double>& dSlackPos,
                        std::vector<double>& dSlackNeg);

  void computeAverageIterate(std::vector<double>& ax_avg,
                             std::vector<double>& aty_avg);

  // Step updates
  std::vector<double> updateX(const std::vector<double>& x,
                              const std::vector<double>& aty,
                              double primal_step);
  std::vector<double> updateY(const std::vector<double>& y,
                              const std::vector<double>& ax,
                              const std::vector<double>& ax_next,
                              double dual_step);
  void updateIteratesFixed();
  void updateIteratesAdaptive();
  bool updateIteratesMalitskyPock(bool first_malitsky_iteration);

  // Restart helpers
  void computeStepSizeRatio(PrimalDualParams& working_params);
  void applyHalpernAveraging(std::vector<double>& x, std::vector<double>& y,
                             std::vector<double>& ax, std::vector<double>& aty);
  void updateAverageIterates(const std::vector<double>& x,
                             const std::vector<double>& y, HighsInt inner_iter);
  void performHalpernPdhgStep(bool is_major, int k_offset);
  void updatePrimalWeightAtRestart(const SolverResults& results);

  // Feasibility calculations
  double computePrimalFeasibility(const std::vector<double>& Ax_vector);
  double computeDualFeasibility(const std::vector<double>& ATy_vector,
                                std::vector<double>& dSlackPos,
                                std::vector<double>& dSlackNeg);
  void computeDualSlacks(const std::vector<double>& dualResidual,
                         std::vector<double>& dSlackPos,
                         std::vector<double>& dSlackNeg);
  double computeDualObjective(const std::vector<double>& y,
                              const std::vector<double>& dSlackPos,
                              const std::vector<double>& dSlackNeg);
  std::vector<double> computeLambda(const std::vector<double>& y,
                                    const std::vector<double>& ATy_vector);
  std::tuple<double, double, double, double, double> computeDualityGap(
      const std::vector<double>& x, const std::vector<double>& y,
      const std::vector<double>& lambda);

  // Other utilities
  void printConstraintInfo();
  bool checkNumericalStability(const std::vector<double>& delta_x,
                               const std::vector<double>& delta_y);
  double computeMovement(const std::vector<double>& delta_primal,
                         const std::vector<double>& delta_dual);
  double computeNonlinearity(const std::vector<double>& delta_primal,
                             const std::vector<double>& delta_aty);
  double computeFixedPointError();

  // --- Data Members ---
  HighsLp lp_;
  const HighsLp* original_lp_ = nullptr;
  HighsLp unscaled_processed_lp_;

  PrimalDualParams params_;
  StepSizeConfig stepsize_;
  Logger logger_;
  HighsLogOptions log_options_;
  SolverResults results_;

  Scaling scaling_;
  RestartScheme restart_scheme_;

  // Problem dimensions & metadata
  HighsInt original_num_col_ = 0;
  HighsInt num_eq_rows_ = 0;
  HighsInt sense_origin_ = 1;
  double unscaled_rhs_norm_ = 0.0;
  double unscaled_c_norm_ = 0.0;
  std::vector<HighsBool> is_equality_row_;
  std::vector<HighsInt> constraint_new_idx_;
  std::vector<ConstraintType> constraint_types_;

  // Iteration State Vectors
  std::vector<double> x_current_, y_current_;
  std::vector<double> x_next_, y_next_;
  std::vector<double> x_avg_, y_avg_;
  std::vector<double> x_sum_, y_sum_;
  std::vector<double> x_at_last_restart_, y_at_last_restart_;
  std::vector<double> reflected_x_, reflected_y_;  // For over-relaxed Halpern
  std::vector<double> x_anchor_, y_anchor_;        // For Halpern

  // Caches
  std::vector<double> Ax_cache_, ATy_cache_;
  std::vector<double> Ax_next_, ATy_next_;
  // Average Caches
  std::vector<double> Ax_avg_, ATy_avg_;
  // Residual Caches
  std::vector<double> K_times_x_diff_;
  std::vector<double> dSlackPos_, dSlackNeg_;
  std::vector<double> dSlackPosAvg_, dSlackNegAvg_;
  std::vector<double> halpern_dual_slack_next_;
  bool halpern_dual_slack_next_valid_ = false;
  double initial_fpe_ = 0.0;  // For Halpern restart
  double fpe_ = 0.0;

  // Scalars
  HighsInt final_iter_count_ = 0;
  HighsInt num_rejected_steps_ = 0;
  HighsInt halpern_iteration_ = 0;
  double sum_weights_ = 0.0;
  double current_eta_ = 0.0;
  double ratio_last_two_step_sizes_ = 1.0;

  // PID
  double primal_weight_ = 1.0;
  double best_primal_weight_ = 1.0;
  double primal_weight_error_sum_ = 0.0;
  double primal_weight_last_error_ = 0.0;
  double best_primal_dual_residual_gap_ =
      std::numeric_limits<double>::infinity();

  HighsTimer* highs_timer_p_;

#if PDLP_PROFILE
  HipdlpTimer hipdlp_timer_;
  HighsTimerClock hipdlp_clocks_;
  void hipdlpTimerStart(const HighsInt hipdlp_clock);
  void hipdlpTimerStop(const HighsInt hipdlp_clock);
#endif

#if defined(CUPDLP_GPU) || defined(HIPDLP_GPU)
  HighsInt a_num_rows_ = 0;
  HighsInt a_num_cols_ = 0;
  HighsInt a_nnz_ = 0;
  double sum_weights_gpu_ = 0.0;
  // --- GPU Specifics ---
  gpuStream_t gpu_stream_ = nullptr;
  gpuSparseHandle_t cusparse_handle_ = nullptr;
  gpuBlasHandle_t cublas_handle_ = nullptr;

  // Matrix descriptors
  gpuSpMatDescr_t mat_a_csr_ = nullptr;
  gpuSpMatDescr_t mat_a_T_csr_ = nullptr;

  // Device Pointers
  HighsInt *d_a_row_ptr_ = nullptr, *d_a_col_ind_ = nullptr;
  double* d_a_val_ = nullptr;
  HighsInt *d_at_row_ptr_ = nullptr, *d_at_col_ind_ = nullptr;
  double* d_at_val_ = nullptr;
  // Packed into a single device allocation to avoid 3 separate H→D transfers
  // per iteration. Layout: [primal_step, dual_step, halpern_iteration (as
  // int)].
  struct GpuStepParams {
    double primal_step;
    double dual_step;
    int halpern_iteration;
  };
  GpuStepParams* d_step_params_ = nullptr;
  // Convenience device pointers aliasing into d_step_params_
  int* d_halpern_iteration_ = nullptr;
  double* d_primal_step_size_ = nullptr;
  double* d_dual_step_size_ = nullptr;

  // GPU Vectors (Device memory)
  double *d_x_current_ = nullptr, *d_y_current_ = nullptr;
  double *d_x_next_ = nullptr, *d_y_next_ = nullptr;
  double *d_x_avg_ = nullptr, *d_y_avg_ = nullptr;
  double *d_x_sum_ = nullptr, *d_y_sum_ = nullptr;
  double *d_ax_current_ = nullptr, *d_aty_current_ = nullptr;
  double *d_ax_next_ = nullptr, *d_aty_next_ = nullptr;
  double *d_ax_avg_ = nullptr, *d_aty_avg_ = nullptr;

  double *d_x_at_last_restart_ = nullptr, *d_y_at_last_restart_ = nullptr;
  double *d_x_anchor_ = nullptr, *d_y_anchor_ = nullptr;
  double *d_pdhg_primal_ = nullptr, *d_pdhg_dual_ = nullptr;
  double* d_delta_x_ = nullptr;
  double* d_delta_y_ = nullptr;
  double* d_AT_delta_y_ = nullptr;

  double *d_col_cost_ = nullptr, *d_col_lower_ = nullptr,
         *d_col_upper_ = nullptr;
  double *d_row_lower_ = nullptr, *d_col_scale_ = nullptr,
         *d_row_scale_ = nullptr;
  bool* d_is_equality_row_ = nullptr;

  // Scratch / Output buffers
  double* d_convergence_results_ = nullptr;
  double *d_dSlackPos_ = nullptr, *d_dSlackNeg_ = nullptr;
  double *d_dSlackPosAvg_ = nullptr, *d_dSlackNegAvg_ = nullptr;
  double *d_buffer_ = nullptr, *d_buffer2_ = nullptr;  // General purpose
  double *d_x_temp_diff_norm_result_ = nullptr,
         *d_y_temp_diff_norm_result_ = nullptr;

  // SpMV buffers
  void* d_spmv_buffer_ax_ = nullptr;
  size_t spmv_buffer_size_ax_ = 0;
  void* d_spmv_buffer_aty_ = nullptr;
  size_t spmv_buffer_size_aty_ = 0;

  // Vector Descriptors
  gpuDnVecDescr_t vec_x_desc_ = nullptr, vec_y_desc_ = nullptr;
  gpuDnVecDescr_t vec_ax_desc_ = nullptr, vec_aty_desc_ = nullptr;

  // GPU Methods
  void setupGpu();
  void cleanupGpu();
  void linalgGpuAx(const double* d_x_in, double* d_ax_out);
  void linalgGpuATy(const double* d_y_in, double* d_aty_out);
  bool checkConvergenceGpu(const HighsInt iter, const double* d_x,
                           const double* d_y, const double* d_ax,
                           const double* d_aty, double epsilon,
                           SolverResults& results, const char* type,
                           double* d_slackPos_out, double* d_slackNeg_out,
                           bool use_halpern_slack);
  void computeStepSizeRatioGpu(PrimalDualParams& working_params);
  void updateAverageIteratesGpu(HighsInt inner_iter);
  void computeAverageIterateGpu();
  double computeMovementGpu(const double* d_x_new, const double* d_x_old,
                            const double* d_y_new, const double* d_y_old);
  double computeNonlinearityGpu(const double* d_x_new, const double* d_x_old,
                                const double* d_aty_new,
                                const double* d_aty_old);
  double computeDiffNormCuBLAS(const double* d_a, const double* d_b, int n);
  void performHalpernPdhgStepGpu(bool is_major, int k_offset);
  double computeFixedPointErrorGpu();

  // Kernel Wrappers
  void launchKernelUpdateX(double primal_step);
  void launchKernelUpdateY(double dual_step);
  void launchKernelUpdateAverages(double weight);
  void launchKernelScaleVector(double* d_out, const double* d_in, double scale,
                               HighsInt n);
#endif
};

#endif  // PDHG_HPP
