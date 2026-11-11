/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/defs.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_DEFS_HPP
#define PDLP_HIPDLP_DEFS_HPP

#include <cmath>
#include <vector>

#include "Highs.h"

enum class Device { CPU, GPU };

enum class ScalingMethod { NONE, RUIZ, POCK_CHAMBOLLE, L2_NORM, COMBINED };

enum class RestartStrategy {
  NO_RESTART,
  FIXED_RESTART,
  ADAPTIVE_RESTART,
  HALPERN_RESTART
};

enum class StepSizeStrategy { FIXED, ADAPTIVE, MALITSKY_POCK, PID };

enum class PostSolveRetcode {
  OK = 0,
  INVALID_SOLUTION = 1,
  SCALING_ERROR = 2,
  DIMENSION_MISMATCH = 3,
  NUMERICAL_ERROR = 4,
  CONSTRAINT_VIOLATION = 5
};

struct StepSizeConfig {
  double primal_step;
  double dual_step;
  double beta;
  double power_method_lambda;
  HighsInt step_size_iter = 0;  // nStepSizeIter
};

struct MalitskyPockParams {
  double step_size_interpolation = 0.5;  // Between 0 and 1
  double step_size_downscaling_factor = 0.7;
  double linesearch_contraction_factor = 0.99;
  void initialise();
};

struct AdaptiveLinesearchParams {
  double step_size_reduction_exponent = 0.3;
  double step_size_growth_exponent = 0.6;
  void initialise();
};

struct PrimalDualParams {
  double eta;
  double omega;
  double tolerance;
  HighsInt max_iterations;
  Device device_type;
  double time_limit = 3600.0;

  // Restart parameters
  RestartStrategy restart_strategy;
  HighsInt fixed_restart_interval;

  bool use_halpern_restart = false;
  double halpern_gamma = 1.0;  // 0: standard Halpern, 1: full reflection, in
                               // between: over relaxation

  // Scaling parameters
  bool use_ruiz_scaling = false;
  bool use_pc_scaling = false;
  bool use_l2_scaling = false;

  // Ruiz scaling parameters
  HighsInt ruiz_iterations = 10;
  double ruiz_norm = INFINITY;

  // Pock-Chambolle scaling parameters
  double pc_alpha = 1.0;

  // Step sizes strategy
  StepSizeStrategy step_size_strategy = StepSizeStrategy::FIXED;
  double k_p = 0.99;
  double k_i = 0.01;
  double k_d = 0.0;
  double i_smooth = 0.3;

  MalitskyPockParams malitsky_pock_params;
  AdaptiveLinesearchParams adaptive_linesearch_params;
  HighsLogOptions log_options_;
  void initialise();
};

struct PdlpIterate {
  // Primary variables
  std::vector<double> x;
  std::vector<double> y;

  // Cached matrix-vector products
  mutable std::vector<double> Ax;
  mutable std::vector<double> Aty;
  mutable bool Ax_valid = false;
  mutable bool Aty_valid = false;

  // Constructors
  PdlpIterate() = default;
  PdlpIterate(HighsInt num_cols, HighsInt num_rows)
      : x(num_cols, 0.0),
        y(num_rows, 0.0),
        Ax(num_rows, 0.0),
        Aty(num_cols, 0.0) {}
  PdlpIterate(const std::vector<double>& x_init,
              const std::vector<double>& y_init)
      : x(x_init), y(y_init), Ax(y_init.size(), 0.0), Aty(x_init.size(), 0.0) {}

  // Arithmetic operations
  // z = alpha * this + beta * other
  void LinearCombination(const PdlpIterate& other, double alpha, double beta);

  // this = this + alpha * other
  void AddScaled(const PdlpIterate& other, double alpha);

  // this = alpha * this
  void Scale(double alpha);

  // Copy operations
  void CopyFrom(const PdlpIterate& other);
  PdlpIterate& operator=(const PdlpIterate& other);

  // Norms and metrics
  double PrimalNorm() const;  // ||x||_2
  double DualNorm() const;    // ||y||_2
  double WeightedNorm(
      double omega) const;  // sqrt(||x||_2^2 + omega^2 * ||y||_2^2)

  // Distance metrics
  double Distance(const PdlpIterate& other, double omega = 1.0) const;

  // Matrix-vector product management
  void ComputeAx(const HighsLp& lp) const;
  void ComputeATy(const HighsLp& lp) const;
  void InvalidateProducts();  // Call when x or y change

  const std::vector<double>& GetAx(const HighsLp& lp) const;
  const std::vector<double>& GetATy(const HighsLp& lp) const;

  // For block-structured problems
  struct BlockStructure {
    std::vector<HighsInt> x_block_sizes;
    std::vector<HighsInt> y_block_sizes;
    // std::vector<std::vector<double>> x_blocks;  // Future: for block problems
    // std::vector<std::vector<double>> y_blocks;
  };

  // Optional: block structure for future extensions
  std::unique_ptr<BlockStructure> block_structure = nullptr;

 private:
  void EnsureAxComputed(const HighsLp& lp) const;
  void EnsureATyComputed(const HighsLp& lp) const;
};

namespace pdlp_iterate_ops {
// Compute z_new = z_old -
};

struct DetailedTimings {
  double total_time = 0.0;
  double iterate_update_time = 0.0;
  double matrix_multiply_time = 0.0;  // Ax and ATy
  double convergence_check_time = 0.0;
  double restart_check_time = 0.0;
  double average_iterate_time = 0.0;
  double projection_time = 0.0;
  double step_size_adjustment_time = 0.0;
  double other_time = 0.0;

  void print(const std::string& solver_name,
             const HighsLogOptions& log_options) const {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "\n=== %s Detailed Timings ===\n", solver_name.c_str());
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Total time:              %6.2f s\n", total_time);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Iterate update:          %6.2f s (%3.0f)\n",
                 iterate_update_time, iterate_update_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "  - Matrix multiply:     %6.2f s (%3.0f)\n",
                 matrix_multiply_time, matrix_multiply_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "  - Projection:          %6.2f s (%3.0f)\n", projection_time,
                 projection_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "  - Step size adjust:    %6.2f s (%3.0f)\n",
                 step_size_adjustment_time,
                 step_size_adjustment_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Convergence check:       %6.2f s (%3.0f)\n",
                 convergence_check_time,
                 convergence_check_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Restart check:           %6.2f s (%3.0f)\n",
                 restart_check_time, restart_check_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Average iterate comp:    %6.2f s (%3.0f)\n",
                 average_iterate_time, average_iterate_time / total_time * 100);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Other:                   %6.2f s (%3.0f)\n", other_time,
                 other_time / total_time * 100);
  }
};

#endif
