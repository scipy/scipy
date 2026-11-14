/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkControl.cpp
 * @brief
 */

#include "simplex/HEkk.h"

void HEkk::initialiseControl() {
  // Copy tolerances from options
  info_.allow_dual_steepest_edge_to_devex_switch =
      options_->simplex_dual_edge_weight_strategy ==
      kSimplexEdgeWeightStrategyChoose;
  info_.dual_steepest_edge_weight_log_error_threshold =
      options_->dual_steepest_edge_weight_log_error_threshold;
  info_.control_iteration_count0 = iteration_count_;
  // Initialise the densities
  info_.col_aq_density = 0;
  info_.row_ep_density = 0;
  info_.row_ap_density = 0;
  info_.row_DSE_density = 0;
  info_.col_steepest_edge_density = 0;
  info_.col_basic_feasibility_change_density = 0;
  info_.row_basic_feasibility_change_density = 0;
  info_.col_BFRT_density = 0;
  info_.primal_col_density = 0;
  // Set the row_dual_density to 1 since it's assumed all costs are at
  // least perturbed from zero, if not initially nonzero
  info_.dual_col_density = 1;
  // Initialise the data used to determine the switch from DSE to
  // Devex
  info_.costly_DSE_frequency = 0;
  info_.num_costly_DSE_iteration = 0;
  info_.costly_DSE_measure = 0;
  info_.average_log_low_DSE_weight_error = 0;
  info_.average_log_high_DSE_weight_error = 0;
}

void HEkk::updateOperationResultDensity(const double local_density,
                                        double& density) const {
  density = (1 - kRunningAverageMultiplier) * density +
            kRunningAverageMultiplier * local_density;
}

void HEkk::assessDSEWeightError(const double computed_edge_weight,
                                const double updated_edge_weight) {
  // Compute the (relative) dual steepest edge weight error for
  // analysis and debugging
  edge_weight_error_ = std::fabs(updated_edge_weight - computed_edge_weight) /
                       std::max(1.0, computed_edge_weight);
  if (edge_weight_error_ > options_->dual_steepest_edge_weight_error_tolerance)
    highsLogDev(options_->log_options, HighsLogType::kInfo,
                "Dual steepest edge weight error is %g\n", edge_weight_error_);
  // Compute the relative deviation in the updated weight compared
  // with the computed weight
  double weight_relative_deviation;
  if (updated_edge_weight < computed_edge_weight) {
    // Updated weight is low
    weight_relative_deviation = computed_edge_weight / updated_edge_weight;
    info_.average_log_low_DSE_weight_error =
        0.99 * info_.average_log_low_DSE_weight_error +
        0.01 * log(weight_relative_deviation);
  } else {
    // Updated weight is correct or high
    weight_relative_deviation = updated_edge_weight / computed_edge_weight;
    info_.average_log_high_DSE_weight_error =
        0.99 * info_.average_log_high_DSE_weight_error +
        0.01 * log(weight_relative_deviation);
  }
}

bool HEkk::switchToDevex() {
  // Parameters controlling switch from DSE to Devex on cost
  const double kCostlyDseMeasureLimit = 1000.0;
  const double kCostlyDseMinimumDensity = 0.01;
  const double kCostlyDseFractionNumTotalIterationBeforeSwitch = 0.1;
  const double kCostlyDseFractionNumCostlyDseIterationBeforeSwitch = 0.05;
  bool switch_to_devex = false;
  // Firstly consider switching on the basis of NLA cost
  double costly_DSE_measure_denominator;
  costly_DSE_measure_denominator = max(
      max(info_.row_ep_density, info_.col_aq_density), info_.row_ap_density);
  if (costly_DSE_measure_denominator > 0) {
    info_.costly_DSE_measure =
        info_.row_DSE_density / costly_DSE_measure_denominator;
    info_.costly_DSE_measure =
        info_.costly_DSE_measure * info_.costly_DSE_measure;
  } else {
    info_.costly_DSE_measure = 0;
  }
  bool costly_DSE_iteration =
      info_.costly_DSE_measure > kCostlyDseMeasureLimit &&
      info_.row_DSE_density > kCostlyDseMinimumDensity;
  info_.costly_DSE_frequency =
      (1 - kRunningAverageMultiplier) * info_.costly_DSE_frequency;
  if (costly_DSE_iteration) {
    info_.num_costly_DSE_iteration++;
    info_.costly_DSE_frequency += kRunningAverageMultiplier * 1.0;
    // What if non-dual iterations have been performed: need to think about this
    HighsInt local_iteration_count =
        iteration_count_ - info_.control_iteration_count0;
    HighsInt local_num_tot = lp_.num_col_ + lp_.num_row_;
    // Switch to Devex if at least 5% of the (at least) 0.1NumTot iterations
    // have been costly
    switch_to_devex =
        info_.allow_dual_steepest_edge_to_devex_switch &&
        (info_.num_costly_DSE_iteration >
         local_iteration_count *
             kCostlyDseFractionNumCostlyDseIterationBeforeSwitch) &&
        (local_iteration_count >
         kCostlyDseFractionNumTotalIterationBeforeSwitch * local_num_tot);

    if (switch_to_devex) {
      highsLogDev(options_->log_options, HighsLogType::kInfo,
                  "Switch from DSE to Devex after %" HIGHSINT_FORMAT
                  " costly DSE iterations of %" HIGHSINT_FORMAT
                  " with "
                  "densities C_Aq = %11.4g; R_Ep = %11.4g; R_Ap = "
                  "%11.4g; DSE = %11.4g\n",
                  info_.num_costly_DSE_iteration, local_iteration_count,
                  info_.col_aq_density, info_.row_ep_density,
                  info_.row_ap_density, info_.row_DSE_density);
    }
  }
  if (!switch_to_devex) {
    // Secondly consider switching on the basis of weight accuracy
    double local_measure = info_.average_log_low_DSE_weight_error +
                           info_.average_log_high_DSE_weight_error;
    double local_threshold =
        info_.dual_steepest_edge_weight_log_error_threshold;
    switch_to_devex = info_.allow_dual_steepest_edge_to_devex_switch &&
                      local_measure > local_threshold;
    if (switch_to_devex) {
      highsLogDev(options_->log_options, HighsLogType::kInfo,
                  "Switch from DSE to Devex with log error measure of %g > "
                  "%g = threshold\n",
                  local_measure, local_threshold);
    }
  }
  return switch_to_devex;
}
