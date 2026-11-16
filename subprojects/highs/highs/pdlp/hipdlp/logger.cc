/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/logger.cc
 * @brief
 */
#include "logger.hpp"

#include <iostream>
#include <map>
#include <sstream>

// Helper to convert enums to strings for pretty printing
template <typename T>
std::string EnumToString(T e, const std::map<T, std::string>& map) {
  auto it = map.find(e);
  return it != map.end() ? it->second : "UNKNOWN";
}

// Logger implementation
void Logger::initialise(const HighsInt log_dev_level,
                        const HighsLogOptions& log_options,
                        HighsTimer* highs_timer_p) {
  if (log_dev_level == kHighsLogDevLevelVerbose) {
    console_level_ = LogLevel::kDebug;
  } else if (log_dev_level == kHighsLogDevLevelDetailed) {
    console_level_ = LogLevel::kVerbose;
  } else if (log_dev_level == kHighsLogDevLevelInfo) {
    console_level_ = LogLevel::kDetailed;
  } else {
    console_level_ = LogLevel::kInfo;  // None;
  }
  log_options_ = log_options;
  highs_timer_p_ = highs_timer_p;
  iteration_stats_count_ = kHighsIInf;
  iteration_stats_time_ = -kHighsInf;
}

void Logger::log(LogLevel level, const std::string& message) const {
  // Now using HiGHS IO
  if (level <= console_level_)
    highsLogUser(log_options_, HighsLogType::kInfo, "%s\n", message.c_str());
}

void Logger::info(const std::string& message) const {
  log(LogLevel::kInfo, message);
}
void Logger::detailed(const std::string& message) const {
  log(LogLevel::kDetailed, message);
}
void Logger::verbose(const std::string& message) const {
  log(LogLevel::kVerbose, message);
}
void Logger::debug(const std::string& message) const {
  log(LogLevel::kDebug, message);
}

void Logger::printHeader() const {
  detailed("Using HiPDLP with developer logging");
}

void Logger::printParams(const PrimalDualParams& params) const {
  detailed("\nSolver Parameters:");
  std::stringstream ss;

  std::map<RestartStrategy, std::string> restart_map = {
      {RestartStrategy::NO_RESTART, "None"},
      {RestartStrategy::FIXED_RESTART, "Fixed"},
      {RestartStrategy::ADAPTIVE_RESTART, "Adaptive"}};
  std::map<StepSizeStrategy, std::string> step_size_map = {
      {StepSizeStrategy::FIXED, "Fixed"},
      {StepSizeStrategy::ADAPTIVE, "Adaptive"},
      {StepSizeStrategy::MALITSKY_POCK, "Malitsky-Pock"},
      {StepSizeStrategy::PID, "PID"}};

  ss << "  - Max Iterations: " << params.max_iterations;
  detailed(ss.str());
  ss.str("");
  ss << "  - Tolerance: " << params.tolerance;
  detailed(ss.str());
  ss.str("");
  ss << "  - Restart Strategy: "
     << EnumToString(params.restart_strategy, restart_map);
  detailed(ss.str());
  ss.str("");

  std::string scaling_method = "None";
  if (params.use_ruiz_scaling && params.use_pc_scaling &&
      params.use_l2_scaling) {
    scaling_method = "Combined (Ruiz + PC + L2)";
  } else if (params.use_ruiz_scaling && params.use_pc_scaling) {
    scaling_method = "Ruiz + Pock-Chambolle";
  } else if (params.use_ruiz_scaling && params.use_l2_scaling) {
    scaling_method = "Ruiz + L2-Norm";
  } else if (params.use_pc_scaling && params.use_l2_scaling) {
    scaling_method = "Pock-Chambolle + L2-Norm";
  } else if (params.use_ruiz_scaling) {
    scaling_method = "Ruiz";
  } else if (params.use_pc_scaling) {
    scaling_method = "Pock-Chambolle";
  } else if (params.use_l2_scaling) {
    scaling_method = "L2-Norm";
  }

  ss << "  - Scaling Method: " << scaling_method;
  detailed(ss.str());
  ss.str("");

  // Optionally, add more details about scaling parameters if any scaling is
  // used
  if (params.use_ruiz_scaling) {
    ss.str("");
    ss << "   * Ruiz iterations: " << params.ruiz_iterations;
    detailed(ss.str());
  }
  if (params.use_pc_scaling) {
    ss.str("");
    ss << "   * PC alpha: " << params.pc_alpha;
    detailed(ss.str());
  }
  ss << "  - Step Size Strategy: "
     << EnumToString(params.step_size_strategy, step_size_map);
  detailed(ss.str());
  detailed("------------------------------------------------------------");
}

void Logger::printIterationHeader() {
  info(
      "     Iter       primal obj         dual obj       pinf       dinf       "
      "gap     pr wt    time");
  this->iteration_stats_count_ = 0;
}

void Logger::printIterationStats(const HighsInt iter,
                                 const SolverResults& results,
                                 const double step_size, const bool forced) {
  // Repeat the header if sufficient iterations have been performed
  if (this->iteration_stats_count_ > kHipdlpIterationStatsHeaderFrequency)
    this->printIterationHeader();
  // Determine whether to log iterations
  bool iteration_log = console_level_ >= LogLevel::kDetailed || forced;
  double time_now = highs_timer_p_->read();
  iteration_log =
      iteration_log ||
      time_now > this->iteration_stats_time_ + kHipdlpIterationStatsFrequency;
  if (!iteration_log) return;

  std::stringstream ss;
  // Using functions to print using streams, taken from HiPO (and
  // originally IPX)
  ss << integer(iter, 9);
  ss << " " << sci(results.primal_obj, 16, 8);
  ss << " " << sci(results.dual_obj, 16, 8);
  ss << " " << sci(results.primal_feasibility, 10, 2);
  ss << " " << sci(results.dual_feasibility, 10, 2);
  ss << " " << sci(results.duality_gap, 9, 2);
  ss << " " << sci(step_size, 9, 2);
  ss << " " << fix(time_now, 7, 1);
  info(ss.str());
  this->iteration_stats_time_ = time_now;
  this->iteration_stats_count_++;
}

void Logger::printSummary(const SolverResults& results, HighsInt total_iter,
                          double total_time) const {
  detailed("\n-------------------- Solver Summary --------------------");
  std::stringstream ss;

  std::map<TerminationStatus, std::string> term_map = {
      {TerminationStatus::NOTSET, "Not set"},
      {TerminationStatus::OPTIMAL, "Optimal"},
      {TerminationStatus::INFEASIBLE, "Infeasible"},
      {TerminationStatus::UNBOUNDED, "Unbounded"},
      {TerminationStatus::MAXITER, "Maxiter"},
      {TerminationStatus::TIMEOUT, "Timeout"},
      {TerminationStatus::ERROR, "Error"}};

  ss << "  - Termination Status: " << EnumToString(results.term_code, term_map);
  detailed(ss.str());
  ss.str("");
  ss << "  - Total Iterations: " << total_iter;
  detailed(ss.str());
  ss.str("");
  ss << "  - Total Time: " << std::fixed << std::setprecision(3) << total_time
     << " seconds";
  detailed(ss.str());
  ss.str("");
  ss << "  - Final Primal Objective: " << std::scientific
     << std::setprecision(6) << results.primal_obj;
  detailed(ss.str());
  ss.str("");
  ss << "  - Final Duality Gap: " << std::scientific << std::setprecision(6)
     << results.duality_gap;
  detailed(ss.str());
  ss.str("");
  ss << "  - Final Primal Feasibility: " << std::scientific
     << std::setprecision(6) << results.primal_feasibility;
  detailed(ss.str());
  ss.str("");
  ss << "  - Final Dual Feasibility: " << std::scientific
     << std::setprecision(6) << results.dual_feasibility;
  detailed(ss.str());
  detailed("------------------------------------------------------------");
}

std::string format(const double d, HighsInt width, HighsInt prec,
                   std::ios_base::fmtflags floatfield) {
  std::ostringstream s;
  s.precision(prec);
  s.width(width);
  s.setf(floatfield, std::ios_base::floatfield);
  s << d;
  return s.str();
}

std::string integer(const HighsInt i, HighsInt width) {
  std::ostringstream s;
  s.width(width);
  s << i;
  return s.str();
}
