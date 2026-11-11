/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/logger.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_LOGGER_HPP
#define PDLP_HIPDLP_LOGGER_HPP

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "restart.hpp"  // For PrimalDualParams
#include "solver_results.hpp"

// Log verbosity level
enum class LogLevel {
  kNone,      // No output
  kInfo,      // Standard output: summary, termination, major events
  kDetailed,  // Developer output
  kVerbose,   // Detailed output: iteration-level info
  kDebug      // Verbose + debug info for developers
};

class Logger {
 public:
  void initialise(const HighsInt log_dev_level,
                  const HighsLogOptions& log_options,
                  HighsTimer* highs_timer_p);
  LogLevel getLogLevel() const { return console_level_; }
  // Logging methods for different levels
  void info(const std::string& message) const;
  void detailed(const std::string& message) const;
  void verbose(const std::string& message) const;
  void debug(const std::string& message) const;

  // Formatted printing functions
  void printHeader() const;
  void printParams(const PrimalDualParams& params) const;
  void printIterationHeader();
  void printIterationStats(const HighsInt iter,
                           const SolverResults& current_results,
                           const double current_eta, const bool forced = false);
  void printSummary(const SolverResults& results, HighsInt total_iter,
                    double total_time) const;
  void setLogOptions(const HighsLogOptions& log_options) {
    log_options_ = log_options;
  }

 private:
  void log(LogLevel level, const std::string& message) const;
  LogLevel console_level_;
  HighsLogOptions log_options_;
  HighsTimer* highs_timer_p_;
  HighsInt iteration_stats_count_;
  double iteration_stats_time_;
  const HighsInt kHipdlpIterationStatsHeaderFrequency = 50;
  const double kHipdlpIterationStatsFrequency = 5.0;
};

// Functions to print using streams, taken from HiPO (and originally
// IPX)
std::string format(const double d, HighsInt width, HighsInt prec,
                   std::ios_base::fmtflags floatfield);
std::string integer(const HighsInt i, HighsInt width = 0);
inline std::string sci(const double d, HighsInt width, HighsInt prec) {
  return format(d, width, prec, std::ios_base::scientific);
}
inline std::string fix(const double d, HighsInt width, HighsInt prec) {
  return format(d, width, prec, std::ios_base::fixed);
}

#endif  // PDLP_HIPDLP_LOGGER_HPP
