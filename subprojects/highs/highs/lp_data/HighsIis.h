/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsIis.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSIIS_H_
#define LP_DATA_HIGHSIIS_H_

#include "model/HighsModel.h"

const bool kIisDevReport = false;

enum IisBoundStatus : int {
  kIisBoundStatusDropped = -1,
  kIisBoundStatusNull,   // 0
  kIisBoundStatusFree,   // 1
  kIisBoundStatusLower,  // 2
  kIisBoundStatusUpper,  // 3
  kIisBoundStatusBoxed   // 4
};

enum IisModelStatus {
  kIisModelStatusFeasible = -1,
  kIisModelStatusUnknown,     // 0
  kIisModelStatusTimeLimit,   // 1
  kIisModelStatusReducible,   // 2
  kIisModelStatusIrreducible  // 3
};

struct HighsIisInfo {
  HighsInt num_lp_solved = 0;
  HighsInt sum_simplex_iteration_counts = 0;
  HighsInt min_simplex_iteration_count = 0;
  HighsInt max_simplex_iteration_count = 0;
  double sum_simplex_times = 0.0;
  double min_simplex_time = 0.0;
  double max_simplex_time = 0.0;
  // IIS logging state (persist across calls)
  double iis_last_disptime = -kHighsInf;
  HighsInt iis_num_disp_lines = 0;

  void clear() {
    num_lp_solved = 0;
    sum_simplex_iteration_counts = 0;
    min_simplex_iteration_count = 0;
    max_simplex_iteration_count = 0;
    sum_simplex_times = 0.0;
    min_simplex_time = 0.0;
    max_simplex_time = 0.0;
    // Reset IIS logging state
    iis_last_disptime = -kHighsInf;
    iis_num_disp_lines = 0;
  }

  void update(const double simplex_time, const HighsInt simplex_iterations) {
    if (num_lp_solved == 0) {
      min_simplex_iteration_count = kHighsIInf;
      min_simplex_time = kHighsInf;
    }
    num_lp_solved += 1;
    sum_simplex_times += simplex_time;
    min_simplex_time = std::min(simplex_time, min_simplex_time);
    max_simplex_time = std::max(simplex_time, max_simplex_time);
    sum_simplex_iteration_counts += simplex_iterations;
    min_simplex_iteration_count =
        std::min(simplex_iterations, min_simplex_iteration_count);
    max_simplex_iteration_count =
        std::max(simplex_iterations, max_simplex_iteration_count);
  }
};

class HighsIis {
 public:
  HighsIis() {}

  void clearData();
  void clear();
  void clearLogInfo();
  void invalid(const HighsLp& lp);
  std::string iisBoundStatusToString(HighsInt bound_status) const;
  std::string iisModelStatusToString(HighsInt model_status) const;
  void report(const std::string& message, const HighsLp& lp) const;
  void reportIteration(const HighsOptions& options, const HighsInt iter,
                       const HighsInt num_rows_remaining, const bool force);
  void reportFinal(const HighsOptions& options) const;
  void addCol(const HighsInt col, const HighsInt status = kIisBoundStatusNull);
  void addRow(const HighsInt row, const HighsInt status = kIisBoundStatusNull);
  void removeCol(const HighsInt col);
  void removeRow(const HighsInt row);
  HighsStatus deduce(const HighsLp& lp, const HighsOptions& options,
                     const HighsCallback& callback, const HighsBasis& basis);
  void setLp(const HighsLp& lp);
  HighsInt nonIsStatus() const;
  void setStatus(const HighsLp& lp);
  HighsInt determineBoundStatus(const double lower, const double upper,
                                const bool is_row) const;

  HighsStatus compute(const HighsLp& lp, const HighsOptions& options,
                      const HighsCallback& callback,
                      const HighsBasis* basis = nullptr);
  void processBoundRelaxation(Highs& highs, const bool row_deletion,
                              const bool drop_lower, const HighsInt iX,
                              double& lower, double& upper,
                              IisModelStatus& iis_status,
                              HighsStatus& search_return_status);

  bool trivial(const HighsLp& lp, const HighsOptions& options);
  bool rowValueBounds(const HighsLp& lp, const HighsOptions& options);

  bool indexStatusOkReturn(const bool return_value) const {
    return return_value;
  }
  bool lpDataOkReturn(const bool return_value) const { return return_value; }
  bool lpOkReturn(const bool return_value) const { return return_value; }

  bool indexStatusOk(const HighsLp& lp) const;
  bool lpDataOk(const HighsLp& lp, const HighsOptions& options) const;
  bool lpOk(const HighsOptions& options) const;

  // Data members
  bool valid_ = false;
  HighsInt status_;
  HighsInt strategy_ = kIisStrategyMin;
  std::vector<HighsInt> col_index_;
  std::vector<HighsInt> row_index_;
  std::vector<HighsInt> col_bound_;
  std::vector<HighsInt> row_bound_;
  std::vector<HighsInt> col_status_;
  std::vector<HighsInt> row_status_;
  HighsIisInfo info_;
  HighsModel model_;
};

#endif  // LP_DATA_HIGHSIIS_H_
