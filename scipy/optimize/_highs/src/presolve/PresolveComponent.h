/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file PresolveComponent.h
 * @brief The HiGHS class
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_COMPONENT_H_
#define PRESOLVE_PRESOLVE_COMPONENT_H_

#include "presolve/Presolve.h"
#include "util/HighsComponent.h"

// Class defining the Presolve Component to be used in HiGHS.
// What used to be in Presolve.h but allowing for further testing and dev.

// The structure of component is general, of the presolve component - presolve
// specific.

struct PresolveComponentData : public HighsComponentData {
  std::vector<presolve::Presolve> presolve_;
  HighsLp reduced_lp_;

  // todo: make reduced one const.
  HighsSolution reduced_solution_;
  HighsSolution recovered_solution_;
};

// HighsComponentInfo is a placeholder for details we want to query from outside
// of HiGHS like execution information.
struct PresolveComponentInfo : public HighsComponentInfo {
  int n_rows_removed = 0;
  int n_cols_removed = 0;
  int n_nnz_removed = 0;

  double presolve_time = 0;
  double postsolve_time = 0;

  double solve_time = 0;
  double cleanup_time = 0;
};

// HighsComponentOptions is a placeholder for options specific to this component
struct PresolveComponentOptions : public HighsComponentOptions {
  bool is_valid = false;
  // presolve options later when needed.
  bool presolve_on = true;
  std::vector<presolve::Presolver> order;
};

class PresolveComponent : public HighsComponent {
 public:
  void clear() override;

  HighsStatus init(const HighsLp& lp, HighsTimer& timer);

  HighsPresolveStatus run();

  HighsLp& getReducedProblem() { return data_.reduced_lp_; }

  HighsStatus setOptions(const HighsOptions& options);

  void setBasisInfo(const std::vector<HighsBasisStatus>& pass_col_status,
                    const std::vector<HighsBasisStatus>& pass_row_status);

  void negateReducedLpColDuals(bool reduced);
  void negateReducedLpCost();

  bool has_run_ = false;

  PresolveComponentInfo info_;
  PresolveComponentData data_;
  PresolveComponentOptions options_;

  HighsPresolveStatus presolve_status_ = HighsPresolveStatus::NotPresolved;
  HighsPostsolveStatus postsolve_status_ = HighsPostsolveStatus::NotPresolved;
};

#endif