/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsComponent.h
 * @brief The HiGHS class
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHS_COMPONENT_H_
#define UTIL_HIGHS_COMPONENT_H_

#include "lp_data/HighsOptions.h"

// HighsComponentData is a placeholder for structs which we will keep after
// run() is done, internally.
struct HighsComponentData {
  bool is_valid = false;
};

// HighsComponentInfo is a placeholder for details we want to query from outside
// of HiGHS like execution information.
struct HighsComponentInfo {
  bool is_valid = false;
};

// HighsComponentOptions is a placeholder for options specific to this component
struct HighsComponentOptions {
  bool is_valid = false;
};

class HighsComponent {
 public:
  virtual void clear() = 0;
  HighsStatus run();
  HighsStatus setOptions(const HighsOptions& options);

  const HighsComponentInfo& getInfo() { return info_; }
  const HighsComponentData& getData() { return data_; }
  const HighsComponentOptions& getOptions() { return options_; }

  virtual ~HighsComponent() = default;

 private:
  bool has_run_ = false;

  HighsComponentInfo info_;
  HighsComponentData data_;
  HighsComponentOptions options_;
};

#endif
