/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_LP_H_
#define LP_DATA_HIGHS_LP_H_

#include <string>
#include <vector>

#include "lp_data/HConst.h"

class HighsLp;

class HighsLp {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  ObjSense sense_ = ObjSense::MINIMIZE;
  double offset_ = 0;

  std::string model_name_ = "";
  std::string lp_name_ = "";

  std::vector<std::string> col_names_;
  std::vector<std::string> row_names_;

  std::vector<int> integrality_;

  bool operator==(const HighsLp& lp);
  bool equalButForNames(const HighsLp& lp);
  void clear();
};

#endif
