/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/KktChStep.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef TEST_KKTCH2_H_
#define TEST_KKTCH2_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "test/DevKkt.h"

namespace presolve {

namespace dev_kkt_check {

class KktCheck;

class KktChStep {
 public:
  KktChStep() {}
  virtual ~KktChStep() {}

  std::vector<double> RcolCost;
  std::vector<double> RcolLower;
  std::vector<double> RcolUpper;
  std::vector<double> RrowLower;
  std::vector<double> RrowUpper;

  int print = 1;

  std::stack<std::vector<std::pair<int, double> > > rLowers;
  std::stack<std::vector<std::pair<int, double> > > rUppers;
  std::stack<std::vector<std::pair<int, double> > > cLowers;
  std::stack<std::vector<std::pair<int, double> > > cUppers;
  std::stack<std::vector<std::pair<int, double> > > costs;

  // full matrix
  void setBoundsCostRHS(const std::vector<double>& colUpper_,
                        const std::vector<double>& colLower_,
                        const std::vector<double>& cost,
                        const std::vector<double>& rowLower_,
                        const std::vector<double>& rowUpper_);
  void addChange(int type, int row, int col, double valC, double dualC,
                 double dualR);
  void addCost(int col, double value);

  dev_kkt_check::State initState(
      const int numCol_, const int numRow_, const std::vector<int>& Astart_,
      const std::vector<int>& Aend_, const std::vector<int>& Aindex_,
      const std::vector<double>& Avalue_, const std::vector<int>& ARstart_,
      const std::vector<int>& ARindex_, const std::vector<double>& ARvalue_,
      const std::vector<int>& flagCol_, const std::vector<int>& flagRow_,
      const std::vector<double>& colValue_, const std::vector<double>& colDual_,
      const std::vector<double>& rowValue_, const std::vector<double>& rowDual_,
      const std::vector<HighsBasisStatus>& col_status_,
      const std::vector<HighsBasisStatus>& row_status_);
};

}  // namespace dev_kkt_check

}  // namespace presolve
#endif /* TEST_KKTCHSTEP_H_ */
