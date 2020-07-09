/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/DevKkt.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef TEST_DEV_KKT_H_
#define TEST_DEV_KKT_H_

#include <map>
#include <vector>

#include "lp_data/HConst.h"

namespace presolve {
namespace dev_kkt_check {

struct State {
  State(const int numCol_, const int numRow_, const std::vector<int>& Astart_,
        const std::vector<int>& Aend_, const std::vector<int>& Aindex_,
        const std::vector<double>& Avalue_, const std::vector<int>& ARstart_,
        const std::vector<int>& ARindex_, const std::vector<double>& ARvalue_,
        const std::vector<double>& colCost_,
        const std::vector<double>& colLower_,
        const std::vector<double>& colUpper_,
        const std::vector<double>& rowLower_,
        const std::vector<double>& rowUpper_, const std::vector<int>& flagCol_,
        const std::vector<int>& flagRow_, const std::vector<double>& colValue_,
        const std::vector<double>& colDual_,
        const std::vector<double>& rowValue_,
        const std::vector<double>& rowDual_,
        const std::vector<HighsBasisStatus>& col_status_,
        const std::vector<HighsBasisStatus>& row_status_)
      : numCol(numCol_),
        numRow(numRow_),
        Astart(Astart_),
        Aend(Aend_),
        Aindex(Aindex_),
        Avalue(Avalue_),
        ARstart(ARstart_),
        ARindex(ARindex_),
        ARvalue(ARvalue_),
        colCost(colCost_),
        colLower(colLower_),
        colUpper(colUpper_),
        rowLower(rowLower_),
        rowUpper(rowUpper_),
        flagCol(flagCol_),
        flagRow(flagRow_),
        colValue(colValue_),
        colDual(colDual_),
        rowValue(rowValue_),
        rowDual(rowDual_),
        col_status(col_status_),
        row_status(row_status_) {}

  const int numCol;
  const int numRow;

  const std::vector<int>& Astart;
  const std::vector<int>& Aend;
  const std::vector<int>& Aindex;
  const std::vector<double>& Avalue;

  const std::vector<int>& ARstart;
  const std::vector<int>& ARindex;
  const std::vector<double>& ARvalue;

  const std::vector<double>& colCost;
  const std::vector<double>& colLower;
  const std::vector<double>& colUpper;
  const std::vector<double>& rowLower;
  const std::vector<double>& rowUpper;

  const std::vector<int>& flagCol;
  const std::vector<int>& flagRow;

  // solution
  const std::vector<double>& colValue;
  const std::vector<double>& colDual;
  const std::vector<double>& rowValue;
  const std::vector<double>& rowDual;

  // basis
  const std::vector<HighsBasisStatus>& col_status;
  const std::vector<HighsBasisStatus>& row_status;
};

enum class KktCondition {
  kColBounds,
  kPrimalFeasibility,
  kDualFeasibility,
  kComplementarySlackness,
  kStationarityOfLagrangian,
  kBasicFeasibleSolution,
  kUnset,
};

struct KktConditionDetails {
  KktConditionDetails() {}
  KktConditionDetails(KktCondition type_) : type(type_) {}

  KktCondition type = KktCondition::kUnset;
  double max_violation = 0.0;
  double sum_violation_2 = 0.0;
  int checked = 0;
  int violated = 0;
};

struct KktInfo {
  std::map<KktCondition, KktConditionDetails> rules;
  bool pass_col_bounds = false;
  bool pass_primal_feas_matrix = false;
  bool pass_dual_feas = false;
  bool pass_st_of_L = false;
  bool pass_comp_slackness = false;
  bool pass_bfs = false;
};

KktInfo initInfo();

bool checkKkt(const State& state, KktInfo info);

void checkPrimalBounds(const State& state, KktConditionDetails& details);
void checkPrimalFeasMatrix(const State& state, KktConditionDetails& details);
void checkDualFeasibility(const State& state, KktConditionDetails& details);
void checkComplementarySlackness(const State& state,
                                 KktConditionDetails& details);
void checkStationarityOfLagrangian(const State& state,
                                   KktConditionDetails& details);
void checkBasicFeasibleSolution(const State& state,
                                KktConditionDetails& details);

}  // namespace dev_kkt_check
}  // namespace presolve

#endif /* TEST_KKTCHSTEP_H_ */
