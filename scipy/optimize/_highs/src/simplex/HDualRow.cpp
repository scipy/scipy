/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDualRow.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HDualRow.h"

#include <cassert>
#include <iostream>

#include "lp_data/HConst.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HSimplex.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HVector.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsSort.h"

using std::make_pair;
using std::pair;
using std::set;

void HDualRow::setupSlice(int size) {
  workSize = size;
  workMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  workDual = &workHMO.simplex_info_.workDual_[0];
  workRange = &workHMO.simplex_info_.workRange_[0];
  work_devex_index = &workHMO.simplex_info_.devex_index_[0];

  // Allocate spaces
  packCount = 0;
  packIndex.resize(workSize);
  packValue.resize(workSize);

  workCount = 0;
  workData.resize(workSize);
  analysis = &workHMO.simplex_analysis_;
}

void HDualRow::setup() {
  // Setup common vectors
  const int numTot = workHMO.simplex_lp_.numCol_ + workHMO.simplex_lp_.numRow_;
  setupSlice(numTot);
  workNumTotPermutation = &workHMO.simplex_info_.numTotPermutation_[0];

  // deleteFreelist() is being called in Phase 1 and Phase 2 since
  // it's in updatePivots(), but create_Freelist() is only called in
  // Phase 2. Hence freeList is not initialised when freeList.empty()
  // is used in deleteFreelist(), clear freeList now.
  freeList.clear();
}

void HDualRow::clear() {
  packCount = 0;
  workCount = 0;
}

void HDualRow::chooseMakepack(const HVector* row, const int offset) {
  /**
   * Pack the indices and values for the row
   *
   * Offset of numCol is used when packing row_ep
   */
  const int rowCount = row->count;
  const int* rowIndex = &row->index[0];
  const double* rowArray = &row->array[0];

  for (int i = 0; i < rowCount; i++) {
    const int index = rowIndex[i];
    const double value = rowArray[index];
    packIndex[packCount] = index + offset;
    packValue[packCount++] = value;
  }
}

void HDualRow::choosePossible() {
  /**
   * Determine the possible variables - candidates for CHUZC
   * TODO: Check with Qi what this is doing
   */
  const double Ta = workHMO.simplex_info_.update_count < 10
                        ? 1e-9
                        : workHMO.simplex_info_.update_count < 20 ? 3e-8 : 1e-6;
  const double Td = workHMO.scaled_solution_params_.dual_feasibility_tolerance;
  const int sourceOut = workDelta < 0 ? -1 : 1;
  workTheta = HIGHS_CONST_INF;
  workCount = 0;
  for (int i = 0; i < packCount; i++) {
    const int iCol = packIndex[i];
    const int move = workMove[iCol];
    const double alpha = packValue[i] * sourceOut * move;
    if (alpha > Ta) {
      workData[workCount++] = make_pair(iCol, alpha);
      const double relax = workDual[iCol] * move + Td;
      if (workTheta * alpha > relax) workTheta = relax / alpha;
    }
  }
}

void HDualRow::chooseJoinpack(const HDualRow* otherRow) {
  /**
   * Join pack of possible candidates in this row with possible
   * candidates in otherRow
   */
  const int otherCount = otherRow->workCount;
  const pair<int, double>* otherData = &otherRow->workData[0];
  copy(otherData, otherData + otherCount, &workData[workCount]);
  workCount = workCount + otherCount;
  workTheta = min(workTheta, otherRow->workTheta);
}

bool HDualRow::chooseFinal() {
  /**
   * Chooses the entering variable via BFRT and EXPAND
   *
   * It will
   * (1) reduce the candidates as a small collection
   * (2) choose by BFRT by going over break points
   * (3) choose final by alpha
   * (4) determine final flip variables
   */

  // 1. Reduce by large step BFRT
  analysis->simplexTimerStart(Chuzc2Clock);
  int fullCount = workCount;
  workCount = 0;
  double totalChange = 0;
  const double totalDelta = fabs(workDelta);
  double selectTheta = 10 * workTheta + 1e-7;
  for (;;) {
    for (int i = workCount; i < fullCount; i++) {
      int iCol = workData[i].first;
      double alpha = workData[i].second;
      double tight = workMove[iCol] * workDual[iCol];
      if (alpha * selectTheta >= tight) {
        swap(workData[workCount++], workData[i]);
        totalChange += workRange[iCol] * alpha;
      }
    }
    selectTheta *= 10;
    if (totalChange >= totalDelta || workCount == fullCount) break;
  }
  analysis->simplexTimerStop(Chuzc2Clock);

  // 2. Choose by small step BFRT

  bool use_quad_sort = false;
  bool use_heap_sort = false;
  const int dual_chuzc_sort_strategy =
      workHMO.options_.dual_chuzc_sort_strategy;
  if (dual_chuzc_sort_strategy == SIMPLEX_DUAL_CHUZC_STRATEGY_CHOOSE) {  // 0
    // Use the quadratic cost sort for smaller values of workCount,
    // otherwise use the heap-based sort
    use_quad_sort = workCount < 100;
    use_heap_sort = !use_quad_sort;
  } else if (dual_chuzc_sort_strategy ==
             SIMPLEX_DUAL_CHUZC_STRATEGY_QUAD) {  // 1
    // Use the quadratic cost sort
    use_quad_sort = true;
  } else if (dual_chuzc_sort_strategy ==
             SIMPLEX_DUAL_CHUZC_STRATEGY_HEAP) {  // 2
    // Use the heap-based sort
    use_heap_sort = true;
  } else if (dual_chuzc_sort_strategy ==
             SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH) {  // 3
    // Use the both sorts - for debugging
    use_quad_sort = true;
    use_heap_sort = true;
  }
  // Ensure that at least one sort is used!
  assert(use_heap_sort || use_quad_sort);

  if (use_heap_sort) {
    // Take a copy of workData and workCount for the independent
    // heap-based code
    original_workData = workData;
    alt_workCount = workCount;
  }
  analysis->simplexTimerStart(Chuzc3Clock);
  if (use_quad_sort) {
    // Use the O(n^2) quadratic sort for the candidates
    analysis->simplexTimerStart(Chuzc3a0Clock);
    bool choose_ok = chooseFinalWorkGroupQuad();
    analysis->simplexTimerStop(Chuzc3a0Clock);
    if (!choose_ok) {
      analysis->simplexTimerStop(Chuzc3Clock);
      return true;
    }
  }
  if (use_heap_sort) {
    // Use the O(n log n) heap sort for the candidates
    analysis->simplexTimerStart(Chuzc3a1Clock);
    chooseFinalWorkGroupHeap();
    analysis->simplexTimerStop(Chuzc3a1Clock);
  }
  // 3. Choose large alpha
  analysis->simplexTimerStart(Chuzc3bClock);
  int breakIndex;
  int breakGroup;
  int alt_breakIndex;
  int alt_breakGroup;
  if (use_quad_sort)
    chooseFinalLargeAlpha(breakIndex, breakGroup, workCount, workData,
                          workGroup);
  if (use_heap_sort)
    chooseFinalLargeAlpha(alt_breakIndex, alt_breakGroup, alt_workCount,
                          sorted_workData, alt_workGroup);
  analysis->simplexTimerStop(Chuzc3bClock);

  if (!use_quad_sort) {
    // If the quadratic sort is not being used, revert to the heap
    // sort results
    breakIndex = alt_breakIndex;
    breakGroup = alt_breakGroup;
  }
  analysis->simplexTimerStart(Chuzc3cClock);

  int sourceOut = workDelta < 0 ? -1 : 1;
  if (use_quad_sort) {
    workPivot = workData[breakIndex].first;
    workAlpha = workData[breakIndex].second * sourceOut * workMove[workPivot];
  } else {
    workPivot = sorted_workData[breakIndex].first;
    workAlpha =
        sorted_workData[breakIndex].second * sourceOut * workMove[workPivot];
  }
  if (workDual[workPivot] * workMove[workPivot] > 0) {
    workTheta = workDual[workPivot] / workAlpha;
  } else {
    workTheta = 0;
  }

  analysis->simplexTimerStop(Chuzc3cClock);

  if (use_quad_sort && use_heap_sort)
    debugDualChuzcWorkDataAndGroup(
        workHMO, workDelta, workTheta, workCount, alt_workCount, breakIndex,
        alt_breakIndex, workData, sorted_workData, workGroup, alt_workGroup);

  analysis->simplexTimerStart(Chuzc3dClock);

  // 4. Determine BFRT flip index: flip all
  fullCount = breakIndex;
  workCount = 0;
  if (use_quad_sort) {
    for (int i = 0; i < workGroup[breakGroup]; i++) {
      const int iCol = workData[i].first;
      const int move = workMove[iCol];
      workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
    }
  } else {
    for (int i = 0; i < alt_workGroup[breakGroup]; i++) {
      const int iCol = sorted_workData[i].first;
      const int move = workMove[iCol];
      workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
    }
  }
  if (workTheta == 0) workCount = 0;
  analysis->simplexTimerStop(Chuzc3dClock);

  analysis->simplexTimerStart(Chuzc3eClock);
  /*
  if (!use_quad_sort) {
    for (int i = 0; i < workCount; i++) workData[i] = sorted_workData[i];
  }
  */
  sort(workData.begin(), workData.begin() + workCount);
  analysis->simplexTimerStop(Chuzc3eClock);
  analysis->simplexTimerStop(Chuzc3Clock);
  return false;
}

bool HDualRow::chooseFinalWorkGroupQuad() {
  const double Td = workHMO.scaled_solution_params_.dual_feasibility_tolerance;
  int fullCount = workCount;
  workCount = 0;
  double totalChange = initial_total_change;
  double selectTheta = workTheta;
  const double totalDelta = fabs(workDelta);
  workGroup.clear();
  workGroup.push_back(0);
  int prev_workCount = workCount;
  double prev_remainTheta = initial_remain_theta;
  double prev_selectTheta = selectTheta;
  int debug_num_loop = 0;

  while (selectTheta < max_select_theta) {
    double remainTheta = initial_remain_theta;
    debug_num_loop++;
    int debug_loop_ln = 0;
    for (int i = workCount; i < fullCount; i++) {
      int iCol = workData[i].first;
      double value = workData[i].second;
      double dual = workMove[iCol] * workDual[iCol];
      // Tight satisfy
      if (dual <= selectTheta * value) {
        swap(workData[workCount++], workData[i]);
        totalChange += value * (workRange[iCol]);
      } else if (dual + Td < remainTheta * value) {
        remainTheta = (dual + Td) / value;
      }
      debug_loop_ln++;
    }
    workGroup.push_back(workCount);

    // Update selectTheta with the value of remainTheta;
    selectTheta = remainTheta;
    // Check for no change in this loop - to prevent infinite loop
    if ((workCount == prev_workCount) && (prev_selectTheta == selectTheta) &&
        (prev_remainTheta == remainTheta)) {
      debugDualChuzcFail(workHMO.options_, workCount, workData, workDual,
                         selectTheta, remainTheta);
      return false;
    }
    // Record the initial values of workCount, remainTheta and selectTheta for
    // the next pass through the loop - to check for infinite loop condition
    prev_workCount = workCount;
    prev_remainTheta = remainTheta;
    prev_selectTheta = selectTheta;
    if (totalChange >= totalDelta || workCount == fullCount) break;
  }
  return true;
}

bool HDualRow::chooseFinalWorkGroupHeap() {
  const double Td = workHMO.scaled_solution_params_.dual_feasibility_tolerance;
  int fullCount = alt_workCount;
  double totalChange = initial_total_change;
  double selectTheta = workTheta;
  const double totalDelta = fabs(workDelta);
  int heap_num_en = 0;
  std::vector<int> heap_i;
  std::vector<double> heap_v;
  heap_i.resize(fullCount + 1);
  heap_v.resize(fullCount + 1);
  for (int i = 0; i < fullCount; i++) {
    int iCol = original_workData[i].first;
    double value = original_workData[i].second;
    double dual = workMove[iCol] * workDual[iCol];
    double ratio = dual / value;
    if (ratio < max_select_theta) {
      heap_num_en++;
      heap_i[heap_num_en] = i;
      heap_v[heap_num_en] = ratio;
    }
  }
  maxheapsort(&heap_v[0], &heap_i[0], heap_num_en);

  alt_workCount = 0;
  alt_workGroup.clear();
  alt_workGroup.push_back(alt_workCount);
  int this_group_first_entry = alt_workCount;
  sorted_workData.resize(heap_num_en);
  for (int en = 1; en <= heap_num_en; en++) {
    int i = heap_i[en];
    int iCol = original_workData[i].first;
    double value = original_workData[i].second;
    double dual = workMove[iCol] * workDual[iCol];
    if (dual > selectTheta * value) {
      // Breakpoint is in the next group, so record the pointer to its
      // first entry
      alt_workGroup.push_back(alt_workCount);
      this_group_first_entry = alt_workCount;
      selectTheta = (dual + Td) / value;
      // End loop if all permitted groups have been identified
      if (totalChange >= totalDelta) break;
    }
    // Store the breakpoint
    sorted_workData[alt_workCount].first = iCol;
    sorted_workData[alt_workCount].second = value;
    totalChange += value * (workRange[iCol]);
    alt_workCount++;
  }
  if (alt_workCount > this_group_first_entry)
    alt_workGroup.push_back(alt_workCount);
  return true;
}

void HDualRow::chooseFinalLargeAlpha(
    int& breakIndex, int& breakGroup, int pass_workCount,
    const std::vector<std::pair<int, double>>& pass_workData,
    const std::vector<int>& pass_workGroup) {
  double finalCompare = 0;
  for (int i = 0; i < pass_workCount; i++)
    finalCompare = max(finalCompare, pass_workData[i].second);
  finalCompare = min(0.1 * finalCompare, 1.0);
  int countGroup = pass_workGroup.size() - 1;
  breakGroup = -1;
  breakIndex = -1;
  for (int iGroup = countGroup - 1; iGroup >= 0; iGroup--) {
    double dMaxFinal = 0;
    int iMaxFinal = -1;
    for (int i = pass_workGroup[iGroup]; i < pass_workGroup[iGroup + 1]; i++) {
      if (dMaxFinal < pass_workData[i].second) {
        dMaxFinal = pass_workData[i].second;
        iMaxFinal = i;
      } else if (dMaxFinal == pass_workData[i].second) {
        int jCol = pass_workData[iMaxFinal].first;
        int iCol = pass_workData[i].first;
        if (workNumTotPermutation[iCol] < workNumTotPermutation[jCol]) {
          iMaxFinal = i;
        }
      }
    }

    if (pass_workData[iMaxFinal].second > finalCompare) {
      breakIndex = iMaxFinal;
      breakGroup = iGroup;
      break;
    }
  }
}

void HDualRow::updateFlip(HVector* bfrtColumn) {
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  double dual_objective_value_change = 0;
  bfrtColumn->clear();
  for (int i = 0; i < workCount; i++) {
    const int iCol = workData[i].first;
    const double change = workData[i].second;
    double local_dual_objective_change = change * workDual[iCol];
    local_dual_objective_change *= workHMO.scale_.cost_;
    dual_objective_value_change += local_dual_objective_change;
    flip_bound(workHMO, iCol);
    workHMO.matrix_.collect_aj(*bfrtColumn, iCol, change);
  }
  workHMO.simplex_info_.updated_dual_objective_value +=
      dual_objective_value_change;
}

void HDualRow::updateDual(double theta) {
  analysis->simplexTimerStart(UpdateDualClock);
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  double dual_objective_value_change = 0;
  for (int i = 0; i < packCount; i++) {
    workDual[packIndex[i]] -= theta * packValue[i];
    // Identify the change to the dual objective
    int iCol = packIndex[i];
    const double delta_dual = theta * packValue[i];
    const double local_value = workHMO.simplex_info_.workValue_[iCol];
    double local_dual_objective_change =
        workHMO.simplex_basis_.nonbasicFlag_[iCol] *
        (-local_value * delta_dual);
    local_dual_objective_change *= workHMO.scale_.cost_;
    dual_objective_value_change += local_dual_objective_change;
  }
  workHMO.simplex_info_.updated_dual_objective_value +=
      dual_objective_value_change;
  analysis->simplexTimerStop(UpdateDualClock);
}

void HDualRow::createFreelist() {
  freeList.clear();
  for (int i = 0; i < workHMO.simplex_lp_.numCol_ + workHMO.simplex_lp_.numRow_;
       i++) {
    if (workHMO.simplex_basis_.nonbasicFlag_[i] &&
        highs_isInfinity(-workHMO.simplex_info_.workLower_[i]) &&
        highs_isInfinity(workHMO.simplex_info_.workUpper_[i]))
      freeList.insert(i);
  }
  debugFreeListNumEntries(workHMO, freeList);
}

void HDualRow::createFreemove(HVector* row_ep) {
  // TODO: Check with Qi what this is doing and why it's expensive
  if (!freeList.empty()) {
    double Ta = workHMO.simplex_info_.update_count < 10
                    ? 1e-9
                    : workHMO.simplex_info_.update_count < 20 ? 3e-8 : 1e-6;
    int sourceOut = workDelta < 0 ? -1 : 1;
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) {
      int iCol = *sit;
      assert(iCol < workHMO.simplex_lp_.numCol_);
      double alpha = workHMO.matrix_.compute_dot(*row_ep, iCol);
      if (fabs(alpha) > Ta) {
        if (alpha * sourceOut > 0)
          workHMO.simplex_basis_.nonbasicMove_[iCol] = 1;
        else
          workHMO.simplex_basis_.nonbasicMove_[iCol] = -1;
      }
    }
  }
}
void HDualRow::deleteFreemove() {
  if (!freeList.empty()) {
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) {
      int iCol = *sit;
      assert(iCol < workHMO.simplex_lp_.numCol_);
      workHMO.simplex_basis_.nonbasicMove_[iCol] = 0;
    }
  }
}

void HDualRow::deleteFreelist(int iColumn) {
  if (!freeList.empty()) {
    if (freeList.count(iColumn)) freeList.erase(iColumn);
  }
}

void HDualRow::computeDevexWeight(const int slice) {
  const bool rp_computed_edge_weight = false;
  computed_edge_weight = 0;
  for (int el_n = 0; el_n < packCount; el_n++) {
    int vr_n = packIndex[el_n];
    if (!workHMO.simplex_basis_.nonbasicFlag_[vr_n]) {
      //      printf("Basic variable %d in packIndex is skipped\n", vr_n);
      continue;
    }
    double pv = work_devex_index[vr_n] * packValue[el_n];
    if (pv) {
      computed_edge_weight += pv * pv;
    }
  }
  if (rp_computed_edge_weight) {
    if (slice >= 0)
      printf(
          "HDualRow::computeDevexWeight: Slice %1d; computed_edge_weight = "
          "%11.4g\n",
          slice, computed_edge_weight);
  }
}
