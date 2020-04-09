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
#include "simplex/HVector.h"
#include "simplex/SimplexTimer.h"

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

  // delete_Freelist() is being called in Phase 1 and Phase 2 since
  // it's in updatePivots(), but create_Freelist() is only called in
  // Phase 2. Hence freeList and freeListSize are not initialised when
  // freeList.empty() is used to identify that freeListSize should be
  // tested for zero. Suddenly freeListSize is 1212631365 rather than
  // zero when uninitialised, triggering a warning. So, let's set
  // clear freeList and set freeListSize = 0.
  freeList.clear();
  freeListSize = 0;
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

#ifdef HiGHSDEV
  bool rp_Choose_final = false;
  //   rp_Choose_final = true;
#endif
  // 1. Reduce by large step BFRT
  analysis->simplexTimerStart(Chuzc2Clock);
  int fullCount = workCount;
  workCount = 0;
  double totalChange = 0;
  double totalDelta = fabs(workDelta);
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

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 1\n");
#endif
  // 2. Choose by small step BFRT
  analysis->simplexTimerStart(Chuzc3Clock);
  const double Td = workHMO.scaled_solution_params_.dual_feasibility_tolerance;
  fullCount = workCount;
  workCount = 0;
  totalChange = 1e-12;
  selectTheta = workTheta;
  workGroup.clear();
  workGroup.push_back(0);
  const double iz_remainTheta = 1e100;
  int prev_workCount = workCount;
  double prev_remainTheta = iz_remainTheta;
  double prev_selectTheta = selectTheta;
  while (selectTheta < 1e18) {
    double remainTheta = iz_remainTheta;
#ifdef HiGHSDEV
    if (rp_Choose_final)
      printf(
          "Performing choose_final 2; selectTheta = %11.4g; workCount=%d; "
          "fullCount=%d\n",
          selectTheta, workCount, fullCount);
#endif
    for (int i = workCount; i < fullCount; i++) {
      int iCol = workData[i].first;
      double value = workData[i].second;
      double dual = workMove[iCol] * workDual[iCol];
#ifdef HiGHSDEV
      if (rp_Choose_final)
        printf("iCol=%4d; v=%11.4g; d=%11.4g |", iCol, value, dual);
#endif
        // Tight satisfy
#ifdef HiGHSDEV
      if (rp_Choose_final)
        printf(" %11.4g = dual ?<=? sTh * v = %11.4g; workCount=%2d", dual,
               selectTheta * value, workCount);
#endif
      if (dual <= selectTheta * value) {
        swap(workData[workCount++], workData[i]);
        totalChange += value * (workRange[iCol]);
      } else if (dual + Td < remainTheta * value) {
        remainTheta = (dual + Td) / value;
      }
#ifdef HiGHSDEV
      if (rp_Choose_final)
        printf(": totCg=%11.4g; rmTh=%11.4g\n", totalChange, remainTheta);
#endif
    }
    workGroup.push_back(workCount);
    // Update selectTheta with the value of remainTheta;
    selectTheta = remainTheta;
    // Check for no change in this loop - to prevent infinite loop
    if ((workCount == prev_workCount) && (prev_selectTheta == selectTheta) &&
        (prev_remainTheta == remainTheta)) {
#ifdef HiGHSDEV
      printf("In choose_final: No change in loop 2 so return error\n");
      double workDataNorm = 0;
      double dualNorm = 0;
      for (int i = 0; i < workCount; i++) {
        int iCol = workData[i].first;
        double value = workData[i].second;
        workDataNorm += value * value;
        value = workDual[iCol];
        dualNorm += value * value;
      }
      workDataNorm += sqrt(workDataNorm);
      dualNorm += sqrt(dualNorm);
      printf("   workCount = %d; selectTheta=%g; remainTheta=%g\n", workCount,
             selectTheta, remainTheta);
      printf("workDataNorm = %g; dualNorm = %g\n", workDataNorm, dualNorm);
#endif
      analysis->simplexTimerStop(Chuzc3Clock);
      return true;
    }
    // Record the initial values of workCount, remainTheta and selectTheta for
    // the next pass through the loop
    prev_workCount = workCount;
    prev_remainTheta = remainTheta;
    prev_selectTheta = selectTheta;
    if (totalChange >= totalDelta || workCount == fullCount) break;
  }

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 2\n");
#endif
  // 3. Choose large alpha
  double finalCompare = 0;
  for (int i = 0; i < workCount; i++)
    finalCompare = max(finalCompare, workData[i].second);
  finalCompare = min(0.1 * finalCompare, 1.0);
  int countGroup = workGroup.size() - 1;
  int breakGroup = -1;
  int breakIndex = -1;
  for (int iGroup = countGroup - 1; iGroup >= 0; iGroup--) {
    double dMaxFinal = 0;
    int iMaxFinal = -1;
    for (int i = workGroup[iGroup]; i < workGroup[iGroup + 1]; i++) {
      if (dMaxFinal < workData[i].second) {
        dMaxFinal = workData[i].second;
        iMaxFinal = i;
      } else if (dMaxFinal == workData[i].second) {
        int jCol = workData[iMaxFinal].first;
        int iCol = workData[i].first;
        if (workNumTotPermutation[iCol] < workNumTotPermutation[jCol]) {
          iMaxFinal = i;
        }
      }
    }

    if (workData[iMaxFinal].second > finalCompare) {
      breakIndex = iMaxFinal;
      breakGroup = iGroup;
      break;
    }
  }

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 3\n");
#endif
  int sourceOut = workDelta < 0 ? -1 : 1;
  workPivot = workData[breakIndex].first;
  workAlpha = workData[breakIndex].second * sourceOut * workMove[workPivot];
  if (workDual[workPivot] * workMove[workPivot] > 0) {
    workTheta = workDual[workPivot] / workAlpha;
  } else {
    workTheta = 0;
  }

  // 4. Determine BFRT flip index: flip all
  fullCount = breakIndex;
  workCount = 0;
  for (int i = 0; i < workGroup[breakGroup]; i++) {
    const int iCol = workData[i].first;
    const int move = workMove[iCol];
    workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
  }
  if (workTheta == 0) workCount = 0;
  sort(workData.begin(), workData.begin() + workCount);
  analysis->simplexTimerStop(Chuzc3Clock);
#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 4\n");
#endif
  return false;
}

void HDualRow::updateFlip(HVector* bfrtColumn) {
  //  checkDualObjectiveValue("Before update_flip");
  double* workDual = &workHMO.simplex_info_.workDual_[0];  //
  //  double *workLower = &workHMO.simplex_info_.workLower_[0];
  //  double *workUpper = &workHMO.simplex_info_.workUpper_[0];
  //  double *workValue = &workHMO.simplex_info_.workValue_[0];
  double dual_objective_value_change = 0;
  bfrtColumn->clear();
  for (int i = 0; i < workCount; i++) {
    const int iCol = workData[i].first;
    const double change = workData[i].second;

    double lcdual_objective_value_change = change * workDual[iCol];
    //    printf("%6d: [%11.4g, %11.4g, %11.4g], (%11.4g) DlObj = %11.4g
    //    dual_objective_value_change = %11.4g\n",
    //	   iCol, workLower[iCol], workValue[iCol], workUpper[iCol], change,
    // lcdual_objective_value_change, dual_objective_value_change);
    dual_objective_value_change += lcdual_objective_value_change;
    flip_bound(workHMO, iCol);  // workModel->flipBound(iCol);
    workHMO.matrix_.collect_aj(*bfrtColumn, iCol, change);
  }
  workHMO.simplex_info_.updated_dual_objective_value +=
      dual_objective_value_change;
  //  &workHMO.>checkDualObjectiveValue("After  update_flip");
}

void HDualRow::updateDual(double theta) {
  //  &workHMO.>checkDualObjectiveValue("Before update_dual");
  analysis->simplexTimerStart(UpdateDualClock);
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  for (int i = 0; i < packCount; i++) {
    workDual[packIndex[i]] -= theta * packValue[i];
    // Identify the change to the dual objective
    int iCol = packIndex[i];
    double dlDual = theta * packValue[i];
    double iColWorkValue = workHMO.simplex_info_.workValue_[iCol];
    double dlDuObj =
        workHMO.simplex_basis_.nonbasicFlag_[iCol] * (-iColWorkValue * dlDual);
    dlDuObj *= workHMO.scale_.cost_;
    workHMO.simplex_info_.updated_dual_objective_value += dlDuObj;
  }
  analysis->simplexTimerStop(UpdateDualClock);
}

void HDualRow::createFreelist() {
  freeList.clear();
  const int* nonbasicFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  int ckFreeListSize = 0;
  const int numTot = workHMO.simplex_lp_.numCol_ + workHMO.simplex_lp_.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (nonbasicFlag[i] && workRange[i] > 1.5 * HIGHS_CONST_INF) {
      freeList.insert(i);
      ckFreeListSize++;
    }
  }
  if (freeList.size() > 0) {
    //  int freeListSa = *freeList.begin();
    //  int freeListE = *freeList.end();
    freeListSize = *freeList.end();
    if (freeListSize != ckFreeListSize) {
      printf("!! STRANGE: freeListSize != ckFreeListSize\n");
    }
    // const int numTot = workHMO.simplex_lp_.numCol_ +
    // workHMO.simplex_lp_.numRow_;
    //  printf("Create Freelist %d:%d has size %d (%3d%%)\n", freeListSa,
    //  freeListE, freeListSize, 100*freeListSize/numTot);
  }
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
    //  int freeListSa = *freeList.begin();
    //  int freeListE = *freeList.end();
    int ckFreeListSize = 0;
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) ckFreeListSize++;
    freeListSize = *freeList.end();
    if (freeListSize != ckFreeListSize) {
      printf("!! STRANGE: freeListSize != ckFreeListSize\n");
    }
    // const int numTot = workHMO.simplex_lp_.numCol_ +
    // workHMO.simplex_lp_.numRow_;
    //  printf("Update Freelist %d:%d has size %d (%3d%%)\n", freeListSa,
    //  freeListE, freeListSize, 100*freeListSize/numTot); if
    //  (freeList.empty()) {
    //    printf("Empty  Freelist\n");
    //  } else {
    //    printf("\n");
    //  }
  } else {
    if (freeListSize > 0)
      printf("!! STRANGE: Empty Freelist has size %d\n", freeListSize);
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
