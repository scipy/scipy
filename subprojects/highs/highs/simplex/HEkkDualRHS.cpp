/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkDualRHS.cpp
 * @brief
 */
#include "simplex/HEkkDualRHS.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>

#include "../extern/pdqsort/pdqsort.h"
#include "simplex/SimplexTimer.h"

using std::fill_n;
using std::make_pair;
using std::min;
using std::nth_element;
using std::pair;

void HEkkDualRHS::setup() {
  const HighsInt numRow = ekk_instance_.lp_.num_row_;
  workMark.resize(numRow);
  workIndex.resize(numRow);
  work_infeasibility.resize(numRow);
  partNum = 0;
  partSwitch = 0;
  analysis = &ekk_instance_.analysis_;
}

void HEkkDualRHS::chooseNormal(HighsInt* chIndex) {
  // Moved the following to the top to avoid starting the clock for a trivial
  // call. NB Must still call HighsInt to maintain sequence of random numbers
  // for code reproducibility!! Never mind if we're not timing the random number
  // call!!
  // HighsInt random = ekk_instance_.random_.integer();
  if (workCount == 0) {
    *chIndex = -1;
    return;
  }

  // Since chooseNormal calls itself, only start the clock if it's not
  // currently running
  bool keep_timer_running = analysis->simplexTimerRunning(ChuzrDualClock);
  //      timer.clock_start[info.clock_[ChuzrDualClock]] < 0;
  if (!keep_timer_running) {
    analysis->simplexTimerStart(ChuzrDualClock);
  }

  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  if (workCount < 0) {
    // DENSE mode
    const HighsInt numRow = -workCount;
    HighsInt randomStart = ekk_instance_.random_.integer(numRow);
    double bestMerit = 0;
    HighsInt bestIndex = -1;
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? numRow : randomStart;
      for (HighsInt iRow = start; iRow < end; iRow++) {
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = edge_weight[iRow];
          //	  printf("Dense: Row %4" HIGHSINT_FORMAT " weight = %g\n", iRow,
          // myWeight);
          if (bestMerit * myWeight < myInfeas) {
            bestMerit = myInfeas / myWeight;
            bestIndex = iRow;
          }
        }
      }
    }
    *chIndex = bestIndex;
  } else {
    // SPARSE mode
    HighsInt randomStart = ekk_instance_.random_.integer(workCount);
    double bestMerit = 0;
    HighsInt bestIndex = -1;
    std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? workCount : randomStart;
      for (HighsInt i = start; i < end; i++) {
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = edge_weight[iRow];
          /*
          const double myMerit = myInfeas / myWeight;
          printf("CHUZR: iRow = %6" HIGHSINT_FORMAT "; Infeas = %11.4g; Weight =
          %11.4g; Merit = %11.4g\n", iRow, myInfeas, myWeight, myMerit);
          */
          if (bestMerit * myWeight < myInfeas) {
            bestMerit = myInfeas / myWeight;
            bestIndex = iRow;
          }
        }
      }
    }

    HighsInt createListAgain = 0;
    if (bestIndex == -1) {
      createListAgain = workCutoff > 0;
    } else if (bestMerit <= workCutoff * 0.99) {
      createListAgain = 1;
    }
    if (createListAgain) {
      createInfeasList(0);
      chooseNormal(&bestIndex);
    }
    *chIndex = bestIndex;
  }
  // Since chooseNormal calls itself, only stop the clock if it's not currently
  // running
  if (!keep_timer_running) analysis->simplexTimerStop(ChuzrDualClock);
}

void HEkkDualRHS::chooseMultiGlobal(HighsInt* chIndex, HighsInt* chCount,
                                    HighsInt chLimit) {
  analysis->simplexTimerStart(ChuzrDualClock);

  for (HighsInt i = 0; i < chLimit; i++) chIndex[i] = -1;

  const HighsUInt chooseCHECK = chLimit * 2;
  vector<pair<double, HighsInt>> setP;
  setP.reserve(chooseCHECK);

  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  if (workCount < 0) {
    // DENSE mode
    const HighsInt numRow = -workCount;
    HighsInt randomStart = ekk_instance_.random_.integer(numRow);
    double cutoffMerit = 0;
    // Now
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? numRow : randomStart;
      for (HighsInt iRow = start; iRow < end; iRow++) {
        // Was
        //    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
        // Continue
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = edge_weight[iRow];
          if (cutoffMerit * myWeight < myInfeas) {
            // Save
            setP.push_back(make_pair(-myInfeas / myWeight, iRow));
            // Shrink
            if (setP.size() >= chooseCHECK) {
              pdqsort(setP.begin(), setP.end());
              setP.resize(chLimit);
              cutoffMerit = -setP.back().first;
            }
          }
        }
      }
    }
  } else {
    // SPARSE Mode
    HighsInt randomStart;
    if (workCount) {
      randomStart = ekk_instance_.random_.integer(workCount);
    } else {
      // workCount = 0
      randomStart = 0;
    }
    double cutoffMerit = 0;
    // Now
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? workCount : randomStart;
      for (HighsInt i = start; i < end; i++) {
        // Was
        //    for (HighsInt i = 0; i < workCount; i++) {
        // Continue
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = edge_weight[iRow];
          /*
          const double myMerit = myInfeas / myWeight;
          printf("CHUZR: iRow = %6" HIGHSINT_FORMAT "; Infeas = %11.4g; Weight =
          %11.4g; Merit = %11.4g\n", iRow, myInfeas, myWeight, myMerit);
          */
          if (cutoffMerit * myWeight < myInfeas) {
            // Save
            setP.push_back(make_pair(-myInfeas / myWeight, iRow));
            // Shrink
            if (setP.size() >= chooseCHECK) {
              pdqsort(setP.begin(), setP.end());
              setP.resize(chLimit);
              cutoffMerit = -setP.back().first;
            }
          }
        }
      }
    }
  }

  // Store the setP
  pdqsort(setP.begin(), setP.end());
  if ((HighsInt)(setP.size()) > chLimit) setP.resize(chLimit);
  *chCount = static_cast<HighsInt>(setP.size());
  for (size_t i = 0; i < setP.size(); i++) chIndex[i] = setP[i].second;
  analysis->simplexTimerStop(ChuzrDualClock);
}

void HEkkDualRHS::chooseMultiHyperGraphAuto(HighsInt* chIndex,
                                            HighsInt* chCount,
                                            HighsInt chLimit) {
  // Automatically decide to use partition or not
  if (partSwitch)
    chooseMultiHyperGraphPart(chIndex, chCount, chLimit);
  else
    chooseMultiGlobal(chIndex, chCount, chLimit);
}

void HEkkDualRHS::chooseMultiHyperGraphPart(HighsInt* chIndex,
                                            HighsInt* chCount,
                                            HighsInt chLimit) {
  analysis->simplexTimerStart(ChuzrDualClock);

  // Force to use partition method, unless doesn't exist
  if (partNum != chLimit) {
    chooseMultiGlobal(chIndex, chCount, chLimit);
    partSwitch = 0;
    analysis->simplexTimerStop(ChuzrDualClock);
    return;
  }

  // Initialise
  for (HighsInt i = 0; i < chLimit; i++) chIndex[i] = -1;
  *chCount = 0;

  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  if (workCount < 0) {
    // DENSE mode
    const HighsInt numRow = -workCount;
    HighsInt randomStart = ekk_instance_.random_.integer(numRow);
    vector<double> bestMerit(chLimit, 0);
    vector<HighsInt> bestIndex(chLimit, -1);
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? numRow : randomStart;
      for (HighsInt iRow = start; iRow < end; iRow++) {
        if (work_infeasibility[iRow] > kHighsZero) {
          HighsInt iPart = workPartition[iRow];
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = edge_weight[iRow];
          if (bestMerit[iPart] * myWeight < myInfeas) {
            bestMerit[iPart] = myInfeas / myWeight;
            bestIndex[iPart] = iRow;
          }
        }
      }
    }
    HighsInt count = 0;
    for (HighsInt i = 0; i < chLimit; i++) {
      if (bestIndex[i] != -1) {
        chIndex[count++] = bestIndex[i];
      }
    }
    *chCount = count;
  } else {
    // SPARSE mode
    if (workCount == 0) {
      analysis->simplexTimerStop(ChuzrDualClock);
      return;
    }

    HighsInt randomStart = ekk_instance_.random_.integer(workCount);
    vector<double> bestMerit(chLimit, 0);
    vector<HighsInt> bestIndex(chLimit, -1);
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? workCount : randomStart;
      for (HighsInt i = start; i < end; i++) {
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > kHighsZero) {
          HighsInt iPart = workPartition[iRow];
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = edge_weight[iRow];
          if (bestMerit[iPart] * myWeight < myInfeas) {
            bestMerit[iPart] = myInfeas / myWeight;
            bestIndex[iPart] = iRow;
          }
        }
      }
    }
    HighsInt count = 0;
    for (HighsInt i = 0; i < chLimit; i++) {
      if (bestIndex[i] != -1) {
        chIndex[count++] = bestIndex[i];
      }
    }
    *chCount = count;
  }

  analysis->simplexTimerStop(ChuzrDualClock);
}

bool HEkkDualRHS::updatePrimal(HVector* column, double theta) {
  analysis->simplexTimerStart(UpdatePrimalClock);

  const HighsInt numRow = ekk_instance_.lp_.num_row_;
  const HighsInt columnCount = column->count;
  const HighsInt* variable_index = column->index.data();
  const double* columnArray = column->array.data();

  const double* baseLower = ekk_instance_.info_.baseLower_.data();
  const double* baseUpper = ekk_instance_.info_.baseUpper_.data();
  const double Tp = ekk_instance_.options_->primal_feasibility_tolerance;
  double* baseValue = ekk_instance_.info_.baseValue_.data();

  HighsInt num_excessive_primal = 0;
  bool updatePrimal_inDense = columnCount < 0 || columnCount > 0.4 * numRow;

  const HighsInt to_entry = updatePrimal_inDense ? numRow : columnCount;
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow =
        updatePrimal_inDense ? iEntry : variable_index[iEntry];
    baseValue[iRow] -= theta * columnArray[iRow];
    // @primal_infeasibility calculation
    const double value = baseValue[iRow];
    const double lower = baseLower[iRow];
    const double upper = baseUpper[iRow];
    double primal_infeasibility = 0;
    if (value < lower - Tp) {
      primal_infeasibility = lower - value;
    } else if (value > upper + Tp) {
      primal_infeasibility = value - upper;
    }
    if (ekk_instance_.info_.store_squared_primal_infeasibility)
      work_infeasibility[iRow] = primal_infeasibility * primal_infeasibility;
    else
      work_infeasibility[iRow] = fabs(primal_infeasibility);
    if (baseValue[iRow] <= -kExcessivePrimalValue ||
        baseValue[iRow] >= kExcessivePrimalValue)
      num_excessive_primal++;
  }
  analysis->simplexTimerStop(UpdatePrimalClock);
  // Flag detection of excessive values in return
  if (num_excessive_primal) return false;
  return true;
}

void HEkkDualRHS::updatePivots(const HighsInt iRow, const double value) {
  // Update the primal value for the row (iRow) where the basis change
  // has occurred, and set the corresponding primal infeasibility
  // value in work_infeasibility
  //
  const double Tp = ekk_instance_.options_->primal_feasibility_tolerance;
  ekk_instance_.info_.baseValue_[iRow] = value;
  // @primal_infeasibility calculation
  const double lower = ekk_instance_.info_.baseLower_[iRow];
  const double upper = ekk_instance_.info_.baseUpper_[iRow];
  double primal_infeasibility = 0;
  if (value < lower - Tp) {
    primal_infeasibility = lower - value;
  } else if (value > upper + Tp) {
    primal_infeasibility = value - upper;
  }
  if (ekk_instance_.info_.store_squared_primal_infeasibility)
    work_infeasibility[iRow] = primal_infeasibility * primal_infeasibility;
  else
    work_infeasibility[iRow] = fabs(primal_infeasibility);
}

void HEkkDualRHS::updateInfeasList(HVector* column) {
  const HighsInt columnCount = column->count;
  const HighsInt* variable_index = column->index.data();

  // DENSE mode: disabled
  if (workCount < 0) return;

  analysis->simplexTimerStart(UpdatePrimalClock);

  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  if (workCutoff <= 0) {
    // The regular sparse way
    for (HighsInt i = 0; i < columnCount; i++) {
      HighsInt iRow = variable_index[i];
      if (!workMark[iRow]) {
        if (work_infeasibility[iRow]) {
          workIndex[workCount++] = iRow;
          workMark[iRow] = true;
        }
      }
    }
  } else {
    // The hyper sparse way
    for (HighsInt i = 0; i < columnCount; i++) {
      HighsInt iRow = variable_index[i];
      if (!workMark[iRow]) {
        if (work_infeasibility[iRow] > edge_weight[iRow] * workCutoff) {
          workIndex[workCount++] = iRow;
          workMark[iRow] = true;
        }
      }
    }
  }

  analysis->simplexTimerStop(UpdatePrimalClock);
}

void HEkkDualRHS::createArrayOfPrimalInfeasibilities() {
  HighsInt numRow = ekk_instance_.lp_.num_row_;
  const double* baseValue = ekk_instance_.info_.baseValue_.data();
  const double* baseLower = ekk_instance_.info_.baseLower_.data();
  const double* baseUpper = ekk_instance_.info_.baseUpper_.data();
  const double Tp = ekk_instance_.options_->primal_feasibility_tolerance;
  for (HighsInt i = 0; i < numRow; i++) {
    // @primal_infeasibility calculation
    const double value = baseValue[i];
    const double lower = baseLower[i];
    const double upper = baseUpper[i];
    double primal_infeasibility = 0;
    if (value < lower - Tp) {
      primal_infeasibility = lower - value;
    } else if (value > upper + Tp) {
      primal_infeasibility = value - upper;
    }
    if (ekk_instance_.info_.store_squared_primal_infeasibility)
      work_infeasibility[i] = primal_infeasibility * primal_infeasibility;
    else
      work_infeasibility[i] = fabs(primal_infeasibility);
  }
}

void HEkkDualRHS::createInfeasList(double columnDensity) {
  HighsInt numRow = ekk_instance_.lp_.num_row_;
  double* dwork = ekk_instance_.scattered_dual_edge_weight_.data();

  // 1. Build the full list
  fill_n(workMark.data(), numRow, false);
  workCount = 0;
  workCutoff = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    if (work_infeasibility[iRow]) {
      workMark[iRow] = true;
      workIndex[workCount++] = iRow;
    }
  }

  // 2. See if it worth to try to go sparse
  //    (Many candidates, really sparse RHS)
  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  if (workCount > max(numRow * 0.01, 500.0) && columnDensity < 0.05) {
    HighsInt icutoff = max(workCount * 0.001, 500.0);
    double maxMerit = 0;
    for (HighsInt iRow = 0, iPut = 0; iRow < numRow; iRow++)
      if (workMark[iRow]) {
        double myMerit = work_infeasibility[iRow] / edge_weight[iRow];
        if (maxMerit < myMerit) maxMerit = myMerit;
        dwork[iPut++] = -myMerit;
      }
    nth_element(dwork, dwork + icutoff, dwork + workCount);
    double cutMerit = -dwork[icutoff];
    workCutoff = min(maxMerit * 0.99999, cutMerit * 1.00001);

    // Create again
    fill_n(workMark.data(), numRow, false);
    workCount = 0;
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      if (work_infeasibility[iRow] >= edge_weight[iRow] * workCutoff) {
        workIndex[workCount++] = iRow;
        workMark[iRow] = true;
      }
    }

    // Reduce by drop smaller
    if (workCount > icutoff * 1.5) {
      // Firstly take up "icutoff" number of elements
      HighsInt fullCount = workCount;
      workCount = icutoff;
      for (HighsInt i = icutoff; i < fullCount; i++) {
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > edge_weight[iRow] * cutMerit) {
          workIndex[workCount++] = iRow;
        } else {
          workMark[iRow] = false;
        }
      }
    }
  }

  // 3. If there are still too many candidates: disable them
  if (workCount > 0.2 * numRow) {
    workCount = -numRow;
    workCutoff = 0;
  }
}

void HEkkDualRHS::assessOptimality() {
  HighsInt num_work_infeasibilities = 0;
  double max_work_infeasibility = 0;
  const HighsInt num_row = ekk_instance_.lp_.num_row_;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (work_infeasibility[iRow] > kHighsZero) {
      num_work_infeasibilities++;
      max_work_infeasibility =
          std::max(work_infeasibility[iRow], max_work_infeasibility);
    }
  }
  ekk_instance_.computeSimplexPrimalInfeasible();
  const HighsInt num_primal_infeasibilities =
      ekk_instance_.info_.num_primal_infeasibilities;
  const double max_primal_infeasibility =
      ekk_instance_.info_.max_primal_infeasibility;
  const bool regular_report = false;
  if (regular_report ||
      (num_work_infeasibilities && !num_primal_infeasibilities)) {
    printf(
        "assessOptimality: %6d rows; workCount = %4d (%6.4f) "
        "num / max infeasibilities: work = %4d / %11.4g; simplex = %4d / "
        "%11.4g: %s\n",
        (int)num_row, (int)workCount,
        workCount > 0 ? (1.0 * workCount) / num_row : 0,
        (int)num_work_infeasibilities, max_work_infeasibility,
        (int)num_primal_infeasibilities, max_primal_infeasibility,
        num_primal_infeasibilities == 0 ? "Optimal" : "");
    if (num_work_infeasibilities && !num_primal_infeasibilities) {
      printf("assessOptimality: call %d; tick %d; iter %d\n",
             (int)ekk_instance_.debug_solve_call_num_,
             (int)ekk_instance_.debug_initial_build_synthetic_tick_,
             (int)ekk_instance_.iteration_count_);
    }
  }
}
