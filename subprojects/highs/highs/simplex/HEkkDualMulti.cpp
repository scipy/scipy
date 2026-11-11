/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkDualMulti.cpp
 * @brief
 */
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <set>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "parallel/HighsParallel.h"
#include "simplex/HEkkDual.h"
#include "simplex/SimplexTimer.h"

using std::cout;
using std::endl;

void HEkkDual::iterateMulti() {
  slice_PRICE = 1;

  // Report candidate
  majorChooseRow();
  minorChooseRow();
  if (row_out == kNoRowChosen) {
    rebuild_reason = kRebuildReasonPossiblyOptimal;
    return;
  }

  // Assign the slice_row_ep, skip if possible
  if (1.0 * multi_finish[multi_nFinish].row_ep->count / solver_num_row < 0.01)
    slice_PRICE = 0;

  if (slice_PRICE) {
    // #pragma omp parallel
    // #pragma omp single
    chooseColumnSlice(multi_finish[multi_nFinish].row_ep);
  } else {
    chooseColumn(multi_finish[multi_nFinish].row_ep);
  }
  // If we failed.
  if (rebuild_reason) {
    if (multi_nFinish) {
      majorUpdate();
    } else {
      highsLogDev(
          ekk_instance_.options_->log_options, HighsLogType::kWarning,
          "PAMI skipping majorUpdate() due to multi_nFinish = %" HIGHSINT_FORMAT
          "; "
          "rebuild_reason = %" HIGHSINT_FORMAT "\n",
          multi_nFinish, rebuild_reason);
    }
    return;
  }

  minorUpdate();
  majorUpdate();
}

void HEkkDual::majorChooseRow() {
  /**
   * 0. Initial check to see if we need to do it again
   */
  if (ekk_instance_.info_.update_count == 0) multi_chooseAgain = 1;
  if (!multi_chooseAgain) return;
  multi_chooseAgain = 0;
  multi_iteration++;
  /**
   * Major loop:
   *     repeat 1-5, until we found a good sets of choices
   */
  std::vector<HighsInt> choiceIndex(multi_num, 0);
  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  for (;;) {
    // 1. Multiple CHUZR
    HighsInt initialCount = 0;

    // Call the hyper-graph method, but partSwitch=0 so just uses
    // choose_multi_global
    dualRHS.chooseMultiHyperGraphAuto(choiceIndex.data(), &initialCount,
                                      multi_num);
    // dualRHS.chooseMultiGlobal(choiceIndex.data(), &initialCount, multi_num);
    if (initialCount == 0 && dualRHS.workCutoff == 0) {
      // OPTIMAL
      return;
    }

    // 2. Shrink the size by cutoff
    HighsInt choiceCount = 0;
    for (HighsInt i = 0; i < initialCount; i++) {
      HighsInt iRow = choiceIndex[i];
      if (dualRHS.work_infeasibility[iRow] / edge_weight[iRow] >=
          dualRHS.workCutoff) {
        choiceIndex[choiceCount++] = iRow;
      }
    }

    if (initialCount == 0 || choiceCount <= initialCount / 3) {
      // Need to do the list again
      dualRHS.createInfeasList(ekk_instance_.info_.col_aq_density);
      continue;
    }

    // 3. Store the choiceIndex to buffer
    for (HighsInt ich = 0; ich < multi_num; ich++)
      multi_choice[ich].row_out = kNoRowChosen;
    for (HighsInt ich = 0; ich < choiceCount; ich++)
      multi_choice[ich].row_out = choiceIndex[ich];

    // 4. Parallel BTRAN and compute weight
    majorChooseRowBtran();

    // 5. Update row densities
    for (HighsInt ich = 0; ich < multi_num; ich++) {
      if (multi_choice[ich].row_out >= 0) {
        const double local_row_ep_density =
            (double)multi_choice[ich].row_ep.count / solver_num_row;
        ekk_instance_.updateOperationResultDensity(
            local_row_ep_density, ekk_instance_.info_.row_ep_density);
      }
    }

    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
      // 6. Check updated and computed weight - just for dual steepest edge
      HighsInt countWrongEdWt = 0;
      for (HighsInt i = 0; i < multi_num; i++) {
        const HighsInt iRow = multi_choice[i].row_out;
        if (iRow < 0) continue;
        double updated_edge_weight = edge_weight[iRow];
        computed_edge_weight = edge_weight[iRow] = multi_choice[i].infeasEdWt;
        //      if (updated_edge_weight < 0.25 * computed_edge_weight) {
        if (!acceptDualSteepestEdgeWeight(updated_edge_weight)) {
          multi_choice[i].row_out = kNoRowChosen;
          countWrongEdWt++;
        }
      }
      if (countWrongEdWt <= choiceCount / 3) break;
    } else {
      // No checking required if not using dual steepest edge so break
      break;
    }
  }

  // 6. Take other info associated with choices
  multi_chosen = 0;
  const double kPamiCutoff = 0.95;
  for (HighsInt i = 0; i < multi_num; i++) {
    const HighsInt iRow = multi_choice[i].row_out;
    if (iRow < 0) continue;
    multi_chosen++;
    // Other info
    multi_choice[i].baseValue = baseValue[iRow];
    multi_choice[i].baseLower = baseLower[iRow];
    multi_choice[i].baseUpper = baseUpper[iRow];
    multi_choice[i].infeasValue = dualRHS.work_infeasibility[iRow];
    multi_choice[i].infeasEdWt = edge_weight[iRow];
    multi_choice[i].infeasLimit =
        dualRHS.work_infeasibility[iRow] / edge_weight[iRow];
    multi_choice[i].infeasLimit *= kPamiCutoff;
  }

  // 6. Finish count
  multi_nFinish = 0;
}

void HEkkDual::majorChooseRowBtran() {
  analysis->simplexTimerStart(BtranClock);

  // 4.1. Prepare BTRAN buffer
  HighsInt multi_ntasks = 0;
  HighsInt multi_iRow[kSimplexConcurrencyLimit];
  HighsInt multi_iwhich[kSimplexConcurrencyLimit];
  double multi_EdWt[kSimplexConcurrencyLimit];
  HVector_ptr multi_vector[kSimplexConcurrencyLimit];
  for (HighsInt ich = 0; ich < multi_num; ich++) {
    if (multi_choice[ich].row_out >= 0) {
      multi_iRow[multi_ntasks] = multi_choice[ich].row_out;
      multi_vector[multi_ntasks] = &multi_choice[ich].row_ep;
      multi_iwhich[multi_ntasks] = ich;
      multi_ntasks++;
    }
  }

  if (analysis->analyse_simplex_summary_data) {
    for (HighsInt i = 0; i < multi_ntasks; i++)
      analysis->operationRecordBefore(kSimplexNlaBtranEp, 1,
                                      ekk_instance_.info_.row_ep_density);
  }
  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  // 4.2 Perform BTRAN
  // #pragma omp parallel for schedule(static, 1)
  // printf("start %d tasks for btran\n", multi_ntasks);
  // std::vector<HighsInt> tmp(multi_ntasks);
  highs::parallel::for_each(0, multi_ntasks, [&](HighsInt start, HighsInt end) {
    for (HighsInt i = start; i < end; i++) {
      // printf("worker %d runs task %i\n", highs::parallel::thread_num(),
      // i);
      // tmp[i] = highs::parallel::thread_num();
      const HighsInt iRow = multi_iRow[i];
      HVector_ptr work_ep = multi_vector[i];
      work_ep->clear();
      work_ep->count = 1;
      work_ep->index[0] = iRow;
      work_ep->array[iRow] = 1;
      work_ep->packFlag = true;
      HighsTimerClock* factor_timer_clock_pointer =
          analysis->getThreadFactorTimerClockPointer();
      ekk_instance_.simplex_nla_.btran(*work_ep,
                                       ekk_instance_.info_.row_ep_density,
                                       factor_timer_clock_pointer);
      if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
        // For Dual steepest edge we know the exact weight as the 2-norm of
        // work_ep
        multi_EdWt[i] = work_ep->norm2();
      } else {
        // For Devex (and Dantzig) we take the updated edge weight
        multi_EdWt[i] = edge_weight[iRow];
      }
    }
  });

  // printf("multi_ntask task schedule:");
  // for (HighsInt i = 0; i < multi_ntasks; ++i) {
  //  printf(" %d->w%d", i, tmp[i]);
  //}
  // printf("\n");

  if (analysis->analyse_simplex_summary_data) {
    for (HighsInt i = 0; i < multi_ntasks; i++)
      analysis->operationRecordAfter(kSimplexNlaBtranEp,
                                     multi_vector[i]->count);
  }
  // 4.3 Put back edge weights: the edge weights for the chosen rows
  // are stored in multi_choice[*].infeasEdWt
  for (HighsInt i = 0; i < multi_ntasks; i++)
    multi_choice[multi_iwhich[i]].infeasEdWt = multi_EdWt[i];

  analysis->simplexTimerStop(BtranClock);
}

void HEkkDual::minorChooseRow() {
  /**
   * 1. Find which to go out
   *        Because we had other checking code
   *        We know current best is OK to be used
   */
  multi_iChoice = -1;
  double bestMerit = 0;
  for (HighsInt ich = 0; ich < multi_num; ich++) {
    const HighsInt iRow = multi_choice[ich].row_out;
    if (iRow < 0) continue;
    double infeasValue = multi_choice[ich].infeasValue;
    double infeasEdWt = multi_choice[ich].infeasEdWt;
    double infeasMerit = infeasValue / infeasEdWt;
    if (bestMerit < infeasMerit) {
      bestMerit = infeasMerit;
      multi_iChoice = ich;
    }
  }

  /**
   * 2. Obtain other info for current sub-optimization choice
   */
  row_out = kNoRowChosen;
  if (multi_iChoice != -1) {
    MChoice* workChoice = &multi_choice[multi_iChoice];

    // Assign useful variables
    row_out = workChoice->row_out;
    variable_out = ekk_instance_.basis_.basicIndex_[row_out];
    double valueOut = workChoice->baseValue;
    double lowerOut = workChoice->baseLower;
    double upperOut = workChoice->baseUpper;
    delta_primal = valueOut - (valueOut < lowerOut ? lowerOut : upperOut);
    move_out = delta_primal < 0 ? -1 : 1;

    // Assign buffers
    MFinish* finish = &multi_finish[multi_nFinish];
    finish->row_out = row_out;
    finish->variable_out = variable_out;
    finish->row_ep = &workChoice->row_ep;
    finish->col_aq = &workChoice->col_aq;
    finish->col_BFRT = &workChoice->col_BFRT;
    // Save the edge weight - over-written later when using Devex
    finish->EdWt = workChoice->infeasEdWt;

    // Disable current row
    workChoice->row_out = kNoRowChosen;
  }
}

void HEkkDual::minorUpdate() {
  // Minor update - store roll back data
  MFinish* finish = &multi_finish[multi_nFinish];
  finish->move_in = ekk_instance_.basis_.nonbasicMove_[variable_in];
  finish->shiftOut = ekk_instance_.info_.workShift_[variable_out];
  finish->flipList.clear();
  for (HighsInt i = 0; i < dualRow.workCount; i++)
    finish->flipList.push_back(dualRow.workData[i].first);

  // Minor update - key parts
  minorUpdateDual();
  minorUpdatePrimal();
  minorUpdatePivots();
  minorUpdateRows();
  if (minor_new_devex_framework) {
    /*
    printf("Iter %7" HIGHSINT_FORMAT " (Major %7" HIGHSINT_FORMAT "): Minor new
    Devex framework\n", ekk_instance_.iteration_count_,
           multi_iteration);
    */
    minorInitialiseDevexFramework();
  }
  multi_nFinish++;
  // Analyse the iteration: possibly report; possibly switch strategy
  iterationAnalysisMinor();

  // Minor update - check for the next iteration
  HighsInt countRemain = 0;
  for (HighsInt i = 0; i < multi_num; i++) {
    HighsInt iRow = multi_choice[i].row_out;
    if (iRow < 0) continue;
    double myInfeas = multi_choice[i].infeasValue;
    double myWeight = multi_choice[i].infeasEdWt;
    countRemain += (myInfeas / myWeight > multi_choice[i].infeasLimit);
  }
  if (countRemain == 0) multi_chooseAgain = 1;
}

void HEkkDual::minorUpdateDual() {
  /**
   * 1. Update the dual solution
   *    XXX Data parallel (depends on the ap partition before)
   */
  if (theta_dual == 0) {
    shiftCost(variable_in, -workDual[variable_in]);
  } else {
    dualRow.updateDual(theta_dual);
    if (slice_PRICE) {
      for (HighsInt i = 0; i < slice_num; i++)
        slice_dualRow[i].updateDual(theta_dual);
    }
  }
  workDual[variable_in] = 0;
  workDual[variable_out] = -theta_dual;
  shiftBack(variable_out);

  /**
   * 2. Apply global bound flip
   */
  dualRow.updateFlip(multi_finish[multi_nFinish].col_BFRT);

  /**
   * 3. Apply local bound flips
   */
  for (HighsInt ich = 0; ich < multi_num; ich++) {
    if (ich == multi_iChoice || multi_choice[ich].row_out >= 0) {
      HVector* this_ep = &multi_choice[ich].row_ep;
      for (HighsInt i = 0; i < dualRow.workCount; i++) {
        double dot = a_matrix->computeDot(*this_ep, dualRow.workData[i].first);
        multi_choice[ich].baseValue -= dualRow.workData[i].second * dot;
      }
    }
  }
}

void HEkkDual::minorUpdatePrimal() {
  MChoice* choice = &multi_choice[multi_iChoice];
  MFinish* finish = &multi_finish[multi_nFinish];
  double valueOut = choice->baseValue;
  double lowerOut = choice->baseLower;
  double upperOut = choice->baseUpper;
  if (delta_primal < 0) {
    theta_primal = (valueOut - lowerOut) / alpha_row;
    finish->basicBound = lowerOut;
  }
  if (delta_primal > 0) {
    theta_primal = (valueOut - upperOut) / alpha_row;
    finish->basicBound = upperOut;
  }
  finish->theta_primal = theta_primal;

  if (edge_weight_mode == EdgeWeightMode::kDevex && !new_devex_framework) {
    std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
    assert(row_out >= 0);
    if (row_out < 0)
      printf("ERROR: row_out = %" HIGHSINT_FORMAT " in minorUpdatePrimal\n",
             row_out);
    const double updated_edge_weight = edge_weight[row_out];
    new_devex_framework = newDevexFramework(updated_edge_weight);
    minor_new_devex_framework = new_devex_framework;
    // Transform the edge weight of the pivotal row according to the
    // simplex update
    double new_pivotal_edge_weight =
        computed_edge_weight / (alpha_row * alpha_row);
    new_pivotal_edge_weight = max(1.0, new_pivotal_edge_weight);
    // Store the Devex weight of the leaving row now - OK since it's
    // stored in finish->EdWt and the updated weights are stored in
    // multi_choice[*].infeasEdWt
    finish->EdWt = new_pivotal_edge_weight;
  }

  /**
   * 5. Update the other primal value
   *    By the pivot (theta_primal)
   */
  for (HighsInt ich = 0; ich < multi_num; ich++) {
    if (multi_choice[ich].row_out >= 0) {
      HVector* this_ep = &multi_choice[ich].row_ep;
      double dot = a_matrix->computeDot(*this_ep, variable_in);
      multi_choice[ich].baseValue -= theta_primal * dot;
      double value = multi_choice[ich].baseValue;
      double lower = multi_choice[ich].baseLower;
      double upper = multi_choice[ich].baseUpper;
      double infeas = 0;
      if (value < lower - Tp) infeas = value - lower;
      if (value > upper + Tp) infeas = value - upper;
      infeas *= infeas;
      multi_choice[ich].infeasValue = infeas;
      if (edge_weight_mode == EdgeWeightMode::kDevex) {
        // Update the other Devex weights
        const double new_pivotal_edge_weight = finish->EdWt;
        double aa_iRow = dot;
        multi_choice[ich].infeasEdWt =
            max(multi_choice[ich].infeasEdWt,
                new_pivotal_edge_weight * aa_iRow * aa_iRow);
      }
    }
  }
}
void HEkkDual::minorUpdatePivots() {
  MFinish* finish = &multi_finish[multi_nFinish];
  ekk_instance_.updatePivots(variable_in, row_out, move_out);
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    // Transform the edge weight of the pivotal row according to the
    // simplex update
    finish->EdWt /= (alpha_row * alpha_row);
  }
  finish->basicValue =
      ekk_instance_.info_.workValue_[variable_in] + theta_primal;
  ekk_instance_.updateMatrix(variable_in, variable_out);
  finish->variable_in = variable_in;
  finish->alpha_row = alpha_row;
  // numericalTrouble is not set in minor iterations, only in
  // majorUpdate, so set it to an illegal value so that its
  // distribution is not updated
  numericalTrouble = -1;
  ekk_instance_.iteration_count_++;
}

void HEkkDual::minorUpdateRows() {
  analysis->simplexTimerStart(UpdateRowClock);
  const HVector* Row = multi_finish[multi_nFinish].row_ep;
  HighsInt updateRows_inDense =
      (Row->count < 0) || (Row->count > 0.1 * solver_num_row);
  if (updateRows_inDense) {
    HighsInt multi_nTasks = 0;
    HighsInt multi_iwhich[kSimplexConcurrencyLimit];
    double multi_xpivot[kSimplexConcurrencyLimit];
    HVector_ptr multi_vector[kSimplexConcurrencyLimit];

    /*
     * Dense mode
     *  1. Find which ones to do and the pivotX
     *  2. Do all of them in task parallel
     */

    // Collect tasks
    for (HighsInt ich = 0; ich < multi_num; ich++) {
      if (multi_choice[ich].row_out >= 0) {
        HVector* next_ep = &multi_choice[ich].row_ep;
        double pivotX = a_matrix->computeDot(*next_ep, variable_in);
        if (fabs(pivotX) < kHighsTiny) continue;
        multi_vector[multi_nTasks] = next_ep;
        multi_xpivot[multi_nTasks] = -pivotX / alpha_row;
        multi_iwhich[multi_nTasks] = ich;
        multi_nTasks++;
      }
    }

    // Perform tasks
    // #pragma omp parallel for schedule(dynamic)
    // printf("minorUpdatesRows: starting %d tasks\n", multi_nTasks);
    highs::parallel::for_each(
        0, multi_nTasks, [&](HighsInt start, HighsInt end) {
          for (HighsInt i = start; i < end; i++) {
            HVector_ptr nextEp = multi_vector[i];
            const double xpivot = multi_xpivot[i];
            nextEp->saxpy(xpivot, Row);
            nextEp->tight();
            if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
              multi_xpivot[i] = nextEp->norm2();
            }
          }
        });

    // Put weight back
    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
      for (HighsInt i = 0; i < multi_nTasks; i++)
        multi_choice[multi_iwhich[i]].infeasEdWt = multi_xpivot[i];
    }
  } else {
    // Sparse mode: just do it sequentially
    for (HighsInt ich = 0; ich < multi_num; ich++) {
      if (multi_choice[ich].row_out >= 0) {
        HVector* next_ep = &multi_choice[ich].row_ep;
        double pivotX = a_matrix->computeDot(*next_ep, variable_in);
        if (fabs(pivotX) < kHighsTiny) continue;
        next_ep->saxpy(-pivotX / alpha_row, Row);
        next_ep->tight();
        if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
          multi_choice[ich].infeasEdWt = next_ep->norm2();
        }
      }
    }
  }
  analysis->simplexTimerStop(UpdateRowClock);
}

void HEkkDual::minorInitialiseDevexFramework() {
  // Set the local Devex weights to 1
  for (HighsInt i = 0; i < multi_num; i++) {
    multi_choice[i].infeasEdWt = 1.0;
  }
  minor_new_devex_framework = false;
}

void HEkkDual::majorUpdate() {
  /**
   * 0. See if it's ready to perform a major update
   */
  if (rebuild_reason) multi_chooseAgain = 1;
  if (!multi_chooseAgain) return;

  // Major update - FTRANs
  majorUpdateFtranPrepare();
  majorUpdateFtranParallel();
  majorUpdateFtranFinal();

  // Major update - check for roll back
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* iFinish = &multi_finish[iFn];
    HVector* iColumn = iFinish->col_aq;
    HighsInt iRow_Out = iFinish->row_out;

    // Use the two pivot values to identify numerical trouble
    if (ekk_instance_.reinvertOnNumericalTrouble(
            "HEkkDual::majorUpdate", numericalTrouble, iColumn->array[iRow_Out],
            iFinish->alpha_row, kMultiNumericalTroubleTolerance)) {
      // HighsInt startUpdate = ekk_instance_.info_.update_count -
      // multi_nFinish;
      rebuild_reason = kRebuildReasonPossiblySingularBasis;
      // if (startUpdate > 0) {
      majorRollback();
      return;
      // }
    }
  }

  // Major update - primal and factor
  majorUpdatePrimal();
  majorUpdateFactor();
  if (new_devex_framework) initialiseDevexFramework();
  iterationAnalysisMajor();
}

void HEkkDual::majorUpdateFtranPrepare() {
  // Prepare FTRAN BFRT buffer
  col_BFRT.clear();
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HVector* Vec = finish->col_BFRT;
    a_matrix->collectAj(*Vec, finish->variable_in, finish->theta_primal);

    // Update this buffer by previous Row_ep
    for (HighsInt jFn = iFn - 1; jFn >= 0; jFn--) {
      MFinish* jFinish = &multi_finish[jFn];
      double* jRow_epArray = jFinish->row_ep->array.data();
      double pivotX = 0;
      for (HighsInt k = 0; k < Vec->count; k++) {
        HighsInt iRow = Vec->index[k];
        pivotX += Vec->array[iRow] * jRow_epArray[iRow];
      }
      if (fabs(pivotX) > kHighsTiny) {
        pivotX /= jFinish->alpha_row;
        a_matrix->collectAj(*Vec, jFinish->variable_in, -pivotX);
        a_matrix->collectAj(*Vec, jFinish->variable_out, pivotX);
      }
    }
    col_BFRT.saxpy(1.0, Vec);
  }

  // Prepare regular FTRAN buffer
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* iFinish = &multi_finish[iFn];
    HVector* iColumn = iFinish->col_aq;
    iColumn->clear();
    iColumn->packFlag = true;
    a_matrix->collectAj(*iColumn, iFinish->variable_in, 1);
  }
}

void HEkkDual::majorUpdateFtranParallel() {
  analysis->simplexTimerStart(FtranMixParClock);

  // Prepare buffers
  HighsInt multi_ntasks = 0;
  double multi_density[kSimplexConcurrencyLimit * 2 + 1];
  HVector_ptr multi_vector[kSimplexConcurrencyLimit * 2 + 1];
  // BFRT first
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordBefore(kSimplexNlaFtranBfrt, col_BFRT.count,
                                    ekk_instance_.info_.col_aq_density);
  multi_density[multi_ntasks] = ekk_instance_.info_.col_aq_density;
  multi_vector[multi_ntasks] = &col_BFRT;
  multi_ntasks++;
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    // Then DSE
    for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
      if (analysis->analyse_simplex_summary_data)
        analysis->operationRecordBefore(kSimplexNlaFtranDse,
                                        multi_finish[iFn].row_ep->count,
                                        ekk_instance_.info_.row_DSE_density);
      multi_density[multi_ntasks] = ekk_instance_.info_.row_DSE_density;
      multi_vector[multi_ntasks] = multi_finish[iFn].row_ep;
      multi_ntasks++;
    }
  }
  // Then Column
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    if (analysis->analyse_simplex_summary_data)
      analysis->operationRecordBefore(kSimplexNlaFtran,
                                      multi_finish[iFn].col_aq->count,
                                      ekk_instance_.info_.col_aq_density);
    multi_density[multi_ntasks] = ekk_instance_.info_.col_aq_density;
    multi_vector[multi_ntasks] = multi_finish[iFn].col_aq;
    multi_ntasks++;
  }

  // Perform FTRAN
  // #pragma omp parallel for schedule(dynamic, 1)
  // printf("majorUpdateFtranParallel: starting %d tasks\n", multi_ntasks);
  highs::parallel::for_each(0, multi_ntasks, [&](HighsInt start, HighsInt end) {
    for (HighsInt i = start; i < end; i++) {
      HVector_ptr rhs = multi_vector[i];
      double density = multi_density[i];
      HighsTimerClock* factor_timer_clock_pointer =
          analysis->getThreadFactorTimerClockPointer();
      ekk_instance_.simplex_nla_.ftran(*rhs, density,
                                       factor_timer_clock_pointer);
    }
  });

  // Update ticks
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HVector* Col = finish->col_aq;
    HVector* Row = finish->row_ep;
    ekk_instance_.total_synthetic_tick_ += Col->synthetic_tick;
    ekk_instance_.total_synthetic_tick_ += Row->synthetic_tick;
  }

  // Update rates
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaFtranBfrt, col_BFRT.count);
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HVector* Col = finish->col_aq;
    HVector* Row = finish->row_ep;
    const double local_col_aq_density = (double)Col->count / solver_num_row;
    ekk_instance_.updateOperationResultDensity(
        local_col_aq_density, ekk_instance_.info_.col_aq_density);
    if (analysis->analyse_simplex_summary_data)
      analysis->operationRecordAfter(kSimplexNlaFtran, Col->count);
    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
      const double local_row_DSE_density = (double)Row->count / solver_num_row;
      ekk_instance_.updateOperationResultDensity(
          local_row_DSE_density, ekk_instance_.info_.row_DSE_density);
      if (analysis->analyse_simplex_summary_data)
        analysis->operationRecordAfter(kSimplexNlaFtranDse, Row->count);
    }
  }
  analysis->simplexTimerStop(FtranMixParClock);
}

void HEkkDual::majorUpdateFtranFinal() {
  analysis->simplexTimerStart(FtranMixFinalClock);
  HighsInt updateFTRAN_inDense = dualRHS.workCount < 0;
  if (updateFTRAN_inDense) {
    for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
      multi_finish[iFn].col_aq->count = -1;
      multi_finish[iFn].row_ep->count = -1;
      double* myCol = multi_finish[iFn].col_aq->array.data();
      double* myRow = multi_finish[iFn].row_ep->array.data();
      for (HighsInt jFn = 0; jFn < iFn; jFn++) {
        HighsInt pivotRow = multi_finish[jFn].row_out;
        const double pivotAlpha = multi_finish[jFn].alpha_row;
        const double* pivotArray = multi_finish[jFn].col_aq->array.data();
        double pivotX1 = myCol[pivotRow];
        double pivotX2 = myRow[pivotRow];

        // The FTRAN regular buffer
        if (fabs(pivotX1) > kHighsTiny) {
          const double pivot = pivotX1 / pivotAlpha;
          // #pragma omp parallel for
          highs::parallel::for_each(
              0, solver_num_row,
              [&](HighsInt start, HighsInt end) {
                for (HighsInt i = start; i < end; i++)
                  myCol[i] -= pivot * pivotArray[i];
              },
              100);
          myCol[pivotRow] = pivot;
        }
        // The FTRAN-DSE buffer
        if (fabs(pivotX2) > kHighsTiny) {
          const double pivot = pivotX2 / pivotAlpha;
          // #pragma omp parallel for
          highs::parallel::for_each(
              0, solver_num_row,
              [&](HighsInt start, HighsInt end) {
                for (HighsInt i = start; i < end; i++)
                  myRow[i] -= pivot * pivotArray[i];
              },
              100);
          myRow[pivotRow] = pivot;
        }
      }
    }
  } else {
    for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
      MFinish* finish = &multi_finish[iFn];
      HVector* Col = finish->col_aq;
      HVector* Row = finish->row_ep;
      for (HighsInt jFn = 0; jFn < iFn; jFn++) {
        MFinish* jFinish = &multi_finish[jFn];
        HighsInt pivotRow = jFinish->row_out;
        double pivotX1 = Col->array[pivotRow];
        // The FTRAN regular buffer
        if (fabs(pivotX1) > kHighsTiny) {
          pivotX1 /= jFinish->alpha_row;
          Col->saxpy(-pivotX1, jFinish->col_aq);
          Col->array[pivotRow] = pivotX1;
        }
        // The FTRAN-DSE buffer
        double pivotX2 = Row->array[pivotRow];
        if (fabs(pivotX2) > kHighsTiny) {
          pivotX2 /= jFinish->alpha_row;
          Row->saxpy(-pivotX2, jFinish->col_aq);
          Row->array[pivotRow] = pivotX2;
        }
      }
    }
  }
  analysis->simplexTimerStop(FtranMixFinalClock);
}

void HEkkDual::majorUpdatePrimal() {
  const bool updatePrimal_inDense = dualRHS.workCount < 0;
  if (updatePrimal_inDense) {
    // Dense update of primal values, infeasibility list and
    // non-pivotal edge weights
    const double* mixArray = col_BFRT.array.data();
    double* local_work_infeasibility = dualRHS.work_infeasibility.data();
    // #pragma omp parallel for schedule(static)
    highs::parallel::for_each(
        0, solver_num_row,
        [&](HighsInt start, HighsInt end) {
          for (HighsInt iRow = start; iRow < end; iRow++) {
            baseValue[iRow] -= mixArray[iRow];
            const double value = baseValue[iRow];
            const double less = baseLower[iRow] - value;
            const double more = value - baseUpper[iRow];
            double infeas = less > Tp ? less : (more > Tp ? more : 0);
            if (ekk_instance_.info_.store_squared_primal_infeasibility)
              local_work_infeasibility[iRow] = infeas * infeas;
            else
              local_work_infeasibility[iRow] = fabs(infeas);
          }
        },
        100);

    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge ||
        (edge_weight_mode == EdgeWeightMode::kDevex && !new_devex_framework)) {
      // Dense update of any edge weights (except weights for pivotal rows)
      std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
      for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
        // multi_finish[iFn].EdWt has already been transformed to correspond to
        // the new basis
        const double new_pivotal_edge_weight = multi_finish[iFn].EdWt;
        const double* colArray = multi_finish[iFn].col_aq->array.data();
        if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
          // Update steepest edge weights
          const double* dseArray = multi_finish[iFn].row_ep->array.data();
          const double Kai = -2 / multi_finish[iFn].alpha_row;
          // #pragma omp parallel for schedule(static)
          highs::parallel::for_each(
              0, solver_num_row,
              [&](HighsInt start, HighsInt end) {
                for (HighsInt iRow = start; iRow < end; iRow++) {
                  const double aa_iRow = colArray[iRow];
                  edge_weight[iRow] +=
                      aa_iRow * (new_pivotal_edge_weight * aa_iRow +
                                 Kai * dseArray[iRow]);
                  edge_weight[iRow] =
                      max(kMinDualSteepestEdgeWeight, edge_weight[iRow]);
                }
              },
              100);
        } else {
          // Update Devex weights
          for (HighsInt iRow = 0; iRow < solver_num_row; iRow++) {
            const double aa_iRow = colArray[iRow];
            edge_weight[iRow] = max(
                edge_weight[iRow], new_pivotal_edge_weight * aa_iRow * aa_iRow);
          }
        }
      }
    }
  } else {
    // Sparse update of primal values, infeasibility list and
    // non-pivotal edge weights
    dualRHS.updatePrimal(&col_BFRT, 1);
    dualRHS.updateInfeasList(&col_BFRT);

    // Sparse update of any edge weights and infeasList. Weights for
    // rows pivotal in this set of MI are updated, but this is based
    // on the previous updated weights not the computed pivotal weight
    // that's known. For the rows pivotal in this set of MI, the
    // weights will be over-written in the next section of code.
    for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
      MFinish* finish = &multi_finish[iFn];
      HVector* Col = finish->col_aq;
      const double new_pivotal_edge_weight = finish->EdWt;
      if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
        // Update steepest edge weights
        HVector* Row = finish->row_ep;
        double Kai = -2 / finish->alpha_row;
        ekk_instance_.updateDualSteepestEdgeWeights(
            finish->row_out, finish->variable_in, Col, new_pivotal_edge_weight,
            Kai, Row->array.data());
      } else if (edge_weight_mode == EdgeWeightMode::kDevex &&
                 !new_devex_framework) {
        // Update Devex weights
        ekk_instance_.updateDualDevexWeights(Col, new_pivotal_edge_weight);
      }
      dualRHS.updateInfeasList(Col);
    }
  }

  // Update primal value for the rows pivotal in this set of MI
  for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HighsInt iRow = finish->row_out;
    double value = baseValue[iRow] - finish->basicBound + finish->basicValue;
    dualRHS.updatePivots(iRow, value);
  }

  // Update any edge weights for the rows pivotal in this set of MI.
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge ||
      (edge_weight_mode == EdgeWeightMode::kDevex && !new_devex_framework)) {
    // Update weights for the pivots using the computed values.
    std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
    for (HighsInt iFn = 0; iFn < multi_nFinish; iFn++) {
      const HighsInt iRow = multi_finish[iFn].row_out;
      const double new_pivotal_edge_weight = multi_finish[iFn].EdWt;
      const double* colArray = multi_finish[iFn].col_aq->array.data();
      // The weight for this pivot is known, but weights for rows
      // pivotal earlier need to be updated
      if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
        // Steepest edge
        const double* dseArray = multi_finish[iFn].row_ep->array.data();
        double Kai = -2 / multi_finish[iFn].alpha_row;
        for (HighsInt jFn = 0; jFn < iFn; jFn++) {
          HighsInt jRow = multi_finish[jFn].row_out;
          double value = colArray[jRow];
          edge_weight[jRow] +=
              value * (new_pivotal_edge_weight * value + Kai * dseArray[jRow]);
          edge_weight[jRow] =
              max(kMinDualSteepestEdgeWeight, edge_weight[jRow]);
        }
        edge_weight[iRow] = new_pivotal_edge_weight;
      } else {
        // Devex
        for (HighsInt jFn = 0; jFn < iFn; jFn++) {
          HighsInt jRow = multi_finish[jFn].row_out;
          const double value = colArray[jRow];
          edge_weight[jRow] =
              max(edge_weight[jRow], new_pivotal_edge_weight * value * value);
        }
        edge_weight[iRow] = new_pivotal_edge_weight;
        num_devex_iterations++;
      }
    }
  }
  checkNonUnitWeightError("999");
}

void HEkkDual::majorUpdateFactor() {
  /**
   * 9. Update the factor by CFT
   */
  HighsInt* iRows = new HighsInt[multi_nFinish];
  for (HighsInt iCh = 0; iCh < multi_nFinish - 1; iCh++) {
    multi_finish[iCh].row_ep->next = multi_finish[iCh + 1].row_ep;
    multi_finish[iCh].col_aq->next = multi_finish[iCh + 1].col_aq;
    iRows[iCh] = multi_finish[iCh].row_out;
  }
  iRows[multi_nFinish - 1] = multi_finish[multi_nFinish - 1].row_out;
  if (multi_nFinish > 0)
    ekk_instance_.updateFactor(multi_finish[0].col_aq, multi_finish[0].row_ep,
                               iRows, &rebuild_reason);

  // Determine whether to reinvert based on the synthetic clock
  const double use_build_synthetic_tick =
      ekk_instance_.build_synthetic_tick_ * kMultiBuildSyntheticTickMu;
  const bool reinvert_syntheticClock =
      ekk_instance_.total_synthetic_tick_ >= use_build_synthetic_tick;
  const bool performed_min_updates =
      ekk_instance_.info_.update_count >=
      kMultiSyntheticTickReinversionMinUpdateCount;
  if (reinvert_syntheticClock && performed_min_updates)
    rebuild_reason = kRebuildReasonSyntheticClockSaysInvert;

  delete[] iRows;
}

void HEkkDual::majorRollback() {
  for (HighsInt iFn = multi_nFinish - 1; iFn >= 0; iFn--) {
    MFinish* finish = &multi_finish[iFn];

    // 1. Roll back pivot
    ekk_instance_.basis_.nonbasicMove_[finish->variable_in] = finish->move_in;
    ekk_instance_.basis_.nonbasicFlag_[finish->variable_in] = 1;
    ekk_instance_.basis_.nonbasicMove_[finish->variable_out] = 0;
    ekk_instance_.basis_.nonbasicFlag_[finish->variable_out] = 0;
    ekk_instance_.basis_.basicIndex_[finish->row_out] = finish->variable_out;

    // 2. Roll back matrix
    ekk_instance_.updateMatrix(finish->variable_out, finish->variable_in);

    // 3. Roll back flips
    for (HighsInt var : finish->flipList) {
      ekk_instance_.flipBound(var);
    }

    // 4. Roll back cost
    ekk_instance_.info_.workShift_[finish->variable_in] = 0;
    ekk_instance_.info_.workShift_[finish->variable_out] = finish->shiftOut;

    // 5. The iteration count
    ekk_instance_.iteration_count_--;
  }
}

bool HEkkDual::checkNonUnitWeightError(const std::string& message) const {
  bool error_found = false;
  if (edge_weight_mode == EdgeWeightMode::kDantzig) {
    std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
    double unit_wt_error = 0;
    for (HighsInt iRow = 0; iRow < solver_num_row; iRow++) {
      unit_wt_error += fabs(edge_weight[iRow] - 1.0);
    }
    error_found = unit_wt_error > 1e-4;
    if (error_found)
      printf("Non-unit Edge weight error of %g: %s\n", unit_wt_error,
             message.c_str());
  }
  return error_found;
}

void HEkkDual::iterationAnalysisMinorData() {
  analysis->multi_iteration_count = multi_iteration;
  analysis->multi_chosen = multi_chosen;
  analysis->multi_finished = multi_nFinish;
}

void HEkkDual::iterationAnalysisMinor() {
  // Possibly report on the iteration
  // PAMI uses alpha_row but serial solver uses alpha
  alpha_col = alpha_row;
  iterationAnalysisData();
  iterationAnalysisMinorData();
  analysis->iterationReport();
  if (analysis->analyse_simplex_summary_data) analysis->iterationRecord();
}

void HEkkDual::iterationAnalysisMajorData() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  analysis->numerical_trouble = numericalTrouble;
  analysis->min_concurrency = info.min_concurrency;
  analysis->num_concurrency = info.num_concurrency;
  analysis->max_concurrency = info.max_concurrency;
}

void HEkkDual::iterationAnalysisMajor() {
  iterationAnalysisMajorData();
  // Possibly switch from DSE to Devex
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    bool switch_to_devex = false;
    //    switch_to_devex = analysis->switchToDevex();
    switch_to_devex = ekk_instance_.switchToDevex();
    if (switch_to_devex) {
      edge_weight_mode = EdgeWeightMode::kDevex;
      // Set up the Devex framework
      initialiseDevexFramework();
    }
  }
  if (analysis->analyse_simplex_summary_data) {
    analysis->iterationRecord();
    analysis->iterationRecordMajor();
  }
}
