/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HFactorDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HFactorDebug.h"

#include "simplex/HVector.h"
#include "util/HighsRandom.h"

const double solve_large_error = 1e-12;
const double solve_excessive_error = sqrt(solve_large_error);

const double inverse_large_error = 1e-12;
const double inverse_excessive_error = sqrt(inverse_large_error);

HighsDebugStatus debugCheckInvert(const HighsOptions& options,
                                  const HFactor& factor) {
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;

  return_status = HighsDebugStatus::OK;
  const int numRow = factor.numRow;
  const int numCol = factor.numCol;
  const int* Astart = factor.getAstart();
  const int* Aindex = factor.getAindex();
  const double* Avalue = factor.getAvalue();
  const int* baseIndex = factor.getBaseIndex();

  HVector column;
  HVector rhs;
  column.setup(numRow);
  rhs.setup(numRow);
  double rhsDensity = 1;

  // Solve for a random solution
  HighsRandom random;
  column.clear();
  rhs.clear();
  column.count = -1;
  for (int iRow = 0; iRow < numRow; iRow++) {
    rhs.index[rhs.count++] = iRow;
    double value = random.fraction();
    column.array[iRow] = value;
    int iCol = baseIndex[iRow];
    if (iCol < numCol) {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        int index = Aindex[k];
        rhs.array[index] += value * Avalue[k];
      }
    } else {
      int index = iCol - numCol;
      rhs.array[index] += value;
    }
  }
  factor.ftran(rhs, rhsDensity);
  double solve_error_norm = 0;
  for (int iRow = 0; iRow < numRow; iRow++) {
    double solve_error = fabs(rhs.array[iRow] - column.array[iRow]);
    solve_error_norm = std::max(solve_error, solve_error_norm);
  }
  std::string value_adjective;
  int report_level;
  return_status = HighsDebugStatus::OK;

  if (solve_error_norm) {
    if (solve_error_norm > solve_excessive_error) {
      value_adjective = "Excessive";
      report_level = ML_ALWAYS;
      return_status = HighsDebugStatus::ERROR;
    } else if (solve_error_norm > solve_large_error) {
      value_adjective = "Large";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else {
      value_adjective = "Small";
      report_level = ML_VERBOSE;
    }
    HighsPrintMessage(
        options.output, options.message_level, report_level,
        "CheckINVERT:   %-9s (%9.4g) norm for random solution solve error\n",
        value_adjective.c_str(), solve_error_norm);
  }

  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_EXPENSIVE)
    return return_status;

  double columnDensity = 0;
  double inverse_error_norm = 0;
  for (int iRow = 0; iRow < numRow; iRow++) {
    int iCol = baseIndex[iRow];
    column.clear();
    column.packFlag = true;
    if (iCol < numCol) {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        int index = Aindex[k];
        column.array[index] = Avalue[k];
        column.index[column.count++] = index;
      }
    } else {
      int index = iCol - numCol;
      column.array[index] = 1.0;
      column.index[column.count++] = index;
    }
    factor.ftran(column, columnDensity);
    double inverse_column_error_norm = 0;
    for (int lc_iRow = 0; lc_iRow < numRow; lc_iRow++) {
      double value = column.array[lc_iRow];
      double ckValue;
      if (lc_iRow == iRow) {
        ckValue = 1;
      } else {
        ckValue = 0;
      }
      double inverse_error = fabs(value - ckValue);
      inverse_column_error_norm =
          std::max(inverse_error, inverse_column_error_norm);
    }
    inverse_error_norm =
        std::max(inverse_column_error_norm, inverse_error_norm);
  }
  if (inverse_error_norm) {
    if (inverse_error_norm > inverse_excessive_error) {
      value_adjective = "Excessive";
      report_level = ML_ALWAYS;
      return_status = HighsDebugStatus::ERROR;
    } else if (inverse_error_norm > inverse_large_error) {
      value_adjective = "Large";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else {
      value_adjective = "Small";
      report_level = ML_VERBOSE;
    }
    HighsPrintMessage(options.output, options.message_level, report_level,
                      "CheckINVERT:   %-9s (%9.4g) norm for inverse error\n",
                      value_adjective.c_str(), inverse_error_norm);
  }

  return return_status;
}

void debugReportRankDeficiency(const int call_id, const int highs_debug_level,
                               FILE* output, const int message_level,
                               const int numRow, const vector<int>& permute,
                               const vector<int>& iwork, const int* baseIndex,
                               const int rank_deficiency,
                               const vector<int>& noPvR,
                               const vector<int>& noPvC) {
  if (highs_debug_level == HIGHS_DEBUG_LEVEL_NONE) return;
  if (call_id == 0) {
    if (numRow > 123) return;
    printf("buildRankDeficiency0:");
    printf("\nIndex  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", i);
    printf("\nPerm   ");
    for (int i = 0; i < numRow; i++) printf(" %2d", permute[i]);
    printf("\nIwork  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", iwork[i]);
    printf("\nBaseI  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", baseIndex[i]);
    printf("\n");
  } else if (call_id == 1) {
    if (rank_deficiency > 100) return;
    printf("buildRankDeficiency1:");
    printf("\nIndex  ");
    for (int i = 0; i < rank_deficiency; i++) printf(" %2d", i);
    printf("\nnoPvR  ");
    for (int i = 0; i < rank_deficiency; i++) printf(" %2d", noPvR[i]);
    printf("\nnoPvC  ");
    for (int i = 0; i < rank_deficiency; i++) printf(" %2d", noPvC[i]);
    if (numRow > 123) {
      printf("\nIndex  ");
      for (int i = 0; i < numRow; i++) printf(" %2d", i);
      printf("\nIwork  ");
      for (int i = 0; i < numRow; i++) printf(" %2d", iwork[i]);
    }
    printf("\n");
  } else if (call_id == 2) {
    if (numRow > 123) return;
    printf("buildRankDeficiency2:");
    printf("\nIndex  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", i);
    printf("\nPerm   ");
    for (int i = 0; i < numRow; i++) printf(" %2d", permute[i]);
    printf("\n");
  }
}

void debugReportRankDeficientASM(
    const int highs_debug_level, FILE* output, const int message_level,
    const int numRow, const vector<int>& MCstart, const vector<int>& MCcountA,
    const vector<int>& MCindex, const vector<double>& MCvalue,
    const vector<int>& iwork, const int rank_deficiency,
    const vector<int>& noPvC, const vector<int>& noPvR) {
  if (highs_debug_level == HIGHS_DEBUG_LEVEL_NONE) return;
  if (rank_deficiency > 10) return;
  double* ASM;
  ASM = (double*)malloc(sizeof(double) * rank_deficiency * rank_deficiency);
  for (int i = 0; i < rank_deficiency; i++) {
    for (int j = 0; j < rank_deficiency; j++) {
      ASM[i + j * rank_deficiency] = 0;
    }
  }
  for (int j = 0; j < rank_deficiency; j++) {
    int ASMcol = noPvC[j];
    int start = MCstart[ASMcol];
    int end = start + MCcountA[ASMcol];
    for (int en = start; en < end; en++) {
      int ASMrow = MCindex[en];
      int i = -iwork[ASMrow] - 1;
      if (i < 0 || i >= rank_deficiency) {
        printf("STRANGE: 0 > i = %d || %d = i >= rank_deficiency = %d\n", i, i,
               rank_deficiency);
      } else {
        if (noPvR[i] != ASMrow) {
          printf("STRANGE: %d = noPvR[i] != ASMrow = %d\n", noPvR[i], ASMrow);
        }
        printf("Setting ASM(%2d, %2d) = %11.4g\n", i, j, MCvalue[en]);
        ASM[i + j * rank_deficiency] = MCvalue[en];
      }
    }
  }
  printf("ASM:                    ");
  for (int j = 0; j < rank_deficiency; j++) printf(" %11d", j);
  printf("\n                        ");
  for (int j = 0; j < rank_deficiency; j++) printf(" %11d", noPvC[j]);
  printf("\n                        ");
  for (int j = 0; j < rank_deficiency; j++) printf("------------");
  printf("\n");
  for (int i = 0; i < rank_deficiency; i++) {
    printf("%11d %11d|", i, noPvR[i]);
    for (int j = 0; j < rank_deficiency; j++) {
      printf(" %11.4g", ASM[i + j * rank_deficiency]);
    }
    printf("\n");
  }
  free(ASM);
}

void debugReportMarkSingC(const int call_id, const int highs_debug_level,
                          FILE* output, const int message_level,
                          const int numRow, const vector<int>& iwork,
                          const int* baseIndex) {
  if (highs_debug_level == HIGHS_DEBUG_LEVEL_NONE) return;
  if (numRow > 123) return;
  if (call_id == 0) {
    printf("\nMarkSingC1");
    printf("\nIndex  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", i);
    printf("\niwork  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", iwork[i]);
    printf("\nBaseI  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", baseIndex[i]);
  } else if (call_id == 1) {
    printf("\nMarkSingC2");
    printf("\nIndex  ");
    for (int i = 0; i < numRow; i++) printf(" %2d", i);
    printf("\nNwBaseI");
    for (int i = 0; i < numRow; i++) printf(" %2d", baseIndex[i]);
    printf("\n");
  }
}

void debugLogRankDeficiency(const int highs_debug_level, FILE* output,
                            const int message_level, const int rank_deficiency,
                            const int basis_matrix_num_el,
                            const int invert_num_el, const int& kernel_dim,
                            const int kernel_num_el, const int nwork) {
  if (highs_debug_level == HIGHS_DEBUG_LEVEL_NONE) return;
  if (!rank_deficiency) return;
  printf(  // highsPrintMessage(output, message_level, ML_DETAILED,
      "Rank deficiency %1d: basis_matrix (%d el); INVERT (%d el); kernel (%d "
      "dim; %d el): nwork = %d\n",
      rank_deficiency, basis_matrix_num_el, invert_num_el, kernel_dim,
      kernel_num_el, nwork);
}
