/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HFactorDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HFACTORDEBUG_H_
#define SIMPLEX_HFACTORDEBUG_H_

#include "lp_data/HighsOptions.h"
#include "simplex/HFactor.h"

HighsDebugStatus debugCheckInvert(const HighsOptions& options,
                                  const HFactor& factor);

void debugReportRankDeficiency(const int call_id, const int highs_debug_level,
                               FILE* output, const int message_level,
                               const int numRow, const vector<int>& permute,
                               const vector<int>& iwork, const int* baseIndex,
                               const int rank_deficiency,
                               const vector<int>& noPvR,
                               const vector<int>& noPvC);

void debugReportRankDeficientASM(
    const int highs_debug_level, FILE* output, const int message_level,
    const int numRow, const vector<int>& MCstart, const vector<int>& MCcountA,
    const vector<int>& MCindex, const vector<double>& MCvalue,
    const vector<int>& iwork, const int rank_deficiency,
    const vector<int>& noPvC, const vector<int>& noPvR);

void debugReportMarkSingC(const int call_id, const int highs_debug_level,
                          FILE* output, const int message_level,
                          const int numRow, const vector<int>& iwork,
                          const int* baseIndex);

void debugLogRankDeficiency(const int highs_debug_level, FILE* output,
                            const int message_level, const int rank_deficiency,
                            const int basis_matrix_num_el,
                            const int invert_num_el, const int& kernel_dim,
                            const int kernel_num_el, const int nwork);
#endif  // SIMPLEX_HFACTORDEBUG_H_
