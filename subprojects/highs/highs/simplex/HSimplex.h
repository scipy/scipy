/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.h
 * @brief
 */
#ifndef SIMPLEX_HSIMPLEX_H_
#define SIMPLEX_HSIMPLEX_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                               HighsInt XnumNewCol);
void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               HighsInt XnumNewCol);

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                            HighsInt XnumNewRow);
void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis,
                            HighsInt XnumNewRow);

void getUnscaledInfeasibilities(const HighsOptions& options,
                                const HighsScale& scale,
                                const SimplexBasis& basis,
                                const HighsSimplexInfo& info,
                                HighsInfo& highs_info);

void setSolutionStatus(HighsInfo& highs_info);
// SCALE:

void scaleSimplexCost(const HighsOptions& options, HighsLp& lp,
                      double& cost_scale);
void unscaleSimplexCost(HighsLp& lp, double cost_scale);

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis);

#endif  // SIMPLEX_HSIMPLEX_H_
