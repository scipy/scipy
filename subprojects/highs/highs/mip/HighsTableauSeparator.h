/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsTableauSeparator.h
 * @brief Class for separating cuts from the LP tableaux rows
 *
 */

#ifndef MIP_HIGHS_TABLEAU_SEPARATOR_H_
#define MIP_HIGHS_TABLEAU_SEPARATOR_H_

#include "mip/HighsSeparator.h"

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsTableauSeparator : public HighsSeparator {
 private:
  int64_t numTries;

 public:
  void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                          HighsLpAggregator& lpAggregator,
                          HighsTransformedLp& transLp,
                          HighsCutPool& cutpool) override;

  HighsTableauSeparator(const HighsMipSolver& mipsolver)
      : HighsSeparator(mipsolver, kTableauSepaString), numTries(0) {}
};

#endif
