/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsModkSeparator.h
 * @brief Class for separating maximally violated mod-k MIR cuts.
 *
 * Contrary to mod-k CG cuts as described in the literature, continuous
 * variables are allowed to appear in the rows used for separation. In case an
 * LP row is already an integral row it is included into the congruence system
 * in the same way as for mod-k CG cuts. Should the LP row contain continuous
 * variables that have a non-zero solution value after bound substitution, then
 * it is discarded, as it can not participate in a maximally violated mod-K MIR
 * cut.
 *
 * If a row contains continuous variables that sit at zero after bound
 * substitution, then those rows are included in the congruence system, as the
 * presence of such variables does not reduce the cuts violation when applying
 * the MIR procedure. In order to handle their presence the row must simply be
 * scaled, such that all integer variables that have a non-zero solution value
 * after bound substitution, as well as the right hand side value, attain an
 * integral value. If we succeed in finding such a scale that is not too large,
 * the resulting row might get a non-zero weight in the solution of the
 * congruence system. The aggregated row therefore can contain continuous
 * variables. These variables, however, all sit at zero in the current LP
 * solution. Using the weights from the solution of the congruence system all
 * integer variables with non-zero solution value will attain a coefficient that
 * is divisible by k, and the integral right hand side value will have a
 * remainder of k - 1 when dividing by k. All other variables do not contribute
 * to the activity of the cut in this LP solution, hence applying the MIR
 * procedure will yield a cut that is violated by (k-1)/k. However, we prefer to
 * generate inequalities with superadditive lifting from the aggregated row
 * whenever all integer variables are bounded.
 *
 */

#ifndef MIP_HIGHS_MODK_SEPARATOR_H_
#define MIP_HIGHS_MODK_SEPARATOR_H_

#include <vector>

#include "mip/HighsSeparator.h"

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsModkSeparator : public HighsSeparator {
 public:
  void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                          HighsLpAggregator& lpAggregator,
                          HighsTransformedLp& transLp,
                          HighsCutPool& cutpool) override;

  HighsModkSeparator(const HighsMipSolver& mipsolver)
      : HighsSeparator(mipsolver, kModKSepaString) {}
};

#endif
