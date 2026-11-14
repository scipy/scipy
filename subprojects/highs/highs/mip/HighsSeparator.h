/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsSeparator.h
 * @brief Base class for separators
 *
 */

#ifndef MIP_HIGHS_SEPARATOR_H_
#define MIP_HIGHS_SEPARATOR_H_

#include <string>

#include "util/HighsType.h"

const std::string kImplboundSepaString = "Separation: Implied bounds";
const std::string kCliqueSepaString = "Separation: Clique";
const std::string kTableauSepaString = "Separation: Tableau";
const std::string kPathAggrSepaString = "Separation: Path aggregation";
const std::string kModKSepaString = "Separation: Mod-k";
const std::string kMachineSchedSepaString =
    "Separation: Single Machine Scheduling";

class HighsLpRelaxation;
class HighsTransformedLp;
class HighsCutPool;
class HighsLpAggregator;
class HighsMipSolver;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsSeparator {
 private:
  HighsInt numCutsFound;
  HighsInt numCalls;
  int clockIndex;

 public:
  HighsSeparator(const HighsMipSolver& mipsolver, const std::string& name);

  virtual void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                  HighsLpAggregator& lpAggregator,
                                  HighsTransformedLp& transLp,
                                  HighsCutPool& cutpool) = 0;

  void run(HighsLpRelaxation& lpRelaxation, HighsLpAggregator& lpAggregator,
           HighsTransformedLp& transLp, HighsCutPool& cutpool);

  HighsInt getNumCutsFound() const { return numCutsFound; }

  HighsInt getNumCalls() const { return numCalls; }

  HighsInt getClockIndex() const { return clockIndex; }

  virtual ~HighsSeparator() {}
};

#endif
