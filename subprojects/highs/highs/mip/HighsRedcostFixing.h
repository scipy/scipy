/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsRedcostFixing.h
 * @brief reduced cost fixing using the current cutoff bound
 */

#ifndef HIGHS_REDCOST_FIXING_H_
#define HIGHS_REDCOST_FIXING_H_

#include <map>
#include <vector>

#include "HighsConflictPool.h"
#include "mip/HighsDomainChange.h"

class HighsDomain;
class HighsMipSolver;
class HighsLpRelaxation;

class HighsRedcostFixing {
  std::vector<std::multimap<double, HighsInt>> lurkingColUpper;
  std::vector<std::multimap<double, HighsInt>> lurkingColLower;

 public:
  std::vector<std::pair<double, HighsDomainChange>> getLurkingBounds(
      const HighsMipSolver& mipsolver, const HighsDomain& globaldom) const;

  void propagateRootRedcost(const HighsMipSolver& mipsolver);

  static void propagateRedCost(const HighsMipSolver& mipsolver,
                               HighsDomain& localdomain, HighsDomain& globaldom,
                               const HighsLpRelaxation& lp,
                               HighsConflictPool& conflictpool,
                               HighsPseudocost& pseudocost, double upper_limit);

  void addRootRedcost(const HighsMipSolver& mipsolver,
                      const std::vector<double>& lpredcost, double lpobjective);
};

#endif
