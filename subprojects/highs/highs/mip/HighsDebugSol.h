/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsDebugSol.h
 * @brief Debug solution for MIP solver
 */

#ifndef HIGHS_DEBUG_SOL_H_
#define HIGHS_DEBUG_SOL_H_

class HighsDomain;
class HighsMipSolver;
class HighsLp;

#include <set>
#include <vector>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomain.h"

#ifdef HIGHS_DEBUGSOL

#include <unordered_map>

struct HighsDebugSol {
  const HighsMipSolver* mipsolver;
  double debugSolObjective;
  std::vector<double> debugOrigSolution;
  std::vector<double> debugSolution;
  bool debugSolActive;
  std::unordered_map<const HighsDomain*, std::multiset<HighsDomainChange>>
      conflictingBounds;

  HighsDebugSol(HighsMipSolver& mipsolver);

  void newIncumbentFound();

  void activate();

  void shrink(const std::vector<HighsInt>& newColIndex);

  void registerDomain(const HighsDomain& domain);

  void boundChangeAdded(const HighsDomain& domain,
                        const HighsDomainChange& domchg,
                        bool branching = false);

  void boundChangeRemoved(const HighsDomain& domain,
                          const HighsDomainChange& domchg);

  void resetDomain(const HighsDomain& domain);

  void nodePruned(const HighsDomain& localdomain, double cutoff = kHighsInf);

  void checkCut(const HighsInt* Rindex, const double* Rvalue, HighsInt Rlen,
                double rhs);

  void checkRow(const HighsInt* Rindex, const double* Rvalue, HighsInt Rlen,
                double Rlower, double Rupper);

  void checkRowAggregation(const HighsLp& lp, const HighsInt* Rindex,
                           const double* Rvalue, HighsInt Rlen);

  void checkClique(const HighsCliqueTable::CliqueVar* clq, HighsInt clqlen);

  void checkVub(HighsInt col, HighsInt vubcol, double vubcoef,
                double vubconstant) const;

  void checkVlb(HighsInt col, HighsInt vlbcol, double vlbcoef,
                double vlbconstant) const;

  void checkConflictReasonFrontier(
      const std::set<HighsDomain::ConflictSet::LocalDomChg>& reasonSideFrontier,
      const std::vector<HighsDomainChange>& domchgstack) const;

  void checkConflictReconvergenceFrontier(
      const std::set<HighsDomain::ConflictSet::LocalDomChg>&
          reconvergenceFrontier,
      const HighsDomain::ConflictSet::LocalDomChg& reconvDomchgPos,
      const std::vector<HighsDomainChange>& domchgstack) const;
};

#else
struct HighsDebugSol {
  HighsDebugSol(HighsMipSolver&) {}

  void newIncumbentFound() const {}

  void activate() const {}

  void shrink(const std::vector<HighsInt>&) const {}

  void registerDomain(const HighsDomain&) const {}

  void boundChangeAdded(const HighsDomain&, const HighsDomainChange&,
                        bool = false) const {}

  void boundChangeRemoved(const HighsDomain&, const HighsDomainChange&) const {}

  void resetDomain(const HighsDomain&) const {}

  void nodePruned(const HighsDomain&, double = kHighsInf) const {}

  void checkCut(const HighsInt*, const double*, HighsInt, double) const {}

  void checkRow(const HighsInt*, const double*, HighsInt, double,
                double) const {}

  void checkRowAggregation(const HighsLp&, const HighsInt*, const double*,
                           HighsInt) const {}

  void checkClique(const HighsCliqueTable::CliqueVar*, HighsInt) const {}

  void checkVub(HighsInt, HighsInt, double, double) const {}

  void checkVlb(HighsInt, HighsInt, double, double) const {}

  void checkConflictReasonFrontier(
      const std::set<HighsDomain::ConflictSet::LocalDomChg>&,
      const std::vector<HighsDomainChange>&) const {}

  void checkConflictReconvergenceFrontier(
      const std::set<HighsDomain::ConflictSet::LocalDomChg>&,
      const HighsDomain::ConflictSet::LocalDomChg&,
      const std::vector<HighsDomainChange>&) const {}
};
#endif

#endif
