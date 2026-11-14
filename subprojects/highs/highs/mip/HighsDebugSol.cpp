/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsDebugSol.h"

#ifdef HIGHS_DEBUGSOL

#include <fstream>

#include "io/FilereaderMps.h"
#include "lp_data/HighsLpUtils.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"

HighsDebugSol::HighsDebugSol(HighsMipSolver& mipsolver)
    : debugSolActive(false) {
  this->mipsolver = &mipsolver;
  debugSolObjective = -kHighsInf;
  debugSolActive = false;
}

void HighsDebugSol::activate() {
  if (!mipsolver->submip &&
      debugSolObjective <= mipsolver->mipdata_->upper_limit &&
      !mipsolver->options_mip_->mip_debug_solution_file.empty()) {
    highsLogDev(mipsolver->options_mip_->log_options, HighsLogType::kInfo,
                "reading debug solution file %s\n",
                mipsolver->options_mip_->mip_debug_solution_file.c_str());
    std::ifstream file(mipsolver->options_mip_->mip_debug_solution_file);
    if (file) {
      std::string varname;
      double varval;
      std::string line;
      bool incolsection = false;
      std::map<std::string, int> nametoidx;

      for (HighsInt i = 0; i != mipsolver->orig_model_->num_col_; ++i)
        nametoidx[mipsolver->orig_model_->col_names_[i]] = i;

      debugOrigSolution.resize(mipsolver->orig_model_->num_col_, 0.0);
      while (std::getline(file, line)) {
        // Check for start and stop markers
        if (line.find("# Columns") != std::string::npos) {
          incolsection = true;
          continue;
        }
        if (line.find("# Rows") != std::string::npos) {
          break;
        }
        if (!incolsection) continue;

        std::istringstream linestream(line);
        linestream >> varname;
        auto it = nametoidx.find(varname);
        if (it != nametoidx.end()) {
          linestream >> varval;
          debugOrigSolution[it->second] = varval;
          highsLogDev(mipsolver->options_mip_->log_options, HighsLogType::kInfo,
                      "%s = %g\n", varname.c_str(), varval);
        }
      }
      debugSolution = debugOrigSolution;

      HighsCDouble debugsolobj = 0.0;
      for (HighsInt i = 0; i != mipsolver->orig_model_->num_col_; ++i)
        debugsolobj +=
            static_cast<HighsCDouble>(mipsolver->orig_model_->col_cost_[i]) *
            debugOrigSolution[i];

      debugSolObjective =
          static_cast<double>(debugsolobj + mipsolver->orig_model_->offset_);
      debugSolActive = true;
      printf("debug sol active\n");
      registerDomain(mipsolver->mipdata_->getDomain());
    } else {
      highsLogUser(mipsolver->options_mip_->log_options, HighsLogType::kWarning,
                   "debug solution: could not open file '%s'\n",
                   mipsolver->options_mip_->mip_debug_solution_file.c_str());
      HighsModel model;
      model.lp_ = *mipsolver->model_;
      model.lp_.col_names_.clear();
      model.lp_.row_names_.clear();
      model.lp_.col_lower_ = mipsolver->mipdata_->getDomain().col_lower_;
      model.lp_.col_upper_ = mipsolver->mipdata_->getDomain().col_upper_;
      FilereaderMps().writeModelToFile(*mipsolver->options_mip_,
                                       "debug_mip.mps", model);
    }
  }
}

void HighsDebugSol::shrink(const std::vector<HighsInt>& newColIndex) {
  if (!debugSolActive) return;

  HighsInt oldNumCol = debugSolution.size();
  for (HighsInt i = 0; i != oldNumCol; ++i)
    if (newColIndex[i] != -1) debugSolution[newColIndex[i]] = debugSolution[i];

  debugSolution.resize(mipsolver->model_->num_col_);
  HighsCDouble debugsolobj = 0.0;
  for (HighsInt i = 0; i != mipsolver->model_->num_col_; ++i)
    debugsolobj += static_cast<HighsCDouble>(mipsolver->model_->col_cost_[i]) *
                   debugSolution[i];

  debugSolObjective = static_cast<double>(debugsolobj);

  conflictingBounds.clear();
}

void HighsDebugSol::registerDomain(const HighsDomain& domain) {
  conflictingBounds.emplace(&domain, std::multiset<HighsDomainChange>());

  if (!debugSolActive) return;

  for (HighsInt i = 0; i != mipsolver->numCol(); ++i) {
    assert(domain.col_lower_[i] <=
           debugSolution[i] + mipsolver->mipdata_->feastol);
    assert(domain.col_upper_[i] >=
           debugSolution[i] - mipsolver->mipdata_->feastol);
  }
}

void HighsDebugSol::newIncumbentFound() {
  if (debugSolActive && debugSolObjective > mipsolver->mipdata_->upper_limit) {
    printf("debug sol inactive\n");
    debugSolActive = false;
  }
}

void HighsDebugSol::boundChangeAdded(const HighsDomain& domain,
                                     const HighsDomainChange& domchg,
                                     bool branching) {
  if (!debugSolActive) return;

  if (conflictingBounds.count(&domain) == 0) return;

  if (domchg.boundtype == HighsBoundType::kLower) {
    if (domchg.boundval <=
        debugSolution[domchg.column] + mipsolver->mipdata_->feastol)
      return;
  } else {
    if (domchg.boundval >=
        debugSolution[domchg.column] - mipsolver->mipdata_->feastol)
      return;
  }

  if (branching || !conflictingBounds[&domain].empty()) {
    conflictingBounds[&domain].insert(domchg);
    return;
  }

  assert(false);
}

void HighsDebugSol::boundChangeRemoved(const HighsDomain& domain,
                                       const HighsDomainChange& domchg) {
  if (!debugSolActive) return;

  if (conflictingBounds.count(&domain) == 0) return;

  auto i = conflictingBounds[&domain].find(domchg);
  if (i != conflictingBounds[&domain].end())
    conflictingBounds[&domain].erase(i);
}

void HighsDebugSol::checkCut(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double rhs) {
  if (!debugSolActive) return;

  HighsCDouble violation = -rhs;

  for (HighsInt i = 0; i != Rlen; ++i)
    violation += debugSolution[Rindex[i]] * Rvalue[i];

  assert(violation <= mipsolver->mipdata_->feastol);
}

void HighsDebugSol::checkRowAggregation(const HighsLp& lp,
                                        const HighsInt* Rindex,
                                        const double* Rvalue, HighsInt Rlen) {
  if (!debugSolActive) return;
  HighsCDouble violation = 0.0;

  HighsSolution dbgSol;
  dbgSol.dual_valid = false;
  dbgSol.value_valid = true;
  dbgSol.col_value = debugSolution;
  calculateRowValuesQuad(lp, dbgSol);
  for (HighsInt i = 0; i < Rlen; ++i) {
    if (Rindex[i] < lp.num_col_)
      violation += dbgSol.col_value[Rindex[i]] * Rvalue[i];
    else
      violation += dbgSol.row_value[Rindex[i] - lp.num_col_] * Rvalue[i];
  }

  double viol = fabs(double(violation));

  assert(viol <= mipsolver->mipdata_->feastol);
}

void HighsDebugSol::checkRow(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rlower, double Rupper) {
  if (!debugSolActive) return;

  HighsCDouble activity = 0;

  for (HighsInt i = 0; i != Rlen; ++i)
    activity += debugSolution[Rindex[i]] * Rvalue[i];

  assert(activity - mipsolver->mipdata_->feastol <= Rupper);
  assert(activity + mipsolver->mipdata_->feastol >= Rlower);
}

void HighsDebugSol::resetDomain(const HighsDomain& domain) {
  if (conflictingBounds.count(&domain) == 0) return;

  conflictingBounds[&domain].clear();
}

void HighsDebugSol::nodePruned(const HighsDomain& localdomain, double cutoff) {
  if (!debugSolActive || debugSolObjective > cutoff) return;

  if (conflictingBounds.count(&localdomain) == 0) return;

  assert(!conflictingBounds[&localdomain].empty());
}

void HighsDebugSol::checkClique(const HighsCliqueTable::CliqueVar* clq,
                                HighsInt clqlen) {
  if (!debugSolActive) return;

  HighsInt violation = -1;

  for (HighsInt i = 0; i != clqlen; ++i)
    violation += (HighsInt)(clq[i].weight(debugSolution) + 0.5);

  assert(violation <= 0);
}

void HighsDebugSol::checkVub(HighsInt col, HighsInt vubcol, double vubcoef,
                             double vubconstant) const {
  if (!debugSolActive || std::abs(vubcoef) == kHighsInf) return;

  assert(debugSolution[col] <= debugSolution[vubcol] * vubcoef + vubconstant +
                                   mipsolver->mipdata_->feastol);
}

void HighsDebugSol::checkVlb(HighsInt col, HighsInt vlbcol, double vlbcoef,
                             double vlbconstant) const {
  if (!debugSolActive || std::abs(vlbcoef) == kHighsInf) return;

  assert(debugSolution[col] >= debugSolution[vlbcol] * vlbcoef + vlbconstant -
                                   mipsolver->mipdata_->feastol);
}

void HighsDebugSol::checkConflictReasonFrontier(
    const std::set<HighsDomain::ConflictSet::LocalDomChg>& reasonSideFrontier,
    const std::vector<HighsDomainChange>& domchgstack) const {
  if (!debugSolActive) return;

  HighsInt numActiveBoundChgs = 0;
  for (const HighsDomain::ConflictSet::LocalDomChg& domchg :
       reasonSideFrontier) {
    HighsInt col = domchg.domchg.column;

    if (domchg.domchg.boundtype == HighsBoundType::kLower) {
      if (debugSolution[col] >=
          domchg.domchg.boundval - mipsolver->mipdata_->feastol)
        ++numActiveBoundChgs;
    } else {
      if (debugSolution[col] <=
          domchg.domchg.boundval + mipsolver->mipdata_->feastol)
        ++numActiveBoundChgs;
    }
  }

  assert(numActiveBoundChgs < (HighsInt)reasonSideFrontier.size());
}

void HighsDebugSol::checkConflictReconvergenceFrontier(
    const std::set<HighsDomain::ConflictSet::LocalDomChg>&
        reconvergenceFrontier,
    const HighsDomain::ConflictSet::LocalDomChg& reconvDomchg,
    const std::vector<HighsDomainChange>& domchgstack) const {
  if (!debugSolActive) return;

  HighsInt numActiveBoundChgs = 0;
  for (const HighsDomain::ConflictSet::LocalDomChg& domchg :
       reconvergenceFrontier) {
    HighsInt col = domchg.domchg.column;

    if (domchg.domchg.boundtype == HighsBoundType::kLower) {
      if (debugSolution[col] >= domchg.domchg.boundval) ++numActiveBoundChgs;
    } else {
      if (debugSolution[col] <= domchg.domchg.boundval) ++numActiveBoundChgs;
    }
  }

  auto reconvChg = mipsolver->mipdata_->getDomain().flip(reconvDomchg.domchg);

  if (reconvChg.boundtype == HighsBoundType::kLower) {
    if (debugSolution[reconvChg.column] >= reconvChg.boundval)
      ++numActiveBoundChgs;
  } else {
    if (debugSolution[reconvChg.column] <= reconvChg.boundval)
      ++numActiveBoundChgs;
  }

  assert(numActiveBoundChgs <= (HighsInt)reconvergenceFrontier.size());
}

#endif
