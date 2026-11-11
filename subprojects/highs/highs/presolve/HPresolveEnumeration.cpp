/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::enumerateSolutions(
    HighsPostsolveStack& postsolve_stack) {
  // enumerate all solutions for pure binary constraints with a small number of
  // variables
  HighsProfilingScope clock(*mipsolver->profiling_,
                            kMipClockEnumerationPresolve);

  // prepare probing
  bool firstCall = false;
  PresolveProbingScope probingFlag(mipsolver->mipdata_->cliquetable);
  HPRESOLVE_CHECKED_CALL(prepareProbing(postsolve_stack, firstCall));

  HighsDomain& domain = mipsolver->mipdata_->getDomain();
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;

  // candidate row for the enumeration
  struct CandidateRow {
    // scores of the row (negated, so that the best rows come first)
    double score;
    double score2;
    // random tie-breaker
    HighsInt tieBreaker;
    // row index, -1 if the row is removed
    HighsInt row;
    uint64_t signature;
    HighsInt numBinaries;
    // number of binaries minus the number of distinct signature bits
    HighsInt numCoincidences;

    // the row index is unique, so the remaining fields do not affect the order
    bool operator<(const CandidateRow& other) const {
      return std::make_tuple(score, score2, tieBreaker, row) <
             std::make_tuple(other.score, other.score2, other.tieBreaker,
                             other.row);
    }
  };

  // maximum size of a row and maximum number of rows that will be checked
  const size_t maxRowSize = 8;
  const HighsInt maxNumRowsChecked = 400;
  const size_t maxNumSolutions = 1 << maxRowSize;
  // maximum percentage of overlap
  const size_t maxPercentageRowOverlap = 50;
  // maximum number of consecutive fails
  const HighsInt maxNumFails = 6;

  // lambda for checking binary rows
  auto getBinaryRow = [&](HighsInt row, std::vector<HighsInt>& binvars,
                          size_t& numnzs) {
    numnzs = 0;
    for (HighsInt j = mipsolver->mipdata_->ARstart_[row];
         j < mipsolver->mipdata_->ARstart_[row + 1]; j++) {
      // get index
      HighsInt col = mipsolver->mipdata_->ARindex_[j];
      // skip fixed variables
      if (domain.isFixed(col)) continue;
      // skip row if there are non-binary variables or maximum number of
      // elements is reached
      if (!domain.isBinary(col) || numnzs >= maxRowSize) return false;
      // store binary variable
      binvars[numnzs++] = col;
    }
    if (numnzs == 0) return false;
    pdqsort(binvars.begin(), binvars.begin() + numnzs);
    return true;
  };

  // lambda for computing row score
  auto computeRowScore = [&](const std::vector<HighsInt>& binvars,
                             size_t numnzs) {
    int64_t score = 0;
    HighsInt score2 = 0;
    for (size_t i = 0; i < numnzs; i++) {
      auto probingScore = computeProbingScore(binvars[i]);
      score += probingScore.first;
      score2 += probingScore.second;
    }
    return std::make_pair<double, double>(score / static_cast<double>(numnzs),
                                          score2 / static_cast<double>(numnzs));
  };

  // lambda for computing row signature
  auto computeRowSignature = [&](const std::vector<HighsInt>& binvars,
                                 size_t numnzs) {
    uint64_t signature = 0;
    for (size_t i = 0; i < numnzs; i++) {
      HighsInt colHashedPos = (HighsHashHelpers::hash(binvars[i]) >> 58);
      assert(colHashedPos < 64);
      signature |= uint64_t{1} << colHashedPos;
    }
    return signature;
  };

  // lambda for computing the smallest overlap for which a row with numnzs
  // binaries is too similar to a row that is at least as long (see
  // removeSimilarRows)
  auto requiredRowOverlap = [&](size_t numnzs) {
    return ((maxPercentageRowOverlap + 1) * numnzs + 99) / 100;
  };

  // lambda for checking with the signatures whether two rows may overlap too
  // much. a column in both rows sets the same bit in both signatures, so the
  // common bits are at least the overlap minus the bits lost to coinciding
  // columns in either row. the check never rejects similar rows; whether
  // they are similar is decided by computeRowOverlap
  auto mayBeSimilar = [&](const CandidateRow& row, const CandidateRow& row2) {
    HighsInt numCoincidences =
        std::min(row.numCoincidences, row2.numCoincidences);
    HighsInt numCommonBits =
        HighsHashHelpers::popcnt(row.signature & row2.signature);
    return static_cast<size_t>(numCommonBits + numCoincidences) >=
           requiredRowOverlap(static_cast<size_t>(
               std::min(row.numBinaries, row2.numBinaries)));
  };

  // lambda for computing overlap of two binary rows
  auto computeRowOverlap = [&](const std::vector<HighsInt>& binvars,
                               const std::vector<HighsInt>& binvars2,
                               size_t numnzs, size_t numnzs2) {
    size_t overlap = 0;
    size_t ir = 0;
    size_t ir2 = 0;
    while (ir < numnzs && ir2 < numnzs2) {
      if (binvars[ir] < binvars2[ir2])
        ir++;
      else if (binvars[ir] > binvars2[ir2])
        ir2++;
      else {
        ir++;
        ir2++;
        overlap++;
      }
    }
    return overlap;
  };

  // lambda for compiling candidate rows
  auto compileRows = [&](std::vector<CandidateRow>& rows) {
    std::vector<HighsInt> binvars(maxRowSize);
    size_t numnzs = 0;
    HighsRandom random(options->random_seed);
    for (HighsInt i = 0; i < mipsolver->numRow(); i++) {
      // skip redundant rows
      if (domain.isRedundantRow(i)) continue;
      // skip non-binary rows
      if (!getBinaryRow(i, binvars, numnzs)) continue;
      // compute row score
      auto score = computeRowScore(binvars, numnzs);
      // compute signature and add row to vector
      uint64_t signature = computeRowSignature(binvars, numnzs);
      HighsInt numCoincidences =
          static_cast<HighsInt>(numnzs) - HighsHashHelpers::popcnt(signature);
      rows.push_back(CandidateRow{-score.first, -score.second, random.integer(),
                                  i, signature, static_cast<HighsInt>(numnzs),
                                  numCoincidences});
    }
  };

  // lambda for removing similar rows
  auto removeSimilarRows = [&](std::vector<CandidateRow>& rows) {
    if (rows.size() <= 1) return;
    std::vector<HighsInt> binvars(maxRowSize);
    std::vector<HighsInt> binvars2(maxRowSize);
    size_t numnzs;
    size_t numnzs2;
    HighsInt numRowsAccepted = 0;
    HighsInt numRowsRemoved = 0;
    size_t numComparisons = 0;
    for (size_t i = 0; i < rows.size() - 1; i++) {
      // get row index and skip removed rows
      HighsInt r = rows[i].row;
      if (r == -1) continue;
      // check if maximum number of rows is reached
      if ((++numRowsAccepted) >= maxNumRowsChecked) break;
      // get indices of binary variables in the row
      getBinaryRow(r, binvars, numnzs);
      // initialise counters
      HighsInt numRowsActive = 0;
      HighsInt oldNumRowsRemoved = numRowsRemoved;
      for (size_t ii = i + 1; ii < rows.size(); ii++) {
        // get row index and skip removed rows
        HighsInt& r2 = rows[ii].row;
        if (r2 == -1) continue;
        // compare signatures to see if there may be too much overlap
        numRowsActive++;
        if (!mayBeSimilar(rows[i], rows[ii])) continue;
        // get indices of binary variables in the row
        getBinaryRow(r2, binvars2, numnzs2);
        // check if there is too much overlap (compute overlap coefficient)
        numComparisons++;
        size_t overlap = computeRowOverlap(binvars, binvars2, numnzs, numnzs2);
        if ((100 * overlap) / std::min(numnzs, numnzs2) >
            maxPercentageRowOverlap) {
          // mark row for removal
          numRowsRemoved++;
          r2 = -1;
        }
      }
      // stop iterating if at most one of the remaining rows is active (not
      // deleted)
      if (numRowsActive - numRowsRemoved + oldNumRowsRemoved <= 1) break;
    }
    if (numRowsRemoved > 0)
      rows.erase(
          std::remove_if(rows.begin(), rows.end(),
                         [](const CandidateRow& p) { return p.row == -1; }),
          rows.end());
    if (rows.size() > static_cast<size_t>(maxNumRowsChecked))
      rows.resize(maxNumRowsChecked);
  };

  // vector of rows
  std::vector<CandidateRow> rows;
  rows.reserve(model->num_row_);

  // compile rows and sort them
  compileRows(rows);
  pdqsort(rows.begin(), rows.end());

  // remove similar rows
  removeSimilarRows(rows);

  // vectors for storing branching decisions and solutions
  struct branch {
    size_t numDomainChanges;
    size_t numChangedCols;
  };
  std::vector<std::array<HighsInt, maxNumSolutions>> solutions(maxRowSize);
  std::vector<HighsInt> vars(maxRowSize);
  std::vector<branch> branches(maxRowSize);
  std::vector<HighsInt> worstCaseBounds(model->num_col_);
  std::vector<double> worstCaseLowerBound(model->num_col_, kHighsInf);
  std::vector<double> worstCaseUpperBound(model->num_col_, -kHighsInf);
  std::vector<double> col_lower(domain.col_lower_);
  std::vector<double> col_upper(domain.col_upper_);

  // lambda for finding a variable to branch on
  auto findBranchVar = [&](size_t numVars) {
    // find variable for branching
    for (size_t i = 0; i < numVars; i++)
      if (!domain.isFixed(vars[i])) return vars[i];
    return HighsInt{-1};
  };

  // lambda for branching (just performs initial lower branch)
  auto doBranch = [&](size_t numVars, HighsInt& numBranches) {
    // find variable for branching
    HighsInt branchvar = findBranchVar(numVars);
    assert(branchvar >= 0);

    // branch downwards
    branches[++numBranches] = {domain.getDomainChangeStack().size(),
                               domain.getChangedCols().size()};
    domain.changeBound(HighsBoundType::kUpper, branchvar, 0);
  };

  // lambda for backtracking
  auto doBacktrack = [&](HighsInt& numBranches) {
    while (numBranches >= 0) {
      // get column index
      const auto& domchg =
          domain.getDomainChangeStack()[branches[numBranches].numDomainChanges];
      HighsInt col = domchg.column;
      HighsBoundType bndtype = domchg.boundtype;
      // backtrack
      domain.backtrack();
      domain.clearChangedCols(
          static_cast<HighsInt>(branches[numBranches].numChangedCols));
      if (bndtype == HighsBoundType::kUpper) {
        // branch upwards
        domain.changeBound(HighsBoundType::kLower, col, 1);
        break;
      } else {
        // remove branch
        branches[numBranches--] = {0, 0};
      }
    }
    // check if enumeration is complete
    return (numBranches >= 0);
  };

  // lambda for checking whether the values of two binary variables are
  // identical in all feasible solutions
  auto identicalVars = [&](size_t numSolutions, size_t index1, size_t index2) {
    for (size_t sol = 0; sol < numSolutions; sol++) {
      if (solutions[index1][sol] != solutions[index2][sol]) return false;
    }
    return true;
  };

  // lambda for checking whether the values of two binary variables are
  // complementary in all feasible solutions
  auto complementaryVars = [&](size_t numSolutions, size_t index1,
                               size_t index2) {
    for (size_t sol = 0; sol < numSolutions; sol++) {
      if (solutions[index1][sol] != 1 - solutions[index2][sol]) return false;
    }
    return true;
  };

  auto handleInfeasibility = [&](bool infeasible) {
    if (infeasible) return Result::kPrimalInfeasible;
    return Result::kOk;
  };

  auto removeWorstCaseBounds = [&](size_t pos, size_t& numWorstCaseBounds) {
    worstCaseLowerBound[worstCaseBounds[pos]] = kHighsInf;
    worstCaseUpperBound[worstCaseBounds[pos]] = -kHighsInf;
    worstCaseBounds[pos] = worstCaseBounds[numWorstCaseBounds - 1];
    worstCaseBounds[numWorstCaseBounds - 1] = 0;
    numWorstCaseBounds--;
  };

  auto updateWorstCaseBounds = [&](HighsInt col) {
    // update worst-case bounds
    worstCaseLowerBound[col] =
        std::min(worstCaseLowerBound[col], domain.col_lower_[col]);
    worstCaseUpperBound[col] =
        std::max(worstCaseUpperBound[col], domain.col_upper_[col]);
    // check if worst-case bounds are not tighter than global bounds
    return (worstCaseLowerBound[col] <= col_lower[col] &&
            worstCaseUpperBound[col] >= col_upper[col]);
  };

  auto handleSolution =
      [&](size_t numVars, size_t& numSolutions, size_t& numWorstCaseBounds,
          size_t& minNumActiveCols, size_t& maxNumActiveCols) {
        // propagate
        domain.propagate();
        if (domain.infeasible()) return;
        // handling of worst-case bounds
        if (numSolutions == 0) {
          // initialize
          for (HighsInt col : domain.getChangedCols()) {
            worstCaseBounds[numWorstCaseBounds++] = col;
            updateWorstCaseBounds(col);
          }
        } else {
          size_t i = 0;
          while (i < numWorstCaseBounds) {
            HighsInt col = worstCaseBounds[i];
            if (!domain.isChangedCol(col)) {
              // no bound changes for this variable -> reset worst-case
              // bounds and remove variable
              removeWorstCaseBounds(i, numWorstCaseBounds);
            } else {
              // update worst-case bounds
              if (updateWorstCaseBounds(col))
                removeWorstCaseBounds(i, numWorstCaseBounds);
              else
                i++;
            }
          }
        }
        // store solution and compute minimum and maximum number of active
        // variables (i.e. those having solution value of 1)
        size_t numActiveCols = 0;
        for (size_t i = 0; i < numVars; i++) {
          HighsInt solValue = domain.col_lower_[vars[i]] == 0.0 ? 0 : 1;
          solutions[i][numSolutions] = solValue;
          if (solValue != 0) numActiveCols++;
        }
        minNumActiveCols = std::min(minNumActiveCols, numActiveCols);
        maxNumActiveCols = std::max(maxNumActiveCols, numActiveCols);
        numSolutions++;
      };

  // loop over candidate rows
  HighsInt numCliquesFound = 0;
  HighsInt numFails = 0;
  for (const auto& r : rows) {
    // get row index
    HighsInt row = r.row;
    // skip redundant rows
    if (domain.isRedundantRow(row)) continue;
    // check row
    size_t numVars = 0;
    if (!getBinaryRow(row, vars, numVars)) continue;

    // main loop
    HighsInt numBranches = -1;
    size_t numWorstCaseBounds = 0;
    size_t numSolutions = 0;
    size_t minNumActiveCols = numVars;
    size_t maxNumActiveCols = 0;
    while (true) {
      bool backtrack = domain.infeasible();
      if (!backtrack) {
        backtrack = findBranchVar(numVars) < 0;
        if (backtrack)
          handleSolution(numVars, numSolutions, numWorstCaseBounds,
                         minNumActiveCols, maxNumActiveCols);
      }
      // branch or backtrack
      if (!backtrack)
        doBranch(numVars, numBranches);
      else if (!doBacktrack(numBranches))
        break;
    }

    // no solutions -> infeasible
    HPRESOLVE_CHECKED_CALL(handleInfeasibility(numSolutions == 0));

    // store current number of bound changes etc.
    size_t oldNumChangedCols = domain.getChangedCols().size();
    HighsInt oldNumCliques = cliquetable.numCliques();
    size_t oldNumSubstitutions = cliquetable.getSubstitutions().size();

    // check if all variables form a clique
    if (maxNumActiveCols == 1 || minNumActiveCols == numVars - 1) {
      numCliquesFound++;
      std::vector<HighsCliqueTable::CliqueVar> clique(numVars);
      for (size_t i = 0; i < numVars; i++)
        clique[i] =
            HighsCliqueTable::CliqueVar(vars[i], maxNumActiveCols == 1 ? 1 : 0);
      cliquetable.addClique(*mipsolver, clique.data(),
                            static_cast<HighsInt>(numVars),
                            minNumActiveCols == maxNumActiveCols);
      HPRESOLVE_CHECKED_CALL(handleInfeasibility(domain.infeasible()));
    }

    // analyse worst-case bounds
    for (size_t i = 0; i < numWorstCaseBounds; i++) {
      HighsInt col = worstCaseBounds[i];
      if (worstCaseLowerBound[col] > domain.col_lower_[col]) {
        // tighten lower bound
        domain.changeBound(HighsBoundType::kLower, col,
                           worstCaseLowerBound[col],
                           HighsDomain::Reason::unspecified());
        HPRESOLVE_CHECKED_CALL(handleInfeasibility(domain.infeasible()));
      }
      if (worstCaseUpperBound[col] < domain.col_upper_[col]) {
        // tighten upper bound
        domain.changeBound(HighsBoundType::kUpper, col,
                           worstCaseUpperBound[col],
                           HighsDomain::Reason::unspecified());
        HPRESOLVE_CHECKED_CALL(handleInfeasibility(domain.infeasible()));
      }
      // clean up
      worstCaseLowerBound[col] = kHighsInf;
      worstCaseUpperBound[col] = -kHighsInf;
      worstCaseBounds[i] = 0;
    }

    for (size_t i = 0; i < numVars - 1; i++) {
      // get column index
      HighsInt col = vars[i];
      for (size_t ii = i + 1; ii < numVars; ii++) {
        // get column index
        HighsInt col2 = vars[ii];
        // check if two binary variables take identical or complementary
        // values in all feasible solutions
        if (identicalVars(numSolutions, i, ii)) {
          // add clique (1 - x_1) + x_2 = 1 to clique table
          numCliquesFound++;
          std::array<HighsCliqueTable::CliqueVar, 2> clique;
          clique[0] = HighsCliqueTable::CliqueVar(col, 0);
          clique[1] = HighsCliqueTable::CliqueVar(col2, 1);
          cliquetable.addClique(*mipsolver, clique.data(), 2, true);
          HPRESOLVE_CHECKED_CALL(handleInfeasibility(domain.infeasible()));
        } else if (complementaryVars(numSolutions, i, ii)) {
          // add clique (1 - x_1) + (1 - x_2) = 1 to clique table
          numCliquesFound++;
          std::array<HighsCliqueTable::CliqueVar, 2> clique;
          clique[0] = HighsCliqueTable::CliqueVar(col, 0);
          clique[1] = HighsCliqueTable::CliqueVar(col2, 0);
          cliquetable.addClique(*mipsolver, clique.data(), 2, true);
          HPRESOLVE_CHECKED_CALL(handleInfeasibility(domain.infeasible()));
        }
      }
    }

    // update bounds
    size_t numChangedCols = domain.getChangedCols().size();
    for (HighsInt col : domain.getChangedCols()) {
      col_lower[col] = domain.col_lower_[col];
      col_upper[col] = domain.col_upper_[col];
    }

    // clear changed cols
    domain.clearChangedCols();

    // check if solution enumeration failed
    if (numChangedCols != oldNumChangedCols ||
        cliquetable.numCliques() != oldNumCliques ||
        cliquetable.getSubstitutions().size() != oldNumSubstitutions)
      numFails = 0;
    else {
      numFails++;
      if (numFails > maxNumFails) break;
    }
  }

  // finalise probing
  HighsInt numVarsFixed = 0;
  HighsInt numBndsTightened = 0;
  HighsInt numVarsSubstituted = 0;
  HighsInt liftedNonzeros = 0;
  HPRESOLVE_CHECKED_CALL(finaliseProbing(postsolve_stack, firstCall,
                                         numVarsFixed, numBndsTightened,
                                         numVarsSubstituted, liftedNonzeros));

  if (numVarsFixed > 0 || numBndsTightened > 0 || numVarsSubstituted > 0 ||
      numDeletedRows > 0 || liftedNonzeros > 0)
    highsLogDev(
        options->log_options, HighsLogType::kInfo,
        "Enumeration presolve: %d deleted columns, %d tightened bounds, "
        "%d substitutions, %d deleted rows, "
        "%d lifted nonzeros\n",
        static_cast<int>(numVarsFixed), static_cast<int>(numBndsTightened),
        static_cast<int>(numVarsSubstituted), static_cast<int>(numDeletedRows),
        static_cast<int>(liftedNonzeros));

  return checkLimits(postsolve_stack);
}

}  // namespace presolve
