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

std::pair<int64_t, HighsInt> HPresolve::computeProbingScore(
    HighsInt col) const {
  HighsInt implicsUp =
      mipsolver->mipdata_->cliquetable.getNumImplications(col, 1);
  HighsInt implicsDown =
      mipsolver->mipdata_->cliquetable.getNumImplications(col, 0);
  return std::make_pair(
      std::min(int64_t{5000}, static_cast<int64_t>(implicsUp) * implicsDown) /
          (int64_t{1} + static_cast<int64_t>(numProbes[col])),
      std::min(HighsInt{100}, implicsUp + implicsDown));
}

HPresolve::Result HPresolve::runProbing(HighsPostsolveStack& postsolve_stack) {
  HighsProfilingScope clock(*mipsolver->profiling_, kMipClockProbingPresolve);
  probingEarlyAbort = false;

  HighsInt oldNumProbed = numProbed;

  // prepare probing
  bool firstCall = false;
  PresolveProbingScope probingFlag(mipsolver->mipdata_->cliquetable);
  HPRESOLVE_CHECKED_CALL(prepareProbing(postsolve_stack, firstCall));

  HighsDomain& domain = mipsolver->mipdata_->getDomain();
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;

  // store binary variables in vector with their number of implications on
  // other binaries
  std::vector<std::tuple<int64_t, HighsInt, HighsInt, HighsInt>> binaries;

  if (!cliquetable.isFull()) {
    binaries.reserve(model->num_col_);
    HighsRandom random(options->random_seed);
    for (HighsInt i = 0; i != model->num_col_; ++i) {
      if (domain.isBinary(i)) {
        auto probingScore = computeProbingScore(i);
        binaries.emplace_back(-probingScore.first, -probingScore.second,
                              random.integer(), i);
      }
    }
  }
  if (!binaries.empty()) {
    // sort variables with many implications on other binaries first
    pdqsort(binaries.begin(), binaries.end());

    size_t numChangedCols = 0;
    while (domain.getChangedCols().size() != numChangedCols) {
      if (domain.isFixed(domain.getChangedCols()[numChangedCols++]))
        ++probingNumDelCol;
    }

    HighsInt numCliquesStart = cliquetable.numCliques();
    HighsInt numImplicsStart = implications.getNumImplications();
    HighsInt numDelStart = probingNumDelCol;

    auto calcNumDel = [&]() {
      return probingNumDelCol - numDelStart +
             static_cast<HighsInt>(implications.substitutions.size() +
                                   cliquetable.getSubstitutions().size());
    };

    HighsInt numDel = calcNumDel();

    int64_t splayContingent =
        cliquetable.numNeighbourhoodQueries +
        std::max(mipsolver->submip ? HighsInt{0} : HighsInt{100000},
                 10 * numNonzeros());
    HighsInt numFail = 0;

    // lambda to check if model has enough continuous variables to perform
    // lifting for probing
    auto modelHasPercentageContVars = [&](size_t percentage) {
      size_t num_cols = 0, num_cont_cols = 0;
      for (size_t col = 0; col < colsize.size(); col++) {
        if (colDeleted[col]) continue;
        num_cols++;
        if (model->integrality_[col] == HighsVarType::kContinuous)
          num_cont_cols++;
      }
      return size_t{100} * num_cont_cols >= percentage * num_cols;
    };

    // collect up to 10 lifting opportunities per row
    const size_t maxNumLiftOpps = std::max(
        size_t{100000}, size_t{10} * static_cast<size_t>(model->num_row_));

    // only search for lifting opportunities if at least 2 percent of the
    // variables in the problem are continuous
    size_t numLiftOpps = 0;
    if (mipsolver->options_mip_->mip_lifting_for_probing != -1 &&
        modelHasPercentageContVars(size_t{2})) {
      // store lifting opportunities
      implications.storeLiftingOpportunity = [&](HighsInt row, HighsInt col,
                                                 HighsInt val, double coef) {
        if (coef == 0.0) return;
        // find lifting opportunities for row
        auto& htree = liftingOpportunities[row];
        // add element
        auto insertresult = htree.insert_or_get(std::make_pair(col, val), coef);
        assert(insertresult.second);
        numLiftOpps++;
      };
    }

    // Set up logging for probing
    const bool silent = silentLog();
    HighsInt iBin = -1;
    HighsInt iBin_probed = -1;
    HighsInt num_binary = binaries.size();
    double tt = this->timer->read();
    double tt0 = tt;
    double log_tt = tt0;
    HighsInt log_iBin_probed = iBin_probed;
    auto probingLog = [&]() {
      if (silent || options->timeless_log) return;
      // Ensure that enough time has elapsed, and at least on binary
      // has been probed
      const double log_tt_interval = 5.0;
      if (tt > log_tt + log_tt_interval && iBin_probed > log_iBin_probed) {
        // Get the average rate from the start of probing
        assert(iBin_probed > 0);
        double rate0 = (tt - tt0) / double(iBin_probed);
        // Get the average rate since last logging
        HighsInt dl_iBin_probed = iBin_probed - log_iBin_probed;
        assert(dl_iBin_probed > 0);
        double rate1 = (tt - log_tt) / double(dl_iBin_probed);
        // Assess the time for probing based on the greater rate
        double rate = std::max(rate0, rate1);
        std::string rate_str =
            " (rate " + highsTimeToString(1e3 * rate) + "/ms";
        double expected_probing_finish_time =
            tt + rate * (num_binary - iBin_probed);
        std::string expected_probing_finish_time_str =
            " => expected probing finish time " +
            highsTimeSecondToString(expected_probing_finish_time) + ")";
        std::string time_str = highsTimeSecondToString(tt);
        highsLogUser(options->log_options, HighsLogType::kInfo,
                     "   Considered %d / %d binaries; %d probed %s%s %s\n",
                     int(iBin), int(num_binary), int(iBin_probed),
                     rate_str.c_str(), expected_probing_finish_time_str.c_str(),
                     time_str.c_str());
        log_tt = tt;
        log_iBin_probed = iBin_probed;
      }
    };

    const bool dualFixProbingEnabled =
        allow_rule_[kPresolveRuleDualFixProbing] && !mipsolver->submip;
    if (dualFixProbingEnabled) {
      domain.getDualFixProbingPropagation().recomputeLocks();
    }

    for (const auto& binvar : binaries) {
      // Count the binaries considered
      iBin++;

      HighsInt i = std::get<3>(binvar);
      if (cliquetable.getSubstitution(i) != nullptr || !domain.isBinary(i))
        continue;

      // Count the binaries probed
      iBin_probed++;

      // Check for timeout
      tt = this->timer->read();
      if (tt > options->time_limit) {
        highsLogUser(
            options->log_options, HighsLogType::kInfo,
            "Time limit reached in probing: "
            "consider not using probing by setting option "
            "presolve_rule_off to 2^%-d = %d\n",
            int(kPresolveRuleProbing),
            int(std::pow(int(2), static_cast<int>(kPresolveRuleProbing))));
        return Result::kStopped;
      }

      // Possibly log probing
      probingLog();

      bool tightenLimits = (numProbed - oldNumProbed) >= 2500;

      // when a large percentage of columns have been deleted, stop this round
      // of probing
      // if (numDel > std::max(model->num_col_ * 0.2, 1000.)) break;
      if (!tightenLimits) {
        probingEarlyAbort =
            numDel >
            std::max(HighsInt{1000}, (model->num_row_ + model->num_col_) / 20);
      } else {
        probingEarlyAbort =
            numDel >
            std::min(HighsInt{1000}, (model->num_row_ + model->num_col_) / 20);
      }
      if (probingEarlyAbort) break;

      // break in case of too many new implications to not spent ages in
      // probing
      if (cliquetable.isFull() ||
          cliquetable.numCliques() - numCliquesStart >
              std::max(HighsInt{1000000}, 2 * numNonzeros()) ||
          implications.getNumImplications() - numImplicsStart >
              std::max(HighsInt{1000000}, 2 * numNonzeros()))
        break;

      // if (numProbed % 10 == 0)
      //   printf(
      //       "numprobed=%d  numDel=%d  newcliques=%d "
      //       "numNeighbourhoodQueries=%ld  "
      //       "splayContingent=%ld\n",
      //       numProbed, numDel, cliquetable.numCliques() - numCliquesStart,
      //       cliquetable.numNeighbourhoodQueries, splayContingent);
      if (cliquetable.numNeighbourhoodQueries > splayContingent) break;

      if (probingContingent - numProbed < 0) break;

      HighsInt numBoundChgs = 0;
      HighsInt numNewCliques = -cliquetable.numCliques();
      domain.setDualFixProbingActive(dualFixProbingEnabled);
      const bool probing_result = implications.runProbing(i, numBoundChgs);
      domain.setDualFixProbingActive(false);
      if (!probing_result) continue;
      probingContingent += numBoundChgs;
      numNewCliques += cliquetable.numCliques();
      numNewCliques = std::max(numNewCliques, HighsInt{0});
      while (domain.getChangedCols().size() != numChangedCols) {
        if (domain.isFixed(domain.getChangedCols()[numChangedCols++]))
          ++probingNumDelCol;
      }
      HighsInt newNumDel = calcNumDel();

      if (newNumDel > numDel) {
        probingContingent += numDel;
        if (!mipsolver->submip) {
          splayContingent += 100 * (newNumDel + numDelStart);
          splayContingent += 1000 * numNewCliques;
        }
        numDel = newNumDel;
        numFail = 0;
      } else if (mipsolver->submip || numNewCliques == 0) {
        splayContingent -= (tightenLimits ? 250 : 100) * numFail;
        ++numFail;
      } else {
        splayContingent += 1000 * numNewCliques;
        numFail = 0;
      }

      ++numProbed;
      numProbes[i] += 1;

      // Stop collecting lifting opportunities if maximum is reached
      if (numLiftOpps >= maxNumLiftOpps)
        implications.storeLiftingOpportunity = nullptr;

      if (domain.infeasible()) return Result::kPrimalInfeasible;
    }

    // finalise probing
    HighsInt numVarsFixed = 0;
    HighsInt numBndsTightened = 0;
    HighsInt numVarsSubstituted = 0;
    HighsInt liftedNonzeros = 0;
    HPRESOLVE_CHECKED_CALL(finaliseProbing(postsolve_stack, firstCall,
                                           numVarsFixed, numBndsTightened,
                                           numVarsSubstituted, liftedNonzeros));

    probingNumDelCol += numVarsSubstituted;

    highsLogDev(options->log_options, HighsLogType::kInfo,
                "%" HIGHSINT_FORMAT " probing evaluations: %" HIGHSINT_FORMAT
                " deleted rows, %" HIGHSINT_FORMAT
                " deleted "
                "columns, %" HIGHSINT_FORMAT " lifted nonzeros\n",
                numProbed - oldNumProbed, numDeletedRows, numDeletedCols,
                liftedNonzeros);

    // lifting for probing
    if (mipsolver->options_mip_->mip_lifting_for_probing != -1) {
      // only perform lifting if probing did not modify the problem so far
      if (numDeletedRows == 0 && numDeletedCols == 0 && liftedNonzeros == 0)
        HPRESOLVE_CHECKED_CALL(liftingForProbing(postsolve_stack));
      // clear lifting opportunities
      liftingOpportunities.clear();
      implications.storeLiftingOpportunity = nullptr;
    }
  }

  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::liftingForProbing(
    HighsPostsolveStack& postsolve_stack) {
  // this method implements lifting for probing as described by Achterberg et
  // al. (2019) Presolve Reductions in Mixed Integer Programming. INFORMS
  // Journal on Computing 32(2):473-506.
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  const HighsDomain& domain = mipsolver->mipdata_->getDomain();

  // collect best lifting opportunity for each row in a vector
  typedef std::pair<HighsCliqueTable::CliqueVar, double> liftingvar;
  typedef std::tuple<HighsInt, std::vector<liftingvar>, double, HighsInt>
      liftingdata;
  std::vector<liftingdata> liftingtable;
  liftingtable.reserve(liftingOpportunities.size());

  // remember overall best score
  double bestscoretotal = -kHighsInf;

  // is lifting allowed to add non-zeros?
  const bool fillallowed = mipsolver->options_mip_->mip_lifting_for_probing > 0;

  // store clique variables and coefficients in a map
  auto comp = [](const HighsCliqueTable::CliqueVar& c1,
                 const HighsCliqueTable::CliqueVar& c2) {
    return c1.col < c2.col || (c1.col == c2.col && c1.val < c2.val);
  };
  std::map<HighsCliqueTable::CliqueVar, std::pair<double, HighsInt>,
           decltype(comp)>
      coefficients(comp);

  // consider lifting opportunities
  size_t numrowsremoved = 0;
  for (const auto& elm : liftingOpportunities) {
    // get row index and skip deleted rows
    HighsInt row = elm.first;
    if (rowDeleted[row]) continue;

    // do not add non-zeros to dense rows
    const bool dense =
        rowsize[row] >
        std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);

    // iterate over elements in hash tree
    const auto& htree = elm.second;
    bool isredundant = false;
    htree.for_each([&](const std::pair<HighsInt, HighsInt>& data, double coef) {
      HighsInt col = data.first;
      HighsInt val = data.second;
      HighsInt pos = findNonzero(row, col);
      isredundant = isredundant || htree.contains(std::make_pair(col, 1 - val));
      if (!dense && (fillallowed || pos != -1) && !colDeleted[col] &&
          !domain.isFixed(col))
        coefficients[HighsCliqueTable::CliqueVar{col, val}] = {coef, pos};
    });

    // remove redundant rows
    if (isredundant) {
      numrowsremoved++;
      postsolve_stack.redundantRow(row);
      removeRow(row);
      coefficients.clear();
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
      continue;
    }

    // skip rows with empty coefficient map
    if (coefficients.empty()) continue;

    // vector to hold best clique
    std::vector<liftingvar> bestclique;
    double bestscore = -kHighsInf;
    HighsInt bestnfill = 0;

    // lambda for computing coefficient difference
    auto computeCoeffDiff = [&](double newvalue, HighsInt nzpos) {
      return std::fabs(newvalue - (nzpos == -1 ? 0 : Avalue[nzpos]));
    };

    // store candidates in a vector
    std::vector<HighsCliqueTable::CliqueVar> candidates;
    candidates.reserve(coefficients.size());
    for (const auto& elm : coefficients) {
      candidates.push_back(elm.first);
      // initialize best clique
      double score = computeCoeffDiff(elm.second.first, elm.second.second);
      if (score > bestscore) {
        bestscore = score;
        bestnfill = elm.second.second == -1 ? 1 : 0;
        bestclique = {std::make_pair(elm.first, elm.second.first)};
      }
    }

    if (candidates.size() > 1) {
      // compute cliques
      auto cliques =
          cliquetable.computeMaximalCliques(candidates, primal_feastol);

      // identify clique with highest score
      for (const auto& clique : cliques) {
        HighsCDouble score = 0;
        HighsInt nfill = 0;
        for (const auto& cliquevar : clique) {
          score += computeCoeffDiff(coefficients[cliquevar].first,
                                    coefficients[cliquevar].second);
          if (coefficients[cliquevar].second == -1) nfill++;
        }
        if (score > bestscore) {
          bestscore = static_cast<double>(score);
          bestnfill = nfill;
          bestclique.clear();
          bestclique.reserve(clique.size());
          for (const auto& cliquevar : clique) {
            bestclique.emplace_back(cliquevar, coefficients[cliquevar].first);
          }
        }
      }
    }

    // store best clique
    liftingtable.emplace_back(row, bestclique, bestscore, bestnfill);
    bestscoretotal = std::max(bestscoretotal, bestscore);
    coefficients.clear();
  }

  // lambda for computing score
  auto computeOverallScore = [&](double score, HighsInt numelms,
                                 HighsInt numfillin) {
    const double weight = 0.5;
    return weight * (score / bestscoretotal) +
           (1 - weight) * static_cast<double>(numelms - numfillin) /
               static_cast<double>(numelms);
  };

  // sort according to score
  pdqsort(
      liftingtable.begin(), liftingtable.end(),
      [&](const liftingdata& opp1, const liftingdata& opp2) {
        double score1 = computeOverallScore(
            std::get<2>(opp1), static_cast<HighsInt>(std::get<1>(opp1).size()),
            std::get<3>(opp1));
        double score2 = computeOverallScore(
            std::get<2>(opp2), static_cast<HighsInt>(std::get<1>(opp2).size()),
            std::get<3>(opp2));
        return (score1 == score2 ? std::get<0>(opp1) < std::get<0>(opp2)
                                 : score1 > score2);
      });

  // perform actual lifting
  size_t nfill = 0;
  size_t nmod = 0;
  size_t numrowsmodified = 0;
  const size_t maxnfill = std::max(10 * liftingtable.size(),
                                   static_cast<size_t>(numNonzeros()) / 100);
  for (const auto& lifting : liftingtable) {
    // get clique
    HighsInt row = std::get<0>(lifting);
    const auto& bestclique = std::get<1>(lifting);

    // check against max. fill-in
    size_t newfill = static_cast<size_t>(std::get<3>(lifting));
    if (nfill + newfill > maxnfill) break;
    nfill += newfill;
    nmod += bestclique.size() - newfill;

    // update matrix
    HighsCDouble update = 0.0;
    for (const auto& elm : bestclique) {
      // get data
      const auto& cliquevar = std::get<0>(elm);
      const double& coeff = std::get<1>(elm);
      // add non-zero to matrix
      addToMatrix(row, cliquevar.col, coeff);
      // compute term to update left-hand / right-hand side
      if (cliquevar.val == 0) update += coeff;
    }

    // update left-hand / right-hand sides
    numrowsmodified++;
    addToRowLower(row, update);
    addToRowUpper(row, update);
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Lifting for probing removed %d and modified %d row(s), added %d "
              "new and modified %d existing nonzero(s)\n",
              static_cast<int>(numrowsremoved),
              static_cast<int>(numrowsmodified), static_cast<int>(nfill),
              static_cast<int>(nmod));

  return Result::kOk;
}

}  // namespace presolve
