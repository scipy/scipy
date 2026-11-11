/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsCutPool.h"

#include <cassert>
#include <numeric>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"

static uint64_t compute_cut_hash(const HighsInt* Rindex, const double* Rvalue,
                                 double maxabscoef, const HighsInt Rlen) {
  std::vector<uint32_t> valueHashCodes(Rlen);

  double scale = 1.0 / maxabscoef;
  for (HighsInt i = 0; i < Rlen; ++i)
    valueHashCodes[i] = HighsHashHelpers::double_hash_code(scale * Rvalue[i]);

  return HighsHashHelpers::vector_hash(Rindex, Rlen) ^
         (HighsHashHelpers::vector_hash(valueHashCodes.data(), Rlen) >> 32);
}

#if 0
static void printCut(const HighsInt* Rindex, const double* Rvalue, HighsInt Rlen,
                     double rhs) {
  for (HighsInt i = 0; i != Rlen; ++i) {
    if (Rvalue[i] > 0)
      printf("+%g<x%" HIGHSINT_FORMAT "> ", Rvalue[i], Rindex[i]);
    else
      printf("-%g<x%" HIGHSINT_FORMAT "> ", -Rvalue[i], Rindex[i]);
  }

  printf("<= %g\n", rhs);
}
#endif

bool HighsCutPool::isDuplicate(uint64_t hash, double norm,
                               const HighsInt* Rindex, const double* Rvalue,
                               HighsInt Rlen, double rhs) {
  auto range = hashToCutMap.equal_range(hash);
  const double* ARvalue = matrix_.getARvalue();
  const HighsInt* ARindex = matrix_.getARindex();

  for (auto it = range.first; it != range.second; ++it) {
    HighsInt rowindex = it->second;
    HighsInt start = matrix_.getRowStart(rowindex);
    HighsInt end = matrix_.getRowEnd(rowindex);

    if (end - start != Rlen) continue;
    if (std::equal(Rindex, Rindex + Rlen, &ARindex[start])) {
      double dotprod = 0.0;

      for (HighsInt i = 0; i != Rlen; ++i)
        dotprod += Rvalue[i] * ARvalue[start + i];

      double parallelism = double(dotprod) * rownormalization_[rowindex] * norm;

      // printf("\n\ncuts with same support and parallelism %g:\n",
      // parallelism); printf("CUT1: "); printCut(Rindex, Rvalue, Rlen, rhs);
      // printf("CUT2: ");
      // printCut(Rindex, ARvalue + start, Rlen, rhs_[rowindex]);
      // printf("\n");

      if (parallelism >= 1 - 1e-6) return true;

      //{
      //  if (ages_[rowindex] >= 0) {
      //    matrix_.replaceRowValues(rowindex, Rvalue);
      //    return rowindex;
      //  } else
      //    return -2;
      //}
    }
  }

  return false;
}

double HighsCutPool::getParallelism(HighsInt row1, HighsInt row2) const {
  HighsInt i1 = matrix_.getRowStart(row1);
  const HighsInt end1 = matrix_.getRowEnd(row1);

  HighsInt i2 = matrix_.getRowStart(row2);
  const HighsInt end2 = matrix_.getRowEnd(row2);

  const HighsInt* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();

  double dotprod = 0.0;
  while (i1 != end1 && i2 != end2) {
    HighsInt col1 = ARindex[i1];
    HighsInt col2 = ARindex[i2];

    if (col1 < col2)
      ++i1;
    else if (col2 < col1)
      ++i2;
    else {
      dotprod += ARvalue[i1] * ARvalue[i2];
      ++i1;
      ++i2;
    }
  }

  return dotprod * rownormalization_[row1] * rownormalization_[row2];
}

double HighsCutPool::getParallelism(HighsInt row1, HighsInt row2,
                                    const HighsCutPool& pool2) const {
  HighsInt i1 = matrix_.getRowStart(row1);
  const HighsInt end1 = matrix_.getRowEnd(row1);

  HighsInt i2 = pool2.matrix_.getRowStart(row2);
  const HighsInt end2 = pool2.matrix_.getRowEnd(row2);

  const HighsInt* ARindex1 = matrix_.getARindex();
  const double* ARvalue1 = matrix_.getARvalue();
  const HighsInt* ARindex2 = pool2.matrix_.getARindex();
  const double* ARvalue2 = pool2.matrix_.getARvalue();

  double dotprod = 0.0;
  while (i1 != end1 && i2 != end2) {
    HighsInt col1 = ARindex1[i1];
    HighsInt col2 = ARindex2[i2];

    if (col1 < col2)
      ++i1;
    else if (col2 < col1)
      ++i2;
    else {
      dotprod += ARvalue1[i1] * ARvalue2[i2];
      ++i1;
      ++i2;
    }
  }

  return dotprod * rownormalization_[row1] * pool2.rownormalization_[row2];
}

void HighsCutPool::lpCutRemoved(HighsInt cut, bool thread_safe) {
  const HighsInt n = numLps_[cut].fetch_add(-1, std::memory_order_relaxed);
  if (thread_safe || n > 1) return;
  if (matrix_.columnsLinked(cut)) {
    propRows.erase(std::make_pair(ages_[cut], cut));
    propRows.emplace(1, cut);
  }
  ages_[cut] = 1;
  --numLpCuts;
  ++ageDistribution[1];
}

void HighsCutPool::performAging() {
  HighsInt cutIndexEnd = matrix_.getNumRows();

  HighsInt agelim = agelim_;
  HighsInt numAvailableCuts = getNumAvailableCuts();
  while (agelim > 5 && numAvailableCuts > softlimit_) {
    numAvailableCuts -= ageDistribution[agelim];
    --agelim;
  }

  for (HighsInt i = 0; i != cutIndexEnd; ++i) {
    // Catch buffered changes (should only occur in parallel case)
    if (numLps_[i] > 0 && ages_[i] >= 0) {
      // Cut has been added to the LP, but age changes haven't been made
      --ageDistribution[ages_[i]];
      if (matrix_.columnsLinked(i)) {
        propRows.erase(std::make_pair(ages_[i], i));
        propRows.emplace(-1, i);
      }
      ages_[i] = -1;
      ++numLpCuts;
      ageResetWhileLocked_[i].store(0, std::memory_order_relaxed);
    } else if (numLps_[i] == 0 && ages_[i] == -1 && rhs_[i] != kHighsInf) {
      // Cut was removed from the LP, but age changes haven't been made
      if (matrix_.columnsLinked(i)) {
        propRows.erase(std::make_pair(ages_[i], i));
        propRows.emplace(1, i);
      }
      ages_[i] = 1;
      --numLpCuts;
      ++ageDistribution[1];
      ageResetWhileLocked_[i].store(0, std::memory_order_relaxed);
      continue;
    } else if (ageResetWhileLocked_[i].load(std::memory_order_relaxed) == 1) {
      resetAge(i);
    }
    ageResetWhileLocked_[i].store(0, std::memory_order_relaxed);
    if (ages_[i] < 0) continue;

    bool isPropagated = matrix_.columnsLinked(i);
    if (isPropagated) propRows.erase(std::make_pair(ages_[i], i));
    ageDistribution[ages_[i]] -= 1;
    ages_[i] += 1;

    if (ages_[i] > agelim) {
      for (HighsDomain::CutpoolPropagation* propagationdomain :
           propagationDomains)
        propagationdomain->cutDeleted(i);

      if (isPropagated) {
        --numPropRows;
        numPropNzs -= getRowLength(i);
      }

      matrix_.removeRow(i);
      ages_[i] = -1;
      rhs_[i] = kHighsInf;
      hasSynced_[i] = false;
    } else {
      if (isPropagated) propRows.emplace(ages_[i], i);
      ageDistribution[ages_[i]] += 1;
    }
  }

  assert((HighsInt)propRows.size() == numPropRows);
}

void HighsCutPool::separate(const std::vector<double>& sol,
                            const HighsDomain& domain, HighsCutSet& cutset,
                            double feastol,
                            const std::deque<HighsCutPool>& cutpools,
                            bool thread_safe) {
  HighsInt nrows = matrix_.getNumRows();
  const HighsInt* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();

  std::vector<std::pair<double, HighsInt>> efficacious_cuts;

  HighsInt agelim = agelim_;

  HighsInt numCuts = getNumCuts() - numLpCuts;
  while (agelim > 1 && numCuts > softlimit_) {
    numCuts -= ageDistribution[agelim];
    --agelim;
  }

  for (HighsInt i = 0; i < nrows; ++i) {
    // cuts with an age of -1 are already in the LP and are therefore skipped
    // Warning: Parallel case tries to add cuts already in current LP.
    // Inefficient. Not sure what happens if added twice.
    // The cut shouldn't have enough violation to be added though.
    if (ages_[i] < 0) continue;

    HighsInt start = matrix_.getRowStart(i);
    HighsInt end = matrix_.getRowEnd(i);

    double viol(-rhs_[i]);

    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double solval = sol[col];

      viol += ARvalue[j] * solval;
    }

    // if the cut is not violated more than feasibility tolerance
    // we skip it and increase its age, otherwise we reset its age
    bool isPropagated = matrix_.columnsLinked(i);
    if (!thread_safe) {
      ageDistribution[ages_[i]] -= 1;
      if (isPropagated) {
        propRows.erase(std::make_pair(ages_[i], i));
      }
    }
    if (double(viol) <= feastol) {
      if (thread_safe) continue;
      ++ages_[i];
      if (ages_[i] >= agelim) {
        uint64_t h = compute_cut_hash(&ARindex[start], &ARvalue[start],
                                      maxabscoef_[i], end - start);

        for (HighsDomain::CutpoolPropagation* propagationdomain :
             propagationDomains)
          propagationdomain->cutDeleted(i);

        if (isPropagated) {
          --numPropRows;
          numPropNzs -= getRowLength(i);
        }

        matrix_.removeRow(i);
        ages_[i] = -1;
        rhs_[i] = kHighsInf;
        ageResetWhileLocked_[i].store(0, std::memory_order_relaxed);
        hasSynced_[i] = false;
        auto range = hashToCutMap.equal_range(h);

        for (auto it = range.first; it != range.second; ++it) {
          if (it->second == i) {
            hashToCutMap.erase(it);
            break;
          }
        }
      } else {
        if (isPropagated) propRows.emplace(ages_[i], i);
        ageDistribution[ages_[i]] += 1;
      }
      continue;
    }

    // compute the norm only for those entries that do not sit at their minimal
    // activity in the current solution this avoids the phenomenon that the
    // traditional efficacy gets weaker for stronger cuts E.g. when considering
    // a clique cut which has additional entries whose value in the current
    // solution is 0 then the efficacy gets lower for each such entry even
    // though the cut dominates the clique cut where all those entries are
    // relaxed out.
    HighsCDouble rownorm = 0.0;
    HighsInt numActiveNzs = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double solval = sol[col];
      if (ARvalue[j] > 0) {
        if (solval > domain.col_lower_[col] + feastol) {
          rownorm += ARvalue[j] * ARvalue[j];
          numActiveNzs += 1;
        }
      } else {
        if (solval < domain.col_upper_[col] - feastol) {
          rownorm += ARvalue[j] * ARvalue[j];
          numActiveNzs += 1;
        }
      }
    }

    if (!thread_safe) {
      ages_[i] = 0;
      ++ageDistribution[0];
      if (isPropagated) propRows.emplace(ages_[i], i);
    }
    double score = viol / (numActiveNzs * sqrt(double(rownorm)));

    efficacious_cuts.emplace_back(score, i);
  }
  assert((HighsInt)propRows.size() == numPropRows);
  if (efficacious_cuts.empty()) return;

  pdqsort(efficacious_cuts.begin(), efficacious_cuts.end(),
          [&efficacious_cuts](const std::pair<double, HighsInt>& a,
                              const std::pair<double, HighsInt>& b) {
            if (a.first > b.first) return true;
            if (a.first < b.first) return false;
            return std::make_pair(
                       HighsHashHelpers::hash((uint64_t(a.second) << 32) +
                                              efficacious_cuts.size()),
                       a.second) >
                   std::make_pair(
                       HighsHashHelpers::hash((uint64_t(b.second) << 32) +
                                              efficacious_cuts.size()),
                       b.second);
          });

  double bestObservedScoreCopy = bestObservedScore;
  double& bestObservedScore_ =
      thread_safe ? bestObservedScoreCopy : bestObservedScore;
  bestObservedScore_ = std::max(efficacious_cuts[0].first, bestObservedScore);
  double minScoreFactorCopy = minScoreFactor;
  double& minScoreFactor_ = thread_safe ? minScoreFactorCopy : minScoreFactor;
  double minScore = minScoreFactor_ * bestObservedScore_;

  HighsInt numefficacious =
      std::upper_bound(efficacious_cuts.begin(), efficacious_cuts.end(),
                       minScore,
                       [](double mscore, std::pair<double, HighsInt> const& c) {
                         return mscore > c.first;
                       }) -
      efficacious_cuts.begin();

  HighsInt lowerThreshold = efficacious_cuts.size() / 20;
  HighsInt upperThreshold = efficacious_cuts.size() - 1;

  if (numefficacious <= lowerThreshold) {
    numefficacious = std::max(efficacious_cuts.size() / 2, size_t{1});
    minScoreFactor_ =
        efficacious_cuts[numefficacious - 1].first / bestObservedScore_;
  } else if (numefficacious > upperThreshold) {
    minScoreFactor_ =
        efficacious_cuts[upperThreshold].first / bestObservedScore_;
  }

  efficacious_cuts.resize(numefficacious);

  HighsInt orignumcuts = cutset.numCuts();
  HighsInt origselectednnz = cutset.ARindex_.size();
  HighsInt selectednnz = origselectednnz;

  for (const std::pair<double, HighsInt>& p : efficacious_cuts) {
    bool discard = false;
    double maxpar = 0.1;
    for (HighsInt i = 0; i != static_cast<HighsInt>(cutset.cutindices.size());
         ++i) {
      if (cutset.cutpools[i] == index_) {
        if (getParallelism(cutset.cutindices[i], p.second) > maxpar) {
          discard = true;
          break;
        }
      } else {
        // Warning: This assumes the cuts in the pool are not changing during,
        // this query, i.e., the worker's pool and the global pool.
        // Currently safe, but doesn't generalise to all designs.
        if (getParallelism(p.second, cutset.cutindices[i],
                           cutpools[cutset.cutpools[i]]) > maxpar) {
          discard = true;
          break;
        }
      }
    }

    if (discard) continue;

    numLps_[p.second].fetch_add(1, std::memory_order_relaxed);
    if (!thread_safe) {
      --ageDistribution[ages_[p.second]];
      ++numLpCuts;
      if (matrix_.columnsLinked(p.second)) {
        propRows.erase(std::make_pair(ages_[p.second], p.second));
        propRows.emplace(-1, p.second);
      }
      ages_[p.second] = -1;
    }
    cutset.cutindices.push_back(p.second);
    cutset.cutpools.push_back(index_);
    selectednnz += matrix_.getRowEnd(p.second) - matrix_.getRowStart(p.second);
  }

  cutset.resize(selectednnz);

  assert(int(cutset.ARvalue_.size()) == selectednnz);
  assert(int(cutset.ARindex_.size()) == selectednnz);

  HighsInt offset = origselectednnz;
  for (HighsInt i = orignumcuts; i != cutset.numCuts(); ++i) {
    cutset.ARstart_[i] = offset;
    HighsInt cut = cutset.cutindices[i];
    HighsInt start = matrix_.getRowStart(cut);
    HighsInt end = matrix_.getRowEnd(cut);
    cutset.upper_[i] = rhs_[cut];

    for (HighsInt j = start; j != end; ++j) {
      assert(offset < selectednnz);
      cutset.ARvalue_[offset] = ARvalue[j];
      cutset.ARindex_[offset] = ARindex[j];
      ++offset;
    }
  }

  assert((HighsInt)propRows.size() == numPropRows);
  cutset.ARstart_[cutset.numCuts()] = offset;
}

void HighsCutPool::separateLpCutsAfterRestart(HighsCutSet& cutset) {
  // should only be called after a restart with a fresh row matrix right now
  assert(matrix_.getNumDelRows() == 0);
  HighsInt numcuts = matrix_.getNumRows();

  cutset.cutindices.resize(numcuts);
  cutset.cutpools.resize(numcuts, index_);
  std::iota(cutset.cutindices.begin(), cutset.cutindices.end(), 0);
  cutset.resize(matrix_.nonzeroCapacity());

  HighsInt offset = 0;
  const HighsInt* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();
  for (HighsInt i = 0; i != cutset.numCuts(); ++i) {
    --ageDistribution[ages_[i]];
    ++numLpCuts;
    if (matrix_.columnsLinked(i)) {
      propRows.erase(std::make_pair(ages_[i], i));
      propRows.emplace(-1, i);
    }
    numLps_[i] = 1;
    ages_[i] = -1;
    cutset.ARstart_[i] = offset;
    HighsInt cut = cutset.cutindices[i];
    HighsInt start = matrix_.getRowStart(cut);
    HighsInt end = matrix_.getRowEnd(cut);
    cutset.upper_[i] = rhs_[cut];

    for (HighsInt j = start; j != end; ++j) {
      assert(offset < (HighsInt)matrix_.nonzeroCapacity());
      cutset.ARvalue_[offset] = ARvalue[j];
      cutset.ARindex_[offset] = ARindex[j];
      ++offset;
    }
  }

  cutset.ARstart_[cutset.numCuts()] = offset;

  assert((HighsInt)propRows.size() == numPropRows);
}

HighsInt HighsCutPool::addCut(const HighsMipSolver& mipsolver, HighsInt* Rindex,
                              double* Rvalue, HighsInt Rlen, double rhs,
                              bool integral, bool propagate,
                              bool extractCliques, bool isConflict) {
  mipsolver.mipdata_->debugSolution.checkCut(Rindex, Rvalue, Rlen, rhs);

  sortBuffer.resize(Rlen);

  // compute 1/||a|| for the cut
  // as it is only computed once
  double norm = 0.0;
  double maxabscoef = 0.0;
  for (HighsInt i = 0; i != Rlen; ++i) {
    norm += Rvalue[i] * Rvalue[i];
    maxabscoef = std::max(maxabscoef, std::abs(Rvalue[i]));
    sortBuffer[i].first = Rindex[i];
    sortBuffer[i].second = Rvalue[i];
  }
  pdqsort_branchless(
      sortBuffer.begin(), sortBuffer.end(),
      [](const std::pair<HighsInt, double>& a,
         const std::pair<HighsInt, double>& b) { return a.first < b.first; });
  for (HighsInt i = 0; i != Rlen; ++i) {
    Rindex[i] = sortBuffer[i].first;
    Rvalue[i] = sortBuffer[i].second;
  }
  uint64_t h = compute_cut_hash(Rindex, Rvalue, maxabscoef, Rlen);
  double normalization = 1.0 / double(sqrt(norm));

  // Warning: This global duplicate check assumes the global pool doesn't
  // have cuts added or deleted during time when local pools can add a cut.
  if (this != &mipsolver.mipdata_->getCutPool()) {
    if (mipsolver.mipdata_->getCutPool().isDuplicate(h, normalization, Rindex,
                                                     Rvalue, Rlen, rhs)) {
      return -1;
    }
  }
  if (isDuplicate(h, normalization, Rindex, Rvalue, Rlen, rhs)) return -1;

  // if (Rlen > 0.15 * matrix_.numCols())
  //   printf("cut with len %d not propagated\n", Rlen);
  if (propagate) {
    HighsInt newPropNzs = numPropNzs + Rlen;

    double avgModelNzs = mipsolver.numNonzero() / (double)mipsolver.numRow();

    double newAvgPropNzs = newPropNzs / (double)(numPropRows + 1);

    constexpr double alpha = 2.0;
    if (isConflict) {
      // for conflicts we allow an increased average propagation density
      if (newAvgPropNzs > std::max(alpha * avgModelNzs, minDensityLim)) {
        propagate = false;
      } else {
        ++numPropRows;
        numPropNzs = newPropNzs;
      }
    } else {
      // for cuts we do not want to accept any dense cuts and don't use the
      // average but its actual length
      if (Rlen >= std::max(alpha * avgModelNzs, minDensityLim)) {
        propagate = false;
      } else {
        ++numPropRows;
        numPropNzs = newPropNzs;
      }
    }
  }

  // if we have more than twice the number of nonzeros of the model in use for
  // propagation we stop propagating the rows with the highest age
  HighsInt propRowExcessNzs = numPropNzs - 2 * mipsolver.numNonzero();
  if (propRowExcessNzs > 0) {
    auto it = propRows.rbegin();

    while (propRowExcessNzs > 0 && it != propRows.rend()) {
      HighsInt len = getRowLength(it->second);
      propRowExcessNzs -= len;
      numPropNzs -= len;
      --numPropRows;
      ++it;
    }

    for (auto i = propRows.rbegin(); i != it; ++i) {
      HighsInt row = i->second;
      matrix_.unlinkColumns(row);
      for (HighsDomain::CutpoolPropagation* propagationdomain :
           propagationDomains)
        propagationdomain->cutDeleted(row, true);
    }

    propRows.erase(it.base(), propRows.end());
  }

  // if no such cut exists we append the new cut
  HighsInt rowindex = matrix_.addRow(Rindex, Rvalue, Rlen, propagate);
  hashToCutMap.emplace(h, rowindex);

  if (rowindex == int(rhs_.size())) {
    rhs_.resize(rowindex + 1);
    ages_.resize(rowindex + 1);
    numLps_.resize(rowindex + 1);
    ageResetWhileLocked_.resize(rowindex + 1);
    hasSynced_.resize(rowindex + 1);
    rownormalization_.resize(rowindex + 1);
    maxabscoef_.resize(rowindex + 1);
    rowintegral.resize(rowindex + 1);
  }

  // set the right hand side and reset the age
  rhs_[rowindex] = rhs;
  ages_[rowindex] = std::max(HighsInt{0}, agelim_ - 5);
  ++ageDistribution[ages_[rowindex]];
  rowintegral[rowindex] = integral;
  numLps_[rowindex] = 0;
  ageResetWhileLocked_[rowindex].store(0, std::memory_order_relaxed);
  hasSynced_[rowindex] = false;
  if (propagate) propRows.emplace(ages_[rowindex], rowindex);
  assert((HighsInt)propRows.size() == numPropRows);

  rownormalization_[rowindex] = normalization;
  maxabscoef_[rowindex] = maxabscoef;

  // printf("density: %.2f%%\n", 100.0 * Rlen / (double)matrix_.numCols());
  for (HighsDomain::CutpoolPropagation* propagationdomain : propagationDomains)
    propagationdomain->cutAdded(rowindex, propagate);

  if (extractCliques && this == &mipsolver.mipdata_->getCutPool()) {
    // if this is the global cutpool extract cliques from the cut
    if (Rlen <= 100)
      mipsolver.mipdata_->cliquetable.extractCliquesFromCut(mipsolver, Rindex,
                                                            Rvalue, Rlen, rhs);
  }

  return rowindex;
}

void HighsCutPool::syncCutPool(const HighsMipSolver& mipsolver,
                               HighsCutPool& syncpool) {
  HighsInt cutIndexEnd = matrix_.getNumRows();
  std::vector<HighsInt> idxs;
  std::vector<double> vals;

  for (HighsInt i = 0; i != cutIndexEnd; ++i) {
    // Only sync cuts in the LP that are not already synced
    if ((numLps_[i] > 0 ||
         ageResetWhileLocked_[i].load(std::memory_order_relaxed) == 1) &&
        !hasSynced_[i]) {
      HighsInt Rlen;
      const HighsInt* Rindex;
      const double* Rvalue;
      getCut(i, Rlen, Rindex, Rvalue);
      // copy cut into something mutable (addCut reorders so can't take const)
      idxs.assign(Rindex, Rindex + Rlen);
      vals.assign(Rvalue, Rvalue + Rlen);
      syncpool.addCut(mipsolver, idxs.data(), vals.data(), Rlen, rhs_[i],
                      rowintegral[i]);
      hasSynced_[i] = true;
    }
  }

  assert((HighsInt)propRows.size() == numPropRows);
}