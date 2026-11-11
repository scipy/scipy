/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsDomain.h"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsConflictPool.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsMipSolverData.h"

static HighsCDouble activityContributionMin(double coef, const double& lb,
                                            const double& ub) {
  if (coef < 0) {
    if (ub == kHighsInf) return -kHighsInf;

    return coef * static_cast<HighsCDouble>(ub);
  } else {
    if (lb == -kHighsInf) return -kHighsInf;

    return coef * static_cast<HighsCDouble>(lb);
  }
}

static HighsCDouble activityContributionMax(double coef, const double& lb,
                                            const double& ub) {
  if (coef < 0) {
    if (lb == -kHighsInf) return kHighsInf;

    return coef * static_cast<HighsCDouble>(lb);
  } else {
    if (ub == kHighsInf) return kHighsInf;

    return coef * static_cast<HighsCDouble>(ub);
  }
}

static HighsCDouble computeDelta(double val, double oldbound, double newbound,
                                 double inf, HighsInt& numinfs) {
  if (oldbound == inf) {
    --numinfs;
    return static_cast<HighsCDouble>(newbound) * val;
  } else if (newbound == inf) {
    ++numinfs;
    return static_cast<HighsCDouble>(-oldbound) * val;
  } else {
    return (static_cast<HighsCDouble>(newbound) -
            static_cast<HighsCDouble>(oldbound)) *
           val;
  }
}

static inline double boundRange(double upper_bound, double lower_bound,
                                double tolerance, HighsVarType var_type) {
  // the range is infinite if one of the bounds is infinite
  if (lower_bound <= -kHighsInf || upper_bound >= kHighsInf) return kHighsInf;
  double range = upper_bound - lower_bound;
  return range - (var_type == HighsVarType::kContinuous
                      ? std::max(0.3 * range, 1000.0 * tolerance)
                      : tolerance);
}

HighsDomain::HighsDomain(HighsMipSolver& mipsolver) : mipsolver(&mipsolver) {
  col_lower_ = mipsolver.model_->col_lower_;
  col_upper_ = mipsolver.model_->col_upper_;
  colLowerPos_.assign(mipsolver.numCol(), -1);
  colUpperPos_.assign(mipsolver.numCol(), -1);
  changedcolsflags_.resize(mipsolver.numCol());
  changedcols_.reserve(mipsolver.numCol());
  infeasible_reason = Reason::unspecified();
  infeasible_ = false;
  dualFixProbingPropagation.domain = this;
  dualFixProbingPropagation.mipsolver = &mipsolver;
}

void HighsDomain::addCutpool(HighsCutPool& cutpool) {
  HighsInt cutpoolindex = static_cast<HighsInt>(cutpoolpropagation.size());
  cutpoolpropagation.emplace_back(cutpoolindex, this, cutpool);
}

void HighsDomain::addConflictPool(HighsConflictPool& conflictPool) {
  HighsInt conflictPoolIndex =
      static_cast<HighsInt>(conflictPoolPropagation.size());
  conflictPoolPropagation.emplace_back(conflictPoolIndex, this, conflictPool);
}

void HighsDomain::ConflictPoolPropagation::linkWatchedLiteral(
    HighsInt linkPos) {
  assert(watchedLiterals_[linkPos].domchg.column != -1);
  HighsInt& head =
      watchedLiterals_[linkPos].domchg.boundtype == HighsBoundType::kLower
          ? colLowerWatched_[watchedLiterals_[linkPos].domchg.column]
          : colUpperWatched_[watchedLiterals_[linkPos].domchg.column];

  watchedLiterals_[linkPos].prev = -1;
  watchedLiterals_[linkPos].next = head;
  if (head != -1) watchedLiterals_[head].prev = linkPos;
  head = linkPos;
}

void HighsDomain::ConflictPoolPropagation::unlinkWatchedLiteral(
    HighsInt linkPos) {
  if (watchedLiterals_[linkPos].domchg.column == -1) return;

  HighsInt& head =
      watchedLiterals_[linkPos].domchg.boundtype == HighsBoundType::kLower
          ? colLowerWatched_[watchedLiterals_[linkPos].domchg.column]
          : colUpperWatched_[watchedLiterals_[linkPos].domchg.column];
  watchedLiterals_[linkPos].domchg.column = -1;
  HighsInt prev = watchedLiterals_[linkPos].prev;
  HighsInt next = watchedLiterals_[linkPos].next;
  if (prev != -1)
    watchedLiterals_[prev].next = next;
  else
    head = next;

  if (next != -1) watchedLiterals_[next].prev = prev;
}

HighsDomain::ConflictPoolPropagation::ConflictPoolPropagation(
    HighsInt conflictpoolindex, HighsDomain* domain,
    HighsConflictPool& conflictpool_)
    : conflictpoolindex(conflictpoolindex),
      domain(domain),
      conflictpool_(&conflictpool_) {
  colLowerWatched_.resize(domain->mipsolver->numCol(), -1);
  colUpperWatched_.resize(domain->mipsolver->numCol(), -1);
  conflictpool_.addPropagationDomain(this);
}

HighsDomain::ConflictPoolPropagation::ConflictPoolPropagation(
    const ConflictPoolPropagation& other)
    : conflictpoolindex(other.conflictpoolindex),
      domain(other.domain),
      conflictpool_(other.conflictpool_),
      colLowerWatched_(other.colLowerWatched_),
      colUpperWatched_(other.colUpperWatched_),
      conflictFlag_(other.conflictFlag_),
      propagateConflictInds_(other.propagateConflictInds_),
      watchedLiterals_(other.watchedLiterals_) {
  if (!domain->mipsolver->mipdata_->parallelLockActive() ||
      conflictpool_ != &domain->mipsolver->mipdata_->getConflictPool())
    conflictpool_->addPropagationDomain(this);
}

HighsDomain::ConflictPoolPropagation&
HighsDomain::ConflictPoolPropagation::operator=(
    const ConflictPoolPropagation& other) {
  if (this == &other) return *this;
  if (conflictpool_) conflictpool_->removePropagationDomain(this);
  conflictpoolindex = other.conflictpoolindex;
  domain = other.domain;
  conflictpool_ = other.conflictpool_;
  colLowerWatched_ = other.colLowerWatched_;
  colUpperWatched_ = other.colUpperWatched_;
  conflictFlag_ = other.conflictFlag_;
  propagateConflictInds_ = other.propagateConflictInds_;
  watchedLiterals_ = other.watchedLiterals_;
  assert(!domain->mipsolver->mipdata_->parallelLockActive());
  if (conflictpool_) conflictpool_->addPropagationDomain(this);
  return *this;
}

HighsDomain::ConflictPoolPropagation::~ConflictPoolPropagation() {
  conflictpool_->removePropagationDomain(this);
}

void HighsDomain::ConflictPoolPropagation::conflictDeleted(HighsInt conflict) {
  assert(conflict < static_cast<HighsInt>(conflictFlag_.size()));
  conflictFlag_[conflict] |= 8;
  unlinkWatchedLiteral(2 * conflict);
  unlinkWatchedLiteral(2 * conflict + 1);
}

void HighsDomain::ConflictPoolPropagation::conflictAdded(HighsInt conflict) {
  HighsInt start = conflictpool_->getConflictRanges()[conflict].first;
  HighsInt end = conflictpool_->getConflictRanges()[conflict].second;
  const std::vector<HighsDomainChange>& conflictEntries =
      conflictpool_->getConflictEntryVector();

  if (static_cast<HighsInt>(conflictFlag_.size()) <= conflict) {
    watchedLiterals_.resize(2 * conflict + 2);
    conflictFlag_.resize(conflict + 1);
  }

  HighsInt numWatched = 0;
  for (HighsInt i = start; i != end; ++i) {
    if (domain->isActive(conflictEntries[i])) continue;
    HighsInt watchPos = 2 * conflict + numWatched;
    watchedLiterals_[watchPos].domchg = conflictEntries[i];
    linkWatchedLiteral(watchPos);
    if (++numWatched == 2) break;
  }
  switch (numWatched) {
    case 0: {
      std::array<std::pair<HighsInt, HighsInt>, 2> latestActive;
      HighsInt numActive = 0;
      for (HighsInt i = start; i != end; ++i) {
        HighsInt pos = conflictEntries[i].boundtype == HighsBoundType::kLower
                           ? domain->colLowerPos_[conflictEntries[i].column]
                           : domain->colUpperPos_[conflictEntries[i].column];
        switch (numActive) {
          case 0:
            latestActive[0] = std::make_pair(pos, i);
            numActive = 1;
            break;
          case 1:
            latestActive[1] = std::make_pair(pos, i);
            numActive = 2;
            if (latestActive[0].first < latestActive[1].first)
              std::swap(latestActive[0], latestActive[1]);
            break;
          case 2:
            if (pos > latestActive[1].first) {
              latestActive[1] = std::make_pair(pos, i);
              if (latestActive[0].first < latestActive[1].first)
                std::swap(latestActive[0], latestActive[1]);
            }
        }
      }
      for (HighsInt i = 0; i < numActive; ++i) {
        HighsInt watchPos = 2 * conflict + i;
        watchedLiterals_[watchPos].domchg =
            conflictEntries[latestActive[i].second];
        linkWatchedLiteral(watchPos);
      }
      break;
    }
    case 1: {
      HighsInt latestActive = -1;
      HighsInt latestPos = -1;

      for (HighsInt i = start; i != end; ++i) {
        // the inactive literal is already watched, so only consider the active
        // ones
        if (!domain->isActive(conflictEntries[i])) continue;
        HighsInt pos = conflictEntries[i].boundtype == HighsBoundType::kLower
                           ? domain->colLowerPos_[conflictEntries[i].column]
                           : domain->colUpperPos_[conflictEntries[i].column];
        if (pos > latestPos) {
          latestActive = i;
          latestPos = pos;
        }
      }
      if (latestActive != -1) {
        HighsInt watchPos = 2 * conflict + 1;
        watchedLiterals_[watchPos].domchg = conflictEntries[latestActive];
        linkWatchedLiteral(watchPos);
      }
      break;
    }
    case 2:
      break;
  }

  conflictFlag_[conflict] = numWatched | (conflictFlag_[conflict] & 4);
  markPropagateConflict(conflict);
}

void HighsDomain::ConflictPoolPropagation::markPropagateConflict(
    HighsInt conflict) {
  if (conflictFlag_[conflict] < 2) {
    propagateConflictInds_.push_back(conflict);
    conflictFlag_[conflict] |= 4;
  }
}

void HighsDomain::ConflictPoolPropagation::updateActivityLbChange(
    HighsInt col, double oldbound, double newbound) {
  assert(!domain->infeasible_);

  for (HighsInt i = colLowerWatched_[col]; i != -1;
       i = watchedLiterals_[i].next) {
    HighsInt conflict = i >> 1;

    const HighsDomainChange& domchg = watchedLiterals_[i].domchg;
    uint8_t numInactiveDelta =
        (domchg.boundval > newbound) - (domchg.boundval > oldbound);
    if (numInactiveDelta != 0) {
      conflictFlag_[conflict] += numInactiveDelta;
      markPropagateConflict(conflict);
    }
  }
}

void HighsDomain::ConflictPoolPropagation::updateActivityUbChange(
    HighsInt col, double oldbound, double newbound) {
  assert(!domain->infeasible_);

  for (HighsInt i = colUpperWatched_[col]; i != -1;
       i = watchedLiterals_[i].next) {
    HighsInt conflict = i >> 1;

    const HighsDomainChange& domchg = watchedLiterals_[i].domchg;
    uint8_t numInactiveDelta =
        (domchg.boundval < newbound) - (domchg.boundval < oldbound);
    if (numInactiveDelta != 0) {
      conflictFlag_[conflict] += numInactiveDelta;
      markPropagateConflict(conflict);
    }
  }
}

void HighsDomain::ConflictPoolPropagation::propagateConflict(
    HighsInt conflict) {
  // remove propagate flag, but keep watched and deleted information
  conflictFlag_[conflict] &= (3 | 8);
  // if two inactive literals are watched or conflict has been deleted skip
  if (conflictFlag_[conflict] >= 2) return;

  if (domain->infeasible_) return;

  const std::vector<HighsDomainChange>& entries =
      conflictpool_->getConflictEntryVector();
  HighsInt start = conflictpool_->getConflictRanges()[conflict].first;
  if (start == -1) {
    unlinkWatchedLiteral(2 * conflict);
    unlinkWatchedLiteral(2 * conflict + 1);
    return;
  }
  HighsInt end = conflictpool_->getConflictRanges()[conflict].second;

  WatchedLiteral* watched = watchedLiterals_.data() + 2 * conflict;

  std::array<HighsInt, 2> inactive;
  uint8_t numInactive = 0;
  for (HighsInt i = start; i != end; ++i) {
    if (domain->isActive(entries[i])) continue;

    inactive[numInactive++] = i;
    if (numInactive == 2) break;
  }

  conflictFlag_[conflict] = numInactive;

  switch (numInactive) {
    case 0:
      assert(!domain->infeasible_);
      domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
      domain->infeasible_ = true;
      domain->infeasible_reason = Reason::cut(
          domain->cutpoolpropagation.size() + conflictpoolindex, conflict);
      domain->infeasible_pos = domain->domchgstack_.size();
      conflictpool_->resetAge(conflict);
      // printf("conflict propagation found infeasibility\n");
      break;
    case 1: {
      HighsDomainChange domchg = domain->flip(entries[inactive[0]]);
      if (!domain->isActive(domchg)) {
        domain->changeBound(
            domain->flip(entries[inactive[0]]),
            Reason::cut(domain->cutpoolpropagation.size() + conflictpoolindex,
                        conflict));
        conflictpool_->resetAge(conflict);
      }
      // printf("conflict propagation found bound change\n");
      break;
    }
    case 2: {
      if (watched[0].domchg != entries[inactive[0]]) {
        unlinkWatchedLiteral(2 * conflict);
        watched[0].domchg = entries[inactive[0]];
        linkWatchedLiteral(2 * conflict);
      }

      if (watched[1].domchg != entries[inactive[1]]) {
        unlinkWatchedLiteral(2 * conflict + 1);
        watched[1].domchg = entries[inactive[1]];
        linkWatchedLiteral(2 * conflict + 1);
      }

      return;
    }
  }
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(HighsInt cutpoolindex,
                                                    HighsDomain* domain,
                                                    HighsCutPool& cutpool_)
    : cutpoolindex(cutpoolindex), domain(domain), cutpool(&cutpool_) {
  cutpool->addPropagationDomain(this);
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(
    const CutpoolPropagation& other)
    : cutpoolindex(other.cutpoolindex),
      domain(other.domain),
      cutpool(other.cutpool),
      activitycuts_(other.activitycuts_),
      activitycutsinf_(other.activitycutsinf_),
      propagatecutflags_(other.propagatecutflags_),
      propagatecutinds_(other.propagatecutinds_),
      capacityThreshold_(other.capacityThreshold_) {
  if (!domain->mipsolver->mipdata_->parallelLockActive() ||
      cutpool != &domain->mipsolver->mipdata_->getCutPool())
    cutpool->addPropagationDomain(this);
}

// Warning: When a domain is copy-assigned, e.g., in `resetLocalDomain`,
// with the line `localdom = getDomain()`, then it is going to notify
// all cut / conflict pools that the original was propagating.
// This would be non-deterministic in the order in which the global
// pool gets notified. Currently, such code is only run in serial.

HighsDomain::CutpoolPropagation& HighsDomain::CutpoolPropagation::operator=(
    const CutpoolPropagation& other) {
  if (this == &other) return *this;
  if (cutpool) cutpool->removePropagationDomain(this);
  cutpoolindex = other.cutpoolindex;
  domain = other.domain;
  cutpool = other.cutpool;
  activitycuts_ = other.activitycuts_;
  activitycutsinf_ = other.activitycutsinf_;
  propagatecutflags_ = other.propagatecutflags_;
  propagatecutinds_ = other.propagatecutinds_;
  capacityThreshold_ = other.capacityThreshold_;
  assert(!domain->mipsolver->mipdata_->parallelLockActive());
  if (cutpool) cutpool->addPropagationDomain(this);
  return *this;
}

HighsDomain::CutpoolPropagation::~CutpoolPropagation() {
  cutpool->removePropagationDomain(this);
}

void HighsDomain::CutpoolPropagation::recomputeCapacityThreshold(HighsInt cut) {
  HighsInt start = cutpool->getMatrix().getRowStart(cut);
  HighsInt end = cutpool->getMatrix().getRowEnd(cut);
  const HighsInt* arindex = cutpool->getMatrix().getARindex();
  const double* arvalue = cutpool->getMatrix().getARvalue();
  capacityThreshold_[cut] = -domain->feastol();
  for (HighsInt i = start; i < end; ++i) {
    if (domain->col_upper_[arindex[i]] == domain->col_lower_[arindex[i]])
      continue;

    double threshold =
        std::fabs(arvalue[i]) * boundRange(domain->col_upper_[arindex[i]],
                                           domain->col_lower_[arindex[i]],
                                           domain->feastol(),
                                           domain->variableType(arindex[i]));

    capacityThreshold_[cut] =
        std::max({capacityThreshold_[cut], threshold, domain->feastol()});
  }
}

void HighsDomain::CutpoolPropagation::cutAdded(HighsInt cut, bool propagate) {
  if (!propagate) {
    if (domain != &domain->mipsolver->mipdata_->getDomain()) return;
    HighsInt start = cutpool->getMatrix().getRowStart(cut);
    HighsInt end = cutpool->getMatrix().getRowEnd(cut);
    const HighsInt* arindex = cutpool->getMatrix().getARindex();
    const double* arvalue = cutpool->getMatrix().getARvalue();

    if (static_cast<HighsInt>(activitycuts_.size()) <= cut) {
      activitycuts_.resize(cut + 1);
      activitycutsinf_.resize(cut + 1);
      propagatecutflags_.resize(cut + 1, 2);
      capacityThreshold_.resize(cut + 1);
    }

    propagatecutflags_[cut] &= ~uint8_t{2};
    domain->computeMinActivity(start, end, arindex, arvalue,
                               activitycutsinf_[cut], activitycuts_[cut]);
  } else {
    HighsInt start = cutpool->getMatrix().getRowStart(cut);
    HighsInt end = cutpool->getMatrix().getRowEnd(cut);
    const HighsInt* arindex = cutpool->getMatrix().getARindex();
    const double* arvalue = cutpool->getMatrix().getARvalue();

    if (static_cast<HighsInt>(activitycuts_.size()) <= cut) {
      activitycuts_.resize(cut + 1);
      activitycutsinf_.resize(cut + 1);
      propagatecutflags_.resize(cut + 1, 2);
      capacityThreshold_.resize(cut + 1);
    }

    propagatecutflags_[cut] &= ~uint8_t{2};
    domain->computeMinActivity(start, end, arindex, arvalue,
                               activitycutsinf_[cut], activitycuts_[cut]);

    recomputeCapacityThreshold(cut);
    markPropagateCut(cut);
  }
}

void HighsDomain::CutpoolPropagation::cutDeleted(
    HighsInt cut, bool deletedOnlyForPropagation) {
  if (deletedOnlyForPropagation &&
      domain == &domain->mipsolver->mipdata_->getDomain()) {
    assert(domain->branchPos_.empty());
    return;
  }

  if (cut < static_cast<HighsInt>(propagatecutflags_.size()))
    propagatecutflags_[cut] |= 2;
}

void HighsDomain::CutpoolPropagation::markPropagateCut(HighsInt cut) {
  if (!propagatecutflags_[cut] &&
      (activitycutsinf_[cut] == 1 ||
       (cutpool->getRhs()[cut] - activitycuts_[cut] <=
        capacityThreshold_[cut]))) {
    propagatecutinds_.push_back(cut);
    propagatecutflags_[cut] |= 1;
  }
}

void HighsDomain::CutpoolPropagation::updateActivityLbChange(
    HighsInt col, double oldbound, double newbound, bool threshold,
    bool activity, bool infeasdomain) {
  if (!infeasdomain) assert(!domain->infeasible_);

  if (newbound < oldbound && threshold) {
    cutpool->getMatrix().forEachNegativeColumnEntry(
        col, [&](HighsInt row, double val) {
          domain->updateThresholdLbChange(col, newbound, val,
                                          capacityThreshold_[row]);
          return true;
        });
  }

  if (!activity) return;

  if (!infeasdomain) {
    cutpool->getMatrix().forEachPositiveColumnEntry(
        col, [&](HighsInt row, double val) {
          assert(val > 0);
          HighsCDouble deltamin = computeDelta(
              val, oldbound, newbound, -kHighsInf, activitycutsinf_[row]);
          activitycuts_[row] += deltamin;

          if (deltamin <= 0) {
            domain->updateThresholdLbChange(col, newbound, val,
                                            capacityThreshold_[row]);
            return true;
          }

          if (activitycutsinf_[row] == 0 &&
              activitycuts_[row] - cutpool->getRhs()[row] > domain->feastol()) {
            domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
            domain->infeasible_ = true;
            domain->infeasible_pos = domain->domchgstack_.size();
            domain->infeasible_reason = Reason::cut(cutpoolindex, row);
            return false;
          }

          markPropagateCut(row);

          return true;
        });
  }

  if (domain->infeasible_) {
    // assert(domain->infeasible_reason.type == cutpoolindex);
    assert(domain->infeasible_reason.index >= 0);
    std::swap(oldbound, newbound);
    cutpool->getMatrix().forEachPositiveColumnEntry(
        col, [&](HighsInt row, double val) {
          assert(val > 0);
          activitycuts_[row] += computeDelta(val, oldbound, newbound,
                                             -kHighsInf, activitycutsinf_[row]);

          if (domain->infeasible_reason.index == row &&
              domain->infeasible_reason.type == cutpoolindex)
            return false;

          return true;
        });
  }
}

void HighsDomain::CutpoolPropagation::updateActivityUbChange(
    HighsInt col, double oldbound, double newbound, bool threshold,
    bool activity, bool infeasdomain) {
  if (!infeasdomain) assert(!domain->infeasible_);

  if (newbound > oldbound && threshold) {
    cutpool->getMatrix().forEachPositiveColumnEntry(
        col, [&](HighsInt row, double val) {
          domain->updateThresholdUbChange(col, newbound, val,
                                          capacityThreshold_[row]);
          return true;
        });
  }

  if (!activity) return;

  if (!infeasdomain) {
    cutpool->getMatrix().forEachNegativeColumnEntry(
        col, [&](HighsInt row, double val) {
          assert(val < 0);
          HighsCDouble deltamin = computeDelta(
              val, oldbound, newbound, kHighsInf, activitycutsinf_[row]);

          activitycuts_[row] += deltamin;

          if (deltamin <= 0) {
            domain->updateThresholdUbChange(col, newbound, val,
                                            capacityThreshold_[row]);
            return true;
          }

          if (activitycutsinf_[row] == 0 &&
              activitycuts_[row] - cutpool->getRhs()[row] > domain->feastol()) {
            domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
            domain->infeasible_ = true;
            domain->infeasible_pos = domain->domchgstack_.size();
            domain->infeasible_reason = Reason::cut(cutpoolindex, row);
            return false;
          }

          markPropagateCut(row);

          return true;
        });
  }

  if (domain->infeasible_) {
    // assert(domain->infeasible_reason.type == cutpoolindex);
    assert(domain->infeasible_reason.index >= 0);
    std::swap(oldbound, newbound);
    cutpool->getMatrix().forEachNegativeColumnEntry(
        col, [&](HighsInt row, double val) {
          assert(val < 0);
          activitycuts_[row] += computeDelta(val, oldbound, newbound, kHighsInf,
                                             activitycutsinf_[row]);

          if (domain->infeasible_reason.index == row &&
              domain->infeasible_reason.type == cutpoolindex)
            return false;

          return true;
        });
  }
}

void HighsDomain::DualFixProbingPropagation::recomputeLocks() {
  mipsolver = domain->mipsolver;
  redundantRowFlags_.assign(2 * mipsolver->numRow(), false);
  redundantRowInds_.clear();
  redundantRowInds_.reserve(2 * mipsolver->numRow());
  zeroCostDirections_.assign(mipsolver->numCol(), FixUndecided);
  fixedZeroCostColumns_.clear();
  fixedZeroCostColumns_.reserve(mipsolver->numCol());

  applyingZeroCostFixings_ = false;
  previousRedundantRowSize = 0;
  numGlobalRedundantRows_ = 0;
  numCachedGlobalRedundantRows_ = 0;

  colLowerLocksOriginal_.assign(mipsolver->numCol(), 0);
  colUpperLocksOriginal_.assign(mipsolver->numCol(), 0);
  colLowerReducedNumLocks_.assign(mipsolver->numCol(), 0);
  colUpperReducedNumLocks_.assign(mipsolver->numCol(), 0);

  candidateFixedCols_.clear();
  candidateFixedCols_.reserve(mipsolver->numCol());
  candidateColFixedFlags_.assign(mipsolver->numCol(), false);
  clearColNumReducedLocks_.clear();
  clearColNumReducedLocks_.reserve(mipsolver->numCol());
  globalCandidateFixedCols_.clear();
  globalCandidateFixedCols_.reserve(mipsolver->numCol());

  // compute the locks for each variable
  const HighsLp* model = mipsolver->model_;
  for (HighsInt col = 0; col < model->a_matrix_.num_col_; col++) {
    for (HighsInt k = model->a_matrix_.start_[col];
         k < model->a_matrix_.start_[col + 1]; k++) {
      const HighsInt row = model->a_matrix_.index_[k];
      const double val = model->a_matrix_.value_[k];
      const double lhs = model->row_lower_[row];
      const double rhs = model->row_upper_[row];
      if ((val > 0 && rhs != kHighsInf) || (val < 0 && lhs != -kHighsInf))
        colUpperLocksOriginal_[col]++;
      if ((val > 0 && lhs != -kHighsInf) || (val < 0 && rhs != kHighsInf))
        colLowerLocksOriginal_[col]++;
    }
  }

  // Check for any pre-redundant rows
  for (HighsInt row = 0; row != mipsolver->numRow(); row++) {
    updateRhsRedundant(row);
    updateLhsRedundant(row);
  }
}

void HighsDomain::DualFixProbingPropagation::cacheGlobalRedundantRows() {
  // Go through rows made redundant via domain changes made while
  // dualFixProbingPropagation was not enabled
  for (; numCachedGlobalRedundantRows_ < numGlobalRedundantRows_;
       ++numCachedGlobalRedundantRows_) {
    const RowSide& side = redundantRowInds_[numCachedGlobalRedundantRows_];
    const HighsInt start = mipsolver->mipdata_->ARstart_[side.row];
    const HighsInt end = mipsolver->mipdata_->ARstart_[side.row + 1];
    for (HighsInt i = start; i < end; ++i) {
      const double val = mipsolver->mipdata_->ARvalue_[i];
      if (val == 0.0) continue;
      const HighsInt col = mipsolver->mipdata_->ARindex_[i];
      const bool isUpperLock = (val > 0) == side.isRhs;
      HighsInt& locks = isUpperLock ? colUpperLocksOriginal_[col]
                                    : colLowerLocksOriginal_[col];
      --locks;
      const double cost = mipsolver->model_->col_cost_[col];
      if (locks == 0 &&
          ((isUpperLock && cost <= 0) || (!isUpperLock && cost >= 0))) {
        globalCandidateFixedCols_.push_back(col);
      }
    }
  }

  // Schedule candidate fixes. These will naturally get fixed in the
  // global domain over time.
  HighsInt numCandidates = 0;
  for (HighsInt col : globalCandidateFixedCols_) {
    if (domain->isFixed(col) || candidateColFixedFlags_[col]) continue;
    globalCandidateFixedCols_[numCandidates++] = col;
    candidateColFixedFlags_[col] = true;
    candidateFixedCols_.push_back(col);
  }
  globalCandidateFixedCols_.resize(numCandidates);
}

void HighsDomain::DualFixProbingPropagation::updateRhsRedundant(HighsInt row) {
  RowSide idx{row, true};
  if (domain->activitymaxinf_[row] != 0 || redundantRowFlags_[idx] ||
      mipsolver->model_->row_upper_[row] == kHighsInf)
    return;

  if (domain->isRowUpperRedundant(row)) {
    redundantRowInds_.push_back(idx);
    redundantRowFlags_[idx] = 1;
    if (!isEnabled()) ++numGlobalRedundantRows_;
  }
}

void HighsDomain::DualFixProbingPropagation::updateLhsRedundant(HighsInt row) {
  RowSide idx{row, false};
  if (domain->activitymininf_[row] != 0 || redundantRowFlags_[idx] ||
      mipsolver->model_->row_lower_[row] == -kHighsInf)
    return;

  if (domain->isRowLowerRedundant(row)) {
    redundantRowInds_.push_back(idx);
    redundantRowFlags_[idx] = 1;
    if (!isEnabled()) ++numGlobalRedundantRows_;
  }
}

void HighsDomain::DualFixProbingPropagation::propagate() {
  if (!isActive()) return;

  auto addCandidateFixing = [&](const HighsInt col) {
    if (!candidateColFixedFlags_[col]) {
      candidateFixedCols_.push_back(col);
      candidateColFixedFlags_[col] = true;
    }
  };

  auto collectZeroCostFixing = [&](const HighsInt col,
                                   const DualFixProbingFixDirection direction) {
    fixedZeroCostColumns_.emplace_back(FixedZeroCostColumn{col, direction});
  };

  for (HighsInt i = previousRedundantRowSize;
       i != static_cast<HighsInt>(redundantRowInds_.size()); ++i) {
    const RowSide& idx = redundantRowInds_[i];
    const HighsInt row = idx.row;
    const bool isLhs = !idx.isRhs;
    const HighsInt start = mipsolver->mipdata_->ARstart_[row];
    const HighsInt end = mipsolver->mipdata_->ARstart_[row + 1];
    for (HighsInt j = start; j < end; ++j) {
      const HighsInt col = mipsolver->mipdata_->ARindex_[j];
      if (domain->isFixed(col)) continue;
      const double val = mipsolver->mipdata_->ARvalue_[j];
      const double cost = mipsolver->model_->col_cost_[col];
      if (val != 0.0) {
        // if LHS: Positive val removes a lower lock, negative an upper
        // if RHS: Negative val removes a lower lock, positive an upper
        const bool reducesLowerLock = (val > 0) == isLhs;
        if (reducesLowerLock && cost >= 0) {
          if (colLowerReducedNumLocks_[col] == 0 &&
              colUpperReducedNumLocks_[col] == 0) {
            clearColNumReducedLocks_.push_back(col);
          }
          ++colLowerReducedNumLocks_[col];
          if (colLowerReducedNumLocks_[col] == colLowerLocksOriginal_[col]) {
            addCandidateFixing(col);
          }
        } else if (!reducesLowerLock && cost <= 0) {
          if (colLowerReducedNumLocks_[col] == 0 &&
              colUpperReducedNumLocks_[col] == 0) {
            clearColNumReducedLocks_.push_back(col);
          }
          ++colUpperReducedNumLocks_[col];
          if (colUpperReducedNumLocks_[col] == colUpperLocksOriginal_[col]) {
            addCandidateFixing(col);
          }
        }
      }
    }
  }

  previousRedundantRowSize = static_cast<HighsInt>(redundantRowInds_.size());

  for (const HighsInt col : candidateFixedCols_) {
    if (domain->isFixed(col)) continue;
    const bool canBeFixedToLower =
        ableToFixToLb(col) &&
        colLowerReducedNumLocks_[col] == colLowerLocksOriginal_[col];
    const bool canBeFixedToUpper =
        ableToFixToUb(col) &&
        colUpperReducedNumLocks_[col] == colUpperLocksOriginal_[col];
    if (!canBeFixedToLower && !canBeFixedToUpper) continue;
    const double cost = mipsolver->model_->col_cost_[col];
    if (cost == 0) {
      DualFixProbingFixDirection direction = zeroCostDirections_[col];
      if (direction == FixUndecided) {
        direction = canBeFixedToLower ? FixLowerBound : FixUpperBound;
      }
      if (direction == FixLowerBound && canBeFixedToLower) {
        collectZeroCostFixing(col, FixLowerBound);
      } else if (direction == FixUpperBound && canBeFixedToUpper) {
        collectZeroCostFixing(col, FixUpperBound);
      }
    } else {
      if (canBeFixedToLower) {
        domain->changeBound(HighsBoundType::kUpper, col,
                            domain->col_lower_[col], Reason::unspecified());
      } else if (canBeFixedToUpper) {
        domain->changeBound(HighsBoundType::kLower, col,
                            domain->col_upper_[col], Reason::unspecified());
      }
      if (domain->infeasible()) break;
    }
  }

  for (const auto x : candidateFixedCols_) {
    candidateColFixedFlags_[x] = false;
  }
  candidateFixedCols_.clear();
}

void HighsDomain::DualFixProbingPropagation::propagateZeroCosts() {
  if (fixedZeroCostColumns_.empty()) return;

  if (storeLiftingOpportunity != nullptr) {
    storeLiftingOpportunity();
    storeLiftingOpportunity = nullptr;
  }

  applyingZeroCostFixings_ = true;
  if (zeroCostStartPos_ == kHighsIInf)
    setZeroCostFixingPosition(
        static_cast<HighsInt>(domain->getDomainChangeStack().size()));

  for (const FixedZeroCostColumn& fixing : fixedZeroCostColumns_) {
    if (domain->isFixed(fixing.col)) continue;
    assert(zeroCostDirections_[fixing.col] == FixUndecided ||
           zeroCostDirections_[fixing.col] == fixing.direction);
    if (fixing.direction == FixLowerBound) {
      domain->changeBound(HighsBoundType::kUpper, fixing.col,
                          domain->col_lower_[fixing.col],
                          Reason::unspecified());
    } else {
      domain->changeBound(HighsBoundType::kLower, fixing.col,
                          domain->col_upper_[fixing.col],
                          Reason::unspecified());
    }
    if (!domain->infeasible()) {
      zeroCostDirections_[fixing.col] = fixing.direction;
    } else {
      break;
    }
  }

  fixedZeroCostColumns_.clear();
}

namespace highs {
template <>
struct RbTreeTraits<
    HighsDomain::ObjectivePropagation::ObjectiveContributionTree> {
  using KeyType = std::pair<double, HighsInt>;
  using LinkType = HighsInt;
};
}  // namespace highs

class HighsDomain::ObjectivePropagation::ObjectiveContributionTree
    : public highs::CacheMinRbTree<ObjectiveContributionTree> {
  std::vector<ObjectiveContribution>& nodes;

 public:
  ObjectiveContributionTree(ObjectivePropagation* objProp, HighsInt partition)
      : highs::CacheMinRbTree<ObjectiveContributionTree>(
            objProp->contributionPartitionSets[partition].first,
            objProp->contributionPartitionSets[partition].second),
        nodes(objProp->objectiveLowerContributions) {}

  highs::RbTreeLinks<HighsInt>& getRbTreeLinks(HighsInt node) {
    return nodes[node].links;
  }

  const highs::RbTreeLinks<HighsInt>& getRbTreeLinks(HighsInt node) const {
    return nodes[node].links;
  }

  std::pair<double, HighsInt> getKey(HighsInt node) const {
    return std::make_pair(-nodes[node].contribution, nodes[node].col);
  }
};

HighsDomain::ObjectivePropagation::ObjectivePropagation(HighsDomain* domain)
    : domain(domain),
      objFunc(&domain->mipsolver->mipdata_->objectiveFunction),
      cost(domain->mipsolver->model_->col_cost_.data()) {
  const auto& objNonzeros = objFunc->getObjectiveNonzeros();
  const auto& partitionStarts = objFunc->getCliquePartitionStarts();

  HighsInt numPartitions = objFunc->getNumCliquePartitions();

  if (numPartitions != 0) {
    propagationConsBuffer = objFunc->getObjectiveValuesPacked();
    partitionCliqueData.resize(objFunc->getNumCliquePartitions());
  }

  isPropagated = false;
  capacityThreshold = kHighsInf;
  objectiveLower = 0.0;
  numInfObjLower = 0;
  objectiveLowerContributions.resize(partitionStarts[numPartitions]);
  contributionPartitionSets.resize(numPartitions,
                                   std::make_pair(HighsInt{-1}, HighsInt{-1}));

  // for each clique partition first set up all columns as contributing with
  // their largest possible value and then remove the largest contribution of
  // each clique partition. For these columns we do not need to check for finite
  // bounds as they must be binary. We do not assume, however, that they must be
  // unfixed but could be fixed to 0 or 1. We maintain a heap that fulfills the
  // invariant that the largest column may be not fixed to its value with the
  // larger contribution.
  for (HighsInt i = 0; i < numPartitions; ++i) {
    ObjectiveContributionTree contributionTree(this, i);
    partitionCliqueData[i].rhs = 1;
    for (HighsInt j = partitionStarts[i]; j < partitionStarts[i + 1]; ++j) {
      HighsInt col = objNonzeros[j];
      objectiveLowerContributions[j].col = col;
      objectiveLowerContributions[j].partition = i;
      if (cost[col] > 0.0) {
        objectiveLower += cost[col];
        objectiveLowerContributions[j].contribution = cost[col];
        partitionCliqueData[i].rhs -= 1;
        if (domain->col_lower_[col] == 0.0) contributionTree.link(j);
      } else {
        objectiveLowerContributions[j].contribution = -cost[col];
        if (domain->col_upper_[col] == 1.0) contributionTree.link(j);
      }
    }

    HighsInt worstPos = contributionTree.first();
    if (worstPos != -1)
      objectiveLower -= objectiveLowerContributions[worstPos].contribution;
  }

  // add contribution of remaining objective nonzeros
  const HighsInt numObjNz = static_cast<HighsInt>(objNonzeros.size());
  for (HighsInt i = partitionStarts[numPartitions]; i < numObjNz; ++i) {
    HighsInt col = objNonzeros[i];
    if (cost[col] > 0.0) {
      if (domain->col_lower_[col] == -kHighsInf)
        ++numInfObjLower;
      else
        objectiveLower +=
            static_cast<HighsCDouble>(domain->col_lower_[col]) * cost[col];
    } else {
      if (domain->col_upper_[col] == kHighsInf)
        ++numInfObjLower;
      else
        objectiveLower +=
            static_cast<HighsCDouble>(domain->col_upper_[col]) * cost[col];
    }
  }

  recomputeCapacityThreshold();
  debugCheckObjectiveLower();
}

void HighsDomain::ObjectivePropagation::getPropagationConstraint(
    HighsInt domchgStackSize, const double*& vals, const HighsInt*& inds,
    HighsInt& len, double& rhs, HighsInt domchgCol) {
  const HighsInt numPartitions = objFunc->getNumCliquePartitions();
  inds = objFunc->getObjectiveNonzeros().data();
  len = static_cast<HighsInt>(objFunc->getObjectiveNonzeros().size());
  if (numPartitions == 0) {
    vals = objFunc->getObjectiveValuesPacked().data();
    rhs = domain->mipsolver->mipdata_->upper_limit;
    return;
  }
  const auto& partitionStarts = objFunc->getCliquePartitionStarts();

  HighsCDouble tmpRhs = domain->mipsolver->mipdata_->upper_limit;
  for (HighsInt i = 0; i < numPartitions; ++i) {
    HighsInt start = partitionStarts[i];
    HighsInt end = partitionStarts[i + 1];
    double largest = 0.0;
    for (HighsInt j = start; j < end; ++j) {
      HighsInt c = inds[j];
      HighsInt pos;
      // skip the column we might want to explain a bound change for and take
      // the second largest column instead.
      if (c == domchgCol) continue;
      if (cost[c] > 0) {
        double lb = domain->getColLowerPos(c, domchgStackSize, pos);
        if (lb < 1) largest = std::max(largest, cost[c]);
      } else {
        double ub = domain->getColUpperPos(c, domchgStackSize, pos);
        if (ub > 0) largest = std::max(largest, -cost[c]);
      }
    }

    tmpRhs += static_cast<HighsCDouble>(largest) * partitionCliqueData[i].rhs;
    if (partitionCliqueData[i].multiplier != largest) {
      partitionCliqueData[i].multiplier = largest;
      const auto& packedObjVals = objFunc->getObjectiveValuesPacked();
      for (HighsInt j = start; j < end; ++j)
        propagationConsBuffer[j] =
            packedObjVals[j] - std::copysign(largest, packedObjVals[j]);
    }
  }

  vals = propagationConsBuffer.data();
  rhs = static_cast<double>(tmpRhs);
}

void HighsDomain::ObjectivePropagation::recomputeCapacityThreshold() {
  const auto& partitionStarts = objFunc->getCliquePartitionStarts();
  HighsInt numPartitions = objFunc->getNumCliquePartitions();

  capacityThreshold = -domain->feastol();
  for (HighsInt i = 0; i < numPartitions; ++i) {
    ObjectiveContributionTree contributionTree(this, i);
    HighsInt worstPos = contributionTree.first();
    if (worstPos == -1) continue;
    if (domain->isFixed(objectiveLowerContributions[worstPos].col)) continue;

    double contribution = objectiveLowerContributions[worstPos].contribution;
    HighsInt bestPos = contributionTree.last();
    if (bestPos != worstPos)
      contribution -= objectiveLowerContributions[bestPos].contribution;

    capacityThreshold =
        std::max(capacityThreshold, contribution * (1.0 - domain->feastol()));
  }

  const auto& objNonzeros = objFunc->getObjectiveNonzeros();
  const HighsInt numObjNzs = static_cast<HighsInt>(objNonzeros.size());
  for (HighsInt i = partitionStarts[numPartitions]; i < numObjNzs; ++i) {
    HighsInt col = objNonzeros[i];

    capacityThreshold = std::max(
        capacityThreshold,
        std::fabs(cost[col]) *
            boundRange(domain->col_upper_[col], domain->col_lower_[col],
                       domain->feastol(), domain->variableType(col)));
  }
}

void HighsDomain::ObjectivePropagation::updateActivityLbChange(
    HighsInt col, double oldbound, double newbound) {
  if (cost[col] <= 0.0) {
    if (cost[col] != 0.0 && newbound < oldbound) {
      capacityThreshold =
          std::max(capacityThreshold,
                   -cost[col] * boundRange(domain->col_upper_[col], newbound,
                                           domain->feastol(),
                                           domain->variableType(col)));
      isPropagated = false;
    }
    debugCheckObjectiveLower();
    return;
  }

  isPropagated = false;

  HighsInt partitionPos = objFunc->getColCliquePartition(col);
  if (partitionPos == -1) {
    if (oldbound == -kHighsInf)
      --numInfObjLower;
    else
      objectiveLower -= static_cast<HighsCDouble>(oldbound) * cost[col];

    if (newbound == -kHighsInf)
      ++numInfObjLower;
    else
      objectiveLower += static_cast<HighsCDouble>(newbound) * cost[col];

    debugCheckObjectiveLower();

    if (newbound < oldbound) {
      capacityThreshold = std::max(
          capacityThreshold,
          cost[col] * boundRange(domain->col_upper_[col],
                                 domain->col_lower_[col], domain->feastol(),
                                 domain->variableType(col)));
    } else if (numInfObjLower == 0 &&
               objectiveLower > domain->mipsolver->mipdata_->upper_limit) {
      domain->infeasible_ = true;
      domain->infeasible_pos = domain->domchgstack_.size();
      domain->infeasible_reason = Reason::objective();
      updateActivityLbChange(col, newbound, oldbound);
    }
  } else {
    if (newbound == 0.0) {
      assert(oldbound == 1.0);
      // binary lower bound of variable in clique partition is relaxed to 0
      HighsInt partition = objectiveLowerContributions[partitionPos].partition;
      ObjectiveContributionTree contributionTree(this, partition);
      HighsInt currFirst = contributionTree.first();

      contributionTree.link(partitionPos);

      double oldContribution = 0.0;
      if (currFirst != -1)
        oldContribution = objectiveLowerContributions[currFirst].contribution;

      if (partitionPos == contributionTree.first() &&
          objectiveLowerContributions[partitionPos].contribution !=
              oldContribution) {
        objectiveLower += oldContribution;
        objectiveLower -=
            objectiveLowerContributions[partitionPos].contribution;

        // update the capacity threshold with the difference of the new highest
        // contribution position to the lowest contribution as the column with
        // the lowest contribution can be fixed to its bound that yields the
        // highest objective value.
        HighsInt bestPos = contributionTree.last();
        double delta = objectiveLowerContributions[partitionPos].contribution;
        if (bestPos != partitionPos)
          delta -= objectiveLowerContributions[bestPos].contribution;
        capacityThreshold =
            std::max(delta * (1.0 - domain->feastol()), capacityThreshold);
      } else {
        // the new linked column could be the one with the new lowest
        // contribution so update the capacity threshold to ensure propagation
        // runs when it can be fixed to the bound that yields the highest
        // objective value
        capacityThreshold =
            std::max((oldContribution -
                      objectiveLowerContributions[partitionPos].contribution) *
                         (1.0 - domain->feastol()),
                     capacityThreshold);
      }

      debugCheckObjectiveLower();
    } else {
      // binary lower bound of variable in clique partition is tightened to 1
      assert(oldbound == 0.0);
      assert(newbound == 1.0);

      HighsInt partition = objectiveLowerContributions[partitionPos].partition;
      ObjectiveContributionTree contributionTree(this, partition);
      bool wasFirst = contributionTree.first() == partitionPos;
      if (wasFirst)
        objectiveLower +=
            objectiveLowerContributions[partitionPos].contribution;

      contributionTree.unlink(partitionPos);

      if (wasFirst) {
        HighsInt newWorst = contributionTree.first();
        if (newWorst != -1)
          objectiveLower -= objectiveLowerContributions[newWorst].contribution;
      }

      debugCheckObjectiveLower();

      if (numInfObjLower == 0 &&
          objectiveLower > domain->mipsolver->mipdata_->upper_limit) {
        domain->infeasible_ = true;
        domain->infeasible_pos = domain->domchgstack_.size();
        domain->infeasible_reason = Reason::objective();
        updateActivityLbChange(col, newbound, oldbound);
      }
    }
  }
}

void HighsDomain::ObjectivePropagation::updateActivityUbChange(
    HighsInt col, double oldbound, double newbound) {
  if (cost[col] >= 0.0) {
    if (cost[col] != 0.0 && newbound > oldbound) {
      capacityThreshold = std::max(
          capacityThreshold,
          cost[col] * boundRange(newbound, domain->col_lower_[col],
                                 domain->feastol(), domain->variableType(col)));
      isPropagated = false;
    }
    debugCheckObjectiveLower();
    return;
  }

  isPropagated = false;

  HighsInt partitionPos = objFunc->getColCliquePartition(col);
  if (partitionPos == -1) {
    if (oldbound == kHighsInf)
      --numInfObjLower;
    else
      objectiveLower -= static_cast<HighsCDouble>(oldbound) * cost[col];

    if (newbound == kHighsInf)
      ++numInfObjLower;
    else
      objectiveLower += static_cast<HighsCDouble>(newbound) * cost[col];

    debugCheckObjectiveLower();

    if (newbound > oldbound) {
      capacityThreshold = std::max(
          capacityThreshold,
          -cost[col] * boundRange(domain->col_upper_[col],
                                  domain->col_lower_[col], domain->feastol(),
                                  domain->variableType(col)));
    } else if (numInfObjLower == 0 &&
               objectiveLower > domain->mipsolver->mipdata_->upper_limit) {
      domain->infeasible_ = true;
      domain->infeasible_pos = domain->domchgstack_.size();
      domain->infeasible_reason = Reason::objective();
      updateActivityUbChange(col, newbound, oldbound);
    }
  } else {
    if (newbound == 1.0) {
      assert(oldbound == 0.0);
      // binary upper bound of variable in clique partition is relaxed to 1
      HighsInt partition = objectiveLowerContributions[partitionPos].partition;
      ObjectiveContributionTree contributionTree(this, partition);
      HighsInt currFirst = contributionTree.first();

      contributionTree.link(partitionPos);

      double oldContribution = 0.0;
      if (currFirst != -1)
        oldContribution = objectiveLowerContributions[currFirst].contribution;

      if (partitionPos == contributionTree.first() &&
          objectiveLowerContributions[partitionPos].contribution !=
              oldContribution) {
        objectiveLower += oldContribution;
        objectiveLower -=
            objectiveLowerContributions[partitionPos].contribution;

        // update the capacity threshold with the difference of the new highest
        // contribution position to the lowest contribution as the column with
        // the lowest contribution can be fixed to its bound that yields the
        // highest objective value.
        HighsInt bestPos = contributionTree.last();
        double delta = objectiveLowerContributions[partitionPos].contribution;
        if (bestPos != partitionPos)
          delta -= objectiveLowerContributions[bestPos].contribution;
        capacityThreshold =
            std::max(delta * (1.0 - domain->feastol()), capacityThreshold);
      } else {
        // the new linked column could be the one with the new lowest
        // contribution so update the capacity threshold to ensure propagation
        // runs when it can be fixed to the bound that yields the highest
        // objective value
        capacityThreshold =
            std::max((oldContribution -
                      objectiveLowerContributions[partitionPos].contribution) *
                         (1.0 - domain->feastol()),
                     capacityThreshold);
      }

      debugCheckObjectiveLower();
    } else {
      // binary upper bound of variable in clique partition is tightened to 0
      assert(oldbound == 1.0);
      assert(newbound == 0.0);
      HighsInt partition = objectiveLowerContributions[partitionPos].partition;
      ObjectiveContributionTree contributionTree(this, partition);
      bool wasFirst = contributionTree.first() == partitionPos;
      if (wasFirst)
        objectiveLower +=
            objectiveLowerContributions[partitionPos].contribution;

      contributionTree.unlink(partitionPos);

      if (wasFirst) {
        HighsInt newWorst = contributionTree.first();
        if (newWorst != -1)
          objectiveLower -= objectiveLowerContributions[newWorst].contribution;
      }

      debugCheckObjectiveLower();

      if (numInfObjLower == 0 &&
          objectiveLower > domain->mipsolver->mipdata_->upper_limit) {
        domain->infeasible_ = true;
        domain->infeasible_pos = domain->domchgstack_.size();
        domain->infeasible_reason = Reason::objective();
        updateActivityUbChange(col, newbound, oldbound);
      }
    }
  }
}

bool HighsDomain::ObjectivePropagation::shouldBePropagated() const {
  if (isPropagated) return false;
  if (numInfObjLower > 1) return false;
  if (domain->infeasible_) return false;
  double upperLimit = domain->mipsolver->mipdata_->upper_limit;
  if (upperLimit == kHighsInf) return false;
  if (upperLimit - objectiveLower > capacityThreshold) return false;

  return true;
}

void HighsDomain::ObjectivePropagation::debugCheckObjectiveLower() const {
#ifndef NDEBUG
  if (domain->infeasible_) return;
  HighsCDouble lowerFromScratch = 0.0;
  HighsInt numInf = 0;
  const HighsInt numPartitions = objFunc->getNumCliquePartitions();
  const auto& partitionStarts = objFunc->getCliquePartitionStarts();
  const auto& objNonzeros = objFunc->getObjectiveNonzeros();

  const HighsInt numObjNzs = static_cast<HighsInt>(objNonzeros.size());
  for (HighsInt i = 0; i < numPartitions; ++i) {
    HighsInt start = partitionStarts[i];
    HighsInt end = partitionStarts[i + 1];
    double largest = 0.0;
    for (HighsInt j = start; j < end; ++j) {
      HighsInt c = objNonzeros[j];
      if (cost[c] > 0) {
        lowerFromScratch += cost[c];

        if (domain->col_lower_[c] < 1) largest = std::max(largest, cost[c]);
      } else {
        if (domain->col_upper_[c] > 0) largest = std::max(largest, -cost[c]);
      }
    }
    lowerFromScratch -= largest;
  }

  for (HighsInt i = partitionStarts[numPartitions]; i < numObjNzs; ++i) {
    HighsInt col = objNonzeros[i];
    if (cost[col] > 0) {
      if (domain->col_lower_[col] == 0.0) continue;
      if (domain->col_lower_[col] > -kHighsInf)
        lowerFromScratch +=
            static_cast<HighsCDouble>(domain->col_lower_[col]) * cost[col];
      else
        ++numInf;
    } else {
      if (domain->col_upper_[col] == 0.0) continue;
      if (domain->col_upper_[col] < kHighsInf)
        lowerFromScratch +=
            static_cast<HighsCDouble>(domain->col_upper_[col]) * cost[col];
      else
        ++numInf;
    }
  }
  assert(abs(lowerFromScratch - objectiveLower) <= domain->feastol());
  assert(numInf == numInfObjLower);
#endif
}

void HighsDomain::ObjectivePropagation::propagate() {
  if (!shouldBePropagated()) return;

  debugCheckObjectiveLower();

  const double upperLimit = domain->mipsolver->mipdata_->upper_limit;
  if (numInfObjLower == 0 && objectiveLower > upperLimit) {
    domain->infeasible_ = true;
    domain->infeasible_pos = domain->domchgstack_.size();
    domain->infeasible_reason = Reason::objective();
    return;
  }

  const auto& objNonzeros = objFunc->getObjectiveNonzeros();

  HighsCDouble capacity = upperLimit - objectiveLower;
  if (numInfObjLower == 1) {
    // Scan non-binary columns for infinite bound contribution until the one
    // column that contributes with an infinite bound is found which is the only
    // column that can be propagated
    HighsInt numCol = static_cast<HighsInt>(objNonzeros.size());
    for (HighsInt i = objFunc->getNumBinariesInObjective(); i < numCol; ++i) {
      HighsInt col = objNonzeros[i];
      if ((cost[col] > 0 && domain->col_lower_[col] != -kHighsInf) ||
          (cost[col] < 0 && domain->col_upper_[col] != kHighsInf))
        continue;
      domain->checkChangeBound(
          cost[col] > 0 ? HighsBoundType::kUpper : HighsBoundType::kLower, col,
          capacity / cost[col], Reason::objective());
      break;
    }
  } else {
    HighsInt numPartitions = objFunc->getNumCliquePartitions();
    while (true) {
      HighsInt numBoundChanges = 0;
      for (HighsInt i = 0; i < numPartitions; ++i) {
        ObjectiveContributionTree contributionTree(this, i);
        HighsInt worst = contributionTree.first();
        if (worst == -1) continue;

        double contribution = objectiveLowerContributions[worst].contribution;

        HighsInt secondWorst = contributionTree.successor(worst);
        if (secondWorst != -1)
          contribution -= objectiveLowerContributions[secondWorst].contribution;

        // the upper limit already uses a tolerance, so we can do a hard cutoff
        if (contribution > capacity) {
          HighsInt col = objectiveLowerContributions[worst].col;
          if (cost[col] > 0) {
            if (domain->col_upper_[col] > 0.0) {
              numBoundChanges++;
              domain->changeBound(HighsBoundType::kUpper, col, 0.0,
                                  Reason::objective());
              if (domain->infeasible_) break;
            }
          } else {
            if (domain->col_lower_[col] < 1.0) {
              numBoundChanges++;
              domain->changeBound(HighsBoundType::kLower, col, 1.0,
                                  Reason::objective());
              if (domain->infeasible_) break;
            }
          }
        } else if (secondWorst != -1) {
          // it might be that we can the column with the lowest possible
          // contribution to its bound value that yields the highest objective
          // contribution.
          HighsInt best = contributionTree.last();
          while (best != contributionTree.first()) {
            // difference to the column with the highest contribution is the
            // objective increase when fixing the column to its bound with the
            // lowest objective contribution. Due to the clique information that
            // means the current column with the highest contribution will
            // contribute with its worst bound.
            if (objectiveLowerContributions[contributionTree.first()]
                        .contribution -
                    objectiveLowerContributions[best].contribution >
                capacity) {
              HighsInt col = objectiveLowerContributions[best].col;
              if (cost[col] > 0) {
                assert(domain->col_lower_[col] < 1.0);
                numBoundChanges++;
                domain->changeBound(HighsBoundType::kLower, col, 1.0,
                                    Reason::objective());
                if (domain->infeasible_) break;
              } else {
                assert(domain->col_upper_[col] > 0.0);
                numBoundChanges++;
                domain->changeBound(HighsBoundType::kUpper, col, 0.0,
                                    Reason::objective());
                if (domain->infeasible_) break;
              }
            } else
              break;

            best = contributionTree.last();
          }

          if (domain->infeasible_) break;
        }
      }

      if (domain->infeasible_) break;

      const HighsInt numObjNzs = static_cast<HighsInt>(objNonzeros.size());
      for (HighsInt i = objFunc->getCliquePartitionStarts()[numPartitions];
           i < numObjNzs; ++i) {
        HighsInt col = objNonzeros[i];

        if (cost[col] > 0) {
          if (domain->checkChangeBound(
                  HighsBoundType::kUpper, col,
                  (capacity +
                   static_cast<HighsCDouble>(domain->col_lower_[col]) *
                       cost[col]) /
                      cost[col],
                  Reason::objective()))
            numBoundChanges++;
        } else {
          if (domain->checkChangeBound(
                  HighsBoundType::kLower, col,
                  (capacity +
                   static_cast<HighsCDouble>(domain->col_upper_[col]) *
                       cost[col]) /
                      cost[col],
                  Reason::objective()))
            numBoundChanges++;
        }
        if (domain->infeasible_) break;
      }
      if (domain->infeasible_) break;

      if (numBoundChanges == 0) break;
      capacity = upperLimit - objectiveLower;
    }
  }

  recomputeCapacityThreshold();
  isPropagated = true;
}

void HighsDomain::computeMinActivity(HighsInt start, HighsInt end,
                                     const HighsInt* ARindex,
                                     const double* ARvalue, HighsInt& ninfmin,
                                     HighsCDouble& activitymin) const {
  if (infeasible_) {
    activitymin = 0.0;
    ninfmin = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < static_cast<HighsInt>(col_lower_.size()));

      HighsInt tmp;
      double lb = getColLowerPos(col, infeasible_pos - 1, tmp);
      double ub = getColUpperPos(col, infeasible_pos - 1, tmp);
      HighsCDouble contributionmin = activityContributionMin(val, lb, ub);

      if (contributionmin == -kHighsInf)
        ++ninfmin;
      else
        activitymin += contributionmin;
    }

  } else {
    activitymin = 0.0;
    ninfmin = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < static_cast<HighsInt>(col_lower_.size()));

      HighsCDouble contributionmin =
          activityContributionMin(val, col_lower_[col], col_upper_[col]);

      if (contributionmin == -kHighsInf)
        ++ninfmin;
      else
        activitymin += contributionmin;
    }
  }
  activitymin.renormalize();
}

void HighsDomain::computeMaxActivity(HighsInt start, HighsInt end,
                                     const HighsInt* ARindex,
                                     const double* ARvalue, HighsInt& ninfmax,
                                     HighsCDouble& activitymax) const {
  if (infeasible_) {
    activitymax = 0.0;
    ninfmax = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < static_cast<HighsInt>(col_lower_.size()));

      HighsInt tmp;
      double lb = getColLowerPos(col, infeasible_pos - 1, tmp);
      double ub = getColUpperPos(col, infeasible_pos - 1, tmp);
      HighsCDouble contributionmin = activityContributionMax(val, lb, ub);

      if (contributionmin == kHighsInf)
        ++ninfmax;
      else
        activitymax += contributionmin;
    }
  } else {
    activitymax = 0.0;
    ninfmax = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < static_cast<HighsInt>(col_lower_.size()));

      HighsCDouble contributionmin =
          activityContributionMax(val, col_lower_[col], col_upper_[col]);

      if (contributionmin == kHighsInf)
        ++ninfmax;
      else
        activitymax += contributionmin;
    }
  }
  activitymax.renormalize();
}

double HighsDomain::adjustedUb(HighsInt col, HighsCDouble boundVal,
                               bool& accept) const {
  double bound;

  if (mipsolver->isColIntegral(col)) {
    bound = static_cast<double>(floor(boundVal + feastol()));
    accept = bound < col_upper_[col] &&
             col_upper_[col] - bound > 1000.0 * feastol() * std::fabs(bound);
  } else {
    if (abs(boundVal - col_lower_[col]) <= epsilon())
      bound = col_lower_[col];
    else
      bound = static_cast<double>(boundVal);
    if (col_upper_[col] == kHighsInf)
      accept = true;
    else if (bound + 1000.0 * feastol() < col_upper_[col]) {
      double relativeImprove = col_upper_[col] - bound;
      if (col_lower_[col] != -kHighsInf)
        relativeImprove /= col_upper_[col] - col_lower_[col];
      else
        relativeImprove /=
            std::max(std::fabs(col_upper_[col]), std::fabs(bound));
      accept = relativeImprove >= 0.3;
    } else
      accept = false;
  }

  return bound;
}

double HighsDomain::adjustedLb(HighsInt col, HighsCDouble boundVal,
                               bool& accept) const {
  double bound;

  if (mipsolver->isColIntegral(col)) {
    bound = static_cast<double>(ceil(boundVal - feastol()));
    accept = bound > col_lower_[col] &&
             bound - col_lower_[col] > 1000.0 * feastol() * std::fabs(bound);
  } else {
    if (abs(col_upper_[col] - boundVal) <= epsilon())
      bound = col_upper_[col];
    else
      bound = static_cast<double>(boundVal);
    if (col_lower_[col] == -kHighsInf)
      accept = true;
    else if (bound - 1000.0 * feastol() > col_lower_[col]) {
      double relativeImprove = bound - col_lower_[col];
      if (col_upper_[col] != kHighsInf)
        relativeImprove /= col_upper_[col] - col_lower_[col];
      else
        relativeImprove /=
            std::max(std::fabs(col_lower_[col]), std::fabs(bound));
      accept = relativeImprove >= 0.3;
    } else
      accept = false;
  }

  return bound;
}

HighsInt HighsDomain::propagateRowUpper(const HighsInt* Rindex,
                                        const double* Rvalue, HighsInt Rlen,
                                        double Rupper,
                                        const HighsCDouble& minactivity,
                                        HighsInt ninfmin,
                                        HighsDomainChange* boundchgs) const {
  assert(std::isfinite(static_cast<double>(minactivity)));
  if (ninfmin > 1) return 0;
  HighsInt numchgs = 0;
  for (HighsInt i = 0; i != Rlen; ++i) {
    HighsCDouble minresact;
    HighsCDouble actcontribution = activityContributionMin(
        Rvalue[i], col_lower_[Rindex[i]], col_upper_[Rindex[i]]);
    if (ninfmin == 1) {
      if (actcontribution != -kHighsInf) continue;

      minresact = minactivity;
    } else {
      minresact = minactivity - actcontribution;
    }

    HighsCDouble boundVal = (Rupper - minresact) / Rvalue[i];
    if (abs(boundVal) * kHighsTiny > feastol()) continue;

    if (Rvalue[i] > 0) {
      bool accept;

      double bound = adjustedUb(Rindex[i], boundVal, accept);
      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kUpper};

    } else {
      bool accept;

      double bound = adjustedLb(Rindex[i], boundVal, accept);
      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kLower};
    }
  }

  return numchgs;
}

HighsInt HighsDomain::propagateRowLower(const HighsInt* Rindex,
                                        const double* Rvalue, HighsInt Rlen,
                                        double Rlower,
                                        const HighsCDouble& maxactivity,
                                        HighsInt ninfmax,
                                        HighsDomainChange* boundchgs) const {
  assert(std::isfinite(static_cast<double>(maxactivity)));
  if (ninfmax > 1) return 0;
  HighsInt numchgs = 0;
  for (HighsInt i = 0; i != Rlen; ++i) {
    HighsCDouble maxresact;
    HighsCDouble actcontribution = activityContributionMax(
        Rvalue[i], col_lower_[Rindex[i]], col_upper_[Rindex[i]]);
    if (ninfmax == 1) {
      if (actcontribution != kHighsInf) continue;

      maxresact = maxactivity;
    } else {
      maxresact = maxactivity - actcontribution;
    }

    HighsCDouble boundVal = (Rlower - maxresact) / Rvalue[i];
    if (abs(boundVal) * kHighsTiny > feastol()) continue;

    if (Rvalue[i] < 0) {
      bool accept;

      double bound = adjustedUb(Rindex[i], boundVal, accept);
      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kUpper};
    } else {
      bool accept;

      double bound = adjustedLb(Rindex[i], boundVal, accept);
      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kLower};
    }
  }

  return numchgs;
}

void HighsDomain::updateThresholdLbChange(HighsInt col, double newbound,
                                          double val, double& threshold) const {
  if (newbound != col_upper_[col]) {
    double thresholdNew =
        std::fabs(val) *
        boundRange(col_upper_[col], newbound, feastol(), variableType(col));

    // the new threshold is now the maximum of the new threshold and the current
    // one
    threshold = std::max({threshold, thresholdNew, feastol()});
  }
}

void HighsDomain::updateThresholdUbChange(HighsInt col, double newbound,
                                          double val, double& threshold) const {
  if (newbound != col_lower_[col]) {
    double thresholdNew =
        std::fabs(val) *
        boundRange(newbound, col_lower_[col], feastol(), variableType(col));

    // the new threshold is now the maximum of the new threshold and the current
    // one
    threshold = std::max({threshold, thresholdNew, feastol()});
  }
}

void HighsDomain::updateActivityLbChange(HighsInt col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  HighsInt start = mip->a_matrix_.start_[col];
  HighsInt end = mip->a_matrix_.start_[col + 1];

  assert(!infeasible_);

  if (objProp_.isActive()) {
    objProp_.updateActivityLbChange(col, oldbound, newbound);
    if (infeasible_) return;
  }

  const bool trackRedundancy = newbound > oldbound && getDualFixProbingActive();

  for (HighsInt i = start; i != end; ++i) {
    if (mip->a_matrix_.value_[i] > 0) {
      HighsCDouble deltamin =
          computeDelta(mip->a_matrix_.value_[i], oldbound, newbound, -kHighsInf,
                       activitymininf_[mip->a_matrix_.index_[i]]);
      activitymin_[mip->a_matrix_.index_[i]] += deltamin;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpminact;
        computeMinActivity(
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i]],
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i] + 1],
            mipsolver->mipdata_->ARindex_.data(),
            mipsolver->mipdata_->ARvalue_.data(), tmpinf, tmpminact);
        assert(abs(activitymin_[mip->a_matrix_.index_[i]] - tmpminact) <=
               feastol());
        assert(tmpinf == activitymininf_[mip->a_matrix_.index_[i]]);
      }
#endif
      if (recordRedundantRows_ &&
          !dualFixProbingPropagation.isZeroCostFixingActive() &&
          mip->row_lower_[mip->a_matrix_.index_[i]] != -kHighsInf &&
          mip->row_upper_[mip->a_matrix_.index_[i]] == kHighsInf)
        updateRedundantRows(mip->a_matrix_.index_[i]);

      if (trackRedundancy)
        dualFixProbingPropagation.updateLhsRedundant(mip->a_matrix_.index_[i]);

      if (deltamin <= 0) {
        updateThresholdLbChange(col, newbound, mip->a_matrix_.value_[i],
                                capacityThreshold_[mip->a_matrix_.index_[i]]);
        continue;
      }

      if (isRowUpperInfeasible(mip->a_matrix_.index_[i])) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowUpper(mip->a_matrix_.index_[i]);
        end = i + 1;
        break;
      }

      if (activitymininf_[mip->a_matrix_.index_[i]] <= 1 &&
          !propagateflags_[mip->a_matrix_.index_[i]] &&
          mip->row_upper_[mip->a_matrix_.index_[i]] != kHighsInf)
        markPropagate(mip->a_matrix_.index_[i]);
    } else {
      HighsCDouble deltamax =
          computeDelta(mip->a_matrix_.value_[i], oldbound, newbound, -kHighsInf,
                       activitymaxinf_[mip->a_matrix_.index_[i]]);
      activitymax_[mip->a_matrix_.index_[i]] += deltamax;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpmaxact;
        computeMaxActivity(
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i]],
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i] + 1],
            mipsolver->mipdata_->ARindex_.data(),
            mipsolver->mipdata_->ARvalue_.data(), tmpinf, tmpmaxact);
        assert(abs(activitymax_[mip->a_matrix_.index_[i]] - tmpmaxact) <=
               feastol());
        assert(tmpinf == activitymaxinf_[mip->a_matrix_.index_[i]]);
      }
#endif
      if (recordRedundantRows_ &&
          !dualFixProbingPropagation.isZeroCostFixingActive() &&
          mip->row_lower_[mip->a_matrix_.index_[i]] == -kHighsInf &&
          mip->row_upper_[mip->a_matrix_.index_[i]] != kHighsInf)
        updateRedundantRows(mip->a_matrix_.index_[i]);

      if (trackRedundancy)
        dualFixProbingPropagation.updateRhsRedundant(mip->a_matrix_.index_[i]);

      if (deltamax >= 0) {
        updateThresholdLbChange(col, newbound, mip->a_matrix_.value_[i],
                                capacityThreshold_[mip->a_matrix_.index_[i]]);
        continue;
      }

      if (isRowLowerInfeasible(mip->a_matrix_.index_[i])) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowLower(mip->a_matrix_.index_[i]);
        end = i + 1;
        break;
      }

      if (activitymaxinf_[mip->a_matrix_.index_[i]] <= 1 &&
          !propagateflags_[mip->a_matrix_.index_[i]] &&
          mip->row_lower_[mip->a_matrix_.index_[i]] != -kHighsInf)
        markPropagate(mip->a_matrix_.index_[i]);
    }
  }

  if (infeasible_) {
    assert(infeasible_reason.type == Reason::kModelRowLower ||
           infeasible_reason.type == Reason::kModelRowUpper);
    assert(infeasible_reason.index == mip->a_matrix_.index_[end - 1]);
  } else {
    // Explanation: If cutpoolpropagation[i] returns infeasible,
    // we still need to update the domain threshold values using
    // cutpoolpropagation[j], j > i, and reverse the activity
    // changes made by cutpoolpropagation[j] j < i.
    HighsInt infeascutpool = -1;
    HighsInt ncutpoolprop = static_cast<HighsInt>(cutpoolpropagation.size());
    for (HighsInt i = 0; i != ncutpoolprop; ++i) {
      if (!infeasible_) {
        cutpoolpropagation[i].updateActivityLbChange(col, oldbound, newbound,
                                                     true, true, infeasible_);
        if (infeasible_) infeascutpool = i;
      } else {
        cutpoolpropagation[i].updateActivityLbChange(col, oldbound, newbound,
                                                     true, false, infeasible_);
      }
    }
    for (HighsInt i = 0; i < infeascutpool; ++i) {
      cutpoolpropagation[i].updateActivityLbChange(col, oldbound, newbound,
                                                   false, true, infeasible_);
    }
  }

  if (infeasible_) {
    std::swap(oldbound, newbound);
    for (HighsInt i = start; i != end; ++i) {
      if (mip->a_matrix_.value_[i] > 0) {
        activitymin_[mip->a_matrix_.index_[i]] +=
            computeDelta(mip->a_matrix_.value_[i], oldbound, newbound,
                         -kHighsInf, activitymininf_[mip->a_matrix_.index_[i]]);
      } else {
        activitymax_[mip->a_matrix_.index_[i]] +=
            computeDelta(mip->a_matrix_.value_[i], oldbound, newbound,
                         -kHighsInf, activitymaxinf_[mip->a_matrix_.index_[i]]);
      }
    }

    if (objProp_.isActive()) {
      objProp_.updateActivityLbChange(col, oldbound, newbound);
    }

    return;
  } else {
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.updateActivityLbChange(col, oldbound, newbound);
  }
}

void HighsDomain::updateActivityUbChange(HighsInt col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  HighsInt start = mip->a_matrix_.start_[col];
  HighsInt end = mip->a_matrix_.start_[col + 1];

  assert(!infeasible_);

  if (objProp_.isActive()) {
    objProp_.updateActivityUbChange(col, oldbound, newbound);
    if (infeasible_) return;
  }

  const bool trackRedundancy = newbound < oldbound && getDualFixProbingActive();

  for (HighsInt i = start; i != end; ++i) {
    if (mip->a_matrix_.value_[i] > 0) {
      HighsCDouble deltamax =
          computeDelta(mip->a_matrix_.value_[i], oldbound, newbound, kHighsInf,
                       activitymaxinf_[mip->a_matrix_.index_[i]]);
      activitymax_[mip->a_matrix_.index_[i]] += deltamax;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpmaxact;
        computeMaxActivity(
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i]],
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i] + 1],
            mipsolver->mipdata_->ARindex_.data(),
            mipsolver->mipdata_->ARvalue_.data(), tmpinf, tmpmaxact);
        assert(abs(activitymax_[mip->a_matrix_.index_[i]] - tmpmaxact) <=
               feastol());
        assert(tmpinf == activitymaxinf_[mip->a_matrix_.index_[i]]);
      }
#endif
      if (recordRedundantRows_ &&
          !dualFixProbingPropagation.isZeroCostFixingActive() &&
          mip->row_lower_[mip->a_matrix_.index_[i]] == -kHighsInf &&
          mip->row_upper_[mip->a_matrix_.index_[i]] != kHighsInf)
        updateRedundantRows(mip->a_matrix_.index_[i]);

      if (trackRedundancy)
        dualFixProbingPropagation.updateRhsRedundant(mip->a_matrix_.index_[i]);

      if (deltamax >= 0) {
        updateThresholdUbChange(col, newbound, mip->a_matrix_.value_[i],
                                capacityThreshold_[mip->a_matrix_.index_[i]]);
        continue;
      }

      if (isRowLowerInfeasible(mip->a_matrix_.index_[i])) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowLower(mip->a_matrix_.index_[i]);
        end = i + 1;
        break;
      }

      if (activitymaxinf_[mip->a_matrix_.index_[i]] <= 1 &&
          !propagateflags_[mip->a_matrix_.index_[i]] &&
          mip->row_lower_[mip->a_matrix_.index_[i]] != -kHighsInf) {
        markPropagate(mip->a_matrix_.index_[i]);
        // propagateflags_[mip->a_matrix_.index_[i]] = 1;
        // propagateinds_.push_back(mip->a_matrix_.index_[i]);
      }
    } else {
      HighsCDouble deltamin =
          computeDelta(mip->a_matrix_.value_[i], oldbound, newbound, kHighsInf,
                       activitymininf_[mip->a_matrix_.index_[i]]);
      activitymin_[mip->a_matrix_.index_[i]] += deltamin;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpminact;
        computeMinActivity(
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i]],
            mipsolver->mipdata_->ARstart_[mip->a_matrix_.index_[i] + 1],
            mipsolver->mipdata_->ARindex_.data(),
            mipsolver->mipdata_->ARvalue_.data(), tmpinf, tmpminact);
        assert(abs(activitymin_[mip->a_matrix_.index_[i]] - tmpminact) <=
               feastol());
        assert(tmpinf == activitymininf_[mip->a_matrix_.index_[i]]);
      }
#endif
      if (recordRedundantRows_ &&
          !dualFixProbingPropagation.isZeroCostFixingActive() &&
          mip->row_lower_[mip->a_matrix_.index_[i]] != -kHighsInf &&
          mip->row_upper_[mip->a_matrix_.index_[i]] == kHighsInf)
        updateRedundantRows(mip->a_matrix_.index_[i]);

      if (trackRedundancy)
        dualFixProbingPropagation.updateLhsRedundant(mip->a_matrix_.index_[i]);

      if (deltamin <= 0) {
        updateThresholdUbChange(col, newbound, mip->a_matrix_.value_[i],
                                capacityThreshold_[mip->a_matrix_.index_[i]]);
        continue;
      }

      if (isRowUpperInfeasible(mip->a_matrix_.index_[i])) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowUpper(mip->a_matrix_.index_[i]);
        end = i + 1;
        break;
      }

      if (activitymininf_[mip->a_matrix_.index_[i]] <= 1 &&
          !propagateflags_[mip->a_matrix_.index_[i]] &&
          mip->row_upper_[mip->a_matrix_.index_[i]] != kHighsInf) {
        markPropagate(mip->a_matrix_.index_[i]);
        // propagateflags_[mip->a_matrix_.index_[i]] = 1;
        // propagateinds_.push_back(mip->a_matrix_.index_[i]);
      }
    }
  }

  if (infeasible_) {
    assert(infeasible_reason.type == Reason::kModelRowLower ||
           infeasible_reason.type == Reason::kModelRowUpper);
    assert(infeasible_reason.index == mip->a_matrix_.index_[end - 1]);
  } else {
    // Explanation: If cutpoolpropagation[i] returns infeasible,
    // we still need to update the domain threshold values using
    // cutpoolpropagation[j], j > i, and reverse the activity
    // changes made by cutpoolpropagation[j] j < i.
    HighsInt infeascutpool = -1;
    HighsInt ncutpoolprop = static_cast<HighsInt>(cutpoolpropagation.size());
    for (HighsInt i = 0; i != ncutpoolprop; ++i) {
      if (!infeasible_) {
        cutpoolpropagation[i].updateActivityUbChange(col, oldbound, newbound,
                                                     true, true, infeasible_);
        if (infeasible_) infeascutpool = i;
      } else {
        cutpoolpropagation[i].updateActivityUbChange(col, oldbound, newbound,
                                                     true, false, infeasible_);
      }
    }
    for (HighsInt i = 0; i < infeascutpool; ++i) {
      cutpoolpropagation[i].updateActivityUbChange(col, oldbound, newbound,
                                                   false, true, infeasible_);
    }
  }

  if (infeasible_) {
    std::swap(oldbound, newbound);
    for (HighsInt i = start; i != end; ++i) {
      if (mip->a_matrix_.value_[i] > 0) {
        activitymax_[mip->a_matrix_.index_[i]] +=
            computeDelta(mip->a_matrix_.value_[i], oldbound, newbound,
                         kHighsInf, activitymaxinf_[mip->a_matrix_.index_[i]]);
      } else {
        activitymin_[mip->a_matrix_.index_[i]] +=
            computeDelta(mip->a_matrix_.value_[i], oldbound, newbound,
                         kHighsInf, activitymininf_[mip->a_matrix_.index_[i]]);
      }
    }

    if (objProp_.isActive()) {
      objProp_.updateActivityUbChange(col, oldbound, newbound);
    }

    return;
  } else {
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.updateActivityUbChange(col, oldbound, newbound);
  }
}

void HighsDomain::recomputeCapacityThreshold(HighsInt row) {
  HighsInt start = mipsolver->mipdata_->ARstart_[row];
  HighsInt end = mipsolver->mipdata_->ARstart_[row + 1];

  capacityThreshold_[row] = -feastol();
  for (HighsInt i = start; i < end; ++i) {
    HighsInt col = mipsolver->mipdata_->ARindex_[i];

    if (col_upper_[col] == col_lower_[col]) continue;

    double threshold = std::fabs(mipsolver->mipdata_->ARvalue_[i]) *
                       boundRange(col_upper_[col], col_lower_[col], feastol(),
                                  variableType(col));

    capacityThreshold_[row] =
        std::max({capacityThreshold_[row], threshold, feastol()});
  }
}

void HighsDomain::updateRedundantRows(HighsInt row) {
  if (!isRedundantRow(row)) {
    // row that was found to be redundant should not be non-redundant
    assert(redundantRows_.find(row) == nullptr);
    return;
  }
  // row is redundant
  redundantRows_.insert(row);
}

double HighsDomain::getRedundantRowValue(HighsInt row) const {
  if (mipsolver->rowLower(row) != -kHighsInf) {
    assert(mipsolver->rowUpper(row) == kHighsInf);
    return static_cast<double>(activitymin_[row] - mipsolver->rowLower(row));
  } else {
    assert(mipsolver->rowUpper(row) != kHighsInf);
    return static_cast<double>(activitymax_[row] - mipsolver->rowUpper(row));
  }
}

bool HighsDomain::isRedundantRow(HighsInt row) const {
  return isRowLowerRedundant(row) && isRowUpperRedundant(row);
}

void HighsDomain::markPropagateCut(Reason reason) {
  switch (reason.type) {
    case Reason::kUnknown:
    case Reason::kCliqueTable:
    case Reason::kBranching:
    case Reason::kModelRowLower:
    case Reason::kModelRowUpper:
    case Reason::kConflictingBounds:
    case Reason::kObjective:
      break;
    default:
      assert(reason.type >= 0 &&
             reason.type <
                 static_cast<HighsInt>(cutpoolpropagation.size() +
                                       conflictPoolPropagation.size()));
      if (reason.type < static_cast<HighsInt>(cutpoolpropagation.size()))
        cutpoolpropagation[reason.type].markPropagateCut(reason.index);
      else
        conflictPoolPropagation[reason.type - cutpoolpropagation.size()]
            .markPropagateConflict(reason.index);
  }
}

void HighsDomain::markPropagate(HighsInt row) {
  if (!propagateflags_[row]) {
    bool proplower = !isRowLowerRedundant(row) &&
                     (activitymaxinf_[row] == 1 ||
                      activitymax_[row] - mipsolver->rowLower(row) <=
                          capacityThreshold_[row]);
    bool propupper = !isRowUpperRedundant(row) &&
                     (activitymininf_[row] == 1 ||
                      mipsolver->rowUpper(row) - activitymin_[row] <=
                          capacityThreshold_[row]);

    if (proplower || propupper) {
      propagateinds_.push_back(row);
      propagateflags_[row] = true;
    }
  }
}

void HighsDomain::computeRowActivities() {
  activitymin_.resize(mipsolver->numRow());
  activitymininf_.resize(mipsolver->numRow());
  activitymax_.resize(mipsolver->numRow());
  activitymaxinf_.resize(mipsolver->numRow());
  capacityThreshold_.resize(mipsolver->numRow());
  propagateflags_.resize(mipsolver->numRow());
  propagateinds_.reserve(mipsolver->numRow());

  for (HighsInt i = 0; i != mipsolver->numRow(); ++i) {
    HighsInt start = mipsolver->mipdata_->ARstart_[i];
    HighsInt end = mipsolver->mipdata_->ARstart_[i + 1];

    computeMinActivity(start, end, mipsolver->mipdata_->ARindex_.data(),
                       mipsolver->mipdata_->ARvalue_.data(), activitymininf_[i],
                       activitymin_[i]);
    computeMaxActivity(start, end, mipsolver->mipdata_->ARindex_.data(),
                       mipsolver->mipdata_->ARvalue_.data(), activitymaxinf_[i],
                       activitymax_[i]);

    recomputeCapacityThreshold(i);

    if ((activitymininf_[i] <= 1 && mipsolver->rowUpper(i) != kHighsInf) ||
        (activitymaxinf_[i] <= 1 && mipsolver->rowLower(i) != -kHighsInf)) {
      markPropagate(i);
      // propagateflags_[i] = 1;
      // propagateinds_.push_back(i);
    }
  }
}

double HighsDomain::doChangeBound(const HighsDomainChange& boundchg) {
  double oldbound;

  if (boundchg.boundtype == HighsBoundType::kLower) {
    oldbound = col_lower_[boundchg.column];
    col_lower_[boundchg.column] = boundchg.boundval;
    if (oldbound != boundchg.boundval) {
      if (!infeasible_)
        updateActivityLbChange(boundchg.column, oldbound, boundchg.boundval);

      if (!isChangedCol(boundchg.column)) {
        changedcolsflags_[boundchg.column] = true;
        changedcols_.push_back(boundchg.column);
      }
    }
  } else {
    oldbound = col_upper_[boundchg.column];
    col_upper_[boundchg.column] = boundchg.boundval;
    if (oldbound != boundchg.boundval) {
      if (!infeasible_)
        updateActivityUbChange(boundchg.column, oldbound, boundchg.boundval);

      if (!isChangedCol(boundchg.column)) {
        changedcolsflags_[boundchg.column] = true;
        changedcols_.push_back(boundchg.column);
      }
    }
  }

  return oldbound;
}

void HighsDomain::changeBound(HighsDomainChange boundchg, Reason reason) {
  assert(boundchg.column >= 0);
  assert(boundchg.column < static_cast<HighsInt>(col_upper_.size()));
  // assert(infeasible_ == 0);
  // if (reason.type == Reason::kObjective) {
  //   if (!mipsolver->submip)
  //     printf("objective propagator changed %s bound of column %d to %g\n",
  //            boundchg.boundtype == HighsBoundType::kLower ? "lower" :
  //            "upper", boundchg.column, boundchg.boundval);
  // }

  HighsInt prevPos;
  if (boundchg.boundtype == HighsBoundType::kLower) {
    if (mipsolver->isColIntegral(boundchg.column)) {
      boundchg.boundval = std::ceil(boundchg.boundval - feastol());
    }
    if (boundchg.boundval <= col_lower_[boundchg.column]) {
      if (reason.type != Reason::kBranching) return;
      boundchg.boundval = col_lower_[boundchg.column];
    }
    if (boundchg.boundval > col_upper_[boundchg.column]) {
      if (boundchg.boundval - col_upper_[boundchg.column] > feastol()) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        if (!infeasible_) {
          infeasible_pos = domchgstack_.size();
          infeasible_ = true;
          infeasible_reason = Reason::conflictingBounds(domchgstack_.size());
        }
      } else {
        boundchg.boundval = col_upper_[boundchg.column];
        if (boundchg.boundval == col_lower_[boundchg.column]) return;
      }
    }

    prevPos = colLowerPos_[boundchg.column];
    colLowerPos_[boundchg.column] = domchgstack_.size();
  } else {
    if (mipsolver->isColIntegral(boundchg.column)) {
      boundchg.boundval = std::floor(boundchg.boundval + feastol());
    }
    if (boundchg.boundval >= col_upper_[boundchg.column]) {
      if (reason.type != Reason::kBranching) return;
      boundchg.boundval = col_upper_[boundchg.column];
    }
    if (boundchg.boundval < col_lower_[boundchg.column]) {
      if (col_lower_[boundchg.column] - boundchg.boundval > feastol()) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        if (!infeasible_) {
          infeasible_pos = domchgstack_.size();
          infeasible_ = true;
          infeasible_reason = Reason::conflictingBounds(domchgstack_.size());
        }
      } else {
        boundchg.boundval = col_lower_[boundchg.column];
        if (boundchg.boundval == col_upper_[boundchg.column]) return;
      }
    }

    prevPos = colUpperPos_[boundchg.column];
    colUpperPos_[boundchg.column] = domchgstack_.size();
  }

  mipsolver->mipdata_->debugSolution.boundChangeAdded(
      *this, boundchg, reason.type == Reason::kBranching);

  if (reason.type == Reason::kBranching)
    branchPos_.push_back(domchgstack_.size());

  assert(prevPos < static_cast<HighsInt>(domchgstack_.size()));

  bool binary = isBinary(boundchg.column);

  double oldbound = doChangeBound(boundchg);

  prevboundval_.emplace_back(oldbound, prevPos);
  domchgstack_.push_back(boundchg);
  domchgreason_.push_back(reason);

  if (binary && !infeasible_ && isFixed(boundchg.column)) {
    mipsolver->mipdata_->cliquetable.addImplications(
        *this, boundchg.column, col_lower_[boundchg.column] > 0.5);
    if (!infeasible_) {
      mipsolver->mipdata_->implications.applyImplications(
          *this, boundchg.column, col_lower_[boundchg.column] > 0.5);
    }
  }
}

bool HighsDomain::checkChangeBound(HighsBoundType boundtype, HighsInt col,
                                   HighsCDouble boundval, Reason reason) {
  if (abs(boundval) * kHighsTiny > feastol()) return false;
  bool accept;
  double bound = boundtype == HighsBoundType::kLower
                     ? adjustedLb(col, boundval, accept)
                     : adjustedUb(col, boundval, accept);
  if (!accept) return false;
  changeBound(boundtype, col, bound, reason);
  return true;
}

void HighsDomain::setDomainChangeStack(
    const std::vector<HighsDomainChange>& domchgstack) {
  infeasible_ = false;
  mipsolver->mipdata_->debugSolution.resetDomain(*this);

  if (!domchgstack_.empty()) {
    for (const HighsDomainChange& domchg : domchgstack_) {
      if (domchg.boundtype == HighsBoundType::kLower)
        colLowerPos_[domchg.column] = -1;
      else
        colUpperPos_[domchg.column] = -1;
    }
  }

  prevboundval_.clear();
  domchgstack_.clear();
  domchgreason_.clear();
  branchPos_.clear();
  HighsInt stacksize = static_cast<HighsInt>(domchgstack.size());
  for (HighsInt k = 0; k != stacksize; ++k) {
    if (domchgstack[k].boundtype == HighsBoundType::kLower &&
        domchgstack[k].boundval <= col_lower_[domchgstack[k].column])
      continue;
    if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
        domchgstack[k].boundval >= col_upper_[domchgstack[k].column])
      continue;

    changeBound(domchgstack[k], Reason::unspecified());

    if (infeasible_) break;
  }
}

void HighsDomain::setDomainChangeStack(
    const std::vector<HighsDomainChange>& domchgstack,
    const std::vector<HighsInt>& branchingPositions) {
  infeasible_ = false;
  mipsolver->mipdata_->debugSolution.resetDomain(*this);
  const bool noSymmetry = !mipsolver->mipdata_->detectSymmetries;

  if (!domchgstack_.empty()) {
    for (const HighsDomainChange& domchg : domchgstack_) {
      if (domchg.boundtype == HighsBoundType::kLower)
        colLowerPos_[domchg.column] = -1;
      else
        colUpperPos_[domchg.column] = -1;
    }
  }

  prevboundval_.clear();
  domchgstack_.clear();
  domchgreason_.clear();
  branchPos_.clear();
  HighsInt stacksize = static_cast<HighsInt>(domchgstack.size());
  HighsInt k = 0;
  for (HighsInt branchPos : branchingPositions) {
    for (; k < branchPos; ++k) {
      if (domchgstack[k].boundtype == HighsBoundType::kLower &&
          domchgstack[k].boundval <= col_lower_[domchgstack[k].column])
        continue;
      if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
          domchgstack[k].boundval >= col_upper_[domchgstack[k].column])
        continue;

      changeBound(domchgstack[k], Reason::unspecified());
      if (!infeasible_) propagate();
      if (infeasible_) return;
    }

    if (k == stacksize) return;

    // Remove redundant branching changes unless symmetry is active
    if (noSymmetry) {
      if (domchgstack[k].boundtype == HighsBoundType::kLower) {
        if (domchgstack[k].boundval <= col_lower_[domchgstack[k].column]) {
          if (domchgstack[k].boundval < col_lower_[domchgstack[k].column])
            continue;
          if (colLowerPos_[domchgstack[k].column] != -1) continue;
        }
      } else {
        if (domchgstack[k].boundval >= col_upper_[domchgstack[k].column]) {
          if (domchgstack[k].boundval > col_upper_[domchgstack[k].column])
            continue;
          if (colUpperPos_[domchgstack[k].column] != -1) continue;
        }
      }
    }

    changeBound(domchgstack[k], Reason::branching());
    if (!infeasible_) propagate();
    if (infeasible_) return;
  }

  for (; k < stacksize; ++k) {
    if (domchgstack[k].boundtype == HighsBoundType::kLower &&
        domchgstack[k].boundval <= col_lower_[domchgstack[k].column])
      continue;
    if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
        domchgstack[k].boundval >= col_upper_[domchgstack[k].column])
      continue;

    mipsolver->mipdata_->debugSolution.boundChangeAdded(*this, domchgstack[k],
                                                        true);

    changeBound(domchgstack[k], Reason::unspecified());
    if (!infeasible_) propagate();
    if (infeasible_) break;
  }
}

void HighsDomain::backtrackToGlobal() {
  HighsInt k = static_cast<HighsInt>(domchgstack_.size()) - 1;
  bool old_infeasible = infeasible_;
  Reason old_reason = infeasible_reason;

  if (infeasible_ &&
      infeasible_pos == static_cast<HighsInt>(domchgstack_.size())) {
    assert(old_infeasible);
    assert(k == static_cast<HighsInt>(domchgstack_.size()) - 1);
    infeasible_ = false;
    infeasible_reason = Reason::unspecified();
  }

  while (k >= 0) {
    double prevbound = prevboundval_[k].first;
    HighsInt prevpos = prevboundval_[k].second;
    assert(prevpos < k);

    mipsolver->mipdata_->debugSolution.boundChangeRemoved(*this,
                                                          domchgstack_[k]);

    if (domchgstack_[k].boundtype == HighsBoundType::kLower) {
      assert(colLowerPos_[domchgstack_[k].column] == k);
      colLowerPos_[domchgstack_[k].column] = prevpos;
    } else {
      assert(colUpperPos_[domchgstack_[k].column] == k);
      colUpperPos_[domchgstack_[k].column] = prevpos;
    }

    if (prevbound != domchgstack_[k].boundval) {
      // change back to global bound
      doChangeBound(
          {prevbound, domchgstack_[k].column, domchgstack_[k].boundtype});
    }

    if (infeasible_ && infeasible_pos == k) {
      assert(old_infeasible);
      assert(k == static_cast<HighsInt>(domchgstack_.size()) - 1);
      infeasible_ = false;
      infeasible_reason = Reason::unspecified();
    }

    --k;
  }

  if (old_infeasible) {
    markPropagateCut(old_reason);
    infeasible_reason = Reason::unspecified();
    infeasible_ = false;
  }

  HighsInt numreason = static_cast<HighsInt>(domchgreason_.size());
  for (HighsInt i = k + 1; i < numreason; ++i)
    markPropagateCut(domchgreason_[i]);

  domchgstack_.clear();
  prevboundval_.clear();
  domchgreason_.clear();
  branchPos_.clear();
}

HighsDomainChange HighsDomain::backtrack() {
  HighsInt k = static_cast<HighsInt>(domchgstack_.size()) - 1;
  bool old_infeasible = infeasible_;
  Reason old_reason = infeasible_reason;

  if (infeasible_ &&
      infeasible_pos == static_cast<HighsInt>(domchgstack_.size())) {
    assert(old_infeasible);
    assert(k == static_cast<HighsInt>(domchgstack_.size()) - 1);
    infeasible_ = false;
    infeasible_reason = Reason::unspecified();
  }
  while (k >= 0) {
    double prevbound = prevboundval_[k].first;
    HighsInt prevpos = prevboundval_[k].second;
    assert(prevpos < k);

    mipsolver->mipdata_->debugSolution.boundChangeRemoved(*this,
                                                          domchgstack_[k]);

    if (domchgstack_[k].boundtype == HighsBoundType::kLower) {
      assert(colLowerPos_[domchgstack_[k].column] == k);
      colLowerPos_[domchgstack_[k].column] = prevpos;
    } else {
      assert(colUpperPos_[domchgstack_[k].column] == k);
      colUpperPos_[domchgstack_[k].column] = prevpos;
    }

    // change back to global bound
    doChangeBound(
        {prevbound, domchgstack_[k].column, domchgstack_[k].boundtype});

    if (infeasible_ && infeasible_pos == k) {
      assert(old_infeasible);
      assert(k == static_cast<HighsInt>(domchgstack_.size()) - 1);
      infeasible_ = false;
      infeasible_reason = Reason::unspecified();
    }

    if (domchgreason_[k].type == Reason::kBranching) {
      branchPos_.pop_back();
      break;
    }

    --k;
  }

  if (old_infeasible) {
    markPropagateCut(old_reason);
    infeasible_reason = Reason::unspecified();
    infeasible_ = false;
  }

  HighsInt numreason = static_cast<HighsInt>(domchgreason_.size());
  for (HighsInt i = k + 1; i < numreason; ++i)
    markPropagateCut(domchgreason_[i]);

  if (k < 0) {
    domchgstack_.clear();
    prevboundval_.clear();
    domchgreason_.clear();
    branchPos_.clear();
    return HighsDomainChange({0.0, -1, HighsBoundType::kLower});
  }

  HighsDomainChange backtrackboundchg = domchgstack_[k];
  domchgstack_.erase(domchgstack_.begin() + k, domchgstack_.end());
  domchgreason_.resize(k);
  prevboundval_.resize(k);

  return backtrackboundchg;
}

bool HighsDomain::propagate() {
  std::vector<HighsInt> propagateinds;

  auto havePropagationRows = [&]() {
    if (!propagateinds_.empty()) return true;

    if (objProp_.isActive() && objProp_.shouldBePropagated()) return true;

    for (const auto& cutpoolprop : cutpoolpropagation) {
      if (!cutpoolprop.propagatecutinds_.empty()) return true;
    }

    for (const auto& conflictprop : conflictPoolPropagation) {
      if (!conflictprop.propagateConflictInds_.empty()) return true;
    }

    if (!infeasible_ && dualFixProbingPropagation.isActive()) return true;

    return false;
  };

  if (!havePropagationRows()) return false;

  size_t changedboundsize = 2 * mipsolver->mipdata_->ARvalue_.size();

  for (const auto& cutpoolprop : cutpoolpropagation)
    changedboundsize = std::max(
        changedboundsize, cutpoolprop.cutpool->getMatrix().nonzeroCapacity());

  std::unique_ptr<HighsDomainChange[]> changedbounds(
      new HighsDomainChange[changedboundsize]);

  while (havePropagationRows()) {
    if (objProp_.isActive()) objProp_.propagate();

    const HighsInt numConflictPools =
        static_cast<HighsInt>(conflictPoolPropagation.size());
    for (HighsInt conflictPool = 0; conflictPool < numConflictPools;
         ++conflictPool) {
      auto& conflictprop = conflictPoolPropagation[conflictPool];
      while (!conflictprop.propagateConflictInds_.empty()) {
        propagateinds.swap(conflictprop.propagateConflictInds_);

        for (HighsInt conflict : propagateinds)
          conflictprop.propagateConflict(conflict);

        propagateinds.clear();
      }
    }

    if (!propagateinds_.empty()) {
      propagateinds.swap(propagateinds_);

      HighsInt numproprows = static_cast<HighsInt>(propagateinds.size());
      for (HighsInt i = 0; i != numproprows; ++i) {
        HighsInt row = propagateinds[i];
        propagateflags_[row] = false;
      }

      if (!infeasible_) {
        propRowNumChangedBounds_.assign(
            numproprows, std::make_pair(HighsInt{0}, HighsInt{0}));

        auto propagateIndex = [&](HighsInt k) {
          // for (HighsInt k = 0; k != numproprows; ++k) {
          HighsInt i = propagateinds[k];
          HighsInt start = mipsolver->mipdata_->ARstart_[i];
          HighsInt end = mipsolver->mipdata_->ARstart_[i + 1];
          HighsInt Rlen = end - start;
          const HighsInt* Rindex = mipsolver->mipdata_->ARindex_.data() + start;
          const double* Rvalue = mipsolver->mipdata_->ARvalue_.data() + start;
          bool recomputeCapThreshold = false;

          if (!isRowUpperRedundant(i)) {
            // computeMinActivity(start, end, mipsolver->ARstart_.data(),
            // mipsolver->ARvalue_.data(), activitymininf_[i],
            //           activitymin_[i]);
            activitymin_[i].renormalize();
            propRowNumChangedBounds_[k].first = propagateRowUpper(
                Rindex, Rvalue, Rlen, mipsolver->rowUpper(i), activitymin_[i],
                activitymininf_[i], &changedbounds[2 * start]);

            recomputeCapThreshold = true;
          }

          if (!isRowLowerRedundant(i)) {
            // computeMaxActivity(start, end, mipsolver->ARstart_.data(),
            // mipsolver->ARvalue_.data(), activitymaxinf_[i],
            //           activitymax_[i]);
            activitymax_[i].renormalize();
            propRowNumChangedBounds_[k].second = propagateRowLower(
                Rindex, Rvalue, Rlen, mipsolver->rowLower(i), activitymax_[i],
                activitymaxinf_[i],
                &changedbounds[2 * start + propRowNumChangedBounds_[k].first]);

            recomputeCapThreshold = true;
          }

          if (recomputeCapThreshold) recomputeCapacityThreshold(i);
        };

        // printf("numproprows (model): %" HIGHSINT_FORMAT "\n", numproprows);

        for (HighsInt k = 0; k != numproprows; ++k) propagateIndex(k);

        for (HighsInt k = 0; k != numproprows; ++k) {
          HighsInt i = propagateinds[k];

          if (propRowNumChangedBounds_[k].first != 0) {
            HighsInt start = 2 * mipsolver->mipdata_->ARstart_[i];
            HighsInt end = start + propRowNumChangedBounds_[k].first;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::modelRowUpper(i));

            if (infeasible_) break;
          }
          if (propRowNumChangedBounds_[k].second != 0) {
            HighsInt start = 2 * mipsolver->mipdata_->ARstart_[i] +
                             propRowNumChangedBounds_[k].first;
            HighsInt end = start + propRowNumChangedBounds_[k].second;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::modelRowLower(i));

            if (infeasible_) break;
          }
        }
      }

      propagateinds.clear();
    }

    const HighsInt numpools = static_cast<HighsInt>(cutpoolpropagation.size());
    for (HighsInt cutpool = 0; cutpool != numpools; ++cutpool) {
      auto& cutpoolprop = cutpoolpropagation[cutpool];
      if (!cutpoolprop.propagatecutinds_.empty()) {
        propagateinds.swap(cutpoolprop.propagatecutinds_);

        HighsInt numproprows = static_cast<HighsInt>(propagateinds.size());
        for (HighsInt i = 0; i != numproprows; ++i) {
          HighsInt cut = propagateinds[i];
          cutpoolprop.propagatecutflags_[cut] &= 2;
        }

        if (!infeasible_) {
          propRowNumChangedBounds_.assign(
              numproprows, std::make_pair(HighsInt{0}, HighsInt{0}));

          auto propagateIndex = [&](HighsInt k) {
            HighsInt i = propagateinds[k];
            // first check if cut is marked as deleted
            if (cutpoolprop.propagatecutflags_[i] & 2) return;

            HighsInt Rlen;
            const HighsInt* Rindex;
            const double* Rvalue;
            cutpoolprop.cutpool->getCut(i, Rlen, Rindex, Rvalue);
            cutpoolprop.activitycuts_[i].renormalize();

            propRowNumChangedBounds_[k].first = propagateRowUpper(
                Rindex, Rvalue, Rlen, cutpoolprop.cutpool->getRhs()[i],
                cutpoolprop.activitycuts_[i], cutpoolprop.activitycutsinf_[i],
                &changedbounds[cutpoolprop.cutpool->getMatrix().getRowStart(
                    i)]);

            cutpoolprop.recomputeCapacityThreshold(i);
          };

          // printf("numproprows (cuts): %" HIGHSINT_FORMAT "\n", numproprows);

          for (HighsInt k = 0; k != numproprows; ++k) propagateIndex(k);

          for (HighsInt k = 0; k != numproprows; ++k) {
            HighsInt i = propagateinds[k];
            if (propRowNumChangedBounds_[k].first != 0) {
              cutpoolprop.cutpool->resetAge(
                  i, mipsolver->mipdata_->parallelLockActive());
              HighsInt start = cutpoolprop.cutpool->getMatrix().getRowStart(i);
              HighsInt end = start + propRowNumChangedBounds_[k].first;
              for (HighsInt j = start; j != end && !infeasible_; ++j)
                changeBound(changedbounds[j], Reason::cut(cutpool, i));
            }

            if (infeasible_) break;
          }
        }

        propagateinds.clear();
      }
    }

    if (!infeasible_ && dualFixProbingPropagation.isActive()) {
      dualFixProbingPropagation.propagate();
    }
    if (!infeasible_ && dualFixProbingPropagation.isEnabled() &&
        !havePropagationRows()) {
      dualFixProbingPropagation.propagateZeroCosts();
    }
  }

  return true;
}

double HighsDomain::getColLowerPos(HighsInt col, HighsInt stackpos,
                                   HighsInt& pos) const {
  double lb = col_lower_[col];
  pos = colLowerPos_[col];
  while (pos > stackpos || (pos != -1 && prevboundval_[pos].first == lb)) {
    lb = prevboundval_[pos].first;
    pos = prevboundval_[pos].second;
  }
  return lb;
}

double HighsDomain::getColUpperPos(HighsInt col, HighsInt stackpos,
                                   HighsInt& pos) const {
  double ub = col_upper_[col];
  pos = colUpperPos_[col];
  while (pos > stackpos || (pos != -1 && prevboundval_[pos].first == ub)) {
    ub = prevboundval_[pos].first;
    pos = prevboundval_[pos].second;
  }
  return ub;
}

void HighsDomain::conflictAnalysis(HighsConflictPool& conflictPool,
                                   HighsDomain& globaldom,
                                   HighsPseudocost& pseudocost) {
  if (&globaldom == this) return;
  if (globaldom.infeasible() || !infeasible_) return;

  globaldom.propagate();
  if (globaldom.infeasible()) return;

  ConflictSet conflictSet(*this, globaldom);

  conflictSet.conflictAnalysis(conflictPool, pseudocost);
}

void HighsDomain::conflictAnalysis(const HighsInt* proofinds,
                                   const double* proofvals, HighsInt prooflen,
                                   double proofrhs,
                                   HighsConflictPool& conflictPool,
                                   HighsDomain& globaldom,
                                   HighsPseudocost& pseudocost) {
  if (&globaldom == this) return;

  if (globaldom.infeasible()) return;

  globaldom.propagate();
  if (globaldom.infeasible()) return;

  ConflictSet conflictSet(*this, globaldom);
  conflictSet.conflictAnalysis(proofinds, proofvals, prooflen, proofrhs,
                               conflictPool, pseudocost);
}

void HighsDomain::conflictAnalyzeReconvergence(
    const HighsDomainChange& domchg, const HighsInt* proofinds,
    const double* proofvals, HighsInt prooflen, double proofrhs,
    HighsConflictPool& conflictPool, HighsDomain& globaldom,
    HighsPseudocost& pseudocost) {
  if (&globaldom == this) return;

  if (globaldom.infeasible()) return;

  globaldom.propagate();
  if (globaldom.infeasible()) return;

  ConflictSet conflictSet(*this, globaldom);

  HighsInt ninfmin;
  HighsCDouble activitymin;
  globaldom.computeMinActivity(0, prooflen, proofinds, proofvals, ninfmin,
                               activitymin);
  if (ninfmin != 0) return;

  if (!conflictSet.explainBoundChangeLeq(
          conflictSet.reconvergenceFrontier,
          ConflictSet::LocalDomChg{static_cast<HighsInt>(domchgstack_.size()),
                                   domchg},
          proofinds, proofvals, prooflen, proofrhs,
          static_cast<double>(activitymin)))
    return;

  if (10 * conflictSet.resolvedDomainChanges.size() >
      1000 + 3 * mipsolver->mipdata_->integral_cols.size())
    return;

  conflictSet.reconvergenceFrontier.insert(
      conflictSet.resolvedDomainChanges.begin(),
      conflictSet.resolvedDomainChanges.end());

  HighsInt depth = static_cast<HighsInt>(branchPos_.size());

  while (depth > 0) {
    HighsInt branchPos = branchPos_[depth - 1];
    if (domchgstack_[branchPos].boundval != prevboundval_[branchPos].first)
      break;

    --depth;
  }

  conflictSet.resolveDepth(conflictSet.reconvergenceFrontier, depth, 0,
                           pseudocost);
  conflictPool.addReconvergenceCut(*this, conflictSet.reconvergenceFrontier,
                                   domchg);
}

void HighsDomain::tightenCoefficients(HighsInt* inds, double* vals,
                                      HighsInt len, double& rhs) const {
  HighsCDouble maxactivity = 0;

  for (HighsInt i = 0; i != len; ++i) {
    if (vals[i] > 0) {
      if (col_upper_[inds[i]] == kHighsInf) return;

      maxactivity += static_cast<HighsCDouble>(col_upper_[inds[i]]) * vals[i];
    } else {
      if (col_lower_[inds[i]] == -kHighsInf) return;

      maxactivity += static_cast<HighsCDouble>(col_lower_[inds[i]]) * vals[i];
    }
  }

  HighsCDouble maxabscoef = maxactivity - rhs;
  if (maxabscoef > feastol()) {
    HighsCDouble upper = rhs;
    HighsInt tightened = 0;
    for (HighsInt i = 0; i != len; ++i) {
      if (mipsolver->isColContinuous(inds[i])) continue;
      if (vals[i] > maxabscoef) {
        HighsCDouble delta = vals[i] - maxabscoef;
        upper -= delta * col_upper_[inds[i]];
        vals[i] = static_cast<double>(maxabscoef);
        ++tightened;
      } else if (vals[i] < -maxabscoef) {
        HighsCDouble delta = -vals[i] - maxabscoef;
        upper += delta * col_lower_[inds[i]];
        vals[i] = -static_cast<double>(maxabscoef);
        ++tightened;
      }
    }

    if (tightened != 0) {
      // printf("tightened %" HIGHSINT_FORMAT " coefficients, rhs changed from
      // %g to %g\n",
      //       tightened, rhs, static_cast<double>(upper));
      rhs = static_cast<double>(upper);
    }
  }
}

double HighsDomain::getMinCutActivity(const HighsCutPool& cutpool,
                                      HighsInt cut) const {
  for (auto& cutpoolprop : cutpoolpropagation) {
    if (cutpoolprop.cutpool == &cutpool) {
      // assert((cutpoolprop.propagatecutflags_[cut] & 2) == 0);

      return cut < static_cast<HighsInt>(
                       cutpoolprop.propagatecutflags_.size()) &&
                     (cutpoolprop.propagatecutflags_[cut] & 2) == 0 &&
                     cutpoolprop.activitycutsinf_[cut] == 0
                 ? static_cast<double>(cutpoolprop.activitycuts_[cut])
                 : -kHighsInf;
    }
  }

  return -kHighsInf;
}

bool HighsDomain::isFixing(const HighsDomainChange& domchg) const {
  double otherbound = domchg.boundtype == HighsBoundType::kUpper
                          ? col_lower_[domchg.column]
                          : col_upper_[domchg.column];
  return std::fabs(domchg.boundval - otherbound) <= epsilon();
}

HighsDomainChange HighsDomain::flip(const HighsDomainChange& domchg) const {
  if (domchg.boundtype == HighsBoundType::kLower) {
    HighsDomainChange flipped{domchg.boundval - feastol(), domchg.column,
                              HighsBoundType::kUpper};
    if (mipsolver->isColIntegral(domchg.column))
      flipped.boundval = std::floor(flipped.boundval);
    return flipped;
  } else {
    HighsDomainChange flipped{domchg.boundval + feastol(), domchg.column,
                              HighsBoundType::kLower};
    if (mipsolver->isColIntegral(domchg.column))
      flipped.boundval = std::ceil(flipped.boundval);
    return flipped;
  }
}

double HighsDomain::feastol() const { return mipsolver->mipdata_->feastol; }

double HighsDomain::epsilon() const { return mipsolver->mipdata_->epsilon; }

HighsDomain::ConflictSet::ConflictSet(HighsDomain& localdom_,
                                      const HighsDomain& globaldom_)
    : localdom(localdom_),
      globaldom(globaldom_),
      reasonSideFrontier(),
      reconvergenceFrontier(),
      resolveQueue(),
      resolvedDomainChanges() {}

static inline double computePrio(double val, double bound, double globalbound,
                                 int64_t numnodes) {
  return std::fabs(val * (bound - globalbound) *
                   static_cast<double>(1 + numnodes));
}

bool HighsDomain::ConflictSet::explainBoundChangeGeq(
    const std::set<LocalDomChg>& currentFrontier, const LocalDomChg& domchg,
    const HighsInt* inds, const double* vals, HighsInt len, double rhs,
    double maxAct) {
  if (maxAct == kHighsInf || std::abs(maxAct) > kHighsInvMacheps) return false;

  // get the coefficient value of the column for which we want to explain
  // the bound change
  double domchgVal = 0;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  HighsCDouble M = maxAct;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    if (col == domchg.domchg.column) {
      domchgVal = vals[i];
      continue;
    }

    ResolveCandidate cand;
    cand.valuePos = i;

    if (vals[i] > 0) {
      double ub = localdom.getColUpperPos(col, domchg.pos, cand.boundPos);
      if (globaldom.col_upper_[col] <= ub || cand.boundPos == -1) continue;
      auto it =
          currentFrontier.find(LocalDomChg{cand.boundPos, HighsDomainChange()});
      if (it != currentFrontier.end()) {
        cand.baseBound = it->domchg.boundval;
        if (cand.baseBound != globaldom.col_upper_[col])
          M += vals[i] * (static_cast<HighsCDouble>(cand.baseBound) -
                          globaldom.col_upper_[col]);
        if (cand.baseBound <= ub) continue;
      } else
        cand.baseBound = globaldom.col_upper_[col];

      cand.delta = vals[i] * (ub - cand.baseBound);
      cand.prio = computePrio(vals[i], ub, globaldom.col_upper_[col],
                              nodequeue.numNodesDown(col));
    } else {
      double lb = localdom.getColLowerPos(col, domchg.pos, cand.boundPos);
      if (globaldom.col_lower_[col] >= lb || cand.boundPos == -1) continue;
      auto it =
          currentFrontier.find(LocalDomChg{cand.boundPos, HighsDomainChange()});

      if (it != currentFrontier.end()) {
        cand.baseBound = it->domchg.boundval;
        if (cand.baseBound != globaldom.col_lower_[col])
          M += vals[i] * (static_cast<HighsCDouble>(cand.baseBound) -
                          globaldom.col_lower_[col]);
        if (cand.baseBound >= lb) continue;
      } else
        cand.baseBound = globaldom.col_lower_[col];

      cand.delta = vals[i] * (lb - cand.baseBound);
      cand.prio = computePrio(vals[i], lb, globaldom.col_lower_[col],
                              nodequeue.numNodesUp(col));
    }

    resolveBuffer.push_back(cand);
  }

  if (domchgVal == 0) return false;

  pdqsort(resolveBuffer.begin(), resolveBuffer.end());

  // to explain the bound change we start from the bound constraint,
  // multiply it by the columns coefficient in the constraint. Then the
  // bound constraint is a0 * x0 >= b0 and the constraint
  // a0 * x0 + \sum ai * xi >= b. Let the max activity of \sum ai xi be M,
  // then the constraint yields the bound a0 * x0 >= b - M. If M is
  // sufficiently small the bound constraint is implied. Therefore we
  // decrease M by updating it with the stronger local bounds until
  // M <= b - b0 holds.
  double b0 = domchg.domchg.boundval;
  if (!localdom.mipsolver->isColContinuous(domchg.domchg.column)) {
    // in case of an integral variable the bound was rounded and can be
    // relaxed by 1-feastol. We use 1 - 10 * feastol for numerical safety.
    if (domchg.domchg.boundtype == HighsBoundType::kLower)
      b0 -= (1.0 - 10 * localdom.feastol());
    else
      b0 += (1.0 - 10 * localdom.feastol());
  } else {
    // for a continuous variable we relax the bound by epsilon to
    // accommodate for tiny rounding errors
    if (domchg.domchg.boundtype == HighsBoundType::kLower)
      b0 -= localdom.epsilon();
    else
      b0 += localdom.epsilon();
  }

  // now multiply the bound constraint with the coefficient value in the
  // constraint to obtain b0
  b0 *= domchgVal;

  // compute the lower bound of M that is necessary
  double Mupper = rhs - b0;

  // M is the current residual activity initially
  if (domchgVal < 0)
    M -= static_cast<HighsCDouble>(domchgVal) *
         globaldom.col_lower_[domchg.domchg.column];
  else
    M -= static_cast<HighsCDouble>(domchgVal) *
         globaldom.col_upper_[domchg.domchg.column];

  return resolveLinearGeq(M, Mupper, vals);
}

bool HighsDomain::ConflictSet::explainBoundChangeLeq(
    const std::set<LocalDomChg>& currentFrontier, const LocalDomChg& domchg,
    const HighsInt* inds, const double* vals, HighsInt len, double rhs,
    double minAct) {
  if (minAct == -kHighsInf || std::abs(minAct) > kHighsInvMacheps) return false;
  // get the coefficient value of the column for which we want to explain
  // the bound change
  double domchgVal = 0;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  HighsCDouble M = minAct;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    if (col == domchg.domchg.column) {
      domchgVal = vals[i];
      continue;
    }

    ResolveCandidate cand;
    cand.valuePos = i;

    if (vals[i] > 0) {
      double lb = localdom.getColLowerPos(col, domchg.pos, cand.boundPos);
      if (globaldom.col_lower_[col] >= lb || cand.boundPos == -1) continue;
      auto it =
          currentFrontier.find(LocalDomChg{cand.boundPos, HighsDomainChange()});

      if (it != currentFrontier.end()) {
        cand.baseBound = it->domchg.boundval;
        if (cand.baseBound != globaldom.col_lower_[col])
          M += vals[i] * (static_cast<HighsCDouble>(cand.baseBound) -
                          globaldom.col_lower_[col]);
        if (cand.baseBound >= lb) continue;
      } else
        cand.baseBound = globaldom.col_lower_[col];

      cand.delta = vals[i] * (lb - cand.baseBound);
      cand.prio = computePrio(vals[i], lb, globaldom.col_lower_[col],
                              nodequeue.numNodesUp(col));
    } else {
      double ub = localdom.getColUpperPos(col, domchg.pos, cand.boundPos);
      if (globaldom.col_upper_[col] <= ub || cand.boundPos == -1) continue;
      auto it =
          currentFrontier.find(LocalDomChg{cand.boundPos, HighsDomainChange()});
      if (it != currentFrontier.end()) {
        cand.baseBound = it->domchg.boundval;
        if (cand.baseBound != globaldom.col_upper_[col])
          M += vals[i] * (static_cast<HighsCDouble>(cand.baseBound) -
                          globaldom.col_upper_[col]);
        if (cand.baseBound <= ub) continue;
      } else
        cand.baseBound = globaldom.col_upper_[col];

      cand.delta = vals[i] * (ub - cand.baseBound);
      cand.prio = computePrio(vals[i], ub, globaldom.col_upper_[col],
                              nodequeue.numNodesDown(col));
    }

    resolveBuffer.push_back(cand);
  }

  if (domchgVal == 0) return false;

  pdqsort(resolveBuffer.begin(), resolveBuffer.end());

  assert(domchgVal != 0);

  // to explain the bound change we start from the bound constraint,
  // multiply it by the columns coefficient in the constraint. Then the
  // bound constraint is a0 * x0 <= b0 and the constraint
  // a0 * x0 + \sum ai * xi <= b. Let the min activity of \sum ai xi be M,
  // then the constraint yields the bound a0 * x0 <= b - M. If M is
  // sufficiently large the bound constraint is implied. Therefore we
  // increase M by updating it with the stronger local bounds until
  // M >= b - b0 holds.
  double b0 = domchg.domchg.boundval;
  if (!localdom.mipsolver->isColContinuous(domchg.domchg.column)) {
    // in case of an integral variable the bound was rounded and can be
    // relaxed by 1-feastol. We use 1 - 10 * feastol for numerical safety
    if (domchg.domchg.boundtype == HighsBoundType::kLower)
      b0 -= (1.0 - 10 * localdom.feastol());
    else
      b0 += (1.0 - 10 * localdom.feastol());
  } else {
    // for a continuous variable we relax the bound by epsilon to
    // accommodate for tiny rounding errors
    if (domchg.domchg.boundtype == HighsBoundType::kLower)
      b0 -= localdom.epsilon();
    else
      b0 += localdom.epsilon();
  }

  // now multiply the bound constraint with the coefficient value in the
  // constraint to obtain b0
  b0 *= domchgVal;

  // compute the lower bound of M that is necessary
  double Mlower = rhs - b0;

  // M is the global residual activity initially
  if (domchgVal < 0)
    M -= static_cast<HighsCDouble>(domchgVal) *
         globaldom.col_upper_[domchg.domchg.column];
  else
    M -= static_cast<HighsCDouble>(domchgVal) *
         globaldom.col_lower_[domchg.domchg.column];

  return resolveLinearLeq(M, Mlower, vals);
}

bool HighsDomain::ConflictSet::resolveLinearGeq(HighsCDouble M, double Mupper,
                                                const double* vals) {
  resolvedDomainChanges.clear();
  double covered = static_cast<double>(M - Mupper);
  if (covered > 0) {
    for (const auto& reasonDomchg : resolveBuffer) {
      LocalDomChg locdomchg;
      locdomchg.pos = reasonDomchg.boundPos;
      locdomchg.domchg = localdom.domchgstack_[reasonDomchg.boundPos];
      M += reasonDomchg.delta;
      resolvedDomainChanges.push_back(locdomchg);
      assert(resolvedDomainChanges.back().pos >= 0);
      assert(resolvedDomainChanges.back().pos <
             static_cast<HighsInt>(localdom.domchgstack_.size()));
      covered = static_cast<double>(M - Mupper);
      if (covered <= 0) break;
    }

    if (covered > 0) return false;

    if (covered < -localdom.feastol()) {
      // there is room for relaxing bounds / dropping unneeded bound changes
      // from the explanation
      // HighsInt numRelaxed = 0;
      // HighsInt numDropped = 0;
      for (HighsInt k = static_cast<HighsInt>(resolvedDomainChanges.size()) - 1;
           k >= 0; --k) {
        ResolveCandidate& reasonDomchg = resolveBuffer[k];
        LocalDomChg& locdomchg = resolvedDomainChanges[k];
        HighsInt i = reasonDomchg.valuePos;
        HighsInt col = locdomchg.domchg.column;
        if (locdomchg.domchg.boundtype == HighsBoundType::kLower) {
          double lb = locdomchg.domchg.boundval;
          double glb = reasonDomchg.baseBound;
          double relaxLb = static_cast<double>(
              ((Mupper - (M - reasonDomchg.delta)) / vals[i]) + glb);
          if (!localdom.mipsolver->isColContinuous(col))
            relaxLb = std::ceil(relaxLb);

          if (relaxLb - lb >= -localdom.feastol()) continue;

          locdomchg.domchg.boundval = relaxLb;

          // only drop if explanation stays valid without this entry
          if (relaxLb - glb <= localdom.epsilon() &&
              M - reasonDomchg.delta - Mupper <= localdom.feastol()) {
            // domain change can be fully removed from conflict
            HighsInt last =
                static_cast<HighsInt>(resolvedDomainChanges.size()) - 1;
            std::swap(resolvedDomainChanges[last], resolvedDomainChanges[k]);
            resolvedDomainChanges.resize(last);

            M -= reasonDomchg.delta;
            // ++numDropped;
          } else {
            if (!localdom.walkBoundChain(locdomchg.pos, relaxLb, HighsInt{1}))
              return false;

            // bound can be relaxed
            M += vals[i] * (static_cast<HighsCDouble>(relaxLb) - lb);
            // ++numRelaxed;
          }

          covered = static_cast<double>(M - Mupper);
          if (covered >= -localdom.feastol()) break;
        } else {
          double ub = locdomchg.domchg.boundval;
          double gub = reasonDomchg.baseBound;
          double relaxUb = static_cast<double>(
              ((Mupper - (M - reasonDomchg.delta)) / vals[i]) + gub);
          if (!localdom.mipsolver->isColContinuous(col))
            relaxUb = std::floor(relaxUb);

          if (relaxUb - ub <= localdom.feastol()) continue;

          locdomchg.domchg.boundval = relaxUb;

          // only drop if explanation stays valid without this entry
          if (relaxUb - gub >= -localdom.epsilon() &&
              M - reasonDomchg.delta - Mupper <= localdom.feastol()) {
            // domain change can be fully removed from conflict
            HighsInt last =
                static_cast<HighsInt>(resolvedDomainChanges.size()) - 1;
            std::swap(resolvedDomainChanges[last], resolvedDomainChanges[k]);
            resolvedDomainChanges.resize(last);

            M -= reasonDomchg.delta;
            // ++numDropped;
          } else {
            // bound can be relaxed
            if (!localdom.walkBoundChain(locdomchg.pos, relaxUb, HighsInt{-1}))
              return false;

            M += vals[i] * (static_cast<HighsCDouble>(relaxUb) - ub);
            // ++numRelaxed;
          }

          covered = static_cast<double>(M - Mupper);
          if (covered >= -localdom.feastol()) break;
        }
      }

      // if (numRelaxed + numDropped)
      //   printf("relaxed %d and dropped %d of %d resolved domain changes\n",
      //          (int)numRelaxed, (int)numDropped,
      //          (int)resolvedDomainChanges.size());

      assert(covered <= localdom.feastol());
    }
  }

  return true;
}

bool HighsDomain::ConflictSet::resolveLinearLeq(HighsCDouble M, double Mlower,
                                                const double* vals) {
  resolvedDomainChanges.clear();
  double covered = static_cast<double>(M - Mlower);
  if (covered < 0) {
    for (const auto& reasonDomchg : resolveBuffer) {
      LocalDomChg locdomchg;
      locdomchg.pos = reasonDomchg.boundPos;
      locdomchg.domchg = localdom.domchgstack_[reasonDomchg.boundPos];
      M += reasonDomchg.delta;
      resolvedDomainChanges.push_back(locdomchg);
      assert(resolvedDomainChanges.back().pos >= 0);
      assert(resolvedDomainChanges.back().pos <
             static_cast<HighsInt>(localdom.domchgstack_.size()));
      covered = static_cast<double>(M - Mlower);
      if (covered >= 0) break;
    }

    if (covered < 0) {
      // printf("local bounds reach only value of %.12g, need at least
      // %.12g\n",
      //        M, Mlower);
      return false;
    }

    if (covered > localdom.feastol()) {
      // there is room for relaxing bounds / dropping unneeded bound changes
      // from the explanation
      // HighsInt numRelaxed = 0;
      // HighsInt numDropped = 0;
      for (HighsInt k = static_cast<HighsInt>(resolvedDomainChanges.size()) - 1;
           k >= 0; --k) {
        ResolveCandidate& reasonDomchg = resolveBuffer[k];
        LocalDomChg& locdomchg = resolvedDomainChanges[k];
        HighsInt i = reasonDomchg.valuePos;
        HighsInt col = locdomchg.domchg.column;
        if (locdomchg.domchg.boundtype == HighsBoundType::kLower) {
          double lb = locdomchg.domchg.boundval;
          double glb = reasonDomchg.baseBound;
          double relaxLb = static_cast<double>(
              ((Mlower - (M - reasonDomchg.delta)) / vals[i]) + glb);
          if (!localdom.mipsolver->isColContinuous(col))
            relaxLb = std::ceil(relaxLb);

          if (relaxLb - lb >= -localdom.feastol()) continue;

          locdomchg.domchg.boundval = relaxLb;

          // only drop if explanation stays valid without this entry
          if (relaxLb - glb <= localdom.epsilon() &&
              M - reasonDomchg.delta - Mlower >= -localdom.feastol()) {
            // domain change can be fully removed from conflict
            HighsInt last =
                static_cast<HighsInt>(resolvedDomainChanges.size()) - 1;
            std::swap(resolvedDomainChanges[last], resolvedDomainChanges[k]);
            resolvedDomainChanges.resize(last);

            M -= reasonDomchg.delta;
            // ++numDropped;
          } else {
            // bound can be relaxed
            if (!localdom.walkBoundChain(locdomchg.pos, relaxLb, HighsInt{1}))
              return false;

            M += vals[i] * (static_cast<HighsCDouble>(relaxLb) - lb);
            // ++numRelaxed;
          }

          covered = static_cast<double>(M - Mlower);
          if (covered <= localdom.feastol()) break;
        } else {
          double ub = locdomchg.domchg.boundval;
          double gub = reasonDomchg.baseBound;
          double relaxUb = static_cast<double>(
              ((Mlower - (M - reasonDomchg.delta)) / vals[i]) + gub);
          if (!localdom.mipsolver->isColContinuous(col))
            relaxUb = std::floor(relaxUb);

          if (relaxUb - ub <= localdom.feastol()) continue;

          locdomchg.domchg.boundval = relaxUb;

          // only drop if explanation stays valid without this entry
          if (relaxUb - gub >= -localdom.epsilon() &&
              M - reasonDomchg.delta - Mlower >= -localdom.feastol()) {
            // domain change can be fully removed from conflict
            HighsInt last =
                static_cast<HighsInt>(resolvedDomainChanges.size()) - 1;
            std::swap(resolvedDomainChanges[last], resolvedDomainChanges[k]);
            resolvedDomainChanges.resize(last);

            M -= reasonDomchg.delta;
            // ++numDropped;
          } else {
            // bound can be relaxed
            if (!localdom.walkBoundChain(locdomchg.pos, relaxUb, HighsInt{-1}))
              return false;

            M += vals[i] * (static_cast<HighsCDouble>(relaxUb) - ub);
            // ++numRelaxed;
          }

          covered = static_cast<double>(M - Mlower);
          if (covered <= localdom.feastol()) break;
        }
      }

      // if (numRelaxed + numDropped)
      //   printf("relaxed %d and dropped %d of %d resolved domain changes\n",
      //          (int)numRelaxed, (int)numDropped,
      //          (int)resolvedDomainChanges.size());

      assert(covered >= -localdom.feastol());
    }
  }
  return true;
}

bool HighsDomain::ConflictSet::explainInfeasibility() {
  switch (localdom.infeasible_reason.type) {
    case Reason::kUnknown:
    case Reason::kBranching:
      return false;
    case Reason::kConflictingBounds: {
      resolvedDomainChanges.clear();
      HighsInt conflictingBoundPos = localdom.infeasible_reason.index;
      HighsInt col = localdom.domchgstack_[conflictingBoundPos].column;

      resolvedDomainChanges.push_back(LocalDomChg{
          conflictingBoundPos, localdom.domchgstack_[conflictingBoundPos]});

      HighsInt otherBoundPos;
      if (localdom.domchgstack_[conflictingBoundPos].boundtype ==
          HighsBoundType::kLower) {
        double ub =
            localdom.getColUpperPos(col, conflictingBoundPos, otherBoundPos);
        assert(localdom.domchgstack_[conflictingBoundPos].boundval - ub >
               +localdom.feastol());
      } else {
        double lb =
            localdom.getColLowerPos(col, conflictingBoundPos, otherBoundPos);
        assert(localdom.domchgstack_[conflictingBoundPos].boundval - lb <
               -localdom.feastol());
      }
      if (otherBoundPos != -1)
        resolvedDomainChanges.push_back(
            LocalDomChg{otherBoundPos, localdom.domchgstack_[otherBoundPos]});
      return true;
    }
    case Reason::kCliqueTable:
      assert(false);
      return false;
    case Reason::kModelRowLower: {
      HighsInt rowIndex = localdom.infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double maxAct = globaldom.getMaxActivity(rowIndex);

      return explainInfeasibilityGeq(
          inds, vals, len, localdom.mipsolver->rowLower(rowIndex), maxAct);
    }
    case Reason::kModelRowUpper: {
      HighsInt rowIndex = localdom.infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double minAct = globaldom.getMinActivity(rowIndex);

      return explainInfeasibilityLeq(
          inds, vals, len, localdom.mipsolver->rowUpper(rowIndex), minAct);
    }
    case Reason::kObjective: {
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      double rhs;

      localdom.objProp_.getPropagationConstraint(localdom.infeasible_pos, vals,
                                                 inds, len, rhs);

      HighsInt ninfmin;
      HighsCDouble minAct;
      globaldom.computeMinActivity(0, len, inds, vals, ninfmin, minAct);
      // This case should only happen when a globally unbounded column is
      // bounded in the local domain, e.g., by a branching choice in some
      // heuristic.
      if (ninfmin > 0) return false;

      return explainInfeasibilityLeq(inds, vals, len, rhs,
                                     static_cast<double>(minAct));
    }
    default:
      assert(localdom.infeasible_reason.type >= 0);
      assert(localdom.infeasible_reason.type <
             static_cast<HighsInt>(localdom.cutpoolpropagation.size() +
                                   localdom.conflictPoolPropagation.size()));

      if (localdom.infeasible_reason.type <
          static_cast<HighsInt>(localdom.cutpoolpropagation.size())) {
        HighsInt cutpoolIndex = localdom.infeasible_reason.type;
        HighsInt cutIndex = localdom.infeasible_reason.index;

        // retrieve the matrix values of the cut
        HighsInt len;
        const HighsInt* inds;
        const double* vals;
        localdom.cutpoolpropagation[cutpoolIndex].cutpool->getCut(cutIndex, len,
                                                                  inds, vals);

        double minAct = globaldom.getMinCutActivity(
            *localdom.cutpoolpropagation[cutpoolIndex].cutpool, cutIndex);

        return explainInfeasibilityLeq(inds, vals, len,
                                       localdom.cutpoolpropagation[cutpoolIndex]
                                           .cutpool->getRhs()[cutIndex],
                                       minAct);
      } else {
        HighsInt conflictPoolIndex = localdom.infeasible_reason.type -
                                     localdom.cutpoolpropagation.size();
        HighsInt conflictIndex = localdom.infeasible_reason.index;

        if (localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictFlag_[conflictIndex] &
            8)
          return false;

        // retrieve the conflict entries
        auto conflictRange =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictRanges()[conflictIndex];
        const HighsDomainChange* conflict =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictEntryVector()
                .data() +
            conflictRange.first;
        HighsInt len = conflictRange.second - conflictRange.first;

        return explainInfeasibilityConflict(conflict, len);
      }
  }
}

bool HighsDomain::ConflictSet::explainInfeasibilityConflict(
    const HighsDomainChange* conflict, HighsInt len) {
  resolvedDomainChanges.clear();
  for (HighsInt i = 0; i < len; ++i) {
    if (globaldom.isActive(conflict[i])) continue;

    HighsInt pos;
    if (conflict[i].boundtype == HighsBoundType::kLower) {
      double lb = localdom.getColLowerPos(conflict[i].column,
                                          localdom.infeasible_pos, pos);
      if (pos == -1 || lb < conflict[i].boundval) return false;

      if (!localdom.walkBoundChain(pos, conflict[i].boundval, HighsInt{1}))
        return false;
    } else {
      double ub = localdom.getColUpperPos(conflict[i].column,
                                          localdom.infeasible_pos, pos);
      if (pos == -1 || ub > conflict[i].boundval) return false;

      if (!localdom.walkBoundChain(pos, conflict[i].boundval, HighsInt{-1}))
        return false;
    }

    resolvedDomainChanges.push_back(LocalDomChg{pos, conflict[i]});
  }

  return true;
}

bool HighsDomain::ConflictSet::explainInfeasibilityGeq(const HighsInt* inds,
                                                       const double* vals,
                                                       HighsInt len, double rhs,
                                                       double maxAct) {
  if (maxAct == kHighsInf || std::abs(maxAct) > kHighsInvMacheps) return false;

  HighsInt infeasible_pos = kHighsIInf;
  if (localdom.infeasible_) infeasible_pos = localdom.infeasible_pos;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    ResolveCandidate cand;
    cand.valuePos = i;

    if (vals[i] > 0) {
      double ub = localdom.getColUpperPos(col, infeasible_pos, cand.boundPos);
      cand.baseBound = globaldom.col_upper_[col];
      if (cand.baseBound <= ub || cand.boundPos == -1) continue;
      cand.delta = vals[i] * (ub - cand.baseBound);
      cand.prio = computePrio(vals[i], ub, globaldom.col_upper_[col],
                              nodequeue.numNodesDown(col));
    } else {
      double lb = localdom.getColLowerPos(col, infeasible_pos, cand.boundPos);
      cand.baseBound = globaldom.col_lower_[col];
      if (cand.baseBound >= lb || cand.boundPos == -1) continue;
      cand.delta = vals[i] * (lb - cand.baseBound);
      cand.prio = computePrio(vals[i], lb, globaldom.col_lower_[col],
                              nodequeue.numNodesUp(col));
    }

    resolveBuffer.push_back(cand);
  }

  pdqsort(resolveBuffer.begin(), resolveBuffer.end());

  // compute the lower bound of M that is necessary
  double Mupper = rhs - std::max(10.0, std::fabs(rhs)) * localdom.feastol();

  assert(reasonSideFrontier.empty());
  return resolveLinearGeq(maxAct, Mupper, vals);
}

bool HighsDomain::ConflictSet::explainInfeasibilityLeq(const HighsInt* inds,
                                                       const double* vals,
                                                       HighsInt len, double rhs,
                                                       double minAct) {
  if (minAct == -kHighsInf || std::abs(minAct) > kHighsInvMacheps) return false;

  HighsInt infeasible_pos = kHighsIInf;
  if (localdom.infeasible_) infeasible_pos = localdom.infeasible_pos;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    ResolveCandidate cand;
    cand.valuePos = i;

    if (vals[i] > 0) {
      double lb = localdom.getColLowerPos(col, infeasible_pos, cand.boundPos);
      cand.baseBound = globaldom.col_lower_[col];
      if (cand.baseBound >= lb || cand.boundPos == -1) continue;
      cand.delta = vals[i] * (lb - cand.baseBound);
      cand.prio = computePrio(vals[i], lb, globaldom.col_lower_[col],
                              nodequeue.numNodesUp(col));
    } else {
      double ub = localdom.getColUpperPos(col, infeasible_pos, cand.boundPos);
      cand.baseBound = globaldom.col_upper_[col];
      if (cand.baseBound <= ub || cand.boundPos == -1) continue;
      cand.delta = vals[i] * (ub - cand.baseBound);
      cand.prio = computePrio(vals[i], ub, globaldom.col_upper_[col],
                              nodequeue.numNodesDown(col));
    }

    resolveBuffer.push_back(cand);
  }

  pdqsort(resolveBuffer.begin(), resolveBuffer.end());

  // compute the lower bound of M that is necessary
  double Mlower = rhs + std::max(10.0, std::fabs(rhs)) * localdom.feastol();

  return resolveLinearLeq(minAct, Mlower, vals);
}

bool HighsDomain::ConflictSet::explainBoundChange(
    const std::set<LocalDomChg>& currentFrontier, LocalDomChg domchg) {
  switch (localdom.domchgreason_[domchg.pos].type) {
    case Reason::kUnknown:
    case Reason::kBranching:
    case Reason::kConflictingBounds:
      return false;
    case Reason::kCliqueTable: {
      HighsInt col = localdom.domchgreason_[domchg.pos].index >> 1;
      HighsInt val = localdom.domchgreason_[domchg.pos].index & 1;
      resolvedDomainChanges.clear();
      HighsInt boundPos;
      if (val) {
        assert(localdom.colLowerPos_[col] >= 0);
        assert(localdom.colLowerPos_[col] <
               static_cast<HighsInt>(localdom.domchgstack_.size()));

        localdom.getColLowerPos(col, domchg.pos, boundPos);
      } else {
        assert(localdom.colUpperPos_[col] >= 0);
        assert(localdom.colUpperPos_[col] <
               static_cast<HighsInt>(localdom.domchgstack_.size()));

        localdom.getColUpperPos(col, domchg.pos, boundPos);
      }

      if (boundPos != -1)
        resolvedDomainChanges.push_back(
            LocalDomChg{boundPos, localdom.domchgstack_[boundPos]});

      return true;
    }
    case Reason::kModelRowLower: {
      HighsInt rowIndex = localdom.domchgreason_[domchg.pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double maxAct = globaldom.getMaxActivity(rowIndex);

      return explainBoundChangeGeq(currentFrontier, domchg, inds, vals, len,
                                   localdom.mipsolver->rowLower(rowIndex),
                                   maxAct);
    }
    case Reason::kModelRowUpper: {
      HighsInt rowIndex = localdom.domchgreason_[domchg.pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double minAct = globaldom.getMinActivity(rowIndex);

      return explainBoundChangeLeq(currentFrontier, domchg, inds, vals, len,
                                   localdom.mipsolver->rowUpper(rowIndex),
                                   minAct);
    }
    case Reason::kObjective: {
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      double rhs;

      localdom.objProp_.getPropagationConstraint(domchg.pos, vals, inds, len,
                                                 rhs, domchg.domchg.column);

      HighsInt ninfmin;
      HighsCDouble minAct;
      globaldom.computeMinActivity(0, len, inds, vals, ninfmin, minAct);
      // todo: treat case with a single infinite contribution
      if (ninfmin > 0) return false;

      return explainBoundChangeLeq(currentFrontier, domchg, inds, vals, len,
                                   rhs, static_cast<double>(minAct));
    }
    default:
      assert(localdom.domchgreason_[domchg.pos].type >= 0);
      assert(localdom.domchgreason_[domchg.pos].type <
             static_cast<HighsInt>(localdom.cutpoolpropagation.size() +
                                   localdom.conflictPoolPropagation.size()));
      if (localdom.domchgreason_[domchg.pos].type <
          static_cast<HighsInt>(localdom.cutpoolpropagation.size())) {
        HighsInt cutpoolIndex = localdom.domchgreason_[domchg.pos].type;
        HighsInt cutIndex = localdom.domchgreason_[domchg.pos].index;

        // retrieve the matrix values of the cut
        HighsInt len;
        const HighsInt* inds;
        const double* vals;
        localdom.cutpoolpropagation[cutpoolIndex].cutpool->getCut(cutIndex, len,
                                                                  inds, vals);
        double minAct = globaldom.getMinCutActivity(
            *localdom.cutpoolpropagation[cutpoolIndex].cutpool, cutIndex);

        return explainBoundChangeLeq(currentFrontier, domchg, inds, vals, len,
                                     localdom.cutpoolpropagation[cutpoolIndex]
                                         .cutpool->getRhs()[cutIndex],
                                     minAct);
      } else {
        HighsInt conflictPoolIndex = localdom.domchgreason_[domchg.pos].type -
                                     localdom.cutpoolpropagation.size();
        HighsInt conflictIndex = localdom.domchgreason_[domchg.pos].index;

        if (localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictFlag_[conflictIndex] &
            8)
          break;

        // retrieve the conflict entries
        auto conflictRange =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictRanges()[conflictIndex];
        const HighsDomainChange* conflict =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictEntryVector()
                .data() +
            conflictRange.first;
        HighsInt len = conflictRange.second - conflictRange.first;

        return explainBoundChangeConflict(domchg, conflict, len);
      }
  }

  return false;
}

bool HighsDomain::ConflictSet::explainBoundChangeConflict(
    const LocalDomChg& locdomchg, const HighsDomainChange* conflict,
    HighsInt len) {
  resolvedDomainChanges.clear();
  auto domchg = localdom.flip(locdomchg.domchg);
  bool foundDomchg = false;
  for (HighsInt i = 0; i < len; ++i) {
    if (!foundDomchg && conflict[i].column == domchg.column &&
        conflict[i].boundtype == domchg.boundtype) {
      if (conflict[i].boundtype == HighsBoundType::kLower) {
        if (conflict[i].boundval <= domchg.boundval) {
          foundDomchg = true;
          continue;
        }
      } else {
        if (conflict[i].boundval >= domchg.boundval) {
          foundDomchg = true;
          continue;
        }
      }
    }
    if (globaldom.isActive(conflict[i])) continue;

    HighsInt pos;
    if (conflict[i].boundtype == HighsBoundType::kLower) {
      double lb =
          localdom.getColLowerPos(conflict[i].column, locdomchg.pos - 1, pos);

      if (pos == -1 || lb < conflict[i].boundval) return false;

      if (!localdom.walkBoundChain(pos, conflict[i].boundval, HighsInt{1}))
        return false;
    } else {
      double ub =
          localdom.getColUpperPos(conflict[i].column, locdomchg.pos - 1, pos);

      if (pos == -1 || ub > conflict[i].boundval) return false;

      if (!localdom.walkBoundChain(pos, conflict[i].boundval, HighsInt{-1}))
        return false;
    }

    resolvedDomainChanges.push_back(
        LocalDomChg{pos, localdom.domchgstack_[pos]});
  }

  return foundDomchg;
}

void HighsDomain::ConflictSet::pushQueue(
    std::set<LocalDomChg>::iterator domchg) {
  resolveQueue.emplace_back(domchg);
  std::push_heap(
      resolveQueue.begin(), resolveQueue.end(),
      [](const std::set<LocalDomChg>::iterator& a,
         const std::set<LocalDomChg>::iterator& b) { return a->pos < b->pos; });
}

std::set<HighsDomain::ConflictSet::LocalDomChg>::iterator
HighsDomain::ConflictSet::popQueue() {
  assert(!resolveQueue.empty());
  std::pop_heap(
      resolveQueue.begin(), resolveQueue.end(),
      [](const std::set<LocalDomChg>::iterator& a,
         const std::set<LocalDomChg>::iterator& b) { return a->pos < b->pos; });
  std::set<LocalDomChg>::iterator elem = resolveQueue.back();
  resolveQueue.pop_back();
  return elem;
}

void HighsDomain::ConflictSet::clearQueue() { resolveQueue.clear(); }

HighsInt HighsDomain::ConflictSet::queueSize() const {
  return resolveQueue.size();
}

bool HighsDomain::ConflictSet::resolvable(HighsInt domChgPos) const {
  assert(domChgPos >= 0);
  assert(domChgPos < static_cast<HighsInt>(localdom.domchgreason_.size()));
  // printf("domchgPos: %d\n", domChgPos);
  // printf("stacksize: %ld\n", localdom.domchgreason_.size());
  switch (localdom.domchgreason_[domChgPos].type) {
    case Reason::kBranching:
    case Reason::kUnknown:
      return false;
    default:
      return true;
  }
}

HighsInt HighsDomain::ConflictSet::resolveDepth(std::set<LocalDomChg>& frontier,
                                                HighsInt depthLevel,
                                                HighsInt stopSize,
                                                HighsPseudocost& pseudocost,
                                                HighsInt minResolve,
                                                bool increaseConflictScore) {
  clearQueue();
  LocalDomChg startPos =
      LocalDomChg{depthLevel == 0 ? 0 : localdom.branchPos_[depthLevel - 1] + 1,
                  HighsDomainChange()};
  while (depthLevel < static_cast<HighsInt>(localdom.branchPos_.size())) {
    HighsInt branchPos = localdom.branchPos_[depthLevel];
    if (localdom.domchgstack_[branchPos].boundval !=
        localdom.prevboundval_[branchPos].first)
      break;
    // printf("skipping redundant depth\n");
    ++depthLevel;
  }

  auto iterEnd =
      depthLevel == static_cast<HighsInt>(localdom.branchPos_.size())
          ? frontier.end()
          : frontier.upper_bound(LocalDomChg{localdom.branchPos_[depthLevel],
                                             HighsDomainChange()});
  bool empty = true;
  for (auto it = frontier.lower_bound(startPos); it != iterEnd; ++it) {
    assert(it != frontier.end());
    empty = false;
    if (resolvable(it->pos)) pushQueue(it);
  }

  if (empty) return -1;

  HighsInt numResolved = 0;

  while (queueSize() > stopSize ||
         (queueSize() > 0 && numResolved < minResolve)) {
    std::set<LocalDomChg>::iterator pos = popQueue();
    if (!explainBoundChange(frontier, *pos)) continue;

    ++numResolved;
    frontier.erase(pos);
    for (const LocalDomChg& i : resolvedDomainChanges) {
      auto insertResult = insertIntoFrontier(frontier, i);
      if (insertResult.second) {
        if (increaseConflictScore) {
          if (localdom.domchgstack_[i.pos].boundtype == HighsBoundType::kLower)
            pseudocost.increaseConflictScoreUp(
                localdom.domchgstack_[i.pos].column);
          else
            pseudocost.increaseConflictScoreDown(
                localdom.domchgstack_[i.pos].column);
        }
        if (i.pos >= startPos.pos && resolvable(i.pos))
          pushQueue(insertResult.first);
      }
    }
  }

  return numResolved;
}

std::pair<std::set<HighsDomain::ConflictSet::LocalDomChg>::iterator, bool>
HighsDomain::ConflictSet::insertIntoFrontier(std::set<LocalDomChg>& frontier,
                                             const LocalDomChg& locdomchg) {
  // the frontier is ordered by stack position only, so if an entry for the
  // same position already exists, keep the stronger of the two relaxed bounds
  auto insertResult = frontier.insert(locdomchg);
  if (!insertResult.second) {
    double& bound = insertResult.first->domchg.boundval;
    if (locdomchg.domchg.boundtype == HighsBoundType::kLower)
      bound = std::max(bound, locdomchg.domchg.boundval);
    else
      bound = std::min(bound, locdomchg.domchg.boundval);
  }
  return insertResult;
}

HighsInt HighsDomain::ConflictSet::computeCuts(HighsInt depthLevel,
                                               HighsConflictPool& conflictPool,
                                               HighsPseudocost& pseudocost) {
  HighsInt numResolved = resolveDepth(
      reasonSideFrontier, depthLevel, 1, pseudocost,
      depthLevel == static_cast<HighsInt>(localdom.branchPos_.size()) ? 1 : 0,
      true);
  if (numResolved == -1) return -1;
  HighsInt numConflicts = 0;
  if (numResolved > 0) {
    // add conflict cut
    localdom.mipsolver->mipdata_->debugSolution.checkConflictReasonFrontier(
        reasonSideFrontier, localdom.domchgstack_);

    conflictPool.addConflictCut(localdom, reasonSideFrontier);
    ++numConflicts;
  }

  // if the queue size is 1 then we have a resolvable UIP that is not the
  // branch vertex
  if (queueSize() == 1) {
    LocalDomChg uip = *popQueue();
    clearQueue();

    // compute the UIP reconvergence cut
    reconvergenceFrontier.clear();
    reconvergenceFrontier.insert(uip);
    HighsInt numResolved =
        resolveDepth(reconvergenceFrontier, depthLevel, 0, pseudocost);

    if (numResolved > 0 && reconvergenceFrontier.count(uip) == 0) {
      localdom.mipsolver->mipdata_->debugSolution
          .checkConflictReconvergenceFrontier(reconvergenceFrontier, uip,
                                              localdom.domchgstack_);
      conflictPool.addReconvergenceCut(localdom, reconvergenceFrontier,
                                       uip.domchg);
      ++numConflicts;
    }
  }

  return numConflicts;
}

void HighsDomain::ConflictSet::conflictAnalysis(HighsConflictPool& conflictPool,
                                                HighsPseudocost& pseudocost) {
  resolvedDomainChanges.reserve(localdom.domchgstack_.size());

  if (!explainInfeasibility()) return;

  pseudocost.increaseConflictWeight();

  for (const LocalDomChg& locdomchg : resolvedDomainChanges) {
    if (locdomchg.domchg.boundtype == HighsBoundType::kLower) {
      pseudocost.increaseConflictScoreUp(locdomchg.domchg.column);
    } else {
      pseudocost.increaseConflictScoreDown(locdomchg.domchg.column);
    }
  }

  if (10 * resolvedDomainChanges.size() >
      1000 + 3 * localdom.mipsolver->mipdata_->integral_cols.size())
    return;

  for (const LocalDomChg& locdomchg : resolvedDomainChanges)
    insertIntoFrontier(reasonSideFrontier, locdomchg);

  localdom.mipsolver->mipdata_->debugSolution.checkConflictReasonFrontier(
      reasonSideFrontier, localdom.domchgstack_);

  HighsInt numConflicts = 0;
  HighsInt lastDepth = static_cast<HighsInt>(localdom.branchPos_.size());
  // printf("start conflict analysis\n");
  HighsInt currDepth;
  for (currDepth = lastDepth; currDepth >= 0; --currDepth) {
    if (currDepth > 0) {
      // skip redundant branching changes which are just added for symmetry
      // handling
      HighsInt branchpos = localdom.branchPos_[currDepth - 1];
      if (localdom.domchgstack_[branchpos].boundval ==
          localdom.prevboundval_[branchpos].first) {
        --lastDepth;
        continue;
      }
    }
    HighsInt numNewConflicts = computeCuts(currDepth, conflictPool, pseudocost);
    // if the depth level was empty, do not consider it
    if (numNewConflicts == -1) {
      --lastDepth;
      continue;
    }

    numConflicts += numNewConflicts;
    // if no conflict was found in the first non-empty depth level we stop here
    if (numConflicts == 0) break;
    // if no conflict was found in this depth level and all conflicts of the
    // first 5 non-empty depth levels are generated we stop here
    if (lastDepth - currDepth >= 4 && numNewConflicts == 0) break;
  }

  // if we stopped in the highest non-empty depth it means no conflicts where
  // added yet. We want to at least add the current conflict frontier as it
  // means the bound change leading to infeasibility was the last branching
  // itself and hence should have been propagated in the previous depth but was
  // not, e.g. because the threshold for an integral variable was not reached.
  if (currDepth == lastDepth)
    conflictPool.addConflictCut(localdom, reasonSideFrontier);
}

void HighsDomain::ConflictSet::conflictAnalysis(const HighsInt* proofinds,
                                                const double* proofvals,
                                                HighsInt prooflen,
                                                double proofrhs,
                                                HighsConflictPool& conflictPool,
                                                HighsPseudocost& pseudocost) {
  resolvedDomainChanges.reserve(localdom.domchgstack_.size());

  HighsInt ninfmin;
  HighsCDouble activitymin;
  globaldom.computeMinActivity(0, prooflen, proofinds, proofvals, ninfmin,
                               activitymin);
  if (ninfmin != 0) return;

  if (!explainInfeasibilityLeq(proofinds, proofvals, prooflen, proofrhs,
                               static_cast<double>(activitymin)))
    return;

  pseudocost.increaseConflictWeight();
  for (const LocalDomChg& locdomchg : resolvedDomainChanges) {
    if (locdomchg.domchg.boundtype == HighsBoundType::kLower)
      pseudocost.increaseConflictScoreUp(locdomchg.domchg.column);
    else
      pseudocost.increaseConflictScoreDown(locdomchg.domchg.column);
  }

  if (10 * resolvedDomainChanges.size() >
      1000 + 3 * localdom.mipsolver->mipdata_->integral_cols.size())
    return;

  reasonSideFrontier.insert(resolvedDomainChanges.begin(),
                            resolvedDomainChanges.end());

  assert(resolvedDomainChanges.size() == reasonSideFrontier.size());

  localdom.mipsolver->mipdata_->debugSolution.checkConflictReasonFrontier(
      reasonSideFrontier, localdom.domchgstack_);

  HighsInt numConflicts = 0;
  HighsInt lastDepth = static_cast<HighsInt>(localdom.branchPos_.size());
  HighsInt currDepth;
  for (currDepth = lastDepth; currDepth >= 0; --currDepth) {
    if (currDepth > 0) {
      // skip redundant branching changes which are just added for symmetry
      // handling
      HighsInt branchpos = localdom.branchPos_[currDepth - 1];
      if (localdom.domchgstack_[branchpos].boundval ==
          localdom.prevboundval_[branchpos].first) {
        --lastDepth;
        continue;
      }
    }
    HighsInt numNewConflicts = computeCuts(currDepth, conflictPool, pseudocost);
    // if the depth level was empty, do not consider it
    if (numNewConflicts == -1) {
      --lastDepth;
      continue;
    }

    numConflicts += numNewConflicts;
    // if no conflict was found in the first non-empty depth level we stop here
    if (numConflicts == 0) break;
    // if no conflict was found in this depth level and all conflicts of the
    // first 5 non-empty depth levels are generated we stop here
    if (lastDepth - currDepth >= 4 && numNewConflicts == 0) break;
  }

  // if we stopped in the highest non-empty depth it means no conflicts where
  // added yet. We want to at least add the current conflict frontier as it
  // means the bound change leading to infeasibility was the last branching
  // itself and hence should have been propagated in the previous depth but was
  // not, e.g. because the threshold for an integral variable was not reached.
  if (currDepth == lastDepth)
    conflictPool.addConflictCut(localdom, reasonSideFrontier);
}
