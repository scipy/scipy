/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_CONFLICTPOOL_H_
#define HIGHS_CONFLICTPOOL_H_

#include <atomic>
#include <set>
#include <vector>

#include "mip/HighsDomain.h"
#include "util/HighsType.h"

class HighsConflictPool {
 private:
  HighsInt agelim_;
  HighsInt softlimit_;
  bool age_lock_;
  std::vector<HighsInt> ageDistribution_;
  std::vector<int16_t> ages_;
  std::vector<unsigned> modification_;
  std::deque<std::atomic<uint8_t>> ageResetWhileLocked_;

  std::vector<HighsDomainChange> conflictEntries_;
  std::vector<std::pair<HighsInt, HighsInt>> conflictRanges_;

  /// keep an ordered set of free spaces in the row arrays so that they can be
  /// reused efficiently
  std::set<std::pair<HighsInt, HighsInt>> freeSpaces_;

  /// vector of deleted conflicts so that their indices can be reused
  std::vector<HighsInt> deletedConflicts_;

  std::vector<HighsDomain::ConflictPoolPropagation*> propagationDomains;

  std::pair<HighsInt, HighsInt> allocateConflict(HighsInt conflictLen);

 public:
  HighsConflictPool(HighsInt agelim, HighsInt softlimit)
      : agelim_(agelim),
        softlimit_(softlimit),
        age_lock_(false),
        ageDistribution_(),
        ages_(),
        modification_(),
        ageResetWhileLocked_(),
        conflictEntries_(),
        conflictRanges_(),
        freeSpaces_(),
        deletedConflicts_(),
        propagationDomains() {
    ageDistribution_.resize(agelim_ + 1);
  }

  void addConflictCut(const HighsDomain& domain,
                      const std::set<HighsDomain::ConflictSet::LocalDomChg>&
                          reasonSideFrontier);

  void addReconvergenceCut(
      const HighsDomain& domain,
      const std::set<HighsDomain::ConflictSet::LocalDomChg>&
          reconvergenceFrontier,
      const HighsDomainChange& reconvergenceDomchg);

  void removeConflict(HighsInt conflict);

  void performAging(bool thread_safe = false);

  void addConflictFromOtherPool(const HighsDomainChange* conflictEntries,
                                HighsInt conflictLen);

  void syncConflictPool(HighsConflictPool& syncpool);

  void resetAge(HighsInt conflict) {
    if (ages_[conflict] > 0) {
      if (age_lock_) {
        ageResetWhileLocked_[conflict].store(1, std::memory_order_relaxed);
        return;
      }
      ageDistribution_[ages_[conflict]] -= 1;
      ageDistribution_[0] += 1;
      ages_[conflict] = 0;
    }
  }

  void setAgeLimit(HighsInt agelim) {
    agelim_ = agelim;
    ageDistribution_.resize(agelim_ + 1);
  }

  unsigned getModificationCount(HighsInt cut) const {
    return modification_[cut];
  }

  void addPropagationDomain(HighsDomain::ConflictPoolPropagation* domain) {
    propagationDomains.push_back(domain);
  }

  void removePropagationDomain(HighsDomain::ConflictPoolPropagation* domain) {
    for (HighsInt k = propagationDomains.size() - 1; k >= 0; --k) {
      if (propagationDomains[k] == domain) {
        propagationDomains.erase(propagationDomains.begin() + k);
        return;
      }
    }
  }

  const std::vector<HighsDomainChange>& getConflictEntryVector() const {
    return conflictEntries_;
  }

  const std::vector<std::pair<HighsInt, HighsInt>>& getConflictRanges() const {
    return conflictRanges_;
  }

  HighsInt getNumConflicts() const {
    return conflictRanges_.size() - deletedConflicts_.size();
  }

  void setAgeLock(const bool ageLock) { age_lock_ = ageLock; }
};

#endif
