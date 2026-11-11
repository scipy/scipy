/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_DOMAIN_H_
#define HIGHS_DOMAIN_H_

#include <cstdint>
#include <deque>
#include <functional>
#include <memory>
#include <set>
#include <vector>

#include "HighsPseudocost.h"
#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolver.h"
#include "pdqsort/pdqsort.h"
#include "util/HighsCDouble.h"
#include "util/HighsRbTree.h"

class HighsCutPool;
class HighsConflictPool;
class HighsObjectiveFunction;

class HighsDomain {
 public:
  struct Reason {
    HighsInt type;
    HighsInt index;

    enum {
      kBranching = -1,
      kUnknown = -2,
      kModelRowUpper = -3,
      kModelRowLower = -4,
      kCliqueTable = -5,
      kConflictingBounds = -6,
      kObjective = -7,
    };
    static Reason branching() { return Reason{kBranching, 0}; }
    static Reason unspecified() { return Reason{kUnknown, 0}; }
    static Reason cliqueTable(HighsInt col, HighsInt val) {
      return Reason{kCliqueTable, 2 * col + val};
    }
    static Reason modelRowUpper(HighsInt row) {
      return Reason{kModelRowUpper, row};
    }
    static Reason modelRowLower(HighsInt row) {
      return Reason{kModelRowLower, row};
    }
    static Reason cut(HighsInt cutpool, HighsInt cut) {
      return Reason{cutpool, cut};
    }
    static Reason conflictingBounds(HighsInt pos) {
      return Reason{kConflictingBounds, pos};
    }
    static Reason objective() { return Reason{kObjective, 0}; }
  };

  class ConflictSet {
    friend class HighsDomain;
    HighsDomain& localdom;
    const HighsDomain& globaldom;

   public:
    struct LocalDomChg {
      HighsInt pos;
      mutable HighsDomainChange domchg;

      bool operator<(const LocalDomChg& other) const { return pos < other.pos; }
    };

    ConflictSet(HighsDomain& localdom, const HighsDomain& globaldom);

    void conflictAnalysis(HighsConflictPool& conflictPool,
                          HighsPseudocost& pseudocost);
    void conflictAnalysis(const HighsInt* proofinds, const double* proofvals,
                          HighsInt prooflen, double proofrhs,
                          HighsConflictPool& conflictPool,
                          HighsPseudocost& pseudocost);

   private:
    std::set<LocalDomChg> reasonSideFrontier;
    std::set<LocalDomChg> reconvergenceFrontier;
    std::vector<std::set<LocalDomChg>::iterator> resolveQueue;
    std::vector<LocalDomChg> resolvedDomainChanges;

    struct ResolveCandidate {
      double delta;
      double baseBound;
      double prio;
      HighsInt boundPos;
      HighsInt valuePos;

      bool operator<(const ResolveCandidate& other) const {
        if (prio > other.prio) return true;
        if (other.prio > prio) return false;

        return boundPos < other.boundPos;
      }
    };

    std::vector<ResolveCandidate> resolveBuffer;

    void pushQueue(std::set<LocalDomChg>::iterator domchgPos);
    std::set<LocalDomChg>::iterator popQueue();
    void clearQueue();
    HighsInt queueSize() const;
    bool resolvable(HighsInt domChgPos) const;

    static std::pair<std::set<LocalDomChg>::iterator, bool> insertIntoFrontier(
        std::set<LocalDomChg>& frontier, const LocalDomChg& locdomchg);

    HighsInt resolveDepth(std::set<LocalDomChg>& frontier, HighsInt depthLevel,
                          HighsInt stopSize, HighsPseudocost& pseudocost,
                          HighsInt minResolve = 0,
                          bool increaseConflictScore = false);

    HighsInt computeCuts(HighsInt depthLevel, HighsConflictPool& conflictPool,
                         HighsPseudocost& pseudocost);

    bool explainInfeasibility();

    bool explainInfeasibilityConflict(const HighsDomainChange* conflict,
                                      HighsInt len);

    bool explainInfeasibilityLeq(const HighsInt* inds, const double* vals,
                                 HighsInt len, double rhs, double minActivity);

    bool explainInfeasibilityGeq(const HighsInt* inds, const double* vals,
                                 HighsInt len, double rhs, double maxActivity);

    bool explainBoundChange(const std::set<LocalDomChg>& currentFrontier,
                            LocalDomChg domchg);

    // bool explainBoundChange(HighsInt pos) {
    //   return explainBoundChange(LocalDomChg{pos,
    //   localdom.domchgstack_[pos]});
    // }

    bool explainBoundChangeConflict(const LocalDomChg& domchg,
                                    const HighsDomainChange* conflict,
                                    HighsInt len);

    bool explainBoundChangeLeq(const std::set<LocalDomChg>& currentFrontier,
                               const LocalDomChg& domChg, const HighsInt* inds,
                               const double* vals, HighsInt len, double rhs,
                               double minActivity);

    bool explainBoundChangeGeq(const std::set<LocalDomChg>& currentFrontier,
                               const LocalDomChg& domChg, const HighsInt* inds,
                               const double* vals, HighsInt len, double rhs,
                               double maxActivity);

    bool resolveLinearLeq(HighsCDouble M, double Mlower, const double* vals);

    bool resolveLinearGeq(HighsCDouble M, double Mupper, const double* vals);
  };

  struct CutpoolPropagation {
    HighsInt cutpoolindex;
    HighsDomain* domain;
    HighsCutPool* cutpool;
    std::vector<HighsCDouble> activitycuts_;
    std::vector<HighsInt> activitycutsinf_;
    std::vector<uint8_t> propagatecutflags_;
    std::vector<HighsInt> propagatecutinds_;
    std::vector<double> capacityThreshold_;

    CutpoolPropagation(HighsInt cutpoolindex, HighsDomain* domain,
                       HighsCutPool& cutpool);

    CutpoolPropagation(const CutpoolPropagation& other);

    CutpoolPropagation& operator=(const CutpoolPropagation& other);

    ~CutpoolPropagation();

    void recomputeCapacityThreshold(HighsInt cut);

    void cutAdded(HighsInt cut, bool propagate);

    void cutDeleted(HighsInt cut, bool deletedOnlyForPropagation = false);

    void markPropagateCut(HighsInt cut);

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound,
                                bool threshold, bool activity,
                                bool infeasdomain);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound,
                                bool threshold, bool activity,
                                bool infeasdomain);
  };

  struct ConflictPoolPropagation {
    HighsInt conflictpoolindex;
    HighsDomain* domain;
    HighsConflictPool* conflictpool_;
    std::vector<HighsInt> colLowerWatched_;
    std::vector<HighsInt> colUpperWatched_;
    std::vector<uint8_t> conflictFlag_;
    std::vector<HighsInt> propagateConflictInds_;

    struct WatchedLiteral {
      HighsDomainChange domchg = {0.0, -1, HighsBoundType::kLower};
      HighsInt prev = -1;
      HighsInt next = -1;
    };

    std::vector<WatchedLiteral> watchedLiterals_;

    ConflictPoolPropagation(HighsInt conflictpoolindex, HighsDomain* domain,
                            HighsConflictPool& cutpool);

    ConflictPoolPropagation(const ConflictPoolPropagation& other);

    ConflictPoolPropagation& operator=(const ConflictPoolPropagation& other);

    ~ConflictPoolPropagation();

    void linkWatchedLiteral(HighsInt linkPos);

    void unlinkWatchedLiteral(HighsInt linkPos);

    void conflictAdded(HighsInt conflict);

    void conflictDeleted(HighsInt conflict);

    void markPropagateConflict(HighsInt conflict);

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

    void propagateConflict(HighsInt conflict);
  };

  struct DualFixProbingPropagation {
    HighsDomain* domain = nullptr;
    HighsMipSolver* mipsolver = nullptr;

    struct RowSide {
      HighsInt row;
      bool isRhs;
      operator size_t() const { return 2 * row + isRhs; }
    };

    std::vector<HighsBool> redundantRowFlags_;
    std::vector<RowSide> redundantRowInds_;
    HighsInt previousRedundantRowSize = 0;
    HighsInt numGlobalRedundantRows_ = 0;
    HighsInt numCachedGlobalRedundantRows_ = 0;

    // Track direction of zero fixings so we don't store disagreeing results
    enum DualFixProbingFixDirection {
      FixUndecided,
      FixLowerBound,
      FixUpperBound,
    };

    struct FixedZeroCostColumn {
      HighsInt col;
      DualFixProbingFixDirection direction;
    };

    std::vector<DualFixProbingFixDirection> zeroCostDirections_;
    std::vector<FixedZeroCostColumn> fixedZeroCostColumns_;

    bool applyingZeroCostFixings_ = false;
    HighsInt zeroCostStartPos_ = kHighsIInf;
    std::function<void()> storeLiftingOpportunity;

    bool enabled_ = false;

    // Remaining global lower / upper locks
    std::vector<HighsInt> colLowerLocksOriginal_;
    std::vector<HighsInt> colUpperLocksOriginal_;
    // Number of removed locks after propagation.
    std::vector<HighsInt> colLowerReducedNumLocks_;
    std::vector<HighsInt> colUpperReducedNumLocks_;
    std::vector<HighsInt> clearColNumReducedLocks_;

    std::vector<HighsInt> candidateFixedCols_;
    std::vector<HighsBool> candidateColFixedFlags_;
    std::vector<HighsInt> globalCandidateFixedCols_;

    void setEnabled(const bool val) { enabled_ = val; }

    bool isEnabled() const { return enabled_; }

    // active only when new redundant rows are found.
    bool isActive() const {
      return enabled_ && (!candidateFixedCols_.empty() ||
                          static_cast<HighsInt>(redundantRowInds_.size()) >
                              previousRedundantRowSize);
    }

    bool isZeroCostFixingActive() const { return applyingZeroCostFixings_; }

    // mark the position when the first zero-cost variable can be fixed to its
    // lower or upper bound.
    void setZeroCostFixingPosition(HighsInt v) { zeroCostStartPos_ = v; }

    HighsInt getZeroCostFixingPosition() const { return zeroCostStartPos_; }

    bool ableToFixToLb(const HighsInt col) const {
      return domain->col_lower_[col] != -kHighsInf &&
             mipsolver->model_->col_cost_[col] >= 0;
    }

    bool ableToFixToUb(const HighsInt col) const {
      return domain->col_upper_[col] != kHighsInf &&
             mipsolver->model_->col_cost_[col] <= 0;
    }

    void beginProbing() {
      previousRedundantRowSize = numGlobalRedundantRows_;
      assert(static_cast<HighsInt>(redundantRowInds_.size()) ==
             numGlobalRedundantRows_);
      fixedZeroCostColumns_.clear();
      setZeroCostFixingPosition(kHighsIInf);
      applyingZeroCostFixings_ = false;
      setEnabled(true);

      for (const auto col : clearColNumReducedLocks_) {
        colLowerReducedNumLocks_[col] = 0;
        colUpperReducedNumLocks_[col] = 0;
      }
      clearColNumReducedLocks_.clear();
      cacheGlobalRedundantRows();
    }

    void endProbing() {
      setEnabled(false);
      for (HighsInt i = numGlobalRedundantRows_;
           i < static_cast<HighsInt>(redundantRowInds_.size()); ++i) {
        redundantRowFlags_[redundantRowInds_[i]] = false;
      }
      redundantRowInds_.resize(numGlobalRedundantRows_);
      fixedZeroCostColumns_.clear();
      applyingZeroCostFixings_ = false;
      storeLiftingOpportunity = nullptr;
      for (const HighsInt col : candidateFixedCols_) {
        candidateColFixedFlags_[col] = false;
      }
      candidateFixedCols_.clear();
    }

    DualFixProbingPropagation() = default;

    DualFixProbingPropagation(const DualFixProbingPropagation& other) = delete;

    DualFixProbingPropagation& operator=(
        const DualFixProbingPropagation& other) = delete;

    void recomputeLocks();
    void cacheGlobalRedundantRows();
    void updateRhsRedundant(HighsInt row);
    void updateLhsRedundant(HighsInt row);
    void propagate();
    void propagateZeroCosts();
  };

 private:
  struct ObjectivePropagation {
    HighsDomain* domain = nullptr;
    const HighsObjectiveFunction* objFunc;
    const double* cost;
    HighsCDouble objectiveLower;
    HighsInt numInfObjLower;
    double capacityThreshold;
    bool isPropagated;

    struct ObjectiveContribution {
      double contribution;
      HighsInt col;
      HighsInt partition;
      highs::RbTreeLinks<HighsInt> links;
    };

    class ObjectiveContributionTree;

    std::vector<ObjectiveContribution> objectiveLowerContributions;
    std::vector<std::pair<HighsInt, HighsInt>> contributionPartitionSets;
    std::vector<double> propagationConsBuffer;
    struct PartitionCliqueData {
      double multiplier;
      HighsInt rhs;
      bool changed;
    };

    std::vector<PartitionCliqueData> partitionCliqueData;

    ObjectivePropagation() {
      objFunc = nullptr;
      cost = nullptr;
      objectiveLower = 0.0;
      numInfObjLower = 0;
      capacityThreshold = 0.0;
      isPropagated = false;
    }
    ObjectivePropagation(HighsDomain* domain);

    bool isActive() const { return domain != nullptr; }

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

    bool shouldBePropagated() const;

    void propagate();

    void debugCheckObjectiveLower() const;

    // construct the proot constraint at the time when the domain change stack
    // had the given size
    void getPropagationConstraint(HighsInt domchgStackSize, const double*& vals,
                                  const HighsInt*& inds, HighsInt& len,
                                  double& rhs, HighsInt domchgCol = -1);

   private:
    void recomputeCapacityThreshold();
  };

  std::vector<HighsBool> changedcolsflags_;
  std::vector<HighsInt> changedcols_;

  std::vector<std::pair<HighsInt, HighsInt>> propRowNumChangedBounds_;

  std::vector<HighsDomainChange> domchgstack_;
  std::vector<Reason> domchgreason_;
  std::vector<std::pair<double, HighsInt>> prevboundval_;

  std::vector<HighsCDouble> activitymin_;
  std::vector<HighsCDouble> activitymax_;
  std::vector<HighsInt> activitymininf_;
  std::vector<HighsInt> activitymaxinf_;
  std::vector<double> capacityThreshold_;
  std::vector<HighsBool> propagateflags_;
  std::vector<HighsInt> propagateinds_;
  ObjectivePropagation objProp_;

  HighsMipSolver* mipsolver;

 private:
  std::deque<CutpoolPropagation> cutpoolpropagation;
  std::deque<ConflictPoolPropagation> conflictPoolPropagation;
  DualFixProbingPropagation dualFixProbingPropagation;

  bool infeasible_ = false;
  Reason infeasible_reason;
  HighsInt infeasible_pos;

  void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

  void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

  void updateThresholdLbChange(HighsInt col, double newbound, double val,
                               double& threshold) const;

  void updateThresholdUbChange(HighsInt col, double newbound, double val,
                               double& threshold) const;

  void recomputeCapacityThreshold(HighsInt row);

  void updateRedundantRows(HighsInt row);

  double doChangeBound(const HighsDomainChange& boundchg);

  std::vector<HighsInt> colLowerPos_;
  std::vector<HighsInt> colUpperPos_;
  std::vector<HighsInt> branchPos_;
  HighsHashTable<HighsInt> redundantRows_;
  bool recordRedundantRows_ = false;
  bool dualFixProbingActive_ = false;

 public:
  std::vector<double> col_lower_;
  std::vector<double> col_upper_;

  HighsDomain(HighsMipSolver& mipsolver);

  HighsDomain(const HighsDomain& other)
      : changedcolsflags_(other.changedcolsflags_),
        changedcols_(other.changedcols_),
        domchgstack_(other.domchgstack_),
        domchgreason_(other.domchgreason_),
        prevboundval_(other.prevboundval_),
        activitymin_(other.activitymin_),
        activitymax_(other.activitymax_),
        activitymininf_(other.activitymininf_),
        activitymaxinf_(other.activitymaxinf_),
        capacityThreshold_(other.capacityThreshold_),
        propagateflags_(other.propagateflags_),
        propagateinds_(other.propagateinds_),
        objProp_(other.objProp_),
        mipsolver(other.mipsolver),
        cutpoolpropagation(other.cutpoolpropagation),
        conflictPoolPropagation(other.conflictPoolPropagation),
        infeasible_(other.infeasible_),
        infeasible_reason(other.infeasible_reason),
        infeasible_pos(other.infeasible_pos),
        colLowerPos_(other.colLowerPos_),
        colUpperPos_(other.colUpperPos_),
        branchPos_(other.branchPos_),
        col_lower_(other.col_lower_),
        col_upper_(other.col_upper_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.domain = this;
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.domain = this;
    if (objProp_.domain) objProp_.domain = this;
    dualFixProbingPropagation.domain = this;
    dualFixProbingPropagation.mipsolver = mipsolver;
  }

  HighsDomain& operator=(const HighsDomain& other) {
    if (this == &other) return *this;
    changedcolsflags_ = other.changedcolsflags_;
    changedcols_ = other.changedcols_;
    domchgstack_ = other.domchgstack_;
    domchgreason_ = other.domchgreason_;
    prevboundval_ = other.prevboundval_;
    activitymin_ = other.activitymin_;
    activitymax_ = other.activitymax_;
    activitymininf_ = other.activitymininf_;
    activitymaxinf_ = other.activitymaxinf_;
    capacityThreshold_ = other.capacityThreshold_;
    propagateflags_ = other.propagateflags_;
    propagateinds_ = other.propagateinds_;
    objProp_ = other.objProp_;
    mipsolver = other.mipsolver;
    cutpoolpropagation = other.cutpoolpropagation;
    conflictPoolPropagation = other.conflictPoolPropagation;
    infeasible_ = other.infeasible_;
    infeasible_reason = other.infeasible_reason;
    colLowerPos_ = other.colLowerPos_;
    colUpperPos_ = other.colUpperPos_;
    branchPos_ = other.branchPos_;
    col_lower_ = other.col_lower_;
    col_upper_ = other.col_upper_;
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.domain = this;
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.domain = this;
    if (objProp_.domain) objProp_.domain = this;
    dualFixProbingPropagation.domain = this;
    dualFixProbingActive_ = false;
    dualFixProbingPropagation.mipsolver = mipsolver;
    dualFixProbingPropagation.endProbing();
    return *this;
  }

  void computeMinActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const double* ARvalue, HighsInt& ninfmin,
                          HighsCDouble& activitymin) const;

  void computeMaxActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const double* ARvalue, HighsInt& ninfmax,
                          HighsCDouble& activitymax) const;

  double adjustedUb(HighsInt col, HighsCDouble boundVal, bool& accept) const;

  double adjustedLb(HighsInt col, HighsCDouble boundVal, bool& accept) const;

  HighsInt propagateRowUpper(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rupper,
                             const HighsCDouble& minactivity, HighsInt ninfmin,
                             HighsDomainChange* boundchgs) const;

  HighsInt propagateRowLower(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rlower,
                             const HighsCDouble& maxactivity, HighsInt ninfmax,
                             HighsDomainChange* boundchgs) const;

  const std::vector<HighsInt>& getChangedCols() const { return changedcols_; }

  void addCutpool(HighsCutPool& cutpool);

  void addConflictPool(HighsConflictPool& conflictPool);

  void clearChangedCols() {
    for (HighsInt i : changedcols_) changedcolsflags_[i] = false;
    changedcols_.clear();
  }

  void removeContinuousChangedCols() {
    for (HighsInt i : changedcols_)
      changedcolsflags_[i] = mipsolver->isColIntegral(i);

    changedcols_.erase(
        std::remove_if(changedcols_.begin(), changedcols_.end(),
                       [&](HighsInt i) { return !isChangedCol(i); }),
        changedcols_.end());
  }

  void clearChangedCols(size_t start) {
    for (size_t i = start; i != changedcols_.size(); ++i)
      changedcolsflags_[changedcols_[i]] = false;

    changedcols_.resize(start);
  }

  bool isChangedCol(HighsInt col) const { return changedcolsflags_[col]; }

  void markPropagate(HighsInt row);

  bool isActive(const HighsDomainChange& domchg) const {
    return domchg.boundtype == HighsBoundType::kLower
               ? domchg.boundval <= col_lower_[domchg.column]
               : domchg.boundval >= col_upper_[domchg.column];
  }

  void markPropagateCut(Reason reason);

  void setupObjectivePropagation() { objProp_ = ObjectivePropagation(this); }

  void computeRowActivities();

  void markInfeasible(Reason reason = Reason::unspecified()) {
    infeasible_ = true;
    infeasible_pos = domchgstack_.size();
    infeasible_reason = reason;
  }

  bool infeasible() const { return infeasible_; }

  void changeBound(HighsDomainChange boundchg,
                   Reason reason = Reason::branching());

  void changeBound(HighsBoundType boundtype, HighsInt col, double boundval,
                   Reason reason = Reason::branching()) {
    changeBound({boundval, col, boundtype}, reason);
  }

  bool checkChangeBound(HighsBoundType boundtype, HighsInt col,
                        HighsCDouble boundval, Reason reason);

  void fixCol(HighsInt col, double val, Reason reason = Reason::unspecified()) {
    if (kAllowDeveloperAssert) {
      assert(infeasible_ == 0);
    }
    if (col_lower_[col] < val) {
      changeBound({val, col, HighsBoundType::kLower}, reason);
      if (infeasible_ == 0) propagate();
    }

    if (infeasible_ == 0 && col_upper_[col] > val)
      changeBound({val, col, HighsBoundType::kUpper}, reason);
  }

  void backtrackToGlobal();

  HighsDomainChange backtrack();

  const std::vector<HighsInt>& getBranchingPositions() const {
    return branchPos_;
  }

  const std::vector<std::pair<double, HighsInt>>& getPreviousBounds() const {
    return prevboundval_;
  }

  // Walk the prevboundval_ chain until we find a position where the bound
  // value is strictly past the target (in the given direction).
  // direction = +1: walk while value >= target (lower bound relaxation)
  // direction = -1: walk while value <= target (upper bound relaxation)
  // Returns false if the chain is exhausted (pos becomes -1).
  bool walkBoundChain(HighsInt& pos, double target, HighsInt direction) const {
    while (pos != -1 &&
           direction * prevboundval_[pos].first >= direction * target) {
      pos = prevboundval_[pos].second;
    }
    return pos != -1;
  }

  const std::vector<HighsDomainChange>& getDomainChangeStack() const {
    return domchgstack_;
  }

  const std::vector<Reason>& getDomainChangeReason() const {
    return domchgreason_;
  }

  double getObjectiveLowerBound() const {
    if (objProp_.isActive() && objProp_.numInfObjLower == 0)
      return static_cast<double>(objProp_.objectiveLower);

    return -kHighsInf;
  }

  void getCutoffConstraint(const double*& vals, const HighsInt*& inds,
                           HighsInt& len, double& rhs) {
    objProp_.getPropagationConstraint(domchgstack_.size(), vals, inds, len,
                                      rhs);
  }

  HighsInt getNumDomainChanges() const { return domchgstack_.size(); }

  bool colBoundsAreGlobal(HighsInt col) const {
    return colLowerPos_[col] == -1 && colUpperPos_[col] == -1;
  }

  HighsInt getBranchDepth() const { return branchPos_.size(); }

  std::vector<HighsDomainChange> getReducedDomainChangeStack(
      std::vector<HighsInt>& branchingPositions) const {
    std::vector<HighsDomainChange> reducedstack;
    if (domchgstack_.empty()) return reducedstack;
    auto addChangeToReducedStack = [&](const HighsInt i) {
      if (domchgreason_[i].type == Reason::kBranching)
        branchingPositions.push_back(reducedstack.size());
      else {
        HighsInt k = i;
        while (prevboundval_[k].second != -1) {
          k = prevboundval_[k].second;
          if (domchgreason_[k].type == Reason::kBranching) {
            branchingPositions.push_back(reducedstack.size());
            break;
          }
        }
      }

      reducedstack.push_back(domchgstack_[i]);
    };
    // Collect tightest bound on each column in stack order
    const HighsInt stackSize = static_cast<HighsInt>(domchgstack_.size());
    const HighsInt numCols = mipsolver->numCol();
    HighsInt logNumCols = 0;
    for (HighsInt v = numCols - 1; v > 0; v >>= 1) {
      ++logNumCols;
    }
    if ((stackSize - 1) / numCols <= logNumCols) {
      // Stack smaller than n * (1 + ceil(log2(n)))
      reducedstack.reserve(stackSize);
      branchingPositions.reserve(branchPos_.size());
      for (HighsInt i = 0; i < stackSize; i++) {
        // keep only the tightest bound change for each variable
        if ((domchgstack_[i].boundtype == HighsBoundType::kLower &&
             colLowerPos_[domchgstack_[i].column] != i) ||
            (domchgstack_[i].boundtype == HighsBoundType::kUpper &&
             colUpperPos_[domchgstack_[i].column] != i))
          continue;
        addChangeToReducedStack(i);
      }
    } else {
      // Stack is large. Avoid iterating over entire stack to
      // hedge against worst-case performance
      std::vector<HighsInt> positions;
      positions.reserve(2 * numCols);
      for (HighsInt col = 0; col < numCols; col++) {
        if (colLowerPos_[col] != -1) positions.push_back(colLowerPos_[col]);
        if (colUpperPos_[col] != -1) positions.push_back(colUpperPos_[col]);
      }
      pdqsort(positions.begin(), positions.end());
      reducedstack.reserve(positions.size());
      branchingPositions.reserve(std::min(positions.size(), branchPos_.size()));
      for (const HighsInt i : positions) {
        addChangeToReducedStack(i);
      }
    }

    reducedstack.shrink_to_fit();
    return reducedstack;
  }

  void setDomainChangeStack(const std::vector<HighsDomainChange>& domchgstack);

  void setDomainChangeStack(const std::vector<HighsDomainChange>& domchgstack,
                            const std::vector<HighsInt>& branchingPositions);

  bool propagate();

  double getColLowerPos(HighsInt col, HighsInt stackpos, HighsInt& pos) const;

  double getColUpperPos(HighsInt col, HighsInt stackpos, HighsInt& pos) const;

  void conflictAnalysis(HighsConflictPool& conflictPool, HighsDomain& globaldom,
                        HighsPseudocost& pseudocost);

  void conflictAnalysis(const HighsInt* proofinds, const double* proofvals,
                        HighsInt prooflen, double proofrhs,
                        HighsConflictPool& conflictPool, HighsDomain& globaldom,
                        HighsPseudocost& pseudocost);

  void conflictAnalyzeReconvergence(const HighsDomainChange& domchg,
                                    const HighsInt* proofinds,
                                    const double* proofvals, HighsInt prooflen,
                                    double proofrhs,
                                    HighsConflictPool& conflictPool,
                                    HighsDomain& globaldom,
                                    HighsPseudocost& pseudocost);

  void tightenCoefficients(HighsInt* inds, double* vals, HighsInt len,
                           double& rhs) const;

  double getMinActivity(HighsInt row) const {
    return activitymininf_[row] == 0 ? static_cast<double>(activitymin_[row])
                                     : -kHighsInf;
  }

  double getMaxActivity(HighsInt row) const {
    return activitymaxinf_[row] == 0 ? static_cast<double>(activitymax_[row])
                                     : kHighsInf;
  }

  // the differences between activities and row sides are computed in full
  // precision before comparing with the tolerance, since comparing a
  // HighsCDouble with a double rounds it first; this keeps the redundancy and
  // infeasibility checks consistent at the tolerance boundary
  bool isRowLowerRedundant(HighsInt row) const {
    return mipsolver->rowLower(row) == -kHighsInf ||
           (activitymininf_[row] == 0 &&
            mipsolver->rowLower(row) - activitymin_[row] <= feastol());
  }

  bool isRowUpperRedundant(HighsInt row) const {
    return mipsolver->rowUpper(row) == kHighsInf ||
           (activitymaxinf_[row] == 0 &&
            activitymax_[row] - mipsolver->rowUpper(row) <= feastol());
  }

  bool isRowLowerInfeasible(HighsInt row) const {
    return mipsolver->rowLower(row) != -kHighsInf &&
           activitymaxinf_[row] == 0 &&
           mipsolver->rowLower(row) - activitymax_[row] > feastol();
  }

  bool isRowUpperInfeasible(HighsInt row) const {
    return mipsolver->rowUpper(row) != kHighsInf && activitymininf_[row] == 0 &&
           activitymin_[row] - mipsolver->rowUpper(row) > feastol();
  }

  double getMinCutActivity(const HighsCutPool& cutpool, HighsInt cut) const;

  bool isBinary(HighsInt col) const {
    return mipsolver->isColIntegral(col) && col_lower_[col] == 0.0 &&
           col_upper_[col] == 1.0;
  }

  bool isGlobalBinary(HighsInt col) const {
    return mipsolver->isColIntegral(col) &&
           mipsolver->model_->col_lower_[col] == 0.0 &&
           mipsolver->model_->col_upper_[col] == 1.0;
  }

  HighsVarType variableType(HighsInt col) const {
    return mipsolver->variableType(col);
  }

  bool isFixed(HighsInt col) const {
    return col_lower_[col] == col_upper_[col];
  }

  bool isFixedToVal(HighsInt col, HighsInt val) const {
    return col_lower_[col] == val && col_upper_[col] == val;
  }

  bool isFixing(const HighsDomainChange& domchg) const;

  HighsDomainChange flip(const HighsDomainChange& domchg) const;

  double feastol() const;

  double epsilon() const;

  HighsInt numModelNonzeros() const { return mipsolver->numNonzero(); }

  bool inSubmip() const { return mipsolver->submip; }

  void clearRedundantRows() { redundantRows_.clear(); };

  const HighsHashTable<HighsInt>& getRedundantRows() const {
    return redundantRows_;
  };

  double getRedundantRowValue(HighsInt row) const;

  void setRecordRedundantRows(bool val) { recordRedundantRows_ = val; };

  bool isRedundantRow(HighsInt row) const;

  DualFixProbingPropagation& getDualFixProbingPropagation() {
    return dualFixProbingPropagation;
  }

  void setDualFixProbingActive(const bool val) { dualFixProbingActive_ = val; }

  bool getDualFixProbingActive() const { return dualFixProbingActive_; }
};

#endif
