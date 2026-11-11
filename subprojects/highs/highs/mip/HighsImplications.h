/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_IMPLICATIONS_H_
#define HIGHS_IMPLICATIONS_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "mip/HighsDomain.h"
#include "mip/HighsDomainChange.h"
#include "util/HighsHashTree.h"

class HighsCliqueTable;
class HighsLpRelaxation;

class HighsImplications {
  HighsInt nextCleanupCall;

  struct Implication {
    double lb = -kHighsInf;
    double ub = kHighsInf;
  };

  std::vector<HighsDomainChange> implicationsDown;
  std::vector<HighsDomainChange> implicationsUp;
  std::vector<HighsHashTree<HighsInt, Implication>> implications;
  std::vector<HighsHashTree<HighsInt, bool>> reverseImplications;
  std::vector<uint8_t> hasProbed;
  int64_t numImplications;
  int64_t numVarBounds;
  int64_t maxVarBounds;

  struct ImplIdx {
    HighsInt col;
    HighsInt val;
    operator size_t() const {
      assert(col >= 0);
      assert(val == 0 || val == 1);
      return 2 * col + val;
    }
  };

  bool computeImplications(HighsInt col, bool val);

 public:
  struct VarBound {
    double coef;
    double constant;
    HighsInt origin;

    double minValue() const {
      double m = std::min(coef, 0.0);
      if (std::abs(constant) >= kHighsInf || std::abs(m) >= kHighsInf)
        return constant + m;
      return static_cast<double>(static_cast<HighsCDouble>(constant) + m);
    }
    double maxValue() const {
      double m = std::max(coef, 0.0);
      if (std::abs(constant) >= kHighsInf || std::abs(m) >= kHighsInf)
        return constant + m;
      return static_cast<double>(static_cast<HighsCDouble>(constant) + m);
    }
  };

 private:
  std::vector<HighsHashTree<HighsInt, VarBound>> vubs;
  std::vector<HighsHashTree<HighsInt, VarBound>> vlbs;
  struct TentativeFixing {
    enum Direction : uint8_t { Undecided, FixLower, FixUpper };
    Direction downProbe = Undecided;
    Direction upProbe = Undecided;

    bool isUndecided() const {
      return downProbe == Undecided && upProbe == Undecided;
    }

    void record(bool upProbing, HighsBoundType boundtype) {
      Direction& probe = upProbing ? upProbe : downProbe;
      if (probe != Undecided) return;
      probe = boundtype == HighsBoundType::kLower ? FixUpper : FixLower;
    }

    void clear() {
      downProbe = Undecided;
      upProbe = Undecided;
    }
  };
  std::vector<HighsInt> dualFixProbingBinInds_;
  std::vector<TentativeFixing> dualFixProbingBinFlags_;

  std::vector<HighsHashTree<HighsInt, HighsInt>> rowToVarBounds;

 public:
  const HighsMipSolver& mipsolver;
  std::vector<HighsSubstitution> substitutions;
  std::vector<HighsBool> colsubstituted;

  HighsImplications(const HighsMipSolver& mipsolver) : mipsolver(mipsolver) {
    nextCleanupCall = mipsolver.numNonzero();
    numImplications = 0;
    numVarBounds = 0;
    resize(mipsolver.numCol(), mipsolver.numRow());
  }

  std::function<void(HighsInt, HighsInt, HighsInt, double)>
      storeLiftingOpportunity;

  void reset() {
    colsubstituted.clear();
    colsubstituted.shrink_to_fit();
    implications.clear();
    implications.shrink_to_fit();
    hasProbed.clear();
    hasProbed.shrink_to_fit();
    reverseImplications.clear();
    reverseImplications.shrink_to_fit();
    numImplications = 0;
    vubs.clear();
    vubs.shrink_to_fit();
    vlbs.clear();
    vlbs.shrink_to_fit();
    dualFixProbingBinInds_.clear();
    rowToVarBounds.clear();
    rowToVarBounds.shrink_to_fit();
    resize(mipsolver.numCol(), mipsolver.numRow());
    numVarBounds = 0;
    nextCleanupCall = mipsolver.numNonzero();
  }

  void resize(HighsInt ncols, HighsInt nrows) {
    implications.resize(2 * static_cast<size_t>(ncols));
    hasProbed.resize(2 * static_cast<size_t>(ncols));
    reverseImplications.resize(ncols);
    colsubstituted.resize(ncols);
    vubs.resize(ncols);
    vlbs.resize(ncols);
    rowToVarBounds.resize(nrows);
    maxVarBounds = calcMaxVarBounds(ncols);
    dualFixProbingBinInds_.reserve(ncols);
    dualFixProbingBinFlags_.assign(ncols, TentativeFixing{});
  }

  constexpr static int64_t calcMaxVarBounds(HighsInt numcol) {
    return int64_t{5000000} + 10 * static_cast<int64_t>(numcol);
  };

  HighsInt getNumImplications() const {
    return static_cast<HighsInt>(numImplications);
  }

  bool probedBefore(const HighsInt col, const bool val) const {
    return hasProbed[ImplIdx{col, val}];
  }

  void addImplication(ImplIdx idx, HighsInt implCol, Implication implic);

  bool tooManyVarBounds() const { return numVarBounds >= maxVarBounds; }

  void strengthenVarBound(VarBound& vbnd, HighsInt multiplier) const;

  void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant,
              HighsInt origin = -1);

  void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant,
              double colupperbound, bool colisinteger, HighsInt origin = -1);

  void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant,
              HighsInt origin = -1);

  void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant,
              double collowerbound, bool colisinteger, HighsInt origin = -1);

  void columnTransformed(HighsInt col, double scale, double constant) {
    // Update implications affected by transformation
    auto changeImplications = [&](HighsInt binCol, bool) {
      for (HighsInt val = 0; val != 2; val++) {
        Implication* implic = implications[ImplIdx{binCol, val}].find(col);
        if (implic) {
          if (scale < 0) std::swap(implic->lb, implic->ub);
          implic->lb -= constant;
          implic->lb /= scale;
          implic->ub -= constant;
          implic->ub /= scale;
        }
      }
    };
    reverseImplications[col].for_each(changeImplications);

    // Update variable bounds affected by transformation
    if (scale < 0) std::swap(vubs[col], vlbs[col]);

    auto transformVbd = [&](HighsInt, VarBound& vbd) {
      vbd.constant -= constant;
      vbd.constant /= scale;
      vbd.coef /= scale;
    };

    vlbs[col].for_each(transformVbd);
    vubs[col].for_each(transformVbd);

    // Update substitutions affected by transformation
    for (auto& substitution : substitutions) {
      if (substitution.substcol == col) {
        substitution.offset -= constant;
        substitution.offset /= scale;
        substitution.scale /= scale;
      }
    }
  }

  const HighsHashTree<HighsInt, VarBound>& getVlbs(HighsInt col) const {
    return vlbs[col];
  }
  const HighsHashTree<HighsInt, VarBound>& getVubs(HighsInt col) const {
    return vubs[col];
  }

  void rowModified(HighsInt row);

  const HighsHashTree<HighsInt, HighsInt>& getRowVarBounds(HighsInt row) const {
    return rowToVarBounds[row];
  }

  std::pair<HighsInt, VarBound> getBestVub(HighsInt col,
                                           const HighsSolution& lpSolution,
                                           double& bestUb,
                                           const HighsDomain& globaldom) const;

  std::pair<HighsInt, VarBound> getBestVlb(HighsInt col,
                                           const HighsSolution& lpSolution,
                                           double& bestLb,
                                           const HighsDomain& globaldom) const;

  bool runProbing(HighsInt col, HighsInt& numReductions);

  void reindexRows(HighsInt nrows, const std::vector<HighsInt>& newRowIndex);

  void rebuild(HighsInt ncols, const std::vector<HighsInt>& cIndex,
               const std::vector<HighsInt>& rIndex);

  void buildFrom(const HighsImplications& init);

  void separateImpliedBounds(const HighsLpRelaxation& lpRelaxation,
                             const std::vector<double>& sol,
                             HighsCutPool& cutpool, double feastol,
                             HighsDomain& globaldom, bool thread_safe);

  void cleanupVarbounds(HighsInt col);

  bool redundantVlb(const VarBound& vlb, double lb) const;
  bool redundantVub(const VarBound& vub, double ub) const;
  bool tightenVlb(VarBound& vlb, double lb) const;
  bool tightenVub(VarBound& vub, double ub) const;

  void cleanupVlb(HighsInt col, HighsInt vlbCol,
                  HighsImplications::VarBound& vlb, double lb, bool& redundant,
                  bool& infeasible, bool allowBoundChanges = true) const;

  void cleanupVub(HighsInt col, HighsInt vubCol,
                  HighsImplications::VarBound& vub, double ub, bool& redundant,
                  bool& infeasible, bool allowBoundChanges = true) const;

  void applyImplications(HighsDomain& domain, HighsInt col, HighsInt val);

  void recordTentativeCliques(const HighsInt val,
                              const HighsDomainChange& domchg) {
    const HighsInt col = domchg.column;
    TentativeFixing& fixing = dualFixProbingBinFlags_[col];
    if (fixing.isUndecided()) dualFixProbingBinInds_.push_back(col);
    fixing.record(val == 1, domchg.boundtype);
  }

  void clearTentativeCliques() {
    for (const HighsInt col : dualFixProbingBinInds_)
      dualFixProbingBinFlags_[col].clear();
    dualFixProbingBinInds_.clear();
  }
};

#endif
