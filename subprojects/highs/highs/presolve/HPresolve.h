/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolve.h
 * @brief
 */
#ifndef PRESOLVE_HIGHS_PRESOLVE_H_
#define PRESOLVE_HIGHS_PRESOLVE_H_
#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "mip/HighsMipSolver.h"
#include "presolve/HPresolveAnalysis.h"
#include "presolve/HPresolveCliqueTable.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"
#include "util/HighsHashTree.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsMatrixSlice.h"

namespace presolve {

class HighsPostsolveStack;

class HPresolve {
  // pointer to model and options that where presolved
  HighsLp* model;
  const HighsOptions* options;
  HighsTimer* timer;
  HighsMipSolver* mipsolver = nullptr;
  HPresolveCliqueTable presolveCliqueTable;
  double primal_feastol;
  bool may_require_basis_postsolve_;
  bool may_require_dual_postsolve_;
  std::vector<HighsBool> allow_rule_;

  // triplet storage
  std::vector<double> Avalue;
  std::vector<HighsInt> Arow;
  std::vector<HighsInt> Acol;

  // linked list links for column based links for each nonzero
  std::vector<HighsInt> colhead;
  std::vector<HighsInt> Anext;
  std::vector<HighsInt> Aprev;

  // splay tree links for row based iteration and nonzero lookup
  std::vector<HighsInt> rowroot;
  std::vector<HighsInt> ARleft;
  std::vector<HighsInt> ARright;

  // length of rows and columns
  std::vector<HighsInt> rowsize;
  std::vector<HighsInt> rowsizeInteger;
  std::vector<HighsInt> rowsizeImplInt;
  // number of singleton columns in a row
  std::vector<HighsInt> rowsizeSingleton;
  std::vector<HighsInt> colsize;

  // vector to store the nonzero positions of a row
  std::vector<HighsInt> rowpositions;

  // stack to reuse free slots
  std::vector<HighsInt> freeslots;

  // vectors holding implied bounds on primal and dual variables as well as
  // their origins
  std::vector<double> implColLower;
  std::vector<double> implColUpper;
  std::vector<HighsInt> colLowerSource;
  std::vector<HighsInt> colUpperSource;
  std::vector<double> rowDualLower;
  std::vector<double> rowDualUpper;
  std::vector<double> implRowDualLower;
  std::vector<double> implRowDualUpper;
  std::vector<HighsInt> rowDualLowerSource;
  std::vector<HighsInt> rowDualUpperSource;
  std::vector<std::set<HighsInt>> colImplSourceByRow;
  std::vector<std::set<HighsInt>> implRowDualSourceByCol;

  // implied bounds on values of primal and dual rows computed from the bounds
  // of primal and dual variables
  HighsLinearSumBounds impliedRowBounds;
  HighsLinearSumBounds impliedDualRowBounds;

  std::vector<HighsInt> changedRowIndices;
  std::vector<HighsBool> changedRowFlag;
  std::vector<HighsInt> changedColIndices;
  std::vector<HighsBool> changedColFlag;

  std::vector<std::pair<HighsInt, HighsInt>> substitutionOpportunities;

  std::unordered_map<HighsInt,
                     HighsHashTree<std::pair<HighsInt, HighsInt>, double>>
      liftingOpportunities;

  // set with the sizes and indices of equation rows sorted by the size and a
  // vector to access there iterator positions in the set by index for quick
  // removal
  std::set<std::pair<HighsInt, HighsInt>> equations;
  std::vector<std::set<std::pair<HighsInt, HighsInt>>::iterator> eqiters;

  bool shrinkProblemEnabled;
  size_t reductionLimit;
  size_t last_reduction_;

  // vectors storing singleton rows and columns
  std::vector<HighsInt> singletonRows;
  std::vector<HighsInt> singletonColumns;

  // flags to mark rows/columns as deleted
  std::vector<HighsBool> rowDeleted;
  std::vector<HighsBool> colDeleted;

  // flags to skip repeated single-equation handling (dual fixing) on unchanged
  // rows
  std::vector<HighsBool> singleEquationChecked;

  std::vector<uint16_t> numProbes;

  int64_t probingContingent;
  HighsInt probingNumDelCol;
  HighsInt numProbed;

  // counters for number of deleted rows and columns
  HighsInt numDeletedRows;
  HighsInt numDeletedCols;

  // number of rows with singleEquationChecked set, so that marking a changed
  // column can skip clearing the flags of its rows when none is set
  HighsInt numSingleEquationChecked;

  // store old problem sizes to compute percentage reductions in
  // presolve loop
  HighsInt oldNumCol;
  HighsInt oldNumRow;
  bool probingEarlyAbort;

  enum class Result {
    kOk,
    kPrimalInfeasible,
    kDualInfeasible,
    kStopped,
    kOutOfMemory,
    kException
  };

  struct StatusResult {
   private:
    bool my_flag;
    Result my_result;

   public:
    // clang-format off
    StatusResult(bool flag) : my_flag(flag), my_result(Result::kOk) {};
    // clang-format on
    StatusResult(Result result) : my_result(result) {
      my_flag = (result == Result::kOk);
    };
    explicit operator bool() const { return my_flag; };
    explicit operator Result() const { return my_result; };
  };

  HighsPresolveStatus presolve_status_;
  HPresolveAnalysis analysis_;

  // private functions for different shared functionality and matrix
  // modification

  bool reducedToEmpty() const {
    return numDeletedCols == model->num_col_ &&
           numDeletedRows == model->num_row_;
  }

  bool hasPresolveDataStructures() const { return colDeleted.size() > 0; }

  void chooseRules();

  void link(HighsInt pos);

  void unlink(HighsInt pos);

  void markChangedRow(HighsInt row);

  void markChangedCol(HighsInt col);

  bool resetSingleEquationChecked(HighsInt row);

  double getMaxAbsColVal(HighsInt col) const;

  double getMaxAbsRowVal(HighsInt row) const;

  bool checkUpdateColImpliedBounds(HighsInt row, double* rowLower = nullptr,
                                   double* rowUpper = nullptr) const;

  Result updateColImpliedBounds(HighsInt row, HighsInt col, double val);

  Result updateColImpliedBounds(HighsInt row);

  bool checkUpdateRowDualImpliedBounds(HighsInt col,
                                       double* dualRowLower = nullptr,
                                       double* dualRowUpper = nullptr) const;

  void updateRowDualImpliedBounds(HighsInt row, HighsInt col, double val);

  void updateRowDualImpliedBounds(HighsInt col);

  void resetColImpliedBounds(HighsInt col, HighsInt row = -1);

  void resetRowDualImpliedBounds(HighsInt row, HighsInt col = -1);

  void resetColImpliedBoundsDerivedFromRow(HighsInt row);

  void resetRowDualImpliedBoundsDerivedFromCol(HighsInt col);

  void resetRowDualImpliedBoundsIfDualSideLost(HighsInt col, bool isLowerBound);

  void matrixNonZeroChanged(HighsInt row, HighsInt col, double oldCoef,
                            double newCoef, bool rowIsDeleted,
                            bool colIsDeleted);

  void changeRowLower(HighsInt row, double newLower,
                      bool skipRowDualUpdate = false);

  void changeRowUpper(HighsInt row, double newUpper,
                      bool skipRowDualUpdate = false);

  void addToRowLower(HighsInt row, const HighsCDouble& delta);

  void addToRowUpper(HighsInt row, const HighsCDouble& delta);

  bool rowCoefficientsIntegral(HighsInt row, double scale) const;

  bool isBinary(HighsInt col) const;

  bool isImpliedFree(HighsInt col) const;

  bool isDualImpliedFree(HighsInt row) const;

  bool isDualImpliedLeq(HighsInt row) const;

  bool isDualImpliedLeq(HighsInt row, double implUpper) const;

  bool isDualImpliedGeq(HighsInt row) const;

  bool isDualImpliedGeq(HighsInt row, double implLower) const;

  void dualImpliedFreeGetRhsAndRowType(HighsInt row, double& rhs,
                                       HighsPostsolveStack::RowType& rowType,
                                       bool relaxRowDualBounds = false);

  bool isEquation(HighsInt row) const;

  bool isRanged(HighsInt row) const;

  bool isRedundant(HighsInt row, double sumLower, double sumUpper) const;

  bool isRedundant(HighsInt row) const;

  bool yieldsImpliedLowerBound(HighsInt row, double val) const;

  bool yieldsImpliedUpperBound(HighsInt row, double val) const;

  bool isImpliedEquationAtLower(HighsInt row) const;

  bool isImpliedEquationAtUpper(HighsInt row) const;

  StatusResult isImpliedIntegral(HighsInt col);

  StatusResult isImpliedInteger(HighsInt col) const;

  StatusResult convertImpliedInteger(HighsInt col, HighsInt row = -1,
                                     bool skipInputChecks = false);

  bool isLowerImplied(HighsInt col) const;

  bool isLowerImplied(HighsInt col, double implLower) const;

  bool isLowerStrictlyImplied(HighsInt col, double* tolerance = nullptr) const;

  bool isUpperImplied(HighsInt col) const;

  bool isUpperImplied(HighsInt col, double implUpper) const;

  bool isUpperStrictlyImplied(HighsInt col, double* tolerance = nullptr) const;

  bool isStrictlyImpliedForDual(HighsInt col, bool isLowerBound) const;

  HighsInt countFillin(HighsInt row);

  void reinsertEquation(HighsInt row);

  void clearLiftingOpportunities(HighsInt row) {
    auto search = liftingOpportunities.find(row);
    if (search != liftingOpportunities.end()) search->second.clear();
  };

#ifndef NDEBUG
  void debugPrintRow(HighsPostsolveStack& postsolve_stack, HighsInt row);
#endif

  HighsInt findNonzero(HighsInt row, HighsInt col);

  bool okFromCSC(const std::vector<double>& Aval,
                 const std::vector<HighsInt>& Aindex,
                 const std::vector<HighsInt>& Astart);

  bool okFromCSR(const std::vector<double>& ARval,
                 const std::vector<HighsInt>& ARindex,
                 const std::vector<HighsInt>& ARstart);

  void toCSC(std::vector<double>& Aval, std::vector<HighsInt>& Aindex,
             std::vector<HighsInt>& Astart);

  void getRowPositions(HighsInt row,
                       std::vector<HighsInt>& myrowpositions) const;

  void storeRow(HighsInt row);

  HighsTripletPositionSlice getStoredRow() const;

  HighsTripletListSlice getColumnVector(HighsInt col) const;

  HighsTripletTreeSlicePreOrder getRowVector(HighsInt row) const;

  HighsTripletTreeSliceInOrder getSortedRowVector(HighsInt row) const;

  void markRowDeleted(HighsInt row);

  void markColDeleted(HighsInt col);

  Result fixColToLower(HighsPostsolveStack& postsolve_stack, HighsInt col);

  Result fixColToUpper(HighsPostsolveStack& postsolve_stack, HighsInt col);

  Result fixColToZero(HighsPostsolveStack& postsolve_stack, HighsInt col);

  Result transformColumn(HighsPostsolveStack& postsolve_stack, HighsInt col,
                         double scale, double constant);

  void scaleRow(HighsInt row, double scale, bool integral = false);

  void scaleStoredRow(HighsInt row, double scale, bool integral = false);

  void substitute(HighsInt row, HighsInt col, double rhs);

  Result changeColUpper(HighsInt col, double newUpper);

  Result changeColLower(HighsInt col, double newLower);

  Result changeColBounds(HighsInt col, double newLower, double newUpper);

  Result checkColBounds(HighsInt col, bool* isFixed = nullptr);

  void changeRowDualUpper(HighsInt row, double newUpper);

  void changeRowDualLower(HighsInt row, double newLower);

  void changeImplColUpper(HighsInt col, double newUpper, HighsInt originRow);

  void changeImplColLower(HighsInt col, double newLower, HighsInt originRow);

  void changeImplRowDualUpper(HighsInt row, double newUpper,
                              HighsInt originCol);

  void changeImplRowDualLower(HighsInt row, double newLower,
                              HighsInt originCol);

  void moveCutsToPool(HighsPostsolveStack& postsolve_stack);

  Result scaleMIP(HighsPostsolveStack& postsolve_stack);

  Result applyConflictGraphSubstitutions(HighsPostsolveStack& postsolve_stack,
                                         HighsInt& numDelCol);

  Result fastPresolveLoop(HighsPostsolveStack& postsolve_stack);

  Result presolve(HighsPostsolveStack& postsolve_stack);

  Result checkTimeLimit();

  Result checkLimits(HighsPostsolveStack& postsolve_stack);

  void storeCurrentProblemSize();

  double problemSizeReduction() const;

  void computeColBounds(HighsInt col, HighsInt boundCol = -1,
                        double boundColValue = -kHighsInf,
                        HighsInt boundColCoeffPattern = 0,
                        double* lowerBound = nullptr,
                        double* upperBound = nullptr,
                        double* worstCaseLowerBound = nullptr,
                        double* worstCaseUpperBound = nullptr);

  template <typename storageFormat>
  HighsCDouble computeDynamism(const HighsMatrixSlice<storageFormat>& vector) {
    double minAbsCoef = kHighsInf;
    double maxAbsCoef = -kHighsInf;
    for (const auto& nonzero : vector) {
      double absCoef = std::abs(nonzero.value());
      minAbsCoef = std::min(minAbsCoef, absCoef);
      maxAbsCoef = std::max(maxAbsCoef, absCoef);
    }
    return static_cast<HighsCDouble>(maxAbsCoef) /
           static_cast<HighsCDouble>(minAbsCoef);
  }

  bool silentLog() const;

 public:
  // for LP presolve
  void setInput(HighsLp& model_, const HighsOptions& options_,
                const HighsInt presolve_reduction_limit,
                HighsTimer* timer = nullptr);

  // for MIP presolve
  void setInput(HighsMipSolver& mipsolver,
                const HighsInt presolve_reduction_limit);

  bool okSetupPresolveDataStructures();
  void setupSubstitutionOpportunities();

  HighsInt numNonzeros() const { return int(Avalue.size() - freeslots.size()); }

  void shrinkProblem(HighsPostsolveStack& postsolve_stack);

  void addToMatrix(const HighsInt row, const HighsInt col, const double val);

  bool addToMatrix(HighsPostsolveStack& postsolve_stack,
                   const std::vector<double>& row_lower,
                   const std::vector<double>& row_upper,
                   const std::vector<std::vector<HighsInt>>& row_indices,
                   const std::vector<std::vector<double>>& row_values);

  bool addToMatrix(HighsPostsolveStack& postsolve_stack, double row_lower,
                   double row_upper, const std::vector<HighsInt>& row_indices,
                   const std::vector<double>& row_values);

  Result normaliseCliqueRows(HighsPostsolveStack& postsolve_stack);

  Result prepareProbing(HighsPostsolveStack& postsolve_stack, bool& firstCall);

  Result finaliseProbing(HighsPostsolveStack& postsolve_stack, bool firstCall,
                         HighsInt& numVarsFixed, HighsInt& numBndsTightened,
                         HighsInt& numVarsSubstituted,
                         HighsInt& liftedNonZeros);

  std::pair<int64_t, HighsInt> computeProbingScore(HighsInt col) const;

  Result runProbing(HighsPostsolveStack& postsolve_stack);

  Result liftingForProbing(HighsPostsolveStack& postsolve_stack);

  Result dominatedColumns(HighsPostsolveStack& postsolve_stack);

  Result doubletonEq(HighsPostsolveStack& postsolve_stack, HighsInt row,
                     HighsPostsolveStack::RowType rowType);

  Result singletonRow(HighsPostsolveStack& postsolve_stack, HighsInt row);

  Result emptyCol(HighsPostsolveStack& postsolve_stack, HighsInt col);

  Result singletonCol(HighsPostsolveStack& postsolve_stack, HighsInt col,
                      const bool timing = false);

  void substituteFreeCol(HighsPostsolveStack& postsolve_stack, HighsInt row,
                         HighsInt col, bool relaxRowDualBounds = false);

  Result emptyRow(HighsPostsolveStack& postsolve_stack, HighsInt row);

  Result rowPresolve(HighsPostsolveStack& postsolve_stack, HighsInt row);

  Result colPresolve(HighsPostsolveStack& postsolve_stack, HighsInt col,
                     const bool timing = false);

  Result detectDominatedCol(HighsPostsolveStack& postsolve_stack, HighsInt col,
                            bool handleSingletonRows = true);

  void computeLocks(
      HighsInt col, bool considerObjective,
      const std::function<bool(HighsInt, bool, bool)>& lockCallback) const;

  Result dualFixing(HighsPostsolveStack& postsolve_stack, HighsInt col);

  Result singletonColStuffing(HighsPostsolveStack& postsolve_stack,
                              HighsInt col);

  Result zeroCostSingleton(HighsPostsolveStack& postsolve_stack, HighsInt col);

  Result enumerateSolutions(HighsPostsolveStack& postsolve_stack);

  double computeImpliedLowerBound(HighsInt col, HighsInt boundCol = -1,
                                  double boundColValue = kHighsInf,
                                  HighsInt boundColCoeffPattern = 0);

  double computeImpliedUpperBound(HighsInt col, HighsInt boundCol = -1,
                                  double boundColValue = kHighsInf,
                                  HighsInt boundColCoeffPattern = 0);

  double computeWorstCaseLowerBound(HighsInt col, HighsInt boundCol = -1,
                                    double boundColValue = kHighsInf,
                                    HighsInt boundColCoeffPattern = 0);

  double computeWorstCaseUpperBound(HighsInt col, HighsInt boundCol = -1,
                                    double boundColValue = kHighsInf,
                                    HighsInt boundColCoeffPattern = 0);

  Result checkOriginalModelBounds();

  Result initialRowAndColPresolve(HighsPostsolveStack& postsolve_stack);

  HighsModelStatus run(HighsPostsolveStack& postsolve_stack);

  Result substitute(HighsInt substcol, HighsInt staycol, double offset,
                    double scale, HighsInt row = -1);

  Result removeFixedCol(HighsInt col);

  Result removeFixedCol(HighsInt col, double fixval);

  void unlinkRow(HighsInt row);

  void removeRow(HighsInt row);

  Result removeDependentEquations(HighsPostsolveStack& postsolve_stack);

  Result removeDependentFreeCols(HighsPostsolveStack& postsolve_stack);

  Result aggregator(HighsPostsolveStack& postsolve_stack);

  Result fourierMotzkin(HighsPostsolveStack& postsolve_stack,
                        HighsInt& numColsEliminated);

  Result removeRowSingletons(HighsPostsolveStack& postsolve_stack);

  Result presolveColSingletons(HighsPostsolveStack& postsolve_stack);

  Result presolveChangedRows(HighsPostsolveStack& postsolve_stack);

  Result presolveChangedCols(HighsPostsolveStack& postsolve_stack);

  Result removeDoubletonEquations(HighsPostsolveStack& postsolve_stack);

  Result strengthenInequalities(HighsPostsolveStack& postsolve_stack,
                                HighsInt& num_strengthened);

  Result detectImpliedIntegers();

  Result detectParallelRowsAndCols(HighsPostsolveStack& postsolve_stack);

  template <typename RowStorageFormat>
  Result equalityRowAddition(HighsPostsolveStack& postsolve_stack,
                             HighsInt stayrow, HighsInt removerow, double scale,
                             const HighsMatrixSlice<RowStorageFormat>& vector);

  void extractVarBounds(HighsInt row);

  Result updateCliqueTableFixedCol(HighsInt col, double val);

  Result updateCliqueTableSubstituteCol(HighsInt substCol, HighsInt stayCol,
                                        double offset, double scale);

  void aggregateVarBounds();

  Result sparsify(HighsPostsolveStack& postsolve_stack);

  const HighsPresolveLog& getPresolveLog() const {
    return analysis_.presolve_log_;
  }

  HighsPresolveStatus getPresolveStatus() const { return presolve_status_; }

  bool failed() const {
    return presolve_status_ == HighsPresolveStatus::kOutOfMemory ||
           presolve_status_ == HighsPresolveStatus::kException ||
           presolve_status_ == HighsPresolveStatus::kMatrixValueError;
  }

  bool zeroRowActivityFeasible() const;

  HighsInt debugGetCheckCol() const;
  HighsInt debugGetCheckRow() const;

  Result presolveRuleTest(HighsPostsolveStack& postsolve_stack);
  Result presolveRuleTestColStuffing(HighsPostsolveStack& postsolve_stack);
  Result presolveRuleTestParallelRowsAndCols(
      HighsPostsolveStack& postsolve_stack);
  Result presolveRuleTestProbing(HighsPostsolveStack& postsolve_stack);
  Result presolveRuleTestFourierMotzkin(HighsPostsolveStack& postsolve_stack);
  Result presolveRuleTestAggregator(HighsPostsolveStack& postsolve_stack);

  /*
  // Methods defined and used in HPresolveDebug, and only executed if
  // HPresolve::debug is called. This hasn't been used for ages, and
  // this comment is retained in case it's useful in future
  static void debug(const HighsLp& lp, const HighsOptions& options);
  void computeIntermediateMatrix(std::vector<HighsInt>& flagRow,
                                 std::vector<HighsInt>& flagCol,
                                 size_t& numreductions);
  void setReductionLimit(size_t reductionLimit) {
    this->reductionLimit = reductionLimit;
  }
  */
};

}  // namespace presolve
#endif
