/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsPostsolveStack.h
 * @brief Class to hold all information for postsolve and can transform back
 * primal and dual solutions.
 */

#ifndef PRESOLVE_HIGHS_POSTSOLVE_STACK_H_
#define PRESOLVE_HIGHS_POSTSOLVE_STACK_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <sstream>
#include <tuple>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsCDouble.h"
#include "util/HighsDataStack.h"
#include "util/HighsMatrixSlice.h"

// class HighsOptions;
namespace presolve {
class HighsPostsolveStack {
  // now a section of individual classes for each type of each transformation
  // step that requires postsolve starts each class gets as first argument the
  // current stack of ReductionValues and custom arguments that contain the
  // necessary information that is required to undo the transformation. The
  // constructor is responsible for storing all necessary information in class
  // members and the reduction value stack. The class members should be as slim
  // as possible and putting values on the reduction value stack should be
  // preferred, because the classes are stored in a discriminated union and the
  // largest size counts. The classes should implement an undo() function which
  // gets the ReductionValues as argument and can be expected to be called such
  // that the stack is in the state as after the constructor has been called.
  // The undo() call must pop all values from the stack that were added during
  // the constructor call, and should restore primal/dual solution values, as
  // well as the basis status as appropriate.
 public:
  enum class OrigRowType : uint8_t { kOriginal, kCut, kAppended };

  enum class RowType {
    kGeq,
    kLeq,
    kEq,
  };
  struct Nonzero {
    HighsInt index;
    double value;

    Nonzero(HighsInt index_, double value_) : index(index_), value(value_) {}
    Nonzero() = default;
  };

  struct FmeRowHeader {
    HighsInt row;
    double rowLower;
    double rowUpper;
  };

  template <typename RowStorageFormat>
  struct FmeRowData {
    FmeRowHeader header;
    HighsMatrixSlice<RowStorageFormat> rowVec;
  };

  struct FmeStepHeader {
    double colLower;
    double colUpper;
    HighsInt col;
    HighsInt numPlus;
    HighsInt numMinus;
  };

  struct FmeDescendant {
    HighsInt row;
    double scaleFactor;
  };

  struct FmeNewRow {
    HighsInt row;
    HighsInt plusParentIdx;
    HighsInt minusParentIdx;
  };

  struct FmeAncestryEntry {
    HighsInt step;
    HighsInt parentRowIndex;
    double scale;
    bool isMinus;
  };

  struct FmeBlockStep {
    FmeStepHeader header;
    std::vector<FmeNewRow> newRows;
  };

  size_t debug_prev_numreductions = 0;
  double debug_prev_col_lower = 0;
  double debug_prev_col_upper = 0;
  double debug_prev_row_lower = 0;
  double debug_prev_row_upper = 0;

 private:
  /// transform a column x by a linear mapping with a new column x'.
  /// I.e. substitute x = a * x' + b
  struct LinearTransform {
    double scale;
    double constant;
    HighsInt col;

    void undo(const HighsOptions& options, HighsSolution& solution) const;

    void transformToPresolvedSpace(std::vector<double>& primalSol) const;
  };

  struct FourierMotzkinObjCol {
    double offset;
    HighsInt col;

    void transformToPresolvedSpace(const std::vector<Nonzero>& costEntries,
                                   std::vector<double>& primalSol) const;

    void undo(const std::vector<Nonzero>& costEntries,
              HighsSolution& solution) const;
  };

  struct FreeColSubstitution {
    double rhs;
    double colCost;
    HighsInt row;
    HighsInt col;
    RowType rowType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  struct DoubletonEquation {
    double coef;
    double coefSubst;
    double rhs;
    double substLower;
    double substUpper;
    double substCost;
    HighsInt row;
    HighsInt colSubst;
    HighsInt col;
    bool lowerTightened;
    bool upperTightened;
    RowType rowType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct EqualityRowAddition {
    HighsInt row;
    HighsInt addedEqRow;
    double eqRowScale;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& eqRowValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct EqualityRowAdditions {
    HighsInt addedEqRow;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& eqRowValues,
              const std::vector<Nonzero>& targetRows, HighsSolution& solution,
              HighsBasis& basis) const;
  };
  struct SingletonRow {
    double coef;
    HighsInt row;
    HighsInt col;
    bool colLowerTightened;
    bool colUpperTightened;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  // column fixed to lower or upper bound
  struct FixedCol {
    double fixValue;
    double colCost;
    HighsInt col;
    HighsBasisStatus fixType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct RedundantRow {
    HighsInt row;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct ImpliedEquation {
    HighsInt row;
    bool atLower;

    void undo(const std::vector<Nonzero>& rowValues,
              HighsSolution& solution) const;
  };

  struct ForcingRow {
    double side;
    HighsInt row;
    RowType rowType;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct ForcingColumn {
    double colCost;
    double colBound;
    HighsInt col;
    bool atInfiniteUpper;
    bool colIntegral;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& colValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct ForcingColumnRemovedRow {
    double rhs;
    HighsInt row;
    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct DuplicateRow {
    double duplicateRowScale;
    HighsInt duplicateRow;
    HighsInt row;
    bool rowLowerTightened;
    bool rowUpperTightened;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
  };

  struct DuplicateColumn {
    double colScale;
    double colLower;
    double colUpper;
    double duplicateColLower;
    double duplicateColUpper;
    HighsInt col;
    HighsInt duplicateCol;
    bool colIntegral;
    bool duplicateColIntegral;

    void undo(const HighsOptions& options, HighsSolution& solution,
              HighsBasis& basis) const;
    bool okMerge(const double tolerance) const;
    void undoFix(const HighsOptions& options, HighsSolution& solution,
                 const double mergeValue) const;
    void transformToPresolvedSpace(std::vector<double>& primalSol) const;
  };

  struct ZeroObjSingletonContinuousCol {
    double origRowLower;
    double origRowUpper;
    double new_row_lb;
    double new_row_ub;
    double lb;
    double ub;
    double coef;
    HighsInt col;
    HighsInt row;

    void undo(const HighsOptions& options,
              const std::vector<Nonzero>& rowValues, HighsSolution& solution,
              HighsBasis& basis);
  };

  /// tags for reduction
  enum class ReductionType : uint8_t {
    kLinearTransform,
    kFreeColSubstitution,
    kDoubletonEquation,
    kEqualityRowAddition,
    kEqualityRowAdditions,
    kSingletonRow,
    kFixedCol,
    kRedundantRow,
    kImpliedEquation,
    kForcingRow,
    kForcingColumn,
    kForcingColumnRemovedRow,
    kDuplicateRow,
    kDuplicateColumn,
    kZeroObjSingletonContinuousCol,
    kFourierMotzkinBlock,
    kFourierMotzkinObjCol,
  };

  struct FmeParentRow {
    FmeRowHeader header;
    double coef;
    std::vector<Nonzero> entries;
    std::vector<FmeDescendant> descendants;
  };

  struct FmeStepData {
    FmeStepHeader header;
    std::vector<FmeParentRow> plusRows;
    std::vector<FmeParentRow> minusRows;
    std::vector<FmeNewRow> newRows;
  };

  static std::vector<FmeStepData> popFourierMotzkinBlock(HighsDataStack& stack);
  static void undoFourierMotzkinBlock(const std::vector<FmeStepData>& steps,
                                      const HighsOptions& options,
                                      HighsSolution& solution,
                                      HighsBasis& basis);

  HighsDataStack reductionValues;
  std::vector<std::pair<ReductionType, size_t>> reductions;
  std::vector<HighsInt> origColIndex;
  std::vector<HighsInt> origRowIndex;
  std::vector<OrigRowType> origRowType;
  std::vector<HighsBool> linearlyTransformable;

  std::vector<Nonzero> rowValues;
  std::vector<Nonzero> colValues;
  HighsInt origNumCol = -1;
  HighsInt origNumRow = -1;
  HighsInt numAppendedRows = 0;
  HighsInt nextRowIndex = -1;
  HighsInt nextColIndex = -1;

  void reductionAdded(ReductionType type) {
    size_t position = reductionValues.getCurrentDataSize();
    reductions.emplace_back(type, position);
  }

  std::string presolveTypeToString(const ReductionType& type) const {
    switch (type) {
      case ReductionType::kLinearTransform: {
        return "Linear transform";
      }
      case ReductionType::kFreeColSubstitution: {
        return "Free col substitution";
      }
      case ReductionType::kDoubletonEquation: {
        return "Doubleton equation";
      }
      case ReductionType::kEqualityRowAddition: {
        return "Equality row addition";
      }
      case ReductionType::kEqualityRowAdditions: {
        return "Equality row additions";
      }
      case ReductionType::kSingletonRow: {
        return "Singleton row";
      }
      case ReductionType::kFixedCol: {
        return "Fixed col";
      }
      case ReductionType::kRedundantRow: {
        return "Redundant Row";
      }
      case ReductionType::kForcingRow: {
        return "Forcing row";
      }
      case ReductionType::kForcingColumn: {
        return "Forcing column";
      }
      case ReductionType::kForcingColumnRemovedRow: {
        return "Forcing column removed row";
      }
      case ReductionType::kDuplicateRow: {
        return "Duplicate row";
      }
      case ReductionType::kDuplicateColumn: {
        return "Duplicate column";
      }
      case ReductionType::kImpliedEquation: {
        return "Implied equation";
      }
      case ReductionType::kZeroObjSingletonContinuousCol: {
        return "Zero obj singleton continuous col";
      }
      case ReductionType::kFourierMotzkinBlock: {
        return "Fourier-Motzkin block";
      }
      case ReductionType::kFourierMotzkinObjCol: {
        return "Fourier-Motzkin obj col";
      }
      default:
        return "Unknown";
    }
  }

 public:
  const std::vector<HighsInt>& getOrigColIndex() const { return origColIndex; }
  HighsInt getOrigColIndex(HighsInt col) const { return origColIndex[col]; }

  const std::vector<HighsInt>& getOrigRowIndex() const { return origRowIndex; }
  HighsInt getOrigRowIndex(HighsInt row) const { return origRowIndex[row]; }

  bool isOrigCol(HighsInt col) const { return origColIndex[col] < origNumCol; }

  // Returns presolved-space indices of columns from the original model
  std::vector<HighsInt> getOrigCols() const {
    std::vector<HighsInt> cols;
    for (HighsInt i = 0; i < static_cast<HighsInt>(origColIndex.size()); ++i)
      if (isOrigCol(i)) cols.push_back(i);
    return cols;
  }

  bool isOrigRow(HighsInt row) const {
    return origRowType[row] == OrigRowType::kOriginal;
  }

  bool isAppendedRow(HighsInt row) const {
    return origRowType[row] == OrigRowType::kAppended;
  }

  bool isCutRow(HighsInt row) const {
    return origRowType[row] == OrigRowType::kCut;
  }

  void setRowType(HighsInt row, OrigRowType type) { origRowType[row] = type; }

  // Returns presolved-space indices of rows from the original model
  std::vector<HighsInt> getOrigRows() const {
    std::vector<HighsInt> rows;
    for (HighsInt i = 0; i < static_cast<HighsInt>(origRowType.size()); ++i)
      if (isOrigRow(i)) rows.push_back(i);
    return rows;
  }

  // Returns presolved-space indices of rows that are cuts
  std::vector<HighsInt> getCutRows() const {
    std::vector<HighsInt> rows;
    // keep reverse loop to avoid behavior changes
    for (HighsInt i = static_cast<HighsInt>(origRowType.size()) - 1; i >= 0;
         --i)
      if (isCutRow(i)) rows.push_back(i);
    return rows;
  }

  // Returns presolved-space indices of rows that are not cuts (original +
  // appended)
  std::vector<HighsInt> getNonCutRows() const {
    std::vector<HighsInt> rows;
    for (HighsInt i = 0; i < static_cast<HighsInt>(origRowType.size()); ++i)
      if (!isCutRow(i)) rows.push_back(i);
    return rows;
  }

  void appendToModel(HighsInt& numRows, HighsInt numRowsToAppend,
                     OrigRowType rowType) {
    if (numRowsToAppend <= 0) return;
    size_t currNumRow = origRowIndex.size();
    size_t newNumRow = currNumRow + numRowsToAppend;
    origRowIndex.resize(newNumRow);
    origRowType.resize(newNumRow, rowType);
    for (size_t i = currNumRow; i != newNumRow; ++i)
      origRowIndex[i] = nextRowIndex++;
    numRows += numRowsToAppend;
  }

  HighsInt getPresolvedIndex(
      const HighsInt original_index,
      const std::vector<HighsInt>& original_indices) const {
    auto it = std::find(original_indices.begin(), original_indices.end(),
                        original_index);
    // If found, return the index after subtracting the base iterator
    if (it != original_indices.end())
      return static_cast<HighsInt>(it - original_indices.begin());
    return -1;
  }

  HighsInt getPresolvedColumnIndex(const HighsInt original_col_index) const {
    return getPresolvedIndex(original_col_index, origColIndex);
  }

  HighsInt getPresolvedRowIndex(const HighsInt original_row_index) const {
    return getPresolvedIndex(original_row_index, origRowIndex);
  }

  void appendCutsToModel(HighsInt numCuts) {
    appendToModel(origNumRow, numCuts, OrigRowType::kCut);
  }

  void appendRowsToModel(HighsInt numRows) {
    appendToModel(numAppendedRows, numRows, OrigRowType::kAppended);
  }

  void removeCutsFromModel(HighsInt numCuts) {
    if (numCuts <= 0) return;
    origNumRow -= numCuts;
    size_t newSize = 0;
    for (size_t i = 0; i < origRowIndex.size(); ++i) {
      if (origRowType[i] != OrigRowType::kCut) {
        if (i != newSize) {
          origRowIndex[newSize] = origRowIndex[i];
          origRowType[newSize] = origRowType[i];
        }
        ++newSize;
      }
    }
    origRowIndex.resize(newSize);
    origRowType.resize(newSize);
  }

  HighsInt getOrigNumRow() const { return origNumRow; }

  HighsInt getOrigNumCol() const { return origNumCol; }

  HighsInt getNextRowIndex() const { return nextRowIndex; }

  HighsInt getNextColIndex() const { return nextColIndex; }

  void appendColToModel() {
    origColIndex.push_back(nextColIndex++);
    linearlyTransformable.push_back(false);
  }

  void initializeIndexMaps(HighsInt numRow, HighsInt numCol);

  void compressIndexMaps(const std::vector<HighsInt>& newRowIndex,
                         const std::vector<HighsInt>& newColIndex);

  void compressRowIndexMap(const std::vector<HighsInt>& newRowIndex);
  void compressColIndexMap(const std::vector<HighsInt>& newColIndex);

  template <typename T>
  void compressIndexMap(const std::vector<HighsInt>& newIndex,
                        std::vector<T>& origIndex) {
    size_t numEn = origIndex.size();
    for (size_t i = 0; i != newIndex.size(); ++i) {
      if (newIndex[i] == -1)
        --numEn;
      else
        origIndex[newIndex[i]] = origIndex[i];
    }
    origIndex.resize(numEn);
  }

  /// transform a column x by a linear mapping with a new column x'.
  /// I.e. substitute x = scale * x' + constant
  void linearTransform(HighsInt col, double scale, double constant) {
    reductionValues.push(LinearTransform{scale, constant, origColIndex[col]});
    reductionAdded(ReductionType::kLinearTransform);
  }

  template <typename RowStorageFormat, typename ColStorageFormat>
  void freeColSubstitution(HighsInt row, HighsInt col, double rhs,
                           double colCost, RowType rowType,
                           const HighsMatrixSlice<RowStorageFormat>& rowVec,
                           const HighsMatrixSlice<ColStorageFormat>& colVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FreeColSubstitution{rhs, colCost, origRowIndex[row],
                                             origColIndex[col], rowType});
    reductionValues.push(rowValues);
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kFreeColSubstitution);
  }

  template <typename RowStorageFormat>
  void zeroCostSingleton(HighsInt row, HighsInt col, double origRowLower,
                         double origRowUpper, double new_row_lb,
                         double new_row_ub, double lb, double ub, double coef,
                         const HighsMatrixSlice<RowStorageFormat>& rowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(ZeroObjSingletonContinuousCol{
        origRowLower, origRowUpper, new_row_lb, new_row_ub, lb, ub, coef,
        origColIndex[col], origRowIndex[row]});
    reductionValues.push(rowValues);
    reductionAdded(ReductionType::kZeroObjSingletonContinuousCol);
  }

  template <typename RowStorageFormat>
  void fourierMotzkinBlockPushStep(
      HighsInt col, const std::vector<FmeRowData<RowStorageFormat>>& plusRows,
      const std::vector<FmeRowData<RowStorageFormat>>& minusRows) {
    std::vector<FmeRowHeader> plusHeaders;
    std::vector<double> plusCoefs;
    plusHeaders.reserve(plusRows.size());
    plusCoefs.reserve(plusRows.size());
    for (const auto& rd : plusRows) {
      std::vector<Nonzero> translated;
      double coef = 0.0;
      for (const HighsSliceNonzero& nz : rd.rowVec) {
        if (nz.index() == col)
          coef = nz.value();
        else
          translated.push_back({origColIndex[nz.index()], nz.value()});
      }
      reductionValues.push(translated);
      plusCoefs.push_back(coef);
      plusHeaders.push_back({origRowIndex[rd.header.row], rd.header.rowLower,
                             rd.header.rowUpper});
    }
    reductionValues.push(plusCoefs);
    reductionValues.push(plusHeaders);

    std::vector<FmeRowHeader> minusHeaders;
    std::vector<double> minusCoefs;
    minusHeaders.reserve(minusRows.size());
    minusCoefs.reserve(minusRows.size());
    for (const auto& rd : minusRows) {
      std::vector<Nonzero> translated;
      double coef = 0.0;
      for (const HighsSliceNonzero& nz : rd.rowVec) {
        if (nz.index() == col)
          coef = nz.value();
        else
          translated.push_back({origColIndex[nz.index()], nz.value()});
      }
      reductionValues.push(translated);
      minusCoefs.push_back(coef);
      minusHeaders.push_back({origRowIndex[rd.header.row], rd.header.rowLower,
                              rd.header.rowUpper});
    }
    reductionValues.push(minusCoefs);
    reductionValues.push(minusHeaders);
  }

  void fourierMotzkinBlockFinalise(
      const std::vector<FmeBlockStep>& blockSteps,
      const std::unordered_map<HighsInt, std::vector<FmeAncestryEntry>>&
          rowAncestry) {
    HighsInt numSteps = static_cast<HighsInt>(blockSteps.size());

    std::vector<std::vector<std::vector<FmeDescendant>>> plusDescendantsAll(
        numSteps);
    std::vector<std::vector<std::vector<FmeDescendant>>> minusDescendantsAll(
        numSteps);
    for (HighsInt s = 0; s < numSteps; ++s) {
      plusDescendantsAll[s].resize(blockSteps[s].header.numPlus);
      minusDescendantsAll[s].resize(blockSteps[s].header.numMinus);
    }
    for (const auto& entry : rowAncestry) {
      HighsInt row = entry.first;
      HighsInt origRow = origRowIndex[row];
      for (const auto& a : entry.second) {
        if (a.isMinus)
          minusDescendantsAll[a.step][a.parentRowIndex].push_back(
              {origRow, a.scale});
        else
          plusDescendantsAll[a.step][a.parentRowIndex].push_back(
              {origRow, a.scale});
      }
    }

    for (HighsInt s = 0; s < numSteps; ++s) {
      assert(static_cast<HighsInt>(plusDescendantsAll[s].size()) ==
             blockSteps[s].header.numPlus);
      for (HighsInt p = 0; p < blockSteps[s].header.numPlus; ++p)
        reductionValues.push(plusDescendantsAll[s][p]);
      assert(static_cast<HighsInt>(minusDescendantsAll[s].size()) ==
             blockSteps[s].header.numMinus);
      for (HighsInt m = 0; m < blockSteps[s].header.numMinus; ++m)
        reductionValues.push(minusDescendantsAll[s][m]);
    }

    for (HighsInt s = 0; s < numSteps; ++s) {
      std::vector<FmeNewRow> translated;
      translated.reserve(blockSteps[s].newRows.size());
      for (const auto& nr : blockSteps[s].newRows)
        translated.push_back(
            {origRowIndex[nr.row], nr.plusParentIdx, nr.minusParentIdx});
      reductionValues.push(translated);
    }

    for (HighsInt s = 0; s < numSteps; ++s) {
      FmeStepHeader header{
          blockSteps[s].header.colLower, blockSteps[s].header.colUpper,
          origColIndex[blockSteps[s].header.col], blockSteps[s].header.numPlus,
          blockSteps[s].header.numMinus};
      reductionValues.push(header);
    }

    reductionValues.push(numSteps);
    reductionAdded(ReductionType::kFourierMotzkinBlock);
  }

  void fourierMotzkinObjCol(HighsInt col, double offset,
                            const std::vector<Nonzero>& costEntries) {
    reductionValues.push(FourierMotzkinObjCol{offset, origColIndex[col]});
    std::vector<Nonzero> translatedEntries;
    translatedEntries.reserve(costEntries.size());
    for (const Nonzero& entry : costEntries)
      if (entry.index != col)
        translatedEntries.emplace_back(origColIndex[entry.index], entry.value);
    reductionValues.push(translatedEntries);
    reductionAdded(ReductionType::kFourierMotzkinObjCol);
  }

  template <typename ColStorageFormat>
  void doubletonEquation(HighsInt row, HighsInt colSubst, HighsInt col,
                         double coefSubst, double coef, double rhs,
                         double substLower, double substUpper, double substCost,
                         bool lowerTightened, bool upperTightened,
                         RowType rowType,
                         const HighsMatrixSlice<ColStorageFormat>& colVec) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(DoubletonEquation{
        coef, coefSubst, rhs, substLower, substUpper, substCost,
        row == -1 ? -1 : origRowIndex[row], origColIndex[colSubst],
        origColIndex[col], lowerTightened, upperTightened, rowType});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kDoubletonEquation);
  }

  template <typename RowStorageFormat>
  void equalityRowAddition(HighsInt row, HighsInt addedEqRow, double eqRowScale,
                           const HighsMatrixSlice<RowStorageFormat>& eqRowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : eqRowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(EqualityRowAddition{
        origRowIndex[row], origRowIndex[addedEqRow], eqRowScale});
    reductionValues.push(rowValues);
    reductionAdded(ReductionType::kEqualityRowAddition);
  }

  template <typename RowStorageFormat>
  void equalityRowAdditions(HighsInt addedEqRow,
                            const HighsMatrixSlice<RowStorageFormat>& eqRowVec,
                            const std::vector<Nonzero>& targetRows) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : eqRowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(EqualityRowAdditions{origRowIndex[addedEqRow]});
    reductionValues.push(rowValues);
    // On entry, targetRows are in the presolved space, but in
    // postsolve they need to be in the original space. For
    // equalityRowAddition the transformation is simply
    // origRowIndex[row], but for equalityRowAdditions all indices
    // must be transformed in the std::vector<Nonzero> targetRows
    std::vector<Nonzero> orig_targetRows = targetRows;
    for (auto& targetRow : orig_targetRows)
      targetRow.index = origRowIndex[targetRow.index];
    reductionValues.push(orig_targetRows);
    reductionAdded(ReductionType::kEqualityRowAdditions);
  }

  void singletonRow(HighsInt row, HighsInt col, double coef,
                    bool tightenedColLower, bool tightenedColUpper) {
    reductionValues.push(SingletonRow{coef, origRowIndex[row],
                                      origColIndex[col], tightenedColLower,
                                      tightenedColUpper});
    reductionAdded(ReductionType::kSingletonRow);
  }

  template <typename ColStorageFormat>
  void fixedColAtLower(HighsInt col, double fixValue, double colCost,
                       const HighsMatrixSlice<ColStorageFormat>& colVec) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kLower});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kFixedCol);
  }

  template <typename ColStorageFormat>
  void fixedColAtUpper(HighsInt col, double fixValue, double colCost,
                       const HighsMatrixSlice<ColStorageFormat>& colVec) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kUpper});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kFixedCol);
  }

  template <typename ColStorageFormat>
  void fixedColAtZero(HighsInt col, double colCost,
                      const HighsMatrixSlice<ColStorageFormat>& colVec) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(
        FixedCol{0.0, colCost, origColIndex[col], HighsBasisStatus::kZero});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kFixedCol);
  }

  template <typename ColStorageFormat>
  void removedFixedCol(HighsInt col, double fixValue, double colCost,
                       const HighsMatrixSlice<ColStorageFormat>& colVec) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kNonbasic});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kFixedCol);
  }

  void removedModelFixedCol(HighsInt col, double fixValue, double colCost,
                            HighsInt col_nnz, HighsInt* index, double* value) {
    assert(std::isfinite(fixValue));
    colValues.clear();
    for (HighsInt iEl = 0; iEl < col_nnz; iEl++)
      colValues.emplace_back(origRowIndex[index[iEl]], value[iEl]);

    reductionValues.push(FixedCol{fixValue, colCost, origColIndex[col],
                                  HighsBasisStatus::kNonbasic});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kFixedCol);
  }

  void redundantRow(HighsInt row) {
    reductionValues.push(RedundantRow{origRowIndex[row]});
    reductionAdded(ReductionType::kRedundantRow);
  }

  template <typename RowStorageFormat>
  void impliedEquation(HighsInt row, bool atLower,
                       const HighsMatrixSlice<RowStorageFormat>& rowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());
    reductionValues.push(ImpliedEquation{origRowIndex[row], atLower});
    reductionValues.push(rowValues);
    reductionAdded(ReductionType::kImpliedEquation);
  }

  template <typename RowStorageFormat>
  void forcingRow(HighsInt row,
                  const HighsMatrixSlice<RowStorageFormat>& rowVec, double side,
                  RowType rowType) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(ForcingRow{side, origRowIndex[row], rowType});
    reductionValues.push(rowValues);
    reductionAdded(ReductionType::kForcingRow);
  }

  template <typename ColStorageFormat>
  void forcingColumn(HighsInt col,
                     const HighsMatrixSlice<ColStorageFormat>& colVec,
                     double cost, double boundVal, bool atInfiniteUpper,
                     bool colIntegral) {
    colValues.clear();
    for (const HighsSliceNonzero& colVal : colVec)
      colValues.emplace_back(origRowIndex[colVal.index()], colVal.value());

    reductionValues.push(ForcingColumn{cost, boundVal, origColIndex[col],
                                       atInfiniteUpper, colIntegral});
    reductionValues.push(colValues);
    reductionAdded(ReductionType::kForcingColumn);
  }

  template <typename RowStorageFormat>
  void forcingColumnRemovedRow(
      HighsInt forcingCol, HighsInt row, double rhs,
      const HighsMatrixSlice<RowStorageFormat>& rowVec) {
    rowValues.clear();
    for (const HighsSliceNonzero& rowVal : rowVec)
      if (rowVal.index() != forcingCol)
        rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value());

    reductionValues.push(ForcingColumnRemovedRow{rhs, origRowIndex[row]});
    reductionValues.push(rowValues);
    reductionAdded(ReductionType::kForcingColumnRemovedRow);
  }

  void duplicateRow(HighsInt row, bool rowUpperTightened,
                    bool rowLowerTightened, HighsInt duplicateRow,
                    double duplicateRowScale) {
    // The surviving row must not be a cut absorbing a non-cut
    // constraint. detectParallelRowsAndCols ensures this by
    // iterating non-cut rows before cut rows.
    assert(!isCutRow(row) || isCutRow(duplicateRow));
    reductionValues.push(
        DuplicateRow{duplicateRowScale, origRowIndex[duplicateRow],
                     origRowIndex[row], rowLowerTightened, rowUpperTightened});
    reductionAdded(ReductionType::kDuplicateRow);
  }

  bool duplicateColumn(double colScale, double colLower, double colUpper,
                       double duplicateColLower, double duplicateColUpper,
                       HighsInt col, HighsInt duplicateCol, bool colIntegral,
                       bool duplicateColIntegral,
                       const double ok_merge_tolerance) {
    const HighsInt origCol = origColIndex[col];
    const HighsInt origDuplicateCol = origColIndex[duplicateCol];
    DuplicateColumn debug_values = {
        colScale,          colLower,          colUpper,
        duplicateColLower, duplicateColUpper, origCol,
        origDuplicateCol,  colIntegral,       duplicateColIntegral};
    const bool ok_merge = debug_values.okMerge(ok_merge_tolerance);
    const bool prevent_illegal_merge = true;
    if (!ok_merge && prevent_illegal_merge) return false;
    reductionValues.push(debug_values);
    //    reductionValues.push(DuplicateColumn{
    //        colScale, colLower, colUpper, duplicateColLower,
    //        duplicateColUpper, origCol, origDuplicateCol, colIntegral,
    //        duplicateColIntegral});

    reductionAdded(ReductionType::kDuplicateColumn);

    // mark columns as not linearly transformable
    linearlyTransformable[origCol] = false;
    linearlyTransformable[origDuplicateCol] = false;
    return true;
  }

  std::vector<double> getReducedPrimalSolution(
      const std::vector<double>& origPrimalSolution) {
    std::vector<double> reducedSolution = origPrimalSolution;
    reducedSolution.resize(nextColIndex, 0.0);

    for (const std::pair<ReductionType, size_t>& primalColTransformation :
         reductions) {
      switch (primalColTransformation.first) {
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn duplicateColReduction;
          reductionValues.setPosition(primalColTransformation.second);
          reductionValues.pop(duplicateColReduction);
          duplicateColReduction.transformToPresolvedSpace(reducedSolution);
          break;
        }
        case ReductionType::kLinearTransform: {
          reductionValues.setPosition(primalColTransformation.second);
          LinearTransform linearTransform;
          reductionValues.pop(linearTransform);
          linearTransform.transformToPresolvedSpace(reducedSolution);
          break;
        }
        default:
          continue;
      }
    }

    size_t reducedNumCol = origColIndex.size();
    for (size_t i = 0; i < reducedNumCol; ++i)
      reducedSolution[i] = reducedSolution[origColIndex[i]];

    reducedSolution.resize(reducedNumCol);
    return reducedSolution;
  }

  bool isColLinearlyTransformable(HighsInt col) const {
    assert(col >= 0);
    assert(static_cast<size_t>(col) < origColIndex.size());
    return linearlyTransformable[origColIndex[col]];
  }

  template <typename T>
  void undoIterateBackwards(std::vector<T>& values,
                            const std::vector<HighsInt>& index,
                            HighsInt origSize, T zero) {
    values.resize(origSize);
#ifdef DEBUG_EXTRA
    // Fill vector with NaN for debugging purposes
    std::vector<T> valuesNew;
    valuesNew.resize(origSize, std::numeric_limits<T>::signaling_NaN());
    for (size_t i = index.size(); i > 0; --i) {
      assert(static_cast<size_t>(index[i - 1]) >= i - 1);
      valuesNew[index[i - 1]] = values[i - 1];
    }
    std::copy(valuesNew.cbegin(), valuesNew.cend(), values.begin());
#else
    for (size_t i = index.size(); i < static_cast<size_t>(origSize); i++)
      values[i] = zero;
    for (size_t i = index.size(); i > 0; --i) {
      size_t to_i = static_cast<size_t>(index[i - 1]);
      assert(to_i >= i - 1);
      assert(to_i < static_cast<size_t>(origSize));
      values[to_i] = values[i - 1];
      if (to_i > i - 1) values[i - 1] = zero;
    }
#endif
  }

  /// check if vector contains NaN or Inf
  bool containsNanOrInf(const std::vector<double>& v) const {
    return std::find_if(v.cbegin(), v.cend(), [](const double& d) {
             return (std::isnan(d) || std::isinf(d));
           }) != v.cend();
  }

  /// undo presolve steps for primal dual solution and basis
  void undo(const HighsOptions& options, HighsSolution& solution,
            HighsBasis& basis, size_t numReductions = 0,
            const HighsInt report_col = -1, const bool thread_safe = false) {
    HighsDataStack reductionValuesCopy;
    std::vector<Nonzero> colValuesCopy;
    std::vector<Nonzero> rowValuesCopy;
    if (!thread_safe) {
      reductionValues.resetPosition();
    } else {
      reductionValuesCopy = reductionValues;
      reductionValuesCopy.resetPosition();
      colValuesCopy = colValues;
      rowValuesCopy = rowValues;
    }
    HighsDataStack& reductionValues_ =
        thread_safe ? reductionValuesCopy : reductionValues;
    std::vector<Nonzero>& colValues_ = thread_safe ? colValuesCopy : colValues;
    std::vector<Nonzero>& rowValues_ = thread_safe ? rowValuesCopy : rowValues;

    // Verify that undo can be performed
    assert(solution.value_valid);
    bool perform_dual_postsolve = solution.dual_valid;
    bool perform_basis_postsolve = basis.valid;

    // expand solution to original index space
    assert(nextColIndex > 0);
    undoIterateBackwards(solution.col_value, origColIndex, nextColIndex, 0.0);

    assert(nextRowIndex >= 0);
    undoIterateBackwards(solution.row_value, origRowIndex, nextRowIndex, 0.0);

    if (perform_dual_postsolve) {
      // if dual solution is given, expand dual solution and basis to original
      // index space
      undoIterateBackwards(solution.col_dual, origColIndex, nextColIndex, 0.0);

      undoIterateBackwards(solution.row_dual, origRowIndex, nextRowIndex, 0.0);
    }

    if (perform_basis_postsolve) {
      // if basis is given, expand basis status values to original index space
      undoIterateBackwards(basis.col_status, origColIndex, nextColIndex,
                           HighsBasisStatus::kNonbasic);

      undoIterateBackwards(basis.row_status, origRowIndex, nextRowIndex,
                           HighsBasisStatus::kNonbasic);
    }

    /*
    // Code for logging when debugging
    std::stringstream ss;

    // Lambda for logging the whole solution
    auto solutionLogging = [&](const std::string& message) {
      printf("\n%s\n", message.c_str());
      for (HighsInt iCol = 0; iCol < origNumCol; iCol++) {
        ss.str(std::string());
        ss << highsFormatToString("Col %9d value = %11.4g", int(iCol),
                                  solution.col_value[iCol]);
        if (perform_dual_postsolve)
          ss << highsFormatToString("; dual = %11.4g;",
                                    solution.col_dual[iCol]);
        if (perform_basis_postsolve)
          ss << highsFormatToString("; status = %s",
                                    utilBasisStatusToString(basis.col_status[iCol]).full_.c_str());
        printf("%s\n", ss.str().c_str());
      }
      for (HighsInt iRow = 0; iRow < origNumRow; iRow++) {
        ss.str(std::string());
        ss << highsFormatToString("Row %9d value = %11.4g", int(iRow),
                                  solution.row_value[iRow]);
        if (perform_dual_postsolve)
          ss << highsFormatToString("; dual = %11.4g;",
                                    solution.row_dual[iRow]);
        if (perform_basis_postsolve)
          ss << highsFormatToString("; status = %s",
                                    utilBasisStatusToString(basis.row_status[iRow]).full_.c_str());
        printf("%s\n", ss.str().c_str());
      }
    };

    // Initialise to illegal values so that initial values are logged
    double report_col_value = kHighsInf;

    // Lambda for logging the solution for a specific column whenever its value
    changes

    auto reportColLogging = [&](const HighsInt reduction) {
      assert(report_col >= 0);
      if (static_cast<size_t>(report_col) >= solution.col_value.size()) return;
      double col_value = solution.col_value[report_col];
      bool report = col_value != report_col_value;
      if (reduction >= 0) {
        ReductionType type = reductions[reduction].first;
        if (report)
          printf("After reduction %9d (type %s):", int(reduction),
                 presolveTypeToString(type).c_str());
      } else if (reduction == -1) {
        report = true;
        printf("Before undo:                        ");
      } else {
        report = true;
        printf("After last reduction:               ");
      }
      if (!report) return;
      ss.str(std::string());
      ss << highsFormatToString(" Col %7d value = %11.4g",
                                int(report_col), col_value);
      if (perform_dual_postsolve)
        ss << highsFormatToString(", dual = %11.4g",
                                  solution.col_dual[report_col]);
      if (perform_basis_postsolve)
        ss << highsFormatToString(" status = %s",
                                  utilBasisStatusToString(basis.col_status[report_col]).full_.c_str());
      printf("%s\n", ss.str().c_str());
      report_col_value = col_value;
    };

    if (report_col >= 0) reportColLogging(-1);

    size_t check_reduction = kHighsSize_tInf;
    if (reductions.size() == check_reduction)
      solutionLogging("After solving presolved LP");
    */

    // now undo the changes
    for (size_t i = reductions.size(); i > numReductions; --i) {
      if (report_col >= 0)
        printf("Before  reduction %2d (type %2d): col_value[%2d] = %g\n",
               int(i - 1), int(reductions[i - 1].first), int(report_col),
               solution.col_value[report_col]);
      /*
      if (i - 1 == check_reduction) {
        printf("Checking reduction %d\n", int(check_reduction));
        solutionLogging("In reductions loop");
      }
      */
      switch (reductions[i - 1].first) {
        case ReductionType::kLinearTransform: {
          LinearTransform reduction;
          reductionValues_.pop(reduction);
          reduction.undo(options, solution);
          break;
        }
        case ReductionType::kFreeColSubstitution: {
          FreeColSubstitution reduction;
          reductionValues_.pop(colValues_);
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, rowValues_, colValues_, solution, basis);
          break;
        }
        case ReductionType::kDoubletonEquation: {
          DoubletonEquation reduction;
          reductionValues_.pop(colValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, colValues_, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAddition: {
          EqualityRowAddition reduction;
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, rowValues_, solution, basis);
          break;
        }
        case ReductionType::kEqualityRowAdditions: {
          EqualityRowAdditions reduction;
          reductionValues_.pop(colValues_);
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, rowValues_, colValues_, solution, basis);
          break;
        }
        case ReductionType::kSingletonRow: {
          SingletonRow reduction;
          reductionValues_.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kFixedCol: {
          FixedCol reduction;
          reductionValues_.pop(colValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, colValues_, solution, basis);
          break;
        }
        case ReductionType::kRedundantRow: {
          RedundantRow reduction;
          reductionValues_.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kImpliedEquation: {
          ImpliedEquation reduction;
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(rowValues_, solution);
          break;
        }
        case ReductionType::kForcingRow: {
          ForcingRow reduction;
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, rowValues_, solution, basis);
          break;
        }
        case ReductionType::kForcingColumn: {
          ForcingColumn reduction;
          reductionValues_.pop(colValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, colValues_, solution, basis);
          break;
        }
        case ReductionType::kForcingColumnRemovedRow: {
          ForcingColumnRemovedRow reduction;
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, rowValues_, solution, basis);
          break;
        }
        case ReductionType::kDuplicateRow: {
          DuplicateRow reduction;
          reductionValues_.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kDuplicateColumn: {
          DuplicateColumn reduction;
          reductionValues_.pop(reduction);
          reduction.undo(options, solution, basis);
          break;
        }
        case ReductionType::kZeroObjSingletonContinuousCol: {
          ZeroObjSingletonContinuousCol reduction;
          reductionValues_.pop(rowValues_);
          reductionValues_.pop(reduction);
          reduction.undo(options, rowValues_, solution, basis);
          break;
        }
        case ReductionType::kFourierMotzkinBlock: {
          auto steps = popFourierMotzkinBlock(reductionValues_);
          undoFourierMotzkinBlock(steps, options, solution, basis);
          break;
        }
        case ReductionType::kFourierMotzkinObjCol: {
          std::vector<Nonzero> costEntries;
          reductionValues_.pop(costEntries);
          FourierMotzkinObjCol reduction;
          reductionValues_.pop(reduction);
          reduction.undo(costEntries, solution);
          break;
        }
        default:
          printf("Reduction case %d not handled\n",
                 int(reductions[i - 1].first));
          if (kAllowDeveloperAssert) assert(1 == 0);
      }
      //      if (report_col >= 0) reportColLogging(i - 1);
    }
    //    if (report_col >= 0) reportColLogging(-2);

    solution.col_value.resize(origNumCol);
    if (perform_dual_postsolve) solution.col_dual.resize(origNumCol);
    if (perform_basis_postsolve) basis.col_status.resize(origNumCol);

    solution.row_value.resize(origNumRow);
    if (perform_dual_postsolve) solution.row_dual.resize(origNumRow);
    if (perform_basis_postsolve) basis.row_status.resize(origNumRow);

#ifdef DEBUG_EXTRA
    // solution should not contain NaN or Inf
    assert(!containsNanOrInf(solution.col_value));
    // row values are not determined by postsolve
    // assert(!containsNanOrInf(solution.row_value));
    assert(!containsNanOrInf(solution.col_dual));
    assert(!containsNanOrInf(solution.row_dual));
#endif
  }

  /// undo presolve steps for primal solution
  void undoPrimal(const HighsOptions& options, HighsSolution& solution,
                  const HighsInt report_col = -1,
                  const bool thread_safe = false) {
    // Call to reductionValues.resetPosition(); seems unnecessary as
    // it's the first thing done in undo
    if (!thread_safe) reductionValues.resetPosition();
    HighsBasis basis;
    basis.valid = false;
    solution.dual_valid = false;
    undo(options, solution, basis, 0, report_col, thread_safe);
  }

  size_t numReductions() const { return reductions.size(); }
};

}  // namespace presolve

#endif
