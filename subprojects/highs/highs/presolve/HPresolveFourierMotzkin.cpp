/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::fourierMotzkin(
    HighsPostsolveStack& postsolve_stack, HighsInt& numColsEliminated) {
  assert(this->allow_rule_[kPresolveRuleFourierMotzkin]);
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleFourierMotzkin);

  using FmeRow = HighsPostsolveStack::FmeRowData<HighsTripletTreeSlicePreOrder>;
  using FmeAncestryEntry = HighsPostsolveStack::FmeAncestryEntry;
  using FmeBlockStep = HighsPostsolveStack::FmeBlockStep;

  // max. absolute coefficient
  const double maxCoef = 1e3;

  // max. number of consecutive failures (while trying to build the heap)
  const HighsInt maxNumFails = 100;
  // max. size of the heap
  const HighsInt maxHeapSize = 10000;

  // sentinel row indices for variable bounds and objective row
  const HighsInt kUpperBoundRow = -2;
  const HighsInt kLowerBoundRow = -3;
  const HighsInt kObjectiveRow = -4;

  // structs
  struct Heap {
    struct Entry {
      HighsInt col;
      int64_t neRed;
      int64_t mrRed;
    };

    std::vector<Entry> entries;
    std::vector<HighsInt> pos;

    bool empty() const { return entries.empty(); }
    HighsInt size() const { return static_cast<HighsInt>(entries.size()); }
    HighsInt top() const { return entries[0].col; }
    bool contains(HighsInt col) const { return pos[col] != -1; }

    void reset(HighsInt numCol, HighsInt reserveSize) {
      entries.clear();
      entries.reserve(reserveSize);
      pos.assign(numCol, -1);
    }

    void push(HighsInt col, int64_t neRed, int64_t mrRed) {
      pos[col] = size();
      entries.push_back({col, neRed, mrRed});
    }

    void insert(HighsInt col, int64_t neRed, int64_t mrRed) {
      push(col, neRed, mrRed);
      siftUp(pos[col]);
    }

    void remove(HighsInt col) {
      HighsInt p = pos[col];
      if (p == -1) return;
      swap(p, size() - 1);
      pos[col] = -1;
      entries.pop_back();
      siftUp(p);
      siftDown(p);
    }

    void update(HighsInt col, int64_t neRed, int64_t mrRed) {
      HighsInt p = pos[col];
      if (p == -1) return;
      entries[p].neRed = neRed;
      entries[p].mrRed = mrRed;
      siftUp(p);
      siftDown(p);
    }

    void heapify() {
      for (HighsInt i = size() / 2 - 1; i >= 0; --i) siftDown(i);
    }

   private:
    bool better(HighsInt i, HighsInt j) const {
      if (entries[i].neRed != entries[j].neRed)
        return entries[i].neRed > entries[j].neRed;
      return entries[i].mrRed > entries[j].mrRed;
    }

    void swap(HighsInt i, HighsInt j) {
      if (i == j) return;
      std::swap(entries[i], entries[j]);
      pos[entries[i].col] = i;
      pos[entries[j].col] = j;
    }

    void siftUp(HighsInt i) {
      if (i >= size()) return;
      while (i > 0) {
        HighsInt parent = (i - 1) / 2;
        if (!better(i, parent)) break;
        swap(i, parent);
        i = parent;
      }
    }

    void siftDown(HighsInt i) {
      HighsInt n = size();
      if (i >= n) return;
      while (true) {
        HighsInt best = i;
        HighsInt left = 2 * i + 1;
        HighsInt right = 2 * i + 2;
        if (left < n && better(left, best)) best = left;
        if (right < n && better(right, best)) best = right;
        if (best == i) break;
        swap(i, best);
        i = best;
      }
    }
  };

  struct newRowEntry {
    HighsInt col;
    HighsCDouble val;
  };

  struct newRow {
    std::vector<newRowEntry> entries;
    double upper;
    HighsInt plusIndex;
    HighsInt minusIndex;
    double plusScale;
    double minusScale;
  };

  struct NewRowOrigin {
    HighsInt plusRow;
    HighsInt minusRow;
    double plusScale;
    double minusScale;
  };

  auto finalise = [&]() {
    analysis_.logging_on_ = logging_on;
    if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleFourierMotzkin);
    return checkLimits(postsolve_stack);
  };

  auto acceptCoef = [&](double val) {
    double absval = std::abs(val);
    return absval == 0.0 || (absval >= 1.0 / maxCoef && absval <= maxCoef);
  };

  auto isCandidate = [&](HighsInt col) {
    if (colDeleted[col]) return false;
    if (colsize[col] == 0) return false;
    if (col == model->fme_obj_col_) return false;
    if (model->integrality_[col] != HighsVarType::kContinuous) return false;
    if (!acceptCoef(model->col_cost_[col])) return false;
    if (!options->presolve_fm_obj_reformulation && model->col_cost_[col] != 0.0)
      return false;
    for (const auto& nz : getColumnVector(col))
      if (isEquation(nz.index()) || !acceptCoef(nz.value())) return false;
    return true;
  };

  auto computeCandidates = [&](std::vector<HighsInt>& candidates) {
    candidates.clear();
    for (HighsInt col = 0; col < model->num_col_; col++)
      if (isCandidate(col)) candidates.push_back(col);
    pdqsort(candidates.begin(), candidates.end(),
            [&](HighsInt a, HighsInt b) { return colsize[a] < colsize[b]; });
    return !candidates.empty();
  };

  auto checkRows = [&](HighsInt col, const std::vector<HighsInt>& objRowCols,
                       std::vector<HighsInt>& iPlus,
                       std::vector<HighsInt>& iMinus, int64_t& nePlus,
                       int64_t& neMinus) {
    nePlus = 0;
    neMinus = 0;
    iPlus.clear();
    iMinus.clear();
    for (const auto& nz : getColumnVector(col)) {
      HighsInt row = nz.index();
      if (rowDeleted[row]) continue;

      if (isRanged(row)) {
        iPlus.push_back(row);
        nePlus += rowsize[row];
        iMinus.push_back(row);
        neMinus += rowsize[row];
      } else {
        HighsInt direction = model->row_upper_[row] != kHighsInf ? 1 : -1;

        if (direction * nz.value() > 0) {
          iPlus.push_back(row);
          nePlus += rowsize[row];
        } else {
          iMinus.push_back(row);
          neMinus += rowsize[row];
        }
      }
    }

    // include finite variable bounds as singleton rows
    if (model->col_upper_[col] != kHighsInf) {
      iPlus.push_back(kUpperBoundRow);
      nePlus += 1;
    }
    if (model->col_lower_[col] != -kHighsInf) {
      iMinus.push_back(kLowerBoundRow);
      neMinus += 1;
    }

    // simulate the objective constraint row for candidates with nonzero
    // cost when the reformulation has not yet been performed
    if (!objRowCols.empty() && model->col_cost_[col] != 0.0) {
      int64_t objRowSize = static_cast<int64_t>(objRowCols.size());
      if (model->col_cost_[col] > 0.0) {
        iPlus.push_back(kObjectiveRow);
        nePlus += objRowSize;
      } else {
        iMinus.push_back(kObjectiveRow);
        neMinus += objRowSize;
      }
    }
  };

  auto collectAffectedCols = [&](HighsInt col, const std::vector<HighsInt>& set,
                                 const std::vector<HighsInt>& objRowCols,
                                 std::vector<HighsInt>& mark,
                                 std::vector<HighsInt>& otherMark,
                                 std::vector<HighsInt>& affectedCols) {
    for (HighsInt row : set) {
      if (row == kObjectiveRow) {
        for (HighsInt k : objRowCols) {
          if (k == col) continue;
          if (mark[k] == 0 && otherMark[k] == 0) affectedCols.push_back(k);
          mark[k]++;
        }
      } else {
        if (row < 0) continue;
        for (const auto& nz : getRowVector(row)) {
          HighsInt k = nz.index();
          if (k == col) continue;
          if (mark[k] == 0 && otherMark[k] == 0) affectedCols.push_back(k);
          mark[k]++;
        }
      }
    }
  };

  auto checkNonZeros = [&](HighsInt col,
                           const std::vector<HighsInt>& objRowCols,
                           std::vector<HighsInt>& iPlus,
                           std::vector<HighsInt>& iMinus,
                           std::vector<HighsInt>& pPlus,
                           std::vector<HighsInt>& pMinus,
                           std::vector<HighsInt>& affectedCols, int64_t& neRed,
                           int64_t& mrRed) {
    // initialise
    neRed = 0;
    mrRed = 0;

    // check rows
    int64_t nePlus;
    int64_t neMinus;
    checkRows(col, objRowCols, iPlus, iMinus, nePlus, neMinus);

    if (iPlus.size() == 0 || iMinus.size() == 0) {
      // other presolve reductions may handle this case (e.g., implied free
      // column substitution)
      iPlus.clear();
      iMinus.clear();
      return false;
    }

    // take into account other variables present in the rows
    collectAffectedCols(col, iPlus, objRowCols, pPlus, pMinus, affectedCols);
    collectAffectedCols(col, iMinus, objRowCols, pMinus, pPlus, affectedCols);

    // compute correction term
    int64_t correction = 0;
    for (HighsInt k : affectedCols) {
      correction += static_cast<int64_t>(pPlus[k]) * pMinus[k];
      pPlus[k] = 0;
      pMinus[k] = 0;
    }

    int64_t mPlus = static_cast<int64_t>(iPlus.size());
    int64_t mMinus = static_cast<int64_t>(iMinus.size());
    int64_t neOld = nePlus + neMinus;
    // note that we subtract the entries for column 'col' since these are
    // eliminated
    int64_t neNew =
        mPlus * (neMinus - mMinus) + mMinus * (nePlus - mPlus) - correction;
    neRed = neOld - neNew;
    mrRed = mPlus + mMinus - mPlus * mMinus;
    return true;
  };

  auto checkNewRow = [&](const newRow& nr, bool& isRedundant) {
    HighsCDouble impliedLower = 0;
    HighsCDouble impliedUpper = 0;
    bool lowerFinite = true;
    bool upperFinite = true;
    isRedundant = false;
    for (const auto& e : nr.entries) {
      double lb = model->col_lower_[e.col];
      double ub = model->col_upper_[e.col];
      if (e.val > 0) {
        lowerFinite = lowerFinite && lb != -kHighsInf;
        if (lowerFinite) impliedLower += e.val * lb;
        upperFinite = upperFinite && ub != kHighsInf;
        if (upperFinite) impliedUpper += e.val * ub;
      } else {
        lowerFinite = lowerFinite && ub != kHighsInf;
        if (lowerFinite) impliedLower += e.val * ub;
        upperFinite = upperFinite && lb != -kHighsInf;
        if (upperFinite) impliedUpper += e.val * lb;
      }
      if (!lowerFinite && !upperFinite) return Result::kOk;
    }

    double lower = lowerFinite ? static_cast<double>(impliedLower) : -kHighsInf;
    double upper = upperFinite ? static_cast<double>(impliedUpper) : kHighsInf;

    // check for infeasibility
    if (lower > nr.upper + primal_feastol) return Result::kPrimalInfeasible;

    // check for redundancy
    isRedundant = upper <= nr.upper + primal_feastol;

    return Result::kOk;
  };

  auto getRowData = [&](HighsInt row, HighsInt col, HighsInt multiplier,
                        double& absCoef, HighsInt& direction, double& bound) {
    if (row < 0) {
      // artificial lower / upper bound row
      direction = 1;
      absCoef = 1.0;
      bound = multiplier > 0 ? model->col_upper_[col] : -model->col_lower_[col];
    } else {
      HighsInt pPos = findNonzero(row, col);
      assert(pPos != -1);
      direction = multiplier * Avalue[pPos] > 0 ? HighsInt{1} : HighsInt{-1};
      absCoef = std::abs(Avalue[pPos]);
      bound = direction > 0 ? model->row_upper_[row] : -model->row_lower_[row];
    }
  };

  auto collectRowEntries = [&](HighsInt row, HighsInt col, double scale,
                               std::vector<newRowEntry>& newRowEntries,
                               std::vector<HighsInt>& newRowMark) {
    if (row < 0) return;
    for (const auto& nz : getRowVector(row)) {
      if (nz.index() == col) continue;
      double val = scale * nz.value();
      if (newRowMark[nz.index()] == -1) {
        newRowMark[nz.index()] = static_cast<HighsInt>(newRowEntries.size());
        newRowEntries.push_back({nz.index(), val});
      } else {
        newRowEntries[newRowMark[nz.index()]].val += val;
      }
    }
  };

  auto isReduction = [](int64_t neRed, int64_t mrRed) {
    return neRed > 0 || (neRed == 0 && mrRed > 0);
  };

  auto insertOriginals =
      [&](std::set<HighsInt>& rows,
          const std::unordered_map<HighsInt, std::set<HighsInt>>& originals,
          HighsInt row, HighsInt col) {
        if (row == kUpperBoundRow)
          rows.insert(-(2 * col + 1));
        else if (row == kLowerBoundRow)
          rows.insert(-(2 * col + 2));
        else {
          auto it = originals.find(row);
          if (it != originals.end())
            rows.insert(it->second.begin(), it->second.end());
          else
            rows.insert(row);
        }
      };

  auto mergeOriginals =
      [&](std::set<HighsInt>& rows,
          const std::unordered_map<HighsInt, std::set<HighsInt>>& originals,
          HighsInt plusRow, HighsInt minusRow, HighsInt col) {
        rows.clear();
        insertOriginals(rows, originals, plusRow, col);
        insertOriginals(rows, originals, minusRow, col);
      };

  auto cernikovRedundant =
      [&](std::set<HighsInt>& rows,
          const std::unordered_map<HighsInt, std::set<HighsInt>>& originals,
          HighsInt plusRow, HighsInt minusRow, HighsInt col,
          HighsInt numColsElim) {
        mergeOriginals(rows, originals, plusRow, minusRow, col);
        return static_cast<HighsInt>(rows.size()) > numColsElim + 2;
      };

  // reformulate objective as a constraint: min c^T x + offset becomes
  // min z with c^T x - z <= -offset. this allows FME to eliminate
  // continuous columns with nonzero cost.
  auto reformulateObjective = [&]() {
    if (model->fme_obj_col_ != -1) {
      assert(!colDeleted[model->fme_obj_col_]);
      return;
    }

    HighsInt zCol = model->num_col_;
    model->num_col_++;
    model->a_matrix_.num_col_++;

    // extend model vectors
    model->col_cost_.push_back(1.0);
    model->col_lower_.push_back(-kHighsInf);
    model->col_upper_.push_back(kHighsInf);
    model->integrality_.push_back(HighsVarType::kContinuous);
    model->a_matrix_.start_.push_back(model->a_matrix_.start_.back());
    if (model->col_names_.size() > 0) model->col_names_.push_back("fme_obj_z");

    // extend presolve vectors
    colhead.push_back(-1);
    colsize.push_back(0);
    colDeleted.push_back(0);
    implColLower.push_back(-kHighsInf);
    implColUpper.push_back(kHighsInf);
    colLowerSource.push_back(-1);
    colUpperSource.push_back(-1);
    implRowDualSourceByCol.push_back({});
    changedColFlag.push_back(1);
    numProbes.push_back(0);

    // update implied bound structures (pointers may be invalidated by
    // reallocation of column vectors above)
    impliedRowBounds.setBoundArrays(
        model->col_lower_.data(), model->col_upper_.data(), implColLower.data(),
        implColUpper.data(), colLowerSource.data(), colUpperSource.data());
    impliedDualRowBounds.setNumSums(model->num_col_);

    // register in postsolve stack
    postsolve_stack.appendColToModel();

    // build the objective constraint row: c^T x - z <= -offset
    double offset = model->offset_;
    std::vector<HighsInt> objIndices;
    std::vector<double> objValues;
    for (HighsInt j = 0; j < zCol; ++j) {
      if (!colDeleted[j] && model->col_cost_[j] != 0.0) {
        objIndices.push_back(j);
        objValues.push_back(model->col_cost_[j]);
      }
    }
    objIndices.push_back(zCol);
    objValues.push_back(-1.0);

    // zero out original costs and offset
    for (HighsInt j = 0; j < zCol; ++j) model->col_cost_[j] = 0.0;
    model->offset_ = 0.0;

    // add the constraint row to the matrix
    addToMatrix(postsolve_stack, -kHighsInf, -offset, objIndices, objValues);

    // register reduction so getReducedPrimalSolution can compute z
    std::vector<HighsPostsolveStack::Nonzero> costEntries;
    for (size_t k = 0; k < objIndices.size(); k++)
      costEntries.emplace_back(objIndices[k], objValues[k]);
    postsolve_stack.fourierMotzkinObjCol(zCol, offset, costEntries);

    model->fme_obj_col_ = zCol;

    shrinkProblem(postsolve_stack);
  };

  auto collectCandidatesAndBuildHeap =
      [&](std::vector<HighsInt>& candidates, Heap& heap,
          std::vector<HighsInt>& iPlus, std::vector<HighsInt>& iMinus,
          std::vector<HighsInt>& pPlus, std::vector<HighsInt>& pMinus,
          std::vector<HighsInt>& affectedCols,
          const std::vector<HighsInt>& objRowCols) {
        // compute candidates
        if (!computeCandidates(candidates)) return false;
        // set up data structures for heap
        heap.reset(model->num_col_, static_cast<HighsInt>(candidates.size()));
        pPlus.assign(model->num_col_, 0);
        pMinus.assign(model->num_col_, 0);
        iPlus.reserve(model->num_row_);
        iMinus.reserve(model->num_row_);
        affectedCols.reserve(model->num_col_);
        // inspect candidates (with limits)
        HighsInt numFails = 0;
        for (HighsInt col : candidates) {
          int64_t neRed = 0;
          int64_t mrRed = 0;
          bool elimCandidate =
              checkNonZeros(col, objRowCols, iPlus, iMinus, pPlus, pMinus,
                            affectedCols, neRed, mrRed);
          affectedCols.clear();
          if (!elimCandidate || !isReduction(neRed, mrRed)) {
            // count number of failures
            if (++numFails > maxNumFails) break;
            continue;
          }
          // add to heap
          numFails = 0;
          heap.push(col, neRed, mrRed);
          if (heap.size() >= maxHeapSize) break;
        }
        if (heap.empty()) return false;
        // heapify
        heap.heapify();
        return true;
      };

  // find index of a row within a list
  auto findRowIndex = [](HighsInt row,
                         const std::vector<FmeRow>& rows) -> HighsInt {
    for (HighsInt i = 0; i < static_cast<HighsInt>(rows.size()); ++i)
      if (rows[i].header.row == row) return i;
    return -1;
  };

  auto collectRows = [&](const std::vector<HighsInt>& rows) {
    std::vector<FmeRow> result;
    for (HighsInt r : rows) {
      if (r < 0) continue;
      result.push_back(
          {{r, model->row_lower_[r], model->row_upper_[r]}, getRowVector(r)});
    }
    return result;
  };

  auto inheritAncestry =
      [&](std::unordered_map<HighsInt, std::vector<FmeAncestryEntry>>&
              rowAncestry,
          HighsInt newModelRow, HighsInt parentRow, HighsInt parentRowIndex,
          HighsInt stepIndex, double scale, bool isMinus) {
        if (parentRow < 0) return;
        auto it = rowAncestry.find(parentRow);
        if (it != rowAncestry.end()) {
          for (const auto& a : it->second)
            rowAncestry[newModelRow].push_back(
                {a.step, a.parentRowIndex, a.scale * scale, a.isMinus});
        }
        if (parentRowIndex >= 0)
          rowAncestry[newModelRow].push_back(
              {stepIndex, parentRowIndex, scale, isMinus});
      };

  auto printLog = [&](HighsInt colsRemoved, HighsInt rowsRemoved,
                      HighsInt rowsAdded) {
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "Fourier-Motzkin (%s objective reformulation) added "
                "%" HIGHSINT_FORMAT " rows and eliminated %" HIGHSINT_FORMAT
                " rows and %" HIGHSINT_FORMAT " columns\n",
                options->presolve_fm_obj_reformulation ? "with" : "without",
                rowsAdded, rowsRemoved, colsRemoved);
  };

  // workspace vectors
  std::vector<HighsInt> candidates;
  std::vector<HighsInt> iPlus;
  std::vector<HighsInt> iMinus;
  std::vector<HighsInt> pPlus;
  std::vector<HighsInt> pMinus;
  std::vector<HighsInt> affectedCols;

  // indexed max-heap
  Heap heap;

  // precompute the objective row: columns with nonzero cost
  // used to simulate the objective constraint in checkRows before
  // reformulation actually happens
  std::vector<HighsInt> objRowCols;
  if (options->presolve_fm_obj_reformulation && model->fme_obj_col_ == -1) {
    for (HighsInt j = 0; j < model->num_col_; ++j) {
      if (!colDeleted[j] && model->col_cost_[j] != 0.0) objRowCols.push_back(j);
    }
  }

  // compute candidates and build initial heap
  if (!collectCandidatesAndBuildHeap(candidates, heap, iPlus, iMinus, pPlus,
                                     pMinus, affectedCols, objRowCols))
    return finalise();

  // vectors for computing new row entries
  std::vector<newRowEntry> newRowEntries;
  std::vector<HighsInt> newRowMark(model->num_col_, -1);

  // vector for storing new rows
  std::vector<newRow> newRows;

  // workspace for filtering new rows
  std::vector<double> rowLower;
  std::vector<double> rowUpper;
  std::vector<std::vector<HighsInt>> rowIndices;
  std::vector<std::vector<double>> rowValues;
  std::vector<NewRowOrigin> newRowOrigins;

  // vector for saving affected candidates
  std::vector<HighsInt> saveAffectedCols;

  // counters for numbers of eliminations
  numColsEliminated = 0;
  HighsInt numColsEliminatedBlock = 0;
  HighsInt numRowsEliminated = 0;
  HighsInt numRowsAdded = 0;

  // FM block data for postsolve
  std::vector<FmeBlockStep> blockSteps;

  // surviving row to its ancestry (which parent rows it descends from)
  std::unordered_map<HighsInt, std::vector<FmeAncestryEntry>> rowAncestry;

  // distinct original parent rows for each derived row (Cernikov check)
  std::unordered_map<HighsInt, std::set<HighsInt>> rowOriginals;
  std::set<HighsInt> mergedOriginals;

  // main loop: eliminate variables from heap
  while (!heap.empty()) {
    HighsInt col = heap.top();
    heap.remove(col);

    // if this candidate has nonzero cost and objective has not yet been
    // reformulated, perform the reformulation now and rebuild the heap
    if (model->fme_obj_col_ == -1 && model->col_cost_[col] != 0.0) {
      // finalise any in-progress FM block before reformulating, since
      // reformulateObjective pushes other reductions onto the data stack
      if (!blockSteps.empty()) {
        postsolve_stack.fourierMotzkinBlockFinalise(blockSteps, rowAncestry);
        printLog(numColsEliminatedBlock, numRowsEliminated, numRowsAdded);
        blockSteps.clear();
        rowAncestry.clear();
        rowOriginals.clear();
        numColsEliminatedBlock = 0;
        numRowsEliminated = 0;
        numRowsAdded = 0;
      }
      // reformulate objective
      reformulateObjective();
      // clear vector for objective and resize marker
      objRowCols.clear();
      newRowMark.resize(model->num_col_, -1);
      // re-compute candidates and re-build heap
      if (!collectCandidatesAndBuildHeap(candidates, heap, iPlus, iMinus, pPlus,
                                         pMinus, affectedCols, objRowCols))
        return finalise();
      continue;
    }

    // compute affected columns
    int64_t neRed = 0;
    int64_t mrRed = 0;
    bool elimCandidate = checkNonZeros(col, objRowCols, iPlus, iMinus, pPlus,
                                       pMinus, affectedCols, neRed, mrRed);

    // heap data should be up-to-date
    assert(elimCandidate && isReduction(neRed, mrRed));

    HighsInt stepIdx = static_cast<HighsInt>(blockSteps.size());

    // perform elimination: generate new rows
    newRows.clear();
    for (HighsInt pRow : iPlus) {
      double pCoefAbs;
      double pBound;
      HighsInt pDirection;
      getRowData(pRow, col, HighsInt{1}, pCoefAbs, pDirection, pBound);

      for (HighsInt mRow : iMinus) {
        double mCoefAbs;
        double mBound;
        HighsInt mDirection;
        getRowData(mRow, col, HighsInt{-1}, mCoefAbs, mDirection, mBound);

        // scale factor to preserve violation tolerances (see section 4.3):
        double s = (pCoefAbs * mCoefAbs) / (pCoefAbs + mCoefAbs);
        double pScale = s / pCoefAbs;
        double mScale = s / mCoefAbs;

        // collect row entries
        collectRowEntries(pRow, col, pDirection * pScale, newRowEntries,
                          newRowMark);
        collectRowEntries(mRow, col, mDirection * mScale, newRowEntries,
                          newRowMark);

        // reset marker before removing near-zeros
        for (const auto& e : newRowEntries) newRowMark[e.col] = -1;

        // remove near-zero entries
        newRowEntries.erase(
            std::remove_if(newRowEntries.begin(), newRowEntries.end(),
                           [&](const newRowEntry& e) {
                             return abs(e.val) <= options->small_matrix_value;
                           }),
            newRowEntries.end());

        // store new row
        double new_upper =
            static_cast<double>(static_cast<HighsCDouble>(pScale) * pBound +
                                static_cast<HighsCDouble>(mScale) * mBound);
        newRows.push_back({newRowEntries, new_upper, pRow, mRow,
                           pDirection * pScale, mDirection * mScale});

        // clear vector
        newRowEntries.clear();
      }
    }

    // add new rows, filtering out redundant ones
    rowLower.clear();
    rowUpper.clear();
    rowIndices.clear();
    rowValues.clear();
    newRowOrigins.clear();

    for (const auto& nr : newRows) {
      bool redundant = false;
      HPRESOLVE_CHECKED_CALL(checkNewRow(nr, redundant));
      if (redundant) continue;

      // Cernikov redundancy check
      if (cernikovRedundant(mergedOriginals, rowOriginals, nr.plusIndex,
                            nr.minusIndex, col, numColsEliminated))
        continue;

      std::vector<HighsInt> indices;
      std::vector<double> values;
      indices.reserve(nr.entries.size());
      values.reserve(nr.entries.size());
      for (const auto& e : nr.entries) {
        indices.push_back(e.col);
        values.push_back(static_cast<double>(e.val));
      }
      rowLower.push_back(-kHighsInf);
      rowUpper.push_back(nr.upper);
      rowIndices.push_back(std::move(indices));
      rowValues.push_back(std::move(values));
      newRowOrigins.push_back(
          {nr.plusIndex, nr.minusIndex, nr.plusScale, nr.minusScale});
    }

    // serialize row data for postsolve before addToMatrix invalidates slices
    std::vector<FmeRow> plusRows = collectRows(iPlus);
    std::vector<FmeRow> minusRows = collectRows(iMinus);

    // push row data for this elimination step onto the postsolve stack
    postsolve_stack.fourierMotzkinBlockPushStep(col, plusRows, minusRows);

    // save block metadata
    assert(model->col_cost_[col] == 0.0);
    blockSteps.push_back({{model->col_lower_[col], model->col_upper_[col], col,
                           static_cast<HighsInt>(plusRows.size()),
                           static_cast<HighsInt>(minusRows.size())},
                          {}});

    // add new rows to matrix
    HighsInt firstNewRow = model->num_row_;
    if (!addToMatrix(postsolve_stack, rowLower, rowUpper, rowIndices,
                     rowValues))
      return finalise();
    numRowsAdded += static_cast<HighsInt>(rowIndices.size());

    // build FmeNewRow data and ancestry for this step
    auto& stepNewRows = blockSteps.back().newRows;
    stepNewRows.reserve(newRowOrigins.size());
    for (HighsInt k = 0; k < static_cast<HighsInt>(newRowOrigins.size()); ++k) {
      HighsInt newModelRow = firstNewRow + k;
      const auto& origin = newRowOrigins[k];
      HighsInt pIdx = findRowIndex(origin.plusRow, plusRows);
      HighsInt mIdx = findRowIndex(origin.minusRow, minusRows);
      inheritAncestry(rowAncestry, newModelRow, origin.plusRow, pIdx, stepIdx,
                      origin.plusScale, false);
      inheritAncestry(rowAncestry, newModelRow, origin.minusRow, mIdx, stepIdx,
                      origin.minusScale, true);
      mergeOriginals(mergedOriginals, rowOriginals, origin.plusRow,
                     origin.minusRow, col);
      rowOriginals[newModelRow] = mergedOriginals;
      stepNewRows.push_back({newModelRow, pIdx, mIdx});
    }

    // mark column as deleted
    markColDeleted(col);
    ++numColsEliminatedBlock;
    ++numColsEliminated;

    // remove old rows containing col (skip bound rows)
    for (HighsInt rp : iPlus) {
      if (rp < 0) continue;
      rowAncestry.erase(rp);
      rowOriginals.erase(rp);
      removeRow(rp);
      ++numRowsEliminated;
    }
    for (HighsInt rm : iMinus) {
      if (rm < 0) continue;
      rowAncestry.erase(rm);
      rowOriginals.erase(rm);
      if (rowDeleted[rm]) continue;
      removeRow(rm);
      ++numRowsEliminated;
    }

    // update affected candidates in the heap
    saveAffectedCols.swap(affectedCols);
    for (HighsInt k : saveAffectedCols) {
      // check if variable is a candidate
      bool isCandidateCol = isCandidate(k);
      // skip variable if it is not on the heap and no candidate
      if (!heap.contains(k) && !isCandidateCol) continue;
      // check column non-zeros
      int64_t ne = 0;
      int64_t mr = 0;
      bool elimCandidate =
          isCandidateCol && checkNonZeros(k, objRowCols, iPlus, iMinus, pPlus,
                                          pMinus, affectedCols, ne, mr);
      affectedCols.clear();
      if (!elimCandidate || !isReduction(ne, mr)) {
        // no candidate or not beneficial -> remove from heap
        heap.remove(k);
      } else if (!heap.contains(k)) {
        // new candidate -> insert into heap
        heap.insert(k, ne, mr);
      } else {
        // update heap
        heap.update(k, ne, mr);
      }
    }
    saveAffectedCols.clear();

    if (checkLimits(postsolve_stack) != Result::kOk) break;
  }

  if (numColsEliminatedBlock > 0) {
    // finalize the FM block
    postsolve_stack.fourierMotzkinBlockFinalise(blockSteps, rowAncestry);

    // log message
    printLog(numColsEliminatedBlock, numRowsEliminated, numRowsAdded);
  }

  return finalise();
}

}  // namespace presolve
