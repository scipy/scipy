/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsMipSolverData.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::detectParallelRowsAndCols(
    HighsPostsolveStack& postsolve_stack) {
  assert(this->allow_rule_[kPresolveRuleParallelRowsAndCols]);
  const bool logging_on = analysis_.logging_on_;
  if (logging_on)
    analysis_.startPresolveRuleLog(kPresolveRuleParallelRowsAndCols);

  std::vector<std::uint64_t> rowHashes;
  std::vector<std::uint64_t> colHashes;
  std::vector<std::pair<double, HighsInt>> rowMax(rowsize.size());
  std::vector<std::pair<double, HighsInt>> colMax(colsize.size());

  HighsInt nnz = Avalue.size();
  // the initial row hashes are the numbers of non-singleton columns in the rows
  rowHashes.resize(rowsize.size());
  std::transform(rowsize.begin(), rowsize.end(), rowsizeSingleton.begin(),
                 rowHashes.begin(), std::minus<HighsInt>());
  colHashes.assign(colsize.begin(), colsize.end());

  // Step 1: Determine scales for rows and columns
  for (HighsInt i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    assert(!colDeleted[Acol[i]]);
    if (colsize[Acol[i]] == 1) {
      colMax[Acol[i]].first = Avalue[i];
      continue;
    }
    double absVal = std::abs(Avalue[i]);
    double absRowMax = std::abs(rowMax[Arow[i]].first);

    // among the largest values which are equal in tolerance
    // we use the nonzero with the smallest row/column index for the column/row
    // scale so that we ensure that duplicate rows/columns are scaled to have
    // the same sign
    if (absVal >= absRowMax - options->small_matrix_value) {
      // we are greater or equal with tolerances, check if we are either
      // strictly larger or equal with a smaller index and remember the signed
      // nonzero if one of those things is the case
      if (absVal > absRowMax + options->small_matrix_value ||
          Acol[i] < rowMax[Arow[i]].second) {
        rowMax[Arow[i]].first = Avalue[i];
        rowMax[Arow[i]].second = Acol[i];
      }
    }

    double absColMax = std::abs(colMax[Acol[i]].first);
    if (absVal >= absColMax - options->small_matrix_value) {
      if (absVal > absColMax + options->small_matrix_value ||
          Arow[i] < colMax[Acol[i]].second) {
        colMax[Acol[i]].first = Avalue[i];
        colMax[Acol[i]].second = Arow[i];
      }
    }
  }

  // Step 2: Compute hash values for rows and columns excluding singleton
  // columns
  for (HighsInt i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    assert(!rowDeleted[Arow[i]] && !colDeleted[Acol[i]]);
    if (colsize[Acol[i]] == 1) {
      colHashes[Acol[i]] = Arow[i];
    } else {
      HighsHashHelpers::sparse_combine(rowHashes[Arow[i]], Acol[i],
                                       HighsHashHelpers::double_hash_code(
                                           Avalue[i] / rowMax[Arow[i]].first));
      HighsHashHelpers::sparse_combine(colHashes[Acol[i]], Arow[i],
                                       HighsHashHelpers::double_hash_code(
                                           Avalue[i] / colMax[Acol[i]].first));
    }
  }

  // Step 3: Loop over the rows and columns and put them into buckets using the
  // computed hash values. Whenever a bucket already contains a row/column,
  // check if we can apply a (nearly) parallel row reduction or a
  // parallel/dominated column reduction.
  std::unordered_multimap<std::uint64_t, HighsInt> buckets;

  const bool debug_report = false;
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (colDeleted[i]) continue;
    if (colsize[i] == 0) {
      HPRESOLVE_CHECKED_CALL(colPresolve(postsolve_stack, i));
      continue;
    }
    auto it = buckets.find(colHashes[i]);
    decltype(it) last = it;

    HighsInt delCol = -1;
    HighsInt parallelColCandidate = -2;

    while (it != buckets.end() && it->first == colHashes[i]) {
      parallelColCandidate = it->second;
      last = it++;

      // we want to check if the columns are parallel, first rule out
      // hash collisions with different size columns
      if (colsize[i] != colsize[parallelColCandidate]) continue;
      // The columns have the same length. Next we determine whether domination
      // is possible in one of the directions, and if it is we designate the
      // dominating column as column 2. The first thing we check is whether the
      // the objective value of one of the (scaled) columns is strictly better
      // then the objective value of the other column which rules out domination
      // in one direction.

      HighsInt col = -1;
      HighsInt duplicateCol = -1;
      double colScale;

      // helpers for checking dominance between parallel columns which is
      // possible for different cases of the variable types: if col can be
      // increased infinitely in which case duplicateCol can be fixed to its
      // lower bound. duplicateCol can be decreased infinitely in which case col
      // can be fixed to its upper bound. for both cases we exploit that the
      // column that remains unfixed can always compensate for the fixed column.
      // This only holds if the compensating column can compensate exactly for
      // feasible value of the fixed column. In the continuous case this
      // trivially holds. In the case where both variables are integer and the
      // scale is +- 1 this also holds trivially. If the scale is > 1 and both
      // variables are integer, this only holds in one direction. We can apply
      // the reduction due to the following reasoning: Applying the scale to
      // col, means we change its meaning and it is not an integer variable
      // anymore, but a variable that moves on multiples of 1/scale. As we have
      // taken care that the scale is >=1 and integral for two integer
      // variables, the scaled column can always exactly compensate for the
      // other column as it can move by 1/k with k being integer. Hence every
      // kth allowed value is integral and no integral value is skipped. If the
      // compensating column is integral
      bool checkColImplBounds = true;
      bool checkDuplicateColImplBounds = true;
      auto colUpperInf = [&]() {
        if (!checkColImplBounds) return false;
        if (mipsolver == nullptr) {
          // for LP we check strict redundancy of the bounds as otherwise dual
          // postsolve might fail when the bound is used in the optimal solution
          return colScale > 0 ? isUpperStrictlyImplied(col)
                              : isLowerStrictlyImplied(col);
        } else {
          // for MIP we do not need dual postsolve so the reduction is valid if
          // the bound is weakly redundant
          return colScale > 0 ? isUpperImplied(col) : isLowerImplied(col);
        }
      };

      auto colLowerInf = [&]() {
        if (!checkColImplBounds) return false;
        if (mipsolver == nullptr) {
          return colScale > 0 ? isLowerStrictlyImplied(col)
                              : isUpperStrictlyImplied(col);
        } else {
          return colScale > 0 ? isLowerImplied(col) : isUpperImplied(col);
        }
      };

      auto duplicateColUpperInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        if (mipsolver == nullptr) {
          return isUpperStrictlyImplied(duplicateCol);
        } else {
          return isUpperImplied(duplicateCol);
        }
      };

      auto duplicateColLowerInf = [&]() {
        if (!checkDuplicateColImplBounds) return false;
        if (mipsolver == nullptr) {
          return isLowerStrictlyImplied(duplicateCol);
        } else {
          return isLowerImplied(duplicateCol);
        }
      };

      // Now check the if the variable types rule out domination in one
      // direction and already skip the column if that rules out domination in
      // both directions due to the previous check on the objective.
      if (model->integrality_[i] == HighsVarType::kInteger &&
          model->integrality_[parallelColCandidate] == HighsVarType::kInteger) {
        // both variables are integral, hence the scale must be integral
        // therefore first choose the smaller colMax value for col2, then check
        // integrality of colMax[col1] / colMax[col2].
        if (std::abs(colMax[i].first) <
            std::abs(colMax[parallelColCandidate].first)) {
          col = i;
          duplicateCol = parallelColCandidate;
        } else {
          col = parallelColCandidate;
          duplicateCol = i;
        }

        double scaleCand = colMax[duplicateCol].first / colMax[col].first;
        if (fractionality(scaleCand, &colScale) > options->small_matrix_value)
          continue;
        assert(std::abs(colScale) >= 1.0);

        // if the scale is larger than 1, duplicate column cannot compensate for
        // all values of scaled col due to integrality as the scaled column
        // moves on a grid of 1/scale.
        //
        // ToDo: Check whether this is too restrictive
        if (colScale != 1.0) checkDuplicateColImplBounds = false;
      } else if (model->integrality_[i] == HighsVarType::kInteger) {
        col = i;
        duplicateCol = parallelColCandidate;
        colScale = colMax[duplicateCol].first / colMax[col].first;

        // as col is integral and duplicateCol is not col cannot compensate for
        // duplicate col
        checkColImplBounds = false;
      } else {
        col = parallelColCandidate;
        duplicateCol = i;
        colScale = colMax[duplicateCol].first / colMax[col].first;

        // as col might be integral and duplicateCol is not integral. In that
        // case col cannot compensate for duplicate col
        checkColImplBounds =
            model->integrality_[parallelColCandidate] != HighsVarType::kInteger;
      }

      double objDiff = static_cast<double>(
          model->col_cost_[col] * static_cast<HighsCDouble>(colScale) -
          model->col_cost_[duplicateCol]);
      // if (std::abs(objDiff) > options->small_matrix_value) continue;
      constexpr HighsInt kMergeParallelCols = 0;
      constexpr HighsInt kDominanceColToUpper = 1;
      constexpr HighsInt kDominanceColToLower = 2;
      constexpr HighsInt kDominanceDuplicateColToUpper = 3;
      constexpr HighsInt kDominanceDuplicateColToLower = 4;

      HighsInt reductionCase = kMergeParallelCols;
      // now do the case distinctions for dominated columns
      // the cases are a lot simpler due to the helper functions
      // for checking the infinite bounds which automatically
      // incorporate the check for the variable types that allow domination.
      if (objDiff < -options->dual_feasibility_tolerance) {
        // scaled col is better than duplicate col
        if (colUpperInf() && model->col_lower_[duplicateCol] != -kHighsInf)
          reductionCase = kDominanceDuplicateColToLower;
        else if (duplicateColLowerInf() &&
                 (colScale < 0 || model->col_upper_[col] != kHighsInf) &&
                 (colScale > 0 || model->col_lower_[col] != -kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
        else
          continue;
      } else if (objDiff > options->dual_feasibility_tolerance) {
        // duplicate col is better than scaled col
        if (colLowerInf() && model->col_upper_[duplicateCol] != kHighsInf)
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf() &&
                 (colScale < 0 || model->col_lower_[col] != -kHighsInf) &&
                 (colScale > 0 || model->col_upper_[col] != kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else
          continue;
      } else {
        if (colUpperInf() && model->col_lower_[duplicateCol] != -kHighsInf)
          reductionCase = kDominanceDuplicateColToLower;
        else if (colLowerInf() && model->col_upper_[duplicateCol] != kHighsInf)
          reductionCase = kDominanceDuplicateColToUpper;
        else if (duplicateColUpperInf() &&
                 (colScale < 0 || model->col_lower_[col] != -kHighsInf) &&
                 (colScale > 0 || model->col_upper_[col] != kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToLower : kDominanceColToUpper;
        else if (duplicateColLowerInf() &&
                 (colScale < 0 || model->col_upper_[col] != kHighsInf) &&
                 (colScale > 0 || model->col_lower_[col] != -kHighsInf))
          reductionCase =
              colScale > 0 ? kDominanceColToUpper : kDominanceColToLower;
      }
      if (reductionCase == kMergeParallelCols) {
        const bool x_int = model->integrality_[col] == HighsVarType::kInteger;
        const bool y_int =
            model->integrality_[duplicateCol] == HighsVarType::kInteger;
        bool illegal_scale = true;
        if (x_int) {
          // The only possible reduction if the column parallelism check
          // succeeds is to merge the two columns into one. If one column is
          // integral this means we have restrictions on integers and need to
          // check additional conditions to allow the merging of two integer
          // columns, or a continuous column and an integer.
          if (model->integrality_[duplicateCol] != HighsVarType::kInteger) {
            assert(!y_int);
            // only one column is integral which cannot be duplicateCol due to
            // the way we assign the columns above
            //
            // Scale must not exceed 1/(y_u-y_l) in magnitude
            illegal_scale =
                std::abs(colScale * (model->col_upper_[duplicateCol] -
                                     model->col_lower_[duplicateCol])) <
                1.0 - primal_feastol;
            if (!illegal_scale && debug_report)
              printf(
                  "kMergeParallelCols: T-F is %s legal with scale %.4g and "
                  "duplicateCol = [%.4g, %.4g]\n",
                  illegal_scale ? "not" : "   ", colScale,
                  model->col_lower_[duplicateCol],
                  model->col_upper_[duplicateCol]);
          } else {
            // Both columns integer
            assert(x_int && y_int);
            // Scale must be integer and not exceed (x_u-x_l)+1 in magnitude
            const double scale_limit = model->col_upper_[col] -
                                       model->col_lower_[col] + 1 +
                                       primal_feastol;
            illegal_scale = std::fabs(colScale) > scale_limit;
          }
          if (illegal_scale) continue;
        } else {
          // Neither column integer: no problem with
          assert(!x_int && !y_int);
        }
      }

      bool parallel = true;
      // now check whether the coefficients are actually parallel
      for (const HighsSliceNonzero& colNz : getColumnVector(col)) {
        HighsInt duplicateColRowPos = findNonzero(colNz.index(), duplicateCol);
        parallel = duplicateColRowPos != -1;
        if (!parallel) break;

        parallel = abs(Avalue[duplicateColRowPos] -
                       static_cast<HighsCDouble>(colScale) * colNz.value()) <=
                   options->small_matrix_value;
        if (!parallel) break;
      }

      if (!parallel) continue;

      switch (reductionCase) {
        case kDominanceDuplicateColToLower:
          delCol = duplicateCol;
          HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, duplicateCol));
          break;
        case kDominanceDuplicateColToUpper:
          delCol = duplicateCol;
          HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, duplicateCol));
          break;
        case kDominanceColToLower:
          delCol = col;
          HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
          break;
        case kDominanceColToUpper:
          delCol = col;
          HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
          break;
        case kMergeParallelCols:
          const bool ok_merge = postsolve_stack.duplicateColumn(
              colScale, model->col_lower_[col], model->col_upper_[col],
              model->col_lower_[duplicateCol], model->col_upper_[duplicateCol],
              col, duplicateCol,
              model->integrality_[col] == HighsVarType::kInteger,
              model->integrality_[duplicateCol] == HighsVarType::kInteger,
              options->mip_feasibility_tolerance);
          if (!ok_merge) {
            if (debug_report) {
              printf(
                  "HPresolve::detectParallelRowsAndCols Illegal merge "
                  "prevented\n");
            }
            continue;
          }
          // When merging a continuous variable into an integer
          // variable, the integer will become continuous - since any
          // value in its range can be mapped back to an integer and a
          // continuous variable. Hence the number of integer
          // variables in the rows corresponding to the former integer
          // variable reduces.
          //
          // With the opposite - merging an integer variable into a
          // continuous variable - the retained variable is
          // continuous, so no action is required
          bool rowsizeIntReduction =
              model->integrality_[duplicateCol] != HighsVarType::kInteger &&
              model->integrality_[col] == HighsVarType::kInteger;
          if (rowsizeIntReduction)
            model->integrality_[col] = HighsVarType::kContinuous;

          markChangedCol(col);

          // compute bounds of merged variable
          double mergeLower = 0;
          double mergeUpper = 0;
          if (colScale > 0) {
            mergeLower = model->col_lower_[col] +
                         colScale * model->col_lower_[duplicateCol];
            mergeUpper = model->col_upper_[col] +
                         colScale * model->col_upper_[duplicateCol];
          } else {
            mergeLower = model->col_lower_[col] +
                         colScale * model->col_upper_[duplicateCol];
            mergeUpper = model->col_upper_[col] +
                         colScale * model->col_lower_[duplicateCol];
          }

          if (mipsolver != nullptr && mipsolver->mipdata_->cliquesExtracted) {
            presolveCliqueTable.eliminateCol(col);
          }

          // change bounds
          HPRESOLVE_CHECKED_CALL(changeColBounds(col, mergeLower, mergeUpper));

          // mark duplicate column as deleted
          markColDeleted(duplicateCol);
          // remove all nonzeros of duplicateCol
          for (HighsInt coliter = colhead[duplicateCol]; coliter != -1;) {
            assert(Acol[coliter] == duplicateCol);

            HighsInt colpos = coliter;
            HighsInt colrow = Arow[coliter];

            // if an integer column was merged into a continuous one make
            // sure to update the integral rowsize
            if (rowsizeIntReduction) rowsizeInteger[colrow] -= 1;

            coliter = Anext[coliter];

            unlink(colpos);

            reinsertEquation(colrow);
          }
          // set cost to zero
          model->col_cost_[duplicateCol] = 0;
          delCol = duplicateCol;

          // remove implied bounds, since they might in general not be valid
          // anymore
          resetColImpliedBounds(col);

          // if an implicit integer and an integer column were merged, check if
          // merged continuous column is implicit integer after merge
          if (rowsizeIntReduction && model->integrality_[duplicateCol] ==
                                         HighsVarType::kImplicitInteger) {
            StatusResult impliedInteger = isImpliedInteger(col);
            HPRESOLVE_CHECKED_CALL(static_cast<Result>(impliedInteger));
            if (impliedInteger) convertImpliedInteger(col, -1, true);
          }
          break;
      }

      break;
    }

    if (delCol != -1) {
      if (delCol != i) buckets.erase(last);

      // we could have new row singletons since a column was removed. Remove
      // those rows immediately
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    } else {
      buckets.emplace_hint(last, colHashes[i], i);
    }
  }

  buckets.clear();

  // Iterate non-cut rows before cut rows so that when a cut is
  // parallel to a non-cut, the non-cut is already in the bucket as
  // the surviving row and the cut is removed. This prevents a cut
  // from absorbing a non-cut constraint which would cause the
  // constraint to be lost when moveCutsToPool removes it.
  std::vector<HighsInt> rowOrder(model->num_row_);
  std::iota(rowOrder.begin(), rowOrder.end(), 0);
  pdqsort(rowOrder.begin(), rowOrder.end(), [&](HighsInt a, HighsInt b) {
    if (postsolve_stack.isCutRow(a) != postsolve_stack.isCutRow(b))
      return !postsolve_stack.isCutRow(a);
    return a < b;
  });

  for (HighsInt rowIndex = 0; rowIndex != model->num_row_; ++rowIndex) {
    HighsInt i = rowOrder[rowIndex];
    if (rowDeleted[i]) continue;
    if (rowsize[i] <= 1 || (rowsize[i] == 2 && isEquation(i)) ||
        rowsize[i] == rowsizeSingleton[i]) {
      HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, i));
      continue;
    }
    auto it = buckets.find(rowHashes[i]);
    decltype(it) last = it;

    // equalityRowAddition is only called if row i or the candidate row
    // contains singleton columns, i.e. if one of the rows is not exactly
    // parallel. hence, if may_require_basis_postsolve_ is true, row i is
    // skipped here if it contains singleton columns, and a candidate row is
    // skipped below if it contains singleton columns
    if (may_require_basis_postsolve_ && rowsizeSingleton[i] != 0) continue;

    // reductions with row i change or delete row i or the candidate rows, and
    // the rowPresolve called by equalityRowAddition can also change other
    // rows. hence, the hashes and the row maxima computed above only serve as
    // filters, and the parallel check below uses the current rows
    //
    // delRow is the row that the last reduction with row i changed or
    // deleted, and -1 if there was no reduction. delRow = i means that row i
    // has changed or has been deleted
    HighsInt delRow = -1;
    if (it != buckets.end()) storeRow(i);
    while (it != buckets.end() && it->first == rowHashes[i]) {
      HighsInt parallelRowCand = it->second;
      last = it++;

      if (rowDeleted[parallelRowCand]) continue;

      const HighsInt numSingleton = rowsizeSingleton[i];
      const HighsInt numSingletonCandidate = rowsizeSingleton[parallelRowCand];

      if (may_require_basis_postsolve_ && numSingletonCandidate != 0) continue;

      if (rowsize[i] - numSingleton !=
          rowsize[parallelRowCand] - numSingletonCandidate)
        continue;

      if (numSingletonCandidate > 1 || numSingleton > 1) {
        // we only handle the case where the rows have at most one extra
        // singleton except when one row has no extra singleton and is an
        // equation. In that case we sparsify the other row by adding the
        // equation and can subsequently solve it as an individual component as
        // it is a row which only contains singletons
        if ((numSingleton != 0 || !isEquation(i)) &&
            (numSingletonCandidate != 0 || !isEquation(parallelRowCand)))
          continue;
      } else if (numSingletonCandidate != numSingleton) {
        // if only one of the two constraints has an extra singleton,
        // we require at least one of the constraints to be an equation
        // if that is the case we can add that equation to the other row
        // and will make it into either a row singleton or a doubleton equation
        // which is removed afterwards
        if (!isEquation(i) && !isEquation(parallelRowCand)) continue;
      }

      double rowScale = rowMax[parallelRowCand].first / rowMax[i].first;
      // check parallel case
      bool parallel = true;
      for (const HighsSliceNonzero& rowNz : getStoredRow()) {
        if (colsize[rowNz.index()] == 1)  // skip singletons
          continue;
        HighsInt nzPos = findNonzero(parallelRowCand, rowNz.index());
        parallel = nzPos != -1;
        if (!parallel) break;

        parallel = abs(Avalue[nzPos] -
                       static_cast<HighsCDouble>(rowScale) * rowNz.value()) <=
                   options->small_matrix_value;
        if (!parallel) break;
      }
      if (!parallel) continue;

      if (may_require_basis_postsolve_)
        assert(numSingleton == 0 && numSingletonCandidate == 0);

      if (numSingleton == 0 && numSingletonCandidate == 0) {
        bool rowLowerTightened = false;
        bool rowUpperTightened = false;
        double newUpper;
        double newLower;
        if (rowScale > 0) {
          newUpper = model->row_upper_[i] * rowScale;
          newLower = model->row_lower_[i] * rowScale;
        } else {
          newLower = model->row_upper_[i] * rowScale;
          newUpper = model->row_lower_[i] * rowScale;
        }

        if (newUpper < model->row_upper_[parallelRowCand]) {
          if (newUpper < model->row_lower_[parallelRowCand] - primal_feastol)
            return Result::kPrimalInfeasible;

          if (newUpper < model->row_lower_[parallelRowCand])
            newUpper = model->row_lower_[parallelRowCand];

          if (newUpper < model->row_upper_[parallelRowCand]) {
            rowUpperTightened = true;
            if (rowScale > 0) {
              double tmp = rowDualLower[i] / rowScale;
              rowDualLower[i] = rowDualLower[parallelRowCand] * rowScale;
              rowDualLower[parallelRowCand] = tmp;
            } else {
              double tmp = rowDualUpper[i] / rowScale;
              rowDualUpper[i] = rowDualLower[parallelRowCand] * rowScale;
              rowDualLower[parallelRowCand] = tmp;
            }

            changeRowUpper(parallelRowCand, newUpper, true);
          }
        }

        if (newLower > model->row_lower_[parallelRowCand]) {
          if (newLower > model->row_upper_[parallelRowCand] + primal_feastol)
            return Result::kPrimalInfeasible;

          if (newLower > model->row_upper_[parallelRowCand])
            newLower = model->row_upper_[parallelRowCand];

          if (newLower > model->row_lower_[parallelRowCand]) {
            // the rows lower bound is tightened
            // instead of updating the activities of dual constraints, we
            // can simply swap the bounds on the row duals. If the old
            // lower bound on the row dual was finite, the new row dual
            // lower bound is infinite as the new row lower bound must be
            // a finite value. This infinite contribution, was, however,
            // already counted from the parallel row. Therefore by
            // swapping the bounds unlinking the other row will not
            // decrease the infinity counter, but simply remove a bound
            // with zero contribution. For a negative scale we need to
            // swap with the negated upper bound of the row dual of row i.
            rowLowerTightened = true;
            if (rowScale > 0) {
              double tmp = rowDualUpper[i] / rowScale;
              rowDualUpper[i] = rowDualUpper[parallelRowCand] * rowScale;
              rowDualUpper[parallelRowCand] = tmp;
            } else {
              double tmp = rowDualLower[i] / rowScale;
              rowDualLower[i] = rowDualUpper[parallelRowCand] * rowScale;
              rowDualUpper[parallelRowCand] = tmp;
            }

            changeRowLower(parallelRowCand, newLower, true);
          }
        }
        // remove implied bounds, since they might in general not be valid
        // anymore
        resetRowDualImpliedBounds(parallelRowCand);

        postsolve_stack.duplicateRow(parallelRowCand, rowUpperTightened,
                                     rowLowerTightened, i, rowScale);
        delRow = i;
        markRowDeleted(i);
        for (HighsInt rowiter : rowpositions) unlink(rowiter);
        break;
      } else if (isEquation(i)) {
        // row i is equation and parallel (except for singletons)
        // add to the row parallelRowCand
        // printf(
        //    "nearly parallel case with %" HIGHSINT_FORMAT " singletons in eq
        //    row and %" HIGHSINT_FORMAT " " "singletons in other row(eq=%"
        //    HIGHSINT_FORMAT ")\n", numSingleton, numSingletonCandidate,
        //    model->row_lower_[parallelRowCand] ==
        //        model->row_upper_[parallelRowCand]);
        assert(!may_require_basis_postsolve_);
        HPRESOLVE_CHECKED_CALL(equalityRowAddition(
            postsolve_stack, i, parallelRowCand, -rowScale, getStoredRow()));
        delRow = parallelRowCand;
        // the candidate row has changed, remove it from the buckets
        last = buckets.erase(last);
        // equalityRowAddition calls rowPresolve, which can delete or change
        // row i and overwrite rowpositions. stop if row i is deleted or only
        // contains singleton columns
        if (rowDeleted[i] || rowsize[i] == rowsizeSingleton[i]) break;
        storeRow(i);
      } else if (isEquation(parallelRowCand)) {
        // printf(
        //    "nearly parallel case with %" HIGHSINT_FORMAT " singletons in eq
        //    row and %" HIGHSINT_FORMAT " " "singletons in other inequality
        //    row\n", numSingletonCandidate, numSingleton);
        // the row parallelRowCand is an equation; add it to the other row
        assert(!may_require_basis_postsolve_);
        HPRESOLVE_CHECKED_CALL(equalityRowAddition(
            postsolve_stack, parallelRowCand, i,
            -rowMax[i].first / rowMax[parallelRowCand].first,
            getRowVector(parallelRowCand)));
        // row i has changed or has been deleted, stop scanning
        delRow = i;
        break;
      } else {
        assert(numSingleton == 1);
        assert(numSingletonCandidate == 1);

        double rowUpper;
        double rowLower;
        if (rowScale > 0) {
          rowUpper = model->row_upper_[i] * rowScale;
          rowLower = model->row_lower_[i] * rowScale;
        } else {
          rowLower = model->row_upper_[i] * rowScale;
          rowUpper = model->row_lower_[i] * rowScale;
        }
        // todo: two inequalities with one singleton. check whether the rows can
        // be converted to equations by introducing a shared slack variable
        // which is the case if the singletons have similar properties
        // (objective sign, bounds, scaled coefficient) and the scaled right
        // hand sides match. Then the case reduces to adding one equation to the
        // other and substituting one of the singletons due to the resulting
        // doubleton equation.
        //        printf("todo, two inequalities with one additional
        //        singleton\n");
        (void)rowLower;
        (void)rowUpper;
      }
    }

    // reinsert row i only if it is unchanged and contains non-singleton
    // columns. the rowPresolve called by equalityRowAddition may also have
    // deleted row i
    if (delRow != i && !rowDeleted[i] && rowsize[i] != rowsizeSingleton[i])
      buckets.emplace_hint(last, rowHashes[i], i);

    if (delRow != -1) HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  analysis_.logging_on_ = logging_on;
  if (logging_on)
    analysis_.stopPresolveRuleLog(kPresolveRuleParallelRowsAndCols);

  return Result::kOk;
}

template <typename RowStorageFormat>
HPresolve::Result HPresolve::equalityRowAddition(
    HighsPostsolveStack& postsolve_stack, HighsInt stayrow, HighsInt removerow,
    double scale, const HighsMatrixSlice<RowStorageFormat>& rowvector) {
  assert(!may_require_basis_postsolve_);
  // extract non-zero positions
  std::vector<HighsInt> stay_rowpositions;
  getRowPositions(stayrow, stay_rowpositions);

  // update postsolve information
  postsolve_stack.equalityRowAddition(removerow, stayrow, scale, rowvector);

  // iterate over non-zero positions instead of iterating over the
  // HighsMatrixSlice because the latter contains pointers to Acol and Avalue
  // that may be invalidated if these vectors are reallocated
  // (see std::vector::push_back performed in HPresolve::addToMatrix).
  for (HighsInt rowiter : stay_rowpositions) {
    HighsInt pos = findNonzero(removerow, Acol[rowiter]);
    if (pos != -1)
      unlink(pos);  // all common nonzeros are cancelled, as the rows are
                    // parallel
    else            // might introduce a singleton
      addToMatrix(removerow, Acol[rowiter], scale * Avalue[rowiter]);
  }

  addToRowLower(removerow,
                static_cast<HighsCDouble>(scale) * model->row_upper_[stayrow]);
  addToRowUpper(removerow,
                static_cast<HighsCDouble>(scale) * model->row_upper_[stayrow]);
  // equalityRowAddition changes the size of removerow (unlink for
  // common nonzeros, addToMatrix for the rest), so removerow must be
  // re-keyed in the size-ordered equations set, otherwise the next
  // removeDoubletonEquations (in rowPresolve, for example) finds a
  // stale size
  //
  reinsertEquation(removerow);
  //
  // row is now a singleton row, doubleton equation, or a row
  // that contains only singletons and we let the normal row presolve
  // handle the cases
  HPRESOLVE_CHECKED_CALL(rowPresolve(postsolve_stack, removerow));
  return Result::kOk;
}

}  // namespace presolve
