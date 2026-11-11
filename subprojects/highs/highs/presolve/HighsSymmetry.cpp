/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsSymmetry.cpp
 * @brief Facilities for symmetry detection
 * @author Leona Gottwald
 */

#include "presolve/HighsSymmetry.h"

#include <algorithm>
#include <numeric>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomain.h"
#include "parallel/HighsParallel.h"
#include "util/HighsDisjointSets.h"
#include "util/HighsTimer.h"

void HighsSymmetryDetection::removeFixPoints() {
  Gend.resize(numVertices);
  for (HighsInt i = 0; i < numVertices; ++i) {
    Gend[i] =
        std::partition(Gedge.begin() + Gstart[i], Gedge.begin() + Gstart[i + 1],
                       [&](const std::pair<HighsInt, HighsUInt>& edge) {
                         return cellSize(vertexToCell[edge.first]) > 1;
                       }) -
        Gedge.begin();
    assert(Gend[i] >= Gstart[i] && Gend[i] <= Gstart[i + 1]);
  }

  HighsInt unitCellIndex = numVertices;
  currentPartition.erase(
      std::remove_if(currentPartition.begin(), currentPartition.end(),
                     [&](HighsInt vertex) {
                       if (cellSize(vertexToCell[vertex]) == 1) {
                         --unitCellIndex;
                         vertexToCell[vertex] = unitCellIndex;
                         return true;
                       }
                       return false;
                     }),
      currentPartition.end());

  for (HighsInt i = 0; i < numVertices; ++i) {
    if (Gend[i] == Gstart[i + 1]) continue;

    for (HighsInt j = Gend[i]; j < Gstart[i + 1]; ++j)
      Gedge[j].first = vertexToCell[Gedge[j].first];
  }

  if ((HighsInt)currentPartition.size() < numVertices) {
    numVertices = currentPartition.size();
    if (numVertices == 0) {
      numActiveCols = 0;
      return;
    }
    currentPartitionLinks.resize(numVertices);
    cellInRefinementQueue.assign(numVertices, false);
    assert(refinementQueue.empty());
    refinementQueue.clear();
    HighsInt cellStart = 0;
    HighsInt cellNumber = 0;
    for (HighsInt i = 0; i < numVertices; ++i) {
      HighsInt vertex = currentPartition[i];
      // if the cell number is different to the current cell number this is the
      // start of a new cell
      if (cellNumber != vertexToCell[vertex]) {
        // remember the number of this cell to identify its end
        cellNumber = vertexToCell[vertex];
        // set the link of the cell start to point to its end
        currentPartitionLinks[cellStart] = i;
        // remember start of this cell
        cellStart = i;
      }

      // correct the vertexToCell array to store the start index of the
      // cell, not its number
      updateCellMembership(i, cellStart, false);
    }

    // set the column partition link of the last started cell to point past the
    // end
    assert((int)currentPartitionLinks.size() > 0);
    currentPartitionLinks[cellStart] = numVertices;

    numActiveCols =
        std::partition_point(currentPartition.begin(), currentPartition.end(),
                             [&](HighsInt v) { return v < numCol; }) -
        currentPartition.begin();
  } else
    numActiveCols = numCol;
}

void HighsSymmetries::clear() {
  permutationColumns.clear();
  permutations.clear();
  orbitPartition.clear();
  orbitSize.clear();
  columnPosition.clear();
  linkCompressionStack.clear();
  columnToOrbitope.clear();
  orbitopes.clear();
  numPerms = 0;
  numGenerators = 0;
}

void HighsSymmetries::mergeOrbits(HighsInt v1, HighsInt v2,
                                  std::vector<HighsInt>& orbitPartition,
                                  std::vector<HighsInt>& orbitSize,
                                  std::vector<HighsInt>& linkCompressionStack) {
  if (v1 == v2) return;

  const HighsInt orbit1 = getOrbit(v1, orbitPartition, linkCompressionStack);
  const HighsInt orbit2 = getOrbit(v2, orbitPartition, linkCompressionStack);

  if (orbit1 == orbit2) return;

  if (orbitSize[orbit2] < orbitSize[orbit1]) {
    orbitPartition[orbit2] = orbit1;
    orbitSize[orbit1] += orbitSize[orbit2];
  } else {
    orbitPartition[orbit1] = orbit2;
    orbitSize[orbit2] += orbitSize[orbit1];
  }
}

HighsInt HighsSymmetries::getOrbit(
    HighsInt col, std::vector<HighsInt>& orbitPartition,
    std::vector<HighsInt>& linkCompressionStack) {
  HighsInt i = columnPosition[col];
  if (i == -1) return -1;

  HighsInt orbit = orbitPartition[i];
  if (orbit != orbitPartition[orbit]) {
    do {
      linkCompressionStack.push_back(i);
      i = orbit;
      orbit = orbitPartition[orbit];
    } while (orbit != orbitPartition[orbit]);

    do {
      i = linkCompressionStack.back();
      linkCompressionStack.pop_back();
      orbitPartition[i] = orbit;
    } while (!linkCompressionStack.empty());
  }

  return orbit;
}

HighsInt HighsSymmetries::propagateOrbitopes(HighsDomain& domain) const {
  if (columnToOrbitope.size() == 0 || domain.getBranchDepth() == 0) return 0;

  std::set<HighsInt> propagationOrbitopes;
  const auto& domchgstack = domain.getDomainChangeStack();
  for (HighsInt pos : domain.getBranchingPositions()) {
    const HighsInt* orbitope = columnToOrbitope.find(domchgstack[pos].column);
    if (orbitope) propagationOrbitopes.insert(*orbitope);
  }

  HighsInt numFixed = 0;
  for (HighsInt i : propagationOrbitopes) {
    numFixed += orbitopes[i].orbitalFixing(domain);
    if (domain.infeasible()) break;
  }

  // if (numFixed)
  //   printf("orbital fixing for full orbitope found %d fixings\n", numFixed);

  return numFixed;
}

std::shared_ptr<const StabilizerOrbits>
HighsSymmetries::computeStabilizerOrbits(const HighsDomain& localdom,
                                         StabilizerOrbitWorkspace& workspace) {
  const auto& domchgStack = localdom.getDomainChangeStack();
  const auto& branchingPos = localdom.getBranchingPositions();

  StabilizerOrbits stabilizerOrbits;
  stabilizerOrbits.stabilizedCols.reserve(permutationColumns.size());
  for (HighsInt i : branchingPos) {
    HighsInt col = domchgStack[i].column;
    if (columnPosition[col] == -1) continue;

    assert(localdom.variableType(col) != HighsVarType::kContinuous);

    // if we branch a variable upwards it is either binary and branched to one
    // and needs to be stabilized or it is a general integer and needs to be
    // stabilized regardless of branching direction.
    // if we branch downwards we only need to stabilize on this
    // branching column if it is a general integer
    if (domchgStack[i].boundtype == HighsBoundType::kLower ||
        !localdom.isGlobalBinary(col))
      stabilizerOrbits.stabilizedCols.push_back(columnPosition[col]);
  }

  HighsInt permLength = permutationColumns.size();
  workspace.orbitPartition.resize(permLength);
  std::iota(workspace.orbitPartition.begin(), workspace.orbitPartition.end(),
            0);
  workspace.orbitSize.assign(permLength, 1);
  workspace.linkCompressionStack.clear();

  for (HighsInt i = 0; i < numPerms; ++i) {
    const HighsInt* perm = permutations.data() + i * permutationColumns.size();

    bool permRespectsBranchings = true;
    for (HighsInt i : stabilizerOrbits.stabilizedCols) {
      if (permutationColumns[i] != perm[i]) {
        permRespectsBranchings = false;
        break;
      }
    }

    if (!permRespectsBranchings) continue;

    for (HighsInt j = 0; j < permLength; ++j) {
      mergeOrbits(permutationColumns[j], perm[j], workspace.orbitPartition,
                  workspace.orbitSize, workspace.linkCompressionStack);
    }
  }

  stabilizerOrbits.stabilizedCols.clear();

  stabilizerOrbits.orbitCols.reserve(permLength);
  for (HighsInt i = 0; i < permLength; ++i) {
    if (localdom.variableType(permutationColumns[i]) ==
        HighsVarType::kContinuous)
      continue;
    HighsInt orbit = getOrbit(permutationColumns[i], workspace.orbitPartition,
                              workspace.linkCompressionStack);
    if (workspace.orbitSize[orbit] == 1)
      stabilizerOrbits.stabilizedCols.push_back(permutationColumns[i]);
    else if (localdom.isGlobalBinary(permutationColumns[i]))
      stabilizerOrbits.orbitCols.push_back(permutationColumns[i]);
  }
  stabilizerOrbits.symmetries = this;
  pdqsort(stabilizerOrbits.stabilizedCols.begin(),
          stabilizerOrbits.stabilizedCols.end());
  if (!stabilizerOrbits.orbitCols.empty()) {
    pdqsort(stabilizerOrbits.orbitCols.begin(),
            stabilizerOrbits.orbitCols.end(),
            [&](HighsInt col1, HighsInt col2) {
              return getOrbit(col1, workspace.orbitPartition,
                              workspace.linkCompressionStack) <
                     getOrbit(col2, workspace.orbitPartition,
                              workspace.linkCompressionStack);
            });
    HighsInt numOrbitCols = stabilizerOrbits.orbitCols.size();
    stabilizerOrbits.orbitStarts.reserve(numOrbitCols + 1);
    stabilizerOrbits.orbitStarts.push_back(0);

    for (HighsInt i = 1; i < numOrbitCols; ++i) {
      if (getOrbit(stabilizerOrbits.orbitCols[i], workspace.orbitPartition,
                   workspace.linkCompressionStack) !=
          getOrbit(stabilizerOrbits.orbitCols[i - 1], workspace.orbitPartition,
                   workspace.linkCompressionStack))
        stabilizerOrbits.orbitStarts.push_back(i);
    }
    stabilizerOrbits.orbitStarts.push_back(numOrbitCols);
  }

  return std::make_shared<const StabilizerOrbits>(std::move(stabilizerOrbits));
}

HighsInt StabilizerOrbits::orbitalFixing(HighsDomain& domain) const {
  HighsInt numFixed = symmetries->propagateOrbitopes(domain);
  if (domain.infeasible() || orbitCols.empty()) return numFixed;

  HighsInt numOrbits = orbitStarts.size() - 1;
  for (HighsInt i = 0; i < numOrbits; ++i) {
    HighsInt fixcol = -1;
    for (HighsInt j = orbitStarts[i]; j < orbitStarts[i + 1]; ++j) {
      if (domain.isFixed(orbitCols[j])) {
        fixcol = orbitCols[j];
        break;
      }
    }

    if (fixcol != -1) {
      HighsInt oldNumFixed = numFixed;
      auto oldSize = domain.getDomainChangeStack().size();
      if (domain.col_lower_[fixcol] == 1.0) {
        for (HighsInt j = orbitStarts[i]; j < orbitStarts[i + 1]; ++j) {
          if (domain.col_lower_[orbitCols[j]] == 1.0) continue;
          ++numFixed;
          domain.changeBound(HighsBoundType::kLower, orbitCols[j], 1.0,
                             HighsDomain::Reason::unspecified());
          if (domain.infeasible()) return numFixed;
        }
      } else {
        for (HighsInt j = orbitStarts[i]; j < orbitStarts[i + 1]; ++j) {
          if (domain.col_upper_[orbitCols[j]] == 0.0) continue;
          ++numFixed;
          domain.changeBound(HighsBoundType::kUpper, orbitCols[j], 0.0,
                             HighsDomain::Reason::unspecified());
          if (domain.infeasible()) return numFixed;
        }
      }

      HighsInt newFixed = numFixed - oldNumFixed;

      if (newFixed != 0) {
        domain.propagate();
        if (domain.infeasible()) return numFixed;
        if ((HighsInt)(domain.getDomainChangeStack().size() - oldSize) >
            newFixed)
          i = -1;
      }
    }
  }

  return numFixed;
}

bool StabilizerOrbits::isStabilized(HighsInt col) const {
  return symmetries->columnPosition[col] == -1 ||
         std::binary_search(stabilizedCols.begin(), stabilizedCols.end(), col);
}

void HighsOrbitopeMatrix::detectSetPackingRows(HighsCliqueTable& cliquetable,
                                               HighsInt cliqueVal) {
  RowPackingStatus resultStatus =
      (cliqueVal == 1) ? kRowPacking : kRowPackingNegated;

  for (HighsInt j = 1; j < rowLength; ++j) {
    HighsInt* colj1 = &entry(0, j);

    for (HighsInt j0 = 0; j0 < j; ++j0) {
      HighsInt* colj0 = &entry(0, j0);

      for (HighsInt i = 0; i < numRows; ++i) {
        if (rowIsSetPacking[i] != kRowUndetermined) continue;

        HighsInt xj0 = colj0[i];
        HighsInt xj1 = colj1[i];

        auto commonClique =
            cliquetable.findCommonClique({xj0, cliqueVal}, {xj1, cliqueVal});

        if (commonClique.first == nullptr) {
          rowIsSetPacking[i] = kRowNotPacking;
          continue;
        }

        HighsInt overlap = 0;

        for (HighsInt k = 0; k < commonClique.second; ++k) {
          if (commonClique.first[k].val == 1 - cliqueVal) continue;

          HighsInt* cliqueColRow = columnToRow.find(commonClique.first[k].col);
          if (cliqueColRow && *cliqueColRow == i) ++overlap;
        }

        if (overlap == rowLength) {
          rowIsSetPacking[i] = resultStatus;
          ++numSetPackingRows;
          if (numSetPackingRows == numRows) return;
        }
      }
    }
  }
}

void HighsOrbitopeMatrix::determineOrbitopeType(HighsCliqueTable& cliquetable) {
  for (HighsInt j = 0; j < rowLength; ++j) {
    for (HighsInt i = 0; i < numRows; ++i) {
      columnToRow.insert(entry(i, j), i);
    }
  }

  rowIsSetPacking.assign(numRows, kRowUndetermined);
  numSetPackingRows = 0;

  detectSetPackingRows(cliquetable, HighsInt{1});

  if (numSetPackingRows == numRows) return;

  // now for the rows that do not have a set packing structure check
  // if we have such structure when all columns in the row are negated

  for (HighsInt i = 0; i < numRows; ++i) {
    if (rowIsSetPacking[i] == kRowNotPacking)
      rowIsSetPacking[i] = kRowUndetermined;
  }

  detectSetPackingRows(cliquetable, HighsInt{0});
}

HighsInt HighsOrbitopeMatrix::getBranchingColumn(
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    HighsInt col) const {
  const HighsInt* i = columnToRow.find(col);
  if (i && rowIsSetPacking[*i] > kRowNotPacking) {
    for (HighsInt j = 0; j < rowLength; ++j) {
      HighsInt branchCol = entry(*i, j);
      if (branchCol == col) break;
      if (colLower[branchCol] != colUpper[branchCol]) return branchCol;
    }
  }

  return col;
}

HighsInt HighsOrbitopeMatrix::orbitalFixingForPackingOrbitope(
    const std::vector<HighsInt>& rows, HighsDomain& domain) const {
  HighsInt numDynamicRows = rows.size();

  // printf("propagate packing orbitope\n");

  std::vector<HighsInt> firstOneInRow(numDynamicRows, -1);

  for (HighsInt j = 0; j < rowLength; ++j) {
    for (HighsInt i = 0; i < numDynamicRows; ++i) {
      if (firstOneInRow[i] != -1) continue;
      HighsInt r = rows[i];
      HighsInt colrj = entry(r, j);
      if (rowIsSetPacking[r] == kRowPacking) {
        if (domain.col_lower_[colrj] > 0.5) firstOneInRow[i] = j;
      } else {
        assert(rowIsSetPacking[r] == kRowPackingNegated);
        if (domain.col_upper_[colrj] < 0.5) firstOneInRow[i] = j;
      }
    }
  }

  // we start looping over the rows and keep the index j at the column that
  // is the right most possible location for a 1 entry.
  // For the first row this position is 0.

  HighsInt j = 0;
  HighsInt numFixed = 0;

  // Fix row entries from 0 to kEnd - 1 at orbitope column to zero.
  // Returns true if the domain becomes infeasible.
  auto fixColumnsToZero = [&](HighsInt kEnd, HighsInt col_j) {
    for (HighsInt k = 0; k < kEnd; ++k) {
      assert(firstOneInRow[k] < j);
      HighsInt col_kj = entry(rows[k], col_j);
      bool negate_k = rowIsSetPacking[rows[k]] == kRowPackingNegated;
      if (negate_k) {
        if (domain.col_lower_[col_kj] > 0.5) continue;
        domain.changeBound(HighsBoundType::kLower, col_kj, 1.0,
                           HighsDomain::Reason::unspecified());
      } else {
        if (domain.col_upper_[col_kj] < 0.5) continue;
        domain.changeBound(HighsBoundType::kUpper, col_kj, 0.0,
                           HighsDomain::Reason::unspecified());
      }
      ++numFixed;
      if (domain.infeasible()) return true;
    }
    return false;
  };

  for (HighsInt i = 0; i < numDynamicRows; ++i) {
    // at this position we know that the last possible position
    // of a 1-entry in this row is j. If we have a 1 behind j
    // the face of the orbitope intersecting the set of initial fixings is empty
    if (firstOneInRow[i] > j) {
      domain.markInfeasible();
      return numFixed;
    }
    HighsInt col_ij = entry(rows[i], j);

    bool negate_i = rowIsSetPacking[rows[i]] == kRowPackingNegated;

    bool notZeroFixed = negate_i ? domain.col_lower_[col_ij] < 0.5
                                 : domain.col_upper_[col_ij] > 0.5;

    // as long as the entry is fixed to zero
    // the frontier stays at the same column
    // if we encounter an entry that is not fixed to zero
    // we need to proceed with the next column and found a frontier step
    if (notZeroFixed) {
      // found a frontier step. Now we first check for the current column
      // if it can be fixed to 1.
      // For this we check if we find an infeasibility in the case where this
      // column was fixed to zero in which case the frontier would have stayed
      // at j for the following rows.
      HighsInt j0 = j;
      for (HighsInt k = i + 1; k < numDynamicRows; ++k) {
        if (firstOneInRow[k] > j0) {
          if (negate_i)
            domain.changeBound(HighsBoundType::kUpper, col_ij, 0.0,
                               HighsDomain::Reason::unspecified());
          else
            domain.changeBound(HighsBoundType::kLower, col_ij, 1.0,
                               HighsDomain::Reason::unspecified());

          ++numFixed;
          if (domain.infeasible()) return numFixed;
          break;
        }

        bool negate_k = rowIsSetPacking[rows[k]] == kRowPackingNegated;
        bool notZeroFixed = negate_k
                                ? domain.col_lower_[entry(rows[k], j0)] < 0.5
                                : domain.col_upper_[entry(rows[k], j0)] > 0.5;

        if (notZeroFixed) {
          ++j0;
          if (j0 == rowLength) break;
        }
      }

      ++j;
      if (j == rowLength) break;

      if (fixColumnsToZero(i + 1, j)) return numFixed;
    }
  }

  // check if there are more columns that can be completely fixed to zero
  while (++j < rowLength) {
    if (fixColumnsToZero(numDynamicRows, j)) return numFixed;
  }

  if (!domain.infeasible() && numFixed) domain.propagate();

  // if (numFixed)
  //  printf("orbital fixing for packing case fixed %d columns\n", numFixed);

  return numFixed;
}

HighsInt HighsOrbitopeMatrix::orbitalFixingForFullOrbitope(
    const std::vector<HighsInt>& rows, HighsDomain& domain) const {
  HighsInt numDynamicRows = rows.size();
  std::vector<int8_t> Mminimal(numDynamicRows * rowLength, -1);

  for (HighsInt j = 0; j < rowLength; ++j) {
    for (HighsInt i = 0; i < numDynamicRows; ++i) {
      HighsInt r = rows[i];
      HighsInt colij = matrix[r + j * numRows];
      if (domain.col_lower_[colij] == 1.0)
        Mminimal[i + j * numDynamicRows] = 1;
      else if (domain.col_upper_[colij] == 0.0)
        Mminimal[i + j * numDynamicRows] = 0;
    }
  }

  std::vector<int8_t> Mmaximal = Mminimal;

  int8_t* MminimaljLast = Mminimal.data() + numDynamicRows * (rowLength - 1);
  int8_t* MmaximaljFirst = Mmaximal.data();
  for (HighsInt k = 0; k < numDynamicRows; ++k) {
    if (MminimaljLast[k] == -1) MminimaljLast[k] = 0;
    if (MmaximaljFirst[k] == -1) MmaximaljFirst[k] = 1;
  }

  auto i_fixed = [&](const int8_t* colj0, const int8_t* colj1) {
    for (HighsInt i = 0; i < numDynamicRows; ++i) {
      if (colj0[i] != -1 && colj1[i] != -1 && colj0[i] != colj1[i]) return i;
    }

    return numDynamicRows;
  };

  auto i_discr = [&](const int8_t* colj0, const int8_t* colj1, HighsInt i_f) {
    for (HighsInt i = i_f; i >= 0; --i) {
      if (colj0[i] != 0 && colj1[i] != 1) return i;
    }

    return HighsInt{-1};
  };

  // Propagate lexicographic bounds column by column.
  // direction = -1: compute minimal matrix (right to left)
  // direction = +1: compute maximal matrix (left to right)
  auto propagateColumns = [&](int8_t* M, HighsInt jStart, HighsInt jEnd,
                              HighsInt direction, int8_t fixedVal,
                              int8_t fillAbove) {
    for (HighsInt j = jStart; j != jEnd; j += direction) {
      int8_t* colj0 = M + j * numDynamicRows;
      int8_t* colj1 = colj0 - direction * numDynamicRows;

      const int8_t* first = (direction < 0) ? colj0 : colj1;
      const int8_t* second = (direction < 0) ? colj1 : colj0;
      HighsInt i_f = i_fixed(first, second);

      if (i_f == numDynamicRows) {
        for (HighsInt k = 0; k < numDynamicRows; ++k) {
          int8_t isFree = (colj0[k] == -1);
          colj0[k] += (isFree & colj1[k]) + isFree;
        }
      } else {
        HighsInt i_d = i_discr(first, second, i_f);
        if (i_d == -1) {
          domain.markInfeasible();
          return;
        }
        for (HighsInt k = 0; k < i_d; ++k) {
          int8_t isFree = (colj0[k] == -1);
          colj0[k] += (isFree & colj1[k]) + isFree;
        }
        colj0[i_d] = fixedVal;
        for (HighsInt k = i_d + 1; k < numDynamicRows; ++k)
          colj0[k] += fillAbove * (colj0[k] == -1);
      }
    }
  };

  propagateColumns(Mminimal.data(), rowLength - 2, HighsInt{-1}, HighsInt{-1},
                   int8_t{1}, int8_t{1});
  if (domain.infeasible()) return 0;

  propagateColumns(Mmaximal.data(), HighsInt{1}, rowLength, HighsInt{1},
                   int8_t{0}, int8_t{2});
  if (domain.infeasible()) return 0;

  HighsInt numFixed = 0;

  for (HighsInt j = 0; j < rowLength; ++j) {
    const int8_t* colMaximal = Mmaximal.data() + j * numDynamicRows;
    const int8_t* colMinimal = Mminimal.data() + j * numDynamicRows;

    for (HighsInt i = 0; i < numDynamicRows; ++i) {
      if (colMinimal[i] != colMaximal[i]) {
        assert(colMinimal[i] < colMaximal[i]);
        break;
      }

      HighsInt r = rows[i];
      HighsInt colrj = matrix[r + j * numRows];
      if (domain.isFixed(colrj)) continue;

      ++numFixed;
      if (colMinimal[i] == 1)
        domain.changeBound(HighsBoundType::kLower, colrj, 1.0,
                           HighsDomain::Reason::unspecified());
      else
        domain.changeBound(HighsBoundType::kUpper, colrj, 0.0,
                           HighsDomain::Reason::unspecified());
      if (domain.infeasible()) break;
    }
    if (domain.infeasible()) break;
  }

  if (!domain.infeasible()) domain.propagate();

  return numFixed;
}

HighsInt HighsOrbitopeMatrix::orbitalFixing(HighsDomain& domain) const {
  std::vector<HighsInt> rows;
  std::vector<HighsBool> rowUsed(numRows);

  rows.reserve(numRows);

  const auto& branchpos = domain.getBranchingPositions();
  const auto& domchgstack = domain.getDomainChangeStack();

  bool isPacking = true;
  for (HighsInt pos : branchpos) {
    const HighsInt* i = columnToRow.find(domchgstack[pos].column);
    if (i && !rowUsed[*i]) {
      rowUsed[*i] = true;
      isPacking = isPacking && rowIsSetPacking[*i] > kRowNotPacking;
      rows.push_back(*i);
    }
  }

  if (rows.empty()) return 0;

  if (isPacking) return orbitalFixingForPackingOrbitope(rows, domain);

  return orbitalFixingForFullOrbitope(rows, domain);
}

void HighsSymmetryDetection::initializeGroundSet() {
  vertexGroundSet = currentPartition;
  pdqsort(vertexGroundSet.begin(), vertexGroundSet.end());
  vertexPosition.resize(vertexToCell.size(), -1);
  for (HighsInt i = 0; i < numVertices; ++i)
    vertexPosition[vertexGroundSet[i]] = i;

  orbitPartition.resize(numVertices);
  std::iota(orbitPartition.begin(), orbitPartition.end(), 0);
  orbitSize.assign(numVertices, 1);

  automorphisms.resize(numVertices * 64);
  numAutomorphisms = 0;
  currNodeCertificate.reserve(numVertices);
}

bool HighsSymmetryDetection::mergeOrbits(HighsInt v1, HighsInt v2) {
  if (v1 == v2) return false;

  HighsInt orbit1 = getOrbit(v1);
  HighsInt orbit2 = getOrbit(v2);

  if (orbit1 == orbit2) return false;

  if (orbit1 < orbit2) {
    orbitPartition[orbit2] = orbit1;
    orbitSize[orbit1] += orbitSize[orbit2];
  } else {
    orbitPartition[orbit1] = orbit2;
    orbitSize[orbit2] += orbitSize[orbit1];
  }

  return true;
}

HighsInt HighsSymmetryDetection::getOrbit(HighsInt vertex) {
  HighsInt i = vertexPosition[vertex];
  HighsInt orbit = orbitPartition[i];
  if (orbit != orbitPartition[orbit]) {
    do {
      linkCompressionStack.push_back(i);
      i = orbit;
      orbit = orbitPartition[orbit];
    } while (orbit != orbitPartition[orbit]);

    do {
      i = linkCompressionStack.back();
      linkCompressionStack.pop_back();
      orbitPartition[i] = orbit;
    } while (!linkCompressionStack.empty());
  }

  return orbit;
}

void HighsSymmetryDetection::initializeHashValues() {
  for (HighsInt i = 0; i != numVertices; ++i) {
    HighsInt cell = vertexToCell[i];
    for (HighsInt j = Gstart[i]; j != Gend[i]; ++j) {
      HighsHashHelpers::sparse_combine(vertexHash[Gedge[j].first], cell,
                                       Gedge[j].second);
    }
    markCellForRefinement(cell);
  }
}

bool HighsSymmetryDetection::updateCellMembership(HighsInt i, HighsInt cell,
                                                  bool markForRefinement) {
  HighsInt vertex = currentPartition[i];
  if (vertexToCell[vertex] != cell) {
    // set new cell id
    vertexToCell[vertex] = cell;
    if (i != cell) currentPartitionLinks[i] = cell;

    // update hashes of affected rows
    if (markForRefinement) {
      for (HighsInt j = Gstart[vertex]; j != Gend[vertex]; ++j) {
        HighsInt edgeDestinationCell = vertexToCell[Gedge[j].first];
        if (cellSize(edgeDestinationCell) == 1) continue;
        HighsHashHelpers::sparse_combine(vertexHash[Gedge[j].first], cell,
                                         Gedge[j].second);
        markCellForRefinement(edgeDestinationCell);
      }
    }

    return true;
  }

  return false;
}

bool HighsSymmetryDetection::splitCell(HighsInt cell, HighsInt splitPoint) {
  u32 hSplit = static_cast<u32>(getVertexHash(currentPartition[splitPoint]));
  u32 hCell = static_cast<u32>(getVertexHash(currentPartition[cell]));

  u32 certificateVal =
      (HighsHashHelpers::pair_hash<0>(hSplit, hCell) +
       HighsHashHelpers::pair_hash<1>(
           cell, currentPartitionLinks[cell] - splitPoint) +
       HighsHashHelpers::pair_hash<2>(splitPoint, splitPoint - cell)) >>
      32;

  // employ prefix pruning scheme as in bliss
  if (!firstLeaveCertificate.empty()) {
    firstLeavePrefixLen +=
        (firstLeavePrefixLen ==
         static_cast<HighsInt>(currNodeCertificate.size())) *
        (certificateVal == firstLeaveCertificate[currNodeCertificate.size()]);
    bestLeavePrefixLen +=
        (bestLeavePrefixLen ==
         static_cast<HighsInt>(currNodeCertificate.size())) *
        (certificateVal == bestLeaveCertificate[currNodeCertificate.size()]);

    // if the node certificate is not a prefix of the first leave's certificate
    // and it comes lexicographically after the certificate value of the
    // lexicographically smallest leave certificate we prune the node
    if (firstLeavePrefixLen <=
            static_cast<HighsInt>(currNodeCertificate.size()) &&
        bestLeavePrefixLen <=
            static_cast<HighsInt>(currNodeCertificate.size())) {
      u32 diffVal = bestLeavePrefixLen ==
                            static_cast<HighsInt>(currNodeCertificate.size())
                        ? certificateVal
                        : currNodeCertificate[bestLeavePrefixLen];
      if (diffVal > bestLeaveCertificate[bestLeavePrefixLen]) return false;
    }
  }

  currentPartitionLinks[splitPoint] = currentPartitionLinks[cell];
  currentPartitionLinks[cell] = splitPoint;
  cellCreationStack.push_back(splitPoint);
  currNodeCertificate.push_back(certificateVal);

  return true;
}

void HighsSymmetryDetection::markCellForRefinement(HighsInt cell) {
  if (cellSize(cell) == 1 || cellInRefinementQueue[cell]) return;

  cellInRefinementQueue[cell] = true;
  refinementQueue.push_back(cell);
  std::push_heap(refinementQueue.begin(), refinementQueue.end(),
                 std::greater<HighsInt>());
}

void HighsSymmetryDetection::clearRefinementState() {
  for (HighsInt c : refinementQueue) cellInRefinementQueue[c] = false;
  refinementQueue.clear();
  vertexHash.clear();
}

HighsSymmetryDetection::u64 HighsSymmetryDetection::getVertexHash(HighsInt v) {
  const u64* h = vertexHash.find(v);
  if (h) return *h;
  return 0;
}

bool HighsSymmetryDetection::partitionRefinement() {
  while (!refinementQueue.empty()) {
    std::pop_heap(refinementQueue.begin(), refinementQueue.end(),
                  std::greater<HighsInt>());

    HighsInt cellStart = refinementQueue.back();
    HighsInt firstCellStart = cellStart;
    refinementQueue.pop_back();
    cellInRefinementQueue[cellStart] = false;

    if (cellSize(cellStart) == 1) continue;
    HighsInt cellEnd = currentPartitionLinks[cellStart];
    assert(cellEnd >= cellStart);

    const size_t startCreationStackSize = cellCreationStack.size();
    const size_t startCertificateSize = currNodeCertificate.size();

    // first check which vertices do have updated hash values and put them to
    // the end of the partition
    HighsInt refineStart =
        std::partition(
            currentPartition.begin() + cellStart,
            currentPartition.begin() + cellEnd,
            [&](HighsInt v) { return vertexHash.find(v) == nullptr; }) -
        currentPartition.begin();

    // if there are none there is nothing to refine
    if (refineStart == cellEnd) continue;

    // sort the vertices that have updated hash values by their hash values
    pdqsort(currentPartition.begin() + refineStart,
            currentPartition.begin() + cellEnd, [&](HighsInt v1, HighsInt v2) {
              return vertexHash[v1] < vertexHash[v2];
            });

    // if not all vertices have updated hash values directly create the first
    // new cell at the start of the range that we want to refine
    if (refineStart != cellStart) {
      assert(refineStart != cellStart);
      assert(refineStart != cellEnd);
      if (!splitCell(cellStart, refineStart)) {
        // node can be pruned, make sure hash values are cleared and queue is
        // empty
        clearRefinementState();
        return false;
      }
      cellStart = refineStart;
      updateCellMembership(cellStart, cellStart);
    }

    // now update the remaining vertices
    bool prune = false;
    HighsInt i;
    assert(vertexHash.find(currentPartition[cellStart]) != nullptr);
    // store value of first hash
    u64 lastHash = vertexHash[currentPartition[cellStart]];
    for (i = cellStart + 1; i < cellEnd; ++i) {
      HighsInt vertex = currentPartition[i];
      // get this vertex hash value
      u64 hash = vertexHash[vertex];

      if (hash != lastHash) {
        // hash values do not match -> start of new cell
        if (!splitCell(cellStart, i)) {
          // refinement process yielded bad prefix of certificate
          // -> node can be pruned
          prune = true;
          break;
        }
        cellStart = i;
        // remember hash value of this new cell under lastHash
        lastHash = hash;
      }

      // update membership of vertex to new cell
      updateCellMembership(i, cellStart);
    }

    if (prune) {
      // node can be pruned, make sure hash values are cleared and queue is
      // empty
      clearRefinementState();
      currentPartitionLinks[firstCellStart] = cellEnd;

      // undo possibly incomplete changes done to the cells
      for (--i; i >= refineStart; --i)
        updateCellMembership(i, firstCellStart, false);

      cellCreationStack.resize(startCreationStackSize);
      currNodeCertificate.resize(startCertificateSize);
      if (!firstLeaveCertificate.empty()) {
        firstLeavePrefixLen = std::min(
            firstLeavePrefixLen, static_cast<HighsInt>(startCertificateSize));
        bestLeavePrefixLen = std::min(
            bestLeavePrefixLen, static_cast<HighsInt>(startCertificateSize));
      }

      return false;
    }

    assert(currentPartitionLinks[cellStart] == cellEnd);
  }

  vertexHash.clear();

  return true;
}

HighsInt HighsSymmetryDetection::selectTargetCell() const {
  HighsInt i = 0;
  if (nodeStack.size() > 1) i = nodeStack[nodeStack.size() - 2].targetCell;

  while (i < numVertices) {
    if (cellSize(i) > 1) return i;

    ++i;
  }

  return -1;
}

bool HighsSymmetryDetection::checkStoredAutomorphism(HighsInt vertex) const {
  HighsInt numCheck = std::min(numAutomorphisms, (HighsInt)64);

  for (HighsInt i = 0; i < numCheck; ++i) {
    const HighsInt* automorphism = automorphisms.data() + i * numVertices;
    bool automorphismUseful = true;
    for (HighsInt j = nodeStack.size() - 2; j >= firstPathDepth; --j) {
      HighsInt fixPos = vertexPosition[nodeStack[j].lastDistiguished];

      if (automorphism[fixPos] != vertexGroundSet[fixPos]) {
        automorphismUseful = false;
        break;
      }
    }

    if (!automorphismUseful) continue;

    if (automorphism[vertexPosition[vertex]] < vertex) return false;
  }

  return true;
}

bool HighsSymmetryDetection::determineNextToDistinguish() {
  Node& currNode = nodeStack.back();
  distinguishCands.clear();
  std::vector<HighsInt>::iterator cellStart;
  std::vector<HighsInt>::iterator cellEnd;
  cellStart = currentPartition.begin() + currNode.targetCell;
  cellEnd =
      currentPartition.begin() + currentPartitionLinks[currNode.targetCell];

  if (currNode.lastDistiguished == -1) {
    auto nextDistinguishPos = std::min_element(cellStart, cellEnd);
    distinguishCands.push_back(&*nextDistinguishPos);
  } else if ((HighsInt)nodeStack.size() > firstPathDepth) {
    for (auto i = cellStart; i != cellEnd; ++i) {
      if (*i > currNode.lastDistiguished && checkStoredAutomorphism(*i))
        distinguishCands.push_back(&*i);
    }
    if (distinguishCands.empty()) return false;
    auto nextDistinguishPos =
        std::min_element(distinguishCands.begin(), distinguishCands.end(),
                         [](HighsInt* a, HighsInt* b) { return *a < *b; });
    std::swap(*distinguishCands.begin(), *nextDistinguishPos);
    distinguishCands.resize(1);
  } else {
    for (auto i = cellStart; i != cellEnd; ++i) {
      if (*i > currNode.lastDistiguished && vertexGroundSet[getOrbit(*i)] == *i)
        distinguishCands.push_back(&*i);
    }
    if (distinguishCands.empty()) return false;
    auto nextDistinguishPos =
        std::min_element(distinguishCands.begin(), distinguishCands.end(),
                         [](HighsInt* a, HighsInt* b) { return *a < *b; });
    std::swap(*distinguishCands.begin(), *nextDistinguishPos);
    distinguishCands.resize(1);
  }

  return true;
}

bool HighsSymmetryDetection::distinguishVertex(HighsInt targetCell) {
  assert(distinguishCands.size() == 1u);
  HighsInt targetCellEnd = currentPartitionLinks[targetCell];
  HighsInt newCell = targetCellEnd - 1;
  std::swap(*distinguishCands[0], currentPartition[newCell]);
  nodeStack.back().lastDistiguished = currentPartition[newCell];

  if (!splitCell(targetCell, newCell)) return false;

  updateCellMembership(newCell, newCell);

  return true;
}

void HighsSymmetryDetection::backtrack(HighsInt backtrackStackNewEnd,
                                       HighsInt backtrackStackEnd) {
  // we assume that we always backtrack from a leaf node, i.e. a discrete
  // partition therefore we do not need to remember the values of the hash
  // contributions as it is the indentity for each position and all new cells
  // are on the cell creation stack.
  for (HighsInt stackPos = backtrackStackEnd - 1;
       stackPos >= backtrackStackNewEnd; --stackPos) {
    HighsInt cell = cellCreationStack[stackPos];
    // look up the cell start of the preceding cell with link compression
    HighsInt newStart = getCellStart(cell - 1);
    // remember the current end
    HighsInt currEnd = currentPartitionLinks[cell];
    // change the link to point to the start of the preceding cell
    currentPartitionLinks[cell] = newStart;
    // change the link of the start pointer of the preceding cell to point to
    // the end of this cell
    currentPartitionLinks[newStart] = currEnd;
  }
}

void HighsSymmetryDetection::cleanupBacktrack(HighsInt cellCreationStackPos) {
  // the links have been updated. Even though they might still not be fully
  // compressed the cell starts will all point to the correct cell end and the
  // lookup with path compression will give the correct start
  for (HighsInt stackPos = cellCreationStack.size() - 1;
       stackPos >= cellCreationStackPos; --stackPos) {
    HighsInt cell = cellCreationStack[stackPos];

    HighsInt cellStart = getCellStart(cell);
    HighsInt cellEnd = currentPartitionLinks[cellStart];

    for (HighsInt v = cell;
         v < cellEnd && vertexToCell[currentPartition[v]] == cell; ++v)
      updateCellMembership(v, cellStart, false);
  }

  cellCreationStack.resize(cellCreationStackPos);
}

HighsInt HighsSymmetryDetection::getCellStart(HighsInt pos) {
  HighsInt startPos = currentPartitionLinks[pos];
  if (startPos > pos) return pos;
  if (currentPartitionLinks[startPos] < startPos) {
    do {
      linkCompressionStack.push_back(pos);
      pos = startPos;
      startPos = currentPartitionLinks[startPos];
    } while (currentPartitionLinks[startPos] < startPos);

    do {
      currentPartitionLinks[linkCompressionStack.back()] = startPos;
      linkCompressionStack.pop_back();
    } while (!linkCompressionStack.empty());
  }

  return startPos;
}

void HighsSymmetryDetection::createNode() {
  nodeStack.emplace_back();
  nodeStack.back().stackStart = cellCreationStack.size();
  nodeStack.back().certificateEnd = currNodeCertificate.size();
  nodeStack.back().targetCell = -1;
  nodeStack.back().lastDistiguished = -1;
}

struct MatrixColumn {
  uint32_t cost;
  uint32_t lb;
  uint32_t ub;
  uint32_t integral;
  uint32_t len;

  bool operator==(const MatrixColumn& other) const {
    return std::memcmp(this, &other, sizeof(MatrixColumn)) == 0;
  }
};

struct MatrixRow {
  uint32_t lb;
  uint32_t ub;
  uint32_t len;

  bool operator==(const MatrixRow& other) const {
    return std::memcmp(this, &other, sizeof(MatrixRow)) == 0;
  }
};

void HighsSymmetryDetection::loadModelAsGraph(const HighsLp& model,
                                              double epsilon) {
  this->model = &model;
  numCol = model.num_col_;
  numRow = model.num_row_;
  numVertices = numRow + numCol;

  cellInRefinementQueue.resize(numVertices);
  vertexToCell.resize(numVertices);
  refinementQueue.reserve(numVertices);
  currNodeCertificate.reserve(numVertices);

  HighsHashTable<MatrixColumn, HighsInt> columnSet;
  HighsHashTable<MatrixRow, HighsInt> rowSet;
  HighsMatrixColoring coloring(epsilon);
  edgeBuffer.resize(numVertices);
  // set up row and column based incidence matrix
  HighsInt numNz = model.a_matrix_.index_.size();
  Gedge.resize(2 * numNz);
  std::transform(model.a_matrix_.index_.begin(), model.a_matrix_.index_.end(),
                 Gedge.begin(), [&](HighsInt rowIndex) {
                   return std::make_pair(rowIndex + numCol, HighsUInt{0});
                 });

  Gstart.resize(numVertices + 1);
  std::copy(model.a_matrix_.start_.begin(), model.a_matrix_.start_.end(),
            Gstart.begin());

  // set up the column colors and count row sizes
  std::vector<HighsInt> rowSizes(numRow);
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Gstart[i]; j < Gstart[i + 1]; ++j) {
      Gedge[j].second = coloring.color(model.a_matrix_.value_[j]);
      rowSizes[model.a_matrix_.index_[j]] += 1;
    }
  }

  // next set up the row starts using the computed row sizes
  HighsInt offset = numNz;
  for (HighsInt i = 0; i < numRow; ++i) {
    Gstart[numCol + i] = offset;
    offset += rowSizes[i];
  }
  Gstart[numCol + numRow] = offset;

  Gend.assign(Gstart.begin() + 1, Gstart.end());

  // finally add the nonzeros to the row major matrix
  for (HighsInt i = 0; i < numCol; ++i) {
    for (HighsInt j = Gstart[i]; j < Gstart[i + 1]; ++j) {
      HighsInt row = model.a_matrix_.index_[j];
      HighsInt ARpos = Gstart[numCol + row + 1] - rowSizes[row];
      rowSizes[row] -= 1;
      Gedge[ARpos].first = i;
      Gedge[ARpos].second = Gedge[j].second;
    }
  }

  // loop over the columns and assign them a number that is distinct based on
  // their upper/lower bounds, cost, and integrality status. Use the columnSet
  // hash table to look up whether a column with similar properties exists and
  // use the previous number in that case. The number is stored in the
  // colToCell array which is subsequently used to sort an initial column
  // permutation.
  HighsInt indexOffset = numCol + 1;
  for (HighsInt i = 0; i < numCol; ++i) {
    MatrixColumn matrixCol;

    matrixCol.cost = coloring.color(model.col_cost_[i]);
    matrixCol.lb = coloring.color(model.col_lower_[i]);
    matrixCol.ub = coloring.color(model.col_upper_[i]);
    matrixCol.integral = (u32)model.integrality_[i];
    matrixCol.len = Gstart[i + 1] - Gstart[i];

    HighsInt* columnCell = &columnSet[matrixCol];

    if (*columnCell == 0) {
      *columnCell = columnSet.size();
      if (model.col_lower_[i] != 0.0 || model.col_upper_[i] != 1.0 ||
          model.integrality_[i] == HighsVarType::kContinuous)
        *columnCell += indexOffset;
    }

    vertexToCell[i] = *columnCell;
  }

  indexOffset = 2 * numCol + 1;
  for (HighsInt i = 0; i < numRow; ++i) {
    MatrixRow matrixRow;

    matrixRow.lb = coloring.color(model.row_lower_[i]);
    matrixRow.ub = coloring.color(model.row_upper_[i]);
    matrixRow.len = Gstart[numCol + i + 1] - Gstart[numCol + i];

    HighsInt* rowCell = &rowSet[matrixRow];

    if (*rowCell == 0) *rowCell = rowSet.size();

    vertexToCell[numCol + i] = indexOffset + *rowCell;
  }

  // set up the initial partition array, sort by the colToCell value
  // assigned above
  currentPartition.resize(numVertices);
  std::iota(currentPartition.begin(), currentPartition.end(), 0);
  pdqsort(currentPartition.begin(), currentPartition.end(),
          [&](HighsInt v1, HighsInt v2) {
            return vertexToCell[v1] < vertexToCell[v2];
          });

  // now set up partition links and correct the colToCell array to the
  // correct cell index
  currentPartitionLinks.resize(numVertices);
  HighsInt cellStart = 0;
  HighsInt cellNumber = 0;
  for (HighsInt i = 0; i < numVertices; ++i) {
    HighsInt vertex = currentPartition[i];
    // if the cell number is different to the current cell number this is the
    // start of a new cell
    if (cellNumber != vertexToCell[vertex]) {
      // remember the number of this cell to identify its end
      cellNumber = vertexToCell[vertex];
      // set the link of the cell start to point to its end
      currentPartitionLinks[cellStart] = i;
      // remember start of this cell
      cellStart = i;
    }

    // correct the colToCell array to not store the start index of the
    // cell, not its number
    vertexToCell[vertex] = cellStart;
    // set the link of the column to the cellStart
    currentPartitionLinks[i] = cellStart;
  }

  // set the column partition link of the last started cell to point past the
  // end
  currentPartitionLinks[cellStart] = numVertices;
}

HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>>
HighsSymmetryDetection::dumpCurrentGraph() {
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> graphTriplets;

  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt colCell = vertexToCell[i];
    for (HighsInt j = Gstart[i]; j != Gend[i]; ++j)
      graphTriplets.insert(vertexToCell[Gedge[j].first], colCell,
                           Gedge[j].second);
    for (HighsInt j = Gend[i]; j != Gstart[i + 1]; ++j)
      graphTriplets.insert(Gedge[j].first, colCell, Gedge[j].second);
  }

  return graphTriplets;
}

void HighsSymmetryDetection::switchToNextNode(HighsInt backtrackDepth) {
  HighsInt stackEnd = cellCreationStack.size();
  // we need to backtrack the datastructures
  nodeStack.resize(backtrackDepth);
  if (backtrackDepth == 0) return;
  do {
    Node& currNode = nodeStack.back();
    backtrack(currNode.stackStart, stackEnd);
    stackEnd = currNode.stackStart;
    firstPathDepth = std::min((HighsInt)nodeStack.size(), firstPathDepth);
    bestPathDepth = std::min((HighsInt)nodeStack.size(), bestPathDepth);
    firstLeavePrefixLen =
        std::min(currNode.certificateEnd, firstLeavePrefixLen);
    bestLeavePrefixLen = std::min(currNode.certificateEnd, bestLeavePrefixLen);
    currNodeCertificate.resize(currNode.certificateEnd);
    if (!determineNextToDistinguish()) {
      nodeStack.pop_back();
      continue;
    }

    // call cleanup backtrack with the final stackEnd
    // so that all hashes are up to date and the link arrays do not contain
    // chains anymore
    cleanupBacktrack(stackEnd);
    HighsInt targetCell = currNode.targetCell;

    if (!distinguishVertex(targetCell)) {
      // if distinguishing the next vertex fails, it means that its certificate
      // value is lexicographically larger than that of the best leave
      nodeStack.pop_back();
      continue;
    }

    if (!partitionRefinement()) {
      stackEnd = cellCreationStack.size();
      continue;
    }

    createNode();
    break;
  } while (!nodeStack.empty());
}

bool HighsSymmetryDetection::compareCurrentGraph(
    const HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>>& otherGraph,
    HighsInt& wrongCell) const {
  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt colCell = vertexToCell[i];

    for (HighsInt j = Gstart[i]; j != Gend[i]; ++j)
      if (!otherGraph.find(std::make_tuple(vertexToCell[Gedge[j].first],
                                           colCell, Gedge[j].second))) {
        // return which cell does not match in its neighbourhood as this should
        // have been detected with the hashing it can very rarely happen due to
        // a hash collision. In such a case we want to backtrack to the last
        // time where we targeted this particular cell. Otherwise we could spent
        // a long time searching for a matching leave value until every
        // combination is exhausted and for each leave in this subtree the graph
        // comparison will fail on this edge.
        wrongCell = colCell;
        return false;
      }
    for (HighsInt j = Gend[i]; j != Gstart[i + 1]; ++j)
      if (!otherGraph.find(
              std::make_tuple(Gedge[j].first, colCell, Gedge[j].second))) {
        wrongCell = colCell;
        return false;
      }
  }

  return true;
}

bool HighsSymmetryDetection::isFromBinaryColumn(HighsInt pos) const {
  if (pos >= numActiveCols) return false;

  HighsInt col = currentPartition[pos];

  if (model->col_lower_[col] != 0.0 || model->col_upper_[col] != 1.0 ||
      model->integrality_[col] == HighsVarType::kContinuous)
    return false;

  return true;
}

HighsSymmetryDetection::ComponentData
HighsSymmetryDetection::computeComponentData(
    const HighsSymmetries& symmetries) {
  ComponentData componentData;

  componentData.components.reset(numActiveCols);
  componentData.firstUnfixed.assign(symmetries.numPerms, -1);
  componentData.numUnfixed.assign(symmetries.numPerms, 0);
  for (HighsInt i = 0; i < symmetries.numPerms; ++i) {
    const HighsInt* perm = symmetries.permutations.data() + i * numActiveCols;

    for (HighsInt j = 0; j < numActiveCols; ++j) {
      if (perm[j] != vertexGroundSet[j]) {
        HighsInt pos = vertexPosition[perm[j]];
        componentData.numUnfixed[i] += 1;
        if (componentData.firstUnfixed[i] != -1)
          componentData.components.merge(componentData.firstUnfixed[i], pos);
        else
          componentData.firstUnfixed[i] = pos;
      }
    }
  }

  componentData.componentSets.assign(vertexGroundSet.begin(),
                                     vertexGroundSet.begin() + numActiveCols);
  pdqsort(componentData.componentSets.begin(),
          componentData.componentSets.end(), [&](HighsInt u, HighsInt v) {
            HighsInt uComp = componentData.components.getSet(vertexPosition[u]);
            HighsInt vComp = componentData.components.getSet(vertexPosition[v]);
            return std::make_pair(
                       componentData.components.getSetSize(uComp) == 1, uComp) <
                   std::make_pair(
                       componentData.components.getSetSize(vComp) == 1, vComp);
          });

  HighsInt currentComponentStart = -1;
  HighsInt currentComponent = -1;
  HighsHashTable<HighsInt> currComponentOrbits;
  for (HighsInt i = 0; i < numActiveCols; ++i) {
    HighsInt comp = componentData.components.getSet(
        vertexPosition[componentData.componentSets[i]]);
    if (componentData.components.getSetSize(comp) == 1) break;
    if (comp != currentComponent) {
      currentComponent = comp;
      currentComponentStart = i;
      componentData.componentStarts.push_back(currentComponentStart);
      componentData.componentNumber.push_back(currentComponent);
      componentData.componentNumOrbits.emplace_back();
      currComponentOrbits.clear();
    }

    if (currComponentOrbits.insert(getOrbit(componentData.componentSets[i]))) {
      ++componentData.componentNumOrbits.back();
    }
  }

  componentData.permComponents.reserve(symmetries.numPerms);

  for (HighsInt i = 0; i < symmetries.numPerms; ++i) {
    if (componentData.firstUnfixed[i] == -1) continue;
    componentData.permComponents.push_back(i);
  }

  pdqsort(componentData.permComponents.begin(),
          componentData.permComponents.end(), [&](HighsInt i, HighsInt j) {
            HighsInt seti =
                componentData.components.getSet(componentData.firstUnfixed[i]);
            HighsInt setj =
                componentData.components.getSet(componentData.firstUnfixed[j]);
            return std::make_pair(seti, componentData.numUnfixed[i]) <
                   std::make_pair(setj, componentData.numUnfixed[j]);
          });

  currentComponentStart = -1;
  currentComponent = -1;

  HighsInt numUsedPerms = componentData.permComponents.size();

  for (HighsInt i = 0; i < numUsedPerms; ++i) {
    HighsInt p = componentData.permComponents[i];
    HighsInt comp =
        componentData.components.getSet(componentData.firstUnfixed[p]);
    if (comp != currentComponent) {
      currentComponent = comp;
      currentComponentStart = i;
      componentData.permComponentStarts.push_back(currentComponentStart);
    }
  }

  assert(componentData.permComponentStarts.size() ==
         componentData.componentStarts.size());
  componentData.permComponentStarts.push_back(numUsedPerms);

  // HighsInt numComponents = componentData.componentStarts.size();
  // printf("found %d components\n", numComponents);
  componentData.componentStarts.push_back(numActiveCols);

  return componentData;
}

bool HighsSymmetryDetection::isFullOrbitope(const ComponentData& componentData,
                                            HighsInt component,
                                            HighsSymmetries& symmetries) {
  HighsInt componentSize = componentData.componentStarts[component + 1] -
                           componentData.componentStarts[component];
  if (componentSize == 1) return false;

  // check that component acts only on binary variables
  for (HighsInt i = componentData.componentStarts[component];
       i < componentData.componentStarts[component + 1]; ++i) {
    HighsInt col = componentData.componentSets[i];
    if (model->integrality_[col] == HighsVarType::kContinuous ||
        model->col_lower_[col] != 0.0 || model->col_upper_[col] != 1.0)
      return false;
  }

  // check that the number of unfixed variables in the first permutation is even
  HighsInt p0 =
      componentData
          .permComponents[componentData.permComponentStarts[component]];
  if (componentData.numUnfixed[p0] & 1) return false;

  // check that the other permutations in the component have the same number of
  // unfixed variables as the first one
  for (HighsInt k = componentData.permComponentStarts[component] + 1;
       k < componentData.permComponentStarts[component + 1]; ++k) {
    HighsInt p = componentData.permComponents[k];
    if (componentData.numUnfixed[p] != componentData.numUnfixed[p0]) {
      // printf("number of unfixed cols in perm %d: %d\n", p,
      //       componentData.numUnfixed[p]);
      // printf("wrong number of unfixed columns in permutation\n");
      return false;
    }
  }

  // all unfixed variables in the permutation must be part of a two cycle
  // and each two cycle defines one row of the orbitope. Hence the number
  // of unfixed variables in each permutation must be even and the number of
  // rows is the number of unfixed variables divided by two
  HighsInt orbitopeNumRows = componentData.numUnfixed[p0] >> 1;

  // if this component is a full orbitope the component size must be
  // divisible by the number of orbits and the result is the size of each
  // orbit
  HighsInt orbitopeOrbitSize = componentSize / orbitopeNumRows;
  if (orbitopeOrbitSize * orbitopeNumRows != componentSize) {
    // printf("wrong number of orbits (%d orbits for component of size %d)\n",
    //        orbitopeNumRows, componentSize);
    return false;
  }

  // the number of permutations must be n-1 where n is the size of the
  // orbits in the orbitope
  HighsInt componentNumPerms =
      componentData.permComponentStarts[component + 1] -
      componentData.permComponentStarts[component];
  if (componentNumPerms != (orbitopeOrbitSize - 1)) {
    // printf("wrong number of perms\n");
    return false;
  }

  // set up the first two columns of the orbitope matrix based on the first
  // permutation.
  HighsOrbitopeMatrix orbitopeMatrix;
  orbitopeMatrix.matrix.resize(componentSize, -1);
  orbitopeMatrix.numRows = orbitopeNumRows;
  orbitopeMatrix.rowLength = orbitopeOrbitSize;
  assert(componentSize == orbitopeMatrix.numRows * orbitopeMatrix.rowLength);

  const HighsInt* perm = symmetries.permutations.data() + p0 * numActiveCols;
  HighsHashTable<HighsInt> colSet;
  HighsInt m = 0;
  for (HighsInt j = 0; j < numActiveCols; ++j) {
    HighsInt jImagePos = vertexPosition[perm[j]];
    if (jImagePos <= j) continue;
    if (m == orbitopeNumRows) return false;

    // permutation should consist of two cycles
    if (perm[jImagePos] != vertexGroundSet[j]) return false;

    orbitopeMatrix.matrix[m] = vertexGroundSet[j];
    orbitopeMatrix.matrix[orbitopeMatrix.numRows + m] = perm[j];

    // Remember set of variables of the orbtiope matrix. Each variable should
    // occur only once, otherwise the permutation do not work out. Since we

    if (!colSet.insert(vertexGroundSet[j])) return false;
    if (!colSet.insert(perm[j])) return false;
    ++m;
  }

  // printf("set up first two columns of possible orbitope with permutation
  // %d\n",
  //        p0);

  HighsInt numColsAdded = 2;
  bool triedLeftExtension = false;

  while (numColsAdded < orbitopeMatrix.rowLength) {
    if (colSet.size() != numColsAdded * orbitopeMatrix.numRows) return false;

    HighsInt* thisCol = &orbitopeMatrix(0, numColsAdded);
    HighsInt* prevCol = &orbitopeMatrix(0, numColsAdded - 1);

    bool foundCand = false;

    while (true) {
      HighsInt movePos = vertexPosition[prevCol[0]];
      perm = nullptr;

      for (HighsInt k = componentData.permComponentStarts[component] + 1;
           k < componentData.permComponentStarts[component + 1]; ++k) {
        HighsInt p = componentData.permComponents[k];
        perm = symmetries.permutations.data() + p * numActiveCols;

        if (perm[movePos] != vertexGroundSet[movePos] &&
            !colSet.find(perm[movePos])) {
          foundCand = true;
          break;
        }
      }

      if (!foundCand && !triedLeftExtension) {
        // if we fail to find a permutation extending the last column directly
        // after we set up the first two columns we might try to extend column
        // zero instead of the previous one
        prevCol = &orbitopeMatrix(0, 0);
        triedLeftExtension = true;
        continue;
      }

      break;
    }

    if (!foundCand) {
      // printf("did not find next permutation moving col %d\n", prevCol[0]);
      return false;
    }

    for (HighsInt j = 0; j < orbitopeMatrix.numRows; ++j) {
      HighsInt nextVertex = perm[vertexPosition[prevCol[j]]];

      thisCol[j] = nextVertex;

      // check if this is a two cycle
      if (perm[vertexPosition[nextVertex]] != prevCol[j]) return false;

      if (!colSet.insert(thisCol[j])) {
        // printf("col already exists\n");
        return false;
      }
    }

    ++numColsAdded;
  }

  if (colSet.size() != componentSize) {
    // printf("not all columns of component are mapped\n");
    return false;
  }

  for (HighsInt col : orbitopeMatrix.matrix)
    symmetries.columnToOrbitope.insert(col, symmetries.orbitopes.size());

  symmetries.orbitopes.emplace_back(std::move(orbitopeMatrix));

  // printf("component %d is full orbitope: size %d and %d orbits\n", component,
  //        componentSize, componentData.componentNumOrbits[component]);
  return true;
}

bool HighsSymmetryDetection::initializeDetection() {
  initializeHashValues();
  partitionRefinement();
  removeFixPoints();
  if (numActiveCols == 0) return false;
  return true;
}

bool HighsSymmetryDetection::recordAutomorphism(
    const std::vector<HighsInt>& leavePartition, HighsSymmetries& symmetries,
    HighsInt maxPerms, HighsInt pathDepth, HighsInt& backtrackDepth) {
  HighsInt k = (numAutomorphisms++) & 63;
  HighsInt* permutation = automorphisms.data() + k * numVertices;
  for (HighsInt i = 0; i < numVertices; ++i)
    permutation[vertexPosition[currentPartition[i]]] = leavePartition[i];

  bool report = false;
  for (HighsInt i = 0; i < numVertices; ++i) {
    if (mergeOrbits(permutation[i], vertexGroundSet[i]) && i < numActiveCols) {
      assert(permutation[i] < numCol);
      report = true;
    }
  }

  if (report) {
    symmetries.permutations.insert(symmetries.permutations.end(), permutation,
                                   permutation + numActiveCols);
    ++symmetries.numPerms;
    if (symmetries.numPerms == maxPerms) return true;
  }

  backtrackDepth = std::min(backtrackDepth, pathDepth);
  return false;
}

void HighsSymmetryDetection::run(HighsSymmetries& symmetries,
                                 const HighsTimer* timer, double timeLimit) {
  assert(numActiveCols != 0);
  initializeGroundSet();
  currNodeCertificate.clear();
  cellCreationStack.clear();
  createNode();
  HighsInt maxPerms = 64000000 / numActiveCols;
  HighsSplitDeque* workerDeque = HighsTaskExecutor::getThisWorkerDeque();
  while (!nodeStack.empty()) {
    // stop if the task was cancelled, e.g. if the MIP solver returned before
    // finishing symmetry detection, or if the time limit is reached
    workerDeque->checkInterrupt();
    if (timer && timer->read() >= timeLimit) break;
    HighsInt targetCell = selectTargetCell();
    if (targetCell == -1) {
      if (firstLeavePartition.empty()) {
        firstLeavePartition = currentPartition;
        firstLeaveCertificate = currNodeCertificate;
        bestLeaveCertificate = currNodeCertificate;
        firstLeaveGraph = dumpCurrentGraph();
        firstPathDepth = nodeStack.size();
        bestPathDepth = nodeStack.size();
        firstLeavePrefixLen = currNodeCertificate.size();
        bestLeavePrefixLen = currNodeCertificate.size();

        HighsInt backtrackDepth = firstPathDepth - 1;
        while (backtrackDepth > 0 &&
               !isFromBinaryColumn(nodeStack[backtrackDepth - 1].targetCell))
          --backtrackDepth;
        switchToNextNode(backtrackDepth);
      } else {
        HighsInt wrongCell = -1;
        HighsInt backtrackDepth = nodeStack.size() - 1;
        assert(currNodeCertificate.size() == firstLeaveCertificate.size());
        if (firstLeavePrefixLen ==
                static_cast<HighsInt>(currNodeCertificate.size()) ||
            bestLeavePrefixLen ==
                static_cast<HighsInt>(currNodeCertificate.size())) {
          if (firstLeavePrefixLen ==
                  static_cast<HighsInt>(currNodeCertificate.size()) &&
              compareCurrentGraph(firstLeaveGraph, wrongCell)) {
            if (recordAutomorphism(firstLeavePartition, symmetries, maxPerms,
                                   firstPathDepth, backtrackDepth))
              break;
          } else if (!bestLeavePartition.empty() &&
                     bestLeavePrefixLen ==
                         static_cast<HighsInt>(currNodeCertificate.size()) &&
                     compareCurrentGraph(bestLeaveGraph, wrongCell)) {
            if (recordAutomorphism(bestLeavePartition, symmetries, maxPerms,
                                   bestPathDepth, backtrackDepth))
              break;
          } else if (bestLeavePrefixLen <
                         static_cast<HighsInt>(currNodeCertificate.size()) &&
                     currNodeCertificate[bestLeavePrefixLen] >
                         bestLeaveCertificate[bestLeavePrefixLen]) {
            // certificate value is lexicographically above the smallest one
            // seen so far, so we might be able to backtrack to a higher level
            HighsInt possibleBacktrackDepth = firstPathDepth - 1;
            while (nodeStack[possibleBacktrackDepth].certificateEnd <=
                   bestLeavePrefixLen)
              ++possibleBacktrackDepth;

            backtrackDepth = std::min(possibleBacktrackDepth, backtrackDepth);
          } else {
            // This case can be caused by a hash collision which was now
            // detected in the graph comparison call. The graph comparison call
            // will return the cell where the vertex neighbourhood caused a
            // mismatch on the edges. This would have been detected by
            // an exact partition refinement when we targeted that cell the last
            // time, so that is where we can backtrack to.
            HighsInt possibleBacktrackDepth;
            for (possibleBacktrackDepth = backtrackDepth;
                 possibleBacktrackDepth >= 0; --possibleBacktrackDepth) {
              if (nodeStack[possibleBacktrackDepth].targetCell == wrongCell) {
                backtrackDepth = possibleBacktrackDepth;
                break;
              }
            }
          }
        } else {
          // leave must have a lexicographically smaller certificate value
          // than the current best leave, because its prefix length is smaller
          // than the best leaves and it would have been already pruned if
          // it's certificate value was larger unless it is equal to the first
          // leaf nodes certificate value which is caught by the first case
          // of the if condition. Hence, having a lexicographically smaller
          // certificate value than the best leave is the only way to get
          // here.
          assert(bestLeaveCertificate[bestLeavePrefixLen] >
                     currNodeCertificate[bestLeavePrefixLen] &&
                 std::memcmp(bestLeaveCertificate.data(),
                             currNodeCertificate.data(),
                             bestLeavePrefixLen * sizeof(u32)) == 0);
          bestLeaveCertificate = currNodeCertificate;
          bestLeaveGraph = dumpCurrentGraph();
          bestLeavePartition = currentPartition;
          bestPathDepth = nodeStack.size();
          bestLeavePrefixLen = currNodeCertificate.size();
        }

        switchToNextNode(backtrackDepth);
      }
    } else {
      Node& currNode = nodeStack.back();
      currNode.targetCell = targetCell;
      bool success = determineNextToDistinguish();
      assert(success);
      if (!distinguishVertex(targetCell)) {
        switchToNextNode(nodeStack.size() - 1);
        continue;
      }
      if (!partitionRefinement()) {
        switchToNextNode(nodeStack.size());
        continue;
      }

      createNode();
    }
  }

  symmetries.numGenerators = symmetries.numPerms;
  if (symmetries.numPerms > 0) {
    vertexPosition.resize(numCol);

    ComponentData componentData = computeComponentData(symmetries);
    HighsInt numComponents = componentData.numComponents();

    for (HighsInt i = 0; i < numComponents; ++i) {
      if (componentData.componentSize(i) == 1) continue;

      isFullOrbitope(componentData, i, symmetries);
    }

    HighsHashTable<HighsInt> deletedPerms;
    for (HighsInt p = 0; p < symmetries.numPerms; ++p) {
      HighsInt* perm = symmetries.permutations.data() + p * numActiveCols;
      for (HighsInt i = 0; i < numActiveCols; ++i) {
        if (perm[i] != vertexGroundSet[i] &&
            symmetries.columnToOrbitope.find(perm[i])) {
          deletedPerms.insert(p);
          break;
        }
      }
    }

    // check which columns have non-trivial orbits
    HighsInt numFixed = 0;
    for (HighsInt i = 0; i < numActiveCols; ++i) {
      if (orbitSize[getOrbit(vertexGroundSet[i])] == 1 ||
          symmetries.columnToOrbitope.find(vertexGroundSet[i])) {
        vertexPosition[vertexGroundSet[i]] = -1;
        vertexGroundSet[i] = -1;
        numFixed += 1;
      }
    }

    if (numFixed != 0) {
      // now compress symmetries and the groundset to only contain the unfixed
      // columns and columns not handled by full orbitopes
      HighsInt p = 0;
      HighsInt* perms = symmetries.permutations.data();
      HighsInt* permEnd =
          symmetries.permutations.data() + symmetries.numPerms * numActiveCols;
      HighsInt* permOutput = symmetries.permutations.data();
      while (perms != permEnd) {
        if (!deletedPerms.find(p)) {
          for (HighsInt i = 0; i < numActiveCols; ++i) {
            if (vertexGroundSet[i] == -1) continue;

            *permOutput = perms[i];
            ++permOutput;
          }
        } else {
          --symmetries.numPerms;
        }
        perms += numActiveCols;
        ++p;
      }

      HighsInt outPos = 0;
      for (HighsInt i = 0; i < numActiveCols; ++i) {
        if (vertexGroundSet[i] == -1) continue;

        vertexGroundSet[outPos] = vertexGroundSet[i];
        vertexPosition[vertexGroundSet[outPos]] = outPos;
        outPos += 1;
      }

      numActiveCols -= numFixed;
      assert(permOutput == symmetries.permutations.data() +
                               symmetries.numPerms * numActiveCols);
    }

    vertexGroundSet.resize(numActiveCols);
    symmetries.permutationColumns = std::move(vertexGroundSet);
    symmetries.columnPosition = std::move(vertexPosition);
    symmetries.permutations.resize(symmetries.numPerms * numActiveCols);
  }
}
