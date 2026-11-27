/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsMachineSchedSeparator.cpp
 */

#include "mip/HighsMachineSchedSeparator.h"

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

bool HighsMachineSchedSeparator::findSingleMachineScheduleClique(
    std::vector<std::vector<double>>& vals,
    std::vector<std::vector<HighsInt>>& inds, std::vector<double>& rhss,
    const HighsDomain& globaldom, const HighsMipSolver& mipsolver) {
  enum class ArcType {
    kImplicationWhenOne,
    kImplicationWhenZero,
  };
  HighsInt largestDegree = 0;
  HighsInt largestDegreeCol = -1;
  HighsInt numRows = 0;
  std::vector<HighsInt> degrees(mipsolver.numCol());
  const HighsInt maxRows = std::min(HighsInt{50000}, 2 * mipsolver.numRow());
  // keys are (j,i) to values (p_ji, y_ji, implication-when-one-or-zero)
  // entries map: y_ji = val -> s_i >= s_j + p_ji - M * val
  HighsHashTable<std::pair<HighsInt, HighsInt>,
                 std::tuple<double, HighsInt, ArcType>>
      adjacency(maxRows + 2);
  // Used to track binaries that imply the same order in both directions
  HighsHashTable<std::pair<HighsInt, HighsInt>, std::vector<HighsInt>> arcToBin;
  // The keys in this map are triples of column indices representing the arc
  // endpoints (continuous variables) i, j and the binary variable y.
  // The values are bit-fields representing which values of y imply
  // which dependency relationships between i and j, as follows:
  // 1 -> (i,j) y = 0, 2 -> (i,j) y = 1, 4 -> (j,i) y = 0, 8 -> (j,i) y = 1
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsInt>, uint8_t> jobOrder;

  auto addEntry = [&](HighsInt posCol, HighsInt negCol, HighsInt binCol,
                      double p, ArcType t) {
    // My_ji + s_i - s_j >= p_ji
    // y_ji = val -> s_i >= s_j + p_ji - M * val
    // Make an arc from negCol (j) to posCol (i)
    if (p < 0) return;
    const auto it = adjacency.find({negCol, posCol});
    if (it != nullptr) {
      if (std::get<0>(*it) < p) {
        adjacency[{negCol, posCol}] = std::make_tuple(p, binCol, t);
      }
    } else {
      ++degrees[posCol];
      if (degrees[posCol] > largestDegree ||
          (degrees[posCol] == largestDegree && posCol < largestDegreeCol)) {
        largestDegreeCol = posCol;
        largestDegree = degrees[posCol];
      }
      adjacency.insert(std::make_pair(negCol, posCol),
                       std::make_tuple(p, binCol, t));
    }
    // Store binaries to later check if any two binaries
    // imply the same ordering of two jobs.
    HighsInt u = std::min(negCol, posCol);
    HighsInt v = std::max(negCol, posCol);
    if (jobOrder.find({u, v, binCol}) == nullptr) {
      jobOrder[{u, v, binCol}] = 0;
    }
    if (negCol < posCol) {
      arcToBin[{negCol, posCol}].emplace_back(binCol);
      jobOrder[{negCol, posCol, binCol}] |=
          t == ArcType::kImplicationWhenOne ? 2 : 1;
    } else {
      jobOrder[{posCol, negCol, binCol}] |=
          t == ArcType::kImplicationWhenOne ? 8 : 4;
    }
  };

  auto tryAddArc = [&](const HighsInt posCol, const HighsInt negCol,
                       const HighsInt binCol, const double rhs_0,
                       const double rhs_1) {
    if (rhs_0 > 0 || rhs_1 > 0) {
      if (rhs_0 > rhs_1) {
        addEntry(posCol, negCol, binCol, rhs_0, ArcType::kImplicationWhenZero);
      } else {
        addEntry(posCol, negCol, binCol, rhs_1, ArcType::kImplicationWhenOne);
      }
      ++numRows;
    }
  };

  for (HighsInt row = 0; row != mipsolver.numRow(); row++) {
    const double rowLower = mipsolver.model_->row_lower_[row];
    const double rowUpper = mipsolver.model_->row_upper_[row];
    if (rowLower == rowUpper) continue;
    const HighsInt start = mipsolver.mipdata_->ARstart_[row];
    const HighsInt end = mipsolver.mipdata_->ARstart_[row + 1];
    if (end - start != 3) continue;
    bool machineSchedRow = true;
    HighsInt posContCol = -1;
    HighsInt negContCol = -1;
    HighsInt binCol = -1;
    double binCoef = 0;
    for (HighsInt i = start; i != end; i++) {
      HighsInt col = mipsolver.mipdata_->ARindex_[i];
      if (globaldom.col_lower_[col] == -kHighsInf) {
        machineSchedRow = false;
        break;
      }
      if (globaldom.isBinary(col)) {
        if (binCol != -1) {
          machineSchedRow = false;
          break;
        }
        binCol = col;
        binCoef = mipsolver.mipdata_->ARvalue_[i];
      } else if (mipsolver.mipdata_->ARvalue_[i] == -1) {
        if (negContCol != -1) {
          machineSchedRow = false;
          break;
        }
        negContCol = col;
      } else if (mipsolver.mipdata_->ARvalue_[i] == 1) {
        if (posContCol != -1) {
          machineSchedRow = false;
          break;
        }
        posContCol = col;
      } else {
        machineSchedRow = false;
        break;
      }
    }
    if (!machineSchedRow || binCol == -1 || negContCol == -1 ||
        posContCol == -1)
      continue;
    // We want to put the row into form:
    // My_ji + s_i - s_j >= p_ji, p_ji >= 0
    // y_ji = 1 -> s_i >= s_j + p_ji - M
    if (rowUpper != kHighsInf) {
      // Given My_ij + s_i - s_j <= d
      // The row becomes (after multiplying by -1):
      // -My_ij + s_j - s_i >= -d
      // Add implication s_j >= s_i + p_ij + M, p_ij + M > 0 when binCol = 1
      // Add implication s_j >= s_i + p_ij, p_ij > 0 when binCol = 0
      tryAddArc(negContCol, posContCol, binCol, -rowUpper, -rowUpper + binCoef);
    }
    if (rowLower != -kHighsInf) {
      // Given My_ij + s_i - s_j >= d
      // Add implication s_i >= s_j + p_ji - M, p_ji - M > 0 when binCol = 1
      // Add implication s_i >= s_j + p_ji, p_ji > 0 when binCol = 0
      tryAddArc(posContCol, negContCol, binCol, rowLower, rowLower - binCoef);
    }
    if (numRows >= maxRows) break;
  }

  // Extract binary variables that imply the same job ordering
  if (!mipsolver.mipdata_->parallelLockActive()) {
    std::vector<HighsCliqueTable::CliqueVar> clique(2);
    for (const auto& entry : arcToBin) {
      std::pair<HighsInt, HighsInt> arc = entry.key();
      HighsInt baseBinCol = -1;
      bool baseArcImplicationWhenOne = false;
      for (const HighsInt& binCol : entry.value()) {
        bool arcImplicationWhenOne = false;
        bool impliesOrder = false;
        if ((jobOrder[{arc.first, arc.second, binCol}] & 9) == 9) {
          impliesOrder = true;
        } else if ((jobOrder[{arc.first, arc.second, binCol}] & 6) == 6) {
          impliesOrder = true;
          arcImplicationWhenOne = true;
        }
        if (!impliesOrder || binCol == baseBinCol) continue;
        if (baseBinCol == -1) {
          baseBinCol = binCol;
          baseArcImplicationWhenOne = arcImplicationWhenOne;
          continue;
        }
        // Add cliques x1 + ~x2 <= 1 and ~x1 + x2 <= 1 which together imply x1
        // == x2 (depending on signs may also have x1 == ~x2)
        HighsCliqueTable::CliqueVar stayCliqueVar =
            HighsCliqueTable::CliqueVar(baseBinCol, baseArcImplicationWhenOne);
        HighsCliqueTable::CliqueVar substCliqueVar =
            HighsCliqueTable::CliqueVar(binCol, arcImplicationWhenOne);
        clique[0] = stayCliqueVar;
        clique[1] = substCliqueVar.complement();
        mipsolver.mipdata_->cliquetable.addClique(mipsolver, clique.data(), 2);
        clique[0] = stayCliqueVar.complement();
        clique[1] = substCliqueVar;
        mipsolver.mipdata_->cliquetable.addClique(mipsolver, clique.data(), 2);
        if (globaldom.infeasible()) return false;
      }
    }
  }

  // Skip any clique smaller than size 3
  if (numRows <= 5) return false;

  // Greedily search neighbours of largest degree column for a double-sided
  // clique (corresponds to a single machine schedule)
  std::vector<HighsInt> potentialNeighbours;
  potentialNeighbours.reserve(largestDegree);
  for (const auto& entry : adjacency) {
    auto arc = entry.key();
    const HighsInt col = std::get<1>(arc);
    if (col == largestDegreeCol) {
      potentialNeighbours.emplace_back(std::get<0>(arc));
    }
  }
  pdqsort(potentialNeighbours.begin(), potentialNeighbours.end(),
          [&](const HighsInt c1, const HighsInt c2) {
            return degrees[c1] > degrees[c2];
          });

  std::vector<HighsInt> neighbours;
  neighbours.reserve(largestDegree + 1);
  neighbours.emplace_back(largestDegreeCol);
  double releaseDate = globaldom.col_lower_[largestDegreeCol];
  std::vector<double> processingTimes(largestDegree + 1, kHighsInf);
  // Iterate over potential neighbours and check validity
  // Greedily add neighbours if they're valid
  for (HighsInt col : potentialNeighbours) {
    bool valid_neighbour = true;
    for (HighsInt neighbour : neighbours) {
      const auto fromArc = adjacency.find({col, neighbour});
      const auto toArc = adjacency.find({neighbour, col});
      // Need to verify that the to and fromArc exist and are opposites
      if (toArc == nullptr || fromArc == nullptr ||
          std::get<1>(*fromArc) != std::get<1>(*toArc) ||
          std::get<2>(*fromArc) == std::get<2>(*toArc)) {
        valid_neighbour = false;
        break;
      }
    }
    if (!valid_neighbour) continue;
    double minProcessingTime = kHighsInf;
    // Extract the processing times from the arcs
    // Turn pair-wise dependent times p_ji into p_j = min{p_ji : i \in N / j}
    for (size_t i = 0; i != neighbours.size(); ++i) {
      HighsInt neighbour = neighbours[i];
      const auto fromArc = adjacency.find({col, neighbour});
      const auto toArc = adjacency.find({neighbour, col});
      minProcessingTime = std::min(minProcessingTime, std::get<0>(*fromArc));
      processingTimes[i] = std::min(processingTimes[i], std::get<0>(*toArc));
    }
    releaseDate = std::min(globaldom.col_lower_[col], releaseDate);
    processingTimes[neighbours.size()] = minProcessingTime;
    neighbours.emplace_back(col);
  }
  if (neighbours.size() < 3) return false;

  // Now populate the actual inequalities
  vals.resize(neighbours.size());
  inds.resize(neighbours.size());
  rhss.resize(neighbours.size());
  std::vector<HighsInt> binColCutPositions(mipsolver.numCol(), -1);
  for (size_t i = 0; i != neighbours.size(); ++i) {
    vals[i].reserve(neighbours.size());
    inds[i].reserve(neighbours.size());
    rhss[i] -= releaseDate;
    HighsInt col = neighbours[i];
    for (size_t j = 0; j != neighbours.size(); ++j) {
      if (i == j) continue;
      HighsInt neighbour = neighbours[j];
      const auto toArc = adjacency.find({neighbour, col});
      assert(toArc != nullptr);
      const HighsInt binCol = std::get<1>(*toArc);
      double val = processingTimes[j];
      if (std::get<2>(*toArc) == ArcType::kImplicationWhenZero) {
        rhss[i] -= val;
        val = -val;
      }
      HighsInt& pos = binColCutPositions[binCol];
      if (pos == -1) {
        pos = static_cast<HighsInt>(inds[i].size());
        inds[i].emplace_back(binCol);
        vals[i].emplace_back(val);
      } else {
        vals[i][pos] += val;
      }
    }
    // Drop any potential zeros (if the same binary sets the order of many jobs)
    size_t j = 0;
    while (j != inds[i].size()) {
      binColCutPositions[inds[i][j]] = -1;
      if (vals[i][j] == 0.0) {
        std::swap(inds[i][j], inds[i].back());
        std::swap(vals[i][j], vals[i].back());
        inds[i].pop_back();
        vals[i].pop_back();
      } else {
        ++j;
      }
    }
    // Put the job start time on the LHS
    inds[i].emplace_back(col);
    vals[i].emplace_back(-1);
  }

  return true;
}

void HighsMachineSchedSeparator::separateLpSolution(
    HighsLpRelaxation& lpRelaxation, HighsLpAggregator& lpAggregator,
    HighsTransformedLp& transLp, HighsCutPool& cutpool) {
  // Only try to separate once
  if (already_tried) return;
  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  std::vector<std::vector<double>> vals;
  std::vector<std::vector<HighsInt>> inds;
  std::vector<double> rhss;
  has_single_machine_schedule = findSingleMachineScheduleClique(
      vals, inds, rhss, transLp.getGlobaldom(), mip);
  if (!has_single_machine_schedule) {
    already_tried = true;
    return;
  }

  // Load the cuts
  const std::vector<double>& lpSolution = lpRelaxation.getSolution().col_value;
  for (size_t i = 0; i != inds.size(); ++i) {
    double viol = -rhss[i];
    for (size_t j = 0; j != inds[i].size(); j++) {
      viol += lpSolution[inds[i][j]] * vals[i][j];
    }
    if (viol >= 10 * mip.mipdata_->feastol)
      cutpool.addCut(mip, inds[i].data(), vals[i].data(),
                     static_cast<HighsInt>(inds[i].size()), rhss[i]);
  }
  already_tried = true;
}
