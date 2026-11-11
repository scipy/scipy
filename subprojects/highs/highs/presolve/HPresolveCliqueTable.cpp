/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolveCliqueTable.h"

#include <algorithm>
#include <cassert>

#include "../extern/pdqsort/pdqsort.h"

using CliqueVar = HighsCliqueTable::CliqueVar;
using Clique = HighsCliqueTable::Clique;

void HPresolveCliqueTable::rebuild(HighsCliqueTable& table) {
  setCliqueTable(&table);
  colStates.assign(table.colDeleted.size(), ColState{});
}

void HPresolveCliqueTable::checkCompactClique(const HighsInt cliqueId,
                                              const HighsInt threshold,
                                              const HighsInt activeSize,
                                              const HighsInt actualSize,
                                              const bool equality,
                                              const HighsInt origin) {
  const Clique& clique = table->cliques[cliqueId];
  if (activeSize == 2 ||
      clique.numZeroFixed >= std::max(threshold, actualSize >> 1)) {
    shortenedClique.clear();
    shortenedClique.reserve(activeSize);
    for (HighsInt i = clique.start; i != clique.end; ++i) {
      if (!table->colDeleted[table->cliqueentries[i].col])
        shortenedClique.push_back(table->cliqueentries[i]);
    }
    table->removeClique(cliqueId, false);
    table->doAddClique(shortenedClique.data(),
                       static_cast<HighsInt>(shortenedClique.size()), equality,
                       origin);
  }
}

bool HPresolveCliqueTable::fixCol(HighsInt col, bool val,
                                  std::vector<CliqueVar>& impliedFixings) {
  const ColState state = colStates[col];
  if (!state.isActive()) {
    if (state.isEliminated()) return true;
    return state.isFixedTo(val);
  }

  if (table->invertedHashList[2 * col].empty() &&
      table->invertedHashList[2 * col + 1].empty() &&
      table->invertedHashListSizeTwo[2 * col].empty() &&
      table->invertedHashListSizeTwo[2 * col + 1].empty()) {
    colStates[col] = ColState::fixed(val);
    table->colDeleted[col] = true;
    return true;
  }

  fixingQueue.clear();
  fixingQueue.emplace_back(col, val);
  size_t nextFixing = 0;

  auto collectIncidentCliques = [&](const CliqueVar v) {
    incidentCliques.clear();
    table->invertedHashList[v.index()].for_each(
        [&](const HighsInt cliqueId, HighsInt) {
          incidentCliques.push_back(cliqueId);
        });
    table->invertedHashListSizeTwo[v.index()].for_each(
        [&](const HighsInt cliqueId) { incidentCliques.push_back(cliqueId); });
  };

  while (nextFixing != fixingQueue.size()) {
    CliqueVar v = fixingQueue[nextFixing++];
    const ColState state = colStates[v.col];
    if (state.isEliminated() || state.isFixedTo(v.val)) continue;
    if (state.isFixed()) return false;
    colStates[v.col] = ColState::fixed(v.val);
    table->colDeleted[v.col] = true;
    if (!(v.val == val && static_cast<HighsInt>(v.col) == col))
      impliedFixings.emplace_back(v);

    // Fix all other literals in incident cliques to be inactive
    collectIncidentCliques(v);
    for (const HighsInt cliqueId : incidentCliques) {
      if (table->cliques[cliqueId].start == -1) continue;
      for (HighsInt i = table->cliques[cliqueId].start;
           i != table->cliques[cliqueId].end; ++i) {
        CliqueVar v2 = table->cliqueentries[i];
        if (v.col == v2.col) continue;
        const ColState otherState = colStates[v2.col];
        if (otherState.isFixedTo(v2.val)) return false;
        if (otherState.isActive()) {
          fixingQueue.emplace_back(static_cast<HighsUInt>(v2.col), 1 - v2.val);
        }
      }
      table->removeClique(cliqueId);
    }

    // Remove complement-literal from cliques
    collectIncidentCliques(v.complement());
    table->invertedHashList[v.complement().index()].clear();
    table->invertedHashListSizeTwo[v.complement().index()].clear();
    for (const HighsInt cliqueId : incidentCliques) {
      Clique& clique = table->cliques[cliqueId];
      if (clique.start == -1) continue;
      const bool equality = clique.equality;
      const HighsInt origin = clique.origin;
      ++clique.numZeroFixed;
      const HighsInt actualSize = clique.end - clique.start;
      const HighsInt activeSize = clique.numActive();
      if (activeSize <= 1) {
        if (equality) {
          if (activeSize == 0) return false;
          for (HighsInt i = clique.start; i != clique.end; ++i) {
            if (!table->colDeleted[table->cliqueentries[i].col]) {
              fixingQueue.push_back(table->cliqueentries[i]);
              break;
            }
          }
        }
        table->removeClique(cliqueId);
        continue;
      }

      if (!table->inPresolveProbing)
        checkCompactClique(cliqueId, 10, activeSize, actualSize, equality,
                           origin);
    }
  }
  return true;
}

void HPresolveCliqueTable::eliminateCol(const HighsInt col) {
  if (colStates[col].isEliminated()) return;
  colStates[col] = ColState(ColState::kEliminated);
  if (table->colDeleted[col]) return;
  table->colDeleted[col] = true;

  incidentCliques.clear();
  table->invertedHashList[2 * col].for_each(
      [&](const HighsInt cliqueId, HighsInt) {
        incidentCliques.push_back(cliqueId);
      });
  table->invertedHashList[2 * col].clear();
  table->invertedHashListSizeTwo[2 * col].for_each(
      [&](const HighsInt cliqueId) { incidentCliques.push_back(cliqueId); });
  table->invertedHashListSizeTwo[2 * col].clear();
  table->invertedHashList[2 * col + 1].for_each(
      [&](const HighsInt cliqueId, HighsInt) {
        incidentCliques.push_back(cliqueId);
      });
  table->invertedHashList[2 * col + 1].clear();
  table->invertedHashListSizeTwo[2 * col + 1].for_each(
      [&](const HighsInt cliqueId) { incidentCliques.push_back(cliqueId); });
  table->invertedHashListSizeTwo[2 * col + 1].clear();

  pdqsort(incidentCliques.begin(), incidentCliques.end());

  for (const HighsInt cliqueId : incidentCliques) {
    Clique& clique = table->cliques[cliqueId];
    if (clique.start == -1) continue;
    ++clique.numZeroFixed;
    clique.origin = -1;
    clique.equality = false;
    const HighsInt actualSize = clique.end - clique.start;
    const HighsInt activeSize = clique.numActive();
    if (activeSize <= 1) {
      table->removeClique(cliqueId, false);
      continue;
    }
    if (!table->inPresolveProbing)
      checkCompactClique(cliqueId, 10, activeSize, actualSize, false, -1);
  }
}

bool HPresolveCliqueTable::substituteCol(
    const HighsInt substCol, const CliqueVar replacement,
    std::vector<CliqueVar>& impliedFixings) {
  if (substCol == static_cast<HighsInt>(replacement.col)) {
    // Substituting itself is either redundant or infeasible
    return replacement.val == 1;
  }
  if (!colStates[substCol].isActive()) {
    // Don't do anything if col to be substituted is already handled
    return true;
  }
  if (colStates[replacement.col].isEliminated()) {
    // Eliminate substituted column if replacement is eliminated
    eliminateCol(substCol);
    return true;
  }
  if (colStates[replacement.col].isFixed()) {
    // If replacement is fixed then bypass direct substitutions via fixings
    const bool fixVal = colStates[replacement.col].isFixedTo(replacement.val);
    impliedFixings.push_back(CliqueVar(substCol, fixVal));
    if (!fixCol(substCol, fixVal, impliedFixings)) return false;
    return true;
  }
  struct Overlap {
    HighsInt cliqueId;
    HighsInt substPos;
    HighsInt replacePos;
  };
  std::vector<Overlap> overlaps;
  auto collectOverlaps = [&](const CliqueVar v) {
    table->invertedHashList[v.index()].for_each(
        [&](const HighsInt cliqueId, const HighsInt substPos) {
          const HighsInt* replacePos =
              table->invertedHashList[replacement.index()].find(cliqueId);
          if (replacePos == nullptr)
            replacePos =
                table->invertedHashList[replacement.complement().index()].find(
                    cliqueId);
          if (replacePos != nullptr)
            overlaps.push_back({cliqueId, substPos, *replacePos});
        });
    table->invertedHashListSizeTwo[v.index()].for_each(
        [&](const HighsInt cliqueId) {
          HighsInt substPos = table->cliques[cliqueId].start;
          HighsInt replacePos = substPos + 1;
          if (table->cliqueentries[replacePos] == v)
            std::swap(substPos, replacePos);
          if (table->cliqueentries[replacePos].col == replacement.col)
            overlaps.push_back({cliqueId, substPos, replacePos});
        });
  };
  collectOverlaps(CliqueVar(substCol, 0));
  collectOverlaps(CliqueVar(substCol, 1));
  pdqsort(overlaps.begin(), overlaps.end(),
          [](const Overlap& a, const Overlap& b) {
            return a.cliqueId < b.cliqueId;
          });

  std::vector<CliqueVar> potentialFixings;

  for (const Overlap& overlap : overlaps) {
    const HighsInt cliqueId = overlap.cliqueId;
    if (table->cliques[cliqueId].start == -1) continue;
    const HighsInt substPos = overlap.substPos;
    const HighsInt replacePos = overlap.replacePos;
    table->cliques[cliqueId].origin = -1;
    const bool equality = table->cliques[cliqueId].equality;

    const CliqueVar substVar = table->cliqueentries[substPos].val
                                   ? replacement
                                   : replacement.complement();
    const CliqueVar replacementVar = table->cliqueentries[replacePos];

    // If both replacement and its complement will exist in new clique
    // then clique is automatically fulfilled. Set all other literals to 0
    if (substVar.val != replacementVar.val) {
      for (HighsInt i = table->cliques[cliqueId].start;
           i != table->cliques[cliqueId].end; ++i) {
        if (i != substPos && i != replacePos &&
            !table->colDeleted[table->cliqueentries[i].col]) {
          potentialFixings.push_back(table->cliqueentries[i].complement());
        }
      }
      table->removeClique(cliqueId, false);
      continue;
    }

    // If replacement will exist twice in new clique then the literal
    // has to take value 0
    potentialFixings.push_back(substVar.complement());
    shortenedClique.clear();
    shortenedClique.reserve(table->cliques[cliqueId].end -
                            table->cliques[cliqueId].start);
    for (HighsInt i = table->cliques[cliqueId].start;
         i != table->cliques[cliqueId].end; ++i) {
      if (i != substPos && i != replacePos &&
          !table->colDeleted[table->cliqueentries[i].col]) {
        shortenedClique.push_back(table->cliqueentries[i]);
      }
    }
    table->removeClique(cliqueId, false);

    if (shortenedClique.size() >= 2) {
      table->doAddClique(shortenedClique.data(),
                         static_cast<HighsInt>(shortenedClique.size()),
                         equality, -1);
    } else if (equality) {
      if (shortenedClique.empty()) return false;
      potentialFixings.push_back(shortenedClique[0]);
    }
  }

  // Now substitute all entries
  table->replaceLiteral(CliqueVar(substCol, 1), replacement);
  table->replaceLiteral(CliqueVar(substCol, 0), replacement.complement());
  colStates[substCol] = ColState{ColState::kEliminated};
  table->colDeleted[substCol] = true;

  for (CliqueVar v : potentialFixings) {
    const ColState state = colStates[v.col];
    if (state.isEliminated() || state.isFixedTo(v.val)) continue;
    if (state.isFixed()) return false;
    impliedFixings.push_back(v);
    if (!fixCol(static_cast<HighsInt>(v.col), v.val, impliedFixings)) {
      return false;
    }
  }

  return true;
}
