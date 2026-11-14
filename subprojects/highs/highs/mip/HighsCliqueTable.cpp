/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsCliqueTable.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <numeric>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "parallel/HighsCombinable.h"
#include "parallel/HighsParallel.h"
#include "presolve/HighsPostsolveStack.h"
#include "util/HighsSplay.h"

#define ADD_ZERO_WEIGHT_VARS

static std::pair<HighsCliqueTable::CliqueVar, HighsCliqueTable::CliqueVar>
sortedEdge(HighsCliqueTable::CliqueVar v1, HighsCliqueTable::CliqueVar v2) {
  if (v1.col > v2.col) return std::make_pair(v2, v1);

  return std::make_pair(v1, v2);
}

void HighsCliqueTable::addSizeTwoClique(HighsInt cliqueid) {
  HighsInt start = cliques[cliqueid].start;
  if (cliques[cliqueid].end - start != 2) return;
  auto edge = sortedEdge(cliqueentries[start], cliqueentries[start + 1]);
  HighsInt* sizeTwoCliqueId = sizeTwoCliques.find(edge);
  if (sizeTwoCliqueId == nullptr) {
    // first clique on this edge
    sizeTwoCliques.insert(edge, cliqueid);
    sizeTwoLinks[cliqueid] = SizeTwoLinks{-1, -1};
    return;
  }
  // insert the clique after the stored one
  HighsInt prev = *sizeTwoCliqueId;
  HighsInt next = sizeTwoLinks[prev].next;
  sizeTwoLinks[cliqueid] = SizeTwoLinks{prev, next};
  sizeTwoLinks[prev].next = cliqueid;
  if (next != -1) sizeTwoLinks[next].prev = cliqueid;
}

void HighsCliqueTable::removeSizeTwoClique(HighsInt cliqueid) {
  HighsInt start = cliques[cliqueid].start;
  if (cliques[cliqueid].end - start != 2) return;
  HighsInt prev = sizeTwoLinks[cliqueid].prev;
  HighsInt next = sizeTwoLinks[cliqueid].next;
  if (prev != -1) {
    sizeTwoLinks[prev].next = next;
  } else {
    // the clique is the stored one, so store the next one if there is one
    auto edge = sortedEdge(cliqueentries[start], cliqueentries[start + 1]);
    HighsInt* sizeTwoCliqueId = sizeTwoCliques.find(edge);
    assert(sizeTwoCliqueId != nullptr && *sizeTwoCliqueId == cliqueid);
    if (next != -1)
      *sizeTwoCliqueId = next;
    else
      sizeTwoCliques.erase(edge);
  }
  if (next != -1) sizeTwoLinks[next].prev = prev;
  sizeTwoLinks[cliqueid] = SizeTwoLinks{-1, -1};
}

void HighsCliqueTable::unlink(HighsInt pos, HighsInt cliqueid) {
  assert(pos >= 0);
  --numcliquesvar[cliqueentries[pos].index()];

  HighsInt cliquelen = cliques[cliqueid].end - cliques[cliqueid].start;
  if (cliquelen == 2)
    invertedHashListSizeTwo[cliqueentries[pos].index()].erase(cliqueid);
  else
    invertedHashList[cliqueentries[pos].index()].erase(cliqueid);
}

void HighsCliqueTable::link(HighsInt pos, HighsInt cliqueid) {
  assert(pos >= 0);
  assert(!colDeleted[cliqueentries[pos].col]);
  ++numcliquesvar[cliqueentries[pos].index()];

  HighsInt cliquelen = cliques[cliqueid].end - cliques[cliqueid].start;
  if (cliquelen == 2)
    invertedHashListSizeTwo[cliqueentries[pos].index()].insert(cliqueid);
  else
    invertedHashList[cliqueentries[pos].index()].insert(cliqueid, pos);
}

HighsInt HighsCliqueTable::findCommonCliqueId(int64_t& numQueries, CliqueVar v1,
                                              CliqueVar v2) const {
  ++numQueries;
  if (!invertedHashListSizeTwo[v1.index()].empty() &&
      !invertedHashListSizeTwo[v2.index()].empty()) {
    const HighsInt* sizeTwoCliqueId = sizeTwoCliques.find(sortedEdge(v1, v2));
    if (sizeTwoCliqueId != nullptr) return *sizeTwoCliqueId;
  }

  const HighsHashTree<HighsInt, HighsInt>& h1 = invertedHashList[v1.index()];
  const HighsHashTree<HighsInt, HighsInt>& h2 = invertedHashList[v2.index()];

  const HighsHashTableEntry<HighsInt, HighsInt>* commonClique =
      h1.find_common(h2);

  if (commonClique) return commonClique->key();
  return -1;
}

void HighsCliqueTable::resolveSubstitution(CliqueVar& v) const {
  while (colsubstituted[v.col]) {
    Substitution subst = substitutions[colsubstituted[v.col] - 1];
    v = v.val == 1 ? subst.replace : subst.replace.complement();
  }
}

void HighsCliqueTable::resolveSubstitution(HighsInt& col, double& val,
                                           double& offset) const {
  while (colsubstituted[col]) {
    Substitution subst = substitutions[colsubstituted[col] - 1];
    if (subst.replace.val == 0) {
      offset += val;
      val = -val;
    }
    col = subst.replace.col;
  }
}

HighsInt HighsCliqueTable::runCliqueSubsumption(
    const HighsDomain& globaldom, std::vector<CliqueVar>& clique) {
  // remove deleted cols
  clique.erase(
      std::remove_if(clique.begin(), clique.end(),
                     [&](CliqueVar clqvar) { return colDeleted[clqvar.col]; }),
      clique.end());

  // return if clique does not have enough elements
  if (clique.size() <= 2) return 0;

  // callback for removing cliques
  HighsInt nremoved = 0;
  auto removeCliques = [&](HighsInt cliqueid) {
    ++nremoved;
    cliques[cliqueid].origin = kHighsIInf;
    removeClique(cliqueid);
  };

  // check clique
  bool redundant;
  HighsInt dominatingOrigin;
  cliqueSubsumption(clique, redundant, dominatingOrigin, removeCliques);

  if (redundant) clique.clear();

  if (!infeasvertexstack.empty()) {
    clique.erase(
        std::remove_if(clique.begin(), clique.end(),
                       [&](CliqueVar v) { return globaldom.isFixed(v.col); }),
        clique.end());
  }

  return nremoved;
}

void HighsCliqueTable::cliqueSubsumption(
    const std::vector<CliqueVar>& clique, bool& redundant,
    HighsInt& dominatingOrigin,
    const std::function<void(HighsInt)>& removeCliqueCallback) {
  // collect indices of cliques that contain variables from the
  // provided vector
  collectCliques(clique);

  // initialise
  redundant = false;
  dominatingOrigin = kHighsIInf;

  for (HighsInt cliqueid : cliquehitinds) {
    // remember hits and zero out vector
    HighsInt hits = cliquehits[cliqueid];
    cliquehits[cliqueid] = 0;

    if (hits == static_cast<HighsInt>(clique.size())) {
      // clique is redundant
      redundant = true;
      if (cliques[cliqueid].origin != kHighsIInf &&
          cliques[cliqueid].origin != -1)
        dominatingOrigin = cliques[cliqueid].origin;
    } else if (cliques[cliqueid].numActive() == hits) {
      // clique is subset of another clique
      if (cliques[cliqueid].equality) {
        // subset of clique is an equality clique (stored in the clique table),
        // thus the remaining variables can be fixed
        bool sizeTwo = cliques[cliqueid].end - cliques[cliqueid].start == 2;
        for (CliqueVar v : clique) {
          bool vHasClq =
              sizeTwo ? invertedHashListSizeTwo[v.index()].contains(cliqueid)
                      : invertedHashList[v.index()].contains(cliqueid);
          if (!vHasClq) infeasvertexstack.push_back(v);
        }
      } else {
        // clique from clique table contains a subset of the variables, thus it
        // can be removed
        removeCliqueCallback(cliqueid);
      }
    }
  }
  // clear vector of cliques indices
  cliquehitinds.clear();
}

void HighsCliqueTable::collectCliques(const std::vector<CliqueVar>& clique) {
  // resize vector
  if (cliquehits.size() < cliques.size()) cliquehits.resize(cliques.size());

  // collect indices of cliques that contain variables from the
  // provided vector
  for (CliqueVar v : clique) {
    invertedHashList[v.index()].for_each([&](HighsInt cliqueid) {
      if (cliquehits[cliqueid] == 0) cliquehitinds.push_back(cliqueid);
      ++cliquehits[cliqueid];
    });

    invertedHashListSizeTwo[v.index()].for_each([&](HighsInt cliqueid) {
      if (cliquehits[cliqueid] == 0) cliquehitinds.push_back(cliqueid);
      ++cliquehits[cliqueid];
    });
  }
}

void HighsCliqueTable::bronKerboschRecurse(BronKerboschData& data,
                                           HighsInt Plen, const CliqueVar* X,
                                           HighsInt Xlen) const {
  double w = data.wR;

  for (HighsInt i = 0; i != Plen; ++i) w += data.P[i].weight(data.sol);

  if (w < data.minW - data.feastol) return;

  if (Plen == 0 && Xlen == 0) {
    std::vector<CliqueVar> clique = data.R;

    if (data.minW < w - data.feastol) {
      data.maxcliques -= data.cliques.size();
      data.cliques.clear();
      data.minW = w;
    }
    data.cliques.emplace_back(std::move(clique));
    // do not further search for cliques that are violated less than this
    // current clique
    return;
  }

  ++data.ncalls;

  if (data.stop()) return;

  double pivweight = -1.0;
  CliqueVar pivot{0, 0};

  for (HighsInt i = 0; i != Xlen; ++i) {
    if (X[i].weight(data.sol) > pivweight) {
      pivweight = X[i].weight(data.sol);
      pivot = X[i];
      if (pivweight >= 1.0 - data.feastol) break;
    }
  }

  if (pivweight < 1.0 - data.feastol) {
    for (HighsInt i = 0; i != Plen; ++i) {
      if (data.P[i].weight(data.sol) > pivweight) {
        pivweight = data.P[i].weight(data.sol);
        pivot = data.P[i];
        if (pivweight >= 1.0 - data.feastol) break;
      }
    }
  }

  std::vector<CliqueVar> PminusNu;
  PminusNu.reserve(Plen);
  queryNeighbourhood(data.neighbourhoodInds, data.numNeighbourhoodQueries,
                     pivot, data.P.data(), Plen);
  data.neighbourhoodInds.push_back(Plen);
  HighsInt k = 0;
  for (HighsInt i : data.neighbourhoodInds) {
    while (k < i) PminusNu.push_back(data.P[k++]);
    ++k;
  }

  pdqsort(PminusNu.begin(), PminusNu.end(), [&](CliqueVar a, CliqueVar b) {
    return std::make_pair(a.weight(data.sol), a.index()) >
           std::make_pair(b.weight(data.sol), b.index());
  });

  std::vector<CliqueVar> localX;
  localX.insert(localX.end(), X, X + Xlen);

  for (CliqueVar v : PminusNu) {
    HighsInt newPlen = partitionNeighbourhood(data.neighbourhoodInds,
                                              data.numNeighbourhoodQueries, v,
                                              data.P.data(), Plen);
    HighsInt newXlen = partitionNeighbourhood(data.neighbourhoodInds,
                                              data.numNeighbourhoodQueries, v,
                                              localX.data(), localX.size());

    // add v to R, update the weight, and do the recursive call
    data.R.push_back(v);
    double wv = v.weight(data.sol);
    data.wR += wv;
    bronKerboschRecurse(data, newPlen, localX.data(), newXlen);
    if (data.stop()) return;

    // remove v from R restore the weight and continue the loop in this call
    data.R.pop_back();
    data.wR -= wv;

    w -= wv;
    if (w < data.minW) return;
    // find the position of v in the vertices removed from P for the recursive
    // call
    // and also remove it from the set P for this call
    HighsInt vpos = -1;
    for (HighsInt i = newPlen; i != Plen; ++i) {
      if (data.P[i] == v) {
        vpos = i;
        break;
      }
    }

    // do the removal by first swapping it to the end of P and reduce the size
    // of P accordingly
    assert(vpos != -1);

    --Plen;
    std::swap(data.P[vpos], data.P[Plen]);

    localX.push_back(v);
  }
}

#if 0
static void printRow(const HighsDomain& domain, const HighsInt* inds,
                     const double* vals, HighsInt len, double lhs, double rhs) {
  printf("%g <= ", lhs);

  for (HighsInt i = 0; i != len; ++i) {
    char sign = vals[i] > 0 ? '+' : '-';
    char var = domain.isBinary(inds[i]) ? 'x' : 'y';
    printf("%c%g %c%" HIGHSINT_FORMAT " ", sign, std::abs(vals[i]), var, inds[i]);
  }

  printf("<= %g\n", rhs);
}

static void printClique(
    const std::vector<HighsCliqueTable::CliqueVar>& clique) {
  bool first = true;
  for (HighsCliqueTable::CliqueVar v : clique) {
    if (!first) printf("+ ");
    char complemented = v.val == 0 ? '~' : ' ';
    printf("%cx%" HIGHSINT_FORMAT " ", complemented, v.col);
    first = false;
  }

  printf("<= 1\n");
}
#endif

void HighsCliqueTable::doAddClique(const CliqueVar* cliquevars,
                                   HighsInt numcliquevars, bool equality,
                                   HighsInt origin) {
  HighsInt cliqueid;

  if (freeslots.empty()) {
    cliqueid = cliques.size();
    cliques.emplace_back();
    sizeTwoLinks.push_back(SizeTwoLinks{-1, -1});
  } else {
    cliqueid = freeslots.back();
    freeslots.pop_back();
  }

  cliques[cliqueid].equality = equality;
  cliques[cliqueid].origin = origin;

  decltype(freespaces)::iterator it;
  HighsInt maxEnd;
  if (freespaces.empty() ||
      (it = freespaces.lower_bound(
           std::make_pair(numcliquevars, HighsInt{-1}))) == freespaces.end()) {
    cliques[cliqueid].start = cliqueentries.size();
    cliques[cliqueid].end = cliques[cliqueid].start + numcliquevars;
    maxEnd = cliques[cliqueid].end;
    cliqueentries.resize(cliques[cliqueid].end);
  } else {
    auto freespace = *it;
    freespaces.erase(it);

    cliques[cliqueid].start = freespace.second;
    cliques[cliqueid].end = cliques[cliqueid].start + numcliquevars;
    maxEnd = cliques[cliqueid].start + freespace.first;
  }

  cliques[cliqueid].numZeroFixed = 0;

  bool fixtozero = false;
  HighsInt k = cliques[cliqueid].start;
  for (HighsInt i = 0; i != numcliquevars; ++i) {
    CliqueVar v = cliquevars[i];

    resolveSubstitution(v);

    if (fixtozero) {
      infeasvertexstack.push_back(v);
      continue;
    }

    // due to substitutions the variable may occur together with its complement
    // in this clique and we can fix all other variables in the clique to zero:
    //          x + ~x + ... <= 1
    //   <=> x + 1 - x + ... <= 1
    //   <=>             ... <= 0
    bool clqHasVCompl =
        numcliquevars == 2
            ? invertedHashListSizeTwo[v.complement().index()].contains(cliqueid)
            : invertedHashList[v.complement().index()].contains(cliqueid);

    if (clqHasVCompl) {
      fixtozero = true;
      for (HighsInt j = cliques[cliqueid].start; j != k; ++j) {
        if (cliqueentries[j].col != v.col)
          infeasvertexstack.push_back(cliqueentries[j]);
        unlink(j, cliqueid);
      }
      k = cliques[cliqueid].start;
      continue;
    }

    // due to substitutions the variable may occur twice in this clique and
    // we can fix it to zero:  x + x + ... <= 1  <=>  2x <= 1 <=> x <= 0.5 <=>
    // x = 0
    bool inserted;
    if (numcliquevars == 2)
      inserted = invertedHashListSizeTwo[v.index()].insert(cliqueid);
    else
      inserted = invertedHashList[v.index()].insert(cliqueid, k);

    if (!inserted) {
      infeasvertexstack.push_back(v);
      continue;
    }

    cliqueentries[k] = v;
    ++numcliquesvar[v.index()];
    ++k;
  }

  if (maxEnd > k) {
    if (int(cliqueentries.size()) == maxEnd)
      cliqueentries.resize(k);
    else
      freespaces.emplace(maxEnd - k, k);

    if (cliques[cliqueid].end > k) {
      switch (k - cliques[cliqueid].start) {
        case 0:
          // clique empty, so just mark it as deleted
          cliques[cliqueid].start = -1;
          cliques[cliqueid].end = -1;
          freeslots.push_back(cliqueid);
          return;
        case 1:
          // size 1 clique is redundant, so unlink the single linked entry
          // and mark it as deleted
          unlink(cliques[cliqueid].start, cliqueid);
          cliques[cliqueid].start = -1;
          cliques[cliqueid].end = -1;
          freeslots.push_back(cliqueid);
          return;
        case 2:
          // due to substitutions the clique became smaller and is now of size
          // two as a result we need to link it to the size two cliqueset
          // instead of the normal cliqueset
          assert(cliqueid >= 0 &&
                 cliqueid < static_cast<HighsInt>(cliques.size()));
          assert(cliques[cliqueid].start >= 0 &&
                 cliques[cliqueid].start <
                     static_cast<HighsInt>(cliqueentries.size()));
          unlink(cliques[cliqueid].start, cliqueid);
          unlink(cliques[cliqueid].start + 1, cliqueid);

          cliques[cliqueid].end = k;

          link(cliques[cliqueid].start, cliqueid);
          link(cliques[cliqueid].start + 1, cliqueid);
          break;
        default:
          cliques[cliqueid].end = k;
      }
    }
  }

  HighsInt cliqueLen = cliques[cliqueid].end - cliques[cliqueid].start;
  numEntries += cliqueLen;
  addSizeTwoClique(cliqueid);
}
struct ThreadNeighbourhoodQueryData {
  int64_t numQueries;
  std::vector<HighsInt> neighbourhoodInds;
};

void HighsCliqueTable::queryNeighbourhood(
    std::vector<HighsInt>& neighbourhoodInds, int64_t& numQueries, CliqueVar v,
    CliqueVar* q, HighsInt N) const {
  neighbourhoodInds.clear();

  if (numCliques(v) == 0) return;

  if (!allowParallel ||
      numEntries - sizeTwoCliques.size() * 2 < minEntriesForParallelism) {
    for (HighsInt i = 0; i < N; ++i) {
      if (haveCommonClique(numQueries, v, q[i])) neighbourhoodInds.push_back(i);
    }
  } else {
    auto neighbourhoodData =
        makeHighsCombinable<ThreadNeighbourhoodQueryData>([N]() {
          ThreadNeighbourhoodQueryData d;
          d.neighbourhoodInds.reserve(N);
          d.numQueries = 0;
          return d;
        });
    highs::parallel::for_each(
        0, N,
        [this, &neighbourhoodData, v, q](HighsInt start, HighsInt end) {
          ThreadNeighbourhoodQueryData& d = neighbourhoodData.local();
          for (HighsInt i = start; i < end; ++i) {
            if (haveCommonClique(d.numQueries, v, q[i]))
              d.neighbourhoodInds.push_back(i);
          }
        },
        10);

    neighbourhoodData.combine_each([&](ThreadNeighbourhoodQueryData& d) {
      neighbourhoodInds.insert(neighbourhoodInds.end(),
                               d.neighbourhoodInds.begin(),
                               d.neighbourhoodInds.end());
      numQueries += d.numQueries;
    });
    pdqsort(neighbourhoodInds.begin(), neighbourhoodInds.end());
  }
}

HighsInt HighsCliqueTable::partitionNeighbourhood(
    std::vector<HighsInt>& neighbourhoodInds, int64_t& numQueries, CliqueVar v,
    CliqueVar* q, HighsInt N) const {
  queryNeighbourhood(neighbourhoodInds, numQueries, v, q, N);

  for (size_t i = 0; i < neighbourhoodInds.size(); ++i)
    std::swap(q[i], q[neighbourhoodInds[i]]);

  return static_cast<HighsInt>(neighbourhoodInds.size());
}

HighsInt HighsCliqueTable::shrinkToNeighbourhood(
    std::vector<HighsInt>& neighbourhoodInds, int64_t& numQueries, CliqueVar v,
    CliqueVar* q, HighsInt N) {
  queryNeighbourhood(neighbourhoodInds, numQueries, v, q, N);

  for (size_t i = 0; i < neighbourhoodInds.size(); ++i)
    q[i] = q[neighbourhoodInds[i]];

  return static_cast<HighsInt>(neighbourhoodInds.size());
}

bool HighsCliqueTable::fixCol(HighsDomain& globaldom, CliqueVar v,
                              bool doProcessInfeasibleVertices) {
  bool wasfixed = globaldom.isFixed(v.col);
  globaldom.fixCol(v.col, static_cast<double>(1 - v.val));
  if (globaldom.infeasible()) return false;
  if (!wasfixed) {
    ++nfixings;
    infeasvertexstack.push_back(v);
    if (doProcessInfeasibleVertices) processInfeasibleVertices(globaldom);
  }
  return true;
}

bool HighsCliqueTable::processNewEdge(HighsDomain& globaldom, CliqueVar v1,
                                      CliqueVar v2) {
  if (v1.col == v2.col) {
    if (v1.val == v2.val) {
      fixCol(globaldom, v1, true);
      return false;
    }

    return true;
  }

  // invertedEdgeCache.erase(sortedEdge(v1, v2));

  if (haveCommonClique(v1.complement(), v2)) {
    fixCol(globaldom, v2, true);
    return false;
  } else if (haveCommonClique(v2.complement(), v1)) {
    fixCol(globaldom, v1, true);
    return false;
  } else {
    // return if complements are not in a clique
    if (!foundCover(globaldom, v1.complement(), v2.complement())) return false;
    if (globaldom.infeasible()) return true;

    foundCover(globaldom, v1, v2);
    if (globaldom.isFixed(v1.col) || globaldom.isFixed(v2.col) ||
        globaldom.infeasible())
      return true;

    Substitution substitution;
    if (v2.col < v1.col) {
      if (v1.val == 1) v2 = v2.complement();

      substitution.substcol = v1.col;
      substitution.replace = v2;
    } else {
      if (v2.val == 1) v1 = v1.complement();

      substitution.substcol = v2.col;
      substitution.replace = v1;
    }

    recordSubstitution(substitution);

    return true;
  }
}

void HighsCliqueTable::replaceLiteral(const CliqueVar substitutedVar,
                                      const CliqueVar replacementVar) {
  HighsHashTree<HighsInt, HighsInt>& substList =
      invertedHashList[substitutedVar.index()];
  HighsHashTree<HighsInt, HighsInt>& replaceList =
      invertedHashList[replacementVar.index()];
  numcliquesvar[replacementVar.index()] +=
      numcliquesvar[substitutedVar.index()];
  numcliquesvar[substitutedVar.index()] = 0;

  substList.for_each([&](HighsInt cliqueid, HighsInt location) {
    replaceList.insert(cliqueid, location);
    cliqueentries[location] = replacementVar;
  });

  substList.clear();

  HighsHashTree<HighsInt>& substListSizeTwo =
      invertedHashListSizeTwo[substitutedVar.index()];
  HighsHashTree<HighsInt>& replaceListSizeTwo =
      invertedHashListSizeTwo[replacementVar.index()];

  substListSizeTwo.for_each([&](HighsInt cliqueid) {
    HighsInt pos = cliques[cliqueid].start;
    HighsInt otherPos = pos + 1;

    if (cliqueentries[otherPos] == substitutedVar) std::swap(pos, otherPos);

    replaceListSizeTwo.insert(cliqueid);
    removeSizeTwoClique(cliqueid);
    cliqueentries[pos] = replacementVar;
    addSizeTwoClique(cliqueid);
  });

  substListSizeTwo.clear();
}

void HighsCliqueTable::recordSubstitution(Substitution substitution) {
  substitutions.push_back(substitution);
  colsubstituted[substitution.substcol] =
      static_cast<HighsInt>(substitutions.size());

  replaceLiteral(CliqueVar(substitution.substcol, 1), substitution.replace);
  replaceLiteral(CliqueVar(substitution.substcol, 0),
                 substitution.replace.complement());
}

void HighsCliqueTable::addClique(const HighsMipSolver& mipsolver,
                                 CliqueVar* cliquevars, HighsInt numcliquevars,
                                 bool equality, HighsInt origin) {
  HighsDomain& globaldom = mipsolver.mipdata_->getDomain();
  mipsolver.mipdata_->debugSolution.checkClique(cliquevars, numcliquevars);
  const HighsInt maxNumCliqueVars = 100;

  // lambda for complementing clique
  auto complementClique = [&]() {
    for (HighsInt i = 0; i != numcliquevars; ++i) {
      cliquevars[i] = cliquevars[i].complement();
    }
  };

  // lambda for analysing the clique to see if all variables can be fixed
  auto fixAllVarsInClique = [&](bool& hasNewEdge) {
    for (HighsInt i = 0; i != numcliquevars; ++i) {
      if (!globaldom.isFixedToVal(cliquevars[i].col, cliquevars[i].val))
        continue;
      // column is fixed to 1, every other entry can be fixed to zero
      for (HighsInt k = 0; k != numcliquevars; ++k) {
        if (k == i) continue;
        if (!fixCol(globaldom, cliquevars[k])) return false;
      }
      processInfeasibleVertices(globaldom);
      return true;
    }

    if (numcliquevars > maxNumCliqueVars) return false;

    // todo, sort new clique to allow log n lookup of membership in size by
    // binary search

    for (HighsInt i = 0; i < numcliquevars - 1; ++i) {
      if (globaldom.isFixed(cliquevars[i].col)) continue;
      if (numCliques(cliquevars[i]) == 0 &&
          numCliques(cliquevars[i].complement()) == 0) {
        hasNewEdge = true;
        continue;
      }

      for (HighsInt j = i + 1; j < numcliquevars; ++j) {
        if (globaldom.isFixed(cliquevars[j].col)) continue;

        if (haveCommonClique(cliquevars[i], cliquevars[j])) continue;
        // todo: Instead of haveCommonClique use findCommonClique. If the
        // common clique is smaller than this clique check if it is a subset
        // of this clique. If it is a subset remove the clique and iterate the
        // process until either a common clique that is not a subset of this
        // one is found, or no common clique exists anymore in which case we
        // proceed with the code below and set hasNewEdge to true

        hasNewEdge = true;

        bool iscover = processNewEdge(globaldom, cliquevars[i], cliquevars[j]);
        if (globaldom.infeasible()) return false;

        if (!mipsolver.mipdata_->nodequeue.empty()) {
          const auto& v1Nodes =
              cliquevars[i].val == 1
                  ? mipsolver.mipdata_->nodequeue.getUpNodes(cliquevars[i].col)
                  : mipsolver.mipdata_->nodequeue.getDownNodes(
                        cliquevars[i].col);
          const auto& v2Nodes =
              cliquevars[j].val == 1
                  ? mipsolver.mipdata_->nodequeue.getUpNodes(cliquevars[j].col)
                  : mipsolver.mipdata_->nodequeue.getDownNodes(
                        cliquevars[j].col);

          if (!v1Nodes.empty() && !v2Nodes.empty()) {
            // general integer variables can become binary during the search
            // and the cliques might be discovered. That means we need to take
            // care here, since the set of nodes branched upwards or downwards
            // are not necessarily containing domain changes setting the
            // variables to the corresponding clique value but could be
            // redundant bound changes setting the upper bound to u >= 1 or
            // the lower bound to l <= 0.

            // itV1 will point to the first node where v1 is fixed to val and
            // endV1 to the end of the range of such nodes. Same for
            // itV2/endV2 with v2.
            auto itV1 = v1Nodes.lower_bound(std::make_pair(
                static_cast<double>(cliquevars[i].val), kHighsIInf));
            auto endV1 = v1Nodes.upper_bound(std::make_pair(
                static_cast<double>(cliquevars[i].val), kHighsIInf));
            auto itV2 = v2Nodes.lower_bound(std::make_pair(
                static_cast<double>(cliquevars[j].val), kHighsIInf));
            auto endV2 = v2Nodes.upper_bound(std::make_pair(
                static_cast<double>(cliquevars[j].val), kHighsIInf));

            if (itV1 != endV1 && itV2 != endV2 &&
                (itV1->second <= std::prev(endV2)->second ||
                 itV2->second <= std::prev(endV1)->second)) {
              // node ranges overlap, check for nodes that can be pruned
              while (itV1 != endV1 && itV2 != endV2) {
                if (itV1->second < itV2->second) {
                  ++itV1;
                } else if (itV2->second < itV1->second) {
                  ++itV2;
                } else {
                  // if (!mipsolver.submip)
                  //   printf("node %d can be pruned\n", itV2->second);
                  HighsInt prunedNode = itV2->second;
                  ++itV1;
                  ++itV2;
                  mipsolver.mipdata_->pruned_treeweight +=
                      mipsolver.mipdata_->nodequeue.pruneNode(prunedNode);
                }
              }
            }
          }
        }

        if (iscover) {
          for (HighsInt k = 0; k != numcliquevars; ++k) {
            if (k == i || k == j) continue;
            if (!fixCol(globaldom, cliquevars[k])) return false;
          }
          processInfeasibleVertices(globaldom);
          return true;
        }
      }
    }
    return false;
  };

  // lambda for checking the (complemented) clique
  auto checkClique = [&](bool& hasNewEdge, bool complement) {
    if (complement) complementClique();
    bool done = fixAllVarsInClique(hasNewEdge);
    if (complement) complementClique();
    if (done) return true;
    return false;
  };

  // resolve substitutions
  for (HighsInt i = 0; i != numcliquevars; ++i) {
    resolveSubstitution(cliquevars[i]);
  }

  // initialise flag
  bool hasNewEdge = false;

  // check if all variables can be fixed or infeasibility is detected
  if (checkClique(hasNewEdge, false)) return;
  if (globaldom.infeasible()) return;

  if (numcliquevars == 2 && equality) {
    // try complemented clique
    if (checkClique(hasNewEdge, true)) return;
    if (globaldom.infeasible()) return;
  }

  if (!hasNewEdge && origin == kHighsIInf) return;

  CliqueVar* unfixedend =
      std::remove_if(cliquevars, cliquevars + numcliquevars,
                     [&](CliqueVar v) { return globaldom.isFixed(v.col); });
  numcliquevars = unfixedend - cliquevars;
  if (numcliquevars < 2) return;

  doAddClique(cliquevars, numcliquevars, equality, origin);
  processInfeasibleVertices(globaldom);
}

void HighsCliqueTable::removeClique(HighsInt cliqueid, bool recordDeletedRow) {
  if (recordDeletedRow && cliques[cliqueid].origin != kHighsIInf &&
      cliques[cliqueid].origin != -1)
    deletedrows.push_back(cliques[cliqueid].origin);

  HighsInt start = cliques[cliqueid].start;
  assert(start != -1);
  HighsInt end = cliques[cliqueid].end;
  HighsInt len = end - start;
  removeSizeTwoClique(cliqueid);

  for (HighsInt i = start; i != end; ++i) {
    unlink(i, cliqueid);
  }

  freeslots.push_back(cliqueid);
  freespaces.emplace(len, start);

  cliques[cliqueid].start = -1;
  cliques[cliqueid].end = -1;
  numEntries -= len;
}

void HighsCliqueTable::fixLastActiveAndRemove(HighsDomain& globaldom,
                                              HighsInt cliqueid) {
  if (cliques[cliqueid].equality && cliques[cliqueid].numActive() == 1)
    for (HighsInt i = cliques[cliqueid].start; i != cliques[cliqueid].end;
         ++i) {
      if (colDeleted[cliqueentries[i].col]) continue;
      if (!globaldom.isFixed(static_cast<HighsInt>(cliqueentries[i].col))) {
        fixCol(globaldom, cliqueentries[i].complement(), false);
      } else if (globaldom.col_lower_[cliqueentries[i].col] !=
                 cliqueentries[i].val) {
        globaldom.markInfeasible();
      }
      break;
    }
  removeClique(cliqueid);
}

void HighsCliqueTable::extractCliques(
    const HighsMipSolver& mipsolver, std::vector<HighsInt>& inds,
    std::vector<double>& vals, std::vector<int8_t>& complementation, double rhs,
    HighsInt nbin, std::vector<HighsInt>& perm, std::vector<CliqueVar>& clique,
    double feastol) {
  HighsImplications& implics = mipsolver.mipdata_->implications;
  HighsDomain& globaldom = mipsolver.mipdata_->getDomain();

  perm.resize(inds.size());
  std::iota(perm.begin(), perm.end(), 0);

  auto binaryend = std::partition(perm.begin(), perm.end(), [&](HighsInt pos) {
    return globaldom.isBinary(inds[pos]);
  });
  nbin = binaryend - perm.begin();
  HighsInt ntotal = static_cast<HighsInt>(perm.size());

  // if not all variables are binary, we extract variable upper and lower bounds
  // constraints on the non-binary variable for each binary variable in the
  // constraint
  if (nbin < ntotal) {
    for (HighsInt i = 0; i != nbin; ++i) {
      HighsInt bincol = inds[perm[i]];
      HighsCDouble impliedub = static_cast<HighsCDouble>(rhs) - vals[perm[i]];
      if (implics.tooManyVarBounds()) break;
      for (HighsInt j = nbin; j != ntotal; ++j) {
        HighsInt col = inds[perm[j]];
        if (globaldom.isFixed(col)) continue;

        HighsCDouble colub =
            static_cast<HighsCDouble>(globaldom.col_upper_[col]) -
            globaldom.col_lower_[col];
        HighsCDouble implcolub = impliedub / vals[perm[j]];
        if (mipsolver.isColIntegral(col))
          implcolub = std::floor(static_cast<double>(implcolub) +
                                 mipsolver.mipdata_->feastol);

        if (implcolub < colub - feastol) {
          HighsCDouble coef;
          HighsCDouble constant;

          if (complementation[perm[i]] == -1) {
            coef = colub - implcolub;
            constant = implcolub;
          } else {
            coef = implcolub - colub;
            constant = colub;
          }

          if (complementation[perm[j]] == -1) {
            constant -= globaldom.col_upper_[col];
            implics.addVLB(col, bincol, -static_cast<double>(coef),
                           -static_cast<double>(constant));
          } else {
            constant += globaldom.col_lower_[col];
            implics.addVUB(col, bincol, static_cast<double>(coef),
                           static_cast<double>(constant));
          }
        }
      }
    }
  }

  // only one binary means we do have no cliques
  if (nbin <= 1) return;

  pdqsort(perm.begin(), binaryend, [&](HighsInt p1, HighsInt p2) {
    return std::make_pair(vals[p1], p1) > std::make_pair(vals[p2], p2);
  });
  // check if any cliques exists
  if (vals[perm[0]] + vals[perm[1]] <= rhs + feastol) return;

  // check if this is a set packing constraint (or easily transformable
  // into one)
  if (std::abs(vals[perm[0]] - vals[perm[nbin - 1]]) <= feastol &&
      rhs < 2 * vals[perm[nbin - 1]] - feastol) {
    // the coefficients on the binary variables are all equal and the
    // right hand side is strictly below two times the coefficient value.
    // Therefore the constraint can be transformed into a set packing
    // constraint by relaxing out all non-binary variables (if any),
    // dividing by the coefficient value of the binary variables, and then
    // possibly rounding down the right hand side
    clique.clear();

    for (auto j = 0; j != nbin; ++j) {
      HighsInt pos = perm[j];
      if (complementation[pos] == -1)
        clique.emplace_back(inds[pos], 0);
      else
        clique.emplace_back(inds[pos], 1);
    }

    addClique(mipsolver, clique.data(), nbin);
    // printf("extracted this clique:\n");
    // printClique(clique);
    return;
  }

  for (HighsInt k = nbin - 1; k != 0; --k) {
    double mincliqueval = rhs - vals[perm[k]] + feastol;
    auto cliqueend = std::partition_point(
        perm.begin(), perm.begin() + k,
        [&](HighsInt p) { return vals[p] > mincliqueval; });

    // no clique for this variable
    if (cliqueend == perm.begin()) continue;

    clique.clear();

    for (auto j = perm.begin(); j != cliqueend; ++j) {
      HighsInt pos = *j;
      if (complementation[pos] == -1)
        clique.emplace_back(inds[pos], 0);
      else
        clique.emplace_back(inds[pos], 1);
    }

    if (complementation[perm[k]] == -1)
      clique.emplace_back(inds[perm[k]], 0);
    else
      clique.emplace_back(inds[perm[k]], 1);

    // printf("extracted this clique:\n");
    // printClique(clique);

    if (clique.size() >= 2) {
      // if (clique.size() > 2) runCliqueSubsumption(globaldom, clique);
      // runCliqueMerging(globaldom, clique);
      // if (clique.size() >= 2) {
      addClique(mipsolver, clique.data(), static_cast<HighsInt>(clique.size()));
      if (globaldom.infeasible()) return;
      //}
    }

    // further cliques are just subsets of this clique
    if (cliqueend == perm.begin() + k) return;
  }
}

template <bool Sort, typename Comparator>
HighsInt HighsCliqueTable::extendClique(
    std::vector<HighsInt>& neighbourhoodInds, std::vector<CliqueVar>& clqVars,
    HighsInt seedPos, HighsInt extensionEnd, const Comparator& candidateOrder) {
  for (HighsInt i = seedPos; i < extensionEnd; ++i) {
    if (Sort)
      pdqsort_branchless(clqVars.begin() + i, clqVars.begin() + extensionEnd,
                         candidateOrder);
    CliqueVar v = clqVars[i];
    HighsInt extensionStart = i + 1;
    extensionEnd =
        partitionNeighbourhood(neighbourhoodInds, numNeighbourhoodQueries, v,
                               clqVars.data() + extensionStart,
                               extensionEnd - extensionStart) +
        extensionStart;
  }
  return extensionEnd;
}

HighsInt HighsCliqueTable::extendClique(
    std::vector<HighsInt>& neighbourhoodInds, std::vector<CliqueVar>& clqVars,
    HighsInt seedPos, HighsInt extensionEnd) {
  auto noop = [](const CliqueVar&, const CliqueVar&) { return false; };
  return extendClique<false>(neighbourhoodInds, clqVars, seedPos, extensionEnd,
                             noop);
}

void HighsCliqueTable::cliquePartition(std::vector<CliqueVar>& clqVars,
                                       std::vector<HighsInt>& partitionStart) {
  randgen.shuffle(clqVars.data(), clqVars.size());

  std::vector<HighsInt> neighbourhoodInds;
  neighbourhoodInds.reserve(clqVars.size());

  HighsInt numClqVars = clqVars.size();
  partitionStart.clear();
  partitionStart.reserve(clqVars.size());

  HighsInt pos = 0;
  while (pos < numClqVars) {
    partitionStart.push_back(pos);
    pos = extendClique(neighbourhoodInds, clqVars, pos, numClqVars);
  }
  partitionStart.push_back(numClqVars);
}

void HighsCliqueTable::cliquePartition(const std::vector<double>& objective,
                                       std::vector<CliqueVar>& clqVars,
                                       std::vector<HighsInt>& partitionStart) {
  randgen.shuffle(clqVars.data(), clqVars.size());

  pdqsort_branchless(clqVars.begin(), clqVars.end(),
                     [&](CliqueVar v1, CliqueVar v2) {
                       return (2 * v1.val - 1) * objective[v1.col] >
                              (2 * v2.val - 1) * objective[v2.col];
                     });

  std::vector<HighsInt> neighbourhoodInds;
  neighbourhoodInds.reserve(clqVars.size());

  HighsInt numClqVars = clqVars.size();
  partitionStart.clear();
  partitionStart.reserve(clqVars.size());
  HighsInt extensionEnd = numClqVars;
  partitionStart.push_back(0);
  HighsInt lastSwappedIndex = 0;
  for (HighsInt i = 0; i < numClqVars; ++i) {
    if (i == extensionEnd) {
      partitionStart.push_back(i);
      extensionEnd = numClqVars;
      if (lastSwappedIndex >= i)
        pdqsort_branchless(clqVars.begin() + i,
                           clqVars.begin() + lastSwappedIndex + 1,
                           [&](CliqueVar v1, CliqueVar v2) {
                             return (2 * v1.val - 1) * objective[v1.col] >
                                    (2 * v2.val - 1) * objective[v2.col];
                           });
      lastSwappedIndex = 0;
    }
    CliqueVar v = clqVars[i];
    HighsInt extensionStart = i + 1;
    extensionEnd =
        partitionNeighbourhood(neighbourhoodInds, numNeighbourhoodQueries, v,
                               clqVars.data() + extensionStart,
                               extensionEnd - extensionStart) +
        extensionStart;
    if (!neighbourhoodInds.empty())
      lastSwappedIndex =
          std::max(neighbourhoodInds.back() + extensionStart, lastSwappedIndex);
  }

  partitionStart.push_back(numClqVars);
}

void HighsCliqueTable::cliqueCover(std::vector<CliqueVar>& clqVars,
                                   std::vector<std::vector<CliqueVar>>& cover) {
  cover.clear();
  if (clqVars.empty()) return;

  auto candidateOrder = [&](const CliqueVar& a, const CliqueVar& b) {
    return numcliquesvar[a.index()] > numcliquesvar[b.index()];
  };

  std::vector<HighsInt> neighbourhoodInds;
  neighbourhoodInds.reserve(clqVars.size());

  HighsInt numClqVars = clqVars.size();
  HighsInt pos = 0;

  while (pos < numClqVars) {
    HighsInt cliqueEnd = extendClique(neighbourhoodInds, clqVars, pos,
                                      numClqVars, candidateOrder);

    std::vector<CliqueVar> clique(clqVars.begin() + pos,
                                  clqVars.begin() + cliqueEnd);

    for (HighsInt i = 0; i < pos; ++i) {
      if (std::all_of(clique.begin(), clique.end(),
                      [&](const CliqueVar& clqvar) {
                        return haveCommonClique(clqVars[i], clqvar);
                      }))
        clique.push_back(clqVars[i]);
    }

    cover.push_back(std::move(clique));
    pos = cliqueEnd;
  }

  // count variable occurrences
  HighsHashTable<HighsInt, HighsInt> count;
  for (const auto& clique : cover)
    for (const auto& var : clique) count[var.index()]++;

  // remove cliques where every member appears in at least one other clique
  size_t kept = 0;
  for (auto& clique : cover) {
    if (std::all_of(clique.begin(), clique.end(), [&](const CliqueVar& var) {
          return count[var.index()] >= 2;
        })) {
      for (const auto& var : clique) count[var.index()]--;
    } else {
      std::swap(cover[kept], clique);
      kept++;
    }
  }
  cover.resize(kept);
}

bool HighsCliqueTable::foundCover(HighsDomain& globaldom, CliqueVar v1,
                                  CliqueVar v2) {
  HighsInt commonclique = findCommonCliqueId(v1, v2);
  if (commonclique == -1) return false;

  while (commonclique != -1) {
    HighsInt start = cliques[commonclique].start;
    HighsInt end = cliques[commonclique].end;

    for (HighsInt i = start; i != end; ++i) {
      if (cliqueentries[i] == v1 || cliqueentries[i] == v2 ||
          colDeleted[cliqueentries[i].col])
        continue;
      if (!fixCol(globaldom, cliqueentries[i])) return true;
    }

    removeClique(commonclique);
    commonclique = findCommonCliqueId(v1, v2);
  }

  processInfeasibleVertices(globaldom);
  return true;
}

void HighsCliqueTable::extractCliquesFromCut(const HighsMipSolver& mipsolver,
                                             const HighsInt* inds,
                                             const double* vals, HighsInt len,
                                             double rhs) {
  if (isFull()) return;

  HighsImplications& implics = mipsolver.mipdata_->implications;
  HighsDomain& globaldom = mipsolver.mipdata_->getDomain();

  const double feastol = mipsolver.mipdata_->feastol;

  HighsCDouble minact = 0.0;
  HighsInt nbin = 0;
  for (HighsInt i = 0; i != len; ++i) {
    if (globaldom.isBinary(inds[i])) ++nbin;

    if (vals[i] > 0) {
      if (globaldom.col_lower_[inds[i]] == -kHighsInf) return;
      minact +=
          vals[i] * static_cast<HighsCDouble>(globaldom.col_lower_[inds[i]]);
    } else {
      if (globaldom.col_upper_[inds[i]] == kHighsInf) return;
      minact +=
          vals[i] * static_cast<HighsCDouble>(globaldom.col_upper_[inds[i]]);
    }
  }

  if (rhs - minact < 0.0) {
    minact = rhs;
  }

  for (HighsInt i = 0; i != len; ++i) {
    if (mipsolver.isColContinuous(inds[i])) continue;

    HighsCDouble impliedBound = (rhs - minact) / vals[i];
    const double boundTol =
        std::max(globaldom.feastol(), globaldom.epsilon() / std::abs(vals[i]));
    if (vals[i] > 0) {
      const double boundVal = std::floor(static_cast<double>(
          impliedBound + globaldom.col_lower_[inds[i]] + boundTol));
      globaldom.changeBound(HighsBoundType::kUpper, inds[i], boundVal,
                            HighsDomain::Reason::unspecified());
      if (globaldom.infeasible()) return;
    } else {
      const double boundVal = std::ceil(static_cast<double>(
          impliedBound + globaldom.col_upper_[inds[i]] - boundTol));
      globaldom.changeBound(HighsBoundType::kLower, inds[i], boundVal,
                            HighsDomain::Reason::unspecified());
      if (globaldom.infeasible()) return;
    }
  }

  if (nbin <= 1) return;

  std::vector<HighsInt> perm;
  perm.resize(len);
  std::iota(perm.begin(), perm.end(), 0);

  auto binaryend = std::partition(perm.begin(), perm.end(), [&](HighsInt pos) {
    return globaldom.isBinary(inds[pos]);
  });

  nbin = binaryend - perm.begin();

  // if not all variables are binary, we extract variable upper and lower bounds
  // constraints on the non-binary variable for each binary variable in the
  // constraint:
  if (nbin < len) {
    for (HighsInt i = 0; i != nbin; ++i) {
      HighsInt bincol = inds[perm[i]];
      HighsCDouble impliedActivity = rhs - minact - std::abs(vals[perm[i]]);
      for (HighsInt j = nbin; j != len; ++j) {
        HighsInt col = inds[perm[j]];
        if (globaldom.isFixed(col)) continue;

        if (vals[perm[j]] > 0) {
          double implcolub = static_cast<double>(
              (impliedActivity +
               vals[perm[j]] *
                   static_cast<HighsCDouble>(globaldom.col_lower_[col])) /
              vals[perm[j]]);
          if (mipsolver.isColIntegral(col))
            implcolub = std::floor(implcolub + mipsolver.mipdata_->feastol);

          if (implcolub < globaldom.col_upper_[col] - feastol) {
            double coef;
            double constant;
            if (vals[perm[i]] < 0) {
              coef = globaldom.col_upper_[col] - implcolub;
              constant = implcolub;
            } else {
              // make sure that upper bound is not infinite to avoid adding VUB
              // with coefficient '-kHighsInf' and constant 'kHighsInf'
              if (globaldom.col_upper_[col] == kHighsInf) continue;
              coef = implcolub - globaldom.col_upper_[col];
              constant = globaldom.col_upper_[col];
            }
            // printf("extracted VUB from cut: x%" HIGHSINT_FORMAT " <= %g*y%"
            // HIGHSINT_FORMAT " + %g\n", col, coef,
            //        bincol, constant);
            implics.addVUB(col, bincol, coef, constant);
          }
        } else {
          double implcollb = static_cast<double>(
              (impliedActivity +
               vals[perm[j]] *
                   static_cast<HighsCDouble>(globaldom.col_upper_[col])) /
              vals[perm[j]]);
          if (mipsolver.isColIntegral(col))
            implcollb = std::ceil(implcollb - mipsolver.mipdata_->feastol);

          if (implcollb > globaldom.col_lower_[col] + feastol) {
            double coef;
            double constant;
            if (vals[perm[i]] < 0) {
              coef = globaldom.col_lower_[col] - implcollb;
              constant = implcollb;
            } else {
              // make sure that lower bound is not infinite to avoid adding VLB
              // with coefficient 'kHighsInf' and constant '-kHighsInf'
              if (globaldom.col_lower_[col] == -kHighsInf) continue;
              coef = implcollb - globaldom.col_lower_[col];
              constant = globaldom.col_lower_[col];
            }

            // printf("extracted VLB from cut: x%" HIGHSINT_FORMAT " >= %g*y%"
            // HIGHSINT_FORMAT " + %g\n", col, coef,
            //        bincol, constant);
            implics.addVLB(col, bincol, coef, constant);
            // printf("extracted VLB from cut\n");
          }
        }
      }
    }
  }

  // only one binary means we do have no cliques
  if (nbin <= 1) return;

  std::vector<CliqueVar> clique;
  clique.reserve(nbin);

  pdqsort(perm.begin(), binaryend, [&](HighsInt p1, HighsInt p2) {
    return std::make_pair(std::abs(vals[p1]), p1) >
           std::make_pair(std::abs(vals[p2]), p2);
  });
  // check if any cliques exists
  if (std::abs(vals[perm[0]]) + std::abs(vals[perm[1]]) <=
      static_cast<double>(rhs - minact + feastol))
    return;

  HighsInt maxNewEntries =
      std::min(mipsolver.mipdata_->numCliqueEntriesAfterPresolve + 100000 +
                   4 * globaldom.numModelNonzeros(),
               numEntries + 10 * nbin);

  for (HighsInt k = nbin - 1; k != 0 && numEntries < maxNewEntries; --k) {
    double mincliqueval =
        static_cast<double>(rhs - minact - std::abs(vals[perm[k]]) + feastol);
    auto cliqueend = std::partition_point(
        perm.begin(), perm.begin() + k,
        [&](HighsInt p) { return std::abs(vals[p]) > mincliqueval; });

    // no clique for this variable
    if (cliqueend == perm.begin()) continue;

    clique.clear();

    for (auto j = perm.begin(); j != cliqueend; ++j) {
      HighsInt pos = *j;
      if (vals[pos] < 0)
        clique.emplace_back(inds[pos], 0);
      else
        clique.emplace_back(inds[pos], 1);
    }

    if (vals[perm[k]] < 0)
      clique.emplace_back(inds[perm[k]], 0);
    else
      clique.emplace_back(inds[perm[k]], 1);

    // printf("extracted this clique:\n");
    // printClique(clique);
    if (clique.size() >= 2) {
      // printf("extracted clique from cut\n");
      // if (clique.size() > 2) runCliqueSubsumption(globaldom, clique);

      addClique(mipsolver, clique.data(), static_cast<HighsInt>(clique.size()));
      if (globaldom.infeasible() || numEntries >= maxNewEntries) return;
    }

    // further cliques are just subsets of this clique
    if (cliqueend == perm.begin() + k) return;
  }
}

void HighsCliqueTable::extractCliques(HighsMipSolver& mipsolver,
                                      bool transformRows) {
  std::vector<HighsInt> inds;
  std::vector<double> vals;
  std::vector<HighsInt> perm;
  std::vector<int8_t> complementation;
  std::vector<CliqueVar> clique;
  HighsHashTable<HighsInt, double> entries;
  double offset;

  double rhs;

  HighsDomain& globaldom = mipsolver.mipdata_->getDomain();

  for (HighsInt i : mipsolver.mipdata_->postSolveStack.getNonCutRows()) {
    HighsInt start = mipsolver.mipdata_->ARstart_[i];
    HighsInt end = mipsolver.mipdata_->ARstart_[i + 1];

    // catch set packing and partitioning constraints that already have the
    // form of a clique without transformations and add those cliques with
    // the rows being recorded. only <= and = rows are checked because
    // normaliseCliqueRows has already flipped >= rows to <= form.
    if (mipsolver.rowUpper(i) < kHighsInf) {
      bool issetppc = true;
      HighsCDouble fixedTerm = 0.0;
      HighsInt numComp = 0;
      for (HighsInt j = start; j != end; ++j) {
        HighsInt col = mipsolver.mipdata_->ARindex_[j];
        double val = mipsolver.mipdata_->ARvalue_[j];

        // handle fixed non-binary variables
        if (!globaldom.isBinary(col) && globaldom.isFixed(col)) {
          fixedTerm +=
              val * static_cast<HighsCDouble>(globaldom.col_upper_[col]);
          continue;
        }

        // check if we have a set partitioning / packing row
        issetppc = globaldom.isBinary(col) && std::abs(val) == 1.0;
        if (!issetppc) break;

        // count number of complemented binaries (with coefficient -1)
        if (val < 0) numComp++;
      }

      if (issetppc) {
        // subtract fixed term and round right-hand side
        double rhs = static_cast<double>(floor(
            mipsolver.rowUpper(i) - fixedTerm + mipsolver.mipdata_->feastol));

        if (rhs == 1.0 - numComp) {
          // subtract fixed term and round left-hand side if finite
          double lhs = mipsolver.rowLower(i);
          if (lhs > -kHighsInf)
            lhs = static_cast<double>(
                ceil(lhs - fixedTerm - mipsolver.mipdata_->feastol));

          clique.clear();

          for (HighsInt j = start; j != end; ++j) {
            HighsInt col = mipsolver.mipdata_->ARindex_[j];
            double val = mipsolver.mipdata_->ARvalue_[j];
            HighsInt dir = val > 0 ? 1 : 0;

            // skip non-binary variables (fixed, see previous loop) and binaries
            // that are fixed to "inactive" values
            if (!globaldom.isBinary(col) ||
                globaldom.isFixedToVal(col, 1 - dir))
              continue;

            // add to clique
            clique.emplace_back(col, dir);
          }

          // add clique to clique table
          if (clique.size() >= 2) {
            addClique(mipsolver, clique.data(),
                      static_cast<HighsInt>(clique.size()), rhs == lhs, i);
            if (globaldom.infeasible()) return;
          }
          continue;
        }
      }
    }

    if (!transformRows || isFull()) continue;

    offset = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = mipsolver.mipdata_->ARindex_[j];
      double val = mipsolver.mipdata_->ARvalue_[j];

      resolveSubstitution(col, val, offset);
      entries[col] += val;
    }

    auto checkRow = [&](double rhs, HighsInt direction) {
      if (direction * rhs == kHighsInf) return;
      rhs = direction * (rhs - offset);
      inds.clear();
      vals.clear();
      complementation.clear();
      bool freevar = false;
      HighsInt nbin = 0;

      for (const auto& entry : entries) {
        HighsInt col = entry.key();
        double val = direction * entry.value();

        if (std::abs(val) < mipsolver.mipdata_->epsilon) continue;

        if (globaldom.isBinary(col)) ++nbin;

        if (val < 0) {
          freevar = globaldom.col_upper_[col] == kHighsInf;
          if (freevar) break;

          vals.push_back(-val);
          inds.push_back(col);
          complementation.push_back(-1);
          rhs -= val * globaldom.col_upper_[col];
        } else {
          freevar = globaldom.col_lower_[col] == -kHighsInf;
          if (freevar) break;

          vals.push_back(val);
          inds.push_back(col);
          complementation.push_back(1);
          rhs -= val * globaldom.col_lower_[col];
        }
      }

      if (!freevar && nbin != 0) {
        extractCliques(mipsolver, inds, vals, complementation, rhs, nbin, perm,
                       clique, mipsolver.mipdata_->feastol);
        if (globaldom.infeasible()) return;
      }
    };

    checkRow(mipsolver.rowUpper(i), HighsInt{1});
    checkRow(mipsolver.rowLower(i), HighsInt{-1});

    entries.clear();
  }
}

void HighsCliqueTable::extractObjCliques(HighsMipSolver& mipsolver) {
  HighsInt nbin =
      mipsolver.mipdata_->objectiveFunction.getNumBinariesInObjective();
  if (nbin <= 1) return;
  HighsDomain& globaldom = mipsolver.mipdata_->getDomain();
  if (globaldom.getObjectiveLowerBound() == -kHighsInf) return;

  const double* vals;
  const HighsInt* inds;
  HighsInt len;
  double rhs;
  globaldom.getCutoffConstraint(vals, inds, len, rhs);

  std::vector<HighsInt> perm;
  perm.resize(nbin);
  std::iota(perm.begin(), perm.end(), 0);

  auto binaryend = std::partition(perm.begin(), perm.end(), [&](HighsInt pos) {
    return vals[pos] != 0.0 && !globaldom.isFixed(inds[pos]);
  });

  nbin = binaryend - perm.begin();

  // only one binary means we do have no cliques
  if (nbin <= 1) return;

  std::vector<CliqueVar> clique;
  clique.reserve(nbin);

  pdqsort(perm.begin(), binaryend, [&](HighsInt p1, HighsInt p2) {
    return std::make_pair(std::fabs(vals[p1]), p1) >
           std::make_pair(std::fabs(vals[p2]), p2);
  });

  // check if any cliques exists
  HighsCDouble minact;
  HighsInt ninf;
  globaldom.computeMinActivity(0, len, inds, vals, ninf, minact);
  const double feastol = mipsolver.mipdata_->feastol;
  if (std::fabs(vals[perm[0]]) + std::fabs(vals[perm[1]]) <=
      static_cast<double>(rhs - minact + feastol))
    return;

  for (HighsInt k = nbin - 1; k != 0; --k) {
    double mincliqueval =
        static_cast<double>(rhs - minact - std::fabs(vals[perm[k]]) + feastol);
    auto cliqueend = std::partition_point(
        perm.begin(), perm.begin() + k,
        [&](HighsInt p) { return std::abs(vals[p]) > mincliqueval; });

    // no clique for this variable
    if (cliqueend == perm.begin()) continue;

    clique.clear();

    for (auto j = perm.begin(); j != cliqueend; ++j) {
      HighsInt pos = *j;
      if (vals[pos] < 0)
        clique.emplace_back(inds[pos], 0);
      else
        clique.emplace_back(inds[pos], 1);
    }

    if (vals[perm[k]] < 0)
      clique.emplace_back(inds[perm[k]], 0);
    else
      clique.emplace_back(inds[perm[k]], 1);

    // printf("extracted this clique from obj:\n");
    // printClique(clique);
    if (clique.size() >= 2) {
      // printf("extracted clique from obj\n");
      // if (clique.size() > 2) runCliqueSubsumption(globaldom, clique);

      addClique(mipsolver, clique.data(), static_cast<HighsInt>(clique.size()));
      if (globaldom.infeasible()) return;
    }

    // further cliques are just subsets of this clique
    if (cliqueend == perm.begin() + k) return;
  }
}

void HighsCliqueTable::processInfeasibleVertices(HighsDomain& globaldom) {
  while (!infeasvertexstack.empty() && !globaldom.infeasible()) {
    CliqueVar v = infeasvertexstack.back().complement();
    infeasvertexstack.pop_back();

    resolveSubstitution(v);
    bool wasfixed = globaldom.isFixed(v.col);
    globaldom.fixCol(v.col, static_cast<double>(v.val));
    if (globaldom.infeasible()) return;
    if (!wasfixed) ++nfixings;
    if (colDeleted[v.col]) continue;
    colDeleted[v.col] = true;

    HighsHashTree<HighsInt, HighsInt> vHashLists =
        std::move(invertedHashList[v.index()]);
    HighsHashTree<HighsInt> vHashListsSizeTwo =
        std::move(invertedHashListSizeTwo[v.index()]);

    bool infeas = vHashLists.for_each([&](HighsInt cliqueid) {
      HighsInt start = cliques[cliqueid].start;
      HighsInt end = cliques[cliqueid].end;

      for (HighsInt i = start; i != end; ++i) {
        if (cliqueentries[i].col == v.col || colDeleted[cliqueentries[i].col])
          continue;
        if (!fixCol(globaldom, cliqueentries[i])) return true;
      }

      removeClique(cliqueid);
      return false;
    });

    if (infeas) return;

    infeas = vHashListsSizeTwo.for_each([&](HighsInt cliqueid) {
      HighsInt start = cliques[cliqueid].start;
      HighsInt end = cliques[cliqueid].end;

      for (HighsInt i = start; i != end; ++i) {
        if (cliqueentries[i].col == v.col || colDeleted[cliqueentries[i].col])
          continue;
        if (!fixCol(globaldom, cliqueentries[i])) return true;
      }

      removeClique(cliqueid);
      return false;
    });

    if (infeas) return;

    vHashLists = std::move(invertedHashList[v.complement().index()]);
    vHashListsSizeTwo =
        std::move(invertedHashListSizeTwo[v.complement().index()]);

    vHashListsSizeTwo.for_each([&](HighsInt cliqueid) {
      ++cliques[cliqueid].numZeroFixed;
      fixLastActiveAndRemove(globaldom, cliqueid);
    });

    assert(cliquehitinds.empty());
    std::vector<CliqueVar> clq;
    vHashLists.for_each([&](HighsInt cliqueid) {
      // assert(cliqueentries[entry.value()].val == 1 - v.val);

      cliques[cliqueid].numZeroFixed += 1;
      if (cliques[cliqueid].numActive() <= 1) {
        fixLastActiveAndRemove(globaldom, cliqueid);
      } else if (!inPresolveProbing &&
                 cliques[cliqueid].numZeroFixed >=
                     std::max(HighsInt{10}, (cliques[cliqueid].end -
                                             cliques[cliqueid].start) >>
                                                1)) {
        const bool equality = cliques[cliqueid].equality;
        const HighsInt origin = cliques[cliqueid].origin;
        clq.assign(cliqueentries.begin() + cliques[cliqueid].start,
                   cliqueentries.begin() + cliques[cliqueid].end);

        auto computeNumDeleted = [&](const std::vector<CliqueVar>& clq) {
          HighsInt numDel = 0;
          for (CliqueVar x : clq) numDel += colDeleted[x.col];
          return numDel;
        };

        assert(computeNumDeleted(clq) == cliques[cliqueid].numZeroFixed);

        removeClique(cliqueid, false);
        clq.erase(
            std::remove_if(clq.begin(), clq.end(),
                           [&](CliqueVar x) {
                             return (globaldom.isFixedToVal(x.col, 1 - x.val) ||
                                     colDeleted[x.col]);
                           }),
            clq.end());
        if (clq.size() > 1) {
          doAddClique(clq.data(), clq.size(), equality, origin);
        } else if (equality) {
          if (clq.empty() ||
              globaldom.isFixedToVal(static_cast<HighsInt>(clq[0].col),
                                     1 - clq[0].val)) {
            globaldom.markInfeasible();
          } else if (!globaldom.isFixed(static_cast<HighsInt>(clq[0].col))) {
            fixCol(globaldom, clq[0].complement(), false);
          }
        }
      }
    });
  }

  propagateAndCleanup(globaldom);
}

void HighsCliqueTable::propagateAndCleanup(HighsDomain& globaldom) {
  const auto& domchgstack = globaldom.getDomainChangeStack();
  HighsInt start = domchgstack.size();
  globaldom.propagate();
  HighsInt end = domchgstack.size();

  while (!globaldom.infeasible() && start != end) {
    for (HighsInt k = start; k != end; ++k) {
      HighsInt col = domchgstack[k].column;
      if (!globaldom.isFixed(col)) continue;
      if (globaldom.col_lower_[col] != 1.0 && globaldom.col_lower_[col] != 0.0)
        continue;

      HighsInt fixval = static_cast<HighsInt>(globaldom.col_lower_[col]);
      CliqueVar v(col, 1 - fixval);
      if (numCliques(v) != 0) {
        vertexInfeasible(globaldom, col, 1 - fixval);
        if (globaldom.infeasible()) return;
      }
    }
    start = domchgstack.size();
    globaldom.propagate();
    end = domchgstack.size();
  }
}

void HighsCliqueTable::vertexInfeasible(HighsDomain& globaldom, HighsInt col,
                                        HighsInt val) {
  bool wasfixed = globaldom.isFixed(col);
  globaldom.fixCol(col, static_cast<double>(1 - val));
  if (globaldom.infeasible()) return;
  if (!wasfixed) ++nfixings;
  infeasvertexstack.emplace_back(col, val);
  processInfeasibleVertices(globaldom);
}

void HighsCliqueTable::separateCliques(const HighsMipSolver& mipsolver,
                                       const std::vector<double>& sol,
                                       HighsCutPool& cutpool, double feastol,
                                       HighsRandom& randgen,
                                       int64_t& localNumNeighbourhoodQueries) {
  BronKerboschData data(sol);
  data.feastol = feastol;
  data.maxNeighbourhoodQueries = 1000000 +
                                 int64_t{100} * mipsolver.numNonzero() +
                                 mipsolver.mipdata_->total_lp_iterations * 1000;
  if (numNeighbourhoodQueries > data.maxNeighbourhoodQueries) return;
  data.maxNeighbourhoodQueries -= numNeighbourhoodQueries;
  const HighsDomain& globaldom = mipsolver.mipdata_->getDomain();

  for (HighsInt i : mipsolver.mipdata_->integral_cols) {
    if (colsubstituted[i] || colDeleted[i]) continue;
#ifdef ADD_ZERO_WEIGHT_VARS
    if (numCliques(i, 0) != 0) {
      if (CliqueVar(i, 0).weight(sol) > feastol)
        data.P.emplace_back(i, 0);
      else
        data.Z.emplace_back(i, 0);
    }
    if (numCliques(i, 1) != 0) {
      if (CliqueVar(i, 1).weight(sol) > feastol)
        data.P.emplace_back(i, 1);
      else
        data.Z.emplace_back(i, 1);
    }
#else
    if (numcliquesvar[CliqueVar(i, 0).index()] != 0 &&
        CliqueVar(i, 0).weight(sol) > feastol)
      data.P.emplace_back(i, 0);
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0 &&
        CliqueVar(i, 1).weight(sol) > feastol)
      data.P.emplace_back(i, 1);
#endif
  }

  // auto t1 = std::chrono::high_resolution_clock::now();
  bronKerboschRecurse(data, data.P.size(), nullptr, 0);

  // auto t2 = std::chrono::high_resolution_clock::now();

  // printf(
  //     "bron kerbosch: %" HIGHSINT_FORMAT " calls, %" HIGHSINT_FORMAT "
  //     cliques, %ldms\n", data.ncalls, int(data.cliques.size()),
  //     std::chrono::duration_cast<std::chrono::milliseconds>(t2 -
  //     t1).count());

  bool runcliquesubsumption = false;
  std::vector<HighsInt> inds;
  std::vector<double> vals;
  for (std::vector<CliqueVar>& clique : data.cliques) {
#ifdef ADD_ZERO_WEIGHT_VARS
    HighsInt extensionend = static_cast<HighsInt>(data.Z.size());
    for (CliqueVar v : clique) {
      extensionend = partitionNeighbourhood(data.neighbourhoodInds,
                                            data.numNeighbourhoodQueries, v,
                                            data.Z.data(), extensionend);
      if (extensionend == 0) break;
    }

    if (extensionend != 0) {
      randgen.shuffle(data.Z.data(), extensionend);

      for (HighsInt i = 0; i < extensionend; ++i) {
        HighsInt k = i + 1;
        extensionend =
            k + partitionNeighbourhood(data.neighbourhoodInds,
                                       data.numNeighbourhoodQueries, data.Z[i],
                                       data.Z.data() + k, extensionend - k);
      }

      clique.insert(clique.end(), data.Z.begin(),
                    data.Z.begin() + extensionend);
    }
#endif

    double rhs = 1;
    runcliquesubsumption = cliques.size() > 2;
    inds.clear();
    vals.clear();

    for (CliqueVar v : clique) {
      inds.push_back(v.col);
      if (v.val == 0) {
        vals.push_back(-1);
        rhs -= 1;
      } else
        vals.push_back(1);
    }

    rhs = std::floor(rhs + 0.5);

    cutpool.addCut(mipsolver, inds.data(), vals.data(), inds.size(), rhs, true,
                   false, false);
  }

  localNumNeighbourhoodQueries += data.numNeighbourhoodQueries;

  if (runcliquesubsumption && &randgen == &this->randgen) {
    for (std::vector<CliqueVar>& clique : data.cliques) {
      HighsInt nremoved = runCliqueSubsumption(globaldom, clique);

      if (clique.empty()) continue;
      if (nremoved != 0) doAddClique(clique.data(), clique.size(), false, -1);
    }
  }
}

std::vector<std::vector<HighsCliqueTable::CliqueVar>>
HighsCliqueTable::separateCliques(const std::vector<double>& sol,
                                  const HighsDomain& globaldom,
                                  double feastol) {
  exit(0);
#if 0
  BronKerboschData data(sol);
  data.feastol = feastol;

  HighsInt numcols = globaldom.col_lower_.size();
  assert(int(numcliquesvar.size()) == 2 * numcols);
  for (HighsInt i = 0; i != numcols; ++i) {
    if (colsubstituted[i]) continue;

    if (numcliquesvar[CliqueVar(i, 0).index()] != 0 &&
        CliqueVar(i, 0).weight(sol) > feastol)
      data.P.emplace_back(i, 0);
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0 &&
        CliqueVar(i, 1).weight(sol) > feastol)
      data.P.emplace_back(i, 1);
  }

  bronKerboschRecurse(data, data.P.size(), nullptr, 0);

  return std::move(data.cliques);
#endif
}

std::vector<std::vector<HighsCliqueTable::CliqueVar>>
HighsCliqueTable::computeMaximalCliques(const std::vector<CliqueVar>& vars,
                                        double feastol) {
  // return if there are no variables
  if (vars.empty())
    return std::vector<std::vector<HighsCliqueTable::CliqueVar>>{};

  // find max column index in clique variables
  size_t maxcolindex = 0;
  for (const auto& var : vars)
    maxcolindex = std::max(static_cast<size_t>(var.col), maxcolindex);

  // set up data
  std::vector<double> sol;
  sol.resize(maxcolindex + 1);
  for (const auto& var : vars) sol[var.col] = var.val;
  BronKerboschData data(sol);
  data.feastol = feastol;

  for (const auto& var : vars) {
    if (colsubstituted[var.col] || colDeleted[var.col]) continue;
    if (numCliques(var) != 0) data.P.emplace_back(var);
  }

  bronKerboschRecurse(data, data.P.size(), nullptr, 0);

  return std::move(data.cliques);
}

void HighsCliqueTable::addImplications(HighsDomain& domain, HighsInt col,
                                       HighsInt val) {
  CliqueVar v(col, val);

  while (colsubstituted[v.col]) {
    assert(static_cast<HighsInt>(substitutions.size()) >
           colsubstituted[v.col] - 1);
    Substitution subst = substitutions[colsubstituted[v.col] - 1];
    v = v.val == 1 ? subst.replace : subst.replace.complement();
    if (v.val == 1) {
      if (domain.col_lower_[v.col] == 1.0) continue;

      domain.changeBound(HighsBoundType::kLower, v.col, 1.0,
                         HighsDomain::Reason::cliqueTable(col, val));
      if (domain.infeasible()) return;
    } else {
      if (domain.col_upper_[v.col] == 0.0) continue;

      domain.changeBound(HighsBoundType::kUpper, v.col, 0.0,
                         HighsDomain::Reason::cliqueTable(col, val));
      if (domain.infeasible()) return;
    }
  }

  auto doFixings = [&](HighsInt cliqueid) {
    HighsInt start = cliques[cliqueid].start;
    HighsInt end = cliques[cliqueid].end;

    for (HighsInt i = start; i != end; ++i) {
      if (cliqueentries[i].col == v.col || colDeleted[cliqueentries[i].col])
        continue;

      if (cliqueentries[i].val == 1) {
        if (domain.col_upper_[cliqueentries[i].col] == 0.0) continue;

        domain.changeBound(HighsBoundType::kUpper, cliqueentries[i].col, 0.0,
                           HighsDomain::Reason::cliqueTable(col, val));
        if (domain.infeasible()) return true;
      } else {
        if (domain.col_lower_[cliqueentries[i].col] == 1.0) continue;

        domain.changeBound(HighsBoundType::kLower, cliqueentries[i].col, 1.0,
                           HighsDomain::Reason::cliqueTable(col, val));
        if (domain.infeasible()) return true;
      }
    }

    return false;
  };

  bool infeas = invertedHashList[v.index()].for_each(doFixings);

  if (infeas) return;

  invertedHashListSizeTwo[v.index()].for_each(doFixings);
}

void HighsCliqueTable::cleanupFixed(HighsDomain& globaldom) {
  HighsInt numcol = globaldom.col_upper_.size();
  HighsInt oldnfixings = nfixings;
  for (HighsInt i = 0; i != numcol; ++i) {
    if (colDeleted[i] || !globaldom.isFixed(i)) continue;
    if (globaldom.col_lower_[i] != 1.0 && globaldom.col_lower_[i] != 0.0)
      continue;

    HighsInt fixval = static_cast<HighsInt>(globaldom.col_lower_[i]);
    CliqueVar v(i, 1 - fixval);

    vertexInfeasible(globaldom, v.col, v.val);
    if (globaldom.infeasible()) return;
  }

  if (nfixings != oldnfixings) propagateAndCleanup(globaldom);
}

HighsInt HighsCliqueTable::getNumImplications(HighsInt col) const {
  // first count all cliques as one implication, so that cliques of size two
  // are accounted for already
  HighsInt i0 = CliqueVar(col, 0).index();
  HighsInt i1 = CliqueVar(col, 1).index();
  HighsInt numimplics = numcliquesvar[i0] + numcliquesvar[i1];

  // now loop over cliques larger than size two and add the cliquelength - 1 as
  // additional implications
  auto countImplics = [&](HighsInt cliqueid) {
    HighsInt nimplics = cliques[cliqueid].numActive() - 1;
    nimplics *= (1 + cliques[cliqueid].equality);
    numimplics += nimplics - 1;
  };
  invertedHashList[i0].for_each(countImplics);
  invertedHashList[i1].for_each(countImplics);
  return numimplics;
}

HighsInt HighsCliqueTable::getNumImplications(HighsInt col, bool val) const {
  HighsInt iVal = CliqueVar(col, val).index();

  // each size two clique is one implication
  HighsInt numimplics = numcliquesvar[iVal];

  // now loop over cliques larger than size two and add the cliquelength - 1 as
  // additional implications
  invertedHashList[iVal].for_each([&](HighsInt cliqueid) {
    HighsInt nimplics = cliques[cliqueid].numActive() - 1;
    nimplics *= (1 + cliques[cliqueid].equality);
    numimplics += nimplics - 1;
  });

  return numimplics;
}

void HighsCliqueTable::runCliqueMerging(HighsDomain& globaldomain,
                                        std::vector<CliqueVar>& clique,
                                        bool equation) {
  CliqueVar extensionstart;
  HighsInt numcliques = kHighsIInf;
  iscandidate.resize(invertedHashList.size());
  std::vector<HighsInt> neighbourhoodInds;
  neighbourhoodInds.reserve(invertedHashList.size());

  HighsInt initialCliqueSize = clique.size();
  for (HighsInt i = 0; i != initialCliqueSize; ++i) {
    if (globaldomain.isFixed(clique[i].col)) continue;

    HighsInt thisNumClqs = numCliques(clique[i]);
    if (thisNumClqs < numcliques) {
      numcliques = thisNumClqs;
      extensionstart = clique[i];
    }
  }

  if (numcliques == kHighsIInf) {
    clique.clear();
    return;
  }

  for (HighsInt i = 0; i != initialCliqueSize; ++i)
    iscandidate[clique[i].index()] = true;

  auto addCands = [&](HighsInt cliqueid) {
    HighsInt start = cliques[cliqueid].start;
    HighsInt end = cliques[cliqueid].end;

    for (HighsInt i = start; i != end; ++i) {
      if (iscandidate[cliqueentries[i].index()] ||
          globaldomain.isFixed(cliqueentries[i].col))
        continue;

      iscandidate[cliqueentries[i].index()] = true;
      clique.push_back(cliqueentries[i]);
    }
  };

  invertedHashList[extensionstart.index()].for_each(
      [&](HighsInt cliqueid) { addCands(cliqueid); });

  invertedHashListSizeTwo[extensionstart.index()].for_each(
      [&](HighsInt cliqueid) { addCands(cliqueid); });

  HighsInt sizeWithCandidates = clique.size();
  for (HighsInt i = 0; i != sizeWithCandidates; ++i)
    iscandidate[clique[i].index()] = false;

  for (HighsInt i = 0; i != initialCliqueSize &&
                       initialCliqueSize < static_cast<HighsInt>(clique.size());
       ++i) {
    if (clique[i] == extensionstart) continue;

    HighsInt newSize =
        initialCliqueSize +
        shrinkToNeighbourhood(neighbourhoodInds, numNeighbourhoodQueries,
                              clique[i], clique.data() + initialCliqueSize,
                              clique.size() - initialCliqueSize);
    clique.erase(clique.begin() + newSize, clique.end());
  }

  if (static_cast<size_t>(initialCliqueSize) < clique.size()) {
    // todo, shuffle extension vars?
    randgen.shuffle(clique.data() + initialCliqueSize,
                    clique.size() - initialCliqueSize);
    HighsInt i = initialCliqueSize;
    while (i < static_cast<HighsInt>(clique.size())) {
      CliqueVar extvar = clique[i];
      i += 1;

      HighsInt newSize = i + shrinkToNeighbourhood(
                                 neighbourhoodInds, numNeighbourhoodQueries,
                                 extvar, clique.data() + i, clique.size() - i);
      clique.erase(clique.begin() + newSize, clique.end());
    }
  }

  if (equation) {
    for (HighsInt i = initialCliqueSize;
         i < static_cast<HighsInt>(clique.size()); ++i)
      vertexInfeasible(globaldomain, clique[i].col, clique[i].val);
  } else {
    runCliqueSubsumption(globaldomain, clique);

    if (!clique.empty()) {
      clique.erase(std::remove_if(clique.begin(), clique.end(),
                                  [&](CliqueVar v) {
                                    return globaldomain.isFixedToVal(v.col,
                                                                     1 - v.val);
                                  }),
                   clique.end());
    }
  }

  processInfeasibleVertices(globaldomain);
}

void HighsCliqueTable::runCliqueMerging(HighsDomain& globaldomain) {
  std::vector<CliqueVar> extensionvars;
  iscandidate.resize(invertedHashList.size());
  std::vector<HighsInt> neighbourhoodInds;
  neighbourhoodInds.reserve(invertedHashList.size());

  HighsInt numcliqueslots = static_cast<HighsInt>(cliques.size());
  const HighsInt maxNewEntries = numEntries + globaldomain.numModelNonzeros();
  bool haveNonModelCliquesToMerge = false;
  for (HighsInt k = 0; k != numcliqueslots; ++k) {
    if (cliques[k].start == -1) continue;
    if (!cliques[k].equality && cliques[k].origin == kHighsIInf) continue;
    if (cliques[k].origin == -1) {
      haveNonModelCliquesToMerge = true;
      continue;
    }
    HighsInt numclqvars = cliques[k].end - cliques[k].start;
    assert(numclqvars != 0);
    if (numclqvars == 0) continue;

    CliqueVar* clqvars = &cliqueentries[cliques[k].start];

    CliqueVar extensionstart = clqvars[0];
    HighsInt numcliques = numCliques(clqvars[0]);
    for (HighsInt i = 1; i != numclqvars; ++i) {
      HighsInt thisNumClqs = numCliques(clqvars[i]);
      if (thisNumClqs < numcliques) {
        numcliques = thisNumClqs;
        extensionstart = clqvars[i];
      }
    }

    for (HighsInt i = 0; i != numclqvars; ++i)
      iscandidate[clqvars[i].index()] = true;

    auto addCands = [&](HighsInt cliqueid) {
      HighsInt start = cliques[cliqueid].start;
      HighsInt end = cliques[cliqueid].end;

      for (HighsInt i = start; i != end; ++i) {
        if (iscandidate[cliqueentries[i].index()] ||
            globaldomain.isFixed(cliqueentries[i].col))
          continue;

        iscandidate[cliqueentries[i].index()] = true;
        extensionvars.push_back(cliqueentries[i]);
      }
    };

    invertedHashList[extensionstart.index()].for_each(
        [&](HighsInt cliqueid) { addCands(cliqueid); });

    invertedHashListSizeTwo[extensionstart.index()].for_each(
        [&](HighsInt cliqueid) { addCands(cliqueid); });

    for (HighsInt i = 0; i != numclqvars; ++i)
      iscandidate[clqvars[i].index()] = false;
    for (CliqueVar v : extensionvars) iscandidate[v.index()] = false;

    for (HighsInt i = 0; i != numclqvars && !extensionvars.empty(); ++i) {
      if (clqvars[i] == extensionstart) continue;

      HighsInt newSize = shrinkToNeighbourhood(
          neighbourhoodInds, numNeighbourhoodQueries, clqvars[i],
          extensionvars.data(), static_cast<HighsInt>(extensionvars.size()));
      extensionvars.erase(extensionvars.begin() + newSize, extensionvars.end());
    }

    if (!extensionvars.empty()) {
      // todo, shuffle extension vars?
      randgen.shuffle(extensionvars.data(),
                      static_cast<HighsInt>(extensionvars.size()));
      HighsInt i = 0;
      while (i < static_cast<HighsInt>(extensionvars.size())) {
        CliqueVar extvar = extensionvars[i];
        i += 1;

        HighsInt newSize =
            i + shrinkToNeighbourhood(
                    neighbourhoodInds, numNeighbourhoodQueries, extvar,
                    extensionvars.data() + i,
                    static_cast<HighsInt>(extensionvars.size()) - i);
        extensionvars.erase(extensionvars.begin() + newSize,
                            extensionvars.end());
      }
    }

    if (cliques[k].equality) {
      for (CliqueVar v : extensionvars)
        vertexInfeasible(globaldomain, v.col, v.val);
    } else {
      HighsInt originrow = cliques[k].origin;
      cliques[k].origin = kHighsIInf;

      size_t numExtensions = extensionvars.size();
      extensionvars.insert(extensionvars.end(),
                           cliqueentries.begin() + cliques[k].start,
                           cliqueentries.begin() + cliques[k].end);
      extensionvars.erase(
          std::remove_if(
              extensionvars.begin() + numExtensions, extensionvars.end(),
              [&](CliqueVar clqvar) { return colDeleted[clqvar.col]; }),
          extensionvars.end());
      removeClique(k);

      // callback for removing cliques
      auto removeCliques = [&](HighsInt cliqueid) { removeClique(cliqueid); };

      // check extension
      bool redundant;
      HighsInt dominatingOrigin;
      cliqueSubsumption(extensionvars, redundant, dominatingOrigin,
                        removeCliques);

      if (!redundant) {
        for (size_t i = 0; i < numExtensions; ++i)
          cliqueextensions.emplace_back(originrow, extensionvars[i]);

        extensionvars.erase(
            std::remove_if(extensionvars.begin(), extensionvars.end(),
                           [&](CliqueVar v) {
                             return globaldomain.isFixedToVal(v.col, 1 - v.val);
                           }),
            extensionvars.end());

        if (extensionvars.size() > 1)
          doAddClique(extensionvars.data(),
                      static_cast<HighsInt>(extensionvars.size()), false,
                      originrow);
      } else {
        // the extended clique is redundant, check if the row can be removed
        if (dominatingOrigin != kHighsIInf)
          deletedrows.push_back(originrow);
        else {
          // this clique is redundant in the cliquetable but its row is not
          // necessarily. Also there might be rows that have been deleted due to
          // being dominated by this row after adding the lifted entries so they
          // must be added to the cliqueextension vector
          for (size_t i = 0; i < numExtensions; ++i)
            cliqueextensions.emplace_back(originrow, extensionvars[i]);
        }
      }
    }

    extensionvars.clear();
    processInfeasibleVertices(globaldomain);

    if (numEntries >= maxNewEntries) break;
    // printf("nonzeroDelta: %d, maxNonzeroDelta: %d\n", nonzeroDelta,
    // maxNonzeroDelta);
  }

  if (haveNonModelCliquesToMerge) {
    for (HighsInt k = 0; k != numcliqueslots; ++k) {
      if (cliques[k].start == -1) continue;
      if (cliques[k].origin != -1) continue;
      // if (cliques[k].end - cliques[k].start <= 1000) continue;

      // printf("numEntries before: %d\n", numEntries);
      extensionvars.clear();
      extensionvars.insert(extensionvars.end(),
                           cliqueentries.begin() + cliques[k].start,
                           cliqueentries.begin() + cliques[k].end);
      removeClique(k);
      runCliqueMerging(globaldomain, extensionvars);

      if (extensionvars.size() > 1)
        doAddClique(extensionvars.data(),
                    static_cast<HighsInt>(extensionvars.size()));
    }
  }
}

void HighsCliqueTable::rebuild(
    HighsInt ncols, const presolve::HighsPostsolveStack& postSolveStack,
    const HighsDomain& globaldomain,
    const std::vector<HighsInt>& orig2reducedcol,
    const std::vector<HighsInt>& orig2reducedrow) {
  HighsCliqueTable newCliqueTable(ncols);
  newCliqueTable.setinPresolveProbingFlag(inPresolveProbing);
  newCliqueTable.setMinEntriesForParallelism(minEntriesForParallelism);
  for (size_t i = 0; i != cliques.size(); ++i) {
    if (cliques[i].start == -1) continue;
    const HighsInt numActive = cliques[i].numActive();

    for (HighsInt k = cliques[i].start; k != cliques[i].end; ++k) {
      HighsInt col = orig2reducedcol[cliqueentries[k].col];

      if (colDeleted[cliqueentries[k].col] || col == -1 ||
          !globaldomain.isBinary(col) ||
          !postSolveStack.isColLinearlyTransformable(col))
        cliqueentries[k].col = kHighsIInf;
      else
        cliqueentries[k].col = col;
    }

    auto newend =
        std::remove_if(cliqueentries.begin() + cliques[i].start,
                       cliqueentries.begin() + cliques[i].end,
                       [](CliqueVar v) { return v.col == kHighsIInf; });
    HighsInt numvars = static_cast<HighsInt>(
        newend - (cliqueentries.begin() + cliques[i].start));
    // Preserve equality if all active entries survive, else relax to <=.
    if (numvars <= 1) continue;

    HighsInt origin = cliques[i].origin != kHighsIInf ? -1 : kHighsIInf;
    newCliqueTable.doAddClique(&cliqueentries[cliques[i].start], numvars,
                               numvars == numActive && cliques[i].equality,
                               origin);
  }

  newCliqueTable.setAllowParallel(allowParallel);
  *this = std::move(newCliqueTable);
}

void HighsCliqueTable::buildFrom(const HighsLp* origModel,
                                 const HighsCliqueTable& init) {
  assert(init.colsubstituted.size() == colsubstituted.size());
  HighsInt ncols = init.colsubstituted.size();
  HighsCliqueTable newCliqueTable(ncols);
  newCliqueTable.setinPresolveProbingFlag(inPresolveProbing);
  newCliqueTable.setMinEntriesForParallelism(minEntriesForParallelism);
  HighsInt ncliques = init.cliques.size();
  std::vector<CliqueVar> clqBuffer;
  clqBuffer.reserve(2 * static_cast<size_t>(origModel->num_col_));
  for (HighsInt i = 0; i != ncliques; ++i) {
    if (init.cliques[i].start == -1) continue;
    if (init.cliques[i].numActive() <= 1) continue;

    clqBuffer.assign(init.cliqueentries.begin() + init.cliques[i].start,
                     init.cliqueentries.begin() + init.cliques[i].end);
    clqBuffer.erase(
        std::remove_if(clqBuffer.begin(), clqBuffer.end(),
                       [origModel, &init](CliqueVar v) {
                         return origModel->col_lower_[v.col] != 0.0 ||
                                origModel->col_upper_[v.col] != 1.0 ||
                                init.colDeleted[v.col];
                       }),
        clqBuffer.end());
    if (clqBuffer.size() <= 1) continue;

    HighsInt origin = init.cliques[i].origin != kHighsIInf ? -1 : kHighsIInf;
    newCliqueTable.doAddClique(clqBuffer.data(),
                               static_cast<HighsInt>(clqBuffer.size()), false,
                               origin);
  }

  newCliqueTable.colsubstituted = init.colsubstituted;
  newCliqueTable.substitutions = init.substitutions;
  // Currently assume buildFrom is always used for sub-mips
  newCliqueTable.setAllowParallel(false);
  *this = std::move(newCliqueTable);
}
