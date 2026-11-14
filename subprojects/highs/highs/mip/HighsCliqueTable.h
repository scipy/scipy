/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_CLIQUE_TABLE_H_
#define HIGHS_CLIQUE_TABLE_H_

#include <cstdint>
#include <set>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsHash.h"
#include "util/HighsHashTree.h"
#include "util/HighsRandom.h"

class HighsCutPool;
class HighsDomain;
class HighsMipSolver;
class HighsLp;

namespace presolve {
class HighsPostsolveStack;
}

class HPresolveCliqueTable;

class HighsCliqueTable {
  friend class HPresolveCliqueTable;

 public:
  struct CliqueVar {
#ifdef HIGHSINT64
    HighsUInt col : 63;
    HighsUInt val : 1;
#else
    HighsUInt col : 31;
    HighsUInt val : 1;
#endif

    HighsInt index() const { return 2 * col + val; }

    double weight(const std::vector<double>& sol) const {
      return val ? sol[col] : 1.0 - sol[col];
    }

    CliqueVar complement() const { return CliqueVar(col, 1 - val); }

    bool operator==(const CliqueVar& other) const {
      return index() == other.index();
    }

    CliqueVar(HighsInt col, HighsInt val) : col(col), val(val) {}
    CliqueVar() = default;
  };
  struct Clique {
    HighsInt start;
    HighsInt end;
    HighsInt origin;
    HighsInt numZeroFixed;
    bool equality;
    HighsInt numActive() const { return end - start - numZeroFixed; }
  };

  struct Substitution {
    HighsInt substcol;
    CliqueVar replace;
  };

 private:
  std::vector<CliqueVar> cliqueentries;

  std::vector<HighsHashTree<HighsInt, HighsInt>> invertedHashList;
  std::vector<HighsHashTree<HighsInt>> invertedHashListSizeTwo;
  HighsHashTable<std::pair<CliqueVar, CliqueVar>, HighsInt> sizeTwoCliques;
  // the clique table may contain several size-two cliques on the same edge;
  // sizeTwoCliques stores one of them, and all of them are kept in a doubly
  // linked list (indexed by clique id, -1 marks the ends)
  struct SizeTwoLinks {
    HighsInt prev;
    HighsInt next;
  };
  std::vector<SizeTwoLinks> sizeTwoLinks;

  std::set<std::pair<HighsInt, HighsInt>> freespaces;
  std::vector<HighsInt> freeslots;
  std::vector<Clique> cliques;
  std::vector<HighsInt> numcliquesvar;
  std::vector<CliqueVar> infeasvertexstack;

  std::vector<HighsInt> colsubstituted;
  std::vector<Substitution> substitutions;
  std::vector<HighsInt> deletedrows;
  std::vector<std::pair<HighsInt, CliqueVar>> cliqueextensions;
  std::vector<HighsBool> iscandidate;
  std::vector<HighsBool> colDeleted;
  std::vector<uint32_t> cliquehits;
  std::vector<HighsInt> cliquehitinds;

  // HighsHashTable<std::pair<CliqueVar, CliqueVar>> invertedEdgeCache;

  HighsRandom randgen;
  HighsInt nfixings;
  HighsInt numEntries;
  HighsInt maxEntries;
  HighsInt minEntriesForParallelism;
  bool inPresolveProbing;
  bool allowParallel;

  void unlink(HighsInt pos, HighsInt cliqueid);

  void link(HighsInt pos, HighsInt cliqueid);

  HighsInt findCommonCliqueId(int64_t& numQueries, CliqueVar v1,
                              CliqueVar v2) const;

  HighsInt findCommonCliqueId(CliqueVar v1, CliqueVar v2) {
    return findCommonCliqueId(numNeighbourhoodQueries, v1, v2);
  }

  HighsInt runCliqueSubsumption(const HighsDomain& globaldom,
                                std::vector<CliqueVar>& clique);

  void cliqueSubsumption(
      const std::vector<CliqueVar>& clique, bool& redundant,
      HighsInt& dominatingOrigin,
      const std::function<void(HighsInt)>& removeCliqueCallback);

  void collectCliques(const std::vector<CliqueVar>& clique);

  struct BronKerboschData {
    const std::vector<double>& sol;
    std::vector<CliqueVar> P;
    std::vector<CliqueVar> R;
    std::vector<CliqueVar> Z;
    std::vector<std::vector<CliqueVar>> cliques;
    std::vector<HighsInt> neighbourhoodInds;

    double wR = 0.0;
    double minW = 1.05;
    double feastol = 1e-6;
    HighsInt ncalls = 0;
    HighsInt maxcalls = 10000;
    HighsInt maxcliques = 100;
    int64_t maxNeighbourhoodQueries = std::numeric_limits<int64_t>::max();
    int64_t numNeighbourhoodQueries = 0;

    bool stop() const {
      return maxcalls == ncalls || int(cliques.size()) == maxcliques ||
             numNeighbourhoodQueries > maxNeighbourhoodQueries;
    }

    BronKerboschData(const std::vector<double>& sol) : sol(sol) {}
  };

  void bronKerboschRecurse(BronKerboschData& data, HighsInt Plen,
                           const CliqueVar* X, HighsInt Xlen) const;

  void extractCliques(const HighsMipSolver& mipsolver,
                      std::vector<HighsInt>& inds, std::vector<double>& vals,
                      std::vector<int8_t>& complementation, double rhs,
                      HighsInt nbin, std::vector<HighsInt>& perm,
                      std::vector<CliqueVar>& clique, double feastol);

  void processInfeasibleVertices(HighsDomain& domain);

  void propagateAndCleanup(HighsDomain& globaldom);

  // add / remove a size-two clique to / from the cliques on its edge; do
  // nothing for other cliques
  void addSizeTwoClique(HighsInt cliqueid);

  void removeSizeTwoClique(HighsInt cliqueid);

  void queryNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                          int64_t& numNeighbourhoodqueries, CliqueVar v,
                          CliqueVar* q, HighsInt N) const;

  void recordSubstitution(Substitution substitution);

  void replaceLiteral(CliqueVar substitutedVar, CliqueVar replacementVar);

  template <bool Sort = true, typename Comparator>
  HighsInt extendClique(std::vector<HighsInt>& neighbourhoodInds,
                        std::vector<CliqueVar>& clqVars, HighsInt seedPos,
                        HighsInt extensionEnd,
                        const Comparator& candidateOrder);

  HighsInt extendClique(std::vector<HighsInt>& neighbourhoodInds,
                        std::vector<CliqueVar>& clqVars, HighsInt seedPos,
                        HighsInt extensionEnd);

 public:
  int64_t numNeighbourhoodQueries;

  HighsCliqueTable(HighsInt ncols) {
    resize(static_cast<size_t>(ncols));
    nfixings = 0;
    numNeighbourhoodQueries = 0;
    numEntries = 0;
    maxEntries = kHighsIInf;
    minEntriesForParallelism = kHighsIInf;
    inPresolveProbing = false;
    allowParallel = true;
  }

  void resize(size_t ncols) {
    invertedHashList.resize(2 * ncols);
    invertedHashListSizeTwo.resize(2 * ncols);
    numcliquesvar.resize(2 * ncols, 0);
    colsubstituted.resize(ncols);
    colDeleted.resize(ncols, false);
  }

  void setinPresolveProbingFlag(const bool inPresolveProbing) {
    this->inPresolveProbing = inPresolveProbing;
  }

  HighsInt getNumEntries() const { return numEntries; }

  HighsRandom& getRandgen() { return randgen; }

  int64_t& getNumNeighbourhoodQueries() { return numNeighbourhoodQueries; }

  HighsInt partitionNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                                  int64_t& numNeighbourhoodqueries, CliqueVar v,
                                  CliqueVar* q, HighsInt N) const;

  HighsInt shrinkToNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                                 int64_t& numNeighbourhoodqueries, CliqueVar v,
                                 CliqueVar* q, HighsInt N);

  bool fixCol(HighsDomain& globaldom, CliqueVar v,
              bool doProcessInfeasibleVertices = false);

  bool processNewEdge(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void doAddClique(const CliqueVar* cliquevars, HighsInt numcliquevars,
                   bool equality = false, HighsInt origin = kHighsIInf);

  void addClique(const HighsMipSolver& mipsolver, CliqueVar* cliquevars,
                 HighsInt numcliquevars, bool equality = false,
                 HighsInt origin = kHighsIInf);

  void removeClique(HighsInt cliqueid, bool recordDeletedRow = true);

  void fixLastActiveAndRemove(HighsDomain& globaldom, HighsInt cliqueid);

  void resolveSubstitution(CliqueVar& v) const;

  void resolveSubstitution(HighsInt& col, double& val, double& rhs) const;

  std::vector<HighsInt>& getDeletedRows() { return deletedrows; }

  const std::vector<HighsInt>& getDeletedRows() const { return deletedrows; }

  std::vector<Substitution>& getSubstitutions() { return substitutions; }

  const std::vector<Substitution>& getSubstitutions() const {
    return substitutions;
  }

  const Substitution* getSubstitution(HighsInt col) const {
    return colsubstituted[col] ? &substitutions[colsubstituted[col] - 1]
                               : nullptr;
  }

  std::vector<std::pair<HighsInt, CliqueVar>>& getCliqueExtensions() {
    return cliqueextensions;
  }

  const std::vector<std::pair<HighsInt, CliqueVar>>& getCliqueExtensions()
      const {
    return cliqueextensions;
  }

  void setMaxEntries(HighsInt numNz) {
    this->maxEntries = 2000000 + 10 * numNz;
  }

  void setMinEntriesForParallelism(HighsInt minEntriesForParallelism) {
    this->minEntriesForParallelism = minEntriesForParallelism;
  }

  bool isFull() const { return numEntries >= maxEntries; }

  HighsInt getNumFixings() const { return nfixings; }

  void cliquePartition(std::vector<CliqueVar>& clqVars,
                       std::vector<HighsInt>& partitionStart);

  void cliquePartition(const std::vector<double>& objective,
                       std::vector<CliqueVar>& clqVars,
                       std::vector<HighsInt>& partitionStart);

  void cliqueCover(std::vector<CliqueVar>& clqVars,
                   std::vector<std::vector<CliqueVar>>& cover);

  bool foundCover(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void extractCliques(HighsMipSolver& mipsolver, bool transformRows = true);

  void extractCliquesFromCut(const HighsMipSolver& mipsolver,
                             const HighsInt* inds, const double* vals,
                             HighsInt len, double rhs);

  void extractObjCliques(HighsMipSolver& mipsolver);

  void vertexInfeasible(HighsDomain& globaldom, HighsInt col, HighsInt val);

  bool haveCommonClique(CliqueVar v1, CliqueVar v2) {
    if (v1.col == v2.col) return false;
    return findCommonCliqueId(v1, v2) != -1;
  }

  bool haveCommonClique(int64_t& numQueries, CliqueVar v1, CliqueVar v2) const {
    if (v1.col == v2.col) return false;
    return findCommonCliqueId(numQueries, v1, v2) != -1;
  }

  std::pair<const CliqueVar*, HighsInt> findCommonClique(CliqueVar v1,
                                                         CliqueVar v2) {
    std::pair<const CliqueVar*, HighsInt> c{nullptr, 0};
    if (v1 == v2) return c;
    HighsInt clq = findCommonCliqueId(v1, v2);
    if (clq == -1) return c;

    c.first = &cliqueentries[cliques[clq].start];
    c.second = cliques[clq].end - cliques[clq].start;
    return c;
  }

  void separateCliques(const HighsMipSolver& mipsolver,
                       const std::vector<double>& sol, HighsCutPool& cutpool,
                       double feastol, HighsRandom& randgen,
                       int64_t& localNumNeighbourhoodQueries);

  std::vector<std::vector<CliqueVar>> separateCliques(
      const std::vector<double>& sol, const HighsDomain& globaldom,
      double feastol);

  std::vector<std::vector<CliqueVar>> computeMaximalCliques(
      const std::vector<CliqueVar>& vars, double feastol);

  void cleanupFixed(HighsDomain& globaldom);

  void addImplications(HighsDomain& domain, HighsInt col, HighsInt val);

  HighsInt getNumImplications(HighsInt col) const;

  HighsInt getNumImplications(HighsInt col, bool val) const;

  void runCliqueMerging(HighsDomain& globaldomain);

  void runCliqueMerging(HighsDomain& globaldomain,
                        std::vector<CliqueVar>& clique, bool equation = false);

  void rebuild(HighsInt ncols,
               const presolve::HighsPostsolveStack& postSolveStack,
               const HighsDomain& globaldomain,
               const std::vector<HighsInt>& cIndex,
               const std::vector<HighsInt>& rIndex);

  void buildFrom(const HighsLp* origModel, const HighsCliqueTable& init);

  HighsInt numCliques() const { return cliques.size() - freeslots.size(); }

  HighsInt numCliques(CliqueVar v) const { return numcliquesvar[v.index()]; }

  HighsInt numCliques(HighsInt col, bool val) const {
    return numcliquesvar[CliqueVar(col, val).index()];
  }

  void setAllowParallel(const bool allowParallel) {
    this->allowParallel = allowParallel;
  }
};

#endif
