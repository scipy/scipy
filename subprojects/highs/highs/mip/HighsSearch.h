/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SEARCH_H_
#define HIGHS_SEARCH_H_

#include <cstdint>
#include <queue>
#include <vector>

#include "mip/HighsConflictPool.h"
#include "mip/HighsDebugSol.h"
#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipWorker.h"
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSeparation.h"
#include "presolve/HighsSymmetry.h"
#include "util/HighsHash.h"

class HighsMipSolver;
class HighsMipWorker;
class HighsImplications;
class HighsCliqueTable;

class HighsSearch {
 public:
  HighsMipWorker& mipworker;

  const HighsMipSolver& mipsolver;
  HighsLpRelaxation* lp;
  HighsDomain localdom;
  HighsPseudocost& pseudocost;
  HighsRandom random;
  int64_t nnodes;
  int64_t nleaves;
  int64_t lpiterations;
  int64_t heurlpiterations;
  int64_t sblpiterations;
  double upper_limit;
  HighsCDouble treeweight;
  std::vector<HighsInt> inds;
  std::vector<double> vals;
  HighsInt depthoffset;
  bool inbranching;
  bool inheuristic;
  bool countTreeWeight;

 public:
  enum class ChildSelectionRule {
    kUp,
    kDown,
    kRootSol,
    kObj,
    kRandom,
    kBestCost,
    kWorstCost,
    kDisjunction,
    kHybridInferenceCost,
  };

  enum class NodeResult {
    kBoundExceeding,
    kDomainInfeasible,
    kLpInfeasible,
    kBranched,
    kSubOptimal,
    kOpen,
  };

 private:
  ChildSelectionRule childselrule;

  struct NodeData {
    double lower_bound;
    double estimate;
    double branching_point;
    // we store the lp objective separately to the lower bound since the lower
    // bound could be above the LP objective when cuts age out or below when the
    // LP is unscaled dual infeasible and it is not set. We still want to use
    // the objective for pseudocost updates and tiebreaking of best bound node
    // selection
    double lp_objective;
    double other_child_lb;
    std::shared_ptr<const HighsBasis> nodeBasis;
    std::shared_ptr<const StabilizerOrbits> stabilizerOrbits;
    HighsDomainChange branchingdecision;
    HighsInt domgchgStackPos;
    uint8_t skipDepthCount;
    uint8_t opensubtrees;

    NodeData(double parentlb = -kHighsInf, double parentestimate = -kHighsInf,
             std::shared_ptr<const HighsBasis> parentBasis = nullptr,
             std::shared_ptr<const StabilizerOrbits> stabilizerOrbits = nullptr)
        : lower_bound(parentlb),
          estimate(parentestimate),
          branching_point(0.0),
          lp_objective(-kHighsInf),
          other_child_lb(parentlb),
          nodeBasis(std::move(parentBasis)),
          stabilizerOrbits(std::move(stabilizerOrbits)),
          branchingdecision{0.0, -1, HighsBoundType::kLower},
          domgchgStackPos(-1),
          skipDepthCount(0),
          opensubtrees(2) {}
  };

  enum ReliableFlags {
    kUpReliable = 1,
    kDownReliable = 2,
    kReliable = kDownReliable | kUpReliable,
  };

  std::vector<double> subrootsol;
  std::vector<NodeData> nodestack;
  StabilizerOrbitWorkspace stabilizerOrbitWorkspace;
  HighsHashTable<HighsInt, int> reliableatnode;

  int branchingVarReliableAtNodeFlags(HighsInt col) const {
    auto it = reliableatnode.find(col);
    if (it == nullptr) return 0;
    return *it;
  }

  bool branchingVarReliableAtNode(HighsInt col) const {
    auto it = reliableatnode.find(col);
    if (it == nullptr || *it != kReliable) return false;

    return true;
  }

  void markBranchingVarUpReliableAtNode(HighsInt col) {
    reliableatnode[col] |= kUpReliable;
  }

  void markBranchingVarDownReliableAtNode(HighsInt col) {
    reliableatnode[col] |= kDownReliable;
  }

  bool orbitsValidInChildNode(const HighsDomainChange& branchChg) const;

  void openChildNode(HighsInt col, double boundval, double branchpoint,
                     HighsBoundType boundtype);

  void stashNodeToProcessed(double lb, double estimate, HighsInt depth);

 public:
  HighsSearch(HighsMipWorker& mipworker, HighsPseudocost& pseudocost);

  void setRINSNeighbourhood(const std::vector<double>& basesol,
                            const std::vector<double>& relaxsol);

  void setRENSNeighbourhood(const std::vector<double>& lpsol);

  double getCutoffBound() const;

  void setLpRelaxation(HighsLpRelaxation* lp) { this->lp = lp; }

  double checkSol(const std::vector<double>& sol, bool& integerfeasible) const;

  void createNewNode();

  void cutoffNode();

  void branchDownwards(HighsInt col, double newub, double branchpoint);

  void branchUpwards(HighsInt col, double newlb, double branchpoint);

  void flipBranchingDecision(NodeData& node);

  void setMinReliable(HighsInt minreliable);

  void setHeuristic(bool inheuristic) {
    this->inheuristic = inheuristic;
    if (inheuristic) childselrule = ChildSelectionRule::kHybridInferenceCost;
  }

  void addBoundExceedingConflict();

  void resetLocalDomain();

  int64_t getHeuristicLpIterations() const;

  int64_t getTotalLpIterations() const;

  int64_t getLocalLpIterations() const;

  int64_t& getLocalNodes();

  int64_t& getLocalLeaves();

  int64_t getStrongBranchingLpIterations() const;

  bool hasNode() const { return !nodestack.empty(); }

  bool currentNodePruned() const { return nodestack.back().opensubtrees == 0; }

  double getCurrentEstimate() const { return nodestack.back().estimate; }

  double getCurrentLowerBound() const { return nodestack.back().lower_bound; }

  HighsInt getCurrentDepth() const { return nodestack.size() + depthoffset; }

  void stashOpenNodes();

  void stashCurrentNode();

  void flushStatistics(HighsMipSolver& mipsolver);

  void resetStatistics();

  void installNode(HighsNodeQueue::OpenNode&& node);

  void addInfeasibleConflict();

  HighsInt selectBranchingCandidate(int64_t maxSbIters, double& downNodeLb,
                                    double& upNodeLb);

  void evalUnreliableBranchCands();

  const NodeData* getParentNodeData() const;

  NodeResult evaluateNode();

  NodeResult branch();

  /// backtrack one level in DFS manner
  bool backtrack(bool recoverBasis = true);

  /// backtrack until the next node (DFS manner) that can be worthwhile
  /// to continue the plunge. Put unpromising nodes to workers processedNodes
  bool backtrackPlunge();

  /// for heuristics. Will discard nodes above targetDepth regardless of their
  /// status
  bool backtrackUntilDepth(HighsInt targetDepth);

  void printDisplayLine(char first, bool header = false);

  NodeResult dive(int64_t nodeLim = std::numeric_limits<int64_t>::max());

  HighsDomain& getLocalDomain() { return localdom; }

  const HighsDomain& getLocalDomain() const { return localdom; }

  HighsPseudocost& getPseudoCost() { return pseudocost; }

  const HighsPseudocost& getPseudoCost() const { return pseudocost; }

  void solveDepthFirst(int64_t maxbacktracks = 1,
                       int64_t maxnodes = std::numeric_limits<int64_t>::max());

  double getFeasTol() const;
  double getUpperLimit() const;
  double getEpsilon() const;

  const std::vector<double>& getRootLpSol() const;
  const std::vector<HighsInt>& getIntegralCols() const;

  HighsDomain& getDomain() const;
  HighsConflictPool& getConflictPool() const;
  HighsCutPool& getCutPool() const;

  const HighsNodeQueue& getNodeQueue() const;

  bool checkLimits(int64_t nodeOffset = 0) const;

  bool checkLocalLimits() const;

  HighsSymmetries& getSymmetries() const;

  bool addIncumbent(const std::vector<double>& sol, double solobj,
                    const int solution_source,
                    const bool print_display_line = true);
};

#endif
