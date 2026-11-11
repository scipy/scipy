/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsNodeQueue.h"

#include <algorithm>
#include <tuple>

#include "lp_data/HConst.h"
#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsSplay.h"

#define ESTIMATE_WEIGHT .5
#define LOWERBOUND_WEIGHT .5

namespace highs {
template <>
struct RbTreeTraits<HighsNodeQueue::NodeLowerRbTree> {
  using KeyType = std::tuple<double, HighsInt, double, int64_t>;
  using LinkType = int64_t;
};

template <>
struct RbTreeTraits<HighsNodeQueue::NodeHybridEstimRbTree> {
  using KeyType = std::tuple<double, HighsInt, int64_t>;
  using LinkType = int64_t;
};

template <>
struct RbTreeTraits<HighsNodeQueue::SuboptimalNodeRbTree> {
  using KeyType = std::pair<double, int64_t>;
  using LinkType = int64_t;
};
}  // namespace highs

using namespace highs;

class HighsNodeQueue::NodeLowerRbTree : public CacheMinRbTree<NodeLowerRbTree> {
  HighsNodeQueue* nodeQueue;

 public:
  NodeLowerRbTree(HighsNodeQueue* nodeQueue)
      : CacheMinRbTree<NodeLowerRbTree>(nodeQueue->lowerRoot,
                                        nodeQueue->lowerMin),
        nodeQueue(nodeQueue) {}

  RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) {
    return nodeQueue->nodes[node].lowerLinks;
  }
  const RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) const {
    return nodeQueue->nodes[node].lowerLinks;
  }
  std::tuple<double, HighsInt, double, int64_t> getKey(int64_t node) const {
    return std::make_tuple(
        nodeQueue->nodes[node].lower_bound,
        static_cast<HighsInt>(nodeQueue->nodes[node].domchgstack.size()),
        nodeQueue->nodes[node].estimate, node);
  }
};

class HighsNodeQueue::NodeHybridEstimRbTree
    : public CacheMinRbTree<NodeHybridEstimRbTree> {
  HighsNodeQueue* nodeQueue;

 public:
  NodeHybridEstimRbTree(HighsNodeQueue* nodeQueue)
      : CacheMinRbTree<NodeHybridEstimRbTree>(nodeQueue->hybridEstimRoot,
                                              nodeQueue->hybridEstimMin),
        nodeQueue(nodeQueue) {}

  RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) {
    return nodeQueue->nodes[node].hybridEstimLinks;
  }
  const RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) const {
    return nodeQueue->nodes[node].hybridEstimLinks;
  }
  std::tuple<double, HighsInt, int64_t> getKey(int64_t node) const {
    constexpr double kLbWeight = 0.5;
    constexpr double kEstimWeight = 0.5;
    return std::make_tuple(kLbWeight * nodeQueue->nodes[node].lower_bound +
                               kEstimWeight * nodeQueue->nodes[node].estimate,
                           -HighsInt(nodeQueue->nodes[node].domchgstack.size()),
                           node);
  }
};

class HighsNodeQueue::SuboptimalNodeRbTree
    : public CacheMinRbTree<SuboptimalNodeRbTree> {
  HighsNodeQueue* nodeQueue;

 public:
  SuboptimalNodeRbTree(HighsNodeQueue* nodeQueue)
      : CacheMinRbTree<SuboptimalNodeRbTree>(nodeQueue->suboptimalRoot,
                                             nodeQueue->suboptimalMin),
        nodeQueue(nodeQueue) {}

  RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) {
    return nodeQueue->nodes[node].lowerLinks;
  }
  const RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) const {
    return nodeQueue->nodes[node].lowerLinks;
  }

  std::pair<double, int64_t> getKey(int64_t node) const {
    return std::make_pair(nodeQueue->nodes[node].lower_bound, node);
  }
};

void HighsNodeQueue::link_estim(int64_t node) {
  assert(node != -1);
  NodeHybridEstimRbTree rbTree(this);
  rbTree.link(node);
}

void HighsNodeQueue::unlink_estim(int64_t node) {
  assert(node != -1);
  NodeHybridEstimRbTree rbTree(this);
  rbTree.unlink(node);
}

void HighsNodeQueue::link_lower(int64_t node) {
  assert(node != -1);
  NodeLowerRbTree rbTree(this);
  rbTree.link(node);
}

void HighsNodeQueue::unlink_lower(int64_t node) {
  assert(node != -1);
  NodeLowerRbTree rbTree(this);
  rbTree.unlink(node);
}

void HighsNodeQueue::link_suboptimal(int64_t node) {
  assert(node != -1);
  SuboptimalNodeRbTree rbTree(this);
  rbTree.link(node);
  ++numSuboptimal;
}

void HighsNodeQueue::unlink_suboptimal(int64_t node) {
  assert(node != -1);
  SuboptimalNodeRbTree rbTree(this);
  rbTree.unlink(node);
  --numSuboptimal;
}

void HighsNodeQueue::link_domchgs(int64_t node) {
  assert(node != -1);
  HighsInt numchgs = nodes[node].domchgstack.size();
  nodes[node].domchglinks.resize(numchgs);

  for (HighsInt i = 0; i != numchgs; ++i) {
    double val = nodes[node].domchgstack[i].boundval;
    HighsInt col = nodes[node].domchgstack[i].column;
    switch (nodes[node].domchgstack[i].boundtype) {
      case HighsBoundType::kLower:
        nodes[node].domchglinks[i] =
            colLowerNodesPtr.get()[col].emplace(val, node).first;
        break;
      case HighsBoundType::kUpper:
        nodes[node].domchglinks[i] =
            colUpperNodesPtr.get()[col].emplace(val, node).first;
    }
  }
}

void HighsNodeQueue::unlink_domchgs(int64_t node) {
  assert(node != -1);
  HighsInt numchgs = nodes[node].domchgstack.size();

  for (HighsInt i = 0; i != numchgs; ++i) {
    HighsInt col = nodes[node].domchgstack[i].column;
    switch (nodes[node].domchgstack[i].boundtype) {
      case HighsBoundType::kLower:
        colLowerNodesPtr.get()[col].erase(nodes[node].domchglinks[i]);
        break;
      case HighsBoundType::kUpper:
        colUpperNodesPtr.get()[col].erase(nodes[node].domchglinks[i]);
    }
  }

  nodes[node].domchglinks.clear();
  nodes[node].domchglinks.shrink_to_fit();
}

double HighsNodeQueue::link(int64_t node) {
  if (nodes[node].lower_bound > optimality_limit) {
    assert(nodes[node].estimate != kHighsInf);
    nodes[node].estimate = kHighsInf;
    link_suboptimal(node);
    link_domchgs(node);
    return std::ldexp(1.0, 1 - nodes[node].depth);
  }

  link_estim(node);
  link_lower(node);
  link_domchgs(node);
  return 0.0;
}

void HighsNodeQueue::unlink(int64_t node) {
  if (nodes[node].estimate == kHighsInf) {
    unlink_suboptimal(node);
  } else {
    unlink_estim(node);
    unlink_lower(node);
  }
  unlink_domchgs(node);
  freeslots.push(node);
}

void HighsNodeQueue::setNumCol(HighsInt numCol) {
  if (this->numCol == numCol) return;
  this->numCol = numCol;
  allocatorState = std::unique_ptr<AllocatorState>(new AllocatorState());

  if (numCol == 0) return;
  colLowerNodesPtr =
      NodeSetArray((NodeSet*)::operator new(sizeof(NodeSet) * numCol));
  colUpperNodesPtr =
      NodeSetArray((NodeSet*)::operator new(sizeof(NodeSet) * numCol));

  NodesetAllocator<std::pair<double, int64_t>> allocator(allocatorState.get());
  for (HighsInt i = 0; i < numCol; ++i) {
    new (colLowerNodesPtr.get() + i) NodeSet(allocator);
    new (colUpperNodesPtr.get() + i) NodeSet(allocator);
  }
}

void HighsNodeQueue::checkGlobalBounds(HighsInt col, double lb, double ub,
                                       double feastol,
                                       HighsCDouble& treeweight) {
  std::set<int64_t> delnodes;

  auto colLowerNodes = colLowerNodesPtr.get();
  auto colUpperNodes = colUpperNodesPtr.get();

  auto prunestart =
      colLowerNodes[col].lower_bound(std::make_pair(ub + feastol, -1));
  for (auto it = prunestart; it != colLowerNodes[col].end(); ++it)
    delnodes.insert(it->second);

  auto pruneend =
      colUpperNodes[col].upper_bound(std::make_pair(lb - feastol, kHighsIInf));
  for (auto it = colUpperNodes[col].begin(); it != pruneend; ++it)
    delnodes.insert(it->second);

  for (int64_t delnode : delnodes) {
    if (nodes[delnode].estimate != kHighsInf)
      treeweight += std::ldexp(1.0, 1 - nodes[delnode].depth);
    unlink(delnode);
  }
}

double HighsNodeQueue::pruneInfeasibleNodes(HighsDomain& globaldomain,
                                            double feastol) {
  size_t numchgs;

  HighsCDouble treeweight = 0.0;

  do {
    if (globaldomain.infeasible()) break;

    numchgs = globaldomain.getDomainChangeStack().size();

    assert(numCol == (HighsInt)globaldomain.col_lower_.size());

    for (HighsInt i = 0; i < numCol; ++i) {
      checkGlobalBounds(i, globaldomain.col_lower_[i],
                        globaldomain.col_upper_[i], feastol, treeweight);
    }

    size_t numopennodes = numNodes();
    if (numopennodes == 0) break;

    auto colLowerNodes = colLowerNodesPtr.get();
    auto colUpperNodes = colUpperNodesPtr.get();

    for (HighsInt i = 0; i < numCol; ++i) {
      if (colLowerNodes[i].size() == numopennodes) {
        double globallb = colLowerNodes[i].begin()->first;
        if (globallb > globaldomain.col_lower_[i]) {
          globaldomain.changeBound(HighsBoundType::kLower, i, globallb,
                                   HighsDomain::Reason::unspecified());
          if (globaldomain.infeasible()) break;
        }
      }

      if (colUpperNodes[i].size() == numopennodes) {
        double globalub = colUpperNodes[i].rbegin()->first;
        if (globalub < globaldomain.col_upper_[i]) {
          globaldomain.changeBound(HighsBoundType::kUpper, i, globalub,
                                   HighsDomain::Reason::unspecified());
          if (globaldomain.infeasible()) break;
        }
      }
    }

    globaldomain.propagate();
  } while (numchgs != globaldomain.getDomainChangeStack().size());

  return double(treeweight);
}

double HighsNodeQueue::pruneNode(int64_t nodeId) {
  double treeweight = nodes[nodeId].estimate != kHighsInf
                          ? std::ldexp(1.0, 1 - nodes[nodeId].depth)
                          : 0.0;
  unlink(nodeId);
  return treeweight;
}

double HighsNodeQueue::performBounding(double upper_limit) {
  NodeLowerRbTree lowerTree(this);

  if (lowerTree.empty()) return 0.0;

  HighsCDouble treeweight = 0.0;

  int64_t maxLbNode = lowerTree.last();
  while (maxLbNode != -1) {
    if (nodes[maxLbNode].lower_bound < upper_limit) break;
    int64_t next = lowerTree.predecessor(maxLbNode);
    treeweight += pruneNode(maxLbNode);
    maxLbNode = next;
  }

  if (optimality_limit < upper_limit) {
    while (maxLbNode != -1) {
      if (nodes[maxLbNode].lower_bound < optimality_limit) break;
      int64_t next = lowerTree.predecessor(maxLbNode);
      assert(nodes[maxLbNode].estimate != kHighsInf);
      unlink_estim(maxLbNode);
      unlink_lower(maxLbNode);
      treeweight += std::ldexp(1.0, 1 - nodes[maxLbNode].depth);
      nodes[maxLbNode].estimate = kHighsInf;
      link_suboptimal(maxLbNode);
      maxLbNode = next;
    }
  }

  if (numSuboptimal) {
    SuboptimalNodeRbTree suboptimalTree(this);
    maxLbNode = suboptimalTree.last();
    while (maxLbNode != -1) {
      if (nodes[maxLbNode].lower_bound < upper_limit) break;
      int64_t next = suboptimalTree.predecessor(maxLbNode);
      unlink(maxLbNode);
      maxLbNode = next;
    }
  }

  return double(treeweight);
}

double HighsNodeQueue::emplaceNode(std::vector<HighsDomainChange>&& domchgs,
                                   std::vector<HighsInt>&& branchPositions,
                                   double lower_bound, double estimate,
                                   HighsInt depth) {
  int64_t pos;

  assert(estimate != kHighsInf);

  if (freeslots.empty()) {
    pos = nodes.size();
    nodes.emplace_back(std::move(domchgs), std::move(branchPositions),
                       lower_bound, estimate, depth);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    nodes[pos] = OpenNode(std::move(domchgs), std::move(branchPositions),
                          lower_bound, estimate, depth);
  }

  assert(nodes[pos].lower_bound == lower_bound);
  assert(nodes[pos].estimate == estimate);
  assert(nodes[pos].depth == depth);

  return link(pos);
}

HighsNodeQueue::OpenNode&& HighsNodeQueue::popBestNode() {
  int64_t bestNode = hybridEstimMin;

  unlink(bestNode);

  return std::move(nodes[bestNode]);
}

HighsNodeQueue::OpenNode&& HighsNodeQueue::popBestBoundNode() {
  int64_t bestBoundNode = lowerMin;

  unlink(bestBoundNode);

  return std::move(nodes[bestBoundNode]);
}

double HighsNodeQueue::getBestLowerBound() const {
  double lb = lowerMin == -1 ? kHighsInf : nodes[lowerMin].lower_bound;

  if (suboptimalMin == -1) return lb;

  return std::min(nodes[suboptimalMin].lower_bound, lb);
}

HighsInt HighsNodeQueue::getBestBoundDomchgStackSize() const {
  HighsInt domchgStackSize = lowerMin == -1
                                 ? kHighsIInf
                                 : HighsInt(nodes[lowerMin].domchgstack.size());
  if (suboptimalMin == -1) return domchgStackSize;

  return std::min(HighsInt(nodes[suboptimalMin].domchgstack.size()),
                  domchgStackSize);
}

void HighsNodeQueue::clear() {
  const bool original = false;
  HighsNodeQueue nodequeue;
  nodequeue.setNumCol(numCol);
  if (original) {
    *this = std::move(nodequeue);
  } else {
    (*this).nodes = std::move(nodequeue.nodes);
    (*this).colLowerNodesPtr = std::move(nodequeue.colLowerNodesPtr);
    (*this).colUpperNodesPtr = std::move(nodequeue.colUpperNodesPtr);
    (*this).freeslots = std::move(nodequeue.freeslots);
    (*this).allocatorState = std::move(nodequeue.allocatorState);
    (*this).lowerRoot = nodequeue.lowerRoot;
    (*this).lowerMin = nodequeue.lowerMin;
    (*this).hybridEstimRoot = nodequeue.hybridEstimRoot;
    (*this).hybridEstimMin = nodequeue.hybridEstimMin;
    (*this).suboptimalRoot = nodequeue.suboptimalRoot;
    (*this).suboptimalMin = nodequeue.suboptimalMin;
    (*this).numSuboptimal = nodequeue.numSuboptimal;
    (*this).optimality_limit = nodequeue.optimality_limit;
    (*this).numCol = nodequeue.numCol;
  }
}
