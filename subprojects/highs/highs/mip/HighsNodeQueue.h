/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_NODE_QUEUE_H_
#define HIGHS_NODE_QUEUE_H_

#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <queue>
#include <set>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsDomainChange.h"
#include "util/HighsCDouble.h"
#include "util/HighsRbTree.h"

class HighsDomain;
class HighsLpRelaxation;

class HighsNodeQueue {
 public:
  template <int S>
  struct ChunkWithSize {
    ChunkWithSize* next;
    typename std::aligned_storage<S, alignof(std::max_align_t)>::type storage;
  };

  using Chunk =
      ChunkWithSize<4096 - offsetof(ChunkWithSize<sizeof(void*)>, storage)>;

  struct AllocatorState {
    void* freeListHead = nullptr;
    char* currChunkStart = nullptr;
    char* currChunkEnd = nullptr;
    Chunk* chunkListHead = nullptr;

    AllocatorState() = default;
    AllocatorState(const AllocatorState&) = delete;
    AllocatorState(AllocatorState&& other)
        : freeListHead(other.freeListHead),
          currChunkStart(other.currChunkStart),
          currChunkEnd(other.currChunkEnd),
          chunkListHead(other.chunkListHead) {
      other.chunkListHead = nullptr;
    }

    AllocatorState& operator=(const AllocatorState&) = delete;

    AllocatorState& operator=(AllocatorState&& other) {
      freeListHead = other.freeListHead;
      currChunkStart = other.currChunkStart;
      currChunkEnd = other.currChunkEnd;
      chunkListHead = other.chunkListHead;

      other.chunkListHead = nullptr;
      return *this;
    }

    ~AllocatorState() {
      while (chunkListHead) {
        Chunk* delChunk = chunkListHead;
        chunkListHead = delChunk->next;
        delete delChunk;
      }
    }
  };

  template <typename T>
  class NodesetAllocator {
    union FreelistNode {
      FreelistNode* next;
      typename std::aligned_storage<sizeof(T), alignof(T)>::type storage;
    };

   public:
    AllocatorState* state;
    using value_type = T;
    using size_type = std::size_t;
    using propagate_on_container_move_assignment = std::true_type;

    NodesetAllocator(AllocatorState* state) : state(state) {}
    NodesetAllocator(const NodesetAllocator& other) noexcept = default;
    template <typename U>
    NodesetAllocator(const NodesetAllocator<U>& other) noexcept
        : state(other.state) {}
    NodesetAllocator(NodesetAllocator&& other) noexcept = default;
    NodesetAllocator& operator=(const NodesetAllocator&) noexcept = default;
    NodesetAllocator& operator=(NodesetAllocator&& other) noexcept = default;
    ~NodesetAllocator() noexcept = default;

    T* allocate(size_type n) {
      if (n == 1) {
        T* ptr = reinterpret_cast<T*>(state->freeListHead);
        if (ptr) {
          state->freeListHead =
              reinterpret_cast<FreelistNode*>(state->freeListHead)->next;
        } else {
          ptr = reinterpret_cast<T*>(state->currChunkStart);
          if (!ptr || state->currChunkStart + sizeof(FreelistNode) >
                          state->currChunkEnd) {
            auto newChunk = new Chunk;
            newChunk->next = state->chunkListHead;
            state->chunkListHead = newChunk;
            state->currChunkStart = reinterpret_cast<char*>(&newChunk->storage);
            state->currChunkEnd =
                state->currChunkStart + sizeof(newChunk->storage);
            ptr = reinterpret_cast<T*>(state->currChunkStart);
            state->currChunkStart += sizeof(FreelistNode);
          } else {
            state->currChunkStart += sizeof(FreelistNode);
          }
        }
        return ptr;
      }

      return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* ptr, size_type n) noexcept {
      if (n == 1) {
        FreelistNode* node = reinterpret_cast<FreelistNode*>(ptr);
        node->next = reinterpret_cast<FreelistNode*>(state->freeListHead);
        state->freeListHead = node;
      } else {
        ::operator delete(ptr);
      }
    }
  };

  using NodeSet = std::set<std::pair<double, int64_t>,
                           std::less<std::pair<double, int64_t>>,
                           NodesetAllocator<std::pair<double, int64_t>>>;

  struct OpenNode {
    std::vector<HighsDomainChange> domchgstack;
    std::vector<HighsInt> branchings;
    std::vector<NodeSet::iterator> domchglinks;
    double lower_bound;
    double estimate;
    HighsInt depth;
    highs::RbTreeLinks<int64_t> lowerLinks;
    highs::RbTreeLinks<int64_t> hybridEstimLinks;

    OpenNode()
        : domchgstack(),
          branchings(),
          domchglinks(),
          lower_bound(-kHighsInf),
          estimate(-kHighsInf),
          depth(0),
          lowerLinks(),
          hybridEstimLinks() {}

    OpenNode(std::vector<HighsDomainChange>&& domchgstack,
             std::vector<HighsInt>&& branchings, double lower_bound,
             double estimate, HighsInt depth)
        : domchgstack(std::move(domchgstack)),
          branchings(std::move(branchings)),
          lower_bound(lower_bound),
          estimate(estimate),
          depth(depth),
          lowerLinks(),
          hybridEstimLinks() {}

    OpenNode& operator=(OpenNode&& other) = default;
    OpenNode(OpenNode&&) = default;

    OpenNode& operator=(const OpenNode& other) = delete;
    OpenNode(const OpenNode&) = delete;
  };

  void checkGlobalBounds(HighsInt col, double lb, double ub, double feastol,
                         HighsCDouble& treeweight);

 private:
  class NodeLowerRbTree;
  class NodeHybridEstimRbTree;
  class SuboptimalNodeRbTree;

  std::unique_ptr<AllocatorState> allocatorState;
  std::vector<OpenNode> nodes;
  std::priority_queue<int64_t, std::vector<int64_t>, std::greater<int64_t>>
      freeslots;

  struct GlobalOperatorDelete {
    template <typename T>
    void operator()(T* x) const {
      ::operator delete(x);
    }
  };
  using NodeSetArray = std::unique_ptr<NodeSet, GlobalOperatorDelete>;
  NodeSetArray colLowerNodesPtr;
  NodeSetArray colUpperNodesPtr;
  int64_t lowerRoot = -1;
  int64_t lowerMin = -1;
  int64_t hybridEstimRoot = -1;
  int64_t hybridEstimMin = -1;
  int64_t suboptimalRoot = -1;
  int64_t suboptimalMin = -1;
  int64_t numSuboptimal = 0;
  double optimality_limit = kHighsInf;
  HighsInt numCol = 0;

  void link_estim(int64_t node);

  void unlink_estim(int64_t node);

  void link_lower(int64_t node);

  void unlink_lower(int64_t node);

  void link_suboptimal(int64_t node);

  void unlink_suboptimal(int64_t node);

  void link_domchgs(int64_t node);

  void unlink_domchgs(int64_t node);

  double link(int64_t node);

  void unlink(int64_t node);

 public:
  void setOptimalityLimit(double optimality_limit) {
    this->optimality_limit = optimality_limit;
  }

  double performBounding(double upper_limit);

  void setNumCol(HighsInt numcol);

  double emplaceNode(std::vector<HighsDomainChange>&& domchgs,
                     std::vector<HighsInt>&& branchings, double lower_bound,
                     double estimate, HighsInt depth);

  OpenNode&& popBestNode();

  OpenNode&& popBestBoundNode();

  int64_t numNodesUp(HighsInt col) const {
    return colLowerNodesPtr.get()[col].size();
  }

  int64_t numNodesDown(HighsInt col) const {
    return colUpperNodesPtr.get()[col].size();
  }

  int64_t numNodesUp(HighsInt col, double val) const {
    assert(numCol > col);
    auto colLowerNodes = colLowerNodesPtr.get();
    auto it = colLowerNodes[col].upper_bound(std::make_pair(val, kHighsIInf));
    if (it == colLowerNodes[col].begin()) return colLowerNodes[col].size();
    return std::distance(it, colLowerNodes[col].end());
  }

  int64_t numNodesDown(HighsInt col, double val) const {
    assert(numCol > col);
    auto colUpperNodes = colUpperNodesPtr.get();
    auto it = colUpperNodes[col].lower_bound(std::make_pair(val, -1));
    if (it == colUpperNodes[col].end()) return colUpperNodes[col].size();
    return std::distance(colUpperNodes[col].begin(), it);
  }

  const NodeSet& getUpNodes(HighsInt col) const {
    return colLowerNodesPtr.get()[col];
  }

  const NodeSet& getDownNodes(HighsInt col) const {
    return colUpperNodesPtr.get()[col];
  }

  double pruneInfeasibleNodes(HighsDomain& globaldomain, double feastol);

  double pruneNode(int64_t nodeId);

  double getBestLowerBound() const;

  HighsInt getBestBoundDomchgStackSize() const;

  void clear();

  int64_t numNodes() const { return nodes.size() - freeslots.size(); }

  int64_t numActiveNodes() const {
    return nodes.size() - freeslots.size() - numSuboptimal;
  }

  bool empty() const { return numActiveNodes() == 0; }
};

template <typename T, typename U>
bool operator==(const HighsNodeQueue::NodesetAllocator<T>&,
                const HighsNodeQueue::NodesetAllocator<U>&) {
  return true;
}

template <typename T, typename U>
bool operator!=(const HighsNodeQueue::NodesetAllocator<T>&,
                const HighsNodeQueue::NodesetAllocator<U>&) {
  return false;
}

#endif
