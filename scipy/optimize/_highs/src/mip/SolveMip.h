/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_SOLVEMIP_H_
#define MIP_SOLVEMIP_H_

#include <cassert>
#include <functional>
#include <memory>
#include <stack>
#include <vector>

#include "lp_data/HConst.h"

struct Node {
  int id;
  int parent_id;
  double parent_objective;
  int level;

  Node();
  Node(int parent, double objective, int index, int depth)
      : id(index),
        parent_id(parent),
        parent_objective(objective),
        level(depth) {
    left_child = nullptr;
    right_child = nullptr;
    branch_col = -1;
  }

  std::vector<int> integer_variables;
  std::vector<double> primal_solution;
  double objective_value;

  // Minimal information about changes. Just col and its bounds for the moment.
  int branch_col;
  std::vector<double> col_lower_bound;
  std::vector<double> col_upper_bound;

  std::unique_ptr<Node> left_child;
  std::unique_ptr<Node> right_child;
};

using NodeIndex = int;
constexpr NodeIndex kNoNodeIndex = -1;
constexpr NodeIndex kNodeIndexError = -2;

class Tree {
 public:
  Tree() {}

  void pushRootNode(Node& node) {
    assert(nodes_.size() == 0);
    std::reference_wrapper<Node> ref(node);
    nodes_.push_back(ref);
  }

  Node& getRootNode() {
    assert(nodes_.size() > 0);
    return nodes_[0];
  }

  bool branch(Node& node);

  Node& next() { return nodes_[nodes_.size() - 1]; }
  void pop() { nodes_.erase(nodes_.end() - 1); }
  bool empty() { return (nodes_.size() == 0); }

  const std::vector<double>& getBestSolution() const { return best_solution_; }

  double getBestObjective() { return best_objective_; }
  double getBestBound(int& best_node);
  int getNumIntegerSolutions() { return num_integer_solutions; }
  int getNumNodesFormed() {
    return 1 + num_nodes;
  }  // Root node plus nodes formed by branching
  int getNumNodesLeft() { return (int)nodes_.size(); }
  void setMipReportLevel(const int mip_report_level_) {
    mip_report_level = mip_report_level_;
  }

 private:
  std::vector<std::reference_wrapper<Node> > nodes_;
  std::vector<double> best_solution_;
  double best_objective_ = HIGHS_CONST_INF;

  NodeIndex chooseBranchingVariable(const Node& node);

  int num_nodes = 0;
  int num_integer_solutions = 0;
  int mip_report_level = 0;
};

#endif
