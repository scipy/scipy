#ifndef FACTORHIGHS_CLIQUE_STACK_H
#define FACTORHIGHS_CLIQUE_STACK_H

#include <stack>
#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

// Class to manage the stack of cliques when processing the elimination tree in
// serial.
//
// The stack is used as follows:
// - use init to initialize the stack, providing the max stack size. If the
//   parameter is correct, there will be no reallocation of stack during its
//   operation.
// - use setup, with the size of the clique that is being computed, to
//   initialize enough elements in the workspace. This also returns a pointer to
//   write the new clique. The workspace is simply a chunk of space at the top
//   of the stack where the new clique is being computed before being pushed in
//   the right position.
// - use getChild to obtain information about the child clique that is at the
//   top of the stack. This also returns a pointer to read the data of the child
//   clique.
// - use popChild when the top child clique is no longer needed, to move on to
//   the next one.
// - use pushWork, with the ID of supernode that is being processed, when all
//   children clique have been assembled, the dense factorization is completed,
//   and the clique in the workspace is ready to be pushed onto the stack.

namespace hipo {

class CliqueStack {
  std::vector<double> stack_;
  double* workspace_;
  Int64 worksize_;
  bool empty_ = true;

  // pairs (sn, size) of supernodes that got pushed
  std::stack<std::pair<Int, Int64>> sn_pushed_{};

  Int64 top_{};

 public:
  void init(Int64 size);
  double* setup(Int64 clique_size, bool& reallocation);
  const double* getChild(Int& child_sn) const;
  void popChild();
  void pushWork(Int sn);
  bool empty() const;
};

}  // namespace hipo

// Example: sn 7 has children 5,6 and the stack looks like this:
//
// | - - - | - - - - - - - | - | - - - - - - - ...
// |   4   |       5       | 6 | t
// | - - - | - - - - - - - | - | - - - - - - - ...
//
// where t points to the top of the stack.
// The information about the supernodes pushed (in sn_pushed_) is:
// (sn 4, size 7), (sn 5, size 15), (sn 6, size 3)
//
// The workspace starts at t and fits in the stack.
// The clique of sn 7 is constructed in the workspace. After a child clique is
// assembled, it is popped from the stack, by moving t back. After all children
// are assembled and the dense factorization is done, the stack looks like this:
//
// | - - - | - - - - - - - - - | - - - - - - - - - - - | - - ...
// |   4   | t      free       |           7           |
// | - - - | - - - - - - - - - | - - - - - - - - - - - | - - ...
//
// The clique 7 can then be pushed onto the stack by copying it at position t.
// At the end, the stack looks like this:
//
// | - - - | - - - - - - - - - - - | - - - - - ...
// |   4   |           7           | t
// | - - - | - - - - - - - - - - - | - - - - - ...
//
// The information in sn_pushed_ is:
// (sn 4, size 7), (sn 7, size 23)

#endif