#include "CliqueStack.h"

#include <cassert>
#include <cstring>

namespace hipo {

void CliqueStack::init(Int64 stack_size) {
  top_ = 0;
  workspace_ = nullptr;
  worksize_ = 0;
  empty_ = false;

  stack_.resize(stack_size, 0.0);
}

double* CliqueStack::setup(Int64 clique_size, bool& reallocation) {
  // Clear workspace

  assert(!workspace_ && !worksize_);
  reallocation = false;

  if (clique_size > 0) {
    // This should not trigger reallocation, because the resize in init is done
    // with the maximum possible size of the stack.
    if (top_ + clique_size > static_cast<Int64>(stack_.size())) {
      reallocation = true;
      stack_.resize(top_ + clique_size, 0.0);
    }

    // accessing stack[top] is valid only if the clique is not empty,
    // otherwise it may be out of bounds.
    workspace_ = &stack_[top_];
    worksize_ = clique_size;

    // initialize workspace to zero
    std::memset(workspace_, 0, worksize_ * sizeof(double));
  }

  return workspace_;
}

bool CliqueStack::empty() const { return empty_; }

const double* CliqueStack::getChild(Int& child_sn) const {
  // Get the top of the stack, in terms of supernode ID of the child and pointer
  // to its data.

  child_sn = sn_pushed_.top().first;
  Int64 child_size = sn_pushed_.top().second;
  const double* child = &stack_[top_ - child_size];

  return child;
}

void CliqueStack::popChild() {
  // Remove top child from the stack

  Int64 child_size = sn_pushed_.top().second;
  sn_pushed_.pop();

  top_ -= child_size;
}

void CliqueStack::pushWork(Int sn) {
  // Put the content of the workspace at the top of the stack

  if (worksize_ > 0) {
    // stack_[top_] has lower address than workspace, so no need to resize.
    std::memmove(&stack_[top_], workspace_, worksize_ * sizeof(double));

    top_ += worksize_;

    // keep track of supernodes pushed
    sn_pushed_.push({sn, worksize_});
  }

  worksize_ = 0;
  workspace_ = nullptr;
}

}  // namespace hipo