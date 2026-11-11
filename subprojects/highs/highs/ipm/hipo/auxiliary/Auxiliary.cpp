#include "Auxiliary.h"

#include <stack>

namespace hipo {

void inversePerm(const std::vector<Int>& perm, std::vector<Int>& iperm) {
  // Given the permutation perm, produce the inverse permutation iperm.
  // perm[i] : i-th entry to use in the new order.
  // iperm[i]: where entry i is located in the new order.

  for (Int i = 0; i < static_cast<Int>(perm.size()); ++i) {
    iperm[perm[i]] = i;
  }
}

void subtreeSize(const std::vector<Int>& parent, std::vector<Int>& sizes) {
  // Compute sizes of subtrees of the tree given by parent

  Int n = parent.size();
  sizes.assign(n, 1);

  for (Int i = 0; i < n; ++i) {
    Int k = parent[i];
    if (k != -1) sizes[k] += sizes[i];
  }
}

void transpose(Int n, Int nz, const Int* ptr, const Int* rows, Int* ptrT,
               Int* rowsT) {
  // Compute the transpose of the matrix and return it in rowsT and ptrT

  std::vector<Int> work(n);

  // count the entries in each row into work
  for (Int i = 0; i < nz; ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  counts2Ptr(n, ptrT, work.data());

  for (Int j = 0; j < n; ++j) {
    for (Int el = ptr[j]; el < ptr[j + 1]; ++el) {
      Int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      Int pos = work[i]++;
      rowsT[pos] = j;
    }
  }
}

void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               std::vector<Int>& ptrT, std::vector<Int>& rowsT) {
  transpose(ptr.size() - 1, rows.size(), ptr.data(), rows.data(), ptrT.data(),
            rowsT.data());
}

void transpose(Int n, Int nz, const Int* ptr, const Int* rows,
               const double* val, Int* ptrT, Int* rowsT, double* valT) {
  // Compute the transpose of the matrix and return it in rowsT, ptrT and valT

  std::vector<Int> work(n);

  // count the entries in each row into work
  for (Int i = 0; i < nz; ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  counts2Ptr(n, ptrT, work.data());

  for (Int j = 0; j < n; ++j) {
    for (Int el = ptr[j]; el < ptr[j + 1]; ++el) {
      Int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      Int pos = work[i]++;
      rowsT[pos] = j;
      valT[pos] = val[el];
    }
  }
}

void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               const std::vector<double>& val, std::vector<Int>& ptrT,
               std::vector<Int>& rowsT, std::vector<double>& valT) {
  transpose(ptr.size() - 1, rows.size(), ptr.data(), rows.data(), val.data(),
            ptrT.data(), rowsT.data(), valT.data());
}

void permuteSym(const std::vector<Int>& iperm, std::vector<Int>& ptr,
                std::vector<Int>& rows, std::vector<double>& val, bool lower) {
  // Symmetric permutation of the lower (upper) triangular matrix M based on
  // inverse permutation iperm. The resulting matrix is lower (upper)
  // triangular, regardless of the input matrix.

  const Int n = ptr.size() - 1;
  std::vector<Int> work(n, 0);
  const bool use_val = !val.empty();

  // go through the columns to count the nonzeros
  for (Int j = 0; j < n; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptr[j]; el < ptr[j + 1]; ++el) {
      const Int i = rows[el];

      // ignore potential entries in upper(lower) triangular part
      if ((lower && i < j) || (!lower && i > j)) continue;

      // get new index of row
      const Int row = iperm[i];

      // if only lower triangular part is used, col is smaller than row
      Int actual_col = lower ? std::min(row, col) : std::max(row, col);
      ++work[actual_col];
    }
  }

  std::vector<Int> new_ptr(n + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  counts2Ptr(n, new_ptr.data(), work.data());

  std::vector<Int> new_rows(new_ptr.back());
  std::vector<double> new_val;
  if (use_val) new_val.resize(new_ptr.back());

  // go through the columns to assign row indices
  for (Int j = 0; j < n; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptr[j]; el < ptr[j + 1]; ++el) {
      const Int i = rows[el];

      // ignore potential entries in upper triangular part
      if ((lower && i < j) || (!lower && i > j)) continue;

      // get new index of row
      const Int row = iperm[i];

      // if only lower triangular part is used, col is smaller than row
      const Int actual_col = lower ? std::min(row, col) : std::max(row, col);
      const Int actual_row = lower ? std::max(row, col) : std::min(row, col);

      Int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
      if (use_val) new_val[pos] = val[el];
    }
  }

  ptr = std::move(new_ptr);
  rows = std::move(new_rows);
  if (use_val) val = std::move(new_val);
}

void childrenLinkedList(const std::vector<Int>& parent, std::vector<Int>& head,
                        std::vector<Int>& next) {
  // Create linked lists of children in elimination tree.
  // parent gives the dependencies of the tree,
  // head[node] is the first child of node,
  // next[head[node]] is the second child,
  // next[next[head[node]]] is the third child...
  // until -1 is reached.

  Int n = parent.size();
  head.assign(n, -1);
  next.assign(n, -1);
  for (Int node = n - 1; node >= 0; --node) {
    if (parent[node] == -1) continue;
    next[node] = head[parent[node]];
    head[parent[node]] = node;
  }
}

void reverseLinkedList(std::vector<Int>& head, std::vector<Int>& next) {
  // Reverse the linked list of children of each node.
  // If a node has children (a -> b -> c -> -1), the reverse list contains
  // children (c -> b -> a -> -1).

  const Int n = head.size();

  for (Int node = 0; node < n; ++node) {
    Int prev_node = -1;
    Int curr_node = head[node];
    Int next_node = -1;

    while (curr_node != -1) {
      next_node = next[curr_node];
      next[curr_node] = prev_node;
      prev_node = curr_node;
      curr_node = next_node;
    }

    head[node] = prev_node;
  }
}

void dfsPostorder(Int node, Int& start, std::vector<Int>& head,
                  const std::vector<Int>& next, std::vector<Int>& order) {
  // Perform depth first search starting from root node and order the nodes
  // starting from the value start. head and next contain the linked list of
  // children.

  std::stack<Int> stack;
  stack.push(node);

  while (!stack.empty()) {
    const Int current = stack.top();
    const Int child = head[current];

    if (child == -1) {
      // no children left to order,
      // remove from the stack and order
      stack.pop();
      order[start++] = current;
    } else {
      // at least one child left to order,
      // add it to the stack and remove it from the list of children
      stack.push(child);
      head[current] = next[child];
    }
  }
}

void processEdge(Int j, Int i, const std::vector<Int>& first,
                 std::vector<Int>& maxfirst, std::vector<Int>& delta,
                 std::vector<Int>& prevleaf, std::vector<Int>& ancestor) {
  // Process edge of skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  // j not a leaf of ith row subtree
  if (i <= j || first[j] <= maxfirst[i]) {
    return;
  }

  // max first[j] so far
  maxfirst[i] = first[j];

  // previous leaf of ith row subtree
  Int jprev = prevleaf[i];

  // A(i,j) is in the skeleton matrix
  delta[j]++;

  if (jprev != -1) {
    // find least common ancestor of jprev and j
    Int q = jprev;
    while (q != ancestor[q]) {
      q = ancestor[q];
    }

    // path compression
    Int sparent;
    for (Int s = jprev; s != q; s = sparent) {
      sparent = ancestor[s];
      ancestor[s] = q;
    }

    // consider overlap
    delta[q]--;
  }

  // previous leaf of ith subtree set to j
  prevleaf[i] = j;
}

Int64 getDiagStart(Int n, Int k, Int nb, Int n_blocks,
                   std::vector<Int64>& start, bool triang) {
  // start position of diagonal blocks for blocked dense formats
  if (n_blocks == 0) return 0;

  start.assign(n_blocks, 0);
  for (Int i = 1; i < n_blocks; ++i) {
    start[i] = start[i - 1] + nb * (n - (i - 1) * nb);
    if (triang) start[i] -= nb * (nb - 1) / 2;
  }

  Int jb = std::min(nb, k - (n_blocks - 1) * nb);
  Int64 result = start.back() + (Int64)(n - (n_blocks - 1) * nb) * jb;
  if (triang) result -= jb * (jb - 1) / 2;
  return result;
}

Int maxDepthTree(const std::vector<Int>& parent) {
  Int max_depth = 0;
  Int n = parent.size();
  std::vector<Int> depth(n, -1);
  for (Int i = 0; i < n; ++i) {
    Int node = i;
    Int value = 1;
    while (node != -1) {
      if (value > depth[node]) {
        depth[node] = value;
      } else
        break;

      ++value;
      node = parent[node];
    }
    if (parent[i] == -1) max_depth = std::max(max_depth, depth[i]);
  }
  return max_depth;
}

void fullFromLower(Int n, Int nz, const Int* ptrL, const Int* rowsL,
                   std::vector<Int>& ptrF, std::vector<Int>& rowsF) {
  // Given a sparse matrix in lower triangular format, build the same matrix in
  // full format, without diagonal entries.

  std::vector<Int> rowsU(nz);
  std::vector<Int> ptrU(n + 1);
  transpose(n, nz, ptrL, rowsL, ptrU.data(), rowsU.data());

  std::vector<Int> work(n);
  for (Int j = 0; j < n; ++j) {
    for (Int el = ptrU[j]; el < ptrU[j + 1]; ++el) {
      const Int i = rowsU[el];
      if (i == j) continue;
      ++work[j];
      ++work[i];
    }
  }

  ptrF.assign(n + 1, 0);
  counts2Ptr(n, ptrF.data(), work.data());
  rowsF.assign(ptrF.back(), 0);
  for (Int j = 0; j < n; ++j) {
    for (Int el = ptrU[j]; el < ptrU[j + 1]; ++el) {
      const Int i = rowsU[el];
      if (i == j) continue;
      rowsF[work[j]++] = i;
      rowsF[work[i]++] = j;
    }
  }
}

double snFlops(double size, double clique_size) {
  return (size + clique_size) * (size + clique_size) * size -
         (size + clique_size) * size * (size + 1) +
         size * (size + 1) * (2 * size + 1) / 6;
}
double snSpops(double clique_size) {
  return clique_size * (clique_size + 1) / 2;
}

Clock::Clock() { start(); }
void Clock::start() { t0 = std::chrono::high_resolution_clock::now(); }
double Clock::stop() const {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> d = t1 - t0;
  return d.count();
}

TempTimer::TempTimer(const char* s) : name{s} {}
void TempTimer::start() { clock.start(); }
void TempTimer::stop() { time += clock.stop(); }
TempTimer::~TempTimer() { printf("Measured %10s: %e\n", name.c_str(), time); }

}  // namespace hipo
