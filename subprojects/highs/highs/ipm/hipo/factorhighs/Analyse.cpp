#include "Analyse.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stack>

#include "DataCollector.h"
#include "FactorHighsSettings.h"
#include "ReturnValues.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "util/HighsDisjointSets.h"

namespace hipo {

const Int64 int32_limit = std::numeric_limits<int32_t>::max();
const Int64 int64_limit = std::numeric_limits<int64_t>::max();

Analyse::Analyse(Int n, Int nz, const Int* rows, const Int* ptr,
                 const Int* signs, const FHoptions& FH_opt,
                 const Logger* logger, DataCollector& data, const Int* perm)
    : FH_opt_{FH_opt}, logger_{logger}, data_{data} {
  // Input the symmetric matrix to be analysed in CSC format.
  // rows contains the row indices.
  // ptr contains the starting points of each column.
  // Only the lower triangular part is used.
  // Diagonal entries must be stored for each column, even if zero.
  // signs contains the sign that each pivot should have.

  n_ = n;
  nz_ = nz;

  rows_lower_ = std::vector<Int>(rows, rows + nz_);
  ptr_lower_ = std::vector<Int>(ptr, ptr + n_ + 1);
  signs_ = std::vector<Int>(signs, signs + n_);
  perm_ = std::vector<Int>(perm, perm + n_);

  // adjust data if one-based indexing is used
  if (FH_opt_.one_indexing) {
    for (Int& i : rows_lower_) --i;
    for (Int& i : ptr_lower_) --i;
    for (Int& i : perm_) --i;
  }

  // Create upper triangular part
  rows_upper_.resize(nz_);
  ptr_upper_.resize(n_ + 1);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Permute the matrix with identical permutation, to extract upper triangular
  // part, if the input is not lower triangular.
  std::vector<Int> id_perm(n_);
  for (Int i = 0; i < n_; ++i) id_perm[i] = i;
  permute(id_perm);

  // actual number of nonzeros of only upper triangular part
  nz_ = ptr_upper_.back();

  // number of nonzeros potentially changed after Permute.
  rows_upper_.resize(nz_);

  // double transpose to sort columns
  ptr_lower_.resize(n_ + 1);
  rows_lower_.resize(nz_);
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  iperm_.resize(n_);
  inversePerm(perm_, iperm_);

  ready_ = true;
}

void Analyse::permute(const std::vector<Int>& iperm) {
  std::vector<double> empty_vals;
  permuteSym(iperm, ptr_upper_, rows_upper_, empty_vals, false);
}

void Analyse::eTree() {
  // Find elimination tree.
  // It works only for upper triangular matrices.
  // The tree is stored in the vector parent:
  //  parent[i] = j
  // means that j is the parent of i in the tree.
  // For the root(s) of the tree, parent[root] = -1.

  parent_.resize(n_);
  std::vector<Int> ancestor(n_);
  Int next{};

  for (Int j = 0; j < n_; ++j) {
    // initialise parent and ancestor, which are still unknown
    parent_[j] = -1;
    ancestor[j] = -1;

    for (Int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      for (Int i = rows_upper_[el]; i != -1 && i < j; i = next) {
        // next is used to move up the tree
        next = ancestor[i];

        // ancestor keeps track of the known part of the tree, to avoid
        // repeating (aka path compression): from j there is a known path to i
        ancestor[i] = j;

        if (next == -1) parent_[i] = j;
      }
    }
  }
}

void Analyse::postorder() {
  // Find a postordering of the elimination tree using depth first search

  postorder_.resize(n_);

  // create linked list of children
  std::vector<Int> head, next;
  childrenLinkedList(parent_, head, next);

  // Execute depth first search only for root node(s)
  Int start{};
  for (Int node = 0; node < n_; ++node) {
    if (parent_[node] == -1) {
      dfsPostorder(node, start, head, next, postorder_);
    }
  }

  // Permute elimination tree based on postorder
  std::vector<Int> ipost(n_);
  inversePerm(postorder_, ipost);
  std::vector<Int> new_parent(n_);
  for (Int i = 0; i < n_; ++i) {
    if (parent_[i] != -1) {
      new_parent[ipost[i]] = ipost[parent_[i]];
    } else {
      new_parent[ipost[i]] = -1;
    }
  }
  parent_ = std::move(new_parent);

  // Permute matrix based on postorder
  permute(ipost);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, postorder_);
  inversePerm(perm_, iperm_);
}

void Analyse::colCount() {
  // Columns count using skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  std::vector<Int> first(n_, -1);
  std::vector<Int> ancestor(n_, -1);
  std::vector<Int> max_first(n_, -1);
  std::vector<Int> prev_leaf(n_, -1);

  col_count_.resize(n_);

  // find first descendant
  for (Int k = 0; k < n_; ++k) {
    Int j = k;
    col_count_[j] = (first[j] == -1) ? 1 : 0;
    while (j != -1 && first[j] == -1) {
      first[j] = k;
      j = parent_[j];
    }
  }

  // each node belongs to a separate set
  for (Int j = 0; j < n_; j++) ancestor[j] = j;

  for (Int k = 0; k < n_; ++k) {
    const Int j = k;

    // if not a root, decrement
    if (parent_[j] != -1) col_count_[parent_[j]]--;

    // process edges of matrix
    for (Int el = ptr_lower_[j]; el < ptr_lower_[j + 1]; ++el) {
      processEdge(j, rows_lower_[el], first, max_first, col_count_, prev_leaf,
                  ancestor);
    }

    if (parent_[j] != -1) ancestor[j] = parent_[j];
  }

  // sum contributions from each child
  for (Int j = 0; j < n_; ++j) {
    if (parent_[j] != -1) {
      col_count_[parent_[j]] += col_count_[j];
    }
  }

  // compute nonzeros of L
  dense_ops_norelax_ = 0.0;
  nz_factor_ = 0;
  for (Int j = 0; j < n_; ++j) {
    nz_factor_ += col_count_[j];
    dense_ops_norelax_ += (double)(col_count_[j] - 1) * (col_count_[j] - 1);
  }
}

void Analyse::fundamentalSupernodes() {
  // Find fundamental supernodes.

  // isSN[i] is true if node i is the start of a fundamental supernode
  std::vector<HighsBool> is_sn(n_, false);

  std::vector<Int> prev_nonz(n_, -1);

  // compute sizes of subtrees
  std::vector<Int> subtree_sizes(n_);
  subtreeSize(parent_, subtree_sizes);

  for (Int j = 0; j < n_; ++j) {
    for (Int el = ptr_lower_[j]; el < ptr_lower_[j + 1]; ++el) {
      const Int i = rows_lower_[el];
      const Int k = prev_nonz[i];

      // mark as fundamental sn, nodes which are leaf of subtrees
      if (k < j - subtree_sizes[j] + 1) {
        is_sn[j] = true;
      }

      // mark as fundamental sn, nodes which have more than one child
      if (parent_[i] != -1 &&
          subtree_sizes[i] + 1 != subtree_sizes[parent_[i]]) {
        is_sn[parent_[i]] = true;
      }

      prev_nonz[i] = j;
    }
  }

  // create information about fundamental supernodes
  sn_belong_.resize(n_);
  Int sn_number = -1;
  for (Int i = 0; i < n_; ++i) {
    // if isSN[i] is true, then node i is the start of a new supernode
    if (is_sn[i]) ++sn_number;

    // mark node i as belonging to the current supernode
    sn_belong_[i] = sn_number;
  }

  // number of supernodes found
  sn_count_ = sn_belong_.back() + 1;

  // sn_start_ contains pointers to the starting node of each supernode
  sn_start_.resize(sn_count_ + 1);
  Int next = 0;
  for (Int i = 0; i < n_; ++i) {
    if (is_sn[i]) {
      sn_start_[next] = i;
      ++next;
    }
  }
  sn_start_[next] = n_;

  // build supernodal elimination tree
  sn_parent_.resize(sn_count_);
  for (Int i = 0; i < sn_count_ - 1; ++i) {
    Int j = parent_[sn_start_[i + 1] - 1];
    if (j != -1) {
      sn_parent_[i] = sn_belong_[j];
    } else {
      sn_parent_[i] = -1;
    }
  }
  sn_parent_.back() = -1;
}

void Analyse::relaxSupernodes() {
  // Choose which strategy to use for supernode amalgamation.
  // Relaxing supernodes is fast to run, so just run again the one that is
  // chosen to set the correct information.

  double flops_1{};
  double spops_1{};
  relaxSnMaxNz(flops_1, spops_1);
  double ops_1 = flops_1 + kSpopsWeightSn * spops_1;

  double flops_2{};
  double spops_2{};
  relaxSnNetOps(flops_2, spops_2);
  double ops_2 = flops_2 + kSpopsWeightSn * spops_2;

  bool total_ops_reduced = ops_2 < ops_1;
  bool spops_reduced_enough = spops_2 < spops_1 * 0.6;

  if (total_ops_reduced && spops_reduced_enough) {
    relaxSnNetOps(flops_2, spops_2);
  } else {
    relaxSnMaxNz(flops_1, spops_1);
  }
}

void Analyse::relaxSnNetOps(double& flops, double& spops) {
  // =================================================
  // Build information about supernodes
  // =================================================
  std::vector<Int> sn_size(sn_count_);
  std::vector<Int> clique_size(sn_count_);
  fake_nz_.assign(sn_count_, 0);
  for (Int i = 0; i < sn_count_; ++i) {
    sn_size[i] = sn_start_[i + 1] - sn_start_[i];
    clique_size[i] = col_count_[sn_start_[i]] - sn_size[i];
  }

  // build linked lists of children
  std::vector<Int> first_child, next_child;
  childrenLinkedList(sn_parent_, first_child, next_child);

  // =================================================
  // Merge supernodes
  // =================================================
  merged_into_.assign(sn_count_, -1);
  merged_sn_ = 0;

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // sort children by increasing size
    std::vector<Int> children;
    Int child = first_child[sn];
    while (child != -1) {
      children.push_back(child);
      child = next_child[child];
    }
    std::sort(children.begin(), children.end(),
              [&](Int a, Int b) { return sn_size[a] < sn_size[b]; });

    for (Int child : children) {
      // how many zero rows would become nonzero
      const Int rows_filled =
          sn_size[sn] + clique_size[sn] - clique_size[child];

      // how many zero entries would become nonzero
      const Int64 nz_added = (Int64)rows_filled * sn_size[child];

      // how many artificial nonzeros would the merged supernode have
      const Int64 total_art_nz = nz_added + fake_nz_[sn] + fake_nz_[child];

      const double flops_added =
          rows_filled * rows_filled * sn_size[child] +
          2 * rows_filled *
              (sn_size[child] * (sn_size[child] + clique_size[child]) -
               sn_size[child] * (sn_size[child] + 1) / 2);

      const double spops_removed =
          clique_size[child] * (clique_size[child] + 1) / 2;

      const double net_ops = flops_added - kSpopsWeightSn * spops_removed;

      if (net_ops < 0) {
        // merge child with parent
        sn_size[sn] += sn_size[child];
        fake_nz_[sn] = total_art_nz;
        ++merged_sn_;
        merged_into_[child] = sn;
      }
    }
  }

  flops = 0.0;
  spops = 0.0;
  for (Int i = 0; i < sn_count_; ++i) {
    if (merged_into_[i] == -1) {
      spops += snSpops(clique_size[i]);
      flops += snFlops(sn_size[i], clique_size[i]);
    }
  }
}

double Analyse::doRelaxSnMaxNz(double& flops, double& spops,
                               Int64 max_artificial_nz) {
  // =================================================
  // Build information about supernodes
  // =================================================
  std::vector<Int> sn_size(sn_count_);
  std::vector<Int> clique_size(sn_count_);
  fake_nz_.assign(sn_count_, 0);
  for (Int i = 0; i < sn_count_; ++i) {
    sn_size[i] = sn_start_[i + 1] - sn_start_[i];
    clique_size[i] = col_count_[sn_start_[i]] - sn_size[i];
    fake_nz_[i] = 0;
  }

  // build linked lists of children
  std::vector<Int> first_child, next_child;
  childrenLinkedList(sn_parent_, first_child, next_child);

  // =================================================
  // Merge supernodes
  // =================================================
  merged_into_.assign(sn_count_, -1);
  merged_sn_ = 0;

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // keep iterating through the children of the supernode, until there's no
    // more child to merge with

    while (true) {
      Int child = first_child[sn];

      // info for first criterion
      Int64 nz_fakenz = int64_limit;
      Int size_fakenz = 0;
      Int child_fakenz = -1;

      while (child != -1) {
        // how many zero rows would become nonzero
        const Int rows_filled =
            sn_size[sn] + clique_size[sn] - clique_size[child];

        // how many zero entries would become nonzero
        const Int64 nz_added = (Int64)rows_filled * sn_size[child];

        // how many artificial nonzeros would the merged supernode have
        const Int64 total_art_nz = nz_added + fake_nz_[sn] + fake_nz_[child];

        // Save child with smallest number of artificial zeros created.
        // Ties are broken based on size of child.
        if (total_art_nz < nz_fakenz ||
            (total_art_nz == nz_fakenz && size_fakenz < sn_size[child])) {
          nz_fakenz = total_art_nz;
          size_fakenz = sn_size[child];
          child_fakenz = child;
        }

        child = next_child[child];
      }

      if (nz_fakenz <= max_artificial_nz) {
        // merging creates fewer nonzeros than the maximum allowed

        // update information of parent
        sn_size[sn] += size_fakenz;
        fake_nz_[sn] = nz_fakenz;

        // count number of merged supernodes
        ++merged_sn_;

        // save information about merging of supernodes
        merged_into_[child_fakenz] = sn;

        // remove child from linked list of children
        child = first_child[sn];
        if (child == child_fakenz) {
          // child_smallest is the first child
          first_child[sn] = next_child[child_fakenz];
        } else {
          while (next_child[child] != child_fakenz) {
            child = next_child[child];
          }
          // now child is the previous child of child_smallest
          next_child[child] = next_child[child_fakenz];
        }

      } else {
        // no more children can be merged with parent
        break;
      }
    }
  }

  flops = 0.0;
  spops = 0.0;
  for (Int i = 0; i < sn_count_; ++i) {
    if (merged_into_[i] == -1) {
      spops += snSpops(clique_size[i]);
      flops += snFlops(sn_size[i], clique_size[i]);
    }
  }

  // compute total number of artificial nonzeros and artificial ops for this
  // value of max_artificial_nz
  double temp_art_nz{};
  double temp_art_ops{};
  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (merged_into_[sn] == -1) {
      temp_art_nz += fake_nz_[sn];
      temp_art_ops += snFlops(sn_size[sn], clique_size[sn]);
    }
  }
  temp_art_ops -= dense_ops_norelax_;

  // if enough fake nz or ops have been added, stop.
  const double ratio_fake = temp_art_ops / (temp_art_ops + dense_ops_norelax_);

  return ratio_fake;
}

void Analyse::relaxSnMaxNz(double& flops, double& spops) {
  // Child which produces smallest number of fake nonzeros is merged if
  // resulting sn has fewer than max_artificial_nz fake nonzeros.
  // Multiple values of max_artificial_nz are tried, chosen with bisection
  // method, until the percentage of artificial nonzeros is in the range [1,2]%.

  Int64 max_artificial_nz = kStartThreshRelax;
  Int64 largest_below = -1;
  Int64 smallest_above = -1;

  double best_dist_ratio = kHighsInf;
  Int64 best_max_art_nz = -1;

  for (Int iter = 0; iter < kMaxIterRelax; ++iter) {
    // relax the supernodes and obtain the ratio of how many new ops have been
    // added with the current value of max_artificial_nz
    const double ratio_fake = doRelaxSnMaxNz(flops, spops, max_artificial_nz);

    // store the best ratio, in case a good ratio is never found
    double dist_ratio_fake = std::min(std::abs(ratio_fake - kLowerRatioRelax),
                                      std::abs(ratio_fake - kUpperRatioRelax));
    if (dist_ratio_fake < best_dist_ratio) {
      best_dist_ratio = dist_ratio_fake;
      best_max_art_nz = max_artificial_nz;
    }

    // try to find ratio in interval [0.01,0.02] using bisection
    if (ratio_fake < kLowerRatioRelax) {
      // ratio too small
      largest_below = max_artificial_nz;
      if (smallest_above == -1) {
        max_artificial_nz *= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else if (ratio_fake > kUpperRatioRelax) {
      // ratio too large
      smallest_above = max_artificial_nz;
      if (largest_below == -1) {
        max_artificial_nz /= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else {
      // good ratio
      return;
    }
  }

  // If reach here, no good ratio was found within kMaxIterRelax
  // To avoid having a catastrophically bad ratio in pathological problems,
  // choose the best ratio found

  doRelaxSnMaxNz(flops, spops, best_max_art_nz);
}

void Analyse::afterRelaxSn() {
  // number of new supernodes
  const Int new_snCount = sn_count_ - merged_sn_;

  // keep track of number of row indices needed for each supernode
  sn_indices_.assign(new_snCount, 0);

  // =================================================
  // Create supernodal permutation
  // =================================================

  // permutation of supernodes needed after merging
  std::vector<Int> sn_perm(sn_count_);

  // number of new sn that includes the old sn
  std::vector<Int> new_id(sn_count_);

  // new sn pointer vector
  std::vector<Int> new_snStart(new_snCount + 1);

  // keep track of the children merged into a given supernode
  std::vector<std::vector<Int>> received_from(sn_count_, std::vector<Int>());

  // index to write into sn_perm
  Int start_perm{};

  // index to write into new_snStart
  Int snStart_ind{};

  // next available number for new sn numbering
  Int next_id{};

  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (merged_into_[sn] > -1) {
      // Current sn was merged into its parent.
      // Save information about which supernode sn was merged into
      received_from[merged_into_[sn]].push_back(sn);
    } else {
      // Current sn was not merged into its parent.
      // It is one of the new sn.

      // Add merged supernodes to the permutation, recursively.

      ++snStart_ind;

      std::stack<Int> toadd;
      toadd.push(sn);

      while (!toadd.empty()) {
        const Int current = toadd.top();

        if (!received_from[current].empty()) {
          for (Int i : received_from[current]) toadd.push(i);
          received_from[current].clear();
        } else {
          toadd.pop();
          sn_perm[start_perm++] = current;
          new_id[current] = next_id;

          // count number of nodes in each new supernode
          new_snStart[snStart_ind] +=
              sn_start_[current + 1] - sn_start_[current];
        }
      }

      // keep track of total number of artificial nonzeros
      artificial_nz_ += fake_nz_[sn];

      // Compute number of indices for new sn.
      // This is equal to the number of columns in the new sn plus the clique
      // size of the original supernode where the children where merged.
      sn_indices_[next_id] = new_snStart[snStart_ind] +
                             col_count_[sn_start_[sn]] - sn_start_[sn + 1] +
                             sn_start_[sn];

      ++next_id;
    }
  }

  // new_snStart contain the number of cols in each new sn.
  // sum them to obtain the sn pointers.
  for (Int i = 0; i < new_snCount; ++i) {
    new_snStart[i + 1] += new_snStart[i];
  }

  // include artificial nonzeros in the nonzeros of the factor
  nz_factor_ += artificial_nz_;

  // compute number of flops needed for the factorisation
  dense_ops_ = 0.0;
  for (Int sn = 0; sn < new_snCount; ++sn) {
    const double colcount_sn = (double)sn_indices_[sn];
    for (Int i = 0; i < new_snStart[sn + 1] - new_snStart[sn]; ++i) {
      dense_ops_ += (colcount_sn - i - 1) * (colcount_sn - i - 1);
    }
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation needed after
  // sn merging.

  // permutation to apply to the existing one
  std::vector<Int> new_perm(n_);

  // index to write into new_perm
  Int start{};

  for (Int i = 0; i < sn_count_; ++i) {
    const Int sn = sn_perm[i];
    for (Int j = sn_start_[sn]; j < sn_start_[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<Int> new_iperm(n_);
  inversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<Int> new_snParent(new_snCount, -1);
  for (Int i = 0; i < sn_count_; ++i) {
    if (sn_parent_[i] == -1) continue;

    const Int ii = new_id[i];
    const Int pp = new_id[sn_parent_[i]];

    if (ii == pp) continue;

    new_snParent[ii] = pp;
  }

  // =================================================
  // Save new information
  // =================================================

  // build new snBelong, i.e., the sn to which each column belongs
  for (Int sn = 0; sn < sn_count_; ++sn) {
    for (Int i = sn_start_[sn]; i < sn_start_[sn + 1]; ++i) {
      sn_belong_[i] = new_id[sn];
    }
  }
  permuteVector(sn_belong_, new_perm);

  permuteVector(col_count_, new_perm);

  // Overwrite previous data
  sn_parent_ = std::move(new_snParent);
  sn_start_ = std::move(new_snStart);
  sn_count_ = new_snCount;

  // Permute matrix based on new permutation
  permute(new_iperm);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, new_perm);
  inversePerm(perm_, iperm_);
}

void Analyse::snPattern() {
  // number of total indices needed
  Int64 indices{};

  for (Int64 i : sn_indices_) indices += i;

  // allocate space for sn pattern
  rows_sn_.resize(indices);
  ptr_sn_.resize(sn_count_ + 1);

  // keep track of visited supernodes
  std::vector<Int> mark(sn_count_, -1);

  // compute column pointers of L
  std::vector<Int64> work(sn_indices_.size());
  for (Int i = 0; i < static_cast<Int>(sn_indices_.size()); ++i)
    work[i] = sn_indices_[i];
  counts2Ptr(work.size(), ptr_sn_.data(), work.data());

  // consider each row
  for (Int i = 0; i < n_; ++i) {
    // for all entries in the row of lower triangle
    for (Int el = ptr_upper_[i]; el < ptr_upper_[i + 1]; ++el) {
      // there is nonzero (i,j)
      const Int j = rows_upper_[el];

      // supernode to which column j belongs to
      Int snj = sn_belong_[j];

      // while supernodes are not yet considered
      while (snj != -1 && mark[snj] != i) {
        // we may end up too far
        if (sn_start_[snj] > i) break;

        // supernode snj is now considered for row i
        mark[snj] = i;

        // there is a nonzero entry in supernode snj at row i
        rows_sn_[work[snj]++] = i;

        // go up the elimination tree
        snj = sn_parent_[snj];
      }
    }
  }
}

void Analyse::relativeIndCols() {
  // Find the relative indices of the original column wrt the frontal matrix of
  // the corresponding supernode

  relind_cols_.resize(nz_);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int64 ptL_start = ptr_sn_[sn];
    const Int64 ptL_end = ptr_sn_[sn + 1];

    // go through the columns of the supernode
    for (Int col = sn_start_[sn]; col < sn_start_[sn + 1]; ++col) {
      // go through original column and supernodal column
      Int ptA = ptr_lower_[col];
      Int64 ptL = ptL_start;

      // offset wrt ptrLower[col]
      Int index{};

      // size of the column of the original matrix
      Int col_size = ptr_lower_[col + 1] - ptr_lower_[col];

      while (ptL < ptL_end) {
        // if found all the relative indices that are needed, stop
        if (index == col_size) {
          break;
        }

        // check if indices coincide
        if (rows_sn_[ptL] == rows_lower_[ptA]) {
          // yes: save relative index and move pointers forward
          // NB: difference fits into Int
          relind_cols_[ptr_lower_[col] + index] = ptL - ptL_start;
          ++index;
          ++ptL;
          ++ptA;
        } else {
          // no: move pointer of L forward
          ++ptL;
        }
      }
    }
  }
}

void Analyse::relativeIndClique() {
  // Find the relative indices of the child clique wrt the frontal matrix of the
  // parent supernode

  relind_clique_.resize(sn_count_);
  consecutive_sums_.resize(sn_count_);

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // if there is no parent, skip supernode
    if (sn_parent_[sn] == -1) continue;

    // number of nodes in the supernode
    const Int sn_size = sn_start_[sn + 1] - sn_start_[sn];

    // column of the first node in the supernode
    const Int j = sn_start_[sn];

    // size of the first column of the supernode
    const Int sn_column_size = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // size of the clique of the supernode
    const Int sn_clique_size = sn_column_size - sn_size;

    // count number of assembly operations during factorise
    sparse_ops_ += snSpops(sn_clique_size);

    relind_clique_[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    Int64 ptr_current = ptr_sn_[sn] + sn_size;

    // iterate through the full column of parent sn
    Int64 ptr_parent = ptr_sn_[sn_parent_[sn]];

    // keep track of start and end of parent sn column
    const Int64 ptr_parent_start = ptr_parent;
    const Int64 ptr_parent_end = ptr_sn_[sn_parent_[sn] + 1];

    // where to write into relind
    Int index{};

    // iterate though the column of the parent sn
    while (ptr_parent < ptr_parent_end) {
      // if found all the relative indices that are needed, stop
      if (index == sn_clique_size) {
        break;
      }

      // check if indices coincide
      if (rows_sn_[ptr_current] == rows_sn_[ptr_parent]) {
        // yes: save relative index and move pointers forward
        // NB: difference fits into Int
        relind_clique_[sn][index] = ptr_parent - ptr_parent_start;
        ++index;
        ++ptr_parent;
        ++ptr_current;
      } else {
        // no: move pointer of parent forward
        ++ptr_parent;
      }
    }

    // Difference between consecutive relative indices.
    // Useful to detect chains of consecutive indices.
    consecutive_sums_[sn].resize(sn_clique_size);
    for (Int i = 0; i < sn_clique_size - 1; ++i) {
      consecutive_sums_[sn][i] =
          relind_clique_[sn][i + 1] - relind_clique_[sn][i];
    }

    // Number of consecutive sums that can be done in one blas call.
    consecutive_sums_[sn].back() = 1;
    for (Int i = sn_clique_size - 2; i >= 0; --i) {
      if (consecutive_sums_[sn][i] > 1) {
        consecutive_sums_[sn][i] = 1;
      } else if (consecutive_sums_[sn][i] == 1) {
        consecutive_sums_[sn][i] = consecutive_sums_[sn][i + 1] + 1;
      } else {
        if (logger_) logger_->printInfo("Error in consecutiveSums\n");
      }
    }
  }
}

void Analyse::computeStorage(Int fr, Int sz, Int64& fr_entries,
                             Int64& cl_entries) const {
  // compute storage required by frontal and clique, based on the format used

  const Int cl = fr - sz;
  const Int nb = FH_opt_.nb;

  Int n_blocks = (sz - 1) / nb + 1;
  std::vector<Int64> temp;
  fr_entries = getDiagStart(fr, sz, nb, n_blocks, temp);

  // clique is stored as a collection of rectangles
  n_blocks = (cl - 1) / nb + 1;
  Int64 schur_size{};
  for (Int j = 0; j < n_blocks; ++j) {
    const Int jb = std::min(nb, cl - j * nb);
    schur_size += (Int64)(cl - j * nb) * jb;
  }
  cl_entries = schur_size;
}

void Analyse::computeCriticalPath() {
  // Compute the critical path within the supernodal elimination tree, and the
  // number of operations along the path. This is the number of operations that
  // need to be done sequentially while doing tree parallelism.

  std::vector<double> critical_ops(sn_count_);

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    critical_ops[sn] = snFlops(sz, fr - sz);
  }

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf nodes
    if (head[sn] == -1) continue;

    double max_ops{};
    Int child = head[sn];
    while (child != -1) {
      // critical_ops of this supernode is max over children of
      // (ops_of_this_sn + critical_ops_of_child)
      max_ops = std::max(max_ops, critical_ops[sn] + critical_ops[child]);
      child = next[child];
    }
    critical_ops[sn] = max_ops;
  }

  for (Int sn = 0; sn < sn_count_; ++sn) {
    critical_ops_ = std::max(critical_ops_, critical_ops[sn]);
  }
}

void Analyse::computeCriticalPathSolve() {
  // Compute the critical path within the task elimination tree, and the
  // number of operations along the path. This is the number of operations that
  // need to be done sequentially while doing tree parallelism.

  std::vector<double> critical_ops(schedule_solve_.count());

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(schedule_solve_.task_parent, head, next);

  ops_solve_ = 0.0;
  critical_ops_solve_ = 0.0;

  for (Int task = 0; task < schedule_solve_.count(); ++task) {
    critical_ops[task] = task_ops_solve_[task];
    ops_solve_ += task_ops_solve_[task];
  }

  for (Int task = 0; task < schedule_solve_.count(); ++task) {
    // leaf task
    if (head[task] == -1) continue;

    double max_ops{};
    Int child = head[task];
    while (child != -1) {
      // critical_ops of this supernode is max over children of
      // (ops_of_this_task + critical_ops_of_child)
      max_ops = std::max(max_ops, critical_ops[task] + critical_ops[child]);
      child = next[child];
    }
    critical_ops[task] = max_ops;
  }

  for (Int task = 0; task < schedule_solve_.count(); ++task) {
    critical_ops_solve_ = std::max(critical_ops_solve_, critical_ops[task]);
  }
}

void Analyse::reorderChildren() {
  std::vector<Int64> clique_entries(sn_count_);
  std::vector<Int64> frontal_entries(sn_count_);
  std::vector<double> storage(sn_count_);
  std::vector<double> storage_factors(sn_count_);

  // initialise data of supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = col_count_[sn_start_[sn]];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);

    // compute number of entries in factors within the subtree
    storage_factors[sn] += frontal_entries[sn];
    if (sn_parent_[sn] != -1)
      storage_factors[sn_parent_[sn]] += storage_factors[sn];
  }

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      storage[sn] = frontal_entries[sn] + clique_entries[sn];
      continue;
    }

    // save children and values to sort
    std::vector<std::pair<Int, double>> children{};
    Int child = head[sn];
    while (child != -1) {
      double value =
          storage[child] - clique_entries[child] - storage_factors[child];
      children.push_back({child, value});
      child = next[child];
    }

    // sort children in decreasing order of the values
    std::sort(children.begin(), children.end(),
              [&](std::pair<Int, double>& a, std::pair<Int, double>& b) {
                return a.second > b.second;
              });

    // modify linked lists with new order of children
    head[sn] = children.front().first;
    for (Int i = 0; i < static_cast<Int>(children.size()) - 1; ++i) {
      next[children[i].first] = children[i + 1].first;
    }
    next[children.back().first] = -1;
  }

  // =================================================
  // Create supernodal permutation
  // =================================================
  // build supernodal permutation with dfs
  std::vector<Int> sn_perm(sn_count_);
  Int start{};
  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (sn_parent_[sn] == -1) dfsPostorder(sn, start, head, next, sn_perm);
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation

  // permutation to apply to the existing one
  std::vector<Int> new_perm(n_);

  // index to write into new_perm
  start = 0;

  for (Int i = 0; i < sn_count_; ++i) {
    const Int sn = sn_perm[i];
    for (Int j = sn_start_[sn]; j < sn_start_[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<Int> new_iperm(n_);
  inversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<Int> isn_perm(sn_count_);
  inversePerm(sn_perm, isn_perm);
  std::vector<Int> new_sn_parent(sn_count_);
  for (Int i = 0; i < sn_count_; ++i) {
    if (sn_parent_[i] != -1) {
      new_sn_parent[isn_perm[i]] = isn_perm[sn_parent_[i]];
    } else {
      new_sn_parent[isn_perm[i]] = -1;
    }
  }

  // =================================================
  // Create new snBelong
  // =================================================

  // build new snBelong, i.e., the sn to which each colum belongs
  for (Int sn = 0; sn < sn_count_; ++sn) {
    for (Int i = sn_start_[sn]; i < sn_start_[sn + 1]; ++i) {
      sn_belong_[i] = isn_perm[sn];
    }
  }
  permuteVector(sn_belong_, new_perm);

  // permute other vectors that may be needed
  permuteVector(col_count_, new_perm);
  permuteVector(sn_indices_, sn_perm);

  // =================================================
  // Create new snStart
  // =================================================
  std::vector<Int> cols_per_sn(sn_count_);

  // compute size of each supernode
  for (Int sn = 0; sn < sn_count_; ++sn) {
    cols_per_sn[sn] = sn_start_[sn + 1] - sn_start_[sn];
  }

  // permute according to new order of supernodes
  permuteVector(cols_per_sn, sn_perm);

  // sum number of columns to obtain pointers
  for (Int i = 0; i < sn_count_ - 1; ++i) {
    cols_per_sn[i + 1] += cols_per_sn[i];
  }

  for (Int i = 0; i < sn_count_; ++i) {
    sn_start_[i + 1] = cols_per_sn[i];
  }

  // =================================================
  // Save new data
  // =================================================

  // Overwrite previous data
  sn_parent_ = std::move(new_sn_parent);

  // Permute matrix based on new permutation
  permute(new_iperm);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, new_perm);
  inversePerm(perm_, iperm_);
}

void Analyse::computeBlockStart() {
  const Int nb = FH_opt_.nb;

  clique_block_start_.resize(sn_count_);
  // compute starting position of each block of columns in the clique, for
  // each supernode
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int sn_size = sn_start_[sn + 1] - sn_start_[sn];
    const Int ldf = ptr_sn_[sn + 1] - ptr_sn_[sn];
    const Int ldc = ldf - sn_size;
    const Int n_blocks = (ldc - 1) / nb + 1;

    Int64 schur_size =
        getDiagStart(ldc, ldc, nb, n_blocks, clique_block_start_[sn]);
    clique_block_start_[sn].push_back(schur_size);
  }
}

Int Analyse::checkOverflow() const {
  // In order to use 32-bit BLAS, any data accessed by BLAS must be addressable
  // using 32-bit integer offset. If BLAS is given a pointer double* A, the
  // distance between the first and last entry of A used by BLAS needs to be
  // smaller than int32_limit. Since the matrices are stored in blocked data
  // structures, and BLAS only uses contiguous data from a given block of
  // columns, we need to impose that:
  //   front_size * min(block_size, sn_size) <= int32_limit

  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int sn_size = sn_start_[sn + 1] - sn_start_[sn];
    const Int front_size = ptr_sn_[sn + 1] - ptr_sn_[sn];

    if ((Int64)front_size * std::min(sn_size, FH_opt_.nb) > int32_limit)
      return 1;
  }

  return 0;
}

void Analyse::computeStackSize() {
  // Compute the minimum size of the stack to process the elimination tree
  // serially.

  std::vector<Int64> clique_entries(sn_count_);
  std::vector<Int64> stack_subtrees(sn_count_);
  Int64 total_frontal{};

  // initialise data of supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    Int64 frontal_entries{};

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries, clique_entries[sn]);

    total_frontal += frontal_entries;
  }

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      stack_subtrees[sn] = clique_entries[sn];
      continue;
    }

    // Compute storage
    // storage is found as max(storage_1,storage_2), where
    // storage_1 = max_j stack_size[j] + \sum_{k up to j-1} clique_entries[k]
    // storage_2 = clique_total_entries (including node itself)

    Int64 clique_partial_entries{};
    Int64 storage_1{};

    Int child = head[sn];
    while (child != -1) {
      Int64 current = stack_subtrees[child] + clique_partial_entries;

      clique_partial_entries += clique_entries[child];
      storage_1 = std::max(storage_1, current);

      child = next[child];
    }

    Int64 storage_2 = clique_partial_entries + clique_entries[sn];

    stack_subtrees[sn] = std::max(storage_1, storage_2);
    max_stack_size_ = std::max(max_stack_size_, stack_subtrees[sn]);
  }

  // minimum storage in serial is equal to the space needed to store the
  // factorisation and the maximum size of the stack. Times 8 to obtain the
  // number of bytes.
  serial_storage_ = (total_frontal + max_stack_size_) * 8;
}

void Analyse::computeTreeScheduleSolve() {
  // compute number of operations for each supernode
  std::vector<double> sn_ops(sn_count_);
  double total_ops = 0;
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];
    const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];
    const double this_sn_dense_ops =
        (double)sz * (sz + 1) / 2 + (double)sz * (fr - sz);
    sn_ops[sn] += this_sn_dense_ops;
    total_ops += this_sn_dense_ops;
  }

  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  const double task_ops_thresh =
      std::max(total_ops * kLargeTaskRelativeThresh, kLargeTaskAbsoluteThres);

  HighsDisjointSets<> sets(sn_count_);
  std::vector<double> child_ops(sn_count_, 0.0);
  std::map<Int, Int> task_numbering;
  Int task_count = 0;

  // Assign supernodes to tasks:
  // if a supernode is part of a task that is large enough, then the task is
  // considered complete. Otherwise, the task is still open and the parent
  // supernode will continue it.
  for (Int sn = 0; sn < sn_count_; ++sn) {
    double this_sn_task_ops = sn_ops[sn];
    Int child = head[sn];
    while (child != -1) {
      this_sn_task_ops += child_ops[child];
      child = next[child];
    }

    const bool task_large = this_sn_task_ops > task_ops_thresh;
    const bool task_root = sn_parent_[sn] == -1;

    if (task_large || task_root) {
      // completed a task
      child_ops[sn] = 0.0;
      task_numbering[sets.getSet(sn)] = task_count;
      task_count++;
    } else {
      // task still incomplete
      child_ops[sn] = this_sn_task_ops;
      sets.merge(sn, sn_parent_[sn]);
    }
  }

  schedule_solve_.sn_per_task.resize(task_count);
  task_ops_solve_.assign(task_count, 0.0);
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int task_id = task_numbering[sets.getSet(sn)];
    schedule_solve_.sn_per_task[task_id].push_back(sn);
    task_ops_solve_[task_id] += sn_ops[sn];
  }

  // Create tree of dependencies among tasks.
  // Since the supernodal tree is postordered, the tree of dependencies among
  // tasks should be automatically postordered as well.
  schedule_solve_.task_parent.assign(task_count, -1);
  for (Int task = 0; task < task_count; ++task) {
    for (Int sn : schedule_solve_.sn_per_task[task]) {
      Int child = head[sn];
      while (child != -1) {
        const Int child_task = task_numbering[sets.getSet(child)];
        if (child_task != task) {
          assert(child_task < task);
          assert(schedule_solve_.task_parent[child_task] == -1);
          schedule_solve_.task_parent[child_task] = task;
        }
        child = next[child];
      }
    }
  }

  schedule_solve_.valid = true;

  // verify that task elimination tree has topological ordering
  for (Int task = 0; task < task_count; ++task) {
    const Int this_parent = schedule_solve_.task_parent[task];
    if (this_parent != -1 && this_parent <= task) {
      schedule_solve_.clear();
      logger_->printInfo("Task tree does not have topological ordering\n");
      break;
    }
  }
}

Int Analyse::run(Symbolic& S) {
  // Perform analyse phase and store the result into the symbolic object S.
  // After Run returns, the Analyse object is not valid.

  if (!ready_) return kRetGeneric;
  if (n_ == 0) return kRetOk;

  HIPO_CLOCK_CREATE;

  HIPO_CLOCK_START(2);
  permute(iperm_);
  eTree();
  postorder();
  HIPO_CLOCK_STOP(2, data_, kTimeAnalyseTree);

  HIPO_CLOCK_START(2);
  colCount();
  HIPO_CLOCK_STOP(2, data_, kTimeAnalyseCount);

  HIPO_CLOCK_START(2);
  fundamentalSupernodes();
  relaxSupernodes();
  afterRelaxSn();
  HIPO_CLOCK_STOP(2, data_, kTimeAnalyseSn);

  HIPO_CLOCK_START(2);
  reorderChildren();
  HIPO_CLOCK_STOP(2, data_, kTimeAnalyseReorder);

  HIPO_CLOCK_START(2);
  snPattern();
  HIPO_CLOCK_STOP(2, data_, kTimeAnalysePattern);

  HIPO_CLOCK_START(2);
  relativeIndCols();
  relativeIndClique();
  HIPO_CLOCK_STOP(2, data_, kTimeAnalyseRelInd);

  computeBlockStart();
  computeCriticalPath();
  computeStackSize();

  computeTreeScheduleSolve();
  computeCriticalPathSolve();

  // move relevant stuff into S
  S.n_ = n_;
  S.sn_ = sn_count_;
  S.nz_ = nz_factor_;
  S.fillin_ = (double)nz_factor_ / nz_;
  S.artificial_nz_ = artificial_nz_;
  S.artificial_ops_ = dense_ops_ - dense_ops_norelax_;
  S.spops_ = sparse_ops_;
  S.critops_ = critical_ops_;
  S.largest_front_ = *std::max_element(sn_indices_.begin(), sn_indices_.end());
  S.serial_storage_ = serial_storage_;
  S.flops_ = dense_ops_;
  S.max_stack_size_ = max_stack_size_;
  S.tree_depth_ = maxDepthTree(sn_parent_);
  S.ops_solve_ = ops_solve_;
  S.critops_solve_ = critical_ops_solve_;

  // compute largest supernode
  std::vector<Int> sn_size(sn_start_.begin() + 1, sn_start_.end());
  for (Int i = sn_count_ - 1; i > 0; --i) sn_size[i] -= sn_size[i - 1];
  S.largest_sn_ = *std::max_element(sn_size.begin(), sn_size.end());

  // build statistics about supernodes size
  for (Int i : sn_size) {
    if (i == 1) S.sn_size_1_++;
    if (i <= 10) S.sn_size_10_++;
    if (i <= 100) S.sn_size_100_++;
  }

  if (checkOverflow()) {
    if (logger_) logger_->printe("Integer overflow in analyse phase\n");
    return kRetIntOverflow;
  }

  // permute signs of pivots
  S.pivot_sign_ = std::move(signs_);
  permuteVector(S.pivot_sign_, perm_);

  S.iperm_ = std::move(iperm_);
  S.rows_ = std::move(rows_sn_);
  S.ptr_ = std::move(ptr_sn_);
  S.sn_parent_ = std::move(sn_parent_);
  S.sn_start_ = std::move(sn_start_);
  S.relind_cols_ = std::move(relind_cols_);
  S.relind_clique_ = std::move(relind_clique_);
  S.consecutive_sums_ = std::move(consecutive_sums_);
  S.clique_block_start_ = std::move(clique_block_start_);
  S.schedule_solve_ = std::move(schedule_solve_);

  S.empty_ = false;

  HIPO_CLOCK_STOP(1, data_, kTimeAnalyse);

  return kRetOk;
}

}  // namespace hipo