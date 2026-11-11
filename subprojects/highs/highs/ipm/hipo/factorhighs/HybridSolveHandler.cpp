#include "HybridSolveHandler.h"

#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "FactorHighsSettings.h"
#include "FormatHandler.h"
#include "Swaps.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

HybridSolveHandler::HybridSolveHandler(
    const Symbolic& S, const std::vector<std::vector<double>>& sn_columns,
    const std::vector<std::vector<Int>>& swaps,
    const std::vector<std::vector<char>>& any_swap,
    const std::vector<std::vector<double>>& pivot_2x2, DataCollector& data,
    const FHoptions& options)
    : SolveHandler(S, sn_columns, data, options),
      swaps_{swaps},
      any_swaps_{any_swap},
      pivot_2x2_{pivot_2x2} {
  childrenLinkedList(S_.schedule().task_parent, first_child_, next_child_);

  bool need_parallel_work =
      options_.parallel_backward || options_.parallel_forward;
  bool need_serial_work =
      !options_.parallel_backward || !options_.parallel_forward;

  // allocate workspace for gemv
  if (need_parallel_work) {
    parallel_gemv_workspace_.resize(S_.schedule().count());
    for (Int task = 0; task < S_.schedule().count(); ++task) {
      Int64 largest_front = 0;
      for (Int sn : S_.schedule().sn_per_task[task]) {
        largest_front = std::max(largest_front, S_.ptr(sn + 1) - S_.ptr(sn));
      }
      parallel_gemv_workspace_[task].resize(largest_front);
    }
  }
  if (need_serial_work) {
    serial_gemv_workspace_.resize(S_.largestFront());
  }
}

void HybridSolveHandler::processForwardSn(Int sn, double* x,
                                          std::vector<double>& work, Int task,
                                          Int end_col_in_task) const {
  // Blas calls: dtrsv, dgemv
  // supernode columns in format FH

  HIPO_CLOCK_CREATE;
  const Int nb = options_.nb;

  // if running in parallel, the updates to supernodes that do not belong to
  // this task have to be deferred to avoid data races. Therefore, they are
  // stored in a buffer and applied serially later.
  const bool defer = task >= 0;

  // leading size of supernode
  const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

  // number of columns in the supernode
  const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

  // first colums of the supernode
  const Int sn_start = S_.snStart(sn);

  // index to access S->rows for this supernode
  const Int64 start_row = S_.ptr(sn);

  // number of blocks of columns
  const Int n_blocks = (sn_size - 1) / nb + 1;

  // index to access snColumns[sn]
  Int64 SnCol_ind{};

  if (sn_size < nb) {
    // Fast solve
    // If supernode is small, avoid making BLAS calls
    const Int jb = sn_size;
    const Int x_start = sn_start;

    const Int* current_swaps = nullptr;
    bool any_swaps_in_block = false;

    if (options_.pivoting) {
      HIPO_CLOCK_START(2);
      any_swaps_in_block = any_swaps_[sn][0];
      if (any_swaps_in_block) {
        current_swaps = swaps_[sn].data();
        // apply swaps to portion of rhs that is affected
        permuteWithSwaps(&x[x_start], current_swaps, jb);
      }
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
    }

    HIPO_CLOCK_START(2);
    for (Int row = 0; row < jb; ++row) {
      for (Int col = 0; col < row; ++col) {
        x[x_start + row] -= sn_columns_[sn][col + jb * row] * x[x_start + col];
      }
    }

    for (Int row = jb; row < ldSn; ++row) {
      if (defer) {
        assert(end_col_in_task >= 0);
        const Int row_to_write = S_.rows(start_row + row);
        if (row_to_write < end_col_in_task) {
          for (Int col = 0; col < jb; ++col) {
            x[S_.rows(start_row + row)] -=
                sn_columns_[sn][col + jb * row] * x[x_start + col];
          }
        } else {
          task_rows_[task].push_back(S_.rows(start_row + row));
          task_vals_[task].push_back(0.0);
          for (Int col = 0; col < jb; ++col)
            task_vals_[task].back() +=
                sn_columns_[sn][col + jb * row] * x[x_start + col];
        }
      } else {
        for (Int col = 0; col < jb; ++col) {
          x[S_.rows(start_row + row)] -=
              sn_columns_[sn][col + jb * row] * x[x_start + col];
        }
      }
    }
    HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense);

    if (any_swaps_in_block && options_.pivoting) {
      HIPO_CLOCK_START(2);
      // apply inverse swaps
      permuteWithSwaps(&x[x_start], current_swaps, jb, true);
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
    }

  } else {
    // go through blocks of columns for this supernode
    for (Int j = 0; j < n_blocks; ++j) {
      // number of columns in the block
      const Int jb = std::min(nb, sn_size - nb * j);

      // number of entries in diagonal part
      const Int diag_entries = jb * jb;

      // index to access vector x
      const Int x_start = sn_start + nb * j;

      const Int* current_swaps = nullptr;
      bool any_swaps_in_block = false;

      if (options_.pivoting) {
        HIPO_CLOCK_START(2);
        any_swaps_in_block = any_swaps_[sn][j];
        if (any_swaps_in_block) {
          current_swaps = &swaps_[sn][nb * j];
          // apply swaps to portion of rhs that is affected
          permuteWithSwaps(&x[x_start], current_swaps, jb);
        }
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
      }

      HIPO_CLOCK_START(2);
      callAndTime_dtrsv('U', 'T', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1, data_);

      SnCol_ind += diag_entries;

      // temporary space for gemv
      const Int gemv_size = ldSn - nb * j - jb;
      if (gemv_size > 0) {
        callAndTime_dgemv('T', jb, gemv_size, 1.0, &sn_columns_[sn][SnCol_ind],
                          jb, &x[x_start], 1, 0.0, work.data(), 1, data_);

        SnCol_ind += jb * gemv_size;
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense);

        HIPO_CLOCK_START(2);
        // scatter solution of gemv
        for (Int i = 0; i < gemv_size; ++i) {
          const Int row = S_.rows(start_row + nb * j + jb + i);

          if (defer) {
            assert(end_col_in_task >= 0);
            if (row < end_col_in_task) {
              x[row] -= work[i];
            } else {
              task_rows_[task].push_back(row);
              task_vals_[task].push_back(work[i]);
            }
          } else {
            x[row] -= work[i];
          }
        }
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_sparse);
      }

      if (any_swaps_in_block && options_.pivoting) {
        HIPO_CLOCK_START(2);
        // apply inverse swaps
        permuteWithSwaps(&x[x_start], current_swaps, jb, true);
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
      }
    }
  }
}

void HybridSolveHandler::processForwardTask(Int task, double* x) const {
  // wait for children to complete
  highs::parallel::TaskGroup tg;
  Int child = first_child_[task];
  while (child != -1) {
    tg.spawn([=, &x]() { processForwardTask(child, x); });
    child = next_child_[child];
  }
  tg.taskWait();

  task_rows_[task].clear();
  task_vals_[task].clear();

  const Int lead_sn = S_.schedule().sn_per_task[task].back();
  Int end_col_in_task = S_.snStart(lead_sn + 1);

  // assembel contributions of children
  child = first_child_[task];
  while (child != -1) {
    for (Int i = 0; i < static_cast<Int>(task_rows_[child].size()); ++i) {
      if (task_rows_[child][i] < end_col_in_task)
        x[task_rows_[child][i]] -= task_vals_[child][i];
      else {
        task_rows_[task].push_back(task_rows_[child][i]);
        task_vals_[task].push_back(task_vals_[child][i]);
      }
    }
    child = next_child_[child];
  }

  for (Int sn : S_.schedule().sn_per_task[task]) {
    processForwardSn(sn, x, parallel_gemv_workspace_[task], task,
                     end_col_in_task);
  }
}

void HybridSolveHandler::forwardSolve(double* x) const {
  if (options_.parallel_forward && S_.schedule().valid) {
    // Hard to parallelise: a sn depends on its children in the tree; multiple
    // children may be writing to the same location in x at the same time.
    // Special care is needed for the writes, involving private buffers.
    task_rows_.resize(S_.schedule().count());
    task_vals_.resize(S_.schedule().count());

    highs::parallel::TaskGroup tg;
    for (Int task = 0; task < S_.schedule().count(); ++task) {
      if (S_.schedule().task_parent[task] == -1) {
        tg.spawn([=, &x]() { processForwardTask(task, x); });
      }
    }
    tg.taskWait();

  } else {
    for (Int sn = 0; sn < S_.sn(); ++sn) {
      processForwardSn(sn, x, serial_gemv_workspace_);
    }
  }
}

void HybridSolveHandler::processBackwardSn(Int sn, double* x,
                                           std::vector<double>& work) const {
  // Blas calls: dtrsv, dgemv
  // supernode columns in format FH

  HIPO_CLOCK_CREATE;
  const Int nb = options_.nb;

  // leading size of supernode
  const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

  // number of columns in the supernode
  const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

  // first colums of the supernode
  const Int sn_start = S_.snStart(sn);

  // index to access S->rows for this supernode
  const Int64 start_row = S_.ptr(sn);

  // number of blocks of columns
  const Int n_blocks = (sn_size - 1) / nb + 1;

  // index to access snColumns[sn]
  // initialised with the total number of entries of snColumns[sn]
  Int64 SnCol_ind = sn_columns_[sn].size() - extra_space_frontal;

  if (sn_size < nb) {
    // Fast solve
    // If supernode is small, avoid making BLAS calls
    const Int jb = sn_size;
    const Int x_start = sn_start;

    const Int* current_swaps = nullptr;
    bool any_swaps_in_block = false;

    if (options_.pivoting) {
      HIPO_CLOCK_START(2);
      any_swaps_in_block = any_swaps_[sn][0];
      if (any_swaps_in_block) {
        current_swaps = swaps_[sn].data();
        // apply swaps to portion of rhs that is affected
        permuteWithSwaps(&x[x_start], current_swaps, jb);
      }
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
    }

    HIPO_CLOCK_START(2);
    for (Int row = ldSn - 1; row >= jb; --row) {
      for (Int col = jb - 1; col >= 0; --col) {
        x[x_start + col] -=
            sn_columns_[sn][col + row * jb] * x[S_.rows(start_row + row)];
      }
    }

    for (Int row = jb - 1; row >= 0; --row) {
      for (Int col = row - 1; col >= 0; --col) {
        x[x_start + col] -= sn_columns_[sn][col + row * jb] * x[x_start + row];
      }
    }
    HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense);

    if (any_swaps_in_block && options_.pivoting) {
      HIPO_CLOCK_START(2);
      // apply inverse swaps
      permuteWithSwaps(&x[x_start], current_swaps, jb, true);
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
    }

  } else {
    // go through blocks of columns for this supernode in reverse order
    for (Int j = n_blocks - 1; j >= 0; --j) {
      // number of columns in the block
      const Int jb = std::min(nb, sn_size - nb * j);

      // number of entries in diagonal part
      const Int diag_entries = jb * jb;

      // index to access vector x
      const Int x_start = sn_start + nb * j;

      const Int* current_swaps = nullptr;
      bool any_swaps_in_block = false;

      if (options_.pivoting) {
        HIPO_CLOCK_START(2);
        any_swaps_in_block = any_swaps_[sn][j];
        if (any_swaps_in_block) {
          current_swaps = &swaps_[sn][nb * j];
          // apply swaps to portion of rhs that is affected
          permuteWithSwaps(&x[x_start], current_swaps, jb);
        }
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
      }

      // temporary space for gemv
      const Int gemv_size = ldSn - nb * j - jb;
      if (gemv_size > 0) {
        HIPO_CLOCK_START(2);
        // scatter entries into y
        for (Int i = 0; i < gemv_size; ++i) {
          const Int row = S_.rows(start_row + nb * j + jb + i);
          work[i] = x[row];
        }
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_sparse);

        HIPO_CLOCK_START(2);
        SnCol_ind -= jb * gemv_size;
        callAndTime_dgemv('N', jb, gemv_size, -1.0, &sn_columns_[sn][SnCol_ind],
                          jb, work.data(), 1, 1.0, &x[x_start], 1, data_);
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense);
      }

      HIPO_CLOCK_START(2);
      SnCol_ind -= diag_entries;
      callAndTime_dtrsv('U', 'N', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1, data_);
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense);

      if (any_swaps_in_block && options_.pivoting) {
        HIPO_CLOCK_START(2);
        // apply inverse swaps
        permuteWithSwaps(&x[x_start], current_swaps, jb, true);
        HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
      }
    }
  }
}

void HybridSolveHandler::processBackwardTask(Int task, double* x) const {
  for (auto rit = S_.schedule().sn_per_task[task].rbegin();
       rit != S_.schedule().sn_per_task[task].rend(); ++rit) {
    const Int sn = *rit;
    processBackwardSn(sn, x, parallel_gemv_workspace_[task]);
  }

  // wait for children to complete
  highs::parallel::TaskGroup tg;
  Int child = first_child_[task];
  while (child != -1) {
    tg.spawn([=, &x]() { processBackwardTask(child, x); });
    child = next_child_[child];
  }
  tg.taskWait();
}

void HybridSolveHandler::backwardSolve(double* x) const {
  if (options_.parallel_backward && S_.schedule().valid) {
    // Easy to parallelise: a sn depends on its ancestors in the tree; the
    // ancestor is the only sn running in a given branch when it writes the
    // update, so no special care needs to be taken for the writes. Respecting
    // the dependencies of the tree is enough.
    highs::parallel::TaskGroup tg;
    for (Int task = 0; task < S_.schedule().count(); ++task) {
      if (S_.schedule().task_parent[task] == -1) {
        tg.spawn([=, &x]() { processBackwardTask(task, x); });
      }
    }
    tg.taskWait();

  } else {
    for (Int sn = S_.sn() - 1; sn >= 0; --sn) {
      processBackwardSn(sn, x, serial_gemv_workspace_);
    }
  }
}

void HybridSolveHandler::processDiagSn(Int sn, double* x) const {
  // supernode columns in format FH

  HIPO_CLOCK_CREATE;
  const Int nb = options_.nb;

  // leading size of supernode
  const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

  // number of columns in the supernode
  const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

  // first colums of the supernode
  const Int sn_start = S_.snStart(sn);

  // number of blocks of columns
  const Int n_blocks = (sn_size - 1) / nb + 1;

  // index to access diagonal part of block
  Int diag_start{};

  // go through blocks of columns for this supernode
  for (Int j = 0; j < n_blocks; ++j) {
    // number of columns in the block
    const Int jb = std::min(nb, sn_size - nb * j);

    const Int* current_swaps = nullptr;
    bool any_swaps_in_block = false;

    if (options_.pivoting) {
      HIPO_CLOCK_START(2);
      any_swaps_in_block = any_swaps_[sn][j];
      if (any_swaps_in_block) {
        current_swaps = &swaps_[sn][nb * j];
        // apply swaps to portion of rhs that is affected
        permuteWithSwaps(&x[sn_start + nb * j], current_swaps, jb);
      }
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
    }

    HIPO_CLOCK_START(2);

    const double* current_2x2 = &pivot_2x2_[sn][nb * j];
    Int step = 1;

    // go through columns of block
    for (Int col = 0; col < jb; col += step) {
      if (current_2x2[col] == 0.0) {
        // 1x1 pivots
        step = 1;
        const double inv_d = sn_columns_[sn][diag_start + col + jb * col];
        x[sn_start + nb * j + col] *= inv_d;
      } else {
        // 2x2 pivots
        step = 2;

        // inverse of 2x2 pivot
        const double i_d1 = sn_columns_[sn][diag_start + col + jb * col];
        const double i_d2 =
            sn_columns_[sn][diag_start + col + 1 + jb * (col + 1)];
        const double i_off = current_2x2[col];

        double x1 = x[sn_start + nb * j + col];
        double x2 = x[sn_start + nb * j + col + 1];

        x[sn_start + nb * j + col] = i_d1 * x1 + i_off * x2;
        x[sn_start + nb * j + col + 1] = i_d2 * x2 + i_off * x1;
      }
    }

    HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_dense);

    if (any_swaps_in_block && options_.pivoting) {
      HIPO_CLOCK_START(2);
      // apply inverse swaps
      permuteWithSwaps(&x[sn_start + nb * j], current_swaps, jb, true);
      HIPO_CLOCK_STOP(2, data_, kTimeSolveSolve_swap);
    }

    // move diag_start forward by number of diagonal entries in block
    diag_start += jb * jb;

    // move diag_start forward by number of sub-diagonal entries in block
    diag_start += (ldSn - nb * j - jb) * jb;
  }
}

void HybridSolveHandler::diagSolve(double* x) const {
  if (options_.parallel_diag) {
    // Embarassingly parallel: each sn reads/writes independent entries in x.
    highs::parallel::for_each(
        0, S_.sn(),
        [&](Int start, Int end) {
          for (Int sn = start; sn < end; ++sn) {
            processDiagSn(sn, x);
          }
        },
        // choose grainsize so that the number of tasks in the for_each loop
        // that execute the function is roughly kParallelDiagTargetNumTasks
        std::ceil((double)S_.size() / kParallelDiagTargetNumTasks));

  } else {
    for (Int sn = 0; sn < S_.sn(); ++sn) {
      processDiagSn(sn, x);
    }
  }
}

void HybridSolveHandler::inertia(Int& pos, Int& neg, Int& zero,
                                 double tol) const {
  pos = 0;
  neg = 0;
  zero = 0;

  const Int nb = options_.nb;

  for (Int sn = 0; sn < S_.sn(); ++sn) {
    // leading size of supernode
    const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    // number of columns in the supernode
    const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

    // number of blocks of columns
    const Int n_blocks = (sn_size - 1) / nb + 1;

    // index to access diagonal part of block
    Int diag_start{};

    // go through blocks of columns for this supernode
    for (Int j = 0; j < n_blocks; ++j) {
      // number of columns in the block
      const Int jb = std::min(nb, sn_size - nb * j);
      const double* current_2x2 = &pivot_2x2_[sn][nb * j];
      Int step = 1;

      // go through columns of block
      for (Int col = 0; col < jb; col += step) {
        if (current_2x2[col] == 0.0) {
          // 1x1 pivots
          step = 1;
          const double inv_d = sn_columns_[sn][diag_start + col + jb * col];
          const double pivot = 1.0 / inv_d;

          if (pivot > tol)
            ++pos;
          else if (pivot < -tol)
            ++neg;
          else
            ++zero;

        } else {
          // 2x2 pivots
          step = 2;

          // inverse of 2x2 pivot
          const double i_d1 = sn_columns_[sn][diag_start + col + jb * col];
          const double i_d2 =
              sn_columns_[sn][diag_start + col + 1 + jb * (col + 1)];
          const double i_off = current_2x2[col];

          // determinant and trace of inverse
          const double i_det = i_d1 * i_d2 - i_off * i_off;
          const double i_trace = i_d1 + i_d2;

          // determinant and trace of 2x2 pivot
          const double det = 1.0 / i_det;
          const double trace = det * i_trace;

          if (std::abs(det) < tol) {
            // det is zero, so at least one pivot is zero
            ++zero;

            if (trace > tol)
              ++pos;
            else if (trace < -tol)
              ++neg;
            else
              ++zero;

          } else if (det > tol) {
            // det is positive, so pivots have same sign
            if (trace > 0)
              pos += 2;
            else
              neg -= 2;
          } else {
            // indefinite 2x2 pivot
            ++pos;
            ++neg;
          }
        }
      }

      // move diag_start forward by number of diagonal entries in block
      diag_start += jb * jb;

      // move diag_start forward by number of sub-diagonal entries in block
      diag_start += (ldSn - nb * j - jb) * jb;
    }
  }
}

}  // namespace hipo
