#ifndef FACTOR_HIGHS_H
#define FACTOR_HIGHS_H

#include "CliqueStack.h"
#include "DataCollector.h"
#include "FactorHighsOptions.h"
#include "Numeric.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/Logger.h"

/*

Direct solver for IPM matrices.

Consider a sparse symmetric matrix M in CSC format.
Only its lower triangular part is used; entries in the upper triangle are
ignored. Diagonal entries must be present for each column, even if the entry is
zero.
The matrix has n rows/cols and nz nonzero entries. It is stored using
three arrays:
- ptr, column pointers, of length n+1;
- rows, row indices, of length nz;
- vals, values, of length nz.

The direct solver uses the following objects:
- Symbolic, to store the symbolic factorization;
- FHsolver, to perform analyse and factorise phases.

Define a right-hand side rhs, which will be overwritten with the solution of
M^{-1} * rhs. Use reorderMetis, reorderAmd, reorderRcm to obtain a fill-reducing
ordering of the matrix, or use a precomputed one. The ordering is stored in the
array perm.

Zero-based indexing is assumed for the sparsity pattern and permutation. Use
setOneIndexing(true) to use one-based indexing.

Define a vector signs that contains the expected sign of each pivot:
-  1 for pivots expected to be positive
- -1 for pivots expected to be negative
-  0 for pivots without an expected sign.
This is used to determine the sign of the regularisation to apply.

Then, the factorization is performed as follows.

    Symbolic S;
    FHsolver FH;
    FH.analyse(S, n, nz rows, ptr, signs, perm);
    FH.factorise(S, n, nz rows, ptr, val);
    FH.solve(x);

Printing to screen is achieved using the interface in auxiliary/Logger.h.
Use setLogger to pass the Logger object to use: FH.setLogger(&logger).
To use printf instead, use FH.setLogger(nullptr,true).
Logging is off by default. Use FH.setLogger(nullptr, false) to switch logging
off.

To add static regularisation when the pivots are selected, use
setRegularisation(reg_p,reg_d) to choose values of primal and dual
regularisation. If regularisation is already added to the matrix, ignore.

The expected sign of a pivot is used to determine which static regularisation to
apply. Pivots with negative expected sign are regularised with -reg_p, pivots
with positive expected sign are regularised with +reg_d. If the sign is unknown,
the computed sign of the pivot is used to determine which regularisation to
apply.

Notice that this sign convention corresponds to a saddle point system with
negative (1,1) block and positive (2,2) block.

Dynamic regularisation may perturb a pivot if it is too small, and it
does so based on the expected sign of the pivot. This procedure uses the
computed sign of the pivot, if the sign is unknown.

Notice that the fill-reducing ordering can be modified during the call to
analyse. The inverse permutation used during the factorisation can be accessed
via the Symbolic object, S.iperm().

*/

namespace hipo {

class FHsolver {
  const Logger* logger_;
  DataCollector data_;
  Numeric N_;
  CliqueStack serial_stack_;
  bool local_logger_ = false;

  FHoptions options_{};

  // Columns of factorisation, stored by supernode.
  // This memory is allocated the first time that it is used. Subsequent
  // factorisations reuse the same memory.
  std::vector<std::vector<double>> sn_columns_;

 public:
  // Print collected data (if any) and terminate DataCollector
  ~FHsolver();

  // Perform analyse phase of matrix with sparsity pattern given by rows and
  // ptr, and store symbolic factorisation in object S.
  // See ReturnValues.h for errors.
  Int analyse(Symbolic& S, Int n, Int nz, const Int* rows, const Int* ptr,
              const Int* signs, const Int* perm);

  // Perform factorise phase of matrix given by rows, ptr, vals, and store
  // numerical factorisation in object N. See ReturnValues.h for errors.
  Int factorise(const Symbolic& S, Int n, Int nz, const Int* rows,
                const Int* ptr, const double* vals);

  // Perform solve phase with rhs given by x, which is overwritten with the
  // solution. Multiple rhs are supported by doing the solves in parallel.
  Int solve(double* x, Int k = 1) const;

  // Perform partial solves with L, D, L^T (including permutation)
  Int forwardSolve(double* x, Int k = 1) const;
  Int diagSolve(double* x, Int k = 1) const;
  Int backwardSolve(double* x, Int k = 1) const;

  // If multiple factorisation are performed, call newIter() before each
  // factorisation. This is used only to collect data for debugging, if
  // expensive data collection is turned on at compile time.
  void newIter();

  // Set values for static regularisation to be added when a pivot is selected.
  // If regularisation is already added to the matrix, ignore.
  void setRegularisation(double reg_p = 0.0, double reg_d = 0.0);

  // Extract the regularisation values used, including static and dynamic.
  void getRegularisation(double* reg);

  // Set block size for dense linear algebra
  void setBlockSize(Int nb);

  // Set pivoting optins, on by default.
  // It uses a static variation of Bunch-Kaufman pivoting, with potential
  // dynamic regularisation. If pivoting is switched off, only static
  // regularisation is applied.
  void setPivoting(bool pivoting);

  void setParallel(bool tree, bool node);
  void setParallelSolve(bool forward, bool backward, bool diag);

  // Pass the Logger object to be used for logging. Alternatively, printf can be
  // used for logging, by passing a nullptr and setting use_printf to true.
  // By default, logging is off.
  void setLogger(const Logger* logger = nullptr, bool use_printf = false);
  const Logger* getLogger() const { return logger_; }

  // Compute number of positive, negative and zero pivots, using tol as
  // tolerance for zero.
  void inertia(Int& pos, Int& neg, Int& zero, double tol = 1e-16) const;

  // Set options for 1-based indexing, off by default.
  // If set to true, the vectors passed to analyse and factorise are assumed to
  // use 1-based indexing, so that all entries are shifted down by 1.
  void setOneIndexing(bool one_indexing);

  void iperm(const Symbolic& S, Int* ip) const;

  // Compute fill-reducing permutation using nested dissection (Metis), approx
  // minimum degree (amd) or reverse Cuthill-McGee (rcm).
  // Set full_matrix_0 = false if the sparsity pattern in rows and ptr
  // corresponds to the lower triangle of the matrix.
  // Set full_matrix_0 = true if the sparsity pattern in rows and ptr
  // corresponds to the full matrix (without the diagonal entries), with 0-based
  // indexing.
  Int reorderMetis(Int n, Int nz, const Int* rows, const Int* ptr, Int* perm,
                   bool full_matrix_0, Int seed) const;
  Int reorderAmd(Int n, Int nz, const Int* rows, const Int* ptr, Int* perm,
                 bool full_matrix_0) const;
  Int reorderRcm(Int n, Int nz, const Int* rows, const Int* ptr, Int* perm,
                 bool full_matrix_0) const;
};

}  // namespace hipo

#endif