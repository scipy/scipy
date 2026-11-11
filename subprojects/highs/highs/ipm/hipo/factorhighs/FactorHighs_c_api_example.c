#include "ipm/hipo/factorhighs/FactorHighs_c_api.h"
#include "math.h"
#include "stdio.h"

/*
  This example factorises the following symmetric matrix

    5
    0     3
    3     2     9
    4     0    -1     8
    0     0     1     0     1

  and solves a linear system, using the FactorHighs C API.
  This file needs to be linked with HiGHS:

  clang FactorHighs_c_api_example.c
  -I/path/to/HiGHS/highs -I/path/to/HiGHS/build
  -L/path/to/HiGHS/build/lib -lhighs
  -Wl,-rpath,/path/to/HiGHS/build/lib

*/

int main() {
  // problem size
  enum { N = 5, NZ = 10, NP1 = N + 1 };

  // define the lower triangle in CSC format, with 0-based indexing
  int ptr[NP1] = {0, 3, 5, 8, 9, 10};
  int rows[NZ] = {0, 2, 3, 1, 2, 2, 3, 4, 3, 4};
  double vals[NZ] = {5, 3, 4, 3, 2, 9, -1, 1, 8, 1};

  /*
  1-based indexing is available as well.
  We define the problem shifted by 1:

    int ptr[n + 1] = {1, 4, 6, 9, 10, 11};
    int rows[nz] = {1, 3, 4, 2, 3, 3, 4, 5, 4, 5};
    double vals[nz] = {5, 3, 4, 3, 2, 9, -1, 1, 8, 1};

  and set the corresponding option:

    FactorHighs_setOneIndexing(FH, 1);
  */

  // matrix is spd, so expect all pivots to be positive
  int signs[N] = {1, 1, 1, 1, 1};

  // rhs and expected solution
  double rhs[N] = {1, 2, 3, 4, 5};
  const double lhs[N] = {0.457627118644068, 1.118644067796610,
                         -0.677966101694915, 0.186440677966102,
                         5.677966101694915};

  // initialise with default number of threads
  const int num_threads = 0;
  int initialise_status = FactorHighs_initialise(num_threads);
  if (initialise_status) {
    printf("HiPO or its dependencies are not available in this build\n");
    return 1;
  }
  void* S = FactorHighs_symbolic_create();
  void* FH = FactorHighs_create();

  // set options
  const int logging_on = 1;
  FactorHighs_setLogging(FH, logging_on);

  const int block_size = 64;
  FactorHighs_setBlockSize(FH, block_size);

  const int pivoting_off = 0;
  FactorHighs_setPivoting(FH, pivoting_off);

  FactorHighs_setRegularisation(FH, 0.0, 0.0);

  // compute ordering with metis
  int perm[N];
  int metis_status = FactorHighs_reorderMetis(FH, N, NZ, rows, ptr, perm);
  if (metis_status) {
    printf("Failed to compute a fill-reducing ordering with Metis\n");
    return 1;
  }

  // perform analyse phase
  int analyse_status =
      FactorHighs_analyse(FH, S, N, NZ, rows, ptr, signs, perm);
  if (analyse_status) {
    printf("Failed to compute symbolic factorisation\n");
    return 1;
  }

  // print extended statistics of symbolic factorisation
  int verbose = 1;
  FactorHighs_symbolic_print(FH, S, verbose);

  // print inverse permutation, that may have been modified by analyse phase
  int iperm[N];
  FactorHighs_iperm(FH, S, iperm);
  printf("\niperm: ");
  for (int i = 0; i < N; ++i) printf("%d ", iperm[i]);
  printf("\n");

  // factorise the matrix
  int factorise_status = FactorHighs_factorise(FH, S, N, NZ, rows, ptr, vals);
  if (factorise_status) {
    printf("Failed to factorise the matrix\n");
    return 1;
  }

  // compute the inertia of the factorisation
  int pos, neg, zero;
  double zero_tolerance = 1e-16;
  FactorHighs_inertia(FH, &pos, &neg, &zero, zero_tolerance);
  printf("\nInertia (+/-/0): %d %d %d\n", pos, neg, zero);

  // triangular solve
  int solve_status = FactorHighs_solve(FH, rhs, 1);
  if (solve_status) {
    printf("Failed to solve the linear system\n");
    return 1;
  }

  /*
  could also perform:

    FactorHighs_forwardSolve(FH, rhs, 1);
    FactorHighs_diagSolve(FH, rhs, 1);
    FactorHighs_backwardSolve(FH, rhs, 1);
  */

  // compute error
  double error = 0.0;
  for (int i = 0; i < N; ++i) error += fabs(rhs[i] - lhs[i]);
  printf("\nError: %e\n", error);

  // terminate
  FactorHighs_symbolic_destroy(S);
  FactorHighs_destroy(FH);
  FactorHighs_terminate();

  return 0;
}
