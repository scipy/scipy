#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
const double inf = HIGHS_CONST_INF;
void reportIssue(const int issue) {
  printf("\n *************\n * Issue %3d *\n *************\n", issue);
}
bool objectiveOk(const double optimal_objective,
                 const double require_optimal_objective) {
  double error = std::fabs(optimal_objective - require_optimal_objective) /
                 std::max(1.0, std::fabs(require_optimal_objective));
  bool error_ok = error < 1e-10;
  if (!error_ok)
    printf("Objective is %g but require %g (error %g)\n", optimal_objective,
           require_optimal_objective, error);
  return error_ok;
}

void reportSolution(Highs& highs) {
  const HighsInfo& info = highs.getHighsInfo();
  if (info.primal_status == PrimalDualStatus::STATUS_FEASIBLE_POINT) {
    const HighsSolution& solution = highs.getSolution();
    printf("Solution\n");
    printf("Col       Value        Dual\n");
    for (int iCol = 0; iCol < highs.getLp().numCol_; iCol++)
      printf("%3d %11.4g %11.4g\n", iCol, solution.col_value[iCol],
             solution.col_dual[iCol]);
    printf("Row       Value        Dual\n");
    for (int iRow = 0; iRow < highs.getLp().numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, solution.row_value[iRow],
             solution.row_dual[iRow]);
  } else {
    printf("info.primal_status = %d\n", info.primal_status);
  }
}

void solve(Highs& highs, std::string presolve, std::string solver,
           const HighsModelStatus require_model_status,
           const double require_optimal_objective = 0) {
  const HighsInfo& info = highs.getHighsInfo();
  HighsStatus status;
  HighsModelStatus model_status;

  status = highs.setHighsOptionValue("solver", solver);
  REQUIRE(status == HighsStatus::OK);

  status = highs.setHighsOptionValue("presolve", presolve);
  REQUIRE(status == HighsStatus::OK);

  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    REQUIRE(
        objectiveOk(info.objective_function_value, require_optimal_objective));
  }

  status = highs.resetHighsOptions();
  REQUIRE(status == HighsStatus::OK);
}

void issue272(Highs& highs) {
  reportIssue(272);
  // This is the FuOR MIP that presolve failed to handle as a maximization
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double optimal_objective = 8.83333333333333;
  lp.numCol_ = 2;
  lp.numRow_ = 2;
  lp.colCost_ = {3, 2};
  lp.colLower_ = {0, 0};
  lp.colUpper_ = {inf, inf};
  lp.rowLower_ = {-inf, -inf};
  lp.rowUpper_ = {23, 10};
  lp.Astart_ = {0, 2, 4};
  lp.Aindex_ = {0, 1, 0, 1};
  lp.Avalue_ = {3, 5, 6, 2};
  lp.sense_ = ObjSense::MAXIMIZE;

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  // Presolve reduces to empty, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue280(Highs& highs) {
  reportIssue(280);
  // This is an easy problem from mckib2 that IPX STILL FAILS to handle
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double optimal_objective = 1;
  lp.numCol_ = 2;
  lp.numRow_ = 2;
  lp.colCost_ = {-1, 1};
  lp.colLower_ = {1, 2};
  lp.colUpper_ = {1, 2};
  lp.rowLower_ = {-inf, 2};
  lp.rowUpper_ = {1, 2};
  lp.Astart_ = {0, 1, 2};
  lp.Aindex_ = {0, 1};
  lp.Avalue_ = {1, 1};

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  // Presolve reduces to empty, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  reportSolution(highs);
  // STILL FAILS!!! Reported to Lukas as issue #1 on IPX
  //    solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue282(Highs& highs) {
  reportIssue(282);
  // This is an easy problem from mckib2 on which presolve+simplex
  // failed to give the correct objective
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double optimal_objective = -18;
  lp.numCol_ = 2;
  lp.numRow_ = 3;
  lp.colCost_ = {-3, -2};
  lp.colLower_ = {0, 0};
  lp.colUpper_ = {inf, inf};
  lp.rowLower_ = {-inf, -inf, -inf};
  lp.rowUpper_ = {10, 8, 4};
  lp.Astart_ = {0, 3, 5};
  lp.Aindex_ = {0, 1, 2, 0, 1};
  lp.Avalue_ = {2, 1, 1, 1, 1};

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  // Presolve reduces to empty, so no real need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue285(Highs& highs) {
  reportIssue(285);
  // This is an infeasible LP for which HiGHS segfaulted after "Problem
  // status detected on presolve: Infeasible"
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status =
      HighsModelStatus::PRIMAL_INFEASIBLE;
  lp.numCol_ = 2;
  lp.numRow_ = 3;
  lp.colCost_ = {-4, 1};
  lp.colLower_ = {2, 0};
  lp.colUpper_ = {2, inf};
  lp.rowLower_ = {-inf, -inf, -inf};
  lp.rowUpper_ = {14, 0, 3};
  lp.Astart_ = {0, 2, 5};
  lp.Aindex_ = {0, 2, 0, 1, 2};
  lp.Avalue_ = {7, 2, -2, 1, -2};

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  // Presolve identifies infeasibility, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
}

void issue295(Highs& highs) {
  reportIssue(295);
  // Simplex solver (without presolve) gets a correct solution, IPX
  // (without presolve) reported the correct objective function value
  // but an inconsistent solution. Both simplex and IPX reported an
  // error when presolve is on.
  //
  // The bug in presolve was due to relying on size to get the number
  // of nonzeros. Correct interpretations of IPX inconsistent solution
  // was added
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double optimal_objective = -2;
  lp.numCol_ = 5;
  lp.numRow_ = 2;
  lp.colCost_ = {0, 0, 0, 1, -1};
  lp.colLower_ = {-inf, -inf, -inf, -1, -1};
  lp.colUpper_ = {inf, inf, inf, 1, 1};
  lp.rowLower_ = {-inf, -inf};
  lp.rowUpper_ = {2, -2};
  lp.Astart_ = {0, 1, 2, 3, 3, 3};
  lp.Aindex_ = {0, 1, 0};
  lp.Avalue_ = {1, 1, 0};

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue306(Highs& highs) {
  reportIssue(306);
  // This is test6690 from mckib2 that gave a small inconsistency in
  // a bound after presolve, causing an error in IPX
  //
  // Resulted in the introduction of cleanBounds after presolve
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double optimal_objective = -1.191;
  lp.numCol_ = 10;
  lp.numRow_ = 6;
  lp.colCost_ = {-1.64, 0.7, 1.8, -1.06, -1.16, 0.26, 2.13, 1.53, 0.66, 0.28};
  lp.colLower_ = {-0.84, -0.97, 0.34, 0.4,   -0.33,
                  -0.74, 0.47,  0.09, -1.45, -0.73};
  lp.colUpper_ = {0.37, 0.02, 2.86, 0.86, 1.18, 0.5, 1.76, 0.17, 0.32, -0.15};
  lp.rowLower_ = {0.9626, -1e+200, -1e+200, -1e+200, -1e+200, -1e+200};
  lp.rowUpper_ = {0.9626, 0.615, 0, 0.172, -0.869, -0.022};
  lp.Astart_ = {0, 0, 1, 2, 5, 5, 6, 7, 9, 10, 12};
  lp.Aindex_ = {4, 4, 0, 1, 3, 0, 4, 1, 5, 0, 1, 4};
  lp.Avalue_ = {-1.22, -0.25, 0.93,  1.18, 0.43,  0.65,
                -2.06, -0.2,  -0.25, 0.83, -0.22, 1.37};

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue316(Highs& highs) {
  reportIssue(316);
  // This is a test problem from matbesancon where maximization failed
  //
  // Resulted in fixes being added to unconstrained LP solver
  HighsStatus status;
  bool bool_status;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double min_optimal_objective = -6;
  const double max_optimal_objective = 12;
  status = highs.clearModel();
  REQUIRE(status == HighsStatus::OK);

  bool_status = highs.addCol(2, -3, 6, 0, NULL, NULL);
  REQUIRE(bool_status);

  // Presolve reduces to empty
  solve(highs, "on", "simplex", require_model_status, min_optimal_objective);
  solve(highs, "off", "simplex", require_model_status, min_optimal_objective);

  bool_status = highs.changeObjectiveSense(ObjSense::MAXIMIZE);
  REQUIRE(bool_status);

  solve(highs, "on", "simplex", require_model_status, max_optimal_objective);
  solve(highs, "off", "simplex", require_model_status, max_optimal_objective);
}

void mpsUnbounded(Highs& highs) {
  // As a maximization, adlittle is unbounded, but a bug in hsol [due
  // to jumping to phase 2 if optimal in phase 1 after clean-up
  // yielded no dual infeasiblities despite the phase 1 objective
  // being negative] resulted in the problem being declared infeasible
  //
  // Resulted in fixes being added to hsol dual
  HighsStatus status;
  bool bool_status;
  const HighsModelStatus require_model_status =
      HighsModelStatus::PRIMAL_UNBOUNDED;

  // Unit test fails for IPX with adlittle_max
  const bool solve_adlittle_max = false;
  std::string model = "adlittle";
  if (!solve_adlittle_max) model = "gas11";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  status = highs.readModel(model_file);
  REQUIRE(status == HighsStatus::OK);

  if (solve_adlittle_max) {
    bool_status = highs.changeObjectiveSense(ObjSense::MAXIMIZE);
    REQUIRE(bool_status);
  }
  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
}

void almostNotUnbounded(Highs& highs) {
  // This problem tests how well HiGHS handles
  // near-unboundedness. None of the LPs is reduced by presolve
  //
  // No
  HighsStatus status;
  HighsLp lp;
  const HighsModelStatus require_model_status0 =
      HighsModelStatus::PRIMAL_UNBOUNDED;
  const HighsModelStatus require_model_status1 = HighsModelStatus::OPTIMAL;
  const HighsModelStatus require_model_status2 = HighsModelStatus::OPTIMAL;
  const double optimal_objective1 = -1;
  const double optimal_objective2 = -3;

  // With epsilon = 1e-7 hsol declares optimality as the primal
  // infeasibility in phase 2 that would lead to unboundedness is
  // sufficiently small
  //
  // With epsilon = 1e-6 hsol goes to phase 2 before detecting
  // unboundedness, because the problem with perturbed costs is dual
  // feasible.
  //
  // With epsilon = 1e-5 hsol identifies unboundedness in phase 1
  // because the problem with perturbed costs is not dual feasible.
  double epsilon = 1e-6;
  lp.numCol_ = 2;
  lp.numRow_ = 3;
  lp.colCost_ = {-1, 1 - epsilon};
  lp.colLower_ = {0, 0};
  lp.colUpper_ = {1e+200, 1e+200};
  lp.rowLower_ = {-1 + epsilon, -1, 3};
  lp.rowUpper_ = {1e+200, 1e+200, 1e+200};
  lp.Astart_ = {0, 3, 6};
  lp.Aindex_ = {0, 1, 2, 0, 1, 2};
  lp.Avalue_ = {1 + epsilon, -1, 1, -1, 1, 1};
  // LP is feasible on [1+alpha, alpha] with objective
  // -1-epsilon*alpha so unbounded

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  //  status = highs.writeModel("epsilon_unbounded.mps"); REQUIRE(status ==
  //  HighsStatus::WARNING);
  solve(highs, "off", "simplex", require_model_status0);
  solve(highs, "off", "ipm", require_model_status0);

  // LP is feasible on [1+alpha, alpha] with objective -1 so optimal,
  // but has open set of optimal solutions
  lp.colCost_ = {-1, 1};
  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);

  solve(highs, "off", "simplex", require_model_status1, optimal_objective1);
  reportSolution(highs);
  solve(highs, "off", "ipm", require_model_status1, optimal_objective1);

  // LP has bounded feasible region with optimal solution
  // [1+2/epsilon, 2/epsilon] and objective
  // value -3
  lp.colCost_[1] = 1 - epsilon;
  lp.rowLower_[0] = -1 - epsilon;
  lp.Avalue_[0] = 1 - epsilon;
  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);

  solve(highs, "off", "simplex", require_model_status2, optimal_objective2);
  reportSolution(highs);
  solve(highs, "off", "ipm", require_model_status2, optimal_objective2);
}
TEST_CASE("test-special-lps", "[TestSpecialLps]") {
  Highs highs;
  issue272(highs);
  issue280(highs);
  issue282(highs);
  issue285(highs);
  issue295(highs);
  issue306(highs);
  issue316(highs);
  mpsUnbounded(highs);
  almostNotUnbounded(highs);
}
