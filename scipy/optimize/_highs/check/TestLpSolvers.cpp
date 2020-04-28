#include "Highs.h"
#include "catch.hpp"

struct IterationCount {
  int simplex;
  int ipm;
  int crossover;
};

void testSolver(Highs& highs, const std::string solver,
                IterationCount& default_iteration_count,
                const int int_simplex_strategy = 0) {
  double default_time_limit;
  int default_simplex_iteration_limit;
  int default_ipm_iteration_limit;
  // Initialise use_simplex_iteration and use_ipm_iteration_limit to
  // stop compiler warning about them being uninitialized
  int use_simplex_iteration_limit = 0;
  int use_ipm_iteration_limit = 0;
  int simplex_iteration_count;
  int ipm_iteration_count;
  int crossover_iteration_count;
  HighsModelStatus model_status;
  HighsStatus return_status;
  const bool perform_timeout_test = false;  // true;  //
  const bool use_simplex = solver == "simplex";

  const HighsInfo& info = highs.getHighsInfo();

  return_status = highs.setHighsOptionValue("solver", solver);
  REQUIRE(return_status == HighsStatus::OK);

  if (use_simplex) {
    SimplexStrategy simplex_strategy =
        static_cast<SimplexStrategy>(int_simplex_strategy);
    if (simplex_strategy == SimplexStrategy::SIMPLEX_STRATEGY_DUAL_TASKS)
      return;
    printf("Simplex strategy %d\n", int_simplex_strategy);
    return_status =
        highs.setHighsOptionValue("simplex_strategy", simplex_strategy);
    REQUIRE(return_status == HighsStatus::OK);
  }

  return_status = highs.getHighsOptionValue("time_limit", default_time_limit);
  REQUIRE(return_status == HighsStatus::OK);

  if (use_simplex) {
    simplex_iteration_count = info.simplex_iteration_count;
    return_status = highs.getHighsOptionValue("simplex_iteration_limit",
                                              default_simplex_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
    // Force HiGHS to start from a logical basis - if this is the
    // second or subsequent call to testSolver
    return_status = highs.setBasis();
    REQUIRE(return_status == HighsStatus::OK);
  } else {
    ipm_iteration_count = info.ipm_iteration_count;
    return_status = highs.getHighsOptionValue("ipm_iteration_limit",
                                              default_ipm_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
    crossover_iteration_count = info.crossover_iteration_count;
  }

  // Vanilla solve: get solution time to calibrate time limit test
  double run_time = highs.getHighsRunTime();
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);
  const double single_solve_run_time = highs.getHighsRunTime() - run_time;

  if (use_simplex) {
    simplex_iteration_count += default_iteration_count.simplex;
    REQUIRE(info.simplex_iteration_count == simplex_iteration_count);
  } else {
    printf("IPM: %d; Crossover: %d\n", info.ipm_iteration_count,
           info.crossover_iteration_count);
    ipm_iteration_count += default_iteration_count.ipm;
    REQUIRE(info.ipm_iteration_count == ipm_iteration_count);
    crossover_iteration_count += default_iteration_count.crossover;
    REQUIRE(info.crossover_iteration_count == crossover_iteration_count);
  }

  // Only perform the time limit test if the solve time is large enough
  const double min_run_time_for_test = 0.001;
  if (perform_timeout_test && single_solve_run_time > min_run_time_for_test) {
    const int ideal_num_solve = 10;
    const double local_time_limit = ideal_num_solve * single_solve_run_time;

    // Solve with time limit
    run_time = highs.getHighsRunTime();
    printf("Current run time is %g\n", run_time);

    double use_time_limit = run_time + local_time_limit;
    return_status = highs.setHighsOptionValue("time_limit", use_time_limit);
    REQUIRE(return_status == HighsStatus::OK);

    const int max_num_solve = 10 * ideal_num_solve;
    int num_solve;
    for (num_solve = 0; num_solve < max_num_solve; num_solve++) {
      if (use_simplex) return_status = highs.setBasis();
      return_status = highs.run();
      if (highs.getModelStatus() == HighsModelStatus::REACHED_TIME_LIMIT) break;
    }
    REQUIRE(num_solve < max_num_solve);
    run_time = highs.getHighsRunTime();
    printf("Current run time is %g: time limit is %g (difference = %g)\n",
           run_time, use_time_limit, run_time - use_time_limit);

    printf("Required %d solves (ideally %d - max %d)\n", num_solve,
           ideal_num_solve, max_num_solve);
  } else {
    printf(
        "Not performed the time limit test since solve time is %g <= %g = "
        "min_run_time_for_test\n",
        single_solve_run_time, min_run_time_for_test);
  }
  return_status = highs.setHighsOptionValue("time_limit", default_time_limit);
  REQUIRE(return_status == HighsStatus::OK);
  if (!use_simplex)
    printf("IPM: %d; Crossover: %d\n", info.ipm_iteration_count,
           info.crossover_iteration_count);

  // Solve with iteration limit
  // First of all check that no iterations are performed if the
  // iteration limit is the current number of iterations
  if (use_simplex) {
    simplex_iteration_count = info.simplex_iteration_count;
    use_simplex_iteration_limit = simplex_iteration_count;

    return_status = highs.setHighsOptionValue("simplex_iteration_limit",
                                              use_simplex_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);

    return_status = highs.setBasis();
    REQUIRE(return_status == HighsStatus::OK);
  } else {
    ipm_iteration_count = info.ipm_iteration_count;
    use_ipm_iteration_limit = ipm_iteration_count;

    return_status = highs.setHighsOptionValue("ipm_iteration_limit",
                                              use_ipm_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
  }

  return_status = highs.run();
  model_status = highs.getModelStatus();
  printf("Returns status = %d; model status = %s\n", (int)return_status,
         highs.highsModelStatusToString(model_status).c_str());
  REQUIRE(return_status == HighsStatus::Warning);
  REQUIRE(model_status == HighsModelStatus::REACHED_ITERATION_LIMIT);

  if (use_simplex) {
    REQUIRE(simplex_iteration_count == info.simplex_iteration_count);
    simplex_iteration_count = info.simplex_iteration_count;
  } else {
    REQUIRE(ipm_iteration_count == info.ipm_iteration_count);
    ipm_iteration_count = info.ipm_iteration_count;
  }

  // Now check that it stops after 10 iterations
  if (use_simplex) {
    const int further_simplex_iterations = 10;
    use_simplex_iteration_limit =
        simplex_iteration_count + further_simplex_iterations;
    printf("Setting simplex_iteration_limit = %d\n",
           use_simplex_iteration_limit);
    return_status = highs.setHighsOptionValue("simplex_iteration_limit",
                                              use_simplex_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
    return_status = highs.setBasis();
    REQUIRE(return_status == HighsStatus::OK);
  } else {
    const int further_ipm_iterations = 5;
    use_ipm_iteration_limit = ipm_iteration_count + further_ipm_iterations;
    printf("Setting ipm_iteration_limit = %d\n", use_ipm_iteration_limit);
    return_status = highs.setHighsOptionValue("ipm_iteration_limit",
                                              use_ipm_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
  }

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::Warning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::REACHED_ITERATION_LIMIT);

  if (use_simplex) {
    simplex_iteration_count = use_simplex_iteration_limit;
    REQUIRE(simplex_iteration_count == info.simplex_iteration_count);
    simplex_iteration_count = info.simplex_iteration_count;
    return_status = highs.setHighsOptionValue("simplex_iteration_limit",
                                              default_simplex_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
  } else {
    ipm_iteration_count = use_ipm_iteration_limit;
    REQUIRE(ipm_iteration_count == info.ipm_iteration_count);
    ipm_iteration_count = info.ipm_iteration_count;
    return_status = highs.setHighsOptionValue("ipm_iteration_limit",
                                              default_ipm_iteration_limit);
    REQUIRE(return_status == HighsStatus::OK);
  }
}

void testSolversSetup(const std::string model,
                      IterationCount& model_iteration_count,
                      vector<int>& simplex_strategy_iteration_count) {
  if (model.compare("adlittle") == 0) {
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::SIMPLEX_STRATEGY_CHOOSE] = 86;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::SIMPLEX_STRATEGY_DUAL_PLAIN] = 86;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::SIMPLEX_STRATEGY_DUAL_TASKS] = 86;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::SIMPLEX_STRATEGY_DUAL_MULTI] = 89;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::SIMPLEX_STRATEGY_PRIMAL] = 101;
    model_iteration_count.ipm = 13;  // 14;
    model_iteration_count.crossover = 0;
  }
}

void testSolvers(Highs& highs, IterationCount& model_iteration_count,
                 const vector<int>& simplex_strategy_iteration_count) {
  bool have_omp = false;
#ifdef OPENMP
  have_omp = true;
#endif
  /*
  int i = (int)SimplexStrategy::SIMPLEX_STRATEGY_PRIMAL;
  model_iteration_count.simplex = simplex_strategy_iteration_count[i];
  testSolver(highs, "simplex", model_iteration_count, i);
  */

  int from_i = (int)SimplexStrategy::SIMPLEX_STRATEGY_MIN;
  int to_i = (int)SimplexStrategy::SIMPLEX_STRATEGY_NUM;
  for (int i = from_i; i < to_i; i++) {
    if (!have_omp) {
      if (i == (int)SimplexStrategy::SIMPLEX_STRATEGY_DUAL_TASKS) continue;
      if (i == (int)SimplexStrategy::SIMPLEX_STRATEGY_DUAL_MULTI) continue;
    }
    model_iteration_count.simplex = simplex_strategy_iteration_count[i];
    testSolver(highs, "simplex", model_iteration_count, i);
  }
  testSolver(highs, "ipm", model_iteration_count);
}

// No commas in test case name.
TEST_CASE("LP-solver", "[highs_lp_solver]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string model = "";
  std::string model_file;
  IterationCount model_iteration_count;
  vector<int> simplex_strategy_iteration_count;
  simplex_strategy_iteration_count.resize(
      (int)SimplexStrategy::SIMPLEX_STRATEGY_NUM);

  HighsOptions options;
  HighsLp lp;
  //  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  Highs highs(options);

  // Read mps
  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  testSolversSetup(model, model_iteration_count,
                   simplex_strategy_iteration_count);

  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.setBasis();
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  testSolvers(highs, model_iteration_count, simplex_strategy_iteration_count);

  // Now check that we can change model within the same Highs instance
  // First reset all the options to their default values
  return_status = highs.resetHighsOptions();
  REQUIRE(return_status == HighsStatus::OK);

  model_file = std::string(HIGHS_DIR) + "/check/instances/etamacro.mps";
  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.setBasis();
  REQUIRE(return_status == HighsStatus::OK);

  int simplex_iteration_count;
  return_status = highs.getHighsInfoValue("simplex_iteration_count",
                                          simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  const HighsInfo& info = highs.getHighsInfo();
  REQUIRE(info.num_dual_infeasibilities == 1);

  simplex_iteration_count += 458;
  REQUIRE(simplex_iteration_count == info.simplex_iteration_count);

  HighsModelStatus model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::NOTSET);

  model_status = highs.getModelStatus(true);
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);
}

TEST_CASE("dual-objective-upper-bound", "[highs_lp_solver]") {
  std::string filename;
  HighsStatus status;
  HighsModelStatus model_status;
  bool bool_status;
  const double min_objective_function_value = -11.6389290663705;
  const double max_objective_function_value = 111.650960689315;
  const double smaller_min_dual_objective_value_upper_bound = -110.0;
  const double larger_min_dual_objective_value_upper_bound = -45.876;
  const double use_max_dual_objective_value_upper_bound = 150.0;
  double save_dual_objective_value_upper_bound;
  HighsOptions options;
  Highs highs(options);
  const HighsInfo& info = highs.getHighsInfo();

  //  status = highs.setHighsOptionValue("message_level", 7);REQUIRE(status ==
  //  HighsStatus::OK);

  double error;
  filename = std::string(HIGHS_DIR) + "/check/instances/e226.mps";
  status = highs.setHighsOptionValue("model_file", filename);
  REQUIRE(status == HighsStatus::OK);

  // Solve vanilla
  printf("\nSolving vanilla LP\n");
  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  error = fabs((info.objective_function_value - min_objective_function_value) /
               min_objective_function_value);
  printf("\nOptimal objective value error = %g\n", error);
  REQUIRE(error < 1e-14);

  // Set dual objective value upper bound after saving the default value
  status = highs.getHighsOptionValue("dual_objective_value_upper_bound",
                                     save_dual_objective_value_upper_bound);
  REQUIRE(status == HighsStatus::OK);

  status =
      highs.setHighsOptionValue("dual_objective_value_upper_bound",
                                larger_min_dual_objective_value_upper_bound);
  REQUIRE(status == HighsStatus::OK);

  // Solve again
  printf(
      "\nSolving LP with presolve and dual objective value upper bound of %g\n",
      larger_min_dual_objective_value_upper_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  // Switch off presolve
  status = highs.setHighsOptionValue("presolve", "off");
  REQUIRE(status == HighsStatus::OK);

  // Solve again
  // This larger dual objective value upper bound is satisfied during phase 2
  printf(
      "\nSolving LP without presolve and larger dual objective value upper "
      "bound of %g\n",
      larger_min_dual_objective_value_upper_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status ==
          HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);

  // Solve again
  // This smaller dual objective value upper bound is satisfied at the start of
  // phase 2
  printf(
      "\nSolving LP without presolve and smaller dual objective value upper "
      "bound of %g\n",
      smaller_min_dual_objective_value_upper_bound);
  status =
      highs.setHighsOptionValue("dual_objective_value_upper_bound",
                                smaller_min_dual_objective_value_upper_bound);
  REQUIRE(status == HighsStatus::OK);

  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status ==
          HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);

  // Solve as maximization and ensure that the dual objective value upper bound
  // isn't used
  bool_status = highs.changeObjectiveSense(ObjSense::MAXIMIZE);
  REQUIRE(bool_status);

  status = highs.setHighsOptionValue("dual_objective_value_upper_bound",
                                     use_max_dual_objective_value_upper_bound);
  REQUIRE(status == HighsStatus::OK);

  // Solve again
  printf(
      "\nSolving LP as maximization without presolve and dual objective value "
      "upper bound of %g\n",
      use_max_dual_objective_value_upper_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  error = fabs((info.objective_function_value - max_objective_function_value) /
               max_objective_function_value);
  printf("\nOptimal objective value error = %g\n", error);
  REQUIRE(error < 1e-14);
}
